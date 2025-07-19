import math
import os
import csv
import numpy as np
from scipy.optimize import minimize

# --- Utility: Compute recombination frequency between marker and QTL position ---
def compute_r(marker, q_pos):
    d = abs(marker["pos"] - q_pos)
    r = 0.5 * (1 - math.exp(-2 * d / 100))  # Haldane's mapping function
    return r

# --- Compute probability QTL inherits allele A ---
def compute_qtl_probabilities(a1, aq, a2, r1, r2):
    prob_table_h1q2 = {
        (0, 0, 0): 0.5 * (1 - r1) * (1 - r2),
        (1, 1, 1): 0.5 * (1 - r1) * (1 - r2),
        (0, 0, 1): 0.5 * (1 - r1) * r2,
        (1, 1, 0): 0.5 * (1 - r1) * r2,
        (0, 1, 1): 0.5 * r1 * (1 - r2),
        (1, 0, 0): 0.5 * r1 * (1 - r2),
        (0, 1, 0): 0.5 * r1 * r2,
        (1, 0, 1): 0.5 * r1 * r2,
    }
    prob_table_h12 = {
        (0, 0): 0.5 * (1 - r1) * (1 - r2) + 0.5 * r1 * r2,
        (1, 1): 0.5 * (1 - r1) * (1 - r2) + 0.5 * r1 * r2,
        (0, 1): 0.5 * (1 - r1) * r2 + 0.5 * (1 - r2) * r1,
        (1, 0): 0.5 * (1 - r1) * r2 + 0.5 * (1 - r2) * r1,
    }
    pc = prob_table_h12.get((a1, a2), 0)
    if pc == 0:
        return 0
    pu = prob_table_h1q2.get((a1, aq, a2), 0)
    return pu / pc

# --- Main likelihood computation for a given interval and trait ---
def compute_likelihoods_for_interval(m1, m2, q_pos, trait_values, sample_names):
    r1 = compute_r(m1, q_pos)
    r2 = compute_r(m2, q_pos)

    T = []
    P0 = []
    P1 = []
    for name, trait in zip(sample_names, trait_values):
        gt1 = m1["samples"].get(name, {}).get("GT", "./.")
        gt2 = m2["samples"].get(name, {}).get("GT", "./.")
        a1 = 0 if gt1 == "0/0" else 1 if gt1 == "1/1" else -1
        a2 = 0 if gt2 == "0/0" else 1 if gt2 == "1/1" else -1
        if trait is None or min(a1, a2) < 0:
            continue
        prob0 = compute_qtl_probabilities(a1, 0, a2, r1, r2)
        prob1 = compute_qtl_probabilities(a1, 1, a2, r1, r2)
        T.append(trait)
        P0.append(prob0)
        P1.append(prob1)

    if len(T) < 3:
        return float("-inf"), float("-inf"), float("-inf"), float("-inf")

    T = np.array(T)
    P0 = np.array(P0)
    P1 = np.array(P1)

    def neg_log_likelihood_sum(mu_dmuq_sigma):
        mu, dmuq, sigma = mu_dmuq_sigma
        b0 = ((T - mu) ** 2) / (2 * sigma ** 2)
        b1 = ((T - (mu + dmuq)) ** 2) / (2 * sigma ** 2)
        likelihoods = P0 * np.exp(-b0) + P1 * np.exp(-b1)
        likelihoods = np.clip(likelihoods, 1e-10, None)
        return -np.sum(np.log(likelihoods))

    def neg_log_likelihood_sum0(sigma_only):
        sigma = sigma_only[0]
        likelihoods = np.exp(-(T ** 2) / (2 * sigma ** 2)) / (math.sqrt(2 * math.pi) * sigma)
        likelihoods = np.clip(likelihoods, 1e-10, None)
        return -np.sum(np.log(likelihoods))

    mu0 = np.mean(T)
    dmuq0 = np.std(T)
    sigma0 = np.std(T)

    res = minimize(neg_log_likelihood_sum, [mu0, dmuq0, sigma0], method="L-BFGS-B",
                   bounds=[(None, None), (None, None), (1e-6, None)])

    res0 = minimize(neg_log_likelihood_sum0, [sigma0], method="L-BFGS-B",
                    bounds=[(1e-6, None)])

    L = -res.fun
    L0 = -res0.fun
    ln_lod = L - L0
    lod_score = ln_lod / math.log(10)
    X2 = 2 * (L0 - L)

    # --- Diagnostics & Sanity Checks ---
    sum_probs = P0 + P1
    if np.any(sum_probs > 1.01) or np.any(sum_probs < 0.99):
        print(f"[WARNING] QTL genotype probabilities not summing to 1 at position {q_pos}")
        print(f"  Min(P0 + P1): {np.min(sum_probs):.4f}, Max(P0 + P1): {np.max(sum_probs):.4f}")

    if L < L0:
        print(f"[NOTE] L < L0 at position {q_pos}, QTL model worse than null.")
        print(f"  L (model): {L:.4f}, L0 (null): {L0:.4f}, LOD: {lod_score:.4f}")

    if not np.isfinite(L) or not np.isfinite(L0):
        print(f"[ERROR] Non-finite likelihood at position {q_pos}")
        print(f"  L: {L}, L0: {L0}, Sample count: {len(T)}")

    if L > 1e4 or L0 > 1e4:
        print(f"[INFO] Very high likelihood at position {q_pos}")
        print(f"  L: {L:.2f}, L0: {L0:.2f}, LOD: {lod_score:.2f}")

    return L, L0, lod_score, X2

# --- Entry point: run SIM on all traits and all chromosomes ---
def run_sim_on_selected_markers(vcf_data, selected_markers, traits, sample_names):
    sim_results = {}
    reports_dir = os.path.join("Results", "SIM", "Reports")
    os.makedirs(reports_dir, exist_ok=True)
    trait_reports = {}

    for chrom, markers in selected_markers.items():
        if chrom not in vcf_data or len(markers) < 2:
            continue

        sim_results[chrom] = {}
        chrom_markers = vcf_data[chrom]["markers"]
        pos_to_marker = {m["pos"]: m for m in chrom_markers}

        for trait_name, trait_dict in traits.items():
            trait_values = [trait_dict.get(name) for name in sample_names]
            trait_results = []

            if trait_name not in trait_reports:
                trait_reports[trait_name] = [("Chromosome", "cM", "Position_bp", "L", "L0", "LOD", "X^2")]

            for i in range(len(markers) - 1):
                m1, m2 = markers[i], markers[i + 1]
                pos1, pos2 = m1["pos"], m2["pos"]
                midpoint = (pos1 + pos2) // 2
                g1, g2 = pos_to_marker.get(pos1), pos_to_marker.get(pos2)
                if not g1 or not g2:
                    continue
                L, L0, lod_score, X2 = compute_likelihoods_for_interval(g1, g2, midpoint, trait_values, sample_names)
                trait_results.append((midpoint, lod_score))
                trait_reports[trait_name].append((chrom, m1["cM"], midpoint, L, L0, lod_score, X2))

            sim_results[chrom][trait_name] = trait_results

    for trait_name, rows in trait_reports.items():
        report_path = os.path.join(reports_dir, f"{trait_name}_report.csv")
        with open(report_path, mode='w', newline='') as f:
            writer = csv.writer(f)
            writer.writerows(rows)

    return sim_results
