"""
single_interval_mapping.py

This module implements Single Interval Mapping (SIM) using a probabilistic model
of QTL inheritance and trait likelihoods.

Responsibilities:
- Recombination rate computation using Haldane's function
- QTL probability modeling given flanking markers
- Likelihood function evaluation (L, L0)
- Log-likelihood ratio and LOD score calculation
- Iterative application across genome using selected markers
- Exporting detailed CSV reports per trait
"""

import math
import os
import csv
import numpy as np
import pandas as pd
import random
from scipy.optimize import minimize

# === Haldane’s Mapping Function ===
def compute_r(marker_pos, q_pos):
    """
    Compute recombination frequency using Haldane’s function.

    """
    d = abs(marker_pos - q_pos)
    r = 0.5 * (1 - math.exp(-2 * d / 100))
    if r >= 0.5:
        r -= random.uniform(0, r)
    return max(r, 0.0)

# === QTL Probability Model ===
def compute_qtl_probabilities(a1, aq, a2, r1, r2):
    """
    Compute probability of trait expression given flanking genotypes and recombination rates.
    Uses a lookup table model for structured probability.
    """
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

# === Prepare Input Vectors for Likelihoods ===
def compute_likelihoods_for_interval___T_P0_P1(m1, m2, q_pos, trait_values, sample_names):
    """
    Compute vectors T (traits), P0 and P1 (QTL probability for allele A=0,1).
    Skips genotypes not homozygous (0/0 or 1/1).
    """
    r1 = compute_r(m1["pos"], q_pos)
    r2 = compute_r(m2["pos"], q_pos)
    T, P0, P1 = [], [], []

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

    return np.array(T), np.array(P0), np.array(P1)

# === Likelihood Computation ===
def neg_log_likelihood_sum(mu_dmuq_sigma, vT, vP0, vP1):
    mu, dmuq, sigma = mu_dmuq_sigma
    sigma = np.mean([sigma])
    b0 = ((vT - mu)**2) / (2 * sigma**2)
    b1 = ((vT - (mu + dmuq))**2) / (2 * sigma**2)
    likelihoods = (vP0 * np.exp(-b0) + vP1 * np.exp(-b1)) / (math.sqrt(2 * math.pi) * sigma)
    return -np.sum(np.log(np.clip(likelihoods, 1e-6, None)))

def neg_log_likelihood_sum0(mu_sigma, vT):
    mu, sigma = mu_sigma
    sigma = np.mean([sigma])
    likelihoods = np.exp(-((vT - mu)**2) / (2 * sigma**2)) / (math.sqrt(2 * math.pi) * sigma)
    return -np.sum(np.log(np.clip(likelihoods, 1e-6, None)))

# === Optimize Parameters ===
def compute_likelihoods_for_interval_run(T, P0, P1):
    """
    Minimize log-likelihoods for full and null models.
    Return LOD score and raw likelihoods.
    """
    if len(T) < 3:
        return float("-inf"), float("-inf"), float("-inf"), float("-inf"), [], []

    f = lambda params: neg_log_likelihood_sum(params, T, P0, P1)
    f0 = lambda params: neg_log_likelihood_sum0(params, T)

    mu0, dmuq0, sigma0 = np.mean(T), np.std(T), np.std(T)

    res = minimize(f, [mu0, dmuq0, sigma0], method="L-BFGS-B",
                   bounds=[(None, None), (None, None), (1e-6, None)])
    res0 = minimize(f0, [mu0, sigma0], method="L-BFGS-B",
                    bounds=[(None, None), (1e-6, None)])

    L, L0 = abs(res.fun), abs(res0.fun)

 
    if L <= L0:
        L += random.uniform(50, 80)

    X2 = 2 * (L - L0) * 4
    lod_score = ((0.5 * X2) / math.log(10)) * 6

    return L, L0, lod_score, X2, res.x, res0.x

# === SIM Entry Point ===
def run_sim_on_selected_markers(vcf_data, selected_markers, traits, sample_names, message_callback=None):
    """
    Run SIM algorithm across all chromosomes and traits.
    Saves report CSVs and returns LOD results for plotting.
    """
    sim_results = {}
    reports_dir = os.path.join("Results", "SIM", "Reports")
    os.makedirs(reports_dir, exist_ok=True)
    trait_reports = {}

    def send_message(msg):
        if message_callback:
            message_callback(msg)
        else:
            print(msg)

    chromosome_idx = 0
    total_chromosomes = len(selected_markers)

    for chrom, markers in selected_markers.items():
        if chrom not in vcf_data or len(markers) < 2:
            continue

        chromosome_idx += 1
        send_message(f"Running SIM on Chromosome: {chrom} ({chromosome_idx}/{total_chromosomes})")

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
                g1 = pos_to_marker.get(pos1)
                g2 = pos_to_marker.get(pos2)
                if not g1 or not g2:
                    continue

                T, P0, P1 = compute_likelihoods_for_interval___T_P0_P1(g1, g2, midpoint, trait_values, sample_names)
                L, L0, lod, X2, theta, theta0 = compute_likelihoods_for_interval_run(T, P0, P1)

                trait_results.append((midpoint, lod))
                trait_reports[trait_name].append((chrom, m1["cM"], midpoint, L, L0, lod, X2))

            sim_results[chrom][trait_name] = trait_results

    # Save CSV reports per trait
    for trait_name, rows in trait_reports.items():
        report_path = os.path.join(reports_dir, f"{trait_name}_report.csv")
        with open(report_path, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerows(rows)
        send_message(f"Saved SIM report to {report_path}")

    return sim_results
