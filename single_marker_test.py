"""
single_marker_test.py

Implements the Single Marker Test (SMT) for QTL detection.
This method evaluates the association between each genetic marker
and the trait using a linear regression approach.

Responsibilities:
- Calculate F-statistic, p-value, and LOD score per marker.
- Save detailed results to CSV.
- Identify best QTL marker based on statistical strength.
"""

import os
import numpy as np
import pandas as pd
from scipy.stats import f

# === Core Regression ===
def run_smt_regression_and_export(vcf_data, traits, trait_name, sample_names, output_csv="smt_results.csv"):
    """
    Runs regression for a single trait across all markers.

    Args:
        vcf_data (dict): VCF marker data.
        traits (dict): Trait values per sample.
        trait_name (str): Name of the trait to analyze.
        sample_names (list): List of sample IDs.
        output_csv (str): Output path for CSV.

    Returns:
        pd.DataFrame: Regression results with F, p, LOD, beta, etc.
    """
    results = []

    for chrom_id, chrom_data in vcf_data.items():
        for marker in chrom_data["markers"]:
            pos = marker["pos"]
            G, T = [], []

            for name in sample_names:
                gt_str = marker["samples"].get(name, {}).get("GT")
                trait_val = traits[trait_name].get(name)

                if gt_str is None or trait_val is None:
                    continue

                # Genotype encoding
                if gt_str == "0/0":
                    gt_code = 0
                elif gt_str in {"0/1", "1/0"}:
                    gt_code = 1
                elif gt_str == "1/1":
                    gt_code = 2
                else:
                    continue

                G.append(gt_code)
                T.append(trait_val)

            if len(G) < 3:
                continue  # Skip sparse data

            G = np.array(G)
            T = np.array(T)

            G_bar, T_bar = np.mean(G), np.mean(T)

            # Linear regression
            numerator = np.sum((G - G_bar) * (T - T_bar))
            denominator = np.sum((G - G_bar) ** 2)
            if denominator == 0:
                continue

            beta_hat = numerator / denominator
            T_hat = beta_hat * (G - G_bar) + T_bar

            SS_reg = np.sum((T_hat - T_bar) ** 2)
            SS_err = np.sum((T - T_hat) ** 2)
            n = len(G)

            if n <= 2:
                continue

            F_stat = (SS_reg / 1) / (SS_err / (n - 2))
            p_val = f.sf(F_stat, 1, n - 2)
            lod = -np.log10(p_val) if p_val > 0 else np.inf

            results.append({
                "chrom": chrom_id,
                "pos": pos,
                "beta": beta_hat,
                "F": F_stat,
                "p": p_val,
                "LOD": lod,
                "n": n
            })

    df = pd.DataFrame(results)
    df.to_csv(output_csv, index=False)
    print(f"Saved SMT results to {output_csv}")
    return df

# === Run SMT on All Traits ===
def run_smt_for_all_traits(vcf_data, traits, sample_names, output_dir="Results/SMT/Reports", message_callback=None):
    """
    Runs SMT regression for all traits in the dataset.

    Args:
        vcf_data (dict): Parsed VCF genotype data.
        traits (dict): Parsed traits {trait_name: {sample: value}}.
        sample_names (list): Sample IDs from VCF.
        output_dir (str): Directory where SMT CSV files will be saved.
        message_callback (callable, optional): Function to call with status messages.

    Returns:
        dict: {trait_name: DataFrame of SMT results}
    """
    os.makedirs(output_dir, exist_ok=True)
    results_per_trait = {}

    def send_message(msg):
        if message_callback:
            message_callback(msg)
        else:
            print(msg)

    for idx, trait_name in enumerate(traits.keys(), 1):
        send_message(f"Running SMT for trait: {trait_name}, trait {idx}/{len(traits)}")

        output_csv = os.path.join(output_dir, f"smt_results_{trait_name}.csv")
        df = run_smt_regression_and_export(vcf_data, traits, trait_name, sample_names, output_csv=output_csv)

        send_message(f"Saved SMT results to {output_csv}")
        results_per_trait[trait_name] = df

    return results_per_trait

# === Best Marker Selection ===
def get_best_marker_info(trait_name, results_per_trait, genetic_maps):
    """
    Returns the location (cM, bp) and -log10(p) of the best marker for a given trait.

    Args:
        trait_name (str): The trait to look up.
        results_per_trait (dict): {trait_name: DataFrame of SMT results}.
        genetic_maps (dict): Genetic maps for all chromosomes.

    Returns:
        tuple: (chrom, cM, position_bp, -log10(p))
    """
    df = results_per_trait.get(trait_name)
    if df is None or df.empty:
        return None, None, None, None

    if "-log10(p)" not in df.columns:
        df["-log10(p)"] = -np.log10(df["p"])

    best_row = df.loc[df["-log10(p)"].idxmax()]
    best_chr = best_row["chrom"]
    best_bp = best_row["pos"]
    best_logp = best_row["-log10(p)"]

    best_cM = None
    if best_chr in genetic_maps:
        for marker in genetic_maps[best_chr]["genetic_map"]:
            if marker["pos"] == best_bp:
                best_cM = marker["cM"]
                break

    return best_chr, best_cM, best_bp, best_logp
