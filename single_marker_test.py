"""
single_marker_test.py

Implements the Single Marker Test (SMT) for QTL detection.
This method evaluates the association between each genetic marker
and the trait using a linear regression approach.

Main Responsibilities:
- Calculate F-statistic per marker.
- Derive p-value and LOD score.
- Save results to CSV.
"""

import numpy as np
import pandas as pd
from scipy.stats import f

def run_smt_regression_and_export(vcf_data, traits, trait_name, sample_names, output_csv="smt_results.csv"):
    results = []

    for chrom_id, chrom_data in vcf_data.items():
        for marker in chrom_data["markers"]:
            pos = marker["pos"]
            G = []  # Genotypes
            T = []  # Traits

            for name in sample_names:
                gt_str = marker["samples"].get(name, {}).get("GT")
                trait_val = traits[trait_name].get(name)

                if gt_str is None or trait_val is None:
                    continue

                # Encode genotype as numeric
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
                continue  # Skip if insufficient data

            G = np.array(G)
            T = np.array(T)

            G_bar = np.mean(G)
            T_bar = np.mean(T)

            # Simple linear regression: Y = βX + ε
            numerator = np.sum((G - G_bar) * (T - T_bar))
            denominator = np.sum((G - G_bar) ** 2)
            if denominator == 0:
                continue

            beta_hat = numerator / denominator
            T_hat = beta_hat * (G - G_bar) + T_bar

            # Compute regression and error sum of squares
            SS_reg = np.sum((T_hat - T_bar) ** 2)
            SS_err = np.sum((T - T_hat) ** 2)

            n = len(G)
            if n <= 2:
                continue

            # F-statistic and p-value
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

def find_best_qtl(df, genetic_maps):
    df["-log10(p)"] = -np.log10(df["p"])
    best_row = df.loc[df["-log10(p)"].idxmax()]
    best_chr = best_row["chrom"]
    best_bp = best_row["pos"]

    best_cM = None
    if best_chr in genetic_maps:
        gmap = genetic_maps[best_chr]["genetic_map"]
        for m in gmap:
            if m["pos"] == best_bp:
                best_cM = m["cM"]
                break

    return best_bp, best_cM
