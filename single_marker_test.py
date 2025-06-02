"""
single_marker_test.py

This module implements the Single Marker Test (SMT) approach for 
QTL detection. It evaluates the association between individual genetic 
markers and the trait of interest using simple statistical tests.

Main Responsibilities:
- Perform SMT on input genotype and phenotype data.
- Compute test statistics (e.g., t-tests or ANOVA).
- Report significant markers and their effect estimates.
"""


import os
from vcf_data_handler import parse_vcf_file, generate_genetic_map,parse_trait_file
from plot_utils import plot_genetic_map, plot_genetic_distance_circles,plot_full_recombination_dots
import matplotlib.pyplot as plt
import pandas as pd

# Set the path to the VCF file (relative to this script)
vcf_path = os.path.join("Data", "cataglyphis.final.DZ (1).vcf")

# Parse the VCF file
vcf_data, sample_names = parse_vcf_file(vcf_path)

# Print number of chromosomes and first 5 chromosome IDs
chrom_ids = list(vcf_data.keys())

genetic_maps = generate_genetic_map(vcf_data, sample_names)

trait_path = os.path.join("Data", "traits.txt")
sample_names, traits = parse_trait_file(trait_path)
print("Available traits:", list(traits.keys()))


import numpy as np

def run_smt_regression_and_export(vcf_data, traits, trait_name, sample_names, output_csv="smt_results.csv"):
    results = []

    for chrom_id, chrom_data in vcf_data.items():
        for marker in chrom_data["markers"]:
            pos = marker["pos"]
            G = []
            T = []

            for name in sample_names:
                gt_str = marker["samples"].get(name, {}).get("GT")
                trait_val = traits[trait_name].get(name)

                if gt_str is None or trait_val is None:
                    continue

                if gt_str == "0/0":
                    gt_code = 0
                elif gt_str in {"0/1", "1/0"}:
                    gt_code = 1
                elif gt_str == "1/1":
                    gt_code = 2
                else:
                    continue

                # Two Vectors, G: Genotype of individual i
                G.append(gt_code)
                # T -> Trait value for individual i
                T.append(trait_val)

            if len(G) < 3:
                continue

            G = np.array(G)
            T = np.array(T)

            # Average of all genotypes, Example: G = [0,2,1,2,0,...] -> G-bar = average of G.
            G_bar = np.mean(G)
            # Average of all phenotypes (traits I check)
            T_bar = np.mean(T)

            # B beta calculation
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

            results.append({
                "chrom": chrom_id,
                "pos": pos,
                "beta": beta_hat,
                "F": F_stat,
                "n": n
            })

    df = pd.DataFrame(results)
    df.to_csv(output_csv, index=False)
    #print(f"Saved SMT results to {output_csv}")

    #unique_chroms = df['chrom'].unique()
    #for chrom in unique_chroms:
    chrom = "chr01"
    chrom_data = df[df['chrom'] == chrom]
    chrom_data = df[df['chrom'] == chrom]
    plt.figure()
    plt.plot(chrom_data['pos'], chrom_data['F'], marker='o')
    plt.title(f"F-statistic across chromosome {chrom}")
    plt.xlabel("Marker Position")
    plt.ylabel("F-statistic")
    plt.grid(True)
    plt.tight_layout()
    plt.show()

    return df

# שים כאן את שם התכונה שלך מתוך traits.txt
trait_to_test = "TC25Me3"  # שנה לשם האמיתי

run_smt_regression_and_export(vcf_data, traits, trait_to_test, sample_names)
