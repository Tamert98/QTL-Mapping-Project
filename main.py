"""
main.py

This is the main script for running the QTL mapping pipeline.
It orchestrates the full analysis by loading the input VCF file,
applying the SMT method, and displaying the results.

Main Responsibilities:
- Load genotype and trait data.
- Generate genetic maps (filtered and unfiltered).
- Run Single Marker Test (SMT) and compute statistics.
- Plot the genetic map and SMT results.
"""

import os
import matplotlib.pyplot as plt
from vcf_data_handler import (
    parse_vcf_file,
    generate_genetic_map,
    generate_genetic_map_filtered,
    parse_trait_file,
)
from plot_utils import (
    plot_genetic_map,
    plot_f_statistic,
    plot_p_values,
    plot_log_p_values,
    plot_combined_qtl_significance,
)
from single_marker_test import run_smt_regression_and_export,find_best_qtl

def main():
    # Input paths
    vcf_path = os.path.join("Data", "cataglyphis.final.DZ (1).vcf")
    trait_path = os.path.join("Data", "traits.txt")
    trait_to_test = "TC25Me3"
    chr_id = "chr03"

    # Load VCF data and trait values
    vcf_data, sample_names = parse_vcf_file(vcf_path)
    sample_names, traits = parse_trait_file(trait_path)

    # Generate genetic maps
    genetic_maps_unfiltered = generate_genetic_map(vcf_data, sample_names)
    genetic_maps_filtered = generate_genetic_map_filtered(vcf_data, sample_names)

    # Run Single Marker Test and save results
    df = run_smt_regression_and_export(vcf_data, traits, trait_to_test, sample_names)

    # Plot results for a specific chromosome
    fig1 = plot_genetic_map(genetic_maps_unfiltered[chr_id]["genetic_map"], chr_id)
    fig2 = plot_genetic_map(genetic_maps_filtered[chr_id]["genetic_map"], chr_id)
    fig3 = plot_f_statistic(df, chr_id)
    fig4 = plot_log_p_values(df, chr_id)

    # Plot overall QTL significance
    fig5 = plot_combined_qtl_significance(df)

    # Find most significant QTL in BP and cM
    smtQtl_Bp, smtQtl_Cm = find_best_qtl(df, genetic_maps_filtered)
    print(f"Most significant QTL: BP = {smtQtl_Bp}, cM = {smtQtl_Cm}")

    # Show all figures
    plt.show()

if __name__ == "__main__":
    main()
