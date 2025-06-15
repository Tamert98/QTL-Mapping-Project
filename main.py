"""
main.py

Main script for running the full QTL mapping pipeline.

Responsibilities:
- Load genotype and trait data.
- Generate genetic maps (filtered and unfiltered).
- Generate pairwise recombination maps.
- Compute global axis limits (for unified plots).
- Run Single Marker Test (SMT) and export QTL results.
- Generate and export:
    * Genetic maps (PDF)
    * Filtered genetic maps (PDF)
    * Heatmaps (PDF)
    * SMT significance plots (PDF)
"""

import os
from vcf_data_handler import (
    parse_vcf_file,
    generate_genetic_map,
    generate_genetic_map_filtered,
    generate_pairwise_genetic_map,
    parse_trait_file,
)
from plot_utils import (
    compute_global_physical_max,
    generate_genetic_map_images_and_pdf,
    generate_filtered_genetic_map_images_and_pdf,
    generate_heatmap_images_and_pdf,
    generate_smt_chromosome_images,
    stitch_smt_chromosome_images_pdf,
    plot_concatenated_qtl_significance,
    generate_concatenated_qtl_pdf,
)
from single_marker_test import (
    run_smt_regression_and_export,
    find_best_qtl,
    run_smt_for_all_traits
)

def main():
    # Input paths and parameters
    vcf_path = os.path.join("Data", "cataglyphis.final.DZ (1).vcf")
    trait_path = os.path.join("Data", "traits.txt")
    trait_to_test = "TC25Me3"
    pvalue_cutoff = 0.05

    # Load data
    vcf_data, sample_names = parse_vcf_file(vcf_path)
    sample_names, traits = parse_trait_file(trait_path)

    # Generate genetic maps
    genetic_maps_unfiltered = generate_genetic_map(vcf_data, sample_names)
    genetic_maps_filtered = generate_genetic_map_filtered(vcf_data, sample_names)

    # Compute global x-axis based on unfiltered maps
    global_xmax = compute_global_physical_max(genetic_maps_unfiltered)

    # Export Genetic Maps with unified X axis
    generate_genetic_map_images_and_pdf(genetic_maps_unfiltered, global_xmax)
    generate_filtered_genetic_map_images_and_pdf(genetic_maps_filtered, global_xmax)

    # Generate pairwise recombination maps for heatmaps
    pairwise_maps = {}
    for chrom_id in genetic_maps_unfiltered:
        pairwise = generate_pairwise_genetic_map(vcf_data, sample_names, chrom_id)
        if pairwise:
            pairwise_maps[chrom_id] = pairwise

    # Export Heatmaps with unified X/Y axis
    generate_heatmap_images_and_pdf(pairwise_maps)

    # Run SMT for selected trait
    df = run_smt_regression_and_export(vcf_data, traits, trait_to_test, sample_names)

    # Generate per-chromosome SMT significance images and combine to PDF
    generate_smt_chromosome_images(df, trait_to_test, significance_level=pvalue_cutoff)
    final_pdf_path = stitch_smt_chromosome_images_pdf(trait_to_test)

    # Plot genome-wide concatenated QTL map for selected trait
    plot_concatenated_qtl_significance(df, trait_to_test, significance_level=pvalue_cutoff)

    # Report best QTL location found
    smtQtl_Bp, smtQtl_Cm = find_best_qtl(df, genetic_maps_unfiltered)
    print(f"Most significant QTL: BP = {smtQtl_Bp}, cM = {smtQtl_Cm}")
    print(f"SMT PDF generated: {final_pdf_path}")

    # Run SMT for all traits and export combined QTL maps
    all_smt_results = run_smt_for_all_traits(vcf_data, traits, sample_names)
    generate_concatenated_qtl_pdf(all_smt_results)

if __name__ == "__main__":
    main()
