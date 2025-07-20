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
    select_evenly_spaced_markers,
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
    plot_lod_curve,
    plot_all_lod_curves
)
from single_marker_test import (
    run_smt_regression_and_export,
    find_best_qtl,
    run_smt_for_all_traits,
    get_best_marker_info
)

from single_interval_mapping import run_sim_on_selected_markers

def main():
    # Input paths and parameters
    vcf_path = os.path.join("Data", "cataglyphis.final.DZ (1).vcf")
    trait_path = os.path.join("Data", "traits.txt")

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

    # Generate evenly spaced markers from filtered genetic maps (3–5 cM apart)
    selected_markers = select_evenly_spaced_markers(genetic_maps_filtered, min_cm_spacing=1.0)

    # Generate pairwise recombination maps for heatmaps
    """pairwise_maps = {}
    for chrom_id in genetic_maps_filtered:
        pairwise = generate_pairwise_genetic_map(vcf_data, sample_names, chrom_id)
        if pairwise:
            pairwise_maps[chrom_id] = pairwise

    # Export Heatmaps with unified X/Y axis
    generate_heatmap_images_and_pdf(pairwise_maps)"""

    # Run SMT for selected trait
    #df = run_smt_regression_and_export(vcf_data, traits, trait_to_test, sample_names)

    # Generate per-chromosome SMT significance images and combine to PDF
    #generate_smt_chromosome_images(df, trait_to_test, significance_level=pvalue_cutoff)
    #final_pdf_path = stitch_smt_chromosome_images_pdf(trait_to_test)

    # Run SMT for all traits and export combined QTL maps
    all_smt_results = run_smt_for_all_traits(vcf_data, traits, sample_names)
    generate_concatenated_qtl_pdf(all_smt_results)
    

    #running the sim algorithim and getting the L,L0,q,LOD score in a list for each tarit 
    """sim_results = run_sim_on_selected_markers(vcf_data, selected_markers, traits, sample_names)
    for chrom, chrom_results in sim_results.items():
        print(f"Chromosome {chrom}:")
        for trait_name, results in chrom_results.items():
            print(f"  Trait: {trait_name}")
            for entry in results[:3]:  # show only first 3 entries for readability
                print(f"    {entry}")
            break  # show only first trait for each chromosome
    # Assuming sim_results already computed
    plot_lod_curve(sim_results, chromosome_id='chr01')
    plot_all_lod_curves(sim_results)"""
    
if __name__ == "__main__":
    main()
