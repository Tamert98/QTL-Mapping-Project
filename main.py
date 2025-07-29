"""
main.py

Main script for running the full QTL mapping pipeline.

Responsibilities:
- Load genotype and trait data.
- Generate genetic maps (filtered and unfiltered).
- Generate pairwise recombination heatmaps.
- Compute global axis limits for unified plots.
- Run Single Marker Test (SMT) and Single Interval Mapping (SIM).
- Export results as images and PDFs:
    * Genetic maps
    * Filtered genetic maps
    * Recombination heatmaps
    * SMT significance plots
    * SMT concatenated QTL map
    * SIM LOD curves
"""

# === Standard Library Imports ===
import os
import math

# === Project Module Imports ===
from vcf_data_handler import (
    parse_vcf_file,
    parse_trait_file,
    generate_combined_genetic_maps,
    generate_combined_genetic_maps_filtered,
    select_evenly_spaced_markers,
    generate_pairwise_genetic_map,
)

from plot_utils import (
    compute_global_physical_max,
    generate_genetic_map_images_and_pdf,
    generate_comparative_genetic_map_images_and_pdf,
    generate_heatmap_images_and_pdf,
    generate_smt_chromosome_images_all_traits,
    stitch_all_smt_chromosome_images,
    generate_concatenated_qtl_pdf,
    plot_all_lod_curves,
)

from single_marker_test import (
    run_smt_for_all_traits,
)

from single_interval_mapping import (
    run_sim_on_selected_markers,
)

# === Main Pipeline ===
def main():
    # === Input File Paths ===
    vcf_path = os.path.join("Data", "cataglyphis.final.DZ (1).vcf")
    trait_path = os.path.join("Data", "traits.txt")

    # === Load Input Data ===
    vcf_data, sample_names = parse_vcf_file(vcf_path)
    sample_names, traits = parse_trait_file(trait_path)

    # === Generate Genetic Maps (Unfiltered + Filtered) ===
    genetic_maps_m0, genetic_maps_comparative = generate_combined_genetic_maps(vcf_data, sample_names)
    genetic_maps_m0_filtered, genetic_maps_comparative_filtered = generate_combined_genetic_maps_filtered(vcf_data, sample_names)

    # === Compute Unified X-Axis (based on m0 filtered map) ===
    global_xmax = compute_global_physical_max(genetic_maps_m0_filtered)

    # === Export Genetic Map Images ===
    generate_genetic_map_images_and_pdf(genetic_maps_m0_filtered, global_xmax)
    generate_comparative_genetic_map_images_and_pdf(
        genetic_maps_comparative,
        genetic_maps_comparative_filtered,
        global_xmax
    )

    # === Run SMT (Single Marker Test) ===
    all_smt_results = run_smt_for_all_traits(vcf_data, traits, sample_names)
    generate_smt_chromosome_images_all_traits(all_smt_results)
    stitch_all_smt_chromosome_images(all_smt_results)
    generate_concatenated_qtl_pdf(all_smt_results)

    # === Select Evenly Spaced Markers (for SIM) ===
    selected_markers = select_evenly_spaced_markers(
        genetic_maps_m0_filtered,
        vcf_data,
        min_cm_spacing=1.0  # Adjustable spacing threshold
    )

    # === Generate Pairwise Recombination Maps (Heatmaps) ===
    pairwise_maps = {}
    for chrom_id in genetic_maps_m0_filtered:
        pairwise = generate_pairwise_genetic_map(vcf_data, sample_names, chrom_id)
        if pairwise:
            pairwise_maps[chrom_id] = pairwise

    generate_heatmap_images_and_pdf(pairwise_maps)

    # === Run SIM (Single Interval Mapping) ===
    sim_results = run_sim_on_selected_markers(vcf_data, selected_markers, traits, sample_names)
    plot_all_lod_curves(sim_results)


# === Entry Point ===
if __name__ == "__main__":
    main()
