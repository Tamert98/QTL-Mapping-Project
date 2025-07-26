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
import math
import os
from vcf_data_handler import (
    parse_vcf_file,
    generate_combined_genetic_maps,
    generate_combined_genetic_maps_filtered,
    generate_pairwise_genetic_map,
    parse_trait_file,
    select_evenly_spaced_markers,
)
from plot_utils import (
    compute_global_physical_max,
    generate_genetic_map_images_and_pdf,
    generate_comparative_genetic_map_images_and_pdf,  # <-- ADDED
    generate_smt_chromosome_images_all_traits,
    stitch_all_smt_chromosome_images,
    generate_heatmap_images_and_pdf,
    plot_concatenated_qtl_significance,
    generate_concatenated_qtl_pdf,
    plot_all_lod_curves
)
from single_marker_test import (
    run_smt_regression_and_export,
    find_best_qtl,
    run_smt_for_all_traits,
    get_best_marker_info
)

from single_interval_mapping import run_sim_on_selected_markers

def compute_genetic_map_by_recombination(selected_markers, sample_names):
    for chrom, markers in selected_markers.items():
        cM_positions = [0.0]  # first marker is at 0 cM
        for i in range(len(markers) - 1):
            m1 = markers[i]
            m2 = markers[i + 1]
            r = compute_r_between_markers(m1, m2, sample_names)
            if r is None or r >= 0.5:
                d_cm = 50  # Max cap
            else:
                d_cm = haldane_r_to_cm(r)
            cM_positions.append(cM_positions[-1] + d_cm)
        # Store back
        for m, cm in zip(markers, cM_positions):
            m["cM"] = cm

def compute_r_between_markers(m1, m2, sample_names):
    n = 0
    recomb = 0
    for name in sample_names:
        gt1 = m1["samples"].get(name, {}).get("GT", "./.")
        gt2 = m2["samples"].get(name, {}).get("GT", "./.")
        if gt1 in ("0/0", "1/1") and gt2 in ("0/0", "1/1"):
            a1 = 0 if gt1 == "0/0" else 1
            a2 = 0 if gt2 == "0/0" else 1
            if a1 != a2:
                recomb += 1
            n += 1
    if n == 0:
        return None
    return recomb / n

def haldane_r_to_cm(r):
    if r >= 0.5:
        return float('inf')
    return -50 * math.log(1 - 2 * r)

def main():
    # Input paths and parameters
    vcf_path = os.path.join("Data", "cataglyphis.final.DZ (1).vcf")
    trait_path = os.path.join("Data", "traits.txt")

    # Load data
    vcf_data, sample_names = parse_vcf_file(vcf_path)
    sample_names, traits = parse_trait_file(trait_path)

    # Generate unfiltered and filtered genetic maps
    genetic_maps_m0, genetic_maps_comparative = generate_combined_genetic_maps(vcf_data, sample_names)
    genetic_maps_m0_filtered, genetic_maps_comparative_filtered = generate_combined_genetic_maps_filtered(vcf_data, sample_names)

    # Compute global x-axis based on m0-based maps
    global_xmax = compute_global_physical_max(genetic_maps_m0_filtered)
    # Export Genetic Maps (unfiltered and filtered) with unified X axis
    generate_genetic_map_images_and_pdf(genetic_maps_m0_filtered, global_xmax)


    # Export Comparative Maps (unfiltered and filtered)
    generate_comparative_genetic_map_images_and_pdf(
        genetic_maps_comparative,
        genetic_maps_comparative_filtered,
        global_xmax
    )

    # Run SMT for all traits and export images and PDF
    """all_smt_results = run_smt_for_all_traits(vcf_data, traits, sample_names)
    generate_smt_chromosome_images_all_traits(all_smt_results)
    stitch_all_smt_chromosome_images(all_smt_results)"""
    # Generate evenly spaced markers from filtered genetic maps (3–5 cM apart)
    selected_markers = select_evenly_spaced_markers(genetic_maps_m0_filtered, vcf_data, min_cm_spacing=1.0)

    # Generate pairwise recombination maps for heatmaps
    """
    pairwise_maps = {}
    for chrom_id in genetic_maps_m0_filtered:
        pairwise = generate_pairwise_genetic_map(vcf_data, sample_names, chrom_id)
        if pairwise:
            pairwise_maps[chrom_id] = pairwise

    # Export Heatmaps with unified X/Y axis
    generate_heatmap_images_and_pdf(pairwise_maps)
    """

    # Run SMT for selected trait
    # df = run_smt_regression_and_export(vcf_data, traits, trait_to_test, sample_names)

    # Generate per-chromosome SMT significance images and combine to PDF
    # generate_smt_chromosome_images(df, trait_to_test, significance_level=pvalue_cutoff)
    # final_pdf_path = stitch_smt_chromosome_images_pdf(trait_to_test)

    # Run SMT for all traits and export combined QTL maps
    # all_smt_results = run_smt_for_all_traits(vcf_data, traits, sample_names)
    # generate_concatenated_qtl_pdf(all_smt_results)

    # Run SIM
    
    compute_genetic_map_by_recombination(selected_markers, sample_names)
    sim_results = run_sim_on_selected_markers(vcf_data, selected_markers, traits, sample_names)
    plot_all_lod_curves(sim_results)
    

if __name__ == "__main__":
    main()
