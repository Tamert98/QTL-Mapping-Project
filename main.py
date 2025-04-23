"""
main.py

This is the main script for running the QTL mapping pipeline. It orchestrates 
the full analysis by loading the input VCF file, applying the SMT and SIM 
methods, and displaying or saving the results.

Main Responsibilities:
- Load genotype data via vcf_parser.
- Load phenotype data from file or hardcoded input.
- Run SMT and SIM analyses.
- Display results or export them to files.

To use:
    python main.py --vcf your_file.vcf --phenotype phenotype.csv
"""

import os
from vcf_data_handler import parse_vcf_file, generate_genetic_map,parse_trait_file
from plot_utils import plot_genetic_map, plot_genetic_distance_circles,plot_full_recombination_dots
import matplotlib.pyplot as plt

def main():
    # Set the path to the VCF file (relative to this script)
    vcf_path = os.path.join("Data", "cataglyphis.final.DZ (1).vcf")

    # Parse the VCF file
    vcf_data, sample_names = parse_vcf_file(vcf_path)

    # Print number of chromosomes and first 5 chromosome IDs
    chrom_ids = list(vcf_data.keys())
    print(f"Number of chromosomes: {len(chrom_ids)}")
    print("First 5 chromosome IDs:", chrom_ids[:5])

    # Print number of markers and first 4 markers in chr01
    chr_id = "chr01"
    if chr_id in vcf_data:
        chr01_markers = vcf_data[chr_id]["markers"]
        print(f"\nNumber of markers in {chr_id}: {len(chr01_markers)}")
        print(f"First 4 markers in {chr_id} (showing 5 samples each):")
        for i, marker in enumerate(chr01_markers[:4]):
            print(f"\nMarker {i+1}: pos = {marker['pos']}, ref = {marker['ref']}, alt = {marker['alt']}, qual = {marker['qual']}")
            print("Samples (first 5):")
            for j, (sample, values) in enumerate(marker["samples"].items()):
                if j >= 5:
                    break
                print(f"  {sample}: {values}")
    else:
        print(f"{chr_id} not found in VCF data.")

    # Generate and print genetic map
    genetic_maps = generate_genetic_map(vcf_data, sample_names)
    if chr_id in genetic_maps:
        print(f"\nFirst 3 elements in genetic map for {chr_id}:")
        for entry in genetic_maps[chr_id]["genetic_map"][:3]:
            print(entry)
    else:
        print(f"{chr_id} not found in genetic maps.")

    trait_path = os.path.join("Data", "traits.txt")
    sample_names, traits = parse_trait_file(trait_path)

    print("Available traits:", list(traits.keys())[:3])
    print("Values for TC25Me3 (first 5 samples):")
    for name in sample_names[:5]:
        print(name, "->", traits["TC25Me3"].get(name))
    
    # Generate plots
    fig1 = plot_genetic_map(genetic_maps[chr_id]["genetic_map"], chr_id)
    fig2 = plot_full_recombination_dots(genetic_maps[chr_id]["genetic_map"], chr_id)
    plt.show()


if __name__ == "__main__":
    main()