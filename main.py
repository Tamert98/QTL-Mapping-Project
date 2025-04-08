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
from vcf_data_handler import parse_vcf_file

def main():
    # Set the path to the VCF file (relative to this script)
    vcf_path = os.path.join("Data", "cataglyphis.final.DZ (1).vcf")  # Replace with actual filename

    # Parse the VCF file
    vcf_data, sample_names = parse_vcf_file(vcf_path)

    # Show info for sanity check
    print(f"Parsed chromosomes: {list(vcf_data.keys())}")
    print(f"Sample names: {sample_names}")

    # Preview chr01 data
    target_chr = "chr01"
    if target_chr in vcf_data:
        print(f"\nData for {target_chr}:")
        print(f"Length: {vcf_data[target_chr]['length']}")
        print(f"Number of markers: {len(vcf_data[target_chr]['markers'])}")
        for marker in vcf_data[target_chr]["markers"][:5]:
            print(marker)
    else:
        print(f"{target_chr} not found in VCF.")

if __name__ == "__main__":
    main()