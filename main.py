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
