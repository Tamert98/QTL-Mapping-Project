"""
vcf_data_handler.py

This module is responsible for reading, parsing, and organizing genetic data 
from a VCF (Variant Call Format) file. It provides functionality to convert 
the data into a structured format suitable for QTL analysis.

Main Responsibilities:
- Load VCF files.
- Extract genotype and marker information.
- Organize data for use in QTL mapping methods (SIM and SMT).

"""

import os

# Set the relative path to the VCF file
vcf_file_path = os.path.join("Data", "cataglyphis.final.DZ (1).vcf")  # Replace with your actual filename

# Data structure to store parsed info
vcf_data = {}
sample_names = []

with open(vcf_file_path, "r") as vcf_file:
    for line in vcf_file:
        line = line.strip()

        # Step 1: Parse contig lines to get chromosome ID and length
        if line.startswith("##contig=<ID="):
            # Extract content between the angled brackets
            content = line[line.find("<")+1 : line.find(">")]
            attrs = dict(
                item.split("=") for item in content.split(",")
            )
            chrom_id = attrs["ID"]
            chrom_length = int(attrs["length"])
            vcf_data[chrom_id] = {
                "length": chrom_length,
                "markers": []
            }

        # Step 2: Get the header line with sample names
        elif line.startswith("#CHROM"):
            header_parts = line.split()
            header_index_map = {name: idx for idx, name in enumerate(header_parts)}
            try:
                format_index = header_index_map["FORMAT"]
                sample_names = header_parts[format_index + 1:]
            except KeyError:
                print("ERROR: 'FORMAT' column not found in header.")
                sample_names = []
    

        # Step 3: Process genotype data lines
        elif not line.startswith("#") and sample_names:
            parts = line.split()
            chrom_id = parts[header_index_map["#CHROM"]]
            pos = int(parts[header_index_map["POS"]])
            ref = parts[header_index_map["REF"]]
            alt = parts[header_index_map["ALT"]]
            qual = float(parts[header_index_map["QUAL"]])
            format_keys = parts[header_index_map["FORMAT"]].split(":")

            sample_values = parts[header_index_map["FORMAT"] + 1:]


            # Build per-sample genotype dictionary
            sample_data = {}
            for name, val in zip(sample_names, sample_values):
                values = val.split(":")
                sample_data[name] = dict(zip(format_keys, values))

            # Add marker to chromosome's marker list
            if chrom_id in vcf_data:
                vcf_data[chrom_id]["markers"].append({
                    "pos": pos,
                    "ref": ref,
                    "alt": alt,
                    "qual": qual,
                    "samples": sample_data
                })

target_chr = "chr01"
if target_chr in vcf_data:
    print(f"\nData for {target_chr}:")
    print(f"Length: {vcf_data[target_chr]['length']}")
    print(f"Number of markers: {len(vcf_data[target_chr]['markers'])}")
    for marker in vcf_data[target_chr]["markers"][:4]:  # Just print first 5 for sanity
        print(marker)
else:
    print(f"{target_chr} not found in VCF data.")