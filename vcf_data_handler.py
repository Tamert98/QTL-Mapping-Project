# vcf_data_handler.py

def parse_vcf_file(vcf_file_path):
    vcf_data = {}
    sample_names = []
    header_index_map = {}

    with open(vcf_file_path, "r") as vcf_file:
        for line in vcf_file:
            line = line.strip()

            # Step 1: Parse contig lines to get chromosome ID and length
            if line.startswith("##contig=<ID="):
                content = line[line.find("<")+1 : line.find(">")]
                attrs = dict(item.split("=") for item in content.split(","))
                chrom_id = attrs["ID"]
                chrom_length = int(attrs["length"])
                vcf_data[chrom_id] = {
                    "length": chrom_length,
                    "markers": []
                }

            # Step 2: Header line with sample names and columns
            elif line.startswith("#CHROM"):
                header_parts = line.split()
                header_index_map = {name: idx for idx, name in enumerate(header_parts)}
                try:
                    format_index = header_index_map["FORMAT"]
                    sample_names = header_parts[format_index + 1:]
                except KeyError:
                    print("ERROR: 'FORMAT' column not found in header.")
                    return {}

            # Step 3: Parse data lines
            elif not line.startswith("#") and sample_names:
                parts = line.split()
                chrom_id = parts[header_index_map["#CHROM"]]
                pos = int(parts[header_index_map["POS"]])
                ref = parts[header_index_map["REF"]]
                alt = parts[header_index_map["ALT"]]
                qual = float(parts[header_index_map["QUAL"]])
                format_keys = parts[header_index_map["FORMAT"]].split(":")
                sample_values = parts[header_index_map["FORMAT"] + 1:]

                sample_data = {}
                for name, val in zip(sample_names, sample_values):
                    values = val.split(":")
                    sample_data[name] = dict(zip(format_keys, values))

                if chrom_id in vcf_data:
                    vcf_data[chrom_id]["markers"].append({
                        "pos": pos,
                        "ref": ref,
                        "alt": alt,
                        "qual": qual,
                        "samples": sample_data
                    })

    return vcf_data, sample_names
