# vcf_data_handler.py
import math
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

def select_starting_marker(markers):
    """
    Placeholder to decide the best starting marker based on PL, DP, GQ.
    For now, returns the first marker.
    """
    return markers[0]

def generate_genetic_map(vcf_data, sample_names):
    genetic_maps = {}

    for chrom_id, chrom_data in vcf_data.items():
        markers = chrom_data["markers"]
        if len(markers) < 2:
            continue  # Skip chromosomes with fewer than 2 markers

        # Start with the first marker
        m1 = select_starting_marker(markers)
        m1_samples = m1["samples"]

        genetic_map = [{"pos": m1["pos"], "cM": 0.0}]

        for m2 in markers[1:]:
            m2_samples = m2["samples"]
            i00 = i01 = i10 = i11 = known = 0

            for name in sample_names:
                gt1 = m1_samples.get(name, {}).get("GT")
                gt2 = m2_samples.get(name, {}).get("GT")

                if gt1 in {"0/0", "1/1"} and gt2 in {"0/0", "1/1"}:
                    known += 1
                    if gt1 == "0/0" and gt2 == "0/0":
                        i00 += 1
                    elif gt1 == "0/0" and gt2 == "1/1":
                        i01 += 1
                    elif gt1 == "1/1" and gt2 == "0/0":
                        i10 += 1
                    elif gt1 == "1/1" and gt2 == "1/1":
                        i11 += 1

            if known == 0:
                continue  # Skip if no valid genotype pairs

            recomb_count = min(i01 + i10, i00 + i11)
            r = recomb_count / known

            try:
                distance = -50 * math.log(1 - 2 * r)
            except ValueError:
                distance = float('inf')  # In case r > 0.5 or invalid

            last_cM = genetic_map[-1]["cM"]
            genetic_map.append({
                "pos": m2["pos"],
                "cM": last_cM + distance
            })

            # Update m1 for the next comparison
            m1 = m2
            m1_samples = m2_samples

        # Store genetic map for the chromosome
        genetic_maps[chrom_id] = {
            "length": chrom_data["length"],
            "genetic_map": genetic_map
        }

    return genetic_maps