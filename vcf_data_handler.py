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
                    raw_sample_names = header_parts[format_index + 1:]
                    sample_names = [name.split(".")[0] for name in raw_sample_names]

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



# Generates a genetic map based on the VCF data and sample names.
# Does not check the validity of the marker data (e.g., allele balance or missingness).
def generate_genetic_map(vcf_data, sample_names):
    genetic_maps = {}

    for chrom_id, chrom_data in vcf_data.items():
        markers = chrom_data["markers"]
        if len(markers) < 2:
            continue

        m0 = markers[0]
        m0_samples = m0["samples"]

        m0_i00 = sum(1 for name in sample_names if m0_samples.get(name, {}).get("GT") == "0/0")
        m0_i11 = sum(1 for name in sample_names if m0_samples.get(name, {}).get("GT") == "1/1")

        genetic_map = [{
            "pos": m0["pos"],
            "cM": 0.0,
            "r": 0.0,
            "i00": m0_i00,
            "i11": m0_i11
        }]

        for m in markers[1:]:
            m_samples = m["samples"]
            i00 = i01 = i10 = i11 = known = 0

            for name in sample_names:
                gt0 = m0_samples.get(name, {}).get("GT")
                gt = m_samples.get(name, {}).get("GT")

                if gt0 in {"0/0", "1/1"} and gt in {"0/0", "1/1"}:
                    known += 1
                    if gt0 == "0/0" and gt == "0/0":
                        i00 += 1
                    elif gt0 == "0/0" and gt == "1/1":
                        i01 += 1
                    elif gt0 == "1/1" and gt == "0/0":
                        i10 += 1
                    elif gt0 == "1/1" and gt == "1/1":
                        i11 += 1

            mi_i00 = sum(1 for name in sample_names if m_samples.get(name, {}).get("GT") == "0/0")
            mi_i11 = sum(1 for name in sample_names if m_samples.get(name, {}).get("GT") == "1/1")

            if known == 0:
                genetic_map.append({
                    "pos": m["pos"],
                    "cM": genetic_map[-1]["cM"],
                    "r": None,
                    "i00": mi_i00,
                    "i11": mi_i11
                })
                continue

            recomb_count = min(i01 + i10, i00 + i11)
            r = recomb_count / known

            try:
                distance = -50 * math.log(1 - 2 * r)
            except ValueError:
                distance = float('inf')

            genetic_map.append({
                "pos": m["pos"],
                "cM": distance,
                "r": r,
                "i00": mi_i00,
                "i11": mi_i11
            })

        genetic_maps[chrom_id] = {
            "length": chrom_data["length"],
            "genetic_map": genetic_map
        }

    return genetic_maps

# Generates a genetic map using only markers that pass quality filtering:
# includes only markers where both 0/0 and 1/1 genotype counts are > 25% of total samples.
# Skips markers that do not meet this threshold, and selects the first valid one as the starting point.
def generate_genetic_map_filtered1(vcf_data, sample_names):
    genetic_maps = {}
    n = len(sample_names)
    threshold = 0.25 * n

    for chrom_id, chrom_data in vcf_data.items():
        markers = chrom_data["markers"]
        if len(markers) < 2:
            continue

        # Find first marker that satisfies the filtering condition
        m0 = None
        for marker in markers:
            samples = marker["samples"]
            i00 = sum(1 for name in sample_names if samples.get(name, {}).get("GT") == "0/0")
            i11 = sum(1 for name in sample_names if samples.get(name, {}).get("GT") == "1/1")
            if i00 > threshold and i11 > threshold:
                m0 = marker
                break

        if m0 is None:
            continue  # No valid m0 found

        m0_samples = m0["samples"]
        m0_i00 = sum(1 for name in sample_names if m0_samples.get(name, {}).get("GT") == "0/0")
        m0_i11 = sum(1 for name in sample_names if m0_samples.get(name, {}).get("GT") == "1/1")

        genetic_map = [{
            "pos": m0["pos"],
            "cM": 0.0,
            "r": 0.0,
            "i00": m0_i00,
            "i11": m0_i11
        }]

        m0_index = markers.index(m0)

        for m in markers[m0_index + 1:]:
            m_samples = m["samples"]
            mi_i00 = sum(1 for name in sample_names if m_samples.get(name, {}).get("GT") == "0/0")
            mi_i11 = sum(1 for name in sample_names if m_samples.get(name, {}).get("GT") == "1/1")

            # Skip markers that do not pass the threshold
            if mi_i00 <= threshold or mi_i11 <= threshold:
                continue

            i00 = i01 = i10 = i11 = known = 0
            for name in sample_names:
                gt0 = m0_samples.get(name, {}).get("GT")
                gt = m_samples.get(name, {}).get("GT")

                if gt0 in {"0/0", "1/1"} and gt in {"0/0", "1/1"}:
                    known += 1
                    if gt0 == "0/0" and gt == "0/0":
                        i00 += 1
                    elif gt0 == "0/0" and gt == "1/1":
                        i01 += 1
                    elif gt0 == "1/1" and gt == "0/0":
                        i10 += 1
                    elif gt0 == "1/1" and gt == "1/1":
                        i11 += 1

            if known == 0:
                genetic_map.append({
                    "pos": m["pos"],
                    "cM": genetic_map[-1]["cM"],
                    "r": None,
                    "i00": mi_i00,
                    "i11": mi_i11
                })
                continue

            recomb_count = min(i01 + i10, i00 + i11)
            r = recomb_count / known

            try:
                distance = -50 * math.log(1 - 2 * r)
            except ValueError:
                distance = float('inf')

            genetic_map.append({
                "pos": m["pos"],
                "cM": distance,
                "r": r,
                "i00": mi_i00,
                "i11": mi_i11
            })

        if len(genetic_map) > 1:
            genetic_maps[chrom_id] = {
                "length": chrom_data["length"],
                "genetic_map": genetic_map
            }

    return genetic_maps


def is_good_genotype(gt_info, min_dp=10, min_gq=20, min_delta_pl=20):
    try:
        dp = int(gt_info.get("DP", 0))
        gq = int(gt_info.get("GQ", 0))
        pl_str = gt_info.get("PL", "")
        pl = list(map(int, pl_str.split(","))) if pl_str else []

        if dp < min_dp or gq < min_gq or len(pl) < 2:
            return False

        sorted_pl = sorted(pl)
        return (sorted_pl[1] - sorted_pl[0]) >= min_delta_pl
    except Exception:
        return False

def generate_genetic_map_filtered(vcf_data, sample_names):
    genetic_maps = {}
    n = len(sample_names)
    threshold = 0.25 * n

    for chrom_id, chrom_data in vcf_data.items():
        markers = chrom_data["markers"]
        if len(markers) < 2:
            continue

        # Find the first valid m0
        m0 = None
        for marker in markers:
            samples = marker["samples"]
            i00 = i11 = 0
            for name in sample_names:
                gt_info = samples.get(name, {})
                if not is_good_genotype(gt_info):
                    continue
                gt = gt_info.get("GT")
                if gt == "0/0":
                    i00 += 1
                elif gt == "1/1":
                    i11 += 1
            if i00 > threshold and i11 > threshold:
                m0 = marker
                break

        if m0 is None:
            continue

        m0_samples = m0["samples"]
        m0_i00 = sum(1 for name in sample_names
                     if is_good_genotype(m0_samples.get(name, {})) and m0_samples.get(name, {}).get("GT") == "0/0")
        m0_i11 = sum(1 for name in sample_names
                     if is_good_genotype(m0_samples.get(name, {})) and m0_samples.get(name, {}).get("GT") == "1/1")

        genetic_map = [{
            "pos": m0["pos"],
            "cM": 0.0,
            "r": 0.0,
            "i00": m0_i00,
            "i11": m0_i11
        }]

        m0_index = markers.index(m0)

        for m in markers[m0_index + 1:]:
            m_samples = m["samples"]

            # Filter counts based on good genotypes
            mi_i00 = sum(1 for name in sample_names
                         if is_good_genotype(m_samples.get(name, {})) and m_samples.get(name, {}).get("GT") == "0/0")
            mi_i11 = sum(1 for name in sample_names
                         if is_good_genotype(m_samples.get(name, {})) and m_samples.get(name, {}).get("GT") == "1/1")

            if mi_i00 <= threshold or mi_i11 <= threshold:
                continue

            i00 = i01 = i10 = i11 = known = 0
            for name in sample_names:
                gt0_info = m0_samples.get(name, {})
                gt_info = m_samples.get(name, {})
                if not (is_good_genotype(gt0_info) and is_good_genotype(gt_info)):
                    continue

                gt0 = gt0_info.get("GT")
                gt = gt_info.get("GT")

                if gt0 in {"0/0", "1/1"} and gt in {"0/0", "1/1"}:
                    known += 1
                    if gt0 == "0/0" and gt == "0/0":
                        i00 += 1
                    elif gt0 == "0/0" and gt == "1/1":
                        i01 += 1
                    elif gt0 == "1/1" and gt == "0/0":
                        i10 += 1
                    elif gt0 == "1/1" and gt == "1/1":
                        i11 += 1

            if known == 0:
                genetic_map.append({
                    "pos": m["pos"],
                    "cM": genetic_map[-1]["cM"],
                    "r": None,
                    "i00": mi_i00,
                    "i11": mi_i11
                })
                continue

            recomb_count = min(i01 + i10, i00 + i11)
            r = recomb_count / known

            try:
                distance = -50 * math.log(1 - 2 * r)
            except ValueError:
                distance = float('inf')

            genetic_map.append({
                "pos": m["pos"],
                "cM": distance,
                "r": r,
                "i00": mi_i00,
                "i11": mi_i11
            })

        if len(genetic_map) > 1:
            genetic_maps[chrom_id] = {
                "length": chrom_data["length"],
                "genetic_map": genetic_map
            }

    return genetic_maps





def is_good_genotype(gt_info, min_dp=10, min_gq=20, min_delta_pl=20):
    try:
        dp = int(gt_info.get("DP", 0))
        gq = int(gt_info.get("GQ", 0))
        pl_str = gt_info.get("PL", "")
        pl = list(map(int, pl_str.split(","))) if pl_str else []

        if dp < min_dp or gq < min_gq or len(pl) < 2:
            return False

        sorted_pl = sorted(pl)
        return (sorted_pl[1] - sorted_pl[0]) >= min_delta_pl
    except Exception:
        return False

def generate_pairwise_genetic_map(vcf_data, sample_names, chrom_id):
    n = len(sample_names)
    threshold = 0.25 * n

    if chrom_id not in vcf_data:
        print(f"Chromosome '{chrom_id}' not found in VCF data.")
        return None

    markers = vcf_data[chrom_id]["markers"]
    if len(markers) < 2:
        print(f"Not enough markers in chromosome '{chrom_id}'.")
        return None

    # Step 1: Filter valid markers by quality
    valid_indices = []
    for idx, marker in enumerate(markers):
        samples = marker["samples"]
        i00 = i11 = 0
        for name in sample_names:
            info = samples.get(name, {})
            if not is_good_genotype(info):
                continue
            gt = info.get("GT")
            if gt == "0/0":
                i00 += 1
            elif gt == "1/1":
                i11 += 1
        if i00 > threshold and i11 > threshold:
            valid_indices.append(idx)

    if len(valid_indices) < 2:
        print(f"Not enough valid markers in chromosome '{chrom_id}' after filtering.")
        return None

    # Step 2: Compute pairwise recombination symmetrically
    pairwise_map = []
    for i in range(len(valid_indices)):
        for j in range(len(valid_indices)):
            if i == j:
                continue  # skip same marker
            idx_i = valid_indices[i]
            idx_j = valid_indices[j]

            mi = markers[idx_i]
            mj = markers[idx_j]

            mi_samples = mi["samples"]
            mj_samples = mj["samples"]

            i00 = i01 = i10 = i11 = known = 0
            for name in sample_names:
                gt_i_info = mi_samples.get(name, {})
                gt_j_info = mj_samples.get(name, {})

                if not (is_good_genotype(gt_i_info) and is_good_genotype(gt_j_info)):
                    continue

                gt_i = gt_i_info.get("GT")
                gt_j = gt_j_info.get("GT")

                if gt_i in {"0/0", "1/1"} and gt_j in {"0/0", "1/1"}:
                    known += 1
                    if gt_i == "0/0" and gt_j == "0/0":
                        i00 += 1
                    elif gt_i == "0/0" and gt_j == "1/1":
                        i01 += 1
                    elif gt_i == "1/1" and gt_j == "0/0":
                        i10 += 1
                    elif gt_i == "1/1" and gt_j == "1/1":
                        i11 += 1

            if known == 0:
                continue

            recomb_count = min(i01 + i10, i00 + i11)
            r = recomb_count / known

            try:
                distance = -50 * math.log(1 - 2 * r)
            except ValueError:
                distance = float('inf')

            # Count good GTs in mj (for reference)
            mj_i00 = sum(1 for name in sample_names
                         if is_good_genotype(mj_samples.get(name, {})) and mj_samples.get(name, {}).get("GT") == "0/0")
            mj_i11 = sum(1 for name in sample_names
                         if is_good_genotype(mj_samples.get(name, {})) and mj_samples.get(name, {}).get("GT") == "1/1")

            pairwise_map.append({
                "pos_i": mi["pos"],
                "pos_j": mj["pos"],
                "cM": distance,
                "r": r,
                "i00": mj_i00,
                "i11": mj_i11
            })

    return {
        "length": vcf_data[chrom_id]["length"],
        "pairwise_map": pairwise_map
    } if pairwise_map else None




def parse_trait_file(trait_file_path):
    """
    Parses a trait file with multiple traits and returns:
    - sample_names: list of cleaned sample names (e.g., 'male-1')
    - traits: dictionary mapping each trait name to a dictionary of {sample_name: value}
              where value is a float or None if unknown ($)
    """
    traits = {}

    with open(trait_file_path, "r") as f:
        lines = f.readlines()

        if len(lines) < 2:
            raise ValueError("Trait file must contain at least two lines.")

        # Line 1: sample names (e.g., male-1.bowtie)
        raw_samples = lines[0].strip().split()
        sample_names = [s.split('.')[0] for s in raw_samples]

        for line in lines[1:]:
            parts = line.strip().split()
            trait_name = parts[0]
            raw_values = parts[1:]

            if len(raw_values) != len(sample_names):
                raise ValueError(f"Trait '{trait_name}' has a mismatched number of values.")

            value_dict = {}
            for name, val in zip(sample_names, raw_values):
                if val == "$":
                    value_dict[name] = None
                else:
                    try:
                        value_dict[name] = float(val)
                    except ValueError:
                        value_dict[name] = None

            traits[trait_name] = value_dict

    return sample_names, traits