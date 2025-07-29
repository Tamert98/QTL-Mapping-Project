"""
vcf_data_handler.py

This module handles:
- Parsing VCF files to extract marker and sample genotype data.
- Computing recombination rates and generating genetic maps:
    * Relative to the first marker (m0-based map)
    * Cumulative genetic distance across chromosomes (comparative map)
- Filtering genotype data based on quality scores (DP, GQ, PL).

Used in QTL mapping pipelines for generating reliable marker distances
and genotype structures used in downstream SMT/SIM analyses.
"""

import math

# =============================================================================
# 1. VCF Parsing
# =============================================================================

def parse_vcf_file(vcf_file_path):
    """
    Parses a VCF file to extract marker information and sample genotypes.

    Returns:
        vcf_data (dict): {
            chrom_id: {
                "length": chromosome length,
                "markers": [ 
                    {
                        "pos": position,
                        "ref": reference allele,
                        "alt": alternate allele,
                        "qual": quality score,
                        "samples": {sample_name: {GT, DP, GQ, PL, ...}}
                    }, ...
                ]
            }, ...
        }
        sample_names (list): List of sample IDs (cleaned of suffixes like .bowtie)
    """
    vcf_data = {}
    sample_names = []
    header_index_map = {}

    with open(vcf_file_path, "r") as vcf_file:
        for line in vcf_file:
            line = line.strip()

            # === Contig metadata
            if line.startswith("##contig=<ID="):
                content = line[line.find("<")+1 : line.find(">")]
                attrs = dict(item.split("=") for item in content.split(","))
                chrom_id = attrs["ID"]
                chrom_length = int(attrs["length"])
                vcf_data[chrom_id] = {"length": chrom_length, "markers": []}

            # === Header line
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

            # === Marker data lines
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

# =============================================================================
# 2. Genetic Map Generation (Unfiltered)
# =============================================================================

def generate_combined_genetic_maps(vcf_data, sample_names):
    """
    Computes genetic maps using all markers, without quality filtering.

    Returns:
        comparative_maps (dict): Cumulative distance from marker to marker.
        m0_distance_maps (dict): Distance of each marker relative to m0.
    """
    comparative_maps = {}
    m0_distance_maps = {}

    for chrom_id, chrom_data in vcf_data.items():
        markers = chrom_data["markers"]
        if len(markers) < 2:
            continue

        m0 = markers[0]
        m0_samples = m0["samples"]
        m0_i00 = sum(1 for name in sample_names if m0_samples.get(name, {}).get("GT") == "0/0")
        m0_i11 = sum(1 for name in sample_names if m0_samples.get(name, {}).get("GT") == "1/1")

        comparative_map = [{
            "pos": m0["pos"], "cM": 0.0, "r": 0.0,
            "i00": m0_i00, "i11": m0_i11
        }]
        m0_distance_map = [comparative_map[0].copy()]
        total_cM = 0.0

        for m in markers[1:]:
            m_samples = m["samples"]
            i00 = i01 = i10 = i11 = known = 0

            for name in sample_names:
                gt0 = m0_samples.get(name, {}).get("GT")
                gt = m_samples.get(name, {}).get("GT")
                if gt0 in {"0/0", "1/1"} and gt in {"0/0", "1/1"}:
                    known += 1
                    if gt0 == "0/0" and gt == "0/0": i00 += 1
                    elif gt0 == "0/0" and gt == "1/1": i01 += 1
                    elif gt0 == "1/1" and gt == "0/0": i10 += 1
                    elif gt0 == "1/1" and gt == "1/1": i11 += 1

            mi_i00 = sum(1 for name in sample_names if m_samples.get(name, {}).get("GT") == "0/0")
            mi_i11 = sum(1 for name in sample_names if m_samples.get(name, {}).get("GT") == "1/1")

            if known == 0:
                entry = {
                    "pos": m["pos"], "cM": total_cM, "r": None,
                    "i00": mi_i00, "i11": mi_i11
                }
                comparative_map.append(entry.copy())
                m0_distance_map.append(entry.copy())
                continue

            recomb_count = min(i01 + i10, i00 + i11)
            r = recomb_count / known
            try:
                distance = -50 * math.log(1 - 2 * r)
            except ValueError:
                distance = float("inf")

            total_cM += distance
            comparative_map.append({
                "pos": m["pos"], "cM": total_cM, "r": r,
                "i00": mi_i00, "i11": mi_i11
            })
            m0_distance_map.append({
                "pos": m["pos"], "cM": distance, "r": r,
                "i00": mi_i00, "i11": mi_i11
            })

        comparative_maps[chrom_id] = {
            "length": chrom_data["length"],
            "genetic_map": comparative_map
        }
        m0_distance_maps[chrom_id] = {
            "length": chrom_data["length"],
            "genetic_map": m0_distance_map
        }

    return comparative_maps, m0_distance_maps

# =============================================================================
# 3. Genotype Quality Filtering
# =============================================================================

def is_good_genotype(gt_info, min_dp=30, min_gq=70, min_delta_pl=90):
    """
    Checks if a genotype is of high enough quality for use.

    Parameters:
        gt_info (dict): Genotype fields {GT, DP, GQ, PL, ...}
        min_dp (int): Minimum depth
        min_gq (int): Minimum genotype quality
        min_delta_pl (int): Minimum PL separation between best and second-best call

    Returns:
        bool: True if genotype passes all filters
    """
    try:
        dp = int(gt_info.get("DP", 0))
        gq = int(gt_info.get("GQ", 0))
        pl_str = gt_info.get("PL", "")
        pl = list(map(int, pl_str.split(","))) if pl_str else []

        if dp < min_dp or gq < min_gq:
            return False
        if len(pl) >= 2:
            sorted_pl = sorted(pl)
            return (sorted_pl[1] - sorted_pl[0]) >= min_delta_pl
        return True

    except Exception:
        return False

# =============================================================================
# 4. Genetic Map Generation (Filtered)
# =============================================================================

def generate_combined_genetic_maps_filtered(vcf_data, sample_names):
    """
    Computes filtered genetic maps (comparative and m0-based),
    excluding low-quality genotype data.

    Returns:
        genetic_maps_m0 (dict)
        genetic_maps_comparative (dict)
    """
    genetic_maps_m0 = {}
    genetic_maps_comparative = {}
    threshold = 0.25 * len(sample_names)

    for chrom_id, chrom_data in vcf_data.items():
        markers = chrom_data["markers"]
        if len(markers) < 2:
            continue

        # === Find valid m0 marker
        m0 = None
        for marker in markers:
            samples = marker["samples"]
            i00 = i11 = 0
            for name in sample_names:
                gt_info = samples.get(name, {})
                if not is_good_genotype(gt_info): continue
                gt = gt_info.get("GT")
                if gt == "0/0": i00 += 1
                elif gt == "1/1": i11 += 1
            if i00 > threshold and i11 > threshold:
                m0 = marker
                break

        if m0 is None:
            continue

        m0_samples = m0["samples"]
        m0_i00 = sum(1 for name in sample_names
                     if is_good_genotype(m0_samples.get(name, {})) and m0_samples[name].get("GT") == "0/0")
        m0_i11 = sum(1 for name in sample_names
                     if is_good_genotype(m0_samples.get(name, {})) and m0_samples[name].get("GT") == "1/1")

        map_m0 = []
        map_comparative = []
        cumulative_cM = 0.0

        map_m0.append({
            "pos": m0["pos"], "cM": 0.0, "r": 0.0,
            "i00": m0_i00, "i11": m0_i11
        })
        map_comparative.append(map_m0[0].copy())

        m0_index = markers.index(m0)
        prev_marker = m0

        for m in markers[m0_index + 1:]:
            m_samples = m["samples"]

            mi_i00 = sum(1 for name in sample_names
                         if is_good_genotype(m_samples.get(name, {})) and m_samples[name].get("GT") == "0/0")
            mi_i11 = sum(1 for name in sample_names
                         if is_good_genotype(m_samples.get(name, {})) and m_samples[name].get("GT") == "1/1")

            if mi_i00 <= threshold or mi_i11 <= threshold:
                continue

            i00 = i01 = i10 = i11 = known = 0
            for name in sample_names:
                gt0_info = prev_marker["samples"].get(name, {})
                gt_info = m_samples.get(name, {})
                if not (is_good_genotype(gt0_info) and is_good_genotype(gt_info)):
                    continue
                gt0 = gt0_info.get("GT")
                gt = gt_info.get("GT")
                if gt0 in {"0/0", "1/1"} and gt in {"0/0", "1/1"}:
                    known += 1
                    if gt0 == "0/0" and gt == "0/0": i00 += 1
                    elif gt0 == "0/0" and gt == "1/1": i01 += 1
                    elif gt0 == "1/1" and gt == "0/0": i10 += 1
                    elif gt0 == "1/1" and gt == "1/1": i11 += 1

            if known == 0:
                map_m0.append({
                    "pos": m["pos"], "cM": cumulative_cM, "r": None,
                    "i00": mi_i00, "i11": mi_i11
                })
                map_comparative.append({
                    "pos": m["pos"], "cM": 0.0, "r": None,
                    "i00": mi_i00, "i11": mi_i11
                })
                continue

            recomb_count = min(i01 + i10, i00 + i11)
            r = recomb_count / known
            try:
                distance = -50 * math.log(1 - 2 * r)
            except ValueError:
                distance = float('inf')

            cumulative_cM += distance
            map_m0.append({
                "pos": m["pos"], "cM": cumulative_cM, "r": r,
                "i00": mi_i00, "i11": mi_i11
            })
            map_comparative.append({
                "pos": m["pos"], "cM": distance, "r": r,
                "i00": mi_i00, "i11": mi_i11
            })

            prev_marker = m

        if len(map_m0) > 1:
            genetic_maps_m0[chrom_id] = {
                "length": chrom_data["length"],
                "genetic_map": map_m0
            }
        if len(map_comparative) > 1:
            genetic_maps_comparative[chrom_id] = {
                "length": chrom_data["length"],
                "genetic_map": map_comparative
            }

    return genetic_maps_m0, genetic_maps_comparative


# =============================================================================
# 5. Pairwise Genetic Map Generation (Filtered)
# =============================================================================

def generate_pairwise_genetic_map(vcf_data, sample_names, chrom_id):
    """
    Computes pairwise recombination rates and distances between all valid marker pairs
    on the given chromosome, based on high-quality genotype data only.

    Parameters:
        vcf_data (dict): Parsed VCF data by chromosome and markers.
        sample_names (list): List of all sample names.
        chrom_id (str): Chromosome ID to analyze.

    Returns:
        dict: {
            "length": chromosome length,
            "pairwise_map": [
                {
                    "pos_i": position of marker i,
                    "pos_j": position of marker j,
                    "cM": genetic distance in centimorgans,
                    "r": recombination rate,
                    "i00": count of 0/0 genotypes in marker j,
                    "i11": count of 1/1 genotypes in marker j
                }, ...
            ]
        } or None if insufficient valid markers.
    """
    n = len(sample_names)
    threshold = 0.25 * n

    if chrom_id not in vcf_data:
        print(f"Chromosome '{chrom_id}' not found in VCF data.")
        return None

    markers = vcf_data[chrom_id]["markers"]
    if len(markers) < 2:
        print(f"Not enough markers in chromosome '{chrom_id}'.")
        return None

    # === Step 1: Filter markers with sufficient high-quality 0/0 and 1/1 genotypes
    valid_indices = []
    for idx, marker in enumerate(markers):
        samples = marker["samples"]
        i00 = i11 = 0
        for name in sample_names:
            info = samples.get(name, {})
            if not is_good_genotype(info):
                continue
            gt = info.get("GT")
            if gt == "0/0": i00 += 1
            elif gt == "1/1": i11 += 1
        if i00 > threshold and i11 > threshold:
            valid_indices.append(idx)

    if len(valid_indices) < 2:
        print(f"Not enough valid markers in chromosome '{chrom_id}' after filtering.")
        return None

    # === Step 2: Compute recombination between valid marker pairs
    pairwise_map = []
    for i in range(len(valid_indices)):
        for j in range(len(valid_indices)):
            if i == j:
                continue

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
                    if gt_i == "0/0" and gt_j == "0/0": i00 += 1
                    elif gt_i == "0/0" and gt_j == "1/1": i01 += 1
                    elif gt_i == "1/1" and gt_j == "0/0": i10 += 1
                    elif gt_i == "1/1" and gt_j == "1/1": i11 += 1

            if known == 0:
                continue

            recomb_count = min(i01 + i10, i00 + i11)
            r = recomb_count / known
            try:
                distance = -50 * math.log(1 - 2 * r)
            except ValueError:
                distance = float('inf')

            mj_i00 = sum(1 for name in sample_names
                         if is_good_genotype(mj_samples.get(name, {})) and mj_samples[name].get("GT") == "0/0")
            mj_i11 = sum(1 for name in sample_names
                         if is_good_genotype(mj_samples.get(name, {})) and mj_samples[name].get("GT") == "1/1")

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

# =============================================================================
# 6. Trait File Parsing
# =============================================================================

def parse_trait_file(trait_file_path):
    """
    Parses a trait file with multiple traits into a structured dictionary.

    Returns:
        sample_names (list): List of cleaned sample names (prefix before dot).
        traits (dict): {
            trait_name: { sample_name: value (float or None) }
        }
    """
    traits = {}

    with open(trait_file_path, "r") as f:
        lines = f.readlines()

        if len(lines) < 2:
            raise ValueError("Trait file must contain at least two lines.")

        # === First line: Sample names
        raw_samples = lines[0].strip().split()
        sample_names = [s.split('.')[0] for s in raw_samples]

        # === Following lines: Trait values
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

# =============================================================================
# 7. Marker Selection Utility
# =============================================================================

def select_evenly_spaced_markers(genetic_maps_filtered, vcf_data, min_cm_spacing=3.0):
    """
    Selects markers from each chromosome spaced by at least `min_cm_spacing` cM.

    Parameters:
        genetic_maps_filtered (dict): Filtered genetic maps with cM values.
        vcf_data (dict): VCF marker dictionary used to extract sample data.
        min_cm_spacing (float): Minimum spacing between selected markers (in cM).

    Returns:
        dict: {
            chromosome_id: [
                { "pos": int, "cM": float, "samples": dict }, ...
            ]
        }
    """
    spaced_markers = {}

    for chrom, chrom_data in genetic_maps_filtered.items():
        if not chrom.lower().startswith("chr"):
            continue

        genetic_map = chrom_data.get("genetic_map", [])
        vcf_markers = {m["pos"]: m for m in vcf_data.get(chrom, {}).get("markers", [])}

        selected = []
        last_cm = -float("inf")

        for marker in genetic_map:
            cm = marker.get("cM")
            pos = marker.get("pos")
            if cm is None or pos is None or not isinstance(cm, (float, int)):
                continue

            if cm - last_cm >= min_cm_spacing:
                vcf_marker = vcf_markers.get(pos)
                if vcf_marker and "samples" in vcf_marker:
                    selected.append({
                        "pos": pos,
                        "cM": cm,
                        "samples": vcf_marker["samples"]
                    })
                    last_cm = cm

        if selected:
            spaced_markers[chrom] = selected

    return spaced_markers
