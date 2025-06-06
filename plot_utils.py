import matplotlib.pyplot as plt
import numpy as np
import random
def plot_genetic_map(genetic_map, chr_id="unknown"):
    """
    Plots the genetic map for a single chromosome.

    Parameters:
    - genetic_map: list of dicts, each with 'pos' and 'cM' keys
    - chr_id: optional chromosome label for the title

    Returns:
    - matplotlib Figure object
    """
    if not genetic_map:
        print(f"No data to plot for {chr_id}.")
        return None

    positions = [entry["pos"] for entry in genetic_map]
    distances = [entry["cM"] for entry in genetic_map]

    fig = plt.figure(figsize=(10, 6))
    plt.plot(positions, distances, marker='o', linestyle='-', linewidth=2)
    plt.title(f"Genetic Map for {chr_id}")
    plt.xlabel("Marker Position (bp)")
    plt.ylabel("Genetic Distance (cM)")
    plt.grid(True)
    plt.tight_layout()
    return fig


def plot_genetic_distance_circles(genetic_map, chr_id="unknown"):
    """
    Samples 100 markers from the genetic map and creates a grayscale heatmap
    using circle intensity to represent relative cM distance between marker pairs.

    X and Y axes are marker positions in base pairs (bp).
    """
    n = len(genetic_map)
    if n == 0:
        print("No markers to sample.")
        return

    num_samples = 80
    step = max(1, n // num_samples)

    sampled = []

    for i in range(num_samples):
        # Choose a random offset within the step-sized chunk
        offset = random.randint(0, step - 1)
        idx = i * step + offset

        if idx >= n:
            break  # Avoid out-of-bounds at the end

        marker = genetic_map[idx]
        sampled.append([marker["pos"], marker["cM"]])

    sampled = np.array(sampled)
    positions = sampled[:, 0]
    distances = sampled[:, 1]
    dmax = distances[-1] if len(distances) > 0 else 1e-6  # Last sample's cM

    # Step 2: Create figure and axes
    fig, ax = plt.subplots(figsize=(10, 10))

    for i in range(len(sampled)):
        for j in range(len(sampled)):
            dij = abs(distances[i] - distances[j])
            p = dij / (dmax + 1e-6)  # Prevent divide-by-zero
            p = min(max(p, 0.05), 0.95)  # Clamp for visible gray range
            color = (p, p, p)

            # Normalize positions to integer grid for spacing
            circle = plt.Circle((positions[j], positions[i]), radius=(positions[-1] - positions[0]) * 0.005, color=color)
            ax.add_patch(circle)

    # Step 3: Set axis limits and labels
    ax.set_xlim(positions.min() - 1e5, positions.max() + 1e5)
    ax.set_ylim(positions.min() - 1e5, positions.max() + 1e5)
    ax.set_xticks(np.linspace(positions.min(), positions.max(), 10))
    ax.set_yticks(np.linspace(positions.min(), positions.max(), 10))
    ax.set_xticklabels([f"{x:.1e}" for x in np.linspace(positions.min(), positions.max(), 10)], rotation=90)
    ax.set_yticklabels([f"{x:.1e}" for x in np.linspace(positions.min(), positions.max(), 10)])
    ax.set_title(f"Pairwise Distance HeatMap for {chr_id} (100 sampled markers)")
    ax.set_xlabel("Marker Position (bp)")
    ax.set_ylabel("Marker Position (bp)")
    ax.set_aspect('equal')
    fig.tight_layout()

    return fig


def plot_full_recombination_heatmap(pairwise_map, chr_id="unknown"):
    """
    Plots a dot for each pair of marker positions on a chromosome using actual pairwise recombination rates.
    Dot color = recombination rate (r), normalized by a fixed maximum of 200.
    """
    if not pairwise_map:
        print("No pairwise recombination data to plot.")
        return None

    # Extract valid entries
    valid_entries = [
        entry for entry in pairwise_map
        if isinstance(entry.get("pos_i"), (int, float)) and
           isinstance(entry.get("pos_j"), (int, float)) and
           isinstance(entry.get("r"), (int, float))
    ]
    if len(valid_entries) < 1:
        print("No valid pairwise recombination entries found.")
        return None

    # Extract positions and values
    x = [entry["pos_i"] for entry in valid_entries]
    y = [entry["pos_j"] for entry in valid_entries]
    r_values = [entry["r"] for entry in valid_entries]

    # Normalize using fixed constant 200
    r_norm = np.clip(np.array(r_values) / 200.0, 0, 1)

    fig, ax = plt.subplots(figsize=(10, 10))
    scatter = ax.scatter(x, y, c=r_norm, cmap="gray", s=4, marker='s')
    ax.set_title(f"Pairwise Recombination Dot Heatmap - {chr_id}")
    ax.set_xlabel("Marker Position (bp)")
    ax.set_ylabel("Marker Position (bp)")
    ax.set_aspect('equal')
    fig.tight_layout()

    return fig

def plot_f_statistic(df, chrom_id):
    """
    Plots F-statistic values across a specified chromosome.

    Args:
        df (pd.DataFrame): DataFrame with SMT results.
        chrom_id (str): Chromosome to plot.

    Returns:
        matplotlib Figure object
    """
    chrom_data = df[df['chrom'] == chrom_id]
    if chrom_data.empty:
        print(f"No data found for chromosome {chrom_id}")
        return None

    fig = plt.figure()
    plt.plot(chrom_data['pos'], chrom_data['F'], marker='o')
    plt.title(f"F-statistic across chromosome {chrom_id}")
    plt.xlabel("Marker Position")
    plt.ylabel("F-statistic")
    plt.grid(True)
    plt.tight_layout()
    return fig


def plot_p_values(df, chrom_id):
    """
    Plots p-values (log scale) across a specified chromosome.

    Args:
        df (pd.DataFrame): DataFrame with SMT results.
        chrom_id (str): Chromosome to plot.

    Returns:
        matplotlib Figure object
    """
    chrom_data = df[df['chrom'] == chrom_id]
    if chrom_data.empty:
        print(f"No data found for chromosome {chrom_id}")
        return None

    fig = plt.figure()
    plt.plot(chrom_data['pos'], chrom_data['p'], marker='o', color='red')
    plt.axhline(y=0.05, color='gray', linestyle='--', label='p = 0.05')
    plt.yscale("log")
    plt.title(f"P-values across chromosome {chrom_id}")
    plt.xlabel("Marker Position")
    plt.ylabel("P-value (log scale)")
    plt.legend()
    plt.grid(True, which='both')
    plt.tight_layout()
    return fig


def plot_log_p_values(df, chrom_id, significance_level=0.05):
    """
    Plots QTL significance (-log10(p-value)) across a specified chromosome.

    Args:
        df (pd.DataFrame): DataFrame with SMT results.
        chrom_id (str): Chromosome to plot.
        significance_level (float): Threshold for significance (default = 0.05).

    Returns:
        matplotlib Figure object
    """
    chrom_data = df[df['chrom'] == chrom_id]
    if chrom_data.empty:
        print(f"No data found for chromosome {chrom_id}")
        return None

    log_p = -np.log10(chrom_data['p'])
    threshold = -np.log10(significance_level)

    fig = plt.figure()
    plt.plot(chrom_data['pos'], log_p, marker='o', color='darkblue')
    plt.axhline(y=threshold, color='red', linestyle='--', label=f'p = {significance_level}')
    plt.title(f"QTL Significance across Chromosome {chrom_id}")
    plt.xlabel("Marker Position (bp)")
    plt.ylabel("QTL Significance (-log10 p-value)")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    return fig

def plot_combined_qtl_significance(df, significance_threshold=0.05):
    df["-log10(p)"] = -np.log10(df["p"])
    df["significant"] = df["p"] < significance_threshold
    significant_df = df[df["significant"]].sort_values(by="pos")

    if significant_df.empty:
        print("No significant markers found.")
        return None

    fig, ax = plt.subplots(figsize=(12, 6))
    ax.plot(significant_df["pos"], significant_df["-log10(p)"], marker="o", linestyle='-', color="blue")
    ax.axhline(-np.log10(significance_threshold), color="red", linestyle="--", label=f"p = {significance_threshold}")
    ax.set_title("QTL Significance Map (based on p-value)")
    ax.set_xlabel("Marker Position (bp)")
    ax.set_ylabel("QTL Significance (higher = stronger)")
    ax.legend()
    ax.grid(True)
    fig.tight_layout()
    return fig