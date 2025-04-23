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
import matplotlib.pyplot as plt
import numpy as np

def plot_full_recombination_dots(genetic_map, chr_id="unknown"):
    """
    Plots a dot for each pair of marker positions on a chromosome.
    Dot color = absolute recombination rate difference (|rᵢ - rⱼ|), normalized.
    Black = similar (close markers), White = high recombination distance.
    """
    if not genetic_map:
        print("No markers to plot.")
        return None

    # Filter valid markers with numeric pos and r
    valid_markers = [
        m for m in genetic_map
        if isinstance(m.get("pos"), (int, float)) and isinstance(m.get("r"), (int, float))
    ]
    if len(valid_markers) < 2:
        print("Not enough valid markers.")
        return None

    positions = np.array([m["pos"] for m in valid_markers])
    r_values = np.array([m["r"] for m in valid_markers])

    if np.all(r_values == r_values[0]):
        print("All recombination rates are equal — skipping plot.")
        return None

    # Pairwise recombination difference matrix
    pos_i, pos_j = np.meshgrid(positions, positions)
    r_i, r_j = np.meshgrid(r_values, r_values)

    r_diff = np.abs(r_i - r_j)
    r_norm = r_diff / np.max(r_diff)  # Normalized to [0, 1]

    # Flatten for plotting
    x = pos_j.flatten()
    y = pos_i.flatten()
    c = r_norm.flatten()

    fig, ax = plt.subplots(figsize=(10, 10))
    ax.scatter(x, y, c=c, cmap="gray", s=2, marker='s')  # grayscale
    ax.set_title(f"Pairwise Recombination Map - {chr_id}")
    ax.set_xlabel("Marker Position (bp)")
    ax.set_ylabel("Marker Position (bp)")
    ax.set_aspect('equal')
    fig.tight_layout()

    return fig
