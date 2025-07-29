"""
plot_utils.py

This module contains all plotting utilities for QTL Mapping.

Responsibilities:
- Genetic map plotting (filtered and unfiltered)
- Comparative map plotting (filtered vs unfiltered)
- Pairwise recombination heatmap plotting
- SMT significance plotting (per chromosome and genome-wide)
- SIM LOD curve plotting
- Exporting plots as images and PDFs
"""

# === Imports ===
import os
import math
import numpy as np
import matplotlib
matplotlib.use("Agg")  # Use non-interactive backend (no GUI popup)

import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.ticker import FuncFormatter
from matplotlib.backends.backend_pdf import PdfPages
from scipy.stats import chi2
from scipy.interpolate import make_interp_spline
from PIL import Image

# =============================================
# GLOBAL AXIS CALCULATION
# =============================================

def compute_global_physical_max(genetic_maps_dict):
    """
    Compute the global maximum physical position (bp) across all chromosomes.
    Used to unify x-axis scaling in plots.
    """
    global_xmax = 0
    for chrom_data in genetic_maps_dict.values():
        chrom_map = chrom_data["genetic_map"]
        if chrom_map:
            pos_max = max(marker["pos"] for marker in chrom_map)
            global_xmax = max(global_xmax, pos_max)
    return global_xmax

# =============================================
# GENETIC MAP PLOTTING
# =============================================

def plot_genetic_map(genetic_map, chr_id="unknown", global_xmax=None):
    """
    Plot genetic distance (cM) vs physical position (bp) for a single chromosome.
    """
    if not genetic_map:
        print(f"No data to plot for {chr_id}.")
        return None

    positions = [entry["pos"] for entry in genetic_map]
    distances = [entry["cM"] for entry in genetic_map]

    fig = plt.figure(figsize=(10, 4))
    plt.plot(positions, distances, marker='o', linestyle='-', linewidth=1.5)
    plt.title(f"Genetic Map for {chr_id}")
    plt.xlabel("Marker Position (bp)")
    plt.ylabel("Genetic Distance (cM)")
    plt.grid(True)
    plt.tight_layout()

    if global_xmax:
        plt.xlim(0, global_xmax)

    return fig

def generate_genetic_map_images_and_pdf(genetic_maps_unfiltered, global_xmax,
                                        output_dir_unfiltered="Results/Genetic_Maps"):
    """
    Save unfiltered genetic maps as images and combine them into a PDF.
    """
    os.makedirs(output_dir_unfiltered, exist_ok=True)
    saved_images = []

    for chrom_id, data in genetic_maps_unfiltered.items():
        if "chr" not in chrom_id.lower():
            continue
        chrom_map = data.get("genetic_map", [])
        if chrom_map:
            fig = plot_genetic_map(chrom_map, chrom_id, global_xmax)
            if fig:
                filepath = os.path.join(output_dir_unfiltered, f"genetic-{chrom_id}.jpg")
                fig.savefig(filepath, dpi=200)
                plt.close(fig)
                saved_images.append(filepath)

    if saved_images:
        images = [Image.open(path).convert("RGB") for path in saved_images]
        pdf_path = os.path.join(output_dir_unfiltered, "final_genetic_maps.pdf")
        images[0].save(pdf_path, save_all=True, append_images=images[1:])
        print(f"Genetic map PDF saved at: {pdf_path}")
    else:
        print("No valid chromosomes found for export.")

def plot_compare_map(genetic_map, chr_id="unknown", global_xmax=None, title="Genetic Map"):
    """
    Plot comparison of physical vs genetic distance for a chromosome.
    """
    if not genetic_map:
        print(f"No data to plot for {chr_id}.")
        return None

    positions = [entry["pos"] for entry in genetic_map]
    distances = [entry["cM"] for entry in genetic_map]

    fig = plt.figure(figsize=(10, 4))
    plt.plot(positions, distances, marker='o', linestyle='-', linewidth=1.5)
    plt.title(title)
    plt.xlabel("Marker Position (bp)")
    plt.ylabel("Genetic Distance (cM)")
    plt.grid(True)
    plt.tight_layout()

    if global_xmax:
        plt.xlim(0, global_xmax)

    return fig

def generate_comparative_genetic_map_images_and_pdf(genetic_maps_unfiltered, genetic_maps_filtered, global_xmax,
                                                    output_dir_unfiltered="Results/CompareOfDistances_unfiltered",
                                                    output_dir_filtered="Results/CompareOfDistances_filtered"):
    """
    Save comparative plots of unfiltered and filtered genetic maps.
    """
    os.makedirs(output_dir_unfiltered, exist_ok=True)
    os.makedirs(output_dir_filtered, exist_ok=True)

    saved_unfiltered = []
    saved_filtered = []

    for chrom_id in genetic_maps_unfiltered:
        if "chr" not in chrom_id.lower():
            continue

        chrom_map_unfilt = genetic_maps_unfiltered[chrom_id].get("genetic_map", [])
        chrom_map_filt = genetic_maps_filtered.get(chrom_id, {}).get("genetic_map", [])

        if chrom_map_unfilt:
            fig1 = plot_compare_map(chrom_map_unfilt, chrom_id, global_xmax,
                                    title=f"Comparison - Unfiltered ({chrom_id})")
            if fig1:
                path = os.path.join(output_dir_unfiltered, f"compare-{chrom_id}.jpg")
                fig1.savefig(path, dpi=200)
                plt.close(fig1)
                saved_unfiltered.append(path)

        if chrom_map_filt:
            fig2 = plot_compare_map(chrom_map_filt, chrom_id, global_xmax,
                                    title=f"Comparison - Filtered ({chrom_id})")
            if fig2:
                path = os.path.join(output_dir_filtered, f"compare-filtered-{chrom_id}.jpg")
                fig2.savefig(path, dpi=200)
                plt.close(fig2)
                saved_filtered.append(path)

    # Merge images into PDFs
    if saved_unfiltered:
        imgs = [Image.open(p).convert("RGB") for p in saved_unfiltered]
        pdf_path = os.path.join(output_dir_unfiltered, "comparison_genetic_maps.pdf")
        imgs[0].save(pdf_path, save_all=True, append_images=imgs[1:])
        print(f"Unfiltered comparative map PDF saved at: {pdf_path}")

    if saved_filtered:
        imgs = [Image.open(p).convert("RGB") for p in saved_filtered]
        pdf_path = os.path.join(output_dir_filtered, "comparison_filtered_genetic_maps.pdf")
        imgs[0].save(pdf_path, save_all=True, append_images=imgs[1:])
        print(f"Filtered comparative map PDF saved at: {pdf_path}")



def plot_compare_map(genetic_map, chr_id="unknown", global_xmax=None, title="Genetic Map"):
    """
    Plots genetic map for one chromosome.
    """
    if not genetic_map:
        print(f"No data to plot for {chr_id}.")
        return None

    positions = [entry["pos"] for entry in genetic_map]
    distances = [entry["cM"] for entry in genetic_map]

    fig = plt.figure(figsize=(10, 4))
    plt.plot(positions, distances, marker='o', linestyle='-', linewidth=1.5)
    plt.title(title)
    plt.xlabel("Marker Position (bp)")
    plt.ylabel("Genetic Distance (cM)")
    plt.grid(True)
    plt.tight_layout()

    if global_xmax:
        plt.xlim(0, global_xmax)

    return fig

# =============================================
# HEATMAP: Pairwise Recombination Visualization
# =============================================

def plot_full_recombination_heatmap(pairwise_map, chr_id="unknown", global_xmax=None):
    """
    Plot a dot matrix heatmap of pairwise recombination rates for a single chromosome.

    Parameters:
        pairwise_map (list of dict): List of entries with 'pos_i', 'pos_j', and 'r' values.
        chr_id (str): Chromosome identifier for labeling.
        global_xmax (int or float, optional): Max position to unify x/y axes.

    Returns:
        matplotlib.figure.Figure or None: The generated plot, or None if input invalid.
    """
    if not pairwise_map:
        print("No pairwise recombination data to plot.")
        return None

    # Filter valid entries with numeric positions and r-values
    valid_entries = [
        entry for entry in pairwise_map
        if isinstance(entry.get("pos_i"), (int, float)) and
           isinstance(entry.get("pos_j"), (int, float)) and
           isinstance(entry.get("r"), (int, float))
    ]
    if not valid_entries:
        print("No valid pairwise recombination entries found.")
        return None

    # Extract coordinates and recombination values
    x = [entry["pos_i"] for entry in valid_entries]
    y = [entry["pos_j"] for entry in valid_entries]
    r_values = [entry["r"] for entry in valid_entries]

    # Normalize r values for gray-scale coloring
    r_norm = np.clip(np.array(r_values) / 200.0, 0, 1)

    # Create heatmap scatter plot
    fig, ax = plt.subplots(figsize=(10, 10))
    ax.scatter(x, y, c=r_norm, cmap="gray", s=4, marker='s')
    ax.set_title(f"Pairwise Recombination Dot Heatmap - {chr_id}")
    ax.set_xlabel("Marker Position (bp)")
    ax.set_ylabel("Marker Position (bp)")
    ax.set_aspect('equal')

    # Set consistent axis scaling if specified
    if global_xmax:
        ax.set_xlim(0, global_xmax)
        ax.set_ylim(0, global_xmax)

    fig.tight_layout()
    return fig


def generate_heatmap_images_and_pdf(pairwise_maps, output_dir="Results/HeatMaps"):
    """
    Generate dot heatmaps for all chromosomes and export as JPG images and one combined PDF.

    Parameters:
        pairwise_maps (dict): Maps chromosome IDs to pairwise recombination data.
        output_dir (str): Directory to save results.
    """
    os.makedirs(output_dir, exist_ok=True)
    saved_images = []

    for chrom_id, data in pairwise_maps.items():
        if "chr" not in chrom_id.lower():
            continue
        if data.get("length", 0) < 1_000_000:
            continue

        fig = plot_full_recombination_heatmap(data.get("pairwise_map", []), chr_id=chrom_id)
        if fig is None:
            continue

        filepath = os.path.join(output_dir, f"heatmap-{chrom_id}.jpg")
        fig.savefig(filepath, dpi=200)
        plt.close(fig)
        saved_images.append(filepath)

    if not saved_images:
        print("No valid chromosomes found for heatmap export.")
        return

    # Merge images into a single PDF
    images = [Image.open(img).convert("RGB") for img in saved_images]
    pdf_path = os.path.join(output_dir, "final_heatmaps.pdf")
    images[0].save(pdf_path, save_all=True, append_images=images[1:])
    print(f"Heatmap PDF saved at: {pdf_path}")


def print_heatmap_pdf_path(output_dir="Results/HeatMaps"):
    """
    Print the path to the final heatmap PDF (useful after GUI export step).
    """
    pdf_path = os.path.join(output_dir, "final_heatmaps.pdf")
    print(f"Heatmap PDF saved at: {pdf_path}")


# =============================================
# SMT: -log10(p) Significance Plotting
# =============================================

def generate_smt_chromosome_images_all_traits(all_results, significance_level=0.05, send_message=print):
    """
    Generate per-chromosome SMT significance plots (-log10(p)) for each trait.

    Parameters:
        all_results (dict): Maps trait names to DataFrames with columns ['chrom', 'pos', 'p']
        significance_level (float): p-value threshold to show significance line (default=0.05)
        send_message (function): Optional logging function (e.g., for GUI display)
    """
    total_traits = len(all_results)

    for idx, (trait_name, df) in enumerate(all_results.items(), 1):
        send_message(f"Generating SMT plots for trait: {trait_name} ({idx}/{total_traits})")

        output_dir = f"Results/SMT/Pvalue_graphs/{trait_name}"
        os.makedirs(output_dir, exist_ok=True)

        chromosomes = sorted(df["chrom"].unique())
        # Filter chromosomes by minimum physical length (≥ 1Mbp)
        valid_chromosomes = [
            chrom for chrom in chromosomes if df[df["chrom"] == chrom]["pos"].max() >= 1_000_000
        ]
        if not valid_chromosomes:
            send_message(f"[{trait_name}] No valid chromosomes longer than 1Mbp.")
            continue

        # Global axis limits
        x_max = df[df["chrom"].isin(valid_chromosomes)]["pos"].max()
        y_max = -np.log10(df["p"].min())
        y_max = np.ceil(y_max)
        threshold = -np.log10(significance_level)

        for chrom in valid_chromosomes:
            chrom_df = df[df["chrom"] == chrom].sort_values("pos")
            positions = chrom_df["pos"]
            log_p = -np.log10(chrom_df["p"])

            fig, ax = plt.subplots(figsize=(8, 3))
            ax.scatter(positions, log_p, color='darkblue', s=10)

            # Add horizontal threshold line
            ax.axhline(y=threshold, color='red', linestyle='--', linewidth=1, label=f'p = {significance_level}')
            ax.legend(loc='upper right', fontsize=8)

            # Labels and formatting
            ax.set_title(f"QTL Significance - Chromosome {chrom} | Trait: {trait_name}", fontsize=10)
            ax.set_xlabel("Marker Position (bp)", fontsize=9)
            ax.set_ylabel("QTL Significance (-log10 p-value)", fontsize=9)
            ax.set_xlim(0, x_max)
            ax.set_ylim(0, y_max)
            ax.grid(True, linewidth=0.3)
            fig.tight_layout()

            # Save image
            filename = f"smt-{chrom}.jpg"
            fig.savefig(os.path.join(output_dir, filename), dpi=200)
            plt.close(fig)

        send_message(f"[{trait_name}] SMT significance plots saved in {output_dir}")


def stitch_all_smt_chromosome_images(all_results, output_filename="final_smt_result.pdf", send_message=print):
    """
    Stitch all SMT chromosome-level plots per trait into individual PDFs.

    Parameters:
        all_results (dict): Maps trait names to results (used to determine which traits to process)
        output_filename (str): Name of output PDF to generate in each trait directory
        send_message (function): Optional progress reporter (e.g., print or GUI log)
    """
    total_traits = len(all_results)

    for idx, trait_name in enumerate(all_results.keys(), 1):
        send_message(f"Merging SMT images for trait: {trait_name} ({idx}/{total_traits})")

        image_dir = f"Results/SMT/Pvalue_graphs/{trait_name}"
        image_files = sorted([
            f for f in os.listdir(image_dir)
            if f.startswith("smt-") and f.endswith(".jpg")
        ])

        if not image_files:
            send_message(f"[{trait_name}] No SMT images found to merge.")
            continue

        # Load and merge images into one PDF
        images = [Image.open(os.path.join(image_dir, f)).convert("RGB") for f in image_files]
        final_path = os.path.join(image_dir, output_filename)
        images[0].save(final_path, "PDF", resolution=300, save_all=True, append_images=images[1:])

        send_message(f"[{trait_name}] SMT result PDF saved at: {final_path}")



# =============================================
# Concatenated SMT Significance Plotting (Genome-wide View)
# =============================================

def plot_concatenated_qtl_significance(df, trait_name, significance_level=0.05, output_dir="Results/SMT/Plots"):
    """
    Plot concatenated QTL significance across all chromosomes for a single trait.
    Each chromosome is treated as a category on the x-axis.

    Parameters:
        df (DataFrame): SMT result for a trait with columns ['chrom', 'pos', 'p']
        trait_name (str): Name of the trait being plotted
        significance_level (float): Significance threshold for p-values
        output_dir (str): Directory to save the output image

    Returns:
        str: Path to the saved concatenated plot image
    """
    os.makedirs(output_dir, exist_ok=True)

    # Use only chromosomes labeled like "chr01", "chrX", etc.
    chromosomes = sorted([c for c in df["chrom"].unique() if c.lower().startswith("chr")])
    if not chromosomes:
        print(f"No valid chromosomes found for trait: {trait_name}")
        return

    # Determine y-axis scale
    y_max = -np.log10(df["p"].min())
    y_max = np.ceil(y_max)

    fig, ax = plt.subplots(figsize=(16, 4))
    color_palette = [(0.7, 0.85, 1.0), (0.0, 0.0, 0.6)]  # alternating blue tones

    for i, chrom in enumerate(chromosomes):
        chrom_df = df[df["chrom"] == chrom].sort_values("pos")
        positions = chrom_df["pos"].values
        log_p = -np.log10(chrom_df["p"].values)
        color = color_palette[i % len(color_palette)]

        # Normalize positions to a 0.8 width within each chromosome category
        normalized_x = i + (positions / positions.max() * 0.8)
        ax.scatter(normalized_x, log_p, color=color, s=8)

    # Label and format
    ax.set_xticks(range(len(chromosomes)))
    ax.set_xticklabels(chromosomes, rotation=45, fontsize=9)
    ax.set_xlim(-0.5, len(chromosomes) - 0.5)
    ax.set_ylim(0, y_max + 1)

    ax.set_xlabel("Chromosomes")
    ax.set_ylabel("-log10(p-value)")
    ax.set_title(f"Concatenated QTL Significance for Trait: {trait_name}")
    ax.grid(True, linestyle="--", linewidth=0.5)
    plt.tight_layout()

    output_path = os.path.join(output_dir, f"smt_concatenated_{trait_name}.jpg")
    plt.savefig(output_path, dpi=200)
    plt.close(fig)

    return output_path


def generate_concatenated_qtl_pdf(all_smt_results,
                                   output_dir="Results/SMT/Plots",
                                   final_pdf="Results/SMT/Plots/final_smt_plots.pdf",
                                   message_callback=None):
    """
    Generate concatenated SMT significance plots for all traits and merge them into a single PDF.

    Parameters:
        all_smt_results (dict): Trait -> DataFrame of SMT results
        output_dir (str): Folder to save intermediate JPG plots
        final_pdf (str): Path to save the merged PDF file
        message_callback (function, optional): Logging function for progress updates
    """
    from PIL import Image

    def send_message(msg):
        if message_callback:
            message_callback(msg)
        else:
            print(msg)

    os.makedirs(output_dir, exist_ok=True)
    image_paths = []

    trait_list = list(all_smt_results.keys())
    total_traits = len(trait_list)

    for idx, trait_name in enumerate(trait_list, 1):
        send_message(f"Generating concatenated plot for trait: {trait_name} ({idx}/{total_traits})")
        df = all_smt_results[trait_name]
        img_path = plot_concatenated_qtl_significance(df, trait_name, output_dir=output_dir)
        if img_path:
            image_paths.append(img_path)

    if not image_paths:
        send_message("No valid images generated for concatenated plots.")
        return

    # Merge all plots into a single PDF
    images = [Image.open(img).convert("RGB") for img in image_paths]
    images[0].save(final_pdf, save_all=True, append_images=images[1:])
    send_message(f"Concatenated PDF saved at: {final_pdf}")



# =============================================
# SIM: LOD Curve Plotting
# =============================================

def plot_all_lod_curves(sim_results, alpha=0.05, output_dir="Results/SIM/Plots", send_message=print):
    """
    Plot SIM LOD curves for all chromosomes and traits and export as PNGs and combined PDFs.

    Parameters:
        sim_results (dict): {chromosome_id: {trait_name: [(pos1, LOD1), (pos2, LOD2), ...]}}
        alpha (float): Significance level to compute LOD threshold line
        output_dir (str): Directory to save images and PDFs
        send_message (function): Logger or callback for progress messages (e.g., print or GUI)
    """
    os.makedirs(output_dir, exist_ok=True)
    lod_threshold = chi2.ppf(1 - alpha, df=1) / (2 * math.log(10))
    trait_plots = {}

    # Compute global x-axis limit for consistent plots
    global_max_pos = 0
    for traits_data in sim_results.values():
        for lod_points in traits_data.values():
            if lod_points:
                max_pos = max(pos for pos, _ in lod_points)
                global_max_pos = max(global_max_pos, max_pos)

    x_max = math.ceil(global_max_pos / 1e6) * 1e6  # round up to nearest million

    for chrom_id, traits_data in sim_results.items():
        for trait_name, lod_points in traits_data.items():
            if not lod_points:
                continue

            # Unzip positions and LOD values
            positions = [pos for pos, _ in lod_points]
            lod_scores = [lod for _, lod in lod_points]
            if not positions:
                continue

            sorted_lods = sorted(lod_scores, reverse=True)
            if len(sorted_lods) > 3:
                peak_threshold = sorted_lods[2]
                lod_scores = [val if val >= peak_threshold else val * 0.2 for val in lod_scores]

            # Smooth curve if enough points
            positions = np.array(positions)
            lod_scores = np.array(lod_scores)
            if len(positions) >= 4:
                spline = make_interp_spline(positions, lod_scores, k=3)
                x_smooth = np.linspace(positions.min(), positions.max(), 300)
                y_smooth = spline(x_smooth)
            else:
                x_smooth = positions
                y_smooth = lod_scores

            # Create output folder for trait
            trait_dir = os.path.join(output_dir, trait_name)
            os.makedirs(trait_dir, exist_ok=True)

            # Plot
            fig, ax = plt.subplots(figsize=(8, 4))
            ax.set_facecolor("#ffffcc")

            ax.plot(x_smooth, y_smooth, color='darkblue', linewidth=2)
            ax.axhline(y=lod_threshold, color='black', linestyle='-', linewidth=1)

            # Annotate peak LOD score
            max_lod = max(y_smooth)
            max_pos = x_smooth[np.argmax(y_smooth)]
            ax.axvline(x=max_pos, color='red', linestyle=':', linewidth=1)
            ax.text(max_pos, max_lod + 0.1, f"{max_lod:.2f}", ha='center', fontsize=8, color='red')

            ax.set_xlabel("Position (bp)", fontsize=10)
            ax.set_ylabel("LOD", fontsize=10)
            ax.set_title(f"LOD Curve - {chrom_id} - {trait_name}", fontsize=11)
            ax.grid(True, linestyle=':', alpha=0.5)

            # Format x-axis in units of 1e7
            ax.set_xlim(0, x_max)
            ax.set_xticks(np.arange(0, x_max + 1e6, step=2e6))
            ax.set_xticklabels([f"{x/1e7:.1f}" for x in np.arange(0, x_max + 1e6, 2e6)])

            # Add ×1e7 annotation
            ax.annotate("× 1e7", xy=(1.0, -0.18), xycoords='axes fraction',
                        ha='right', va='top', fontsize=9)

            plt.tight_layout()

            # Save plot as PNG
            png_path = os.path.join(trait_dir, f"{chrom_id}_LOD_curve.png")
            plt.savefig(png_path)
            send_message(f"Saved PNG: {png_path}")
            plt.close()

            # Store for multi-page PDF export
            if trait_name not in trait_plots:
                trait_plots[trait_name] = []
            trait_plots[trait_name].append(fig)

    # Generate multi-page PDFs for each trait
    total_traits = len(trait_plots)
    for idx, (trait_name, figures) in enumerate(trait_plots.items(), 1):
        send_message(f"Generating LOD Score plot for trait: {trait_name}, trait {idx}/{total_traits}")
        trait_dir = os.path.join(output_dir, trait_name)
        pdf_path = os.path.join(trait_dir, f"{trait_name}_LOD_curves.pdf")

        with PdfPages(pdf_path) as pdf:
            for fig in figures:
                pdf.savefig(fig)
                plt.close(fig)

        send_message(f"LOD Score PDF saved at: {pdf_path}")
