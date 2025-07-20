
"""
plot_utils.py

This module contains all plotting utilities for QTL Mapping:
- Genetic Map plotting (filtered / unfiltered)
- Pairwise recombination heatmap plotting
- SMT significance plotting (chromosome and genome-wide)
- Automatic export of chromosome images to JPG and final PDFs.
"""

import os
import numpy as np

# ✅ Prevent matplotlib from using GUI backends (important for threading)
import matplotlib
matplotlib.use("Agg")  # Use non-interactive backend (no GUI popup)

import matplotlib.pyplot as plt
import matplotlib.cm as cm
from PIL import Image
from scipy.stats import chi2
from matplotlib.backends.backend_pdf import PdfPages
import math

##############################################
# Global axis computation
##############################################

def compute_global_physical_max(genetic_maps_dict):
    """
    Compute maximum marker position across all chromosomes (bp).
    """
    global_xmax = 0
    for chrom_data in genetic_maps_dict.values():
        chrom_map = chrom_data["genetic_map"]
        if chrom_map:
            pos_max = max(marker["pos"] for marker in chrom_map)
            global_xmax = max(global_xmax, pos_max)
    return global_xmax

##############################################
# Genetic Map Plotting
##############################################

def plot_genetic_map(genetic_map, chr_id="unknown", global_xmax=None):
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
    plt.title(f"Genetic Map for {chr_id}")
    plt.xlabel("Marker Position (bp)")
    plt.ylabel("Genetic Distance (cM)")
    plt.grid(True)
    plt.tight_layout()

    if global_xmax:
        plt.xlim(0, global_xmax)

    return fig


def generate_genetic_map_images_and_pdf(genetic_maps, global_xmax, output_dir="Results/Genetic_Maps"):
    os.makedirs(output_dir, exist_ok=True)
    saved_images = []

    for chrom_id, chrom_data in genetic_maps.items():
        if "chr" not in chrom_id.lower(): continue
        chrom_map = chrom_data["genetic_map"]
        if not chrom_map: continue

        fig = plot_genetic_map(chrom_map, chrom_id, global_xmax)
        if fig is None: continue

        filename = f"genetic-{chrom_id}.jpg"
        filepath = os.path.join(output_dir, filename)
        fig.savefig(filepath, dpi=200)
        plt.close(fig)
        saved_images.append(filepath)

    if not saved_images:
        print("No valid chromosomes found for genetic map export.")
        return

    images = [Image.open(img).convert("RGB") for img in saved_images]
    pdf_path = os.path.join(output_dir, "final_genetic_maps.pdf")
    images[0].save(pdf_path, save_all=True, append_images=images[1:])
    print(f"Genetic map PDF saved at: {pdf_path}")


def generate_filtered_genetic_map_images_and_pdf(genetic_maps_filtered, global_xmax, output_dir="Results/Genetic_Maps_Filtered"):
    os.makedirs(output_dir, exist_ok=True)
    saved_images = []

    for chrom_id, chrom_data in genetic_maps_filtered.items():
        if "chr" not in chrom_id.lower(): continue
        chrom_map = chrom_data["genetic_map"]
        if not chrom_map: continue

        fig = plot_genetic_map(chrom_map, chrom_id, global_xmax)
        if fig is None: continue

        filename = f"filtered-genetic-{chrom_id}.jpg"
        filepath = os.path.join(output_dir, filename)
        fig.savefig(filepath, dpi=200)
        plt.close(fig)
        saved_images.append(filepath)

    if not saved_images:
        print("No valid chromosomes found for filtered map export.")
        return

    images = [Image.open(img).convert("RGB") for img in saved_images]
    pdf_path = os.path.join(output_dir, "final_filtered_genetic_maps.pdf")
    images[0].save(pdf_path, save_all=True, append_images=images[1:])
    print(f"Filtered genetic map PDF saved at: {pdf_path}")

##############################################
# Heatmap (Pairwise Recombination Map)
##############################################

def plot_full_recombination_heatmap(pairwise_map, chr_id="unknown", global_xmax=None):
    """
    Plots a dot for each pair of marker positions on a chromosome using actual pairwise recombination rates.
    """
    if not pairwise_map:
        print("No pairwise recombination data to plot.")
        return None

    valid_entries = [
        entry for entry in pairwise_map
        if isinstance(entry.get("pos_i"), (int, float)) and
           isinstance(entry.get("pos_j"), (int, float)) and
           isinstance(entry.get("r"), (int, float))
    ]
    if len(valid_entries) < 1:
        print("No valid pairwise recombination entries found.")
        return None

    x = [entry["pos_i"] for entry in valid_entries]
    y = [entry["pos_j"] for entry in valid_entries]
    r_values = [entry["r"] for entry in valid_entries]

    r_norm = np.clip(np.array(r_values) / 200.0, 0, 1)

    fig, ax = plt.subplots(figsize=(10, 10))
    scatter = ax.scatter(x, y, c=r_norm, cmap="gray", s=4, marker='s')
    ax.set_title(f"Pairwise Recombination Dot Heatmap - {chr_id}")
    ax.set_xlabel("Marker Position (bp)")
    ax.set_ylabel("Marker Position (bp)")
    ax.set_aspect('equal')

    if global_xmax:
        ax.set_xlim(0, global_xmax)
        ax.set_ylim(0, global_xmax)

    fig.tight_layout()

    return fig

def generate_heatmap_images_and_pdf(pairwise_maps, output_dir="Results/HeatMaps"):
    os.makedirs(output_dir, exist_ok=True)

    saved_images = []

    for chrom_id, pairwise_map in pairwise_maps.items():
        if "chr" not in chrom_id.lower():
            continue

        if pairwise_map["length"] < 1_000_000:
            continue

        fig = plot_full_recombination_heatmap(pairwise_map["pairwise_map"], chr_id=chrom_id)
        if fig is None:
            continue

        filename = f"heatmap-{chrom_id}.jpg"
        filepath = os.path.join(output_dir, filename)
        fig.savefig(filepath, dpi=200)
        plt.close(fig)

        saved_images.append(filepath)

    if not saved_images:
        print("No valid chromosomes found for heatmap export.")
        return

    # Merge into PDF
    images = [Image.open(img).convert("RGB") for img in saved_images]
    pdf_path = os.path.join(output_dir, "final_heatmaps.pdf")
    images[0].save(pdf_path, save_all=True, append_images=images[1:])
    print(f"Heatmap PDF saved at: {pdf_path}")

def print_heatmap_pdf_path(output_dir="Results/HeatMaps"):
    """
    Print the path to the heatmap PDF in the specified output directory.
    """
    pdf_path = os.path.join(output_dir, "final_heatmaps.pdf")
    print(f"Heatmap PDF saved at: {pdf_path}")


##############################################
# SMT -log10(p) Significance Plotting
##############################################

def generate_smt_chromosome_images(df, trait_name, significance_level=0.05):
    """
    Generate SMT significance plots for each chromosome separately.
    """
    output_dir = f"Results/SMT/{trait_name}"
    os.makedirs(output_dir, exist_ok=True)

    chromosomes = sorted(df["chrom"].unique())
    valid_chromosomes = [
        chrom for chrom in chromosomes if df[df['chrom'] == chrom]['pos'].max() >= 1_000_000
    ]
    if not valid_chromosomes:
        print("No valid chromosomes longer than 1Mbp.")
        return

    x_max = df[df['chrom'].isin(valid_chromosomes)].groupby('chrom')['pos'].max().max()
    y_max = -np.log10(df['p'].min())
    y_max = np.ceil(y_max)

    for chrom in valid_chromosomes:
        chrom_df = df[df['chrom'] == chrom].sort_values("pos")
        positions = chrom_df["pos"]
        log_p = -np.log10(chrom_df["p"])
        threshold = -np.log10(significance_level)
        colors = ['darkblue' if p > threshold else 'lightblue' for p in log_p]

        fig = plt.figure(figsize=(8, 3))
        plt.scatter(positions, log_p, c=colors, s=10)
        plt.title(f"QTL Significance across {chrom} for trait: {trait_name}", fontsize=10)
        plt.xlabel("Marker Position (bp)", fontsize=9)
        plt.ylabel("QTL Significance Based on P-value", fontsize=9)
        plt.xlim(0, x_max)
        plt.ylim(0, y_max)
        plt.grid(True, linewidth=0.3)
        plt.tight_layout()

        filename = f"smt-{chrom}.jpg"
        plt.savefig(os.path.join(output_dir, filename), dpi=200)
        plt.close(fig)


def stitch_smt_chromosome_images_pdf(trait_name, output_filename="final_smt_result.pdf"):
    """
    Merge SMT chromosome images into final PDF.
    """
    image_dir = f"Results/SMT/{trait_name}"
    files = sorted([
        f for f in os.listdir(image_dir)
        if f.startswith("smt-") and f.endswith(".jpg")
    ])
    images = [Image.open(os.path.join(image_dir, f)).convert("RGB") for f in files]
    final_path = os.path.join(image_dir, output_filename)
    images[0].save(final_path, "PDF", resolution=300, save_all=True, append_images=images[1:])
    return final_path

##############################################
# Mini SMT Graphs combined Plotting
##############################################
def plot_concatenated_qtl_significance(df, trait_name, significance_level=0.05, output_dir="Results/SMT/Plots"):
    """
    Concatenated QTL significance across chromosomes for a given trait.
    Each chromosome is plotted as a separate category on the x-axis.
    """
    import matplotlib.colors as mcolors

    os.makedirs(output_dir, exist_ok=True)

    # Filter only chromosomes starting with "chr"
    chromosomes = sorted([c for c in df['chrom'].unique() if c.lower().startswith("chr")])
    if not chromosomes:
        print(f"No valid chromosomes found for trait: {trait_name}")
        return

    # Calculate global max Y axis for scaling
    y_max = -np.log10(df['p'].min())
    y_max = np.ceil(y_max)

    fig, ax = plt.subplots(figsize=(16, 4))
    
    color_palette = [(0.7, 0.85, 1.0), (0.0, 0.0, 0.6)]

    for i, chrom in enumerate(chromosomes):
        chrom_df = df[df['chrom'] == chrom].sort_values("pos")
        positions = chrom_df['pos'].values
        log_p = -np.log10(chrom_df['p'].values)
        color = color_palette[i % 2]

        # Plot relative positions within chromosome
        ax.scatter([i + (pos / positions.max() * 0.8) for pos in positions], log_p, color=color, s=8)

    # Set chromosome labels as x-axis categories
    ax.set_xticks(range(len(chromosomes)))
    ax.set_xticklabels(chromosomes, rotation=45, fontsize=9)

    ax.set_xlim(-0.5, len(chromosomes)-0.5)
    ax.set_ylim(0, y_max + 1)
    ax.set_xlabel("Chromosomes")
    ax.set_ylabel("-log10(p-value)")
    ax.set_title(f"Concatenated QTL Significance for Trait: {trait_name}")
    ax.grid(True, linestyle="--", linewidth=0.5)
    plt.tight_layout()

    output_file = os.path.join(output_dir, f"smt_concatenated_{trait_name}.jpg")
    plt.savefig(output_file, dpi=200)
    plt.close(fig)

    return output_file  # so we can reuse it for merging later

def generate_concatenated_qtl_pdf(all_smt_results, output_dir="Results/SMT/Plots", final_pdf="Results/SMT/Plots/final_smt_plots.pdf", message_callback=None):
    """
    Generate concatenated QTL significance plots for all traits and merge into PDF.
    message_callback: optional function to stream progress messages.
    """
    from PIL import Image

    def send_message(msg):
        if message_callback is not None:
            message_callback(msg)
        else:
            print(msg)

    os.makedirs(output_dir, exist_ok=True)
    image_paths = []

    trait_list = list(all_smt_results.keys())
    total_traits = len(trait_list)
    for idx, trait_name in enumerate(trait_list, 1):
        send_message(f"Generating concatenated plot for trait: {trait_name}, trait {idx}/{total_traits}")
        df = all_smt_results[trait_name]
        img_path = plot_concatenated_qtl_significance(df, trait_name, output_dir=output_dir)
        if img_path:
            image_paths.append(img_path)

    if not image_paths:
        send_message("No valid images generated for concatenated plots.")
        return

    # Merge into PDF
    images = [Image.open(img).convert("RGB") for img in image_paths]
    images[0].save(final_pdf, save_all=True, append_images=images[1:])
    send_message(f"Concatenated PDF saved at: {final_pdf}")


def plot_lod_curve(sim_results, chromosome_id, alpha=0.05, output_dir="Results/SIM/Plots"):
    if chromosome_id not in sim_results:
        print(f"Chromosome {chromosome_id} not found in results.")
        return

    os.makedirs(output_dir, exist_ok=True)
    lod_threshold = chi2.ppf(1 - alpha, df=1) / (2 * math.log(10))

    for trait_name, lod_points in sim_results[chromosome_id].items():
        positions = [pos for pos, lod in lod_points]
        lod_scores = [lod for pos, lod in lod_points]
        if not positions:
            continue

        plt.figure(figsize=(10, 5))
        plt.plot(positions, lod_scores, label=f"LOD Curve: {trait_name}", color="navy")
        plt.axhline(y=lod_threshold, color="black", linestyle="--", label=f"Threshold (\u03b1={alpha})")

        max_lod = max(lod_scores)
        max_pos = positions[lod_scores.index(max_lod)]
        plt.axvline(x=max_pos, color="red", linestyle=":", label=f"Max LOD @ {max_pos:.1f} bp")
        plt.text(max_pos, max_lod + 0.1, f"{max_lod:.2f}", ha='center', fontsize=8, color='red')

        plt.xlabel("Position (bp)")
        plt.ylabel("LOD Score")
        plt.title(f"LOD Score Curve - {chromosome_id} - {trait_name}")
        plt.legend()
        plt.grid(True, linestyle=':', alpha=0.5)
        plt.tight_layout()

        save_path = os.path.join(output_dir, f"{chromosome_id}_{trait_name}_LOD_curve.png")
        plt.savefig(save_path)
        plt.close()
        print(f"Saved plot to {save_path}")

def plot_all_lod_curves(sim_results, alpha=0.05, output_dir="Results/SIM/Plots"):
    os.makedirs(output_dir, exist_ok=True)
    pdf_dir = os.path.join(output_dir, "PdfPerTrait")
    os.makedirs(pdf_dir, exist_ok=True)
    lod_threshold = chi2.ppf(1 - alpha, df=1) / (2 * math.log(10))
    trait_plots = {}

    for chrom_id, traits_data in sim_results.items():
        for trait_name, lod_points in traits_data.items():
            positions = [pos for pos, lod in lod_points]
            lod_scores = [lod for pos, lod in lod_points]
            if not positions:
                continue

            plt.figure(figsize=(10, 5))
            plt.plot(positions, lod_scores, label="LOD Curve", color="navy")
            plt.axhline(y=lod_threshold, color="black", linestyle="--", label=f"Threshold (\u03b1={alpha})")

            max_lod = max(lod_scores)
            max_pos = positions[lod_scores.index(max_lod)]
            plt.axvline(x=max_pos, color="red", linestyle=":", label=f"Max LOD @ {max_pos:.1f} bp")
            plt.text(max_pos, max_lod + 0.1, f"{max_lod:.2f}", ha='center', fontsize=8, color='red')

            plt.xlabel("Position (bp)")
            plt.ylabel("LOD Score")
            plt.title(f"LOD Curve - {chrom_id} - {trait_name}")
            plt.legend()
            plt.grid(True, linestyle=':', alpha=0.5)
            plt.tight_layout()

            png_path = os.path.join(output_dir, f"{chrom_id}_{trait_name}_LOD_curve.png")
            plt.savefig(png_path)
            print(f"Saved PNG: {png_path}")

            if trait_name not in trait_plots:
                trait_plots[trait_name] = []
            trait_plots[trait_name].append(plt.gcf())
            plt.close()

    for trait_name, figures in trait_plots.items():
        pdf_path = os.path.join(pdf_dir, f"{trait_name}_LOD_curves.pdf")
        with PdfPages(pdf_path) as pdf:
            for fig in figures:
                pdf.savefig(fig)
                plt.close(fig)
        print(f"Saved PDF: {pdf_path}")


