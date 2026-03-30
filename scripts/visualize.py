#!/usr/bin/env python3
"""
LC-MS Visualization Toolkit
===========================
Generates publication-quality figures from LC-MS peak detection results.

Outputs:
  - Base Peak Chromatogram (BPC)
  - Total Ion Chromatogram (TIC)
  - m/z vs Retention Time scatter plot
  - Peak intensity heatmap
  - Feature statistics summary

Usage:
    python visualize.py --input data/features.csv --output results/figures/
"""

import argparse
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")  # Non-interactive backend
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import pandas as pd
import seaborn as sns

# Publication-quality style
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.size": 11,
    "axes.titlesize": 13,
    "axes.labelsize": 11,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "legend.fontsize": 10,
    "figure.dpi": 150,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "axes.spines.top": False,
    "axes.spines.right": False,
})


def plot_mz_vs_rt(df, output_path):
    """Scatter plot: m/z vs RT, colored by intensity."""
    fig, ax = plt.subplots(figsize=(8, 6))
    
    sc = ax.scatter(
        df["rt_med"], df["mz_med"],
        c=np.log10(df["intensity_max"] + 1),
        cmap="viridis",
        s=df["intensity_max"] / df["intensity_max"].max() * 100 + 5,
        alpha=0.7,
        edgecolors="white",
        linewidths=0.5,
    )
    cbar = plt.colorbar(sc, ax=ax, label="log₁₀(intensity)")
    ax.set_xlabel("Retention Time (min)")
    ax.set_ylabel("m/z (Th)")
    ax.set_title("LC-MS Feature Map: m/z vs RT")
    
    cbar.ax.set_ylabel("log₁₀(intensity + 1)", fontsize=10)
    ax.xaxis.set_major_locator(ticker.MultipleLocator(1.0))
    ax.yaxis.set_major_locator(ticker.MultipleLocator(100))
    
    plt.savefig(output_path)
    plt.close()
    print(f"  Saved: {output_path}")


def plot_intensity_histogram(df, output_path):
    """Histogram of feature intensities."""
    fig, axes = plt.subplots(1, 2, figsize=(12, 4))
    
    intensities = df["intensity_max"]
    ax1, ax2 = axes
    
    # Linear scale
    ax1.hist(intensities, bins=50, color="steelblue", edgecolor="white", alpha=0.8)
    ax1.set_xlabel("Peak Intensity (BPC)")
    ax1.set_ylabel("Feature Count")
    ax1.set_title("Feature Intensity Distribution")
    ax1.axvline(intensities.median(), color="red", linestyle="--",
                label=f"Median={intensities.median():.0f}")
    ax1.legend()
    
    # Log scale
    ax2.hist(np.log10(intensities + 1), bins=50, color="darkorange",
             edgecolor="white", alpha=0.8)
    ax2.set_xlabel("log₁₀(Intensity + 1)")
    ax2.set_ylabel("Feature Count")
    ax2.set_title("Feature Intensity Distribution (log scale)")
    ax2.axvline(np.log10(intensities.median() + 1), color="red",
                linestyle="--", label=f"Median={np.log10(intensities.median()+1):.2f}")
    ax2.legend()
    
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()
    print(f"  Saved: {output_path}")


def plot_snr_distribution(df, output_path):
    """SNR distribution of detected features."""
    fig, ax = plt.subplots(figsize=(7, 5))
    
    snr = df["snr_max"] if "snr_max" in df.columns else df["snr"]
    ax.hist(snr, bins=50, color="#2ecc71", edgecolor="white", alpha=0.8)
    ax.set_xlabel("Signal-to-Noise Ratio (SNR)")
    ax.set_ylabel("Feature Count")
    ax.set_title("Feature SNR Distribution")
    ax.axvline(snr.median(), color="red", linestyle="--",
               label=f"Median SNR = {snr.median():.1f}")
    ax.legend()
    
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()
    print(f"  Saved: {output_path}")


def plot_feature_stats(df, output_path):
    """Summary statistics panel."""
    n_features = len(df)
    mz_range = (df["mz_min"].min(), df["mz_max"].max())
    rt_range = (df["rt_min"].min(), df["rt_max"].max())
    intensity_range = (df["intensity_max"].min(), df["intensity_max"].max())
    
    fig, axes = plt.subplots(2, 2, figsize=(12, 9))
    
    # m/z distribution
    ax = axes[0, 0]
    ax.hist(df["mz_med"], bins=50, color="steelblue", edgecolor="white", alpha=0.8)
    ax.set_xlabel("m/z (Th)")
    ax.set_ylabel("Count")
    ax.set_title("m/z Distribution")
    
    # RT distribution
    ax = axes[0, 1]
    ax.hist(df["rt_med"], bins=50, color="darkorange", edgecolor="white", alpha=0.8)
    ax.set_xlabel("Retention Time (min)")
    ax.set_ylabel("Count")
    ax.set_title("RT Distribution")
    
    # Intensity vs m/z
    ax = axes[1, 0]
    sc = ax.scatter(df["mz_med"], np.log10(df["intensity_max"] + 1),
                    c=df["rt_med"], cmap="plasma", alpha=0.6, s=20)
    ax.set_xlabel("m/z (Th)")
    ax.set_ylabel("log₁₀(Intensity + 1)")
    ax.set_title("Intensity vs m/z (colored by RT)")
    plt.colorbar(sc, ax=ax, label="RT (min)")
    
    # Top features table
    ax = axes[1, 1]
    ax.axis("off")
    top_features = df.nlargest(15, "intensity_max")[
        ["feature_id", "mz_med", "rt_med", "intensity_max", "snr_max"]
    ]
    table_data = [
        [str(int(r.feature_id)), f"{r.mz_med:.4f}", f"{r.rt_med:.2f}",
         f"{r.intensity_max:.0f}", f"{r.snr_max:.1f}"]
        for _, r in top_features.iterrows()
    ]
    col_labels = ["ID", "m/z", "RT (min)", "Intensity", "SNR"]
    table = ax.table(
        cellText=table_data,
        colLabels=col_labels,
        loc="center",
        cellLoc="center",
    )
    table.auto_set_font_size(False)
    table.set_fontsize(9)
    table.scale(1.2, 1.5)
    ax.set_title("Top 15 Features by Intensity", pad=20)
    
    # Add summary text
    summary_text = (
        f"Total features: {n_features}\n"
        f"m/z range: {mz_range[0]:.2f} – {mz_range[1]:.2f}\n"
        f"RT range: {rt_range[0]:.2f} – {rt_range[1]:.2f} min\n"
        f"Intensity range: {intensity_range[0]:.0f} – {intensity_range[1]:.0f}"
    )
    fig.text(0.02, 0.02, summary_text, fontsize=9, family="monospace",
             verticalalignment="bottom")
    
    plt.suptitle("LC-MS Feature Detection Summary", fontsize=14, fontweight="bold")
    plt.tight_layout(rect=[0, 0.05, 1, 0.97])
    plt.savefig(output_path)
    plt.close()
    print(f"  Saved: {output_path}")


def plot_bpc_from_mzml(filepath, output_path, n_top=20):
    """Plot Base Peak Chromatogram from mzML file."""
    try:
        import pymzml
    except ImportError:
        sys.stderr.write("pymzml not installed. Skipping BPC plot.\n")
        return
    
    rt_list = []
    bpc_list = []
    
    with pymzml.run.Reader(filepath) as reader:
        for spectrum in reader:
            if spectrum.ms_level != 1:
                continue
            rt_list.append(spectrum.scan_time_in_minutes())
            ints = spectrum.i
            bpc_list.append(np.max(ints) if len(ints) > 0 else 0.0)
    
    if not rt_list:
        print(f"  Warning: No MS1 spectra in {filepath}")
        return
    
    fig, ax = plt.subplots(figsize=(12, 4))
    ax.plot(rt_list, bpc_list, color="steelblue", linewidth=0.8)
    ax.fill_between(rt_list, bpc_list, alpha=0.3, color="steelblue")
    ax.set_xlabel("Retention Time (min)")
    ax.set_ylabel("Base Peak Intensity")
    ax.set_title(f"Base Peak Chromatogram: {Path(filepath).name}")
    ax.xaxis.set_major_locator(ticker.MultipleLocator(1.0))
    
    plt.savefig(output_path)
    plt.close()
    print(f"  Saved: {output_path}")


def plot_heatmap(df, output_path):
    """Heatmap of top features across samples (if multiple sources)."""
    if "source_file" not in df.columns:
        # Single sample — plot m/z vs RT heatmap
        fig, ax = plt.subplots(figsize=(10, 6))
        hb = ax.hexbin(df["rt_med"], df["mz_med"],
                       C=np.log10(df["intensity_max"] + 1),
                       gridsize=40, cmap="YlOrRd", reduce_C_function=np.mean)
        ax.set_xlabel("Retention Time (min)")
        ax.set_ylabel("m/z (Th)")
        ax.set_title("Feature Density: m/z vs RT")
        plt.colorbar(hb, ax=ax, label="log₁₀(Intensity + 1)")
        plt.savefig(output_path)
        plt.close()
    else:
        # Multi-sample heatmap
        pivot = df.pivot_table(
            index="feature_id",
            columns="source_file",
            values="intensity_max",
            aggfunc="max",
            fill_value=0,
        )
        top_n = min(50, len(pivot))
        pivot_top = pivot.nlargest(top_n, pivot.columns[0])
        
        fig, ax = plt.subplots(figsize=(10, 8))
        sns.heatmap(
            np.log10(pivot_top + 1),
            cmap="viridis",
            ax=ax,
            xticklabels=True,
            yticklabels=False,
            cbar_kws={"label": "log₁₀(Intensity + 1)"},
        )
        ax.set_xlabel("Sample")
        ax.set_ylabel("Feature (top 50)")
        ax.set_title("Feature Intensity Heatmap")
        plt.savefig(output_path)
        plt.close()
    
    print(f"  Saved: {output_path}")


def main():
    parser = argparse.ArgumentParser(description="LC-MS visualization")
    parser.add_argument("--input", "-i", required=True,
                        help="Input CSV feature table")
    parser.add_argument("--output", "-o", required=True,
                        help="Output directory for figures")
    parser.add_argument("--mzml", "-m", default=None,
                        help="Optional mzML file for BPC/TIC plots")
    parser.add_argument("--top_n", type=int, default=20,
                        help="Number of top features to label (default: 20)")
    
    args = parser.parse_args()
    
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    df = pd.read_csv(args.input)
    print(f"Loaded {len(df)} features from {args.input}")
    
    print("\nGenerating figures...")
    
    plot_mz_vs_rt(df, output_dir / "mz_vs_rt.png")
    plot_intensity_histogram(df, output_dir / "intensity_distribution.png")
    plot_snr_distribution(df, output_dir / "snr_distribution.png")
    plot_feature_stats(df, output_dir / "feature_summary.png")
    plot_heatmap(df, output_dir / "feature_heatmap.png")
    
    if args.mzml:
        plot_bpc_from_mzml(args.mzml, output_dir / "bpc.png")
    
    print(f"\nAll figures saved to: {output_dir}/")
    print(f"  mz_vs_rt.png             - m/z vs RT scatter plot")
    print(f"  intensity_distribution.png - Feature intensity histograms")
    print(f"  snr_distribution.png     - SNR distribution")
    print(f"  feature_summary.png      - Summary statistics panel")
    print(f"  feature_heatmap.png      - Feature density/sample heatmap")


if __name__ == "__main__":
    main()
