#!/usr/bin/env python3
"""
LC-MS Peak Detection Pipeline
==============================-
Detects peaks (features) from raw LC-MS mzML files using:
  1. TIC/BPC construction from mzML
  2. Chromatographic peak detection via scipy.find_peaks
  3. Feature table generation (CSV)

Usage:
    python peak_detection.py --input data/raw/ --output data/peak_list.csv
"""

import argparse
import glob
import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.signal import find_peaks, savgol_filter

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger("peak_detection")


def parse_mzml(filepath: str, ms_level: int = 1):
    """Parse mzML file and return list of (rt, mz, intensity) tuples."""
    try:
        import pymzml
    except ImportError:
        sys.stderr.write("ERROR: pymzml not installed. Run: pip install pymzml\n")
        sys.exit(1)

    spectra = []
    with pymzml.run.Reader(filepath) as reader:
        for spectrum in reader:
            if spectrum.ms_level != ms_level:
                continue
            rt = spectrum.scan_time_in_minutes()
            mz_array = spectrum.mz
            int_array = spectrum.i
            if len(mz_array) == 0:
                continue
            spectra.append((rt, mz_array, int_array))
    log.info("  Parsed %d spectra from %s", len(spectra), Path(filepath).name)
    return spectra


def build_tic(spectra, use_bpc=False):
    """Build TIC (or BPC) from spectra list.
    
    TIC = sum of all intensities per scan
    BPC = max intensity per scan (base peak)
    """
    rt = []
    tic = []
    bpc = []
    for r, mz, ints in spectra:
        rt.append(r)
        tic.append(np.sum(ints))
        bpc.append(np.max(ints) if len(ints) > 0 else 0.0)
    rt = np.array(rt)
    if use_bpc:
        return rt, np.array(bpc)
    return rt, np.array(tic)


def detect_peaks_in_tic(rt, intensities, min_snr=3.0, min_peak_width=3,
                         noise_threshold=None):
    """Detect peaks in TIC/BPC chromatogram."""
    if noise_threshold is None:
        noise_threshold = np.percentile(intensities[intensities > 0], 10)
    
    # Smooth with Savitzky-Golay filter
    window = min(25, len(intensities) - 1)
    if window < 5:
        window = 5
    if window % 2 == 0:
        window -= 1
    smoothed = savgol_filter(intensities, window, polyorder=3)
    
    # Find peaks
    height_threshold = noise_threshold * min_snr
    peaks, properties = find_peaks(
        smoothed,
        height=height_threshold,
        prominence=noise_threshold * min_snr * 0.5,
        width=min_peak_width,
        distance=min_peak_width,
    )
    log.info("  Found %d peaks in TIC (SNR=%.1f)", len(peaks), min_snr)
    return peaks, rt, smoothed


def extract_feature_at_peak(spectra, peak_idx, rt, window_rt=0.1):
    """Extract all ions around a peak retention time."""
    peak_rt = rt[peak_idx]
    features = []
    for spec_idx, (r, mz, ints) in enumerate(spectra):
        if abs(r - peak_rt) > window_rt:
            continue
        for i, m in enumerate(mz):
            if i < len(ints):
                features.append((m, ints[i], r))
    return features


def detect_features_in_spectra(spectra, rt_ref, 
                                mz_tolerance_ppm=10,
                                intensity_threshold=500,
                                min_snr=3.0):
    """Detect features across all spectra (m/z features at each RT peak)."""
    noise_est = []
    for r, mz, ints in spectra:
        if np.percentile(ints, 10) > 0:
            noise_est.append(np.percentile(ints, 10))
    noise_level = np.median(noise_est) if noise_est else 100

    feature_list = []
    feature_id = 1
    
    for spec_idx, (r, mz, ints) in enumerate(spectra):
        above_threshold = ints > noise_level * min_snr
        if not np.any(above_threshold):
            continue
        for i, above in enumerate(above_threshold):
            if above:
                mz_val = mz[i]
                int_val = ints[i]
                snr = int_val / noise_level if noise_level > 0 else 0
                feature_list.append({
                    "feature_id": feature_id,
                    "mz": round(mz_val, 5),
                    "intensity": round(float(int_val), 1),
                    "rt_min": round(r, 2),
                    "snr": round(snr, 2),
                    "spectrum_idx": spec_idx,
                })
                feature_id += 1
    
    log.info("  Detected %d raw features (intensity > %.0f, SNR > %.1f)",
             len(feature_list), noise_level, min_snr)
    return pd.DataFrame(feature_list)


def cluster_features(df, mz_tolerance_ppm=10, rt_tolerance=0.1):
    """Cluster features with similar m/z and RT (simple window-based clustering)."""
    if df.empty:
        return df
    
    df = df.sort_values(["mz", "rt_min"]).reset_index(drop=True)
    clusters = []
    current_cluster = [df.iloc[0]]
    
    for i in range(1, len(df)):
        row = df.iloc[i]
        prev = current_cluster[-1]
        mz_diff_ppm = abs(row["mz"] - prev["mz"]) / prev["mz"] * 1e6
        rt_diff = abs(row["rt_min"] - prev["rt_min"])
        if mz_diff_ppm < mz_tolerance_ppm and rt_diff < rt_tolerance:
            current_cluster.append(row)
        else:
            clusters.append(current_cluster)
            current_cluster = [row]
    if current_cluster:
        clusters.append(current_cluster)
    
    aggregated = []
    for cluster_id, cluster in enumerate(clusters, 1):
        cdf = pd.DataFrame(cluster)
        aggregated.append({
            "feature_id": cluster_id,
            "mz_min": cdf["mz"].min(),
            "mz_max": cdf["mz"].max(),
            "mz_med": round(cdf["mz"].median(), 5),
            "intensity_max": cdf["intensity"].max(),
            "intensity_sum": round(cdf["intensity"].sum(), 1),
            "rt_min": round(cdf["rt_min"].min(), 2),
            "rt_max": round(cdf["rt_min"].max(), 2),
            "rt_med": round(cdf["rt_min"].median(), 2),
            "snr_max": round(cdf["snr"].max(), 2),
            "n_spectra": len(cdf),
        })
    
    log.info("  Clustered into %d features (PPM tol=%d, RT tol=%.1f min)",
             len(aggregated), mz_tolerance_ppm, rt_tolerance)
    return pd.DataFrame(aggregated)


def process_mzml_file(filepath, args):
    """Process a single mzML file."""
    log.info("Processing: %s", Path(filepath).name)
    spectra = parse_mzml(filepath, ms_level=1)
    if not spectra:
        log.warning("No MS1 spectra found in %s", filepath)
        return None
    
    # Build TIC and detect peaks
    rt, tic = build_tic(spectra, use_bpc=False)
    peak_indices, rt_smooth, tic_smooth = detect_peaks_in_tic(
        rt, tic,
        min_snr=args.min_snr,
        min_peak_width=args.min_peak_width,
    )
    
    # Detect features across all spectra
    df = detect_features_in_spectra(
        spectra, rt,
        intensity_threshold=args.noise_threshold,
        min_snr=args.min_snr,
    )
    
    if df.empty:
        return None
    
    # Cluster features
    df_clustered = cluster_features(
        df,
        mz_tolerance_ppm=args.mz_tolerance_ppm,
        rt_tolerance=args.rt_tolerance,
    )
    
    df_clustered["source_file"] = Path(filepath).name
    log.info("  Final feature count: %d", len(df_clustered))
    return df_clustered


def main():
    parser = argparse.ArgumentParser(
        description="LC-MS peak detection from mzML files"
    )
    parser.add_argument("--input", "-i", required=True,
                        help="Input mzML file or directory")
    parser.add_argument("--output", "-o", required=True,
                        help="Output CSV file for peak list")
    parser.add_argument("--mzml_pattern", default="*.mzML",
                        help="Glob pattern for mzML files in directory mode")
    parser.add_argument("--noise_threshold", type=float, default=500,
                        help="Minimum absolute intensity (default: 500)")
    parser.add_argument("--min_snr", type=float, default=3.0,
                        help="Minimum signal-to-noise ratio (default: 3.0)")
    parser.add_argument("--min_peak_width", type=int, default=3,
                        help="Minimum peak width in spectra (default: 3)")
    parser.add_argument("--mz_tolerance_ppm", type=int, default=10,
                        help="m/z clustering tolerance in PPM (default: 10)")
    parser.add_argument("--rt_tolerance", type=float, default=0.1,
                        help="RT clustering tolerance in minutes (default: 0.1)")
    parser.add_argument("--bpc", action="store_true",
                        help="Use base peak chromatogram instead of TIC")
    parser.add_argument("--verbose", "-v", action="store_true")
    
    args = parser.parse_args()
    
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    
    input_path = Path(args.input)
    mzml_files = []
    
    if input_path.is_file():
        mzml_files = [input_path]
    elif input_path.is_dir():
        mzml_files = sorted(input_path.glob(args.mzml_pattern))
        if not mzml_files:
            log.error("No mzML files found in %s with pattern %s",
                      input_path, args.mzml_pattern)
            sys.exit(1)
    else:
        log.error("Input path not found: %s", input_path)
        sys.exit(1)
    
    log.info("Found %d mzML file(s) to process", len(mzml_files))
    
    all_features = []
    for fpath in mzml_files:
        result = process_mzml_file(str(fpath), args)
        if result is not None:
            all_features.append(result)
    
    if not all_features:
        log.warning("No features detected in any file!")
        sys.exit(0)
    
    combined = pd.concat(all_features, ignore_index=True)
    
    # Filter by minimum intensity
    before = len(combined)
    combined = combined[combined["intensity_max"] >= args.noise_threshold]
    log.info("Filtered %d features below noise threshold (%d retained)",
             before - len(combined), len(combined))
    
    # Save
    combined.to_csv(args.output, index=False)
    log.info("Saved peak list to: %s", args.output)
    log.info("Total features: %d", len(combined))
    print(f"\nPeak detection complete!")
    print(f"  Files processed: {len(mzml_files)}")
    print(f"  Features detected: {len(combined)}")
    print(f"  Output: {args.output}")


if __name__ == "__main__":
    main()
