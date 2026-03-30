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
    """Parse mzML file using Python's xml library and return list of (rt, mz, intensity) tuples."""
    import base64
    import struct
    import xml.etree.ElementTree as ET
    
    try:
        import gzip
        with gzip.open(filepath, 'rt', encoding='utf-8') as f:
            content = f.read()
        root = ET.fromstring(content)
    except Exception:
        root = ET.parse(filepath).getroot()
    
    spectra = []
    
    # Find all spectrum elements (handles namespaces)
    def find_all_spectra(elem):
        results = []
        for child in elem:
            tag_local = child.tag.split('}')[-1] if '}' in child.tag else child.tag
            if tag_local == 'spectrum':
                results.append(child)
            else:
                results.extend(find_all_spectra(child))
        return results
    
    # Find spectrumList
    def find_spectrum_list(elem):
        for child in elem:
            tag_local = child.tag.split('}')[-1] if '}' in child.tag else child.tag
            if tag_local == 'spectrumList':
                return child
            result = find_spectrum_list(child)
            if result is not None:
                return result
        return None
    
    spectrum_list = find_spectrum_list(root)
    if spectrum_list is None:
        log.warning("No spectrumList found in %s", filepath)
        return spectra
    
    for spec_elem in find_all_spectra(spectrum_list):
        # Check ms level
        ms_level_found = 1
        for elem in spec_elem.iter():
            tag_local = elem.tag.split('}')[-1] if '}' in elem.tag else elem.tag
            if tag_local == 'cvParam':
                acc = elem.get('accession', '')
                if acc == 'MS:1000511':
                    ms_level_found = int(elem.get('value', 1))
                    break
        
        if ms_level_found != ms_level:
            continue
        
        # Find RT from scanList
        rt = None
        for scan_list in spec_elem.iter():
            tag_local = scan_list.tag.split('}')[-1] if '}' in scan_list.tag else scan_list.tag
            if tag_local == 'scanList':
                for scan in scan_list:
                    scan_tag = scan.tag.split('}')[-1] if '}' in scan.tag else scan.tag
                    if scan_tag == 'scan':
                        for cv in scan:
                            cv_tag = cv.tag.split('}')[-1] if '}' in cv.tag else cv.tag
                            if cv_tag == 'cvParam' and cv.get('accession') == 'MS:1000016':
                                rt = float(cv.get('value', 0))
                                unit_acc = cv.get('unitAccession', '')
                                if unit_acc == 'UO:0000031':
                                    pass  # already minutes
                                elif unit_acc == 'UO:0000030':
                                    rt = rt / 60.0 if rt else None
                                break
                    if rt is not None:
                        break
            if rt is not None:
                break
        
        if rt is None:
            continue
        
        # Find binary data arrays
        arrays = {}
        for bda_list in spec_elem.iter():
            tag_local = bda_list.tag.split('}')[-1] if '}' in bda_list.tag else bda_list.tag
            if tag_local == 'binaryDataArrayList':
                for bda in bda_list:
                    bda_tag = bda.tag.split('}')[-1] if '}' in bda.tag else bda.tag
                    if bda_tag != 'binaryDataArray':
                        continue
                    
                    array_type = None
                    data_type_bits = 64
                    binary_text = None
                    
                    for child in bda:
                        child_tag = child.tag.split('}')[-1] if '}' in child.tag else child.tag
                        if child_tag == 'cvParam':
                            acc = child.get('accession', '')
                            if acc in ('MS:1000518', 'MS:1000523', 'MS:1000038'):
                                array_type = 'mz'
                            elif acc in ('MS:1000519', 'MS:1000521', 'MS:1000042'):
                                array_type = 'int'
                            elif acc == 'MS:1000516':
                                data_type_bits = 64
                            elif acc == 'MS:1000515':
                                data_type_bits = 32
                        elif child_tag == 'binary':
                            binary_text = child.text.strip() if child.text else ''
                    
                    if array_type and binary_text:
                        binary_data = base64.b64decode(binary_text)
                        if data_type_bits == 64:
                            fmt = f'<{len(binary_data)//8}d'
                        else:
                            fmt = f'<{len(binary_data)//4}f'
                        arrays[array_type] = list(struct.unpack(fmt, binary_data))
        
        if 'mz' in arrays and 'int' in arrays and len(arrays['mz']) > 0:
            spectra.append((rt, arrays['mz'], arrays['int']))
    
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
