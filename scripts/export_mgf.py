#!/usr/bin/env python3
"""
Export LC-MS Feature Table to MGF Format
==========================================
Converts a CSV feature table (from peak_detection.py) to MGF format
for GNPS FBMN (Feature-based Molecular Networking) upload.

Usage:
    python export_mgf.py --input data/features.csv --output data/features.mgf
"""

import argparse
import sys
from pathlib import Path

import pandas as pd


def export_mgf(df, output_path, ms_level=1):
    """Export feature table to MGF format."""
    required_cols = ["feature_id", "mz_med", "intensity_max", "rt_med"]
    missing = [c for c in required_cols if c not in df.columns]
    if missing:
        sys.stderr.write(
            f"ERROR: Missing required columns: {missing}\n"
            f"Available: {list(df.columns)}\n"
        )
        sys.exit(1)
    
    count = 0
    with open(output_path, "w") as out:
        for _, row in df.iterrows():
            feature_id = int(row["feature_id"])
            mz = float(row["mz_med"])
            intensity = float(row["intensity_max"])
            rt_seconds = float(row["rt_med"]) * 60.0
            
            # Write MGF entry
            out.write("BEGIN IONS\n")
            out.write(f"FEATURE_ID={feature_id}\n")
            out.write(f"RTINSECONDS={rt_seconds:.1f}\n")
            out.write(f"PEPMASS={mz:.5f} {intensity:.1f}\n")
            
            # Add placeholder spectrum peaks (if raw spectra not available)
            # In practice, you would include actual fragment peaks here
            # For feature-based networking, we include the monoisotopic m/z
            out.write(f"{mz:.5f} {intensity:.1f}\n")
            
            out.write("END IONS\n\n")
            count += 1
    
    print(f"Exported {count} features to {output_path}")
    return count


def main():
    parser = argparse.ArgumentParser(
        description="Export LC-MS feature table to MGF format for GNPS"
    )
    parser.add_argument("--input", "-i", required=True,
                        help="Input CSV feature table")
    parser.add_argument("--output", "-o", required=True,
                        help="Output MGF file")
    parser.add_argument("--ms_level", type=int, default=1,
                        help="MS level (default: 1)")
    
    args = parser.parse_args()
    
    df = pd.read_csv(args.input)
    print(f"Loaded {len(df)} features from {args.input}")
    
    count = export_mgf(df, args.output, args.ms_level)
    print(f"\nMGF export complete! {count} features written.")


if __name__ == "__main__":
    main()
