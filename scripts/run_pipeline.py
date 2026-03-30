#!/usr/bin/env python3
"""
LC-MS Metabolomics Pipeline Runner
==================================
Runs the complete peak detection → annotation → visualization pipeline.

Usage:
    python run_pipeline.py --input data/raw/ --output results/ --metlin metlin.csv
"""

import argparse
import subprocess
import sys
from pathlib import Path


def run_cmd(cmd, description):
    print(f"\n{'='*60}")
    print(f"STEP: {description}")
    print(f"CMD:  {' '.join(cmd)}")
    print('='*60)
    result = subprocess.run(cmd, capture_output=True, text=True)
    print(result.stdout)
    if result.stderr:
        print("STDERR:", result.stderr[:2000])
    if result.returncode != 0:
        print(f"ERROR: Command failed with exit code {result.returncode}")
        return False
    return True


def main():
    parser = argparse.ArgumentParser(
        description="Run the complete LC-MS metabolomics pipeline"
    )
    parser.add_argument("--input", "-i", required=True,
                        help="Directory containing mzML files")
    parser.add_argument("--output", "-o", required=True,
                        help="Output directory")
    parser.add_argument("--metlin", default=None,
                        help="Optional METLIN CSV for annotation")
    parser.add_argument("--min_snr", type=float, default=3.0)
    parser.add_argument("--ppm_tolerance", type=float, default=10.0)
    parser.add_argument("--skip_visualization", action="store_true")
    
    args = parser.parse_args()
    
    input_dir = Path(args.input)
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    data_dir = output_dir / "data"
    data_dir.mkdir(exist_ok=True)
    fig_dir = output_dir / "figures"
    fig_dir.mkdir(exist_ok=True)
    
    scripts_dir = Path(__file__).parent
    
    # Step 1: Peak detection
    peak_csv = data_dir / "peaks_raw.csv"
    ok = run_cmd([
        sys.executable, str(scripts_dir / "peak_detection.py"),
        "--input", str(input_dir),
        "--output", str(peak_csv),
        "--min_snr", str(args.min_snr),
    ], "Peak Detection")
    if not ok:
        sys.exit(1)
    
    # Step 2: Export to MGF
    mgf_file = data_dir / "features.mgf"
    run_cmd([
        sys.executable, str(scripts_dir / "export_mgf.py"),
        "--input", str(peak_csv),
        "--output", str(mgf_file),
    ], "Export to MGF (GNPS format)")
    
    # Step 3: METLIN annotation (if provided)
    if args.metlin:
        annotated_csv = data_dir / "features_annotated.csv"
        run_cmd([
            sys.executable, str(scripts_dir / "annotate_metlin.py"),
            "--input", str(peak_csv),
            "--metlin_csv", str(args.metlin),
            "--output", str(annotated_csv),
            "--ppm_tolerance", str(args.ppm_tolerance),
        ], "METLIN Annotation")
    else:
        annotated_csv = peak_csv
    
    # Step 4: Visualization
    if not args.skip_visualization:
        run_cmd([
            sys.executable, str(scripts_dir / "visualize.py"),
            "--input", str(annotated_csv),
            "--output", str(fig_dir),
        ], "Visualization")
    
    print(f"\n{'='*60}")
    print("PIPELINE COMPLETE")
    print(f"{'='*60}")
    print(f"Peak list:      {peak_csv}")
    print(f"MGF file:       {mgf_file}")
    if args.metlin:
        print(f"Annotated:      {annotated_csv}")
    print(f"Figures:        {fig_dir}/")
    print("\nNext steps:")
    print(f"  1. Upload {mgf_file} to GNPS FBMN: https://gnps.ucsd.edu/ProteoSAFe/")
    print(f"  2. View figures in {fig_dir}/")


if __name__ == "__main__":
    main()
