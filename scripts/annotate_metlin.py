#!/usr/bin/env python3
"""
METLIN Database Annotation
==========================
Annotates LC-MS features against a local METLIN export (CSV format).

METLIN CSV expected columns:
  exact_mass, mz, name, formula, smiles, inchi, CAS

Usage:
    python annotate_metlin.py \
        --input data/features.csv \
        --metlin_csv /path/to/metlin_export.csv \
        --output results/annotated.csv \
        --ppm_tolerance 10
"""

import argparse
import sys

import numpy as np
import pandas as pd


def annotate_features(features_df, metlin_df, mz_col="mz_med",
                      ppm_tolerance=10, score_threshold=0.5):
    """Annotate features by matching m/z to METLIN database."""
    required = ["mz", "name"]
    missing = [c for c in required if c not in metlin_df.columns]
    if missing:
        sys.stderr.write(f"ERROR: METLIN CSV missing columns: {missing}\n")
        sys.exit(1)
    
    # Handle different column names for m/z
    mass_col = None
    for col in ["exact_mass", "mz", "mass"]:
        if col in metlin_df.columns:
            mass_col = col
            break
    if mass_col is None:
        sys.stderr.write("ERROR: METLIN CSV must have exact_mass, mz, or mass column\n")
        sys.exit(1)
    
    results = []
    metlin_masses = metlin_df[mass_col].values.astype(float)
    
    for _, row in features_df.iterrows():
        mz_obs = float(row[mz_col])
        mz_tol = mz_obs * ppm_tolerance * 1e-6
        mz_min = mz_obs - mz_tol
        mz_max = mz_obs + mz_tol
        
        matches = metlin_df[
            (metlin_masses >= mz_min) & (metlin_masses <= mz_max)
        ].copy()
        
        if matches.empty:
            results.append({
                **row.to_dict(),
                "metlin_name": None,
                "metlin_formula": None,
                "metlin_mass": None,
                "metlin_cas": None,
                "metlin_smiles": None,
                "mz_error_ppm": None,
                "num_matches": 0,
            })
        else:
            # Take best match (smallest m/z error)
            matches["mz_error_ppm"] = abs(
                matches[mass_col] - mz_obs
            ) / mz_obs * 1e6
            best = matches.nsmallest(1, "mz_error_ppm").iloc[0]
            
            results.append({
                **row.to_dict(),
                "metlin_name": best.get("name", None),
                "metlin_formula": best.get("formula", None),
                "metlin_mass": best.get(mass_col, None),
                "metlin_cas": best.get("CAS", None),
                "metlin_smiles": best.get("smiles", None),
                "mz_error_ppm": round(best["mz_error_ppm"], 3),
                "num_matches": len(matches),
            })
    
    return pd.DataFrame(results)


def main():
    parser = argparse.ArgumentParser(
        description="Annotate LC-MS features against METLIN database"
    )
    parser.add_argument("--input", "-i", required=True,
                        help="Input CSV feature table")
    parser.add_argument("--metlin_csv", required=True,
                        help="METLIN export CSV file")
    parser.add_argument("--output", "-o", required=True,
                        help="Output annotated CSV")
    parser.add_argument("--mz_col", default="mz_med",
                        help="m/z column name in input (default: mz_med)")
    parser.add_argument("--ppm_tolerance", type=float, default=10.0,
                        help="m/z tolerance in PPM (default: 10)")
    parser.add_argument("--score_threshold", type=float, default=0.5,
                        help="Minimum matching score (0-1, default: 0.5)")
    
    args = parser.parse_args()
    
    print(f"Loading features from {args.input}...")
    features = pd.read_csv(args.input)
    print(f"  {len(features)} features loaded")
    
    print(f"Loading METLIN database from {args.metlin_csv}...")
    metlin = pd.read_csv(args.metlin_csv)
    print(f"  {len(metlin)} METLIN entries loaded")
    
    print(f"Annotating with PPM tolerance = {args.ppm_tolerance}...")
    annotated = annotate_features(
        features, metlin,
        mz_col=args.mz_col,
        ppm_tolerance=args.ppm_tolerance,
    )
    
    # Summary
    n_annotated = annotated["metlin_name"].notna().sum()
    print(f"\nAnnotation summary:")
    print(f"  Total features: {len(annotated)}")
    print(f"  Annotated: {n_annotated} ({100*n_annotated/len(annotated):.1f}%)")
    print(f"  Unannotated: {len(annotated) - n_annotated}")
    
    annotated.to_csv(args.output, index=False)
    print(f"\nSaved annotated features to: {args.output}")


if __name__ == "__main__":
    main()
