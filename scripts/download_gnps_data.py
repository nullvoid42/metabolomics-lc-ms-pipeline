#!/usr/bin/env python3
"""
Download Public LC-MS Data from GNPS/MassIVE
============================================
Uses the GNPS downloadpublicdata tool or direct USI fetching.

Usage:
    # Download a specific GNPS dataset
    python download_gnps_data.py --dataset MSV000100610 --output data/raw/
    
    # Download via USI (MetaboLights)
    python download_gnps_data.py --usi mzspec:MTBLS2:sample01.mzML --output data/raw/
"""

import argparse
import subprocess
import sys
from pathlib import Path


GNPS_DATASETS = {
    "MSV000100610": {
        "name": "GNPS LC-MS Benchmark Dataset",
        "description": "2000 mzML files, diverse chemical space",
        "url": "https://massive.ucsd.edu/ProteoSAFe/dataset_summary.jsp?task=b3f0b972ca8e4e08a5f1070eed152b74",
    },
    "MSV000084596": {
        "name": "Streptomyces Metabolomics",
        "description": "Streptomyces spp. LC-MS/MS",
        "url": "https://massive.ucsd.edu/ProteoSAFe/dataset.jsp?task=0b7c5c81a7e94bf4b1b88e0c0e1b8f84",
    },
}


def download_gnps_dataset(dataset_id, output_dir, use_tool=True):
    """Download a GNPS dataset using the GNPS public data downloader."""
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    # Clone GNPS download tool if not present
    tool_dir = output_path.parent / "downloadpublicdata"
    if not tool_dir.exists():
        print(f"Cloning GNPS downloadpublicdata tool...")
        try:
            subprocess.run(
                ["git", "clone", "https://github.com/Wang-Bioinformatics-Lab/downloadpublicdata.git",
                 str(tool_dir)],
                check=True,
                capture_output=True,
            )
        except Exception as e:
            print(f"Git clone failed: {e}")
            return False
    else:
        print(f"Using existing download tool at: {tool_dir}")
    
    # Create TSV file for the dataset
    tsv_file = tool_dir / "data" / f"{dataset_id}.tsv"
    tsv_file.parent.mkdir(exist_ok=True)
    with open(tsv_file, "w") as f:
        f.write(f"dataset\t{dataset_id}\n")
    
    # Run downloader
    print(f"Downloading {dataset_id} to {output_path}...")
    print("This may take a while for large datasets...")
    
    try:
        result = subprocess.run(
            [sys.executable, str(tool_dir / "bin" / "download_public_data_usi.py"),
             str(tsv_file), str(output_path), str(tool_dir / "data" / f"{dataset_id}_summary.tsv"),
             "--noconversion"],
            check=True,
            capture_output=True,
            text=True,
            timeout=3600,
        )
        print(result.stdout)
        return True
    except subprocess.TimeoutExpired:
        print("Download timed out (>1 hour). Try downloading a subset.")
        return False
    except Exception as e:
        print(f"Download failed: {e}")
        return False


def download_via_ftp(study_id, output_dir):
    """Download MTBLS study via EBI FTP."""
    import ftplib
    
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    ftp_url = "ftp.ebi.ac.uk"
    remote_dir = f"/pub/databases/metabolights/studies/public/{study_id}/mzML"
    
    print(f"Connecting to {ftp_url}...")
    try:
        with ftplib.FTP(ftp_url) as ftp:
            ftp.login()
            ftp.cwd(remote_dir)
            files = ftp.nlst()
            
            print(f"Found {len(files)} files in {remote_dir}")
            count = 0
            for fname in files:
                if not fname.endswith(".mzML"):
                    continue
                local_path = output_path / fname
                if local_path.exists():
                    continue
                with open(local_path, "wb") as f:
                    ftp.retrbinary(f"RETR {fname}", f.write)
                count += 1
                if count % 10 == 0:
                    print(f"  Downloaded {count}/{len([f for f in files if f.endswith('.mzML')])} files...")
            
            print(f"Download complete: {count} files saved to {output_path}")
            return True
    except Exception as e:
        print(f"FTP download failed: {e}")
        return False


def main():
    parser = argparse.ArgumentParser(
        description="Download public LC-MS data from GNPS/MassIVE/MetaboLights"
    )
    parser.add_argument("--dataset", default=None,
                        help="GNPS/MassIVE dataset ID (e.g., MSV000100610)")
    parser.add_argument("--usi", default=None,
                        help="mzSpec USI (e.g., mzspec:MTBLS2:sample.mzML)")
    parser.add_argument("--output", "-o", required=True,
                        help="Output directory")
    parser.add_argument("--list_datasets", action="store_true",
                        help="List available GNPS datasets")
    parser.add_argument("--mtbls", default=None,
                        help="MetaboLights study ID (e.g., MTBLS2)")
    
    args = parser.parse_args()
    
    if args.list_datasets:
        print("Available GNPS Datasets:")
        for did, info in GNPS_DATASETS.items():
            print(f"  {did}: {info['name']}")
            print(f"    {info['description']}")
        return
    
    if args.mtbls:
        print(f"Downloading MetaboLights study {args.mtbls} via FTP...")
        success = download_via_ftp(args.mtbls, args.output)
        if success:
            print("Done!")
        return
    
    if args.dataset:
        success = download_gnps_dataset(args.dataset, args.output)
        if success:
            print("Dataset download complete!")
        else:
            print("Dataset download failed. Try --list_datasets for alternatives.")
        return
    
    print("Please specify --dataset, --mtbls, or --usi. Use --list_datasets to see options.")
    print("\nExample usage:")
    print("  python download_gnps_data.py --dataset MSV000100610 --output data/raw/")
    print("  python download_gnps_data.py --mtbls MTBLS2 --output data/raw/")


if __name__ == "__main__":
    main()
