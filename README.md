# metabolomics-lc-ms-pipeline

**LC-MS Metabolomics Data Analysis Pipeline**  
Automated peak detection, annotation, and visualization for liquid chromatography–mass spectrometry (LC-MS) data.

---

## Overview

This pipeline processes raw LC-MS data (mzML format) through a complete metabolomics workflow:

```
Raw mzML → Peak Detection → Feature Alignment → Annotation → Visualization
```

### Tools & Dependencies

| Tool | Version | Purpose |
|------|---------|---------|
| Python | ≥3.9 | Core scripting |
| pymzml | ≥3.0 | mzML file parsing |
| scipy | ≥1.10 | Peak detection algorithms |
| matplotlib | ≥3.7 | Visualization |
| pandas | ≥2.0 | Tabular data handling |
| MZmine (optional) | ≥3.0 | Full-featured GUI peak detection |

> **Note:** This pipeline includes a native Python mzML parser and peak detection engine
> that runs without any Java dependencies. The MZmine step is entirely optional
> for users who prefer a GUI workflow.

| Tool | Version | Purpose |
|------|---------|---------|
| Python | ≥3.9 | Core scripting |
| pymzML | ≥3.0 | mzML file parsing |
| scipy | ≥1.10 | Peak detection algorithms |
| matplotlib | ≥3.7 | Visualization |
| pandas | ≥2.0 | Tabular data handling |
| MZmine (optional) | ≥3.0 | Full-featured GUI peak detection |

### Reference Datasets

- **MTBLS2**: *Arabidopsis thaliana* LC-MS metabolomics (MetaboLights)
- **GNPS MassIVE MSV000100610**: 2000 mzML LC-MS/MS files

---

## Pipeline Stages

### Stage 1 — Data Acquisition

Download public LC-MS data or place your own `.mzML` files in `data/raw/`.

```bash
# Example: download MTBLS2 sample via FTP
wget -r -np -nH --cut-dirs=4 \
  ftp://ftp.ebi.ac.uk/pub/databases/metabolights/studies/public/MTBLS2/mzML/ \
  -P data/raw/
```

### Stage 2 — Peak Detection

**Python (recommended for automation):**

```bash
python scripts/peak_detection.py \
  --input data/raw/ \
  --output data/peak_list.csv \
  --mzml_pattern "*.mzML" \
  --noise_threshold 500 \
  --min_snr 3.0 \
  --min_peak_width 3
```

**MZmine (GUI alternative):**
1. Open MZmine → `File → Import Raw Data → mzML`
2. Run **Mass detection** (MS1: centroid or vendor pick)
3. Run **Chromatogram builder** → **Chromatogram deconvolution**
4. Run **Feature alignment** (if multiple samples)
5. Export via `Feature list → Export/Submit to GNPS-FBMN`

### Stage 3 — Feature List Export

Export peak tables for downstream annotation:

```bash
# Export to CSV
python scripts/export_csv.py \
  --input data/peak_list.csv \
  --output data/features.csv

# Export to MGF (for GNPS)
python scripts/export_mgf.py \
  --input data/peak_list.csv \
  --output data/features.mgf
```

### Stage 4 — GNPS Annotation

Upload `features.mgf` to [GNPS FBMN](https://gnps.ucsd.edu/ProteoSAFe/index.jsp?params=%7B%22workflow%22:%22FEATURE-BASED-MOLECULAR-NETWORKING%22%7D).

Or use local METLIN database search:

```bash
python scripts/annotate_metlin.py \
  --input data/features.csv \
  --output results/annotated.csv \
  --metlin_csv /path/to/metlin_export.csv
```

### Stage 5 — Visualization

```bash
python scripts/visualize.py \
  --input data/features.csv \
  --output results/figures/
```

Generates:
- Base peak chromatogram (BPC)
- Total ion chromatogram (TIC)
- Peak intensity heatmap
- m/z vs RT scatter plot
- Feature statistics summary

---

## Project Structure

```
metabolomics-lc-ms-pipeline/
├── README.md
├── LICENSE
├── scripts/
│   ├── peak_detection.py      # LC-MS peak detection engine
│   ├── export_csv.py          # CSV exporter
│   ├── export_mgf.py          # MGF mass spectrometry exporter
│   ├── annotate_metlin.py     # METLIN local annotation
│   ├── visualize.py           # Plotting: BPC, TIC, heatmap
│   └── download_gnps_data.py  # GNPS public data downloader
├── data/
│   ├── raw/                   # Raw .mzML files
│   └── peak_list.csv          # Detected peaks (intermediate)
├── results/
│   ├── figures/               # Generated plots
│   └── annotated.csv          # Annotated feature table
└── docs/
    └── pipeline_overview.md   # Detailed pipeline documentation
```

---

## Quick Start

```bash
# 1. Clone this repo
git clone https://github.com/nullvoid42/metabolomics-lc-ms-pipeline.git
cd metabolomics-lc-ms-pipeline

# 2. Install Python dependencies
pip install -r requirements.txt

# 3. Generate synthetic test data (or download real data)
python scripts/generate_test_mzml.py --output data/raw/test_sample.mzML --n_spectra 500

# 4. Run peak detection
python scripts/peak_detection.py --input data/raw/ --output data/peak_list.csv

# 5. Export to MGF for GNPS
python scripts/export_mgf.py --input data/peak_list.csv --output data/features.mgf

# 6. Visualize results
python scripts/visualize.py --input data/peak_list.csv --output results/figures/

# Or run the full pipeline end-to-end:
python scripts/run_pipeline.py --input data/raw/ --output results/
```

### Server Setup (weibin@106.15.234.172)

```bash
# Install conda environment
export MAMBA_ROOT_PREFIX=/home/weibin/mamba
~/bin/micromamba create -n metabolomics python=3.11 pymzml scipy pandas matplotlib seaborn numpy -y

# Activate and run
~/bin/micromamba activate metabolomics
python scripts/peak_detection.py --input data/raw/ --output data/peak_list.csv
```

---

## Local Installation

```bash
# Clone the repo
git clone https://github.com/nullvoid42/metabolomics-lc-ms-pipeline.git
cd metabolomics-lc-ms-pipeline

# Install Python dependencies
pip install -r requirements.txt

# Verify installation
python scripts/peak_detection.py --help
```

### Peak Table (CSV)
```
feature_id,mz_min,mz_max,rt_min,rt_max,intensity,SNR,mzmed,rtmed
1,120.0452,120.0550,245.1,247.3,15234.5,8.2,120.0501,246.2
2,165.1123,165.1201,312.5,314.8,8921.0,5.1,165.1162,313.6
```

### MGF (Mascot Generic Format)
```
COM=mzmine 2 export
FEATURE_ID=1
RTINSECONDS=1477
PEPMASS=120.0501 15234.5
1.0078 12.4
2.0141 8.2
...
END
```

---

## Citation

If you use this pipeline in your research, please cite:

> Zhou Z-Y, et al. (2026). *metabolomics-lc-ms-pipeline: A reproducible LC-MS metabolomics workflow*. GitHub. https://github.com/nullvoid42/metabolomics-lc-ms-pipeline

---

## License

MIT License
