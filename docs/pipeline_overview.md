# LC-MS Metabolomics Pipeline Overview

## Architecture

```
Raw mzML Files
      │
      ▼
┌─────────────────────────┐
│  Stage 1: Peak Detection │  (peak_detection.py)
│  - pymzml reader        │
│  - TIC/BPC construction │
│  - scipy find_peaks     │
│  - Feature clustering   │
└────────────┬────────────┘
             │
             ▼
┌─────────────────────────┐
│  Stage 2: Export Formats │  (export_mgf.py, export_csv.py)
│  - MGF for GNPS FBMN     │
│  - CSV for downstream    │
└────────────┬────────────┘
             │
       ┌─────┴─────┐
       ▼           ▼
┌──────────┐  ┌──────────────┐
│   GNPS   │  │    METLIN    │
│ FBMN web │  │ Local search │
└──────────┘  └──────────────┘
       │           │
       └─────┬─────┘
             ▼
┌─────────────────────────┐
│  Stage 3: Visualization │  (visualize.py)
│  - m/z vs RT scatter    │
│  - BPC/TIC chromatogram │
│  - Intensity heatmap    │
│  - Feature statistics   │
└─────────────────────────┘
```

## Peak Detection Algorithm

### Step 1: Spectrum Parsing
- Read mzML using `pymzml` (supports zlib/gzip compression)
- Extract MS1 spectra with (rt, mz[], intensity[])

### Step 2: Chromatogram Construction
- **TIC**: Sum of intensities per scan → Total Ion Chromatogram
- **BPC**: Max intensity per scan → Base Peak Chromatogram

### Step 3: Peak Detection in Chromatogram
```python
scipy.signal.find_peaks(
    savgol_filter(tic, window=25, polyorder=3),  # smoothed
    height=noise_threshold * min_snr,
    prominence=noise_threshold * min_snr * 0.5,
    width=min_peak_width,
)
```

### Step 4: Feature Clustering
- Window-based clustering across m/z (PPM tolerance) and RT (minutes)
- Aggregates: min/max/median m/z, RT range, sum/max intensity

## m/z Tolerance for Annotation

| Tolerance | Use Case |
|-----------|----------|
| 5 PPM | High-resolution Q-TOF, Orbitrap |
| 10 PPM | Standard LC-MS |
| 20 PPM | Low-resolution LC-MS |
| 50 PPM | GC-MS |

## GNPS FBMN Workflow

1. Export features as MGF file
2. Go to https://gnps.ucsd.edu/ProteoSAFe/
3. Select workflow: `FEATURE-BASED-MOLECULAR-NETWORKING`
4. Upload `features.mgf`
5. Set parameters:
   - Precursor ion tolerance: 0.02 Da
   - Fragment ion tolerance: 0.02 Da
   - Minimum peak matches: 4
6. Submit and wait (~10-30 min)
7. Download results: network (SIF/GraphML), annotated MS/MS spectra

## Key References

- Wang et al. (2016) Nature Methods 13: 904-906 — GNPS
- Schmid et al. (2023) Nature Biotechnology 41: 841-844 — MZmine 3
- Plus iCalGNN for computational annotation
