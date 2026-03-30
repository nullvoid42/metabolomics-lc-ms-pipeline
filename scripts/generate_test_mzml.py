#!/usr/bin/env python3
"""
Generate a synthetic LC-MS mzML file for testing the peak detection pipeline.
Creates realistic-looking LC-MS data with Gaussian peaks.
"""

import argparse
import base64
import gzip
import random
import struct
import sys
import math
from pathlib import Path


def gaussian(x, mu, sigma, amplitude):
    """Gaussian peak function."""
    return amplitude * math.exp(-0.5 * ((x - mu) / sigma) ** 2)


def generate_synthetic_mzml(output_path, n_spectra=500, seed=42):
    """Generate a synthetic mzML file with realistic LC-MS peaks."""
    random.seed(seed)
    
    # Realistic m/z values for small molecules (50-600 Th)
    # with some characteristic clusters
    peak_centers = [
        (120.08, 2.5, 15000),   # mz=120, RT=2.5 min, intensity=15000
        (165.11, 5.0, 25000),   # mz=165, RT=5 min
        (212.09, 7.5, 18000),   # mz=212, RT=7.5 min
        (280.14, 10.0, 30000),  # mz=280, RT=10 min
        (340.21, 12.5, 22000),  # mz=340, RT=12.5 min
        (445.32, 15.0, 35000), # mz=445, RT=15 min
        (520.18, 17.5, 20000),  # mz=520, RT=17.5 min
        (580.25, 20.0, 28000),  # mz=580, RT=20 min
        (89.05, 3.2, 8000),     # small peak
        (200.15, 8.8, 12000),
        (315.09, 13.5, 16000),
        (455.22, 18.2, 19000),
    ]
    
    rt_range = (0.5, 25.0)  # minutes
    mz_range = (50.0, 650.0)
    rt_sigma = 0.3  # RT peak width in minutes
    
    spectra_xml = []
    for spec_idx in range(n_spectra):
        rt = rt_range[0] + (rt_range[1] - rt_range[0]) * spec_idx / (n_spectra - 1)
        
        # Generate mz/intensity pairs
        mz_list = []
        int_list = []
        
        # Add noise floor
        noise_level = 200
        for mz_val in range(int(mz_range[0]), int(mz_range[1]), 5):
            mz_float = float(mz_val) + random.uniform(-0.5, 0.5)
            intensity = random.uniform(noise_level * 0.5, noise_level * 1.5)
            mz_list.append(mz_float)
            int_list.append(intensity)
        
        # Add Gaussian peaks
        for mz_center, rt_center, amplitude in peak_centers:
            rt_diff = abs(rt - rt_center)
            if rt_diff < 3 * rt_sigma:
                peak_intensity = gaussian(rt, rt_center, rt_sigma, amplitude)
                # Add some mz isotope peaks
                for delta_mz in [0.0, 1.003, 2.006]:
                    mz_peak = mz_center + delta_mz
                    # Isotope ratio (M+1 is ~30% of M, M+2 is ~5%)
                    if delta_mz == 0.0:
                        iso_factor = 1.0
                    elif delta_mz == 1.003:
                        iso_factor = 0.3
                    else:
                        iso_factor = 0.05
                    
                    # Add to spectrum
                    peak_int = peak_intensity * iso_factor
                    mz_list.append(mz_peak)
                    int_list.append(peak_int + random.uniform(0, noise_level * 0.2))
        
        # Sort by m/z
        sorted_pairs = sorted(zip(mz_list, int_list), key=lambda x: x[0])
        mz_list = [p[0] for p in sorted_pairs]
        int_list = [p[1] for p in sorted_pairs]
        
        # Format as binary arrays (64-bit float, little-endian, base64 encoded)
        mz_binary = base64.b64encode(struct.pack(f'<{len(mz_list)}d', *mz_list)).decode()
        int_binary = base64.b64encode(struct.pack(f'<{len(int_list)}d', *int_list)).decode()
        
        spec_xml = f"""    <spectrum index="{spec_idx}" id="scan={spec_idx+1}" spotID="">
      <cvParam cvRef="MS" accession="MS:1000511" name="ms level" value="1"/>
      <cvParam cvRef="MS" accession="MS:1000128" name="centroid spectrum"/>
      <cvParam cvRef="MS" accession="MS:1000294" name="mass spectrum"/>
      <cvParam cvRef="MS" accession="MS:1000528" name="lowest observed m/z" value="{min(mz_list):.4f}"/>
      <cvParam cvRef="MS" accession="MS:1000527" name="highest observed m/z" value="{max(mz_list):.4f}"/>
      <cvParam cvRef="MS" accession="MS:1000504" name="base peak m/z" value="{mz_list[int_list.index(max(int_list))]:.4f}"/>
      <cvParam cvRef="MS" accession="MS:1000505" name="base peak intensity" value="{max(int_list):.1f}"/>
      <cvParam cvRef="MS" accession="MS:1000285" name="total ion current" value="{sum(int_list):.1f}"/>
      <scanList count="1">
        <scan instrumentModeRef="">
          <cvParam cvRef="MS" accession="MS:1000016" name="scan start time" value="{rt:.4f}" unitCvRef="UO" unitAccession="UO:0000031" unitName="minute"/>
        </scan>
      </scanList>
      <binaryDataArrayList count="2">
        <binaryDataArray>
          <cvParam cvRef="MS" accession="MS:1000518" name="m/z array" unitCvRef="MS" unitAccession="MS:1000040" arrayLength="{len(mz_list)}"/>
          <cvParam cvRef="MS" accession="MS:1000516" name="64-bit float" value=""/>
          <cvParam cvRef="MS" accession="MS:1000520" name="no compression" value=""/>
          <binary>{mz_binary}</binary>
        </binaryDataArray>
        <binaryDataArray>
          <cvParam cvRef="MS" accession="MS:1000519" name="intensity array" unitCvRef="MS" unitAccession="MS:1000131" arrayLength="{len(int_list)}"/>
          <cvParam cvRef="MS" accession="MS:1000516" name="64-bit float" value=""/>
          <cvParam cvRef="MS" accession="MS:1000520" name="no compression" value=""/>
          <binary>{int_binary}</binary>
        </binaryDataArray>
      </binaryDataArrayList>
    </spectrum>"""
        spectra_xml.append(spec_xml)
    
    # Build full mzML
    mzml_content = f"""<?xml version="1.0" encoding="UTF-8"?>
<mzML xmlns="http://psi.hupo.org/ms/mzml" version="1.1.0">
  <cvList count="1">
    <cv id="MS" fullName="Proteomics Standards Initiative Mass Spectrometry" uri="http://psidev.info/ms/"/>
  </cvList>
  <instrumentConfigurationList count="1">
    <instrumentConfiguration id="IC1">
      <cvParam cvRef="MS" accession="MS:1000031" name="instrument model"/>
    </instrumentConfiguration>
  </instrumentConfigurationList>
  <softwareList count="1">
    <software id="software" version="1.0">
      <cvParam cvRef="MS" accession="MS:1000799" name="custom unreleased software tool"/>
    </software>
  </softwareList>
  <run defaultInstrumentConfigurationRef="IC1" id="run1">
    <spectrumList count="{n_spectra}">
{"".join(spectra_xml)}
    </spectrumList>
  </run>
</mzML>
"""
    
    if output_path.suffix == ".gz":
        with gzip.open(output_path, "wt", encoding="utf-8") as f:
            f.write(mzml_content)
    else:
        with open(output_path, "w", encoding="utf-8") as f:
            f.write(mzml_content)
    
    print(f"Generated synthetic mzML: {output_path}")
    print(f"  Spectra: {n_spectra}")
    print(f"  Peak centers: {len(peak_centers)}")
    print(f"  RT range: {rt_range[0]:.1f} - {rt_range[1]:.1f} min")
    print(f"  m/z range: {mz_range[0]:.0f} - {mz_range[1]:.0f} Th")


def main():
    parser = argparse.ArgumentParser(description="Generate synthetic mzML test file")
    parser.add_argument("--output", "-o", required=True, help="Output mzML file")
    parser.add_argument("--n_spectra", type=int, default=500, help="Number of spectra")
    parser.add_argument("--seed", type=int, default=42, help="Random seed")
    
    args = parser.parse_args()
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    generate_synthetic_mzml(output_path, args.n_spectra, args.seed)


if __name__ == "__main__":
    main()
