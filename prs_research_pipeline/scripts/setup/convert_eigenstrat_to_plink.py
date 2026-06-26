#!/usr/bin/env python3
"""
Convert EIGENSTRAT format (.geno/.snp/.ind) to PLINK format (.bed/.bim/.fam).

EIGENSTRAT .geno format: packed binary, 4 samples per byte.
  Bits 7-6: sample 4n,   bits 5-4: sample 4n+1,
  Bits 3-2: sample 4n+2, bits 1-0: sample 4n+3
  Value: 0=homozygous ref, 1=heterozygous, 2=homozygous alt, 3=missing

This script extracts specific populations (e.g., Jew_Ashkenazi + European reference)
and creates a PLINK binary dataset with overlapping SNPs matched to 1000G.

Usage:
  python3 convert_eigenstrat_to_plink.py \
      --geno reference/human_origins/HumanOriginsPublic2068.geno \
      --snp reference/human_origins/HumanOriginsPublic2068.snp \
      --ind reference/human_origins/HumanOriginsPublic2068.ind \
      --populations "Jew_Ashkenazi,French,Spanish,Basque,Sardinian,English,Orcadian" \
      --bim reference/1000G_full/1000G_full.bim \
      --out reference/human_origins/european_aj_subset
"""

import sys
import os
import argparse
import struct
from pathlib import Path

import numpy as np
import pandas as pd

# EIGENSTRAT genotype decoding: bits 0-1 → genotype 0, bits 2-3 → genotype 1, etc.
GENO_MAP = np.array([0, 0, 0, 0,   # 0000→0000
                      1, 0, 0, 0,   # 0001→1000
                      0, 1, 0, 0,   # 0010→0100
                      0, 0, 1, 0,   # 0011→0010
                      # ... etc. Full 256-entry lookup
                      ], dtype=np.uint8)

# Build the full 256-entry genotype lookup table
GENO_LOOKUP = np.zeros(256, dtype=np.uint8)
for byte_val in range(256):
    g0 = (byte_val >> 6) & 0x03
    g1 = (byte_val >> 4) & 0x03
    g2 = (byte_val >> 2) & 0x03
    g3 = byte_val & 0x03
    # Pack: g0 in bits 6-7, g1 in 4-5, g2 in 2-3, g3 in 0-1
    # 3 = missing → store as 0 (PLINK missing = 0 in FAM, 3 in BED)
    packed = ((3 if g0 == 3 else g0) << 6) | \
             ((3 if g1 == 3 else g1) << 4) | \
             ((3 if g2 == 3 else g2) << 2) | \
             (3 if g3 == 3 else g3)
    GENO_LOOKUP[byte_val] = packed


def decode_geno_block(block: bytes, n_samples: int) -> np.ndarray:
    """Decode one SNP's genotype block. Vectorized for speed."""
    raw = np.frombuffer(block, dtype=np.uint8)
    # Expand: each byte → 4 genotypes via bit shifts
    g0 = (raw >> 6) & 0x03
    g1 = (raw >> 4) & 0x03
    g2 = (raw >> 2) & 0x03
    g3 = raw & 0x03
    # Interleave: [g0[0], g1[0], g2[0], g3[0], g0[1], ...]
    genos = np.column_stack([g0, g1, g2, g3]).ravel()
    return genos[:n_samples]


def main():
    parser = argparse.ArgumentParser(description="Convert EIGENSTRAT to PLINK")
    parser.add_argument("--geno", required=True, help="EIGENSTRAT .geno file")
    parser.add_argument("--snp", required=True, help="EIGENSTRAT .snp file")
    parser.add_argument("--ind", required=True, help="EIGENSTRAT .ind file")
    parser.add_argument("--populations", default="Jew_Ashkenazi",
                       help="Comma-separated populations to extract")
    parser.add_argument("--bim", help="1000G .bim file for SNP matching")
    parser.add_argument("--out", default="european_aj_subset")
    parser.add_argument("--max-snps", type=int, default=10000,
                       help="Max SNPs to extract (subset for speed; 10K is enough for PCA)")
    parser.add_argument("--verbose", "-v", action="store_true")
    args = parser.parse_args()

    target_pops = [p.strip() for p in args.populations.split(",")]
    print(f"Target populations: {target_pops}")

    # ── Load .ind (sample metadata) ──
    ind = pd.read_csv(args.ind, sep=r"\s+", header=None)
    ind.columns = ["sample_id", "sex", "population"]
    print(f"Total samples: {len(ind)}")

    # Find sample indices for target populations
    sample_mask = ind["population"].isin(target_pops)
    sample_indices = np.where(sample_mask)[0]
    n_extract = len(sample_indices)
    if n_extract == 0:
        print("❌ No samples found for the specified populations!")
        sys.exit(1)

    print(f"Extracting {n_extract} samples from {len(target_pops)} populations:")
    for pop in target_pops:
        count = sample_mask[ind["population"] == pop].sum()
        if count > 0:
            print(f"  {pop}: {count} samples")

    # Create reverse index: original index → new index
    new_idx = {old: new for new, old in enumerate(sample_indices)}

    # ── Load .snp (SNP metadata) ──
    snp = pd.read_csv(args.snp, sep=r"\s+", header=None)
    snp.columns = ["rsid", "chrom", "genetic_pos", "physical_pos", "allele1", "allele2"]
    n_snps = len(snp)
    print(f"Total SNPs: {n_snps}")

    # Match with 1000G BIM if provided (by chromosome:position)
    snp_to_keep = np.ones(n_snps, dtype=bool)
    if args.bim and os.path.exists(args.bim):
        bim = pd.read_csv(args.bim, sep=r"\s+", header=None)
        bim.columns = ["chrom", "variant_id", "genetic_dist", "pos", "a1", "a2"]
        # Build chr:pos set from BIM
        bim["chrom"] = bim["chrom"].astype(str).str.replace("chr", "")
        bim["chr_pos"] = bim["chrom"] + ":" + bim["pos"].astype(str)
        bim_keys = set(bim["chr_pos"])

        # Build chr:pos from SNP file
        snp["chrom"] = snp["chrom"].astype(str)
        snp["chr_pos"] = snp["chrom"] + ":" + snp["physical_pos"].astype(str)
        snp_to_keep = snp["chr_pos"].isin(bim_keys).values
        print(f"SNPs overlapping with 1000G (by chr:pos): {snp_to_keep.sum()} / {n_snps}")

        # Also try matching by rsID for additional overlap
        bim_rsids = set(bim["variant_id"].astype(str))
        snp_rsids = snp["rsid"].astype(str)
        rsid_match = snp_rsids.isin(bim_rsids) & ~snp_to_keep
        if rsid_match.sum() > 0:
            print(f"  Additional rsID matches: {rsid_match.sum()}")
            snp_to_keep = snp_to_keep | rsid_match
    else:
        print("No BIM file — keeping all SNPs")

    keep_indices = np.where(snp_to_keep)[0]
    n_keep = len(keep_indices)

    # Subset if too many SNPs
    if args.max_snps and n_keep > args.max_snps:
        # Take evenly spaced subset
        step = n_keep // args.max_snps
        keep_indices = keep_indices[::step][:args.max_snps]
        n_keep = len(keep_indices)
        print(f"Subset to {n_keep} SNPs for speed (PCA needs only ~10K)")

    if n_keep == 0:
        print("❌ No overlapping SNPs found!")
        sys.exit(1)
    n_all_samples = len(ind)
    bytes_per_snp = (n_all_samples + 3) // 4
    extract_bytes_per_snp = (n_extract + 3) // 4

    # ── Write .bim ──
    bim_out = snp.iloc[keep_indices].copy()
    bim_out = bim_out[["chrom", "rsid", "genetic_pos", "physical_pos", "allele1", "allele2"]]
    bim_out.columns = ["chrom", "rsid", "genetic_dist", "pos", "a1", "a2"]
    bim_out.to_csv(f"{args.out}.bim", sep="\t", header=False, index=False)

    # ── Write .fam ──
    fam_data = []
    for i in sample_indices:
        row = ind.iloc[i]
        fam_data.append([row["sample_id"], row["sample_id"], "0", "0",
                         "1" if row["sex"] == "M" else "2", "-9"])
    pd.DataFrame(fam_data).to_csv(f"{args.out}.fam", sep="\t", header=False, index=False)

    # ── Write population_labels.txt ──
    labels_out = []
    for i in sample_indices:
        row = ind.iloc[i]
        pop = row["population"]
        # Map Human Origins labels to our standard codes
        pop_map = {
            "Jew_Ashkenazi": "AJ",
            "French": "CEU", "English": "GBR", "Orcadian": "GBR",
            "Spanish": "IBS", "Basque": "IBS", "Sardinian": "TSI",
            "Italian_North": "TSI", "Italian_South": "TSI",
        }
        labels_out.append([row["sample_id"], pop_map.get(pop, pop)])
    pd.DataFrame(labels_out).to_csv(f"{args.out}_population_labels.txt", sep="\t", header=False, index=False)

    # ── Convert .geno to .bed ──
    print(f"Converting genotypes ({n_keep} SNPs × {n_extract} samples)...")

    # .bed header: magic bytes 0x6c 0x1b + 0x01 (SNP-major)
    bed_header = bytes([0x6c, 0x1b, 0x01])

    with open(args.geno, "rb") as f_geno, open(f"{args.out}.bed", "wb") as f_bed:
        f_bed.write(bed_header)

        prev_pos = 0
        for snp_idx in keep_indices:
            # Seek to this SNP's data
            byte_offset = snp_idx * bytes_per_snp
            f_geno.seek(byte_offset)
            block = f_geno.read(bytes_per_snp)

            if len(block) < bytes_per_snp:
                print(f"  ⚠ Short read at SNP {snp_idx}")
                # Pad with zeros (missing)
                block = block + b'\xff' * (bytes_per_snp - len(block))

            # Decode all samples
            all_genos = decode_geno_block(block, n_all_samples)

            # Extract only target samples
            extracted = np.array([all_genos[i] for i in sample_indices], dtype=np.uint8)

            # Re-pack: 4 samples per byte
            packed = np.zeros(extract_bytes_per_snp, dtype=np.uint8)
            for j in range(n_extract):
                byte_j = j // 4
                shift = (3 - (j % 4)) * 2  # bits: 6,4,2,0
                val = extracted[j] & 0x03
                packed[byte_j] |= (val << shift)

            f_bed.write(packed.tobytes())

            if args.verbose and (snp_idx + 1) % 50000 == 0:
                pct = (snp_idx + 1) / n_keep * 100
                print(f"  {snp_idx + 1}/{n_keep} SNPs ({pct:.1f}%)")

    print(f"✅ Converted to {args.out}.bed/.bim/.fam")
    print(f"   SNPs: {n_keep}")
    print(f"   Samples: {n_extract}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
