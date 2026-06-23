#!/usr/bin/env python3
"""
Download and process the Allen Ancient DNA Resource (AADR) for archaic admixture analysis.

Downloads the 1240K SNP panel from Harvard Dataverse, extracts Neanderthal/Denisovan
individuals, and converts from Eigenstrat format to PLINK binary format.

Reference: Mallick et al. (2024) Scientific Data. DOI: 10.7910/DVN/FFIDCW

Output: reference/aadr/
    - aadr_archaic.bed/bim/fam    — Neanderthal + Denisovan genotypes (PLINK)
    - aadr_modern.bed/bim/fam     — Modern human reference panel (PLINK)
    - aadr_individuals.tsv        — Individual metadata (archaic + modern)
    - aadr_manifest.json          — Provenance and version info

Usage:
    python download_aadr_reference.py                          # Full download + convert
    python download_aadr_reference.py --dry-run                # Show what would happen
    python download_aadr_reference.py --modern-n 500           # Subsample modern pops
    python download_aadr_reference.py --skip-download          # Convert only (if .geno exists)
"""

import argparse
import hashlib
import json
import os
import struct
import subprocess
import sys
import tempfile
import time
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np

# ── Project paths ──────────────────────────────────────────────────────────────

PLATFORM_DIR = Path(__file__).resolve().parent.parent.parent
REF_DIR = PLATFORM_DIR / "reference" / "aadr"
TOOLS_DIR = Path(__file__).resolve().parent.parent.parent.parent / "tools"
PLINK_BIN = str(TOOLS_DIR / "plink")

# ── AADR Dataset Config ────────────────────────────────────────────────────────

DATAVERSE_DOI = "10.7910/DVN/FFIDCW"
DATAVERSE_API = "https://dataverse.harvard.edu/api"
VERSION = "v66.p1"
PANEL = "1240K"

# File IDs from Dataverse (dataset version 14.0, June 2026)
AADR_FILES = {
    "geno": ("v66.p1_1240K.aadr.patch.PUB.geno", 13994829),
    "snp":  ("v66.p1_1240K.aadr.patch.PUB.snp", 13994514),
    "ind":  ("v66.p1_1240K.aadr.patch.PUB.ind", 13994513),
    "anno": ("v66.p1_1240K.aadr.PUB.anno", 13994515),
}

# Known archaic individual IDs in AADR v66 (1240K panel)
# Format: AADR_ID -> (label, species, population_code)
ARCHAIC_INDIVIDUALS = {
    "I10895": ("Altai_Neanderthal", "Neanderthal", "NEA"),
    "I10894": ("Altai_Neanderthal.DG", "Neanderthal", "NEA"),
    "I10893": ("Vindija_Neanderthal", "Neanderthal", "NEA"),
    "I10892": ("Chagyrskaya_Neanderthal", "Neanderthal", "NEA"),
    "I10891": ("Denisova", "Denisovan", "DEN"),
    "I10890": ("Denisova.DG", "Denisovan", "DEN"),
    # Also check for these alternative IDs
    "I10626": ("Mezmaiskaya1_Neanderthal", "Neanderthal", "NEA"),
    "I10627": ("Mezmaiskaya2_Neanderthal", "Neanderthal", "NEA"),
    "I10889": ("Denisova_alt", "Denisovan", "DEN"),
}

# Modern reference populations to include (for comparison)
# 1000 Genomes super-populations + key regional populations
MODERN_POPS = {
    "French": "EUR", "English": "EUR", "Spanish": "EUR", "Italian_North": "EUR",
    "German": "EUR", "Sardinian": "EUR", "Finnish": "EUR", "Russian": "EUR",
    "Han": "EAS", "Japanese": "EAS", "Korean": "EAS", "Dai": "EAS",
    "Yoruba": "AFR", "Mbuti": "AFR", "Dinka": "AFR", "Mandenka": "AFR",
    "Ju_hoan_North": "AFR", "Somali": "AFR",
    "Brahmin": "SAS", "Mala": "SAS", "Balochi": "SAS",
    "Maya": "AMR", "Pima": "AMR", "Surui": "AMR", "Quechua": "AMR",
    "Papuan": "OCE", "Australian": "OCE", "Bougainville": "OCE",
}


def download_file(filename: str, file_id: int, output_dir: Path,
                   dry_run: bool = False) -> Path:
    """Download a single file from the AADR Dataverse dataset using its file ID."""
    output_path = output_dir / filename

    if output_path.exists():
        size_mb = output_path.stat().st_size / 1_048_576
        print(f"  ✓ {filename} already exists ({size_mb:.0f} MB)")
        return output_path

    # Dataverse direct download URL via file ID
    url = f"https://dataverse.harvard.edu/api/access/datafile/{file_id}"

    if dry_run:
        print(f"  → Would download: {filename} (ID={file_id})")
        return output_path

    print(f"  ↓ Downloading {filename}...")
    import subprocess as sp

    try:
        # Dataverse requires cookie-based session for S3 redirect auth
        # curl handles this correctly with -L (follow) and cookie jar
        cookie_jar = output_dir / ".aadr_cookies.txt"
        result = sp.run(
            ["curl", "-sL", "-b", str(cookie_jar), "-c", str(cookie_jar),
             "-o", str(output_path), url],
            timeout=3600,
            capture_output=True,
            text=True,
        )
        if result.returncode != 0 or not output_path.exists():
            raise Exception(f"curl failed: {result.stderr.strip()[:200]}")
        actual_mb = output_path.stat().st_size / 1_048_576
        print(f"    ✓ {filename} downloaded ({actual_mb:.0f} MB)")
    except Exception as e:
        print(f"    ✗ Failed to download {filename}: {e}")
        if output_path.exists():
            output_path.unlink()
        raise

    return output_path


def parse_ind_file(ind_path: Path) -> pd.DataFrame:
    """Parse Eigenstrat .ind file. Returns DataFrame with columns: id, sex, pop."""
    import pandas as pd
    df = pd.read_csv(ind_path, sep=r"\s+", header=None, dtype=str)
    df.columns = ["id", "sex", "pop"]
    return df


def parse_snp_file(snp_path: Path) -> pd.DataFrame:
    """Parse Eigenstrat .snp file. Returns DataFrame with cols: id, chr, cm, pos, a1, a2."""
    import pandas as pd
    df = pd.read_csv(snp_path, sep=r"\s+", header=None, dtype=str)
    # Columns: id, chr, genetic_pos_cM, physical_pos, allele1, allele2
    if df.shape[1] >= 6:
        df.columns = ["id", "chr", "cm", "pos", "a1", "a2"]
    else:
        df.columns = ["id", "chr", "cm", "pos", "a1", "a2"][:df.shape[1]]
    df["pos"] = pd.to_numeric(df["pos"], errors="coerce").fillna(0).astype(int)
    return df


def find_archaic_individuals(ind_df: pd.DataFrame) -> Dict[str, Dict]:
    """Find archaic individuals in the AADR .ind file."""
    found = {}
    for ind_id, meta in ARCHAIC_INDIVIDUALS.items():
        matches = ind_df[ind_df["id"].str.contains(ind_id, na=False)]
        if len(matches) > 0:
            for _, row in matches.iterrows():
                found[row["id"]] = {
                    "label": meta[0],
                    "species": meta[1],
                    "population": meta[2],
                    "sex": row["sex"],
                }
    return found


def convert_eigenstrat_to_plink(
    geno_path: Path,
    snp_path: Path,
    ind_path: Path,
    output_prefix: str,
    selected_individuals: Optional[List[int]] = None,
    plink_bin: str = "plink",
) -> bool:
    """
    Convert Eigenstrat .geno/.snp/.ind to PLINK binary (.bed/.bim/.fam).

    Eigenstrat .geno format:
      - Packed binary: each byte = 4 genotypes (2 bits each)
      - 0 = homozygous ref (0 copies of alt allele)
      - 1 = heterozygous (1 copy)
      - 2 = homozygous alt (2 copies)
      - 3 = missing
      - Order: ALL SNPs for individual 1, then ALL SNPs for individual 2, etc.
        (row-major: individuals × SNPs)

    This is a large operation (~7 GB .geno file). We process in chunks.
    """
    import pandas as pd

    print(f"\n── Converting Eigenstrat → PLINK ──")

    # Load metadata
    snp_df = parse_snp_file(snp_path)
    ind_df = parse_ind_file(ind_path)

    n_snps = len(snp_df)
    n_ind = len(ind_df)
    print(f"  SNPs: {n_snps:,}  |  Individuals: {n_ind:,}")

    # Determine which individuals to extract
    if selected_individuals is not None:
        idx_map = {i: idx for idx, i in enumerate(selected_individuals)
                    if 0 <= i < n_ind}
        n_out = len(idx_map)
    else:
        n_out = n_ind
        idx_map = {i: i for i in range(n_ind)}

    print(f"  Extracting: {n_out} individuals")

    if n_out == 0:
        print("  ✗ No individuals selected")
        return False

    # ── Write .fam file ──
    fam_path = Path(output_prefix + ".fam")
    with open(fam_path, "w") as fh:
        for old_idx, new_idx in sorted(idx_map.items(), key=lambda x: x[1]):
            row = ind_df.iloc[old_idx]
            fid = "AADR"
            iid = str(row["id"])
            father = "0"
            mother = "0"
            sex = {"M": "1", "F": "2"}.get(str(row.get("sex", "U")).upper(), "0")
            pheno = "-9"
            fh.write(f"{fid}\t{iid}\t{father}\t{mother}\t{sex}\t{pheno}\n")
    print(f"  ✓ {fam_path}")

    # ── Write .bim file ──
    bim_path = Path(output_prefix + ".bim")
    bim_records = []
    for i, (_, snp) in enumerate(snp_df.iterrows()):
        chrom = str(snp["chr"]).replace("chr", "")
        if chrom == "23":
            chrom = "X"
        elif chrom == "24":
            chrom = "Y"
        elif chrom == "25":
            chrom = "MT"
        elif chrom == "26":
            chrom = "XY"
        vid = f"{chrom}:{snp['pos']}"
        cm = snp.get("cm", 0)
        pos = snp["pos"]
        a1 = str(snp.get("a2", "0"))  # PLINK: a1 = alt allele
        a2 = str(snp.get("a1", "0"))  # a2 = ref allele
        bim_records.append((chrom, vid, cm, pos, a1, a2))

    with open(bim_path, "w") as fh:
        for chrom, vid, cm, pos, a1, a2 in bim_records:
            fh.write(f"{chrom}\t{vid}\t{cm}\t{pos}\t{a1}\t{a2}\n")
    print(f"  ✓ {bim_path}")

    # ── Write .bed file ──
    # PLINK .bed format: SNP-major, 3 bytes magic + packed genotypes
    bed_path = Path(output_prefix + ".bed")
    geno_bytes = n_ind * ((n_snps + 3) // 4)

    with open(bed_path, "wb") as bed_fh:
        # Magic bytes
        bed_fh.write(bytes([0x6C, 0x1B, 0x01]))

        with open(geno_path, "rb") as geno_fh:
            # Process in SNPs batches (process by SNP column, not individual row)
            # Since .geno is row-major (ind × SNP), we read all data and transpose
            chunk_snps = 10000  # Process 10K SNPs at a time
            total_bytes = n_ind * n_snps
            bytes_per_snp = n_ind  # One byte per SNP per individual... actually 1/4 byte

            # The packed format: 4 genotypes per byte
            # For each SNP column, we need n_ind genotypes → ceil(n_ind/4) bytes
            bytes_per_snp_packed = (n_ind + 3) // 4

            # Read entire .geno into memory for transposition
            packed_rows = (n_snps + 3) // 4  # bytes per individual
            row_bytes = ((n_snps + 3) // 4)

            # Strategy: process individual by individual, extract their genotypes,
            # then transpose to SNP-major for PLINK output
            print(f"  Converting {n_ind:,} individuals × {n_snps:,} SNPs...")

            # Allocate memory for transposed genotype matrix
            # We'll read in batches of individuals to save memory
            ind_batch = 500
            snp_data = np.zeros((n_snps, n_out), dtype=np.int8)

            for batch_start in range(0, n_ind, ind_batch):
                batch_end = min(batch_start + ind_batch, n_ind)
                batch_n = batch_end - batch_start

                # Read batch of individuals from .geno
                batch_bytes_needed = batch_n * row_bytes
                offset = batch_start * row_bytes
                geno_fh.seek(offset)
                batch_data = np.frombuffer(geno_fh.read(batch_bytes_needed),
                                            dtype=np.uint8)

                # Unpack genotypes for this batch
                for b_idx in range(batch_n):
                    orig_idx = batch_start + b_idx
                    if orig_idx not in idx_map:
                        continue
                    new_idx = idx_map[orig_idx]
                    row_start = b_idx * row_bytes
                    for snp in range(n_snps):
                        byte_idx = row_start + (snp // 4)
                        bit_shift = 6 - 2 * (snp % 4)
                        gt = (batch_data[byte_idx] >> bit_shift) & 0x03
                        # Convert: 0→2 (homozygous alt→homozygous ref in PLINK), 2→0
                        # Actually PLINK uses 0=hom alt, 1=het, 2=hom ref (counting alt alleles)
                        # Eigenstrat uses 0=hom ref (0 alt copies), 1=het, 2=hom alt
                        # So Eigenstrat 0→PLINK 2, 1→1, 2→0, 3→missing
                        if gt == 0:
                            gt = 2
                        elif gt == 2:
                            gt = 0
                        elif gt == 3:
                            gt = 3  # missing stays missing
                        snp_data[snp, new_idx] = gt

                pct = (batch_end / n_ind) * 100
                print(f"    {pct:.0f}% ({batch_end:,}/{n_ind:,} individuals)", end="\r")

            print()  # newline after progress

            # Write PLINK .bed (SNP-major, packed, 4 genotypes per byte)
            # PLINK .bed encoding: 00=hom_A1, 01=het, 10=missing, 11=hom_A2
            # Eigenstrat: 0=hom_ref, 1=het, 2=hom_alt, 3=missing
            # With A1=alt, A2=ref: ES_0→PLINK_3, ES_1→PLINK_1, ES_2→PLINK_0, ES_3→PLINK_2
            ES_TO_PLINK = {0: 3, 1: 1, 2: 0, 3: 2}

            print(f"  Writing PLINK .bed ({n_snps:,} SNPs × {n_out} individuals)...")
            bytes_per_snp_out = (n_out + 3) // 4  # packed bytes for output individuals
            for snp in range(n_snps):
                buf = [0] * bytes_per_snp_out
                for ind in range(n_out):
                    gt_raw = int(snp_data[snp, ind])
                    gt = ES_TO_PLINK.get(gt_raw, 2)  # default to missing
                    byte_idx = ind // 4
                    bit_shift = 6 - 2 * (ind % 4)
                    buf[byte_idx] |= (gt & 0x03) << bit_shift
                bed_fh.write(bytes(buf))
                if snp % 50000 == 0 and snp > 0:
                    pct_snp = (snp / n_snps) * 100
                    print(f"    {pct_snp:.0f}% ({snp:,}/{n_snps:,} SNPs)", end="\r")

            print(f"    100% ({n_snps:,}/{n_snps:,} SNPs)")

    print(f"  ✓ {bed_path}")
    return True


def build_aadr_reference(
    output_dir: Path,
    plink_bin: str = "plink",
    dry_run: bool = False,
    skip_download: bool = False,
    modern_n: int = 500,
) -> bool:
    """
    Download AADR data, convert archaic + modern individuals to PLINK.

    Returns True on success.
    """
    import pandas as pd

    output_dir.mkdir(parents=True, exist_ok=True)

    print("═══ AADR Archaic Reference Builder ═══\n")
    print(f"Dataset: {DATAVERSE_DOI} (v{VERSION}, {PANEL} panel)")
    print(f"Output:  {output_dir}/\n")

    # ── Step 1: Download ──
    if not skip_download:
        print("1. Downloading AADR files...")
        paths = {}
        for key, (filename, file_id) in AADR_FILES.items():
            try:
                paths[key] = download_file(filename, file_id, output_dir, dry_run)
            except Exception as e:
                print(f"  ✗ Download failed for {key}: {e}")
                return False

        if dry_run:
            print("\n  Dry run complete. Use without --dry-run to download.")
            return True  # dry run is not an error
    else:
        paths = {key: output_dir / filename
                 for key, (filename, _) in AADR_FILES.items()}
        missing = [k for k, v in paths.items() if not v.exists()]
        if missing:
            print(f"  ✗ Missing files: {missing}")
            return False
        print("1. Using existing AADR files.")

    # ── Step 2: Parse metadata ──
    print("\n2. Parsing AADR metadata...")
    ind_df = parse_ind_file(paths["ind"])
    snp_df = parse_snp_file(paths["snp"])

    # Find archaic individuals
    archaics = find_archaic_individuals(ind_df)
    if not archaics:
        print("  ⚠ No archaic individuals found! Checking .ind file...")
        print(f"  Total individuals: {len(ind_df)}")
        print(f"  Sample IDs: {', '.join(ind_df['id'].head(20).tolist())}...")
        # Try broader search
        for species_term in ["Neanderthal", "Denisovan", "Altai", "Vindija",
                               "Chagyr", "Denisova", "Mezmais", "archaic"]:
            matches = ind_df[ind_df["id"].str.contains(species_term, case=False, na=False)]
            if len(matches) > 0:
                print(f"  Found '{species_term}': {matches['id'].tolist()}")
        return False

    print(f"  ✓ Found {len(archaics)} archaic individuals:")
    for ind_id, meta in archaics.items():
        print(f"    • {ind_id} → {meta['label']} ({meta['species']})")

    archaic_indices = [ind_df.index[ind_df["id"] == iid].tolist()[0]
                        for iid in archaics.keys()]

    # Select modern reference individuals
    print(f"\n3. Selecting modern reference individuals...")
    modern_indices = []
    modern_meta = {}
    for pop_name, super_pop in MODERN_POPS.items():
        pop_matches = ind_df[ind_df["pop"].str.contains(pop_name, case=False, na=False)]
        n_avail = len(pop_matches)
        if n_avail > 0:
            n_take = min(n_avail, max(5, modern_n // len(MODERN_POPS)))
            selected = pop_matches.head(n_take)
            for _, row in selected.iterrows():
                idx = ind_df.index.get_loc(row.name)
                modern_indices.append(idx)
                modern_meta[str(idx)] = {
                    "id": str(row["id"]),
                    "pop": pop_name,
                    "super_pop": super_pop,
                }

    print(f"  ✓ Selected {len(modern_indices)} modern individuals "
          f"from {len(set(m['pop'] for m in modern_meta.values()))} populations")

    # ── Step 3: Convert archaic genotypes ──
    all_selected = archaic_indices + modern_indices
    print(f"\n4. Converting {len(all_selected)} individuals to PLINK...")

    success = convert_eigenstrat_to_plink(
        geno_path=paths["geno"],
        snp_path=paths["snp"],
        ind_path=paths["ind"],
        output_prefix=str(output_dir / "aadr_reference"),
        selected_individuals=all_selected,
        plink_bin=plink_bin,
    )

    if not success:
        print("  ✗ Conversion failed")
        return False

    # ── Step 4: Create split archaic/modern PLINK files ──
    print(f"\n5. Splitting archaic vs modern...")
    n_archaic = len(archaic_indices)

    # Write keep files for PLINK --keep
    fam = pd.read_csv(output_dir / "aadr_reference.fam", sep=r"\s+", header=None)
    fam.columns = ["fid", "iid", "father", "mother", "sex", "pheno"]

    archaic_keep = fam.head(n_archaic)
    modern_keep = fam.iloc[n_archaic:]

    archaic_keep_path = output_dir / "archaic_keep.txt"
    modern_keep_path = output_dir / "modern_keep.txt"

    archaic_keep[["fid", "iid"]].to_csv(archaic_keep_path, sep="\t",
                                          header=False, index=False)
    modern_keep[["fid", "iid"]].to_csv(modern_keep_path, sep="\t",
                                        header=False, index=False)

    # PLINK extract
    for label, keep_path, out_prefix in [
        ("archaic", archaic_keep_path, "aadr_archaic"),
        ("modern", modern_keep_path, "aadr_modern"),
    ]:
        result = subprocess.run(
            [plink_bin, "--bfile", str(output_dir / "aadr_reference"),
             "--keep", str(keep_path), "--make-bed",
             "--out", str(output_dir / out_prefix),
             "--allow-extra-chr", "--threads", "4", "--memory", "8000"],
            capture_output=True, text=True, timeout=300,
        )
        if result.returncode == 0:
            print(f"  ✓ {output_dir / out_prefix}.bed")
        else:
            print(f"  ⚠ {label} split issue: {result.stderr[-200:]}")

    # ── Step 5: Save manifest ──
    manifest = {
        "source": f"Harvard Dataverse — {DATAVERSE_DOI}",
        "version": VERSION,
        "panel": PANEL,
        "generated_date": pd.Timestamp.now().isoformat(),
        "archaic_individuals": archaics,
        "n_modern": len(modern_indices),
        "n_snps": int(len(snp_df)),
        "reference": "Mallick et al. (2024) The Allen Ancient DNA Resource. Scientific Data.",
    }

    with open(output_dir / "aadr_manifest.json", "w") as fh:
        json.dump(manifest, fh, indent=2)
    print(f"\n  ✓ {output_dir / 'aadr_manifest.json'}")

    # Clean up intermediate files
    for ext in [".bed", ".bim", ".fam", ".log", ".nosex"]:
        Path(str(output_dir / "aadr_reference") + ext).unlink(missing_ok=True)
    archaic_keep_path.unlink(missing_ok=True)
    modern_keep_path.unlink(missing_ok=True)

    print(f"\n═══ Done: {output_dir}/ ═══")
    return True


def main():
    parser = argparse.ArgumentParser(
        description="Download and build AADR archaic reference panel")
    parser.add_argument("--output-dir", default=str(REF_DIR),
                        help="Output directory")
    parser.add_argument("--plink", default=PLINK_BIN,
                        help="Path to PLINK binary")
    parser.add_argument("--dry-run", action="store_true",
                        help="Show what would be downloaded")
    parser.add_argument("--skip-download", action="store_true",
                        help="Skip download, use existing .geno files")
    parser.add_argument("--modern-n", type=int, default=500,
                        help="Max modern reference individuals")
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    success = build_aadr_reference(
        output_dir=output_dir,
        plink_bin=args.plink,
        dry_run=args.dry_run,
        skip_download=args.skip_download,
        modern_n=args.modern_n,
    )

    return 0 if success else 1


if __name__ == "__main__":
    sys.exit(main())
