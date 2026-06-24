#!/usr/bin/env python3
"""
Download and index the Vindija Neanderthal VCF for direct archaic comparison.

Source: Max Planck Institute for Evolutionary Anthropology
URL: http://cdna.eva.mpg.de/neandertal/Vindija/VCF/
Reference: Prüfer et al. (2017) Science. doi:10.1126/science.aao1887

The Vindija 33.19 specimen is a ~30x coverage Neanderthal genome aligned to hg19.
This provides the most accurate archaic comparison available.

Output: reference/vindija/
    - vindija_neanderthal.vcf.gz      — Vindija VCF (hg19, ~2 GB)
    - vindija_neanderthal.vcf.gz.tbi  — Tabix index
    - vindija_manifest.json           — Provenance

Usage:
    python download_vindija_reference.py                  # Download + index
    python download_vindija_reference.py --dry-run        # Show what would happen
"""

import argparse
import json
import subprocess
import sys
from pathlib import Path

import pandas as pd

# ── Paths ──────────────────────────────────────────────────────────────────────

PLATFORM_DIR = Path(__file__).resolve().parent.parent.parent
REF_DIR = PLATFORM_DIR / "reference" / "vindija"
VCF_BASE_URL = "http://ftp.eva.mpg.de/neandertal/Vindija/VCF/Vindija33.19"
VCF_FILES = [f"chr{c}_mq25_mapab100.vcf.gz" for c in
             [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22]]
# Total: ~50 GB for all autosomes. chr22 alone = 575 MB.
TOOLS_DIR = Path(__file__).resolve().parent.parent.parent.parent / "tools"

# ── African 1000G populations for baseline ─────────────────────────────────────
# From 1000G Phase 3: YRI, LWK, GWD, MSL, ESN, ASW, ACB
AFRICAN_POPS_1000G = ["YRI", "LWK", "GWD", "MSL", "ESN", "ASW", "ACB"]


def download_vindija_vcf(output_dir: Path, chromosomes: list = None,
                          dry_run: bool = False) -> list[Path]:
    """Download Vindija Neanderthal VCF files per chromosome. Returns list of paths."""
    output_dir.mkdir(parents=True, exist_ok=True)

    if chromosomes is None:
        chromosomes = ["chr22"]  # Default: just chr22 (575 MB)

    paths = []
    for chr_name in chromosomes:
        filename = f"{chr_name}_mq25_mapab100.vcf.gz"
        url = f"{VCF_BASE_URL}/{filename}"
        vcf_path = output_dir / filename

        if vcf_path.exists():
            size_mb = vcf_path.stat().st_size / 1_048_576
            print(f"  ✓ {filename} already exists ({size_mb:.0f} MB)")
            paths.append(vcf_path)
            continue

        if dry_run:
            print(f"  → Would download {filename} from {VCF_BASE_URL}/")
            paths.append(vcf_path)
            continue

        print(f"  ↓ Downloading {filename}...")
        result = subprocess.run(
            ["curl", "-L", "-o", str(vcf_path), url],
            timeout=7200, capture_output=True, text=True,
        )

        if result.returncode != 0 or vcf_path.stat().st_size < 1000:
            print(f"  ✗ Download failed: {result.stderr[-200:]}")
            vcf_path.unlink(missing_ok=True)
            continue

        size_mb = vcf_path.stat().st_size / 1_048_576
        print(f"    ✓ Downloaded ({size_mb:.0f} MB)")
        paths.append(vcf_path)

    return paths


def index_vcfs(vcf_paths: list[Path]) -> list[Path]:
    """Create tabix indexes for Vindija VCF files."""
    tbi_paths = []
    for vcf_path in vcf_paths:
        tbi_path = Path(str(vcf_path) + ".tbi")
        if tbi_path.exists():
            print(f"  ✓ Index exists: {vcf_path.name}.tbi")
        else:
            print(f"  ↓ Indexing {vcf_path.name}...")
            result = subprocess.run(
                ["tabix", "-p", "vcf", str(vcf_path)],
                timeout=600, capture_output=True, text=True,
            )
            if result.returncode == 0:
                print(f"    ✓ Indexed")
            else:
                print(f"    ✗ Index failed: {result.stderr[:100]}")
        tbi_paths.append(tbi_path)
    return tbi_paths


def main():
    parser = argparse.ArgumentParser(
        description="Download Vindija Neanderthal reference for archaic comparison")
    parser.add_argument("--output-dir", default=str(REF_DIR))
    parser.add_argument("--chromosomes", default="chr22",
                        help="Comma-separated list (e.g. chr1,chr2,chr22). Default: chr22")
    parser.add_argument("--all", action="store_true",
                        help="Download all autosomes (~50 GB)")
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()

    output_dir = Path(args.output_dir)

    if args.all:
        chromosomes = [f"chr{c}" for c in range(1, 23)]
    else:
        chromosomes = [c.strip() for c in args.chromosomes.split(",")]

    total_files = len(chromosomes)
    print("═══ Vindija Neanderthal Reference Download ═══\n")
    print(f"Source:  {VCF_BASE_URL}/")
    print(f"Genome:  Vindija 33.19 (~30x coverage, hg19, MQ≥25)")
    print(f"Chromosomes: {', '.join(chromosomes)}")
    print(f"Output:  {output_dir}/")
    print(f"Paper:   Prüfer et al. (2017) Science\n")

    try:
        print(f"1. Downloading {total_files} chromosome file(s)...")
        vcf_paths = download_vindija_vcf(output_dir, chromosomes, args.dry_run)

        if args.dry_run:
            print(f"\n  Total: {total_files} file(s), ~{total_files * 575 / 22:.0f} MB (chr22 extrapolation)")
            print("  Dry run complete.")
            return 0

        if not vcf_paths:
            print("  ✗ No files downloaded")
            return 1

        print(f"\n2. Indexing {len(vcf_paths)} file(s)...")
        index_vcfs(vcf_paths)

        manifest = {
            "source": "Max Planck Institute for Evolutionary Anthropology",
            "url": VCF_BASE_URL,
            "specimen": "Vindija 33.19",
            "coverage": "~30x",
            "genome_build": "hg19/GRCh37",
            "filter": "MQ≥25, MAPQ≥100",
            "citation": "Prüfer et al. (2017) Science. doi:10.1126/science.aao1887",
            "download_date": pd.Timestamp.now().isoformat(),
            "chromosomes_downloaded": chromosomes,
            "files": [p.name for p in vcf_paths],
        }
        manifest_path = output_dir / "vindija_manifest.json"
        with open(manifest_path, "w") as fh:
            json.dump(manifest, fh, indent=2)
        print(f"\n  ✓ {manifest_path}")

        print(f"\n═══ Done: {output_dir}/ ═══")
        print(f"  To add more chromosomes: --chromosomes chr1,chr2,...")
        print(f"  To download all autosomes: --all (~50 GB)")
        return 0

    except Exception as e:
        print(f"\n  ✗ Failed: {e}")
        return 1


if __name__ == "__main__":
    sys.exit(main())
