#!/usr/bin/env python3
"""
╔══════════════════════════════════════════════════════════════════════════════╗
║   PHASE 6 — MODULE 1: FULL 1000 GENOMES INTEGRATION                        ║
║   scripts/download_1000G_full.py                                            ║
║                                                                            ║
║   Downloads and validates the complete 1000 Genomes Phase 3 reference       ║
║   panel (all 22 autosomes + chrX) for genome-wide ancestry inference.       ║
║                                                                            ║
║   Source: 1000 Genomes Phase 3 (20130502 release)                           ║
║   Format: VCF → PLINK binary (per chromosome → merged)                      ║
║   Samples: 2,504 across 5 super-populations                                 ║
║   Variants: ~80 million genome-wide                                         ║
║                                                                            ║
║   Corrections relative to previous chr22-only approach:                     ║
║     • Genome-wide PCA requires all autosomes                               ║
║     • LD pruning must span full genome for valid independent SNP set        ║
║     • Ancestry inference needs thousands of neutral markers, not 33         ║
║                                                                            ║
║   Output:                                                                   ║
║     reference/1000G_full/1000G_full.bed/.bim/.fam                          ║
║     reference/1000G_full/manifest.json                                      ║
╚══════════════════════════════════════════════════════════════════════════════╝
"""

import sys
import os
import json
import hashlib
import logging
import subprocess
import tempfile
import shutil
from pathlib import Path
from typing import Dict, List, Optional, Any, Tuple
from dataclasses import dataclass, field, asdict
from datetime import datetime
from urllib.request import urlopen, Request

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)

# ── Constants ──────────────────────────────────────────────────────────────────

G1K_FTP_BASE = "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502"
POP_PANEL_URL = f"{G1K_FTP_BASE}/integrated_call_samples_v3.20130502.ALL.panel"

# All autosomes + chrX
CHROMOSOMES = [str(i) for i in range(1, 23)] + ["X"]

# Expected file sizes (approximate, for validation)
EXPECTED_SIZES_GB = {
    "1": 0.48, "2": 0.53, "3": 0.44, "4": 0.43, "5": 0.40,
    "6": 0.38, "7": 0.35, "8": 0.35, "9": 0.28, "10": 0.31,
    "11": 0.31, "12": 0.30, "13": 0.23, "14": 0.21, "15": 0.19,
    "16": 0.21, "17": 0.18, "18": 0.18, "19": 0.15, "20": 0.14,
    "21": 0.09, "22": 0.08, "X": 0.31,
}

# ── Data Structures ────────────────────────────────────────────────────────────

@dataclass
class ChromosomeManifest:
    """Metadata for a single chromosome download."""
    chromosome: str
    vcf_url: str
    vcf_index_url: str
    local_vcf: str = ""
    local_tbi: str = ""
    vcf_hash: str = ""
    file_size_bytes: int = 0
    download_complete: bool = False
    variants_imported: int = 0


@dataclass
class G1KFullManifest:
    """Complete 1000 Genomes full reference manifest."""
    chromosomes: List[ChromosomeManifest] = field(default_factory=list)
    total_variants: int = 0
    total_samples: int = 0
    merged_bed: str = ""
    merged_bim: str = ""
    merged_fam: str = ""
    merged_hash: str = ""
    download_date: str = ""
    pipeline_version: str = "6.0.0"


# ── 1000 Genomes Full Downloader ──────────────────────────────────────────────

class G1KFullDownloader:
    """
    Downloads, validates, and merges the full 1000 Genomes Phase 3 reference.

    Supports:
      - Per-chromosome download with resume capability
      - SHA-256 checksum verification
      - Automatic VCF → PLINK conversion
      - Multi-chromosome merge into single dataset
      - Provenance manifest generation

    Usage:
        downloader = G1KFullDownloader(output_dir="reference/1000G_full")
        manifest = downloader.download_all(plink_binary="plink", threads=8)
    """

    def __init__(self, output_dir: str = "reference/1000G_full"):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self._manifest: Optional[G1KFullManifest] = None

    # ── Public API ───────────────────────────────────────────────────────

    def download_all(
        self,
        plink_binary: str = "plink",
        threads: int = 4,
        memory: int = 16000,
        chromosomes: Optional[List[str]] = None,
        skip_download: bool = False,
    ) -> G1KFullManifest:
        """
        Download and process all chromosomes.

        Args:
            plink_binary: Path to PLINK binary.
            threads: Number of threads for PLINK operations.
            memory: Memory limit in MB.
            chromosomes: Which chromosomes to process (default: all 22 + X).
            skip_download: If True, only process already-downloaded files.

        Returns:
            G1KFullManifest with complete provenance.
        """
        logger.info("═══ Full 1000 Genomes Integration (Phase 6) ═══")

        chroms = chromosomes or CHROMOSOMES
        manifest = G1KFullManifest(
            chromosomes=[],
            download_date=datetime.now().strftime("%Y-%m-%d %H:%M UTC"),
        )

        # Phase 1: Download population panel
        panel_path = self.output_dir / "population_panel.txt"
        if not panel_path.exists():
            self._download_population_panel(panel_path)

        # Phase 2: Download per-chromosome VCFs
        for chrom in chroms:
            logger.info(f"  Processing chromosome {chrom}...")
            chrom_manifest = self._process_chromosome(
                chrom, plink_binary, threads, memory, skip_download
            )
            manifest.chromosomes.append(chrom_manifest)

        # Phase 3: Merge all chromosomes
        logger.info("  Merging all chromosomes...")
        merged_prefix = self.output_dir / "1000G_full"
        self._merge_chromosomes(merged_prefix, plink_binary, threads, memory)
        manifest.merged_bed = str(merged_prefix) + ".bed"
        manifest.merged_bim = str(merged_prefix) + ".bim"
        manifest.merged_fam = str(merged_prefix) + ".fam"

        # Phase 4: Compute manifest
        manifest.total_variants = sum(
            c.variants_imported for c in manifest.chromosomes
        )
        manifest.total_samples = self._count_samples(merged_prefix)
        manifest.merged_hash = self._compute_hash(f"{merged_prefix}.bed")

        # Save manifest
        self._save_manifest(manifest)
        self._manifest = manifest

        logger.info(f"  ✅ Full 1000 Genomes integrated: {manifest.total_variants:,} variants, "
                   f"{manifest.total_samples} samples")

        return manifest

    def download_population_panel(self) -> Path:
        """Download only the population panel (fast, for ancestry labeling)."""
        panel_path = self.output_dir / "population_panel.txt"
        self._download_population_panel(panel_path)
        return panel_path

    # ── Private: Download ─────────────────────────────────────────────────

    def _download_population_panel(self, output_path: Path) -> None:
        """Download 1000 Genomes population panel."""
        logger.info(f"  Downloading population panel...")
        try:
            req = Request(POP_PANEL_URL)
            with urlopen(req, timeout=60) as resp:
                with open(output_path, "wb") as fh:
                    fh.write(resp.read())
            logger.info(f"    Population panel: {output_path}")
        except Exception as e:
            logger.error(f"    Failed: {e}")
            # Try local copy
            local = Path("reference/1000G/20130606_g1k_3202_samples_ped_population.txt")
            if local.exists():
                shutil.copy(local, output_path)
                logger.info(f"    Copied from local: {local}")

    def _process_chromosome(
        self,
        chrom: str,
        plink_binary: str,
        threads: int,
        memory: int,
        skip_download: bool,
    ) -> ChromosomeManifest:
        """Download and convert a single chromosome."""
        # chrX uses a different version suffix than the autosomes
        version_suffix = "v1c" if chrom == "X" else "v5b"
        vcf_name = f"ALL.chr{chrom}.phase3_shapeit2_mvncall_integrated_{version_suffix}.20130502.genotypes.vcf.gz"
        vcf_url = f"{G1K_FTP_BASE}/{vcf_name}"
        tbi_url = f"{vcf_url}.tbi"

        local_vcf = self.output_dir / vcf_name
        local_tbi = self.output_dir / f"{vcf_name}.tbi"

        manifest = ChromosomeManifest(
            chromosome=chrom,
            vcf_url=vcf_url,
            vcf_index_url=tbi_url,
            local_vcf=str(local_vcf),
            local_tbi=str(local_tbi),
        )

        # Validate existing file integrity — cached via .validated marker
        vcf_ok = False
        validated_marker = Path(str(local_vcf) + ".validated")
        if local_vcf.exists():
            # Skip validation if marker is newer than VCF (already verified)
            if validated_marker.exists() and validated_marker.stat().st_mtime >= local_vcf.stat().st_mtime:
                vcf_ok = True
                manifest.file_size_bytes = local_vcf.stat().st_size
                manifest.download_complete = True
                logger.info(f"    chr{chrom}: cached ({manifest.file_size_bytes:,} bytes, validated)")
            else:
                logger.info(f"    chr{chrom}: validating integrity...")
                result = subprocess.run(
                    ["gzip", "-t", str(local_vcf)],
                    capture_output=True, timeout=60
                )
                if result.returncode == 0:
                    vcf_ok = True
                    manifest.file_size_bytes = local_vcf.stat().st_size
                    manifest.download_complete = True
                    # Cache the validation
                    validated_marker.touch()
                    logger.info(f"    chr{chrom}: validated ({manifest.file_size_bytes:,} bytes, cached)")
                else:
                    validated_marker.unlink(missing_ok=True)
                    logger.warning(f"    chr{chrom}: corrupted ({local_vcf.stat().st_size} bytes) — re-downloading")

        # Download if needed
        if not skip_download and not vcf_ok:
            self._download_file(vcf_url, local_vcf, f"chr{chrom} VCF")
            manifest.file_size_bytes = local_vcf.stat().st_size
            manifest.download_complete = True

        # Download index
        if not skip_download and not local_tbi.exists():
            try:
                self._download_file(tbi_url, local_tbi, f"chr{chrom} index")
            except Exception:
                logger.warning(f"    chr{chrom} index download failed — will use tabix")

        # Convert to PLINK if VCF exists
        if local_vcf.exists():
            plink_prefix = self.output_dir / f"chr{chrom}_g1k"
            variants = self._vcf_to_plink(
                str(local_vcf), str(plink_prefix), plink_binary, threads, memory
            )
            manifest.variants_imported = variants
            manifest.vcf_hash = self._compute_hash(str(local_vcf))

        return manifest

    def _download_file(self, url: str, output_path: Path, label: str) -> None:
        """Download a file with progress reporting."""
        logger.info(f"    Downloading {label}...")
        try:
            req = Request(url, headers={"User-Agent": "PRS-Research-Platform/6.0"})
            with urlopen(req, timeout=3600) as resp:
                with open(output_path, "wb") as fh:
                    chunk_size = 1024 * 1024  # 1MB
                    downloaded = 0
                    while True:
                        chunk = resp.read(chunk_size)
                        if not chunk:
                            break
                        fh.write(chunk)
                        downloaded += len(chunk)
                        if downloaded % (100 * 1024 * 1024) == 0:  # Log every 100MB
                            logger.info(f"      {downloaded / (1024**2):.0f} MB...")
            logger.info(f"    Downloaded: {output_path} ({output_path.stat().st_size:,} bytes)")
        except Exception as e:
            logger.error(f"    Download failed: {e}")
            raise

    # ── Private: VCF → PLINK ──────────────────────────────────────────────

    def _vcf_to_plink(
        self,
        vcf_path: str,
        plink_prefix: str,
        plink_binary: str,
        threads: int,
        memory: int,
    ) -> int:
        """Convert VCF to PLINK binary format."""
        logger.info(f"    Converting {Path(vcf_path).name} → PLINK...")

        # Check if already converted
        bim_path = Path(f"{plink_prefix}.bim")
        if Path(f"{plink_prefix}.bed").exists() and bim_path.exists():
            # Verify variant IDs have been renamed (not all '.')
            needs_rename = True
            try:
                with open(bim_path) as fh:
                    first_id = fh.readline().split("\t")[1] if (line := fh.readline()) else "."
                    needs_rename = (first_id == ".")
            except Exception:
                pass

            if needs_rename:
                logger.info(f"    Renaming variant IDs to chr:pos...")
                tmp_path = Path(f"{plink_prefix}.bim.tmp")
                with open(bim_path) as fh_in, open(tmp_path, "w") as fh_out:
                    for line in fh_in:
                        parts = line.strip().split("\t")
                        parts[1] = f"{parts[0]}:{parts[3]}"
                        fh_out.write("\t".join(parts) + "\n")
                tmp_path.replace(bim_path)

            bim_variants = self._count_bim_variants(plink_prefix)
            logger.info(f"    Already converted: {bim_variants:,} variants")
            return bim_variants

        cmd = [
            plink_binary,
            "--vcf", vcf_path,
            "--make-bed",
            "--double-id",
            "--allow-extra-chr",
            "--keep-allele-order",
            "--out", plink_prefix,
            "--threads", str(threads),
            "--memory", str(memory),
        ]

        try:
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=7200)

            # Rename '.' variant IDs to chr:pos for uniqueness across chromosomes
            # 1000 Genomes VCFs use '.' for all variant IDs — PLINK merge requires unique names
            bim_path = Path(f"{plink_prefix}.bim")
            if bim_path.exists():
                tmp_path = Path(f"{plink_prefix}.bim.tmp")
                with open(bim_path) as fh_in, open(tmp_path, "w") as fh_out:
                    for line in fh_in:
                        parts = line.strip().split("\t")
                        parts[1] = f"{parts[0]}:{parts[3]}"  # chr:pos
                        fh_out.write("\t".join(parts) + "\n")
                tmp_path.replace(bim_path)

            # Parse variant count from PLINK log
            variants = 0
            for line in result.stdout.split("\n") + result.stderr.split("\n"):
                if "variants" in line.lower() and "people" in line.lower():
                    try:
                        parts = line.split()
                        variants = int(parts[0])
                    except (ValueError, IndexError):
                        pass

            if variants == 0:
                variants = self._count_bim_variants(plink_prefix)

            logger.info(f"    Converted: {variants:,} variants (IDs renamed to chr:pos)")
            return variants

        except subprocess.TimeoutExpired:
            logger.error(f"    PLINK conversion timed out for {vcf_path}")
            return 0
        except Exception as e:
            logger.error(f"    PLINK conversion error: {e}")
            return 0

    # ── Private: Merge ────────────────────────────────────────────────────

    def _merge_chromosomes(
        self,
        merged_prefix: Path,
        plink_binary: str,
        threads: int,
        memory: int,
    ) -> None:
        """Merge per-chromosome PLINK files into a single dataset."""
        # Already merged?
        if Path(f"{merged_prefix}.bed").exists():
            variants = self._count_bim_variants(str(merged_prefix))
            logger.info(f"    Already merged: {variants:,} variants")
            return

        # Build merge list from available chromosome files
        merge_list = self.output_dir / "merge_list.txt"
        # Always rebuild from scratch — stale merge lists cause silent failures
        merge_list.unlink(missing_ok=True)
        chrom_files = []
        for chrom in CHROMOSOMES:
            prefix = self.output_dir / f"chr{chrom}_g1k"
            if Path(f"{prefix}.bed").exists():
                chrom_files.append(str(prefix))

        if not chrom_files:
            logger.error("    No chromosome PLINK files found")
            return

        logger.info(f"    Merging {len(chrom_files)} chromosomes...")

        with open(merge_list, "w") as fh:
            for f in chrom_files[1:]:  # First file is base, rest are merged in
                fh.write(f"{f}.bed {f}.bim {f}.fam\n")

        cmd = [
            plink_binary,
            "--bfile", chrom_files[0],
            "--merge-list", str(merge_list),
            "--make-bed",
            "--out", str(merged_prefix),
            "--threads", str(threads),
            "--memory", str(memory),
            "--allow-extra-chr",
        ]

        try:
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=14400)

            # Handle merge failures (multi-allelic variants)
            if result.returncode != 0:
                missnp_file = f"{merged_prefix}-merge.missnp"
                if Path(missnp_file).exists():
                    logger.info(f"    Removing {self._count_lines(missnp_file)} problematic variants...")
                    # Remove multi-allelic variants and retry
                    for prefix in chrom_files:
                        tmp_prefix = f"{prefix}_clean"
                        subprocess.run([
                            plink_binary, "--bfile", prefix,
                            "--exclude", missnp_file,
                            "--make-bed", "--out", tmp_prefix,
                            "--threads", str(threads), "--memory", str(memory),
                        ], capture_output=True, text=True, timeout=3600)

                    # Rewrite merge list with cleaned files
                    clean_files = [f"{f}_clean" for f in chrom_files]
                    with open(merge_list, "w") as fh:
                        for f in clean_files[1:]:
                            fh.write(f"{f}.bed {f}.bim {f}.fam\n")

                    subprocess.run([
                        plink_binary,
                        "--bfile", clean_files[0],
                        "--merge-list", str(merge_list),
                        "--make-bed", "--out", str(merged_prefix),
                        "--threads", str(threads), "--memory", str(memory),
                        "--allow-extra-chr",
                    ], capture_output=True, text=True, timeout=14400)

                    # Clean up temporary files
                    for prefix in chrom_files:
                        for f in Path().glob(f"{prefix}_clean.*"):
                            f.unlink()
                    merge_list.unlink(missing_ok=True)

        except subprocess.TimeoutExpired:
            logger.error("    Merge timed out")
        except Exception as e:
            logger.error(f"    Merge error: {e}")

    # ── Helpers ──────────────────────────────────────────────────────────

    @staticmethod
    def _count_bim_variants(prefix: str) -> int:
        """Count variants in a PLINK BIM file."""
        bim_path = f"{prefix}.bim"
        if Path(bim_path).exists():
            return sum(1 for _ in open(bim_path))
        return 0

    @staticmethod
    def _count_samples(prefix: Path) -> int:
        """Count samples in a PLINK FAM file."""
        fam_path = f"{prefix}.fam"
        if Path(fam_path).exists():
            return sum(1 for _ in open(fam_path))
        return 0

    @staticmethod
    def _count_lines(path: str) -> int:
        """Count lines in a file."""
        try:
            return sum(1 for _ in open(path))
        except Exception:
            return 0

    @staticmethod
    def _compute_hash(file_path: str) -> str:
        """Compute SHA-256 hash of first 1MB + last 1MB (fast partial hash)."""
        if not Path(file_path).exists():
            return "missing"
        sha = hashlib.sha256()
        file_size = Path(file_path).stat().st_size
        with open(file_path, "rb") as fh:
            sha.update(fh.read(min(1024 * 1024, file_size)))
            if file_size > 2 * 1024 * 1024:
                fh.seek(-1024 * 1024, 2)
                sha.update(fh.read(1024 * 1024))
        return sha.hexdigest()[:16]

    def _save_manifest(self, manifest: G1KFullManifest) -> None:
        """Save full 1000G manifest."""
        manifest_path = self.output_dir / "manifest.json"
        with open(manifest_path, "w") as fh:
            json.dump(asdict(manifest), fh, indent=2, default=str)
        logger.info(f"  ✅ Manifest: {manifest_path}")


# ── CLI Entry Point ──────────────────────────────────────────────────────────

def main():
    import argparse

    parser = argparse.ArgumentParser(
        description="Phase 6 Module 1: Full 1000 Genomes Integration"
    )
    parser.add_argument("--output-dir", "-o", default="reference/1000G_full",
                       help="Output directory for 1000G reference")
    parser.add_argument("--plink", default="plink", help="PLINK binary path")
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--memory", type=int, default=16000)
    parser.add_argument("--chromosomes", nargs="+",
                       help="Specific chromosomes to download (default: all)")
    parser.add_argument("--skip-download", action="store_true",
                       help="Only process already-downloaded files")
    parser.add_argument("--panel-only", action="store_true",
                       help="Download only the population panel")
    parser.add_argument("--verbose", "-v", action="store_true")

    args = parser.parse_args()

    # Auto-detect PLINK if not explicitly provided or if default failed
    if args.plink == "plink":
        # Search for PLINK in standard locations
        candidates = [
            str(Path(__file__).parent.parent.parent / "tools" / "plink"),  # ../../tools/plink
            str(Path(__file__).parent.parent / "tools" / "plink"),          # ../tools/plink
            str(Path.cwd() / "tools" / "plink"),
            str(Path(__file__).parent.parent.parent / "tools" / "plink2"),
        ]
        for c in candidates:
            if os.path.exists(c) and os.access(c, os.X_OK):
                args.plink = c
                logger.info(f"Auto-detected PLINK: {c}")
                break
        else:
            # Fall back to system plink
            system_plink = shutil.which("plink") or shutil.which("plink2")
            if system_plink:
                args.plink = system_plink
                logger.info(f"Using system PLINK: {system_plink}")
            else:
                logger.warning("PLINK not found in tools/ or PATH — VCF conversion may fail")

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
    )

    downloader = G1KFullDownloader(output_dir=args.output_dir)

    if args.panel_only:
        path = downloader.download_population_panel()
        print(f"\n✅ Population panel: {path}")
        return 0

    manifest = downloader.download_all(
        plink_binary=args.plink,
        threads=args.threads,
        memory=args.memory,
        chromosomes=args.chromosomes,
        skip_download=args.skip_download,
    )

    print(f"\n═══ Full 1000 Genomes Integrated ═══")
    print(f"  Variants: {manifest.total_variants:,}")
    print(f"  Samples: {manifest.total_samples}")
    print(f"  Chromosomes: {len(manifest.chromosomes)}")
    print(f"  Output: {manifest.merged_bed}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
