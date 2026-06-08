#!/usr/bin/env python3
"""
╔══════════════════════════════════════════════════════════════════════════════╗
║   CLINVAR VCF DOWNLOAD — GRCh37 Reference                                    ║
║   scripts/setup/download_clinvar.py                                          ║
║                                                                            ║
║   Downloads the ClinVar germline VCF (GRCh37/hg19) from NCBI FTP.           ║
║   Cached locally at reference/clinvar/ for subsequent annotation runs.      ║
║                                                                            ║
║   Files downloaded:                                                         ║
║     • clinvar.vcf.gz      (~183 MB) — aggregate germline classifications   ║
║     • clinvar.vcf.gz.tbi  (~596 KB) — tabix index                          ║
║     • clinvar.vcf.gz.md5  (~132 B)  — MD5 checksum                         ║
║                                                                            ║
║   Source: ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/                      ║
║   Format: VCFv4.1, ~4.4M variants                                           ║
║                                                                            ║
║   Output:                                                                   ║
║     reference/clinvar/clinvar.vcf.gz                                        ║
║     reference/clinvar/clinvar.vcf.gz.tbi                                    ║
║     reference/clinvar/manifest.json                                         ║
╚══════════════════════════════════════════════════════════════════════════════╝
"""

import sys
import os
import json
import logging
import hashlib
import time
from pathlib import Path
from datetime import datetime, timezone
from urllib.request import urlopen, Request
from urllib.error import HTTPError, URLError
from urllib.parse import urlparse
import shutil

logger = logging.getLogger(__name__)

# ── Constants ──────────────────────────────────────────────────────────────────

CLINVAR_FTP_BASE = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37"

FILES = [
    "clinvar.vcf.gz",
    "clinvar.vcf.gz.tbi",
    "clinvar.vcf.gz.md5",
]

CHUNK_SIZE = 1024 * 1024  # 1 MB for progress reporting


def _progress(msg: str):
    logger.info(f"  {msg}")


def download_file(url: str, dest: Path, retries: int = 3) -> bool:
    """Download a file from URL to dest with retry logic. Returns True on success."""
    dest.parent.mkdir(parents=True, exist_ok=True)

    for attempt in range(retries):
        try:
            _progress(f"Downloading {url.rsplit('/', 1)[-1]} ...")
            req = Request(url)
            with urlopen(req, timeout=120) as resp:
                total = int(resp.headers.get("Content-Length", 0))
                downloaded = 0
                with open(dest, "wb") as fh:
                    while True:
                        chunk = resp.read(CHUNK_SIZE)
                        if not chunk:
                            break
                        fh.write(chunk)
                        downloaded += len(chunk)
                        if total > 0 and downloaded % (10 * CHUNK_SIZE) == 0:
                            pct = downloaded / total * 100
                            _progress(f"  {pct:.0f}% ({downloaded / (1024*1024):.0f}/{total / (1024*1024):.0f} MB)")
            _progress(f"  ✓ {dest.name} ({downloaded / (1024*1024):.0f} MB)")
            return True
        except (HTTPError, URLError, OSError) as e:
            if attempt < retries - 1:
                wait = 2 ** attempt
                logger.warning(f"  Retry {attempt + 1}/{retries} in {wait}s: {e}")
                time.sleep(wait)
            else:
                logger.error(f"  ✗ Failed to download {url}: {e}")
                return False
    return False


def verify_md5(filepath: Path, md5_path: Path) -> bool:
    """Verify file against its MD5 checksum file."""
    if not md5_path.exists():
        logger.warning(f"  No MD5 file found at {md5_path} — skipping verification")
        return True  # No MD5 to check against, assume OK

    expected_md5 = md5_path.read_text().strip().split()[0]

    _progress(f"Verifying MD5 for {filepath.name} ...")
    hasher = hashlib.md5()
    with open(filepath, "rb") as fh:
        while True:
            chunk = fh.read(CHUNK_SIZE)
            if not chunk:
                break
            hasher.update(chunk)
    actual_md5 = hasher.hexdigest()

    if actual_md5.lower() == expected_md5.lower():
        _progress(f"  ✓ MD5 verified: {actual_md5}")
        return True
    else:
        logger.error(f"  ✗ MD5 mismatch! Expected {expected_md5}, got {actual_md5}")
        return False


def count_clinvar_variants(vcf_path: Path) -> int:
    """Count variants in a gzipped VCF (skip header lines)."""
    import gzip
    count = 0
    try:
        with gzip.open(vcf_path, "rt") as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                count += 1
    except Exception as e:
        logger.warning(f"  Could not count variants: {e}")
        return 0
    return count


def download_clinvar(output_dir: str = "reference/clinvar", force: bool = False) -> dict:
    """
    Download ClinVar VCF (GRCh37) reference files.

    Args:
        output_dir: Directory to store ClinVar reference files.
        force: Re-download even if cached files exist.

    Returns:
        Summary dict with download status and file info.
    """
    logger.info("═══ ClinVar VCF Download (GRCh37) ═══")

    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    vcf_path = out / "clinvar.vcf.gz"
    tbi_path = out / "clinvar.vcf.gz.tbi"
    md5_path = out / "clinvar.vcf.gz.md5"
    manifest_path = out / "manifest.json"

    # Check if already cached
    if not force and vcf_path.exists() and tbi_path.exists():
        if verify_md5(vcf_path, md5_path):
            logger.info("  ClinVar VCF already cached and verified — skipping download")
            n_variants = count_clinvar_variants(vcf_path)
            return {
                "downloaded": False,
                "cached": True,
                "n_variants": n_variants,
                "path": str(vcf_path),
                "index": str(tbi_path),
            }
        else:
            logger.warning("  Cached ClinVar VCF failed MD5 check — re-downloading")

    # Download all files
    success = True
    for filename in FILES:
        url = f"{CLINVAR_FTP_BASE}/{filename}"
        dest = out / filename
        if not download_file(url, dest):
            success = False

    if not success:
        logger.error("  Some ClinVar downloads failed. Re-run to retry.")
        return {"downloaded": False, "error": "Download failed"}

    # Verify MD5
    if vcf_path.exists() and md5_path.exists():
        if not verify_md5(vcf_path, md5_path):
            logger.error("  MD5 verification failed — re-run with --force to re-download")
            return {"downloaded": False, "error": "MD5 verification failed"}

    # Count variants
    n_variants = count_clinvar_variants(vcf_path)
    file_size = vcf_path.stat().st_size if vcf_path.exists() else 0

    # Write manifest
    manifest = {
        "source": "ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/",
        "build": "GRCh37",
        "download_date": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "file_size_bytes": file_size,
        "file_size_mb": round(file_size / (1024 * 1024), 1),
        "n_variants": n_variants,
        "md5_verified": True,
        "files": {
            "vcf": "clinvar.vcf.gz",
            "index": "clinvar.vcf.gz.tbi",
            "checksum": "clinvar.vcf.gz.md5",
        },
    }
    manifest_path.write_text(json.dumps(manifest, indent=2))

    logger.info(f"  ✅ ClinVar VCF ready: {n_variants:,} variants, {file_size / (1024*1024):.0f} MB")
    logger.info(f"  Manifest: {manifest_path}")

    return {
        "downloaded": True,
        "cached": False,
        "n_variants": n_variants,
        "path": str(vcf_path),
        "index": str(tbi_path),
        "manifest": str(manifest_path),
    }


def main():
    import argparse

    parser = argparse.ArgumentParser(
        description="Download ClinVar VCF (GRCh37) reference for pathogenic variant annotation"
    )
    parser.add_argument("--output-dir", "-o", default="reference/clinvar",
                        help="Output directory (default: reference/clinvar)")
    parser.add_argument("--force", "-f", action="store_true",
                        help="Force re-download even if cached")
    parser.add_argument("--verbose", "-v", action="store_true",
                        help="Verbose output")
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s" if args.verbose else "%(message)s",
    )

    result = download_clinvar(args.output_dir, args.force)

    if result.get("error"):
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
