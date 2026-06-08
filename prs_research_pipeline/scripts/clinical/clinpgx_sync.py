#!/usr/bin/env python3
"""
╔══════════════════════════════════════════════════════════════════════════════╗
║   ClinPGx DOWNLOAD — Pharmacogenomic Reference Data                          ║
║   scripts/clinical/clinpgx_sync.py                                           ║
║                                                                            ║
║   Downloads pharmacogenomic datasets from ClinPGx (free, no login).         ║
║   Replaces per-variant API calls with bulk ZIP downloads.                   ║
║                                                                            ║
║   Files downloaded:                                                         ║
║     • clinicalVariants.zip (72 KB)  — variant-drug pairs, evidence levels  ║
║     • guidelineAnnotations.json.zip (836 KB) — CPIC clinical guidelines     ║
║     • variants.zip (873 KB) — dbSNP-mapped variants                        ║
║                                                                            ║
║   Source: https://api.clinpgx.org/v1/download/file/data/                    ║
║   License: CC BY-SA 4.0 (ClinPGx/PharmGKB)                                 ║
║   Rate limit: 2 requests/second                                             ║
╚══════════════════════════════════════════════════════════════════════════════╝
"""

import sys
import os
import json
import logging
import zipfile
import csv
import io
import time
from pathlib import Path
from urllib.request import urlopen, Request
from urllib.error import HTTPError, URLError
from datetime import datetime, timezone
from typing import Dict, List, Optional

logger = logging.getLogger(__name__)

DOWNLOAD_BASE = "https://api.clinpgx.org/v1/download/file/data"

# Files to download (URL path → local name)
DATASETS = {
    "clinicalVariants.zip": "Clinical variant-drug pairs with evidence levels (TSV)",
    "guidelineAnnotations.json.zip": "CPIC/DPWG clinical guideline annotations (JSON)",
    "variants.zip": "dbSNP-mapped variants annotated by ClinPGx (TSV)",
}

HEADERS = {
    "User-Agent": "PRS-Research-Pipeline/10.0 (research@localhost)",
    "Accept": "*/*",
}


def download_dataset(name: str, output_dir: Path, retries: int = 3) -> Optional[Path]:
    """Download a dataset from ClinPGx. Returns path to downloaded file."""
    url = f"{DOWNLOAD_BASE}/{name}"
    dest = output_dir / name

    for attempt in range(retries):
        try:
            logger.info(f"  Downloading {name}...")
            req = Request(url, headers=HEADERS)
            with urlopen(req, timeout=120) as resp:
                dest.write_bytes(resp.read())
            size_kb = dest.stat().st_size / 1024
            logger.info(f"    ✓ {size_kb:.0f} KB")
            return dest
        except (HTTPError, URLError) as e:
            if attempt < retries - 1:
                wait = 2 ** attempt
                logger.warning(f"    Retry {attempt + 1}/{retries} in {wait}s: {e}")
                time.sleep(wait)
            else:
                logger.error(f"    ✗ Failed: {e}")
                return None
    return None


def parse_clinical_variants(zip_path: Path) -> List[Dict]:
    """Parse clinicalVariants.zip TSV into list of variant-drug dicts.
    Columns: variant, gene, type, level of evidence, chemicals, phenotypes
    """
    results = []
    with zipfile.ZipFile(zip_path) as zf:
        for name in zf.namelist():
            if name.endswith(".tsv"):
                with zf.open(name) as fh:
                    reader = csv.DictReader(io.TextIOWrapper(fh, encoding="utf-8"), delimiter="\t")
                    for row in reader:
                        results.append({
                            "gene": row.get("gene", ""),
                            "drug": row.get("chemicals", ""),
                            "variant": row.get("variant", ""),  # rsID or star allele
                            "type": row.get("type", ""),
                            "evidence_level": row.get("level of evidence", ""),
                            "phenotypes": row.get("phenotypes", ""),
                        })
    return results


def parse_guidelines(zip_path: Path) -> List[Dict]:
    """Parse guidelineAnnotations.json.zip into list of guideline dicts.
    Each JSON file has: {citations: [...], guideline: {...}} or just the guideline object.
    """
    results = []
    with zipfile.ZipFile(zip_path) as zf:
        for name in zf.namelist():
            if name.endswith(".json"):
                with zf.open(name) as fh:
                    data = json.load(fh)

                # Handle wrapper: {"citations": [...], "guideline": {...}}
                if isinstance(data, dict) and "guideline" in data:
                    g = data["guideline"]
                else:
                    g = data

                if not isinstance(g, dict):
                    continue

                chemicals = g.get("relatedChemicals", [])
                genes = g.get("relatedGenes", [])
                results.append({
                    "name": g.get("name", ""),
                    "source": g.get("source", ""),
                    "drugs": [c.get("name", "") for c in chemicals],
                    "genes": [g2.get("symbol", "") for g2 in genes],
                    "summary_markdown": g.get("summaryMarkdown", ""),
                    "text_markdown": g.get("textMarkdown", ""),
                    "recommendation": json.dumps(g.get("recommendation", "")) if g.get("recommendation") else "",
                    "pmids": [lit.get("pmid", "") for lit in g.get("literature", [])],
                })
    return results


def parse_variants(zip_path: Path) -> Dict[str, Dict]:
    """Parse variants.zip TSV into dict keyed by rsID.
    Columns: Variant ID, Variant Name (rsID), Gene Symbols, Location, ...
    """
    results = {}
    with zipfile.ZipFile(zip_path) as zf:
        for name in zf.namelist():
            if name.endswith(".tsv"):
                with zf.open(name) as fh:
                    reader = csv.DictReader(io.TextIOWrapper(fh, encoding="utf-8"), delimiter="\t")
                    for row in reader:
                        rsid = row.get("Variant Name", "")
                        if rsid and rsid.startswith("rs"):
                            results[rsid] = {
                                "gene": row.get("Gene Symbols", ""),
                                "clinpgx_id": row.get("Variant ID", ""),
                                "location": row.get("Location", ""),
                                "clinical_annotations": row.get("Clinical Annotation count", "0"),
                                "guideline_annotations": row.get("Guideline Annotation count", "0"),
                            }
    return results


def sync_all(output_dir: str = "reference/clinpgx", force: bool = False) -> Dict:
    """Download and parse all ClinPGx datasets."""
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    logger.info("═══ ClinPGx Data Download ═══")
    logger.info(f"  Output: {out}/")

    result = {"datasets": {}, "stats": {}}

    for name, description in DATASETS.items():
        logger.info(f"\n  {description}")
        dest = out / name

        if dest.exists() and not force:
            logger.info(f"    Already cached — skipping")
        else:
            downloaded = download_dataset(name, out)
            if not downloaded:
                logger.error(f"    Failed to download {name}")
                continue
            time.sleep(0.6)  # Rate limit: 2 req/sec

        result["datasets"][name] = str(dest)

    # Parse and summarize
    cv_path = out / "clinicalVariants.zip"
    if cv_path.exists():
        cv_data = parse_clinical_variants(cv_path)
        result["stats"]["clinical_variants"] = len(cv_data)
        result["clinical_variants"] = cv_data
        logger.info(f"\n  Clinical variant-drug pairs: {len(cv_data)}")

        # Show high-evidence examples
        level_1a = [v for v in cv_data if v.get("evidence_level") == "1A"]
        logger.info(f"    Level 1A (highest): {len(level_1a)}")
        for v in level_1a[:5]:
            logger.info(f"      {v['gene']} → {v['drug']} ({v['variant']})")

    gl_path = out / "guidelineAnnotations.json.zip"
    if gl_path.exists():
        gl_data = parse_guidelines(gl_path)
        result["stats"]["guidelines"] = len(gl_data)
        result["guidelines"] = gl_data
        logger.info(f"  CPIC guidelines: {len(gl_data)}")
        for g in gl_data[:3]:
            logger.info(f"    {g['source'].upper()}: {g['name'][:80]}...")

    var_path = out / "variants.zip"
    if var_path.exists():
        var_data = parse_variants(var_path)
        result["stats"]["variants"] = len(var_data)
        result["variants"] = var_data
        logger.info(f"  dbSNP-mapped variants: {len(var_data)}")

    # Save manifest
    manifest = {
        "download_date": datetime.now(timezone.utc).isoformat(),
        "source": "api.clinpgx.org (CC BY-SA 4.0)",
        "datasets": result["stats"],
    }
    manifest_path = out / "clinpgx_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2))

    # Save parsed data as JSON for offline use
    parsed = {
        "clinical_variants": result.get("clinical_variants", []),
        "guidelines": result.get("guidelines", []),
        "variants": result.get("variants", {}),
        "manifest": manifest,
    }
    parsed_path = out / "clinpgx_parsed.json"
    parsed_path.write_text(json.dumps(parsed, indent=2, ensure_ascii=False))
    logger.info(f"\n  ✅ Parsed database: {parsed_path} ({parsed_path.stat().st_size / 1024:.0f} KB)")
    logger.info(f"  ✅ Ready for pharmacogenomic annotation")

    return result


def main():
    import argparse

    parser = argparse.ArgumentParser(description="ClinPGx data download for pharmacogenomic annotation")
    parser.add_argument("--sync", action="store_true", help="Download and parse all ClinPGx datasets")
    parser.add_argument("--force", "-f", action="store_true", help="Force re-download")
    parser.add_argument("--output-dir", "-o", default="reference/clinpgx", help="Output directory")
    parser.add_argument("--verbose", "-v", action="store_true")
    parser.add_argument("--show-stats", action="store_true", help="Show cached stats")
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(message)s",
    )

    if args.show_stats:
        manifest_path = Path(args.output_dir) / "clinpgx_manifest.json"
        if manifest_path.exists():
            logger.info(manifest_path.read_text())
        else:
            logger.warning("No cache found. Run --sync first.")
        return 0

    if args.sync:
        result = sync_all(args.output_dir, args.force)
        return 0 if result.get("datasets") else 1

    parser.print_help()
    return 0


if __name__ == "__main__":
    sys.exit(main())
