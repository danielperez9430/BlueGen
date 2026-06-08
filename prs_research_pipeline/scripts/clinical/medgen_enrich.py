#!/usr/bin/env python3
"""
╔══════════════════════════════════════════════════════════════════════════════╗
║   MEDGEN DISEASE DESCRIPTION ENRICHMENT (LOCAL DATABASE)                     ║
║   scripts/clinical/medgen_enrich.py                                          ║
║                                                                            ║
║   Enriches ClinVar pathogenic variants with human-readable disease          ║
║   descriptions from a local MedGen database download (~10 MB).              ║
║                                                                            ║
║   One-time setup:                                                           ║
║     python scripts/clinical/medgen_enrich.py --download                      ║
║     → Downloads reference/medgen/{NAMES,MGDEF,MGSTY}.RRF.gz                 ║
║                                                                            ║
║   Enrichment (instant after first run):                                     ║
║     python scripts/clinical/medgen_enrich.py --enrich                        ║
║     → Matches disease names against MedGen, adds descriptions               ║
║                                                                            ║
║   Source: NCBI MedGen FTP (free, no license required)                       ║
║   Format: UMLS Rich Release Format (RRF), pipe-delimited                    ║
║   Files: NAMES.RRF.gz (3 MB) + MGDEF.RRF.gz (5 MB) + MGSTY.RRF.gz (1.6 MB) ║
╚══════════════════════════════════════════════════════════════════════════════╝
"""

import sys
import os
import json
import logging
import gzip
import re
import time
from pathlib import Path
from collections import defaultdict
from urllib.request import urlopen, Request
from urllib.error import HTTPError, URLError
from datetime import datetime, timezone

logger = logging.getLogger(__name__)

MEDGEN_FTP = "https://ftp.ncbi.nlm.nih.gov/pub/medgen"
MEDGEN_FILES = ["NAMES.RRF.gz", "MGDEF.RRF.gz", "MGSTY.RRF.gz"]

# ═══════════════════════════════════════════════════════════════════════════════
# DOWNLOAD
# ═══════════════════════════════════════════════════════════════════════════════


def _get_remote_info(url: str) -> dict:
    """Get remote file info (size, last-modified) via HEAD request."""
    try:
        req = Request(url, method="HEAD")
        with urlopen(req, timeout=15) as resp:
            return {
                "size": int(resp.headers.get("Content-Length", 0)),
                "last_modified": resp.headers.get("Last-Modified", ""),
            }
    except Exception:
        return {}


def download_medgen(output_dir: str = "reference/medgen", force: bool = False, max_age_days: int = 30) -> bool:
    """Download MedGen RRF files from NCBI FTP. Skips if cached and not stale."""
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    manifest_path = out / "manifest.json"

    # Check if existing files are stale
    all_exist = all((out / f).exists() for f in MEDGEN_FILES)
    if all_exist and not force:
        # Check manifest age
        if manifest_path.exists():
            try:
                manifest = json.loads(manifest_path.read_text())
                download_date_str = manifest.get("download_date", "")
                if download_date_str:
                    download_date = datetime.fromisoformat(download_date_str)
                    age_days = (datetime.now(timezone.utc) - download_date).days
                    if age_days <= max_age_days:
                        logger.info(f"MedGen database cached ({age_days}d old, max {max_age_days}d) — skipping")
                        return True
                    else:
                        logger.info(f"MedGen database is {age_days}d old (max {max_age_days}d) — checking for updates...")
            except Exception:
                pass

    # Download files
    remote_info = {}
    for filename in MEDGEN_FILES:
        url = f"{MEDGEN_FTP}/{filename}"
        dest = out / filename
        remote = _get_remote_info(url)

        # Skip if remote hasn't changed
        if dest.exists() and not force and remote.get("last_modified"):
            local_mtime = datetime.fromtimestamp(dest.stat().st_mtime, tz=timezone.utc)
            logger.info(f"  {filename}: local {local_mtime.strftime('%Y-%m-%d')} — skipping (unchanged)")
            remote_info[filename] = remote
            continue

        logger.info(f"  Downloading {filename} ({_format_mb(url)} MB)...")
        try:
            req = Request(url)
            with urlopen(req, timeout=120) as resp:
                dest.write_bytes(resp.read())
            logger.info(f"    ✓ {dest.name} ({dest.stat().st_size / 1024:.0f} KB)")
            remote_info[filename] = remote
        except Exception as e:
            logger.error(f"    ✗ Failed: {e}")
            return False

    # Write manifest
    manifest = {
        "source": "ftp.ncbi.nlm.nih.gov/pub/medgen/",
        "download_date": datetime.now(timezone.utc).isoformat(),
        "files": {f: remote_info.get(f, {}) for f in MEDGEN_FILES},
        "update_frequency": "weekly (Wednesdays)",
        "update_command": "python scripts/clinical/medgen_enrich.py --download",
    }
    manifest_path.write_text(json.dumps(manifest, indent=2))
    logger.info(f"  ✅ MedGen database ready ({out}/)")
    logger.info(f"  ℹ️  Next check: run with --download to refresh weekly")
    return True


def check_medgen_update(ref_dir: str = "reference/medgen") -> dict:
    """
    Check if the local MedGen database is up to date.
    Returns dict with status info.
    """
    ref = Path(ref_dir)
    manifest_path = ref / "manifest.json"

    if not manifest_path.exists():
        return {"status": "not_downloaded", "message": "MedGen database not downloaded. Run with --download."}

    try:
        manifest = json.loads(manifest_path.read_text())
        download_date_str = manifest.get("download_date", "")
        if download_date_str:
            download_date = datetime.fromisoformat(download_date_str)
            age_days = (datetime.now(timezone.utc) - download_date).days
            if age_days <= 7:
                return {"status": "current", "age_days": age_days, "download_date": download_date_str[:10]}
            elif age_days <= 30:
                return {"status": "recent", "age_days": age_days, "download_date": download_date_str[:10]}
            else:
                return {"status": "stale", "age_days": age_days, "download_date": download_date_str[:10],
                        "message": f"Database is {age_days} days old. Consider running with --download to update."}
    except Exception:
        return {"status": "error", "message": "Could not read manifest."}

    return {"status": "unknown"}


def _format_mb(url: str) -> str:
    """Get file size in MB from Content-Length header."""
    try:
        with urlopen(Request(url, method="HEAD"), timeout=10) as resp:
            size = int(resp.headers.get("Content-Length", 0))
            return f"{size / (1024*1024):.1f}"
    except Exception:
        return "?"


# ═══════════════════════════════════════════════════════════════════════════════
# DATABASE LOADING
# ═══════════════════════════════════════════════════════════════════════════════


def load_medgen_database(ref_dir: str = "reference/medgen") -> dict:
    """
    Load MedGen RRF files into memory.
    Returns dict with:
      - name_to_cui: {normalized_name: CUI}
      - cui_to_def: {CUI: definition_text}
      - cui_to_name: {CUI: preferred_name}
    """
    ref = Path(ref_dir)
    names_file = ref / "NAMES.RRF.gz"
    def_file = ref / "MGDEF.RRF.gz"
    sty_file = ref / "MGSTY.RRF.gz"

    if not names_file.exists():
        logger.error(f"  MedGen NAMES file not found: {names_file}")
        logger.error(f"  Run: python scripts/clinical/medgen_enrich.py --download")
        logger.error(f"  Or download manually from: ftp.ncbi.nlm.nih.gov/pub/medgen/")
        return {}

    logger.info("  Loading MedGen database...")

    # Step 1: Filter to Disease or Syndrome concepts only
    disease_cuis = set()
    if sty_file.exists():
        with gzip.open(sty_file, "rt", encoding="utf-8", errors="replace") as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                parts = line.strip().split("|")
                if len(parts) >= 4:
                    cui, _, _, sty = parts[0], parts[1], parts[2], parts[3]
                    if "Disease or Syndrome" in sty:
                        disease_cuis.add(cui)
    logger.info(f"    Disease concepts: {len(disease_cuis):,}")

    # Step 2: Load definitions (CUI → definition)
    cui_to_def = {}
    if def_file.exists():
        with gzip.open(def_file, "rt", encoding="utf-8", errors="replace") as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                parts = line.strip().split("|")
                if len(parts) >= 2:
                    cui, definition = parts[0], parts[1]
                    if cui in disease_cuis and definition:
                        cui_to_def[cui] = definition.replace("\r", " ").strip()

    # Filter disease CUIs to only those with definitions
    defined_cuis = set(cui_to_def.keys())
    disease_cuis &= defined_cuis
    logger.info(f"    Defined diseases: {len(disease_cuis):,}")

    # Step 3: Load names (CUI → preferred name)
    cui_to_name = {}
    if names_file.exists():
        with gzip.open(names_file, "rt", encoding="utf-8", errors="replace") as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                parts = line.strip().split("|")
                if len(parts) >= 2:
                    cui, name = parts[0], parts[1]
                    if cui in disease_cuis and name:
                        cui_to_name[cui] = name

    # Step 4: Build name → CUI index (normalized)
    name_to_cui = {}
    for cui, name in cui_to_name.items():
        key = _normalize(name)
        # Keep shortest CUI for each normalized name (most general concept)
        if key not in name_to_cui or len(cui) < len(name_to_cui[key]):
            name_to_cui[key] = cui

    logger.info(f"    Name index: {len(name_to_cui):,} unique disease names")
    logger.info(f"    Definitions: {len(cui_to_def):,}")

    return {
        "name_to_cui": name_to_cui,
        "cui_to_def": cui_to_def,
        "cui_to_name": cui_to_name,
    }


def _normalize(name: str) -> str:
    """Normalize a disease name for matching."""
    return name.lower().replace("_", " ").replace(",", "").replace("'", "").replace("-", " ").strip()


# ═══════════════════════════════════════════════════════════════════════════════
# MATCHING
# ═══════════════════════════════════════════════════════════════════════════════


def _clean_disease_name(raw: str) -> list[str]:
    """
    Clean a ClinVar disease name into searchable forms.
    Returns a list of candidate names to try, ordered by specificity.
    """
    name = raw.replace("_", " ").strip()
    if not name or name.lower() in ("not provided", ".", ""):
        return []

    candidates = [name]  # Original

    # Strip qualifier suffixes
    for suffix in [
        ", susceptibility to",
        ", susceptibility to, 1",
        ", susceptibility to, 2",
        ", susceptibility to, 3",
        ", susceptibility to, 4",
        ", susceptibility to, 5",
        ", resistance to",
        ", type 1", ", type 2", ", type 3", ", type 4",
        " 1", " 2", " 3", " 4", " 5",
    ]:
        if name.endswith(suffix):
            stripped = name[:-len(suffix)].strip()
            if stripped and stripped != name:
                candidates.append(stripped)
            break

    # Try just the first part of compound names: "Disease X, subtype Y" → "Disease X"
    if "," in name:
        first_part = name.split(",")[0].strip()
        if first_part and first_part != name:
            candidates.append(first_part)

    return candidates


def find_disease_definition(
    disease_name: str,
    db: dict,
) -> dict | None:
    """
    Find a MedGen definition for a disease name.
    Returns {name, definition, source} or None.
    """
    if not db:
        return None

    name_to_cui = db.get("name_to_cui", {})
    cui_to_def = db.get("cui_to_def", {})
    cui_to_name = db.get("cui_to_name", {})

    candidates = _clean_disease_name(disease_name)
    seen_cuis = set()

    for candidate in candidates:
        key = _normalize(candidate)
        if key in name_to_cui:
            cui = name_to_cui[key]
            if cui in seen_cuis:
                continue
            seen_cuis.add(cui)

            definition = cui_to_def.get(cui, "")
            preferred_name = cui_to_name.get(cui, candidate)

            if definition:
                if len(definition) > 350:
                    definition = definition[:347] + "..."
                return {
                    "name": preferred_name,
                    "definition": definition,
                    "medgen_cui": cui,
                    "source": "NCBI MedGen",
                }

    # Fallback: fuzzy match — check if any MedGen name contains our disease name
    norm = _normalize(disease_name)
    if len(norm) > 10:  # Only for specific-enough names
        for db_name, cui in name_to_cui.items():
            if cui in seen_cuis:
                continue
            if norm in db_name or db_name in norm:
                definition = cui_to_def.get(cui, "")
                preferred_name = cui_to_name.get(cui, db_name)
                if definition:
                    if len(definition) > 350:
                        definition = definition[:347] + "..."
                    return {
                        "name": preferred_name,
                        "definition": definition,
                        "medgen_cui": cui,
                        "source": "NCBI MedGen (fuzzy match)",
                    }

    return None


# ═══════════════════════════════════════════════════════════════════════════════
# ENRICHMENT
# ═══════════════════════════════════════════════════════════════════════════════


def enrich_clinvar_output(
    clinvar_json: str = "clinvar/clinvar_pathogenic_variants.json",
    ref_dir: str = "reference/medgen",
) -> dict:
    """
    Enrich ClinVar output with MedGen disease descriptions from local database.
    """
    logger.info("═══ MedGen Disease Description Enrichment ═══")

    json_path = Path(clinvar_json)
    if not json_path.exists():
        logger.error(f"  ClinVar JSON not found: {json_path}")
        return {"error": "ClinVar JSON not found"}

    # Check if database needs update
    status = check_medgen_update(ref_dir)
    if status.get("status") == "stale":
        logger.warning(f"  ⚠️  {status.get('message', 'Database may be outdated.')}")
        logger.warning(f"  Run: python scripts/clinical/medgen_enrich.py --download")
    elif status.get("status") == "not_downloaded":
        logger.warning("  MedGen database not downloaded.")
        logger.warning("  Run: python scripts/clinical/medgen_enrich.py --download")

    # Load MedGen database
    db = load_medgen_database(ref_dir)
    if not db:
        logger.error("  MedGen database not loaded — run with --download first")
        return {"error": "MedGen database not available"}

    # Load ClinVar data
    data = json.loads(json_path.read_text())
    variants = data.get("pathogenic_variants", [])
    if not variants:
        logger.info("  No variants to enrich")
        return {"enriched": 0}

    # Collect unique disease names
    unique_diseases = set()
    for v in variants:
        dname = v.get("disease_name", "")
        if dname and dname not in (".", "not_provided"):
            for part in dname.split("|"):
                part = part.strip()
                if part and part.lower() != "not provided":
                    unique_diseases.add(part)

    logger.info(f"  Unique disease names: {len(unique_diseases)}")

    # Build lookup cache: disease_name → description
    cache = {}
    n_found = 0
    for disease in sorted(unique_diseases):
        result = find_disease_definition(disease, db)
        if result:
            cache[disease] = result
            n_found += 1

    logger.info(f"  Matched: {n_found}/{len(unique_diseases)} diseases")

    # Enrich variants
    n_enriched = 0
    for v in variants:
        dname = v.get("disease_name", "")
        if not dname or dname in (".", "not_provided"):
            v["disease_description"] = ""
            v["disease_description_source"] = ""
            continue

        descriptions = []
        for part in dname.split("|"):
            part = part.strip()
            if part and part.lower() != "not provided":
                result = cache.get(part)
                if result and result.get("definition"):
                    descriptions.append(result["definition"])

        if descriptions:
            v["disease_description"] = " | ".join(descriptions[:3])
            v["disease_description_source"] = "NCBI MedGen"
            n_enriched += 1
        else:
            v["disease_description"] = ""
            v["disease_description_source"] = ""

    # Save
    json_path.write_text(json.dumps(data, indent=2, ensure_ascii=False))
    logger.info(f"  ✅ {n_enriched}/{len(variants)} variants enriched")
    logger.info(f"  Saved: {json_path}")

    return {
        "enriched": n_enriched,
        "total_variants": len(variants),
        "unique_diseases": len(unique_diseases),
        "matched_diseases": n_found,
    }


# ═══════════════════════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════════════════════


def main():
    import argparse

    parser = argparse.ArgumentParser(
        description="MedGen disease description enrichment for ClinVar variants"
    )
    parser.add_argument("--download", action="store_true",
                        help="Download/update MedGen RRF files from NCBI FTP")
    parser.add_argument("--force", "-f", action="store_true",
                        help="Force re-download even if cached")
    parser.add_argument("--check-update", action="store_true",
                        help="Check if local MedGen database is up to date")
    parser.add_argument("--max-age", type=int, default=30,
                        help="Max age in days before considering database stale (default: 30)")
    parser.add_argument("--ref-dir", default="reference/medgen",
                        help="MedGen reference directory")
    parser.add_argument("--clinvar-json", default="clinvar/clinvar_pathogenic_variants.json",
                        help="Path to ClinVar pathogenic variants JSON")
    parser.add_argument("--verbose", "-v", action="store_true")
    parser.add_argument("--show-samples", action="store_true",
                        help="Show sample descriptions after enrichment")
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(message)s",
    )

    if args.check_update:
        status = check_medgen_update(args.ref_dir)
        logger.info(json.dumps(status, indent=2))
        return 0 if status.get("status") in ("current", "recent") else 1

    if args.download:
        success = download_medgen(args.ref_dir, force=args.force, max_age_days=args.max_age)
        return 0 if success else 1

    result = enrich_clinvar_output(args.clinvar_json, args.ref_dir)

    if args.show_samples and result.get("enriched", 0) > 0:
        data = json.loads(Path(args.clinvar_json).read_text())
        logger.info("\n═══ Sample Descriptions ═══")
        count = 0
        for v in data.get("pathogenic_variants", []):
            desc = v.get("disease_description", "")
            if desc:
                disease = v.get("disease_name", "?")
                confidence = v.get("confidence_tier", "?")
                logger.info(f"\n  [{confidence}] {disease[:70]}")
                logger.info(f"    {desc[:150]}...")
                count += 1
                if count >= 8:
                    break

    return 0


if __name__ == "__main__":
    sys.exit(main())
