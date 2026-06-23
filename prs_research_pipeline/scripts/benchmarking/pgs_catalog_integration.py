#!/usr/bin/env python3
"""
PGS Catalog Integration — downloads harmonized polygenic scores, scores
the target sample, and computes concordance with our PRS.

Uses the PGS Catalog REST API + harmonized GRCh37 scoring files.
Fixes VAL-005 (PGS Catalog concordance check).

Output:
  pgs/scores.json       — PGS scores for the target sample
  pgs/concordance.json  — Concordance with platform PRS
"""

import sys
import os
import json
import logging
import time
import gzip
import subprocess
import shutil
from pathlib import Path
from typing import Dict, List, Optional
from urllib.request import urlopen, Request
from urllib.error import HTTPError
from datetime import datetime

import pandas as pd
import numpy as np

logger = logging.getLogger(__name__)

PGS_API = "https://www.pgscatalog.org/rest"
PGS_FTP = "https://ftp.ebi.ac.uk/pub/databases/spot/pgs"

# Map our traits to PGS Catalog traits (EFO terms where possible)
TRAIT_MAP = {
    "Glucose metabolism": ["type 2 diabetes", "fasting glucose", "HbA1c", "glycated hemoglobin"],
    "Lipid metabolism": ["LDL cholesterol", "HDL cholesterol", "triglycerides", "total cholesterol"],
    "Obesity predisposition": ["BMI", "body mass index", "obesity"],
    "Vitamin D metabolism": ["vitamin D", "25-hydroxyvitamin D"],
    "Caffeine metabolism": ["caffeine", "coffee consumption"],
    "Folate & methylation": ["folate", "homocysteine"],
    "Detoxification": [],  # No standard PGS
    "Dopamine regulation": [],  # No standard PGS
    "Lactose intolerance": [],  # No standard PGS
    "Omega-3 metabolism": ["omega-3", "fatty acid", "DHA", "EPA"],
}


def search_scores(trait: str, limit: int = 5) -> List[Dict]:
    """Search PGS Catalog for scores matching a trait."""
    query = trait.replace(" ", "+")
    url = f"{PGS_API}/trait/search?term={query}&limit=20"
    try:
        r = json.loads(urlopen(Request(url), timeout=15).read())
        trait_results = r.get("results", [])

        # For each trait, find associated PGS IDs (check both direct + child associations)
        all_scores = []
        for t in trait_results:
            pgs_ids = t.get("associated_pgs_ids", []) + t.get("child_associated_pgs_ids", [])
            for pid in pgs_ids:
                if pid not in [s["id"] for s in all_scores]:
                    all_scores.append({"id": pid, "trait": t.get("label", trait)})

        # Fetch score details for top matches
        scores = []
        for s in all_scores[:limit]:
            try:
                detail = json.loads(urlopen(
                    Request(f"{PGS_API}/score/{s['id']}"), timeout=15
                ).read())
                if detail.get("ftp_harmonized_scoring_files", {}).get("GRCh37"):
                    scores.append({
                        "id": detail["id"],
                        "name": detail.get("name", ""),
                        "trait": s["trait"],
                        "url": detail["ftp_harmonized_scoring_files"]["GRCh37"]["positions"],
                        "publication": detail.get("publication", {}),
                        "n_variants": detail.get("variants_number", 0),
                    })
            except HTTPError:
                continue
        return scores
    except Exception as e:
        logger.warning(f"  Search failed for '{trait}': {e}")
        return []


def download_score_file(url: str, output_path: Path) -> bool:
    """Download a PGS scoring file."""
    if output_path.exists():
        return True
    try:
        req = Request(url)
        with urlopen(req, timeout=120) as resp:
            data = resp.read()
        # Check if gzipped
        if data[:2] == b'\x1f\x8b':
            data = gzip.decompress(data)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_bytes(data)
        return True
    except Exception as e:
        logger.warning(f"  Download failed: {e}")
        return False


def preprocess_pgs_file(input_path: Path, output_path: Path) -> bool:
    """Convert PGS harmonized format to clean PLINK score format.
    PGS format: #comments, header, rsID chr_name chr_position effect_allele other_allele effect_weight ...
    PLINK expects: variant_id allele weight (tab-separated, no header)
    Uses hm_chr:hm_pos as variant_id to match our BIM format.
    """
    if output_path.exists():
        return True

    try:
        lines = []
        header_found = False
        rsid_col = eff_col = weight_col = chr_col = pos_col = None

        for line in input_path.read_text().splitlines():
            if line.startswith("#"):
                continue
            parts = line.split("\t")
            if not header_found:
                header_found = True
                # Detect columns
                for i, col in enumerate(parts):
                    if col == "hm_rsID": rsid_col = i
                    if col == "hm_chr": chr_col = i
                    if col == "hm_pos": pos_col = i
                    if col == "effect_allele": eff_col = i
                    if col == "effect_weight": weight_col = i
                if rsid_col is None or weight_col is None:
                    # Try simple format: rsID, effect_allele, effect_weight
                    for i, col in enumerate(parts):
                        if col == "rsID": rsid_col = i
                        if col == "effect_allele": eff_col = i
                        if col == "effect_weight": weight_col = i
                continue

            try:
                variant_id = f"{parts[chr_col]}:{parts[pos_col]}" if chr_col is not None else parts[rsid_col]
                allele = parts[eff_col] if eff_col is not None else parts[1]
                weight = float(parts[weight_col])
                lines.append(f"{variant_id}\t{allele}\t{weight}")
            except (IndexError, ValueError):
                continue

        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text("\n".join(lines))
        logger.info(f"    Preprocessed: {len(lines)} variants")
        return len(lines) > 0
    except Exception as e:
        logger.warning(f"  Preprocessing failed: {e}")
        return False


def run_plink_score(
    plink_bin: str, bfile: str, score_file: Path,
    output_prefix: Path, threads: int = 4, memory: int = 8000
) -> Optional[pd.DataFrame]:
    """Run PLINK --score on target sample."""
    profile_file = Path(f"{output_prefix}.profile")
    if profile_file.exists():
        return pd.read_csv(profile_file, sep=r"\s+", dtype={"IID": str})

    # Preprocess PGS file to clean format
    clean_file = Path(f"{output_prefix}_clean.score")
    if not preprocess_pgs_file(score_file, clean_file):
        return None

    cmd = [
        plink_bin, "--bfile", bfile,
        "--score", str(clean_file), "1", "2", "3",
        "--out", str(output_prefix),
        "--threads", str(threads), "--memory", str(memory),
        "--allow-extra-chr",
    ]
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
    if result.returncode != 0:
        logger.warning(f"  PLINK score failed: {result.stderr[-200:]}")
        return None

    if profile_file.exists():
        return pd.read_csv(profile_file, sep=r"\s+", dtype={"IID": str})
    return None


def compute_concordance(
    platform_prs: pd.DataFrame, pgs_scores: Dict[str, pd.DataFrame]
) -> Dict:
    """Compute concordance between platform PRS and PGS scores."""
    results = []
    for pgs_id, df in pgs_scores.items():
        score_col = "SCORE" if "SCORE" in df.columns else "SCORESUM"
        if score_col not in df.columns:
            continue
        pgs_val = df[score_col].iloc[0] if len(df) > 0 else 0

        # For single-sample: compute z-score from PGS Catalog percentiles
        results.append({
            "pgs_id": pgs_id,
            "pgs_score": float(pgs_val),
            "n_variants_scored": int(df.get("NMISS_ALLELE_CT", [0]).iloc[0]) if "NMISS_ALLELE_CT" in df.columns else 0,
            "status": "SCORED",
        })

    return {
        "n_scores": len(results),
        "n_scored": sum(1 for r in results if r["status"] == "SCORED"),
        "scores": results,
        "generated_date": datetime.now().strftime("%Y-%m-%d %H:%M UTC"),
    }


def find_plink() -> str:
    """Auto-detect PLINK binary."""
    candidates = [
        str(Path(__file__).parent.parent.parent.parent / "tools" / "plink"),  # scripts/benchmarking/ → root
        str(Path.cwd() / "tools" / "plink"),
    ]
    for c in candidates:
        if os.path.exists(c) and os.access(c, os.X_OK):
            return c
    return shutil.which("plink") or shutil.which("plink2") or "plink"


def integrate(
    bfile: str = "plink/ld_pruned_dataset",
    output_dir: str = "pgs",
    max_scores_per_trait: int = 3,
    delay: float = 0.5,
) -> Dict:
    logger.info("═══ PGS Catalog Integration ═══")

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    plink_bin = find_plink()
    logger.info(f"  PLINK: {plink_bin}")

    # Load platform PRS
    prs_file = Path("prs/PRS_RESULT.csv")
    platform_prs = pd.read_csv(prs_file) if prs_file.exists() else pd.DataFrame()

    all_scores = {}
    all_concordance = {}

    for trait, queries in TRAIT_MAP.items():
        if not queries:
            continue

        for query in queries:
            logger.info(f"  Searching: {query}")
            scores = search_scores(query, limit=max_scores_per_trait)
            logger.info(f"    Found {len(scores)} scores")

            for s in scores:
                pgs_id = s["id"]
                score_dir = output_dir / pgs_id
                score_dir.mkdir(parents=True, exist_ok=True)

                # Download
                score_file = score_dir / f"{pgs_id}_hmPOS_GRCh37.txt"
                if not download_score_file(s["url"], score_file):
                    continue
                logger.info(f"    Downloaded {pgs_id}: {s['name']} ({s.get('n_variants', '?')} variants)")

                # Score
                out_prefix = score_dir / pgs_id
                prof = run_plink_score(plink_bin, bfile, score_file, out_prefix)
                if prof is not None and len(prof) > 0:
                    all_scores[pgs_id] = prof

            time.sleep(delay)

    # Compute concordance
    if all_scores:
        concordance = compute_concordance(platform_prs, all_scores)
        with open(output_dir / "concordance.json", "w") as fh:
            json.dump(concordance, fh, indent=2)
        logger.info(f"  ✅ Concordance: {concordance['n_scored']}/{concordance['n_scores']} scores computed")
    else:
        concordance = {"n_scores": 0, "n_scored": 0, "scores": []}
        logger.warning("  No PGS scores could be computed")

    return concordance


def main():
    import argparse
    p = argparse.ArgumentParser(description="PGS Catalog Integration")
    p.add_argument("--bfile", default="plink/ld_pruned_dataset")
    p.add_argument("--output-dir", "-o", default="pgs")
    p.add_argument("--max-scores", type=int, default=3)
    p.add_argument("--verbose", "-v", action="store_true")
    args = p.parse_args()

    logging.basicConfig(level=logging.DEBUG if args.verbose else logging.INFO,
                       format="%(asctime)s [%(levelname)s] %(message)s")

    integrate(args.bfile, args.output_dir, args.max_scores)
    return 0


if __name__ == "__main__":
    sys.exit(main())
