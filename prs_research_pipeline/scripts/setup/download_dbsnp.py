#!/usr/bin/env python3
"""
dbSNP Annotation Enrichment — queries NCBI Entrez API for each rsID.
No large file downloads needed. Enriches SNP database with clinical
significance, gene names, global MAF, and dbSNP validation status.

Output:
  data/snp_database_dbsnp_enriched.csv
  data/dbsnp_annotations.json  (API cache — subsequent runs instant)
"""

import sys
import os
import json
import logging
import time
from pathlib import Path
from urllib.request import urlopen, Request
from urllib.error import HTTPError, URLError
from datetime import datetime

import pandas as pd

logger = logging.getLogger(__name__)

ENTREZ_EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"


def fetch_snp(rsid: str, retries: int = 3) -> dict | None:
    """Fetch SNP annotation from NCBI Entrez E-utilities."""
    for attempt in range(retries):
        try:
            search_url = f"{ENTREZ_EUTILS}/esearch.fcgi?db=snp&term={rsid}[SNP]&retmode=json"
            with urlopen(Request(search_url), timeout=15) as resp:
                search_data = json.loads(resp.read().decode())

            uid_list = search_data.get("esearchresult", {}).get("idlist", [])
            if not uid_list:
                return None

            # Fetch summaries for all matching UIDs, pick the one matching our rsID
            summary_url = f"{ENTREZ_EUTILS}/esummary.fcgi?db=snp&id={','.join(uid_list)}&retmode=json"
            with urlopen(Request(summary_url), timeout=15) as resp2:
                data = json.loads(resp2.read().decode())

            result_data = data.get("result", {})
            for uid in uid_list:
                entry = result_data.get(uid, {})
                if entry.get("snp_id") == rsid:
                    return entry
            # Fallback: return first result
            return result_data.get(uid_list[0], {})

        except (HTTPError, URLError, json.JSONDecodeError):
            if attempt < retries - 1:
                time.sleep(1.5 ** attempt)
                continue
            return None
        except Exception as e:
            if attempt < retries - 1:
                time.sleep(1.5 ** attempt)
                continue
            logger.warning(f"  {rsid}: {e}")
            return None
    return None


def fetch_clinvar(rsid: str) -> str:
    """Try ClinVar API for clinical significance when SNP API doesn't have it."""
    try:
        url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=clinvar&term={rsid}&retmode=json"
        r = json.loads(urlopen(Request(url), timeout=10).read())
        ids = r.get("esearchresult", {}).get("idlist", [])
        if ids:
            furl = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=clinvar&id={ids[0]}&retmode=json"
            cr = json.loads(urlopen(Request(furl), timeout=10).read())
            return cr.get("result", {}).get(ids[0], {}).get("clinical_significance", "")
    except Exception:
        pass
    return ""


def extract_annotations(data: dict) -> dict:
    """Extract key annotations from Entrez esummary response."""
    return {
        "rsid": data.get("snp_id", ""),
        "uid": data.get("uid", 0),
        "chr_pos": f"{data.get('chr','')}:{data.get('chrpos','')}" if data.get("chr") else "",
        "clinical_significance": data.get("clinical_significance", ""),
        "gene_name": (data.get("genes", [{}]) or [{}])[0].get("name", ""),
        "variant_type": data.get("snp_class", ""),
        "global_maf": data.get("global_maf", ""),
        "validated": data.get("validated", "") == "by-cluster,by-frequency",
    }


def enrich_database(
    snp_db: str = "data/snp_database_annotated.csv",
    output_dir: str = "data",
    delay: float = 0.5,
) -> dict:
    logger.info("═══ dbSNP Annotation Enrichment ═══")

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    db = pd.read_csv(snp_db, dtype=str)
    rsids = sorted(db["rsid"].dropna().unique())
    logger.info(f"  SNPs to query: {len(rsids)}")

    # Load cache
    cache_path = output_dir / "dbsnp_annotations.json"
    cache = json.loads(cache_path.read_text()) if cache_path.exists() else {}
    logger.info(f"  Cache: {len(cache)} previously fetched")

    # Fetch missing
    annotations = {}
    n_new = 0
    for i, rsid in enumerate(rsids):
        if rsid in cache:
            annotations[rsid] = cache[rsid]
            continue
        if (i + 1) % 10 == 0:
            logger.info(f"  Progress: {i+1}/{len(rsids)}...")
        data = fetch_snp(rsid)
        ann = extract_annotations(data) if data else {"rsid": rsid, "status": "not_found"}
        # ClinVar fallback for clinical significance
        if data and not ann.get("clinical_significance"):
            clinvar_sig = fetch_clinvar(rsid)
            if clinvar_sig:
                ann["clinical_significance"] = clinvar_sig
        annotations[rsid] = ann
        n_new += 1
        time.sleep(delay)

    # Save cache
    cache_path.write_text(json.dumps(annotations, indent=2))
    logger.info(f"  Newly fetched: {n_new} | Total: {len(annotations)}")

    # Enrich database
    enrichments = [
        {
            "rsid": rsid,
            "dbsnp_status": ann.get("status", "fetched"),
            "dbsnp_chr_pos": ann.get("chr_pos", ""),
            "clinical_significance": ann.get("clinical_significance", ""),
            "dbsnp_gene": ann.get("gene_name", ""),
            "variant_type": ann.get("variant_type", ""),
            "global_maf": ann.get("global_maf", ""),
            "dbsnp_validated": ann.get("validated", False),
        }
        for rsid, ann in annotations.items()
    ]

    enrich_df = pd.DataFrame(enrichments)
    enriched = pd.concat([db, enrich_df.drop(columns=["rsid"])], axis=1)
    enriched.to_csv(output_dir / "snp_database_dbsnp_enriched.csv", index=False)

    n_clin = int(enrich_df["clinical_significance"].astype(bool).sum())
    n_valid = int(enrich_df["dbsnp_validated"].sum())
    logger.info(f"  ✅ Enriched: {output_dir}/snp_database_dbsnp_enriched.csv")
    logger.info(f"  Results: {n_clin} clinical, {n_valid} validated")

    return {"total": len(rsids), "clinical": n_clin, "validated": n_valid}


def main():
    import argparse
    p = argparse.ArgumentParser(description="dbSNP annotation enrichment")
    p.add_argument("--snp-db", default="data/snp_database_annotated.csv")
    p.add_argument("--output-dir", "-o", default="data")
    p.add_argument("--delay", type=float, default=0.5)
    p.add_argument("--verbose", "-v", action="store_true")
    args = p.parse_args()

    logging.basicConfig(level=logging.DEBUG if args.verbose else logging.INFO,
                       format="%(asctime)s [%(levelname)s] %(message)s")
    enrich_database(args.snp_db, args.output_dir, args.delay)
    return 0


if __name__ == "__main__":
    sys.exit(main())
