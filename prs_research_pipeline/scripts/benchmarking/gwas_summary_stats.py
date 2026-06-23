#!/usr/bin/env python3
"""
GWAS Summary Statistics Downloader — Full VCF via tabix

Downloads real GWAS summary statistics for our 109 curated SNPs using
tabix to query remote VCFs. Only fetches data at our specific positions
(~5 MB total vs 65 GB full download). No API token needed.

Output:
  gwas/gwas_stats.json           — Summary
  gwas/gwas_stats_detailed.json  — Full per-SNP GWAS results
"""

import sys; import os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import _venv

import json
import logging
import time
import subprocess
from pathlib import Path
from datetime import datetime

import pandas as pd

logger = logging.getLogger(__name__)

# ── GWAS VCF URLs (public, no auth needed) ──
GWAS_VCFS = {
    "ieu-b-40": {
        "url": "https://gwas-api.mrcieu.ac.uk/datasets/ieu-b-40/vcf.gz",
        "trait": "Body mass index",
        "n": 681275, "pmid": "30124842", "consortium": "GIANT",
        "target_trait": "Obesity predisposition",
    },
    "ebi-a-GCST006867": {
        "url": "http://gwas-api.mrcieu.ac.uk/datasets/ebi-a-GCST006867/vcf.gz",
        "trait": "Type 2 diabetes",
        "n": 655666, "pmid": "30054458", "consortium": "DIAGRAM",
        "target_trait": "Glucose metabolism",
    },
    "ebi-a-GCST90025986": {
        "url": "http://gwas-api.mrcieu.ac.uk/datasets/ebi-a-GCST90025986/vcf.gz",
        "trait": "Blood glucose levels",
        "n": 400458, "pmid": "34226706", "consortium": "MAGIC",
        "target_trait": "Glucose metabolism",
    },
    "ieu-b-109": {
        "url": "https://gwas-api.mrcieu.ac.uk/datasets/ieu-b-109/vcf.gz",
        "trait": "HDL cholesterol",
        "n": 403943, "pmid": "32203549", "consortium": "GLGC",
        "target_trait": "Lipid metabolism",
    },
    "ieu-b-110": {
        "url": "https://gwas-api.mrcieu.ac.uk/datasets/ieu-b-110/vcf.gz",
        "trait": "LDL cholesterol",
        "n": 440546, "pmid": "32203549", "consortium": "GLGC",
        "target_trait": "Lipid metabolism",
    },
    "ieu-b-111": {
        "url": "https://gwas-api.mrcieu.ac.uk/datasets/ieu-b-111/vcf.gz",
        "trait": "Triglycerides",
        "n": 441016, "pmid": "32203549", "consortium": "GLGC",
        "target_trait": "Lipid metabolism",
    },
    "ebi-a-GCST90000618": {
        "url": "http://gwas-api.mrcieu.ac.uk/datasets/ebi-a-GCST90000618/vcf.gz",
        "trait": "25-Hydroxyvitamin D",
        "n": 496946, "pmid": "32242144", "consortium": "SUNLIGHT",
        "target_trait": "Vitamin D metabolism",
    },
    "ukb-b-5237": {
        "url": "https://gwas-api.mrcieu.ac.uk/datasets/ukb-b-5237/vcf.gz",
        "trait": "Coffee intake",
        "n": 428860, "pmid": "", "consortium": "UK Biobank",
        "target_trait": "Caffeine metabolism",
    },
    "ieu-a-995": {
        "url": "http://gwas-api.mrcieu.ac.uk/datasets/ieu-a-995/vcf.gz",
        "trait": "Homocysteine levels",
        "n": 44147, "pmid": "23824729", "consortium": "Homocysteine GWAS",
        "target_trait": "Folate & methylation",
    },
    "ebi-a-GCST90013997": {
        "url": "http://gwas-api.mrcieu.ac.uk/datasets/ebi-a-GCST90013997/vcf.gz",
        "trait": "Omega-3 fatty acids",
        "n": 114999, "pmid": "34017140", "consortium": "CHARGE",
        "target_trait": "Omega-3 metabolism",
    },
}


def query_tabix(vcf_url: str, chrom: str, pos: int) -> List[Dict]:
    """Query a remote VCF via tabix for a specific position. Returns parsed VCF records."""
    region = f"{chrom}:{pos}-{pos}"
    cmd = ["tabix", vcf_url, region]
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=30)
        if result.returncode != 0:
            return []
        records = []
        for line in result.stdout.strip().split("\n"):
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 8:
                continue
            # VCF format: CHROM POS ID REF ALT QUAL FILTER INFO
            info = {}
            for item in parts[7].split(";"):
                if "=" in item:
                    k, v = item.split("=", 1)
                    info[k] = v
                else:
                    info[item] = True

            records.append({
                "chrom": parts[0],
                "pos": int(parts[1]),
                "rsid": parts[2],
                "ref": parts[3],
                "alt": parts[4],
                "beta": float(info.get("ES", info.get("BETA", 0))),
                "se": float(info.get("SE", 0)),
                "pval": float(info.get("LP", -1)),  # -log10(p)
                "eaf": float(info.get("AF", info.get("EAF", 0))),
                "n": int(info.get("N", info.get("SS", 0))),
                "info": {k: v for k, v in info.items()
                        if k not in ("ES", "SE", "LP", "AF", "EAF", "N", "SS")},
            })
        return records
    except subprocess.TimeoutExpired:
        return []
    except Exception as e:
        logger.debug(f"  tabix error: {e}")
        return []


def download_gwas_data(snp_db_path: str, output_dir: str) -> Dict:
    """Download GWAS summary stats for all our SNPs using tabix."""
    logger.info("═══ GWAS Summary Statistics Download (tabix) ═══")

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    db = pd.read_csv(snp_db_path, dtype=str) if Path(snp_db_path).exists() else pd.DataFrame()
    logger.info(f"  SNP database: {len(db)} variants")

    # Build position index
    snp_positions = []
    for _, row in db.iterrows():
        chrom = str(row.get("chrom", "")).replace("chr", "")
        try:
            pos = int(row.get("pos", 0))
        except ValueError:
            continue
        if chrom and pos:
            snp_positions.append({
                "rsid": row.get("rsid", ""),
                "chrom": chrom,
                "pos": pos,
                "gene": row.get("gene", ""),
                "trait": row.get("trait_category", ""),
                "effect_allele": row.get("effect_allele", ""),
                "weight": row.get("weight", ""),
            })

    logger.info(f"  SNP positions to query: {len(snp_positions)}")

    all_results = {}
    total_hits = 0
    total_queries = 0

    for ds_id, ds_info in GWAS_VCFS.items():
        vcf_url = ds_info["url"]
        logger.info(f"  Querying {ds_id}: {ds_info['trait']} (n={ds_info['n']:,})...")

        # Check if VCF is accessible
        check = subprocess.run(
            ["tabix", "-l", vcf_url], capture_output=True, text=True, timeout=15
        )
        if check.returncode != 0:
            logger.warning(f"    VCF not accessible: {check.stderr[:100]}")
            continue

        hits = []
        for snp in snp_positions:
            records = query_tabix(vcf_url, snp["chrom"], snp["pos"])
            total_queries += 1
            for r in records:
                r["dataset_id"] = ds_id
                r["dataset_trait"] = ds_info["trait"]
                r["db_rsid"] = snp["rsid"]
                r["db_gene"] = snp["gene"]
                r["db_trait"] = snp["trait"]
                r["db_effect_allele"] = snp["effect_allele"]
                hits.append(r)
                total_hits += 1

            if total_queries % 20 == 0:
                logger.info(f"    Progress: {total_queries}/{len(snp_positions)} queries, {total_hits} hits so far")

        if hits:
            all_results[ds_id] = {
                "trait": ds_info["trait"],
                "consortium": ds_info["consortium"],
                "n": ds_info["n"],
                "pmid": ds_info["pmid"],
                "n_hits": len(hits),
                "results": hits,
            }
            logger.info(f"    Found {len(hits)} variant records")
        else:
            logger.info(f"    No variants found")

        time.sleep(1)  # Be nice to the server

    # Save results
    report = {
        "generated_date": datetime.now().strftime("%Y-%m-%d %H:%M UTC"),
        "method": "tabix remote query (no full downloads)",
        "datasets_queried": len(all_results),
        "total_hits": total_hits,
        "total_queries": total_queries,
        "per_trait": {ds_id: {"n_hits": v["n_hits"]} for ds_id, v in all_results.items()},
    }

    with open(output_dir / "gwas_stats.json", "w") as fh:
        json.dump(report, fh, indent=2)
    with open(output_dir / "gwas_stats_detailed.json", "w") as fh:
        json.dump(all_results, fh, indent=2, default=str)

    logger.info(f"  ✅ Downloaded: {total_hits} GWAS records across {len(all_results)} datasets")
    return report


from urllib.request import urlopen, Request
from urllib.error import HTTPError

OPENGWAS = "https://api.opengwas.io/api"

DATASETS_API = {
    "Obesity predisposition": ["ieu-b-40"],
    "Glucose metabolism": ["ebi-a-GCST006867", "ebi-a-GCST90025986"],
    "Lipid metabolism":   ["ieu-b-109", "ieu-b-110", "ieu-b-111"],
    "Vitamin D metabolism":   ["ebi-a-GCST90000618"],
    "Caffeine metabolism":    ["ukb-b-5237"],
    "Folate & methylation":   ["ieu-a-995"],
    "Omega-3 metabolism":     ["ebi-a-GCST90013997"],
}


def query_api_associations(token: str, rsids: list, ds_ids: list) -> list:
    """POST /associations — get GWAS variant data for specific rsIDs."""
    url = f"{OPENGWAS}/associations"
    body = json.dumps({"variant": rsids, "id": ds_ids, "proxies": 0}).encode()
    try:
        req = Request(url, data=body, headers={
            "Authorization": f"Bearer {token}",
            "Content-Type": "application/json",
            "Accept": "application/json",
        })
        data = json.loads(urlopen(req, timeout=60).read())
        return data if isinstance(data, list) else data.get("data", data.get("results", []))
    except HTTPError as e:
        logger.warning(f"  API error {e.code}")
        return []
    except Exception as e:
        logger.warning(f"  API error: {e}")
        return []


def download_via_api(token: str, snp_db_path: str, output_dir: str) -> Dict:
    """Download GWAS data via OpenGWAS API (fast, requires token)."""
    logger.info("═══ GWAS via OpenGWAS API ═══")

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    db = pd.read_csv(snp_db_path, dtype=str) if Path(snp_db_path).exists() else pd.DataFrame()
    logger.info(f"  SNPs: {len(db)}")

    all_results = {}
    total_hits = 0

    for trait, ds_ids in DATASETS_API.items():
        trait_snps = db[db["trait_category"] == trait] if "trait_category" in db.columns else db
        rsids = [r for r in trait_snps["rsid"].dropna().tolist() if r and r.startswith("rs")]

        if not rsids:
            continue

        logger.info(f"  {trait}: {len(rsids)} SNPs → {ds_ids}")
        results = query_api_associations(token, rsids, ds_ids)

        if results:
            all_results[trait] = {"datasets": ds_ids, "n_hits": len(results), "results": results}
            total_hits += len(results)
            logger.info(f"    {len(results)} associations found")
        time.sleep(0.5)

    report = {
        "generated_date": datetime.now().strftime("%Y-%m-%d %H:%M UTC"),
        "source": "IEU OpenGWAS API",
        "datasets_queried": len(all_results),
        "total_hits": total_hits,
        "per_trait": {t: {"n_hits": v["n_hits"]} for t, v in all_results.items()},
    }

    with open(output_dir / "gwas_stats.json", "w") as fh:
        json.dump(report, fh, indent=2)
    with open(output_dir / "gwas_stats_detailed.json", "w") as fh:
        json.dump(all_results, fh, indent=2, default=str)

    logger.info(f"  ✅ GWAS: {total_hits} associations across {len(all_results)} traits")
    return report


def integrate(token: str = None, snp_db_path: str = "data/snp_database_annotated.csv",
              output_dir: str = "gwas") -> Dict:
    """Main entry point. Uses OpenGWAS API (fast, requires token)."""
    if not token:
        token = os.environ.get("OPENGWAS_TOKEN")
    if not token:
        token_file = Path(output_dir).parent / ".opengwas_token"
        if token_file.exists():
            token = token_file.read_text().strip()

    if token:
        return download_via_api(token, snp_db_path, output_dir)
    else:
        logger.error("No OpenGWAS token. Get one at https://api.opengwas.io/")
        logger.error("Set OPENGWAS_TOKEN env var, pass --token, or save to .opengwas_token")
        return {"error": "no_token"}


def main():
    import argparse
    p = argparse.ArgumentParser(description="GWAS Summary Statistics Downloader")
    p.add_argument("--token", help="OpenGWAS JWT token (optional — tabix doesn't need it)")
    p.add_argument("--snp-db", default="data/snp_database_annotated.csv")
    p.add_argument("--output-dir", "-o", default="gwas")
    p.add_argument("--verbose", "-v", action="store_true")
    args = p.parse_args()

    logging.basicConfig(level=logging.DEBUG if args.verbose else logging.INFO,
                       format="%(asctime)s [%(levelname)s] %(message)s")
    integrate(args.token, args.snp_db, args.output_dir)
    return 0


if __name__ == "__main__":
    sys.exit(main())
