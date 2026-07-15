#!/usr/bin/env python3
"""
╔══════════════════════════════════════════════════════════════════════════════╗
║   PHARMACOGENOMIC VARIANT ANNOTATOR (PharmGKB / CPIC)                       ║
║   scripts/clinical/pharmgkb_annotator.py                                     ║
║                                                                            ║
║   Annotates user WGS VCF against a curated database of pharmacogenomic      ║
║   variants from CPIC (Clinical Pharmacogenetics Implementation Consortium)  ║
║   guidelines and PharmGKB.                                                  ║
║                                                                            ║
║   Approach:                                                                 ║
║     1. Load curated PharmGKB variant database (data/pharmgkb_variants.csv)  ║
║     2. For each variant, query user VCF via tabix                           ║
║     3. Determine genotype and copies of effect allele                       ║
║     4. Map to phenotype and drug recommendation                             ║
║     5. Output structured JSON + summary                                     ║
║                                                                            ║
║   Output:                                                                   ║
║     pharmgkb/pharmgkb_drug_report.json                                      ║
╚══════════════════════════════════════════════════════════════════════════════╝
"""

import sys
import os
import json
import logging
import subprocess
import csv
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from collections import defaultdict, Counter
from datetime import datetime, timezone

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from utils.constants import PIPELINE_VERSION

logger = logging.getLogger(__name__)


def load_pharmgkb_db(csv_path: str) -> List[Dict]:
    """Load the curated pharmacogenomic variant database."""
    variants = []
    with open(csv_path, "r") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            row["pos"] = int(row["pos"].strip())
            row["chrom"] = row["chrom"].strip()
            row["ref"] = row["ref"].strip().upper()
            row["alt"] = row["alt"].strip().upper()
            row["effect_allele"] = row["effect_allele"].strip()
            row["rsid"] = row["rsid"].strip()
            variants.append(row)
    return variants


def query_vcf_position(vcf_path: str, chrom: str, pos: int) -> Optional[Dict]:
    """
    Query a user VCF for a specific position via tabix.
    Returns parsed variant dict or None.
    Handles chr prefix mismatch and position range (±5bp).
    """
    # Try with chr prefix (most VCFs use this)
    chrom_str = f"chr{chrom}" if not chrom.startswith("chr") else chrom
    region = f"{chrom_str}:{pos}-{pos}"
    cmd = ["tabix", vcf_path, region]

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=30)
        lines = [l.strip() for l in result.stdout.strip().split("\n") if l.strip()]
        if not lines:
            return None

        # Take the first matching record
        parts = lines[0].split("\t")
        if len(parts) < 10:
            return None

        ref = parts[3].upper()
        alt = parts[4].upper()
        fmt = parts[8].split(":")
        sample = parts[9].split(":")

        # Parse genotype (GT field)
        gt = None
        if "GT" in fmt:
            gt_idx = fmt.index("GT")
            gt = sample[gt_idx] if gt_idx < len(sample) else None

        return {
            "chrom": parts[0],
            "pos": int(parts[1]),
            "ref": ref,
            "alt": alt,
            "gt": gt,
            "qual": parts[5],
            "filter": parts[6],
        }

    except subprocess.TimeoutExpired:
        logger.warning(f"  tabix timeout for {chrom}:{pos}")
        return None
    except Exception as e:
        logger.warning(f"  tabix error for {chrom}:{pos}: {e}")
        return None


def count_effect_alleles(gt: str, ref: str, alt: str, effect_allele: str) -> int:
    """
    Count copies of the effect allele from genotype string.
    GT format: "0/0", "0/1", "1/1", "0|1", etc.
    """
    if not gt or gt in ("./.", ".|.", "./.", "NA", ""):
        return -1  # Unknown

    alleles = [ref, alt]  # 0 = ref, 1 = alt for simple biallelic

    # Handle multi-allelic: effect allele might be any of the ALT alleles
    alts = alt.split(",")
    all_alleles = [ref] + alts

    count = 0
    for allele_idx in gt.replace("|", "/").split("/"):
        try:
            idx = int(allele_idx)
            if idx < len(all_alleles):
                allele = all_alleles[idx]
                if allele.upper() == effect_allele.upper():
                    count += 1
        except (ValueError, IndexError):
            return -1

    return count


def infer_phenotype(count: int, star_phenotype: str) -> str:
    """
    Infer metabolic phenotype from allele count and variant phenotype.
    Maps copy number to clinical phenotype.
    """
    if count < 0:
        return "Unknown (no genotype data)"

    if "non-functional" in star_phenotype.lower() or "deficient" in star_phenotype.lower():
        if count == 2:
            return "Homozygous deficient — severe impact"
        elif count == 1:
            return "Heterozygous — reduced function"
        else:
            return "Normal function"

    if "reduced function" in star_phenotype.lower() or "reduced expression" in star_phenotype.lower():
        if count == 2:
            return "Homozygous — significantly reduced"
        elif count == 1:
            return "Heterozygous — mildly reduced"
        else:
            return "Normal function"

    if "poor metabolizer" in star_phenotype.lower():
        if count == 2:
            return "Poor metabolizer"
        elif count == 1:
            return "Intermediate metabolizer"
        else:
            return "Normal metabolizer"

    if "rapid" in star_phenotype.lower() or "ultrarapid" in star_phenotype.lower():
        if count >= 1:
            return "Rapid/Ultrarapid metabolizer"
        else:
            return "Normal metabolizer"

    if "altered function" in star_phenotype.lower():
        if count >= 1:
            return "Altered function"
        else:
            return "Normal function"

    if "unfavorable" in star_phenotype.lower():
        if count >= 1:
            return "Unfavorable genotype"
        else:
            return "Favorable genotype"

    # Default
    if count >= 1:
        return f"Variant carrier ({count} copies)"
    return "Normal"


def classify_actionability(count: int, star_phenotype: str, cpic_level: str) -> str:
    """
    Classify how actionable a finding is.
    Returns 'critical', 'important', 'informative', or 'normal'.
    """
    if count <= 0:
        return "normal"

    pheno_lower = star_phenotype.lower()

    if "non-functional" in pheno_lower or "deficient" in pheno_lower:
        if count == 2:
            return "critical"
        return "important"

    if "poor metabolizer" in pheno_lower:
        if count == 2:
            return "important"
        return "informative"

    if "rapid" in pheno_lower:
        return "informative"

    if "reduced" in pheno_lower:
        return "important" if count == 2 else "informative"

    if cpic_level == "A":
        return "informative" if count >= 1 else "normal"

    return "informative"


def enrich_db_from_clinpgx(
    pharmgkb_db: List[Dict],
    pipeline_root: Path,
) -> List[Dict]:
    """
    Enrich the variant database with ClinPGx clinical variants.
    Uses clinicalVariants TSV for variant-drug pairs and variants TSV for positions.
    Returns enriched database (original + new variants from ClinPGx).
    """
    clinpgx_path = pipeline_root / "reference" / "clinpgx" / "clinpgx_parsed.json"
    if not clinpgx_path.exists():
        logger.info("  ClinPGx data not found — using curated CSV only")
        return pharmgkb_db

    try:
        clinpgx = json.loads(clinpgx_path.read_text())
    except Exception:
        return pharmgkb_db

    clinical_variants = clinpgx.get("clinical_variants", [])
    variants = clinpgx.get("variants", {})

    if not clinical_variants:
        return pharmgkb_db

    # Filter to high-evidence rsID-only variants not already in the database
    existing_rsids = {v["rsid"] for v in pharmgkb_db if v.get("rsid")}
    high_evidence = {"1A", "1B"}

    n_added = 0
    for cv in clinical_variants:
        rsid = cv.get("variant", "")
        if not rsid or not rsid.startswith("rs") or rsid in existing_rsids:
            continue
        if cv.get("evidence_level") not in high_evidence:
            continue

        # Get position info from variants database
        var_info = variants.get(rsid, {})
        location = var_info.get("location", "")
        gene = cv.get("gene", "") or var_info.get("gene", "")

        # Parse NCBI location format: "NC_000003.12:183917980" → chr:pos
        chrom = ""
        pos = ""
        if location:
            parts = location.split(":")
            if len(parts) >= 2:
                chr_part = parts[0].replace("NC_", "")
                # Convert NCBI accession to chromosome number
                chrom_match = chr_part.split(".")
                if chrom_match:
                    chrom = chrom_match[0]

        if not chrom or not pos:
            continue

        pharmgkb_db.append({
            "gene": gene,
            "rsid": rsid,
            "chrom": chrom,
            "pos": int(pos) if pos else 0,
            "ref": "",
            "alt": "",
            "effect_allele": "",
            "star_allele": "",
            "star_phenotype": cv.get("type", ""),
            "drug": cv.get("drug", ""),
            "drug_class": "",
            "recommendation_en": "",
            "recommendation_es": "",
            "cpic_level": cv.get("evidence_level", "").replace("A", "A").replace("B", "B")[:1],
            "pmid": "",
        })
        existing_rsids.add(rsid)
        n_added += 1

    if n_added > 0:
        logger.info(f"  ClinPGx enriched: +{n_added} high-evidence variants ({len(pharmgkb_db)} total)")
    return pharmgkb_db


def annotate_user_vcf(
    user_vcf: str,
    pharmgkb_db: List[Dict],
    output_dir: str = "pharmgkb",
) -> Dict:
    """
    Annotate user VCF against pharmacogenomic variant database.
    """
    logger.info("═══ Pharmacogenomic Variant Annotation ═══")
    logger.info(f"  Database: {len(pharmgkb_db)} variant-drug associations")

    # Group by position — each position may affect MULTIPLE drugs
    pos_to_variants = {}
    for v in pharmgkb_db:
        key = (v["chrom"], v["pos"], v["ref"], v["alt"])
        if key not in pos_to_variants:
            pos_to_variants[key] = []
        pos_to_variants[key].append(v)

    logger.info(f"  Unique positions to query: {len(pos_to_variants)}")

    # Query each position, report ALL drugs affected
    results = []
    n_found = 0
    n_total = len(pos_to_variants)
    for i, ((chrom, pos, ref, alt), variants_at_pos) in enumerate(pos_to_variants.items()):
        if (i + 1) % 10 == 0:
            logger.info(f"  Progress: {i+1}/{n_total}...")

        vcf_record = query_vcf_position(user_vcf, chrom, pos)

        if vcf_record and vcf_record.get("gt"):
            # For each drug at this position, check the effect allele
            for variant in variants_at_pos:
                effect_allele = variant["effect_allele"]
                if not effect_allele:
                    continue
                count = count_effect_alleles(vcf_record["gt"], ref, alt, effect_allele)

                if count > 0:
                    n_found += 1
                    phenotype = infer_phenotype(count, variant["star_phenotype"])
                    actionability = classify_actionability(count, variant["star_phenotype"], variant["cpic_level"])

                    results.append({
                        "gene": variant["gene"],
                        "rsid": variant["rsid"],
                        "star_allele": variant["star_allele"],
                        "chrom": str(chrom),
                        "pos": pos,
                        "ref": ref,
                        "alt": alt,
                        "effect_allele": effect_allele,
                        "genotype": vcf_record["gt"],
                        "copies": count,
                        "phenotype": phenotype,
                        "star_phenotype": variant["star_phenotype"],
                        "drug": variant["drug"],
                        "drug_class": variant["drug_class"],
                        "recommendation_en": variant["recommendation_en"],
                        "recommendation_es": variant["recommendation_es"],
                        "cpic_level": variant["cpic_level"],
                        "actionability": actionability,
                        "pmid": variant["pmid"],
                    })

    logger.info(f"  ✓ {n_found} actionable pharmacogenomic findings")

    # Group by gene
    by_gene = defaultdict(list)
    for r in results:
        by_gene[r["gene"]].append(r)

    # Group by drug
    by_drug = defaultdict(list)
    for r in results:
        by_drug[r["drug"]].append(r)

    # Group by actionability
    by_actionability = Counter(r["actionability"] for r in results)

    # Summary
    genes_with_hits = sorted(by_gene.keys())
    drugs_affected = sorted(by_drug.keys())

    output = {
        "metadata": {
            "analysis_date": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
            "pipeline_version": PIPELINE_VERSION,
            "user_vcf": user_vcf,
            "pharmgkb_db": "data/pharmgkb_variants.csv",
            "cpic_guidelines": "CPIC Level A/B (Clinical Pharmacogenetics Implementation Consortium). Run clinpgx_sync.py --sync for full guideline texts.",
            "variants_queried": len(pos_to_variants),
            "findings": n_found,
        },
        "pharmacogenomic_findings": sorted(results, key=lambda r: (
            {"critical": 0, "important": 1, "informative": 2, "normal": 3}[r["actionability"]],
            r["gene"],
            r["drug"],
        )),
        "summary": {
            "total_findings": n_found,
            "by_actionability": dict(by_actionability),
            "genes_with_findings": genes_with_hits,
            "drugs_affected": drugs_affected,
            "by_gene": {g: len(by_gene[g]) for g in genes_with_hits},
            "by_drug": {d: len(by_drug[d]) for d in drugs_affected},
        },
    }

    return output


def main():
    import argparse

    parser = argparse.ArgumentParser(
        description="Pharmacogenomic variant annotation from curated CPIC/PharmGKB database"
    )
    parser.add_argument("--vcf", required=True, help="User VCF path (.vcf.gz)")
    parser.add_argument("--pharmgkb-db", default="data/pharmgkb_variants.csv",
                        help="PharmGKB variant database CSV")
    parser.add_argument("--output-dir", "-o", default="pharmgkb",
                        help="Output directory")
    parser.add_argument("--verbose", "-v", action="store_true")
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(message)s",
    )

    # Resolve paths
    pipeline_root = Path(__file__).resolve().parent.parent.parent
    output_dir = pipeline_root / args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    db_path = args.pharmgkb_db
    if not Path(db_path).is_absolute():
        db_path = str(pipeline_root / db_path)

    if not Path(db_path).exists():
        logger.error(f"PharmGKB database not found: {db_path}")
        return 1

    user_vcf = args.vcf if Path(args.vcf).is_absolute() else str(Path.cwd() / args.vcf)

    if not Path(user_vcf).exists():
        logger.error(f"User VCF not found: {user_vcf}")
        return 1

    # Load database
    db = load_pharmgkb_db(db_path)

    # Enrich with ClinPGx clinical variants if available
    db = enrich_db_from_clinpgx(db, pipeline_root)

    # Annotate
    output = annotate_user_vcf(user_vcf, db, str(output_dir))

    # Enrich with ClinPGx guideline texts if available
    clinpgx_path = Path(output_dir).parent / "reference" / "clinpgx" / "clinpgx_parsed.json"
    if not clinpgx_path.exists():
        clinpgx_path = Path(output_dir) / "clinpgx_parsed.json"  # Fallback
    if clinpgx_path.exists():
        logger.info("  Enriching with ClinPGx guideline texts...")
        try:
            clinpgx = json.loads(clinpgx_path.read_text())
            guidelines = clinpgx.get("guidelines", [])
            # Build drug+gene → guideline lookup
            guideline_lookup = {}
            for g in guidelines:
                for drug in g.get("drugs", []):
                    for gene in g.get("genes", []):
                        key = (gene.upper(), drug.lower())
                        guideline_lookup[key] = g

            for finding in output["pharmacogenomic_findings"]:
                key = (finding["gene"].upper(), finding["drug"].lower())
                guide = guideline_lookup.get(key)
                if guide:
                    finding["guideline_name"] = guide.get("name", "")
                    finding["guideline_source"] = guide.get("source", "")
                    finding["guideline_summary"] = guide.get("summary_markdown", "")
                    finding["guideline_text"] = guide.get("text_markdown", "")

            output["metadata"]["clinpgx_enriched"] = True
            output["metadata"]["clinpgx_guidelines"] = len(guidelines)
            logger.info(f"    ✓ Enriched with {len(guidelines)} ClinPGx guidelines")
        except Exception as e:
            logger.warning(f"    Could not enrich: {e}")

    # Write output
    output_path = output_dir / "pharmgkb_drug_report.json"
    output_path.write_text(json.dumps(output, indent=2, ensure_ascii=False))
    logger.info(f"  ✅ {output_path}")

    # Log summary
    summary = output["summary"]
    logger.info(f"\n═══ Pharmacogenomic Summary ═══")
    logger.info(f"  Findings: {summary['total_findings']}")
    logger.info(f"  Critical: {summary['by_actionability'].get('critical', 0)}")
    logger.info(f"  Important: {summary['by_actionability'].get('important', 0)}")
    logger.info(f"  Informative: {summary['by_actionability'].get('informative', 0)}")
    if summary['genes_with_findings']:
        logger.info(f"  Genes: {', '.join(summary['genes_with_findings'])}")
    if summary['drugs_affected']:
        logger.info(f"  Drugs: {', '.join(summary['drugs_affected'])}")

    if summary['total_findings'] > 0:
        logger.info(f"\n  Top findings:")
        for f in output['pharmacogenomic_findings'][:8]:
            icon = {"critical": "🔴", "important": "🟠", "informative": "🟡"}.get(f['actionability'], "⚪")
            logger.info(f"  {icon} {f['gene']:8s} {f['star_allele']:6s} → {f['drug']:20s} ({f['copies']} copies)")
            logger.info(f"     {f['recommendation_en'][:120]}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
