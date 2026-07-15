#!/usr/bin/env python3
"""
╔══════════════════════════════════════════════════════════════════════════════╗
║   CLINVAR PATHOGENIC VARIANT ANNOTATOR                                       ║
║   scripts/clinical/clinvar_annotator.py                                      ║
║                                                                            ║
║   Annotates every variant in a user WGS VCF against the ClinVar germline    ║
║   classification database (GRCh37). Finds variants classified as            ║
║   Pathogenic, Likely_pathogenic, or Pathogenic/Likely_pathogenic.           ║
║                                                                            ║
║   Algorithm (two-pass, position-first):                                     ║
║     1. bcftools query → extract CHROM:POS from user VCF as regions BED      ║
║     2. tabix -R → query ClinVar for all positionally overlapping records    ║
║     3. Python → exact REF/ALT match + CLNSIG classification                 ║
║     4. Filter pathogenic → build summary → write JSON                       ║
║                                                                            ║
║   CLNSIG values considered pathogenic:                                      ║
║     • Pathogenic                                                            ║
║     • Likely_pathogenic                                                     ║
║     • Pathogenic/Likely_pathogenic                                          ║
║     • Pathogenic,_low_penetrance                                            ║
║     • Likely_pathogenic,_low_penetrance                                     ║
║     • established_risk_allele                                               ║
║     • Likely_risk_allele                                                    ║
║     • risk_factor                                                           ║
║     • Affects                                                               ║
║                                                                            ║
║   Excluded: Uncertain_significance, Benign, Likely_benign, drug_response    ║
║             (unless co-occurring with pathogenic), not_provided, other      ║
║                                                                            ║
║   Output:                                                                   ║
║     clinvar/clinvar_pathogenic_variants.json                                ║
╚══════════════════════════════════════════════════════════════════════════════╝
"""

import sys
import os
import json
import logging
import subprocess
import tempfile
import gzip
from pathlib import Path
from datetime import datetime, timezone
from typing import Dict, List, Optional, Set, Tuple
from collections import Counter, defaultdict

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from utils.constants import PIPELINE_VERSION

logger = logging.getLogger(__name__)

# ── Pathogenic CLNSIG values ───────────────────────────────────────────────────

# Primary pathogenic categories
PATHOGENIC_DIRECT = {
    "pathogenic",
    "likely_pathogenic",
    "pathogenic/likely_pathogenic",
}

# Pathogenic with qualifiers
PATHOGENIC_QUALIFIED = {
    "pathogenic,_low_penetrance",
    "likely_pathogenic,_low_penetrance",
}

# Risk alleles / established associations
PATHOGENIC_RISK = {
    "established_risk_allele",
    "likely_risk_allele",
    "risk_factor",
    "affects",
    "protective",       # protective allele has clinical significance
}

# All values considered "pathogenic" for filtering
PATHOGENIC_ALL = PATHOGENIC_DIRECT | PATHOGENIC_QUALIFIED | PATHOGENIC_RISK

# Explicitly excluded values (even if they appear alongside pathogenic)
EXCLUDED_CLNSIG = {
    "benign",
    "likely_benign",
    "benign/likely_benign",
    "uncertain_significance",
    "uncertain_significance/uncertain_risk_allele",
    "not_provided",
    "no_classification_for_the_single_variant",
    "no_classifications_from_unflagged_records",
    "other",
    "drug_response",           # pharmacogenomic, not disease-causing
    "conflicting_classifications_of_pathogenicity",  # too uncertain
    "association",             # too vague
    "vus-high",                # still uncertain significance
    "uncertain_risk_allele",
}


def _strip_chr(chrom: str) -> str:
    """Normalize chromosome name: strip 'chr' prefix for consistent matching."""
    return chrom.replace("chr", "").replace("Chr", "").replace("CHR", "")


def _clean_value(val: str) -> str:
    """Normalize ClinVar INFO value: lowercase, strip whitespace, remove quotes."""
    return val.strip().strip('"').lower()


def _split_clnsig(clnsig_raw: str) -> List[str]:
    """
    Split CLNSIG field into individual classification values.
    ClinVar uses '|' as separator for multiple submissions.
    Also handles ';' and ',' separators found in older releases.
    """
    if not clnsig_raw or clnsig_raw == ".":
        return []
    # Primary separator in current ClinVar VCFs is '|'
    parts = clnsig_raw.replace(",", "|").split("|")
    return [_clean_value(p) for p in parts if p.strip()]


def classify_clnsig(clnsig_raw: str) -> Dict:
    """
    Classify a CLNSIG value string.

    Returns dict with:
      - category: 'Pathogenic' | 'Likely_pathogenic' | 'Pathogenic/Likely_pathogenic'
                  | 'Risk_allele' | 'Conflicting' | 'Benign' | 'Uncertain' | 'Other'
      - is_pathogenic: bool
      - is_pathogenic_or_likely: bool
      - raw_value: original string
      - parsed_values: list of individual classifications
    """
    values = _split_clnsig(clnsig_raw)
    if not values:
        return {
            "category": "Unknown",
            "is_pathogenic": False,
            "is_pathogenic_or_likely": False,
            "raw_value": clnsig_raw,
            "parsed_values": [],
        }

    pathogenic_hits = [v for v in values if v in PATHOGENIC_ALL]
    has_pathogenic = any(v in PATHOGENIC_DIRECT for v in values)
    has_likely_pathogenic = any("likely_pathogenic" in v for v in values)
    has_risk = any(v in PATHOGENIC_RISK for v in values)

    # Determine primary category
    if "pathogenic" in values and "likely_pathogenic" in values:
        category = "Pathogenic/Likely_pathogenic"
    elif "pathogenic" in values or any("pathogenic" in v for v in values if "low_penetrance" in v):
        category = "Pathogenic"
    elif "likely_pathogenic" in values or any("likely_pathogenic" in v for v in values):
        category = "Likely_pathogenic"
    elif has_risk:
        category = "Risk_allele"
    elif "conflicting_classifications_of_pathogenicity" in values:
        category = "Conflicting"
    elif any(v in EXCLUDED_CLNSIG for v in values):
        # Check which excluded category
        for v in values:
            if "benign" in v:
                category = "Benign"
                break
            elif "uncertain" in v:
                category = "Uncertain_significance"
                break
            elif v == "drug_response":
                category = "Drug_response"
                break
        else:
            category = "Other"
    else:
        category = "Other"

    return {
        "category": category,
        "is_pathogenic": len(pathogenic_hits) > 0,
        "is_pathogenic_or_likely": has_pathogenic or has_likely_pathogenic,
        "raw_value": clnsig_raw,
        "parsed_values": values,
    }


def is_pathogenic(clnsig_raw: str) -> bool:
    """Quick check: does this CLNSIG value indicate pathogenicity?"""
    return classify_clnsig(clnsig_raw)["is_pathogenic"]


# ── Confidence tier classification ────────────────────────────────────────────

# CLNREVSTAT → confidence tier mapping
CONFIDENCE_TIER_MAP = {
    # Tier 1 — Highest confidence
    "practice_guideline": "high",
    "reviewed_by_expert_panel": "high",
    # Tier 2 — Good confidence
    "criteria_provided,_multiple_submitters,_no_conflicts": "moderate",
    # Tier 3 — Lower confidence
    "criteria_provided,_single_submitter": "low",
    "criteria_provided,_conflicting_interpretations": "low",
    # Tier 4 — Very low / no evidence criteria stated
    "no_assertion_criteria_provided": "very_low",
    "no_assertion_provided": "very_low",
    "no_classification_for_the_single_variant": "very_low",
    "classified_by_submitter": "very_low",
}

CONFIDENCE_TIER_LABELS = {
    "high": {"en": "High Confidence", "es": "Confianza Alta"},
    "moderate": {"en": "Moderate Confidence", "es": "Confianza Moderada"},
    "low": {"en": "Lower Confidence", "es": "Confianza Baja"},
    "very_low": {"en": "Very Low Confidence", "es": "Confianza Muy Baja"},
}

CONFIDENCE_TIER_SYMBOLS = {
    "high": "🏅",
    "moderate": "✓",
    "low": "⚠️",
    "very_low": "❓",
}


def classify_confidence_tier(clnrevstat: str, clnsig: str = "") -> str:
    """
    Classify a ClinVar variant into a confidence tier based on review status
    and clinical significance.

    Tiers:
      - high: practice guideline or expert panel review
      - moderate: multiple submitters agree, criteria provided
      - low: single submitter with criteria
      - very_low: no assertion criteria provided (risk alleles with no evidence criteria)
    """
    if not clnrevstat:
        return "very_low"

    cleaned = _clean_value(clnrevstat)

    if cleaned in CONFIDENCE_TIER_MAP:
        return CONFIDENCE_TIER_MAP[cleaned]

    # Heuristic: if it contains "practice_guideline" or "expert_panel"
    if "practice_guideline" in cleaned or "expert_panel" in cleaned:
        return "high"
    if "multiple_submitters" in cleaned and "no_conflicts" in cleaned:
        return "moderate"
    if "single_submitter" in cleaned:
        return "low"
    if "no_assertion" in cleaned or "no_classification" in cleaned:
        return "very_low"

    return "very_low"


# ── VCF parsing utilities ──────────────────────────────────────────────────────


def parse_clinvar_info(info_str: str) -> Dict[str, str]:
    """
    Parse ClinVar VCF INFO column into a dict.
    Handles quoted values and multi-value fields.
    """
    result = {}
    for field in info_str.split(";"):
        if "=" in field:
            key, val = field.split("=", 1)
            result[key] = val
        else:
            result[field] = "1"  # Flag fields
    return result


def extract_user_vcf_regions(user_vcf: str, output_bed: str) -> int:
    """
    Extract CHROM\\tPOS\\tPOS\\tREF\\tALT from user VCF as BED-like file.
    Uses bcftools query for speed on large VCFs.
    Returns number of variants extracted.
    """
    logger.info("  Extracting variant positions from user VCF ...")

    # Check VCF exists and has index
    vcf_path = Path(user_vcf)
    if not vcf_path.exists():
        raise FileNotFoundError(f"User VCF not found: {user_vcf}")

    tbi_path = Path(str(user_vcf) + ".tbi")
    if not tbi_path.exists() and not Path(str(user_vcf).replace(".gz", "") + ".tbi").exists():
        logger.info("  Creating tabix index for user VCF ...")
        subprocess.run(
            ["bcftools", "index", "-t", user_vcf],
            check=True, capture_output=True, timeout=600,
        )

    # Extract: CHROM POS POS REF ALT (tab-separated, no header)
    # Normalize chromosome names: strip "chr" prefix for ClinVar compatibility
    cmd = [
        "bash", "-c",
        f"bcftools query -f '%CHROM\\t%POS\\t%POS\\t%REF\\t%ALT\\n' {user_vcf} | sed 's/^chr//i'"
    ]

    try:
        result = subprocess.run(
            cmd, capture_output=True, text=True, timeout=600,
        )
        if result.returncode != 0:
            raise RuntimeError(f"bcftools query failed: {result.stderr}")

        lines = result.stdout.strip().split("\n")
        with open(output_bed, "w") as fh:
            fh.write(result.stdout)

        n_variants = len([l for l in lines if l.strip()])
        logger.info(f"  ✓ {n_variants:,} variant positions extracted")
        return n_variants

    except subprocess.TimeoutExpired:
        raise RuntimeError("bcftools query timed out on user VCF")


def query_clinvar_by_regions(clinvar_vcf: str, regions_file: str, output_vcf: str) -> int:
    """
    Query ClinVar VCF using tabix -R for positional overlap.
    Returns number of overlapping ClinVar records.
    """
    logger.info("  Querying ClinVar for positional overlaps via tabix (~30-60s)...")

    # Ensure ClinVar has index
    clinvar_path = Path(clinvar_vcf)
    if not clinvar_path.exists():
        raise FileNotFoundError(
            f"ClinVar VCF not found: {clinvar_vcf}\n"
            f"Run: python scripts/setup/download_clinvar.py"
        )

    tbi_path = Path(str(clinvar_vcf) + ".tbi")
    if not tbi_path.exists():
        raise FileNotFoundError(
            f"ClinVar index not found: {tbi_path}\n"
            f"Run: python scripts/setup/download_clinvar.py"
        )

    cmd = ["tabix", "-R", regions_file, clinvar_vcf]

    try:
        result = subprocess.run(cmd, capture_output=True, timeout=600)
        if result.returncode != 0:
            stderr = result.stderr.decode("utf-8", errors="replace")
            raise RuntimeError(
                f"tabix query failed: {stderr}\n"
                f"Check that both VCFs use the same chromosome naming (chr1 vs 1)."
            )

        output = result.stdout.decode("utf-8", errors="replace")
        with open(output_vcf, "w") as fh:
            fh.write(output)

        n_overlap = len([l for l in output.split("\n") if l.strip() and not l.startswith("#")])
        logger.info(f"  ✓ {n_overlap:,} positional overlaps with ClinVar")
        return n_overlap

    except subprocess.TimeoutExpired:
        raise RuntimeError(
            "tabix query timed out (>10 min). The ClinVar VCF may be corrupted.\n"
            "Run: python scripts/setup/download_clinvar.py --force to re-download."
        )


def load_user_vcf_variants(user_vcf: str) -> Dict[Tuple, Dict]:
    """
    Load user VCF variants into a dict keyed by (chrom, pos, ref, alt).
    Only loads CHROM, POS, REF, ALT — no genotypes needed.
    Supports multi-allelic sites (splits ALT by comma).
    """
    logger.info("  Loading user VCF variant dictionary ...")
    variants = {}

    cmd = [
        "bcftools", "query",
        "-f", "%CHROM\t%POS\t%REF\t%ALT\n",
        user_vcf,
    ]

    try:
        result = subprocess.run(
            cmd, capture_output=True, text=True, timeout=600,
        )
        if result.returncode != 0:
            raise RuntimeError(f"bcftools query failed: {result.stderr}")

        for line in result.stdout.strip().split("\n"):
            parts = line.strip().split("\t")
            if len(parts) < 4:
                continue
            chrom, pos, ref, alt_str = parts[0], parts[1], parts[2], parts[3]
            chrom = _strip_chr(chrom)  # Normalize: "chr1" → "1"

            # Handle multi-allelic: ALT can be "A,G,T"
            alts = alt_str.split(",")
            for alt in alts:
                if alt and alt != "." and alt != "<*>":
                    key = (chrom, pos, ref.upper(), alt.upper())
                    variants[key] = {
                        "chrom": chrom,
                        "pos": int(pos),
                        "ref": ref,
                        "alt": alt,
                    }

        logger.info(f"  ✓ {len(variants):,} unique variant alleles loaded")
        return variants

    except subprocess.TimeoutExpired:
        raise RuntimeError("bcftools query timed out")


def annotate_clinvar_overlap(
    clinvar_overlap_vcf: str,
    user_variants: Dict[Tuple, Dict],
) -> Dict:
    """
    Parse ClinVar overlap VCF, match on exact REF/ALT against user variants,
    and classify CLNSIG. Returns structured annotation results.
    """
    logger.info("  Matching ClinVar records to user variants ...")

    pathogenic_variants = []
    all_matches = []
    stats = {
        "clinvar_records_processed": 0,
        "positional_overlaps": 0,
        "exact_matches": 0,
    }
    # Separate counters: all matches by category vs pathogenic subset
    all_match_categories = Counter()
    pathogenic_counts = Counter()

    review_status_counts = Counter()
    confidence_tier_counts = Counter()
    gene_counts = Counter()
    chrom_counts = Counter()

    with open(clinvar_overlap_vcf, "r") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            stats["positional_overlaps"] += 1
            parts = line.split("\t")
            if len(parts) < 8:
                continue

            clinvar_chrom = parts[0]
            clinvar_pos = int(parts[1])
            clinvar_ref = parts[3].upper()
            clinvar_alt = parts[4].upper()

            # Try matching each ALT allele
            matched = False
            for alt in clinvar_alt.split(","):
                key = (clinvar_chrom, str(clinvar_pos), clinvar_ref, alt)
                if key in user_variants:
                    matched = True
                    stats["exact_matches"] += 1

                    # Parse INFO fields
                    info = parse_clinvar_info(parts[7])

                    clnsig_raw = info.get("CLNSIG", "")
                    clndn = info.get("CLNDN", "")
                    clndisdb = info.get("CLNDISDB", "")
                    geneinfo = info.get("GENEINFO", "")
                    clnrevstat = info.get("CLNREVSTAT", "")
                    rsid = info.get("RS", "")
                    clnvc = info.get("CLNVC", "")
                    allele_id = info.get("ALLELEID", "")
                    af_exac = info.get("AF_EXAC", "")
                    af_tgp = info.get("AF_TGP", "")
                    origin = info.get("ORIGIN", "")

                    # Classify
                    classification = classify_clnsig(clnsig_raw)
                    all_match_categories[classification["category"]] += 1

                    # Parse gene info
                    genes = []
                    if geneinfo:
                        for g in geneinfo.split("|"):
                            gene_name = g.split(":")[0] if ":" in g else g
                            if gene_name:
                                genes.append(gene_name)

                    confidence_tier = classify_confidence_tier(clnrevstat, clnsig_raw)

                    variant_entry = {
                        "chrom": clinvar_chrom,
                        "pos": clinvar_pos,
                        "ref": clinvar_ref,
                        "alt": alt,
                        "rsid": rsid,
                        "clinical_significance": classification["category"],
                        "clnsig_raw": clnsig_raw,
                        "disease_name": clndn,
                        "disease_dbs": clndisdb,
                        "gene_info": geneinfo,
                        "genes": genes,
                        "review_status": clnrevstat,
                        "confidence_tier": confidence_tier,
                        "variant_type": clnvc,
                        "allele_id": allele_id,
                        "af_exac": af_exac,
                        "af_1000g": af_tgp,
                        "origin": origin,
                    }

                    all_matches.append(variant_entry)

                    if classification["is_pathogenic"]:
                        pathogenic_variants.append(variant_entry)
                        pathogenic_counts[classification["category"]] += 1
                        confidence_tier_counts[confidence_tier] += 1
                        for g in genes:
                            gene_counts[g] += 1
                        chrom_counts[clinvar_chrom] += 1
                        if clnrevstat:
                            review_status_counts[_clean_value(clnrevstat)] += 1

            stats["clinvar_records_processed"] += 1

    # Sort pathogenic variants: Pathogenic first, then Likely_pathogenic, then by chromosome
    severity_order = {
        "Pathogenic": 0,
        "Pathogenic/Likely_pathogenic": 1,
        "Likely_pathogenic": 2,
        "Risk_allele": 3,
    }
    pathogenic_variants.sort(
        key=lambda v: (
            severity_order.get(v["clinical_significance"], 9),
            v["chrom"],
            v["pos"],
        )
    )

    logger.info(f"  ✓ {stats['exact_matches']:,} exact REF/ALT matches with ClinVar")
    logger.info(f"  ✓ {len(pathogenic_variants)} pathogenic/likely pathogenic variants found")

    return {
        "pathogenic_variants": pathogenic_variants,
        "all_matches": all_matches,
        "stats": stats,
        "all_match_categories": dict(all_match_categories),
        "pathogenic_counts": dict(pathogenic_counts),
        "confidence_tier_counts": dict(confidence_tier_counts),
        "review_status_counts": dict(review_status_counts),
        "gene_counts": dict(gene_counts.most_common(50)),
        "chrom_counts": dict(sorted(chrom_counts.items())),
    }


def build_output_json(
    user_vcf: str,
    clinvar_vcf: str,
    user_variant_count: int,
    clinvar_variant_count: int,
    annotation: Dict,
    output_dir: str,
) -> Dict:
    """Build the final output JSON structure."""

    stats = annotation["stats"]
    pathogenic_counts = annotation.get("pathogenic_counts", {})
    all_match_categories = annotation.get("all_match_categories", {})
    confidence_tier_counts = annotation.get("confidence_tier_counts", {})

    # Get ClinVar release date from manifest
    clinvar_release_date = ""
    clinvar_manifest_path = Path(clinvar_vcf).parent / "manifest.json"
    if clinvar_manifest_path.exists():
        try:
            manifest = json.loads(clinvar_manifest_path.read_text())
            clinvar_release_date = manifest.get("download_date", "")
        except Exception:
            pass

    metadata = {
        "analysis_date": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "pipeline_version": PIPELINE_VERSION,
        "user_vcf": user_vcf,
        "clinvar_vcf": clinvar_vcf,
        "user_vcf_total_variants": user_variant_count,
        "clinvar_total_variants": clinvar_variant_count,
        "positional_overlaps": stats["positional_overlaps"],
        "exact_matches": stats["exact_matches"],
        "match_rate": round(stats["exact_matches"] / max(stats["positional_overlaps"], 1) * 100, 1),
        "clinvar_release_date": clinvar_release_date,
    }

    # Get pathogenic variants list (needed for cross-reference and total count)
    pathogenic_variants = annotation["pathogenic_variants"]

    pathogenic_variant_summary = {
        "total_pathogenic": pathogenic_counts.get("Pathogenic", 0),
        "total_likely_pathogenic": pathogenic_counts.get("Likely_pathogenic", 0),
        "total_pathogenic_or_likely_pathogenic": pathogenic_counts.get("Pathogenic/Likely_pathogenic", 0),
        "total_risk_alleles": pathogenic_counts.get("Risk_allele", 0),
        "total_combined": len(pathogenic_variants),
        "total_benign": all_match_categories.get("Benign", 0),
        "total_uncertain": all_match_categories.get("Uncertain_significance", 0),
        "total_conflicting": all_match_categories.get("Conflicting", 0),
        "by_review_status": annotation["review_status_counts"],
        "by_confidence_tier": confidence_tier_counts,
        "high_confidence_count": confidence_tier_counts.get("high", 0) + confidence_tier_counts.get("moderate", 0),
        "by_gene": annotation["gene_counts"],
        "by_chromosome": annotation["chrom_counts"],
    }

    # Cross-reference with curated SNP database
    snp_db_path = Path(output_dir).parent / "data" / "snp_database_annotated.csv"
    curated_rsids = set()
    if snp_db_path.exists():
        import csv
        try:
            with open(snp_db_path) as fh:
                reader = csv.DictReader(fh)
                for row in reader:
                    rsid = row.get("rsid", "").strip()
                    if rsid:
                        curated_rsids.add(rsid)
        except Exception:
            pass

    for v in pathogenic_variants:
        v["in_snp_database"] = v.get("rsid", "") in curated_rsids if v.get("rsid") else False

    if curated_rsids:
        pathogenic_variant_summary["in_snp_database"] = sum(
            1 for v in pathogenic_variants if v.get("in_snp_database")
        )

    return {
        "metadata": metadata,
        "pathogenic_variants": pathogenic_variants,
        "pathogenic_variant_summary": pathogenic_variant_summary,
    }


# ── Main ───────────────────────────────────────────────────────────────────────


def main():
    import argparse

    parser = argparse.ArgumentParser(
        description="ClinVar pathogenic variant annotation — annotates user VCF against ClinVar"
    )
    parser.add_argument("--vcf", required=True,
                        help="User VCF path (.vcf.gz)")
    parser.add_argument("--clinvar-vcf", default="reference/clinvar/clinvar.vcf.gz",
                        help="ClinVar VCF path (default: reference/clinvar/clinvar.vcf.gz)")
    parser.add_argument("--output-dir", "-o", default="clinvar",
                        help="Output directory (default: clinvar/)")
    parser.add_argument("--temp-dir", default=None,
                        help="Temporary directory for intermediate files")
    parser.add_argument("--verbose", "-v", action="store_true",
                        help="Verbose output")
    parser.add_argument("--keep-temp", action="store_true",
                        help="Keep intermediate files for debugging")
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s" if args.verbose else "%(message)s",
    )

    logger.info("═══ ClinVar Pathogenic Variant Annotation ═══")

    # Resolve paths relative to pipeline root
    pipeline_root = Path(__file__).resolve().parent.parent.parent
    output_dir = pipeline_root / args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    clinvar_vcf = args.clinvar_vcf
    if not Path(clinvar_vcf).is_absolute():
        clinvar_vcf = str(pipeline_root / clinvar_vcf)
    user_vcf = args.vcf if Path(args.vcf).is_absolute() else str(Path.cwd() / args.vcf)

    if not Path(clinvar_vcf).exists():
        logger.error(f"ClinVar VCF not found: {clinvar_vcf}")
        logger.error("Run: python scripts/setup/download_clinvar.py")
        return 1

    if not Path(user_vcf).exists():
        logger.error(f"User VCF not found: {user_vcf}")
        return 1

    # Setup temp directory
    if args.temp_dir:
        temp_dir = Path(args.temp_dir)
        temp_dir.mkdir(parents=True, exist_ok=True)
        _cleanup = False
    else:
        temp_dir = Path(tempfile.mkdtemp(prefix="clinvar_work_"))
        _cleanup = True

    try:
        # Step 1: Extract user VCF positions
        regions_file = str(temp_dir / "user_regions.tsv")
        n_user = extract_user_vcf_regions(user_vcf, regions_file)

        if n_user == 0:
            logger.warning("  No variants found in user VCF")
            # Write empty output
            _write_empty_output(output_dir, user_vcf, clinvar_vcf)
            return 0

        # Step 2: Tabix query ClinVar
        overlap_vcf = str(temp_dir / "clinvar_overlap.vcf")
        n_overlap = query_clinvar_by_regions(clinvar_vcf, regions_file, overlap_vcf)

        if n_overlap == 0:
            logger.warning("  No positional overlaps with ClinVar — output will be empty")
            _write_empty_output(output_dir, user_vcf, clinvar_vcf, n_user_variants=n_user)
            return 0

        # Step 3: Load user variants and match
        user_variants = load_user_vcf_variants(user_vcf)

        # Step 4: Annotate and classify
        annotation = annotate_clinvar_overlap(overlap_vcf, user_variants)

        # Count ClinVar variants (from manifest if available, else estimate)
        clinvar_manifest = Path(clinvar_vcf).parent / "manifest.json"
        n_clinvar = 4_400_000  # default estimate
        if clinvar_manifest.exists():
            try:
                n_clinvar = json.loads(clinvar_manifest.read_text()).get("n_variants", n_clinvar)
            except Exception:
                pass

        # Step 5: Build and write output
        output = build_output_json(
            user_vcf, clinvar_vcf, n_user, n_clinvar, annotation, str(output_dir),
        )

        output_path = output_dir / "clinvar_pathogenic_variants.json"
        output_path.write_text(json.dumps(output, indent=2, ensure_ascii=False))
        logger.info(f"  ✅ {output_path}")

        # Summary
        summary = output["pathogenic_variant_summary"]
        total = summary["total_combined"]
        logger.info(f"\n═══ Summary ═══")
        logger.info(f"  Pathogenic:              {summary['total_pathogenic']}")
        logger.info(f"  Likely_pathogenic:       {summary['total_likely_pathogenic']}")
        logger.info(f"  Pathogenic/Likely:       {summary['total_pathogenic_or_likely_pathogenic']}")
        logger.info(f"  Risk alleles:            {summary['total_risk_alleles']}")
        logger.info(f"  ─────────────────────────")
        logger.info(f"  TOTAL pathogenic/likely: {total}")
        logger.info(f"  All matches — Benign:    {summary.get('total_benign', 0):,}")
        logger.info(f"  All matches — Uncertain: {summary.get('total_uncertain', 0):,}")
        logger.info(f"  All matches — Conflicting:{summary.get('total_conflicting', 0):,}")

        if summary.get("by_gene"):
            logger.info(f"\n  Top genes:")
            for gene, count in list(summary["by_gene"].items())[:10]:
                logger.info(f"    {gene}: {count}")

        logger.info(f"\n  Results in: {output_dir}/")

        return 0

    except Exception as e:
        logger.error(f"  ✗ Annotation failed: {e}")
        import traceback
        logger.debug(traceback.format_exc())
        return 1

    finally:
        # Cleanup temp files
        if _cleanup and not args.keep_temp and temp_dir.exists():
            import shutil
            shutil.rmtree(temp_dir, ignore_errors=True)


def _write_empty_output(
    output_dir: Path,
    user_vcf: str,
    clinvar_vcf: str,
    n_user_variants: int = 0,
):
    """Write an empty but valid ClinVar annotation JSON."""
    output = {
        "metadata": {
            "analysis_date": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
            "pipeline_version": PIPELINE_VERSION,
            "user_vcf": user_vcf,
            "clinvar_vcf": clinvar_vcf,
            "user_vcf_total_variants": n_user_variants,
            "clinvar_total_variants": 0,
            "positional_overlaps": 0,
            "exact_matches": 0,
            "match_rate": 0,
            "clinvar_release_date": "",
        },
        "pathogenic_variants": [],
        "pathogenic_variant_summary": {
            "total_pathogenic": 0,
            "total_likely_pathogenic": 0,
            "total_pathogenic_or_likely_pathogenic": 0,
            "total_risk_alleles": 0,
            "total_combined": 0,
            "by_review_status": {},
            "by_confidence_tier": {},
            "high_confidence_count": 0,
            "by_gene": {},
            "by_chromosome": {},
        },
    }
    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / "clinvar_pathogenic_variants.json"
    output_path.write_text(json.dumps(output, indent=2, ensure_ascii=False))
    logger.info(f"  Empty output written to {output_path}")


if __name__ == "__main__":
    sys.exit(main())
