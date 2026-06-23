#!/usr/bin/env python3
"""
PRS Pipeline Test Suite — validates all outputs, catches numpy.bool_ regressions,
verifies JSON integrity, and checks data consistency.

Usage:
  python3 prs.py test              # Run all tests
  python3 prs.py test --quick      # JSON integrity only (fast)
  python3 prs.py test --verbose    # Show all checks
"""

import sys
import os
import json
import logging
from pathlib import Path
from typing import Any, Tuple, List
from datetime import datetime
from collections import defaultdict

logger = logging.getLogger(__name__)

PLATFORM_DIR = Path(__file__).parent.parent.parent  # utils/ → scripts/ → prs_research_pipeline/

# ── Test results tracking ───────────────────────────────────────────────────

PASS, FAIL, WARN = 0, 0, 0
RESULTS = []

def ok(msg: str): global PASS; PASS += 1; RESULTS.append(("✅", msg)); print(f"  ✅ {msg}")
def bad(msg: str): global FAIL; FAIL += 1; RESULTS.append(("❌", msg)); print(f"  ❌ {msg}")
def warn(msg: str): global WARN; WARN += 1; RESULTS.append(("⚠️ ", msg)); print(f"  ⚠️  {msg}")


# ═══════════════════════════════════════════════════════════════════════════════
# TEST 1: JSON integrity — no numpy types, valid JSON, required keys
# ═══════════════════════════════════════════════════════════════════════════════

EXPECTED_JSON_FILES = [
    # (path, required_keys, description)
    ("prs/PRS_RESULT.json", ["sample_id", "prs_entries", "metadata"], "Unified PRS output"),
    ("science/global_validation_report.json", ["overall_score", "checks"], "Global validation"),
    ("science/CONSOLIDATION_MANIFEST.json", ["ssst_sources", "summary"], "SSST manifest"),
    ("science/ANCESTRY_MODEL.json", ["assigned_population", "method"], "Ancestry model"),
    ("science/adversarial_validation_report.json", ["results", "overall_robustness_score"], "Adversarial validation"),
    ("science/failure_mode_map.json", ["failures", "n_failures"], "Failure mode map"),
    ("science/leakage_audit.json", ["checks", "pipeline_safe"], "Leakage audit"),
    ("science/pipeline_gate_check.json", ["checks", "pipeline_can_proceed"], "Pipeline gates"),
    ("science/snp_universe.json", ["unified_snp_count"], "SNP universe"),
    ("FINAL_SCIENTIFIC_SCORE.json", ["scientific_integrity_score", "components"], "Final score"),
    ("PUBLICATION_LOCK.md", None, "Publication lock"),
    ("benchmark/VALIDATION_REPORT.json", ["entries", "validation_summary"], "Benchmark validation"),
    ("benchmark/gwas_consortium_validation.json", ["validations", "consortia"], "GWAS consortium"),
    ("benchmark/portability_report.json", ["populations", "global_bias_index"], "Portability"),
    ("benchmark/quality_delta.json", ["components", "mean_delta"], "Quality delta"),
    ("prs/population_calibration_report.json", ["traits_analyzed", "methodology"], "Calibration report"),
    ("prs/uncertainty_report.json", ["results", "global_uncertainty_score"], "Uncertainty report"),
    ("prs/consistency_check_report.json", ["passed", "gwas_ancestry_match"], "Consistency check"),
    ("reproducibility/run_fingerprint.json", ["run_id", "environment"], "Run fingerprint"),
    ("reports/comprehensive_report_en.html", None, "EN HTML report"),
    ("reports/comprehensive_report_es.html", None, "ES HTML report"),
    ("clinvar/clinvar_pathogenic_variants.json",
     ["metadata", "pathogenic_variants", "pathogenic_variant_summary"],
     "ClinVar pathogenic annotation"),
    ("pharmgkb/pharmgkb_drug_report.json",
     ["metadata", "pharmacogenomic_findings", "summary"],
     "Pharmacogenomic drug report"),
]


def validate_json_type(value: Any, path: str = "$") -> List[str]:
    """Recursively check for numpy types and non-serializable objects."""
    issues = []

    # Check for numpy scalars masquerading as Python types
    type_name = type(value).__module__
    if type_name == "numpy":
        if hasattr(value, "item"):
            issues.append(f"{path}: numpy.{type(value).__name__} — should be Python {type(value.item()).__name__}")
        else:
            issues.append(f"{path}: numpy.{type(value).__name__}")

    # Specific checks for numpy bool (our #1 regression)
    if hasattr(value, "dtype"):
        issues.append(f"{path}: has .dtype={value.dtype} — likely numpy scalar")

    if isinstance(value, dict):
        for k, v in value.items():
            issues.extend(validate_json_type(v, f"{path}.{k}"))
    elif isinstance(value, list):
        for i, v in enumerate(value):
            issues.extend(validate_json_type(v, f"{path}[{i}]"))

    return issues


def test_json_integrity():
    """Validate all JSON outputs are valid and free of numpy types."""
    print("\n── JSON Integrity ──")

    for path, required_keys, desc in EXPECTED_JSON_FILES:
        full_path = PLATFORM_DIR / path

        if not full_path.exists():
            warn(f"{desc}: file missing ({path})")
            continue

        # Skip non-JSON files (HTML reports, etc.)
        if not path.endswith(".json"):
            if full_path.exists():
                ok(f"{desc}: file exists ({full_path.stat().st_size / 1024:.0f} KB)")
            else:
                warn(f"{desc}: file missing ({path})")
            continue

        # Validate valid JSON
        try:
            with open(full_path) as fh:
                data = json.load(fh)
        except json.JSONDecodeError as e:
            bad(f"{desc}: invalid JSON — {e}")
            continue

        # Check for numpy types recursively
        issues = validate_json_type(data, "$")
        if issues:
            bad(f"{desc}: {len(issues)} numpy type(s) found")
            for issue in issues[:5]:
                logger.warning(f"    {issue}")
        else:
            ok(f"{desc}: valid JSON, no numpy types")

        # Check required keys
        if required_keys:
            missing = [k for k in required_keys if k not in data]
            if missing:
                warn(f"{desc}: missing keys: {missing}")

        # Check file is not empty
        if isinstance(data, dict) and not data:
            warn(f"{desc}: empty object")
        if isinstance(data, list) and not data:
            warn(f"{desc}: empty list")


# ═══════════════════════════════════════════════════════════════════════════════
# TEST 2: Data consistency checks
# ═══════════════════════════════════════════════════════════════════════════════

def test_data_consistency():
    """Verify data consistency across pipeline outputs."""
    print("\n── Data Consistency ──")

    # PRS_RESULT consistency
    prs_path = PLATFORM_DIR / "prs/PRS_RESULT.json"
    if prs_path.exists():
        with open(prs_path) as fh:
            prs = json.load(fh)
        entries = prs.get("prs_entries", [])

        if len(entries) == 10:
            ok(f"PRS_RESULT: 10 traits")
        elif len(entries) > 0:
            warn(f"PRS_RESULT: {len(entries)} traits (expected 10)")
        else:
            bad(f"PRS_RESULT: 0 traits — pipeline not run?")

        # Check for all-zero z-scores (synthetic calibration)
        all_zero_z = all(float(e.get("population_zscore", 0)) == 0 for e in entries)
        if all_zero_z:
            warn(f"PRS_RESULT: all z-scores = 0 — synthetic calibration? Run build_reference_distributions")
        else:
            ok(f"PRS_RESULT: real z-scores detected")

        # Check reasonable value ranges
        for e in entries:
            z = float(e.get("population_zscore", 0))
            pctl = float(e.get("population_percentile", 50))
            raw = float(e.get("raw_score", 0))
            if pctl < 0 or pctl > 100:
                bad(f"PRS_RESULT: {e['trait']} percentile out of range: {pctl}")
            if abs(z) > 10:
                warn(f"PRS_RESULT: {e['trait']} z-score extreme: {z}")

    # Validation score consistency
    val_path = PLATFORM_DIR / "science/global_validation_report.json"
    fin_path = PLATFORM_DIR / "FINAL_SCIENTIFIC_SCORE.json"
    if val_path.exists():
        with open(val_path) as fh:
            val = json.load(fh)
        n_checks = val.get("total_checks", 0)
        score = val.get("overall_score", 0)
        if n_checks >= 8:
            ok(f"Validation: {n_checks} dimensions checked")
        else:
            warn(f"Validation: only {n_checks} dimensions")
        if score > 0:
            ok(f"Validation: score {score:.0f}/100")
        else:
            warn(f"Validation: score is 0")

    # Integrity score component check
    if fin_path.exists():
        with open(fin_path) as fh:
            fin = json.load(fh)
        components = fin.get("components", [])
        if len(components) >= 6:
            ok(f"Integrity: {len(components)} components")
        else:
            warn(f"Integrity: only {len(components)} components")
        total = fin.get("scientific_integrity_score", 0)
        if total > 0:
            ok(f"Integrity: score {total:.0f}/100")
        else:
            bad(f"Integrity: score is 0")

    # Failure map
    fm_path = PLATFORM_DIR / "science/failure_mode_map.json"
    if fm_path.exists():
        with open(fm_path) as fh:
            fm = json.load(fh)
        n = fm.get("n_failures", 0)
        if n >= 17:
            ok(f"Failure map: {n} modes covered")
        else:
            warn(f"Failure map: only {n} modes")

    # Leakage audit
    lk_path = PLATFORM_DIR / "science/leakage_audit.json"
    if lk_path.exists():
        with open(lk_path) as fh:
            lk = json.load(fh)
        if lk.get("pipeline_safe"):
            ok(f"Leakage: pipeline SAFE")
        if lk.get("passed", 0) == lk.get("total_checks", 0):
            ok(f"Leakage: all {lk['total_checks']} checks passed")


# ═══════════════════════════════════════════════════════════════════════════════
# TEST 3: Core pipeline artifacts
# ═══════════════════════════════════════════════════════════════════════════════

def test_pipeline_artifacts():
    """Verify critical pipeline output files exist."""
    print("\n── Pipeline Artifacts ──")

    # Genotype processing outputs
    checks = [
        ("plink/cohort.bed", "Stage A: PLINK binary"),
        ("plink/cohort.bim", "Stage A: variant map"),
        ("plink/cohort.fam", "Stage A: sample info"),
        ("qc/qc_filtered.bed", "Stage B: QC'd genotypes"),
        ("plink/ld_pruned_dataset.bed", "Stage C: LD-pruned dataset"),
        ("pca/pca_results.eigenvec", "Stage D: PCA results"),
        ("pca/target_pcs.eigenvec", "Stage D: Target sample PCs"),
        ("prs/PRS_RESULT.json", "Stage F-H: PRS output"),
        ("prs/PRS_RESULT.csv", "Stage F-H: PRS CSV"),
        ("prs/population_calibration_report.json", "Stage H: Calibration"),
        ("prs/uncertainty_report.json", "Stage I: Uncertainty"),
        ("interpretations/interpretations_en.json", "Stage K: EN interpretation"),
        ("interpretations/interpretations_es.json", "Stage K: ES interpretation"),
        ("reproducibility/run_fingerprint.json", "Phase 7: Reproducibility"),
        ("science/assumptions.lock.json", "Phase 7: Scientific lock"),
    ]

    for path, desc in checks:
        if (PLATFORM_DIR / path).exists():
            ok(f"{desc}")
        else:
            warn(f"{desc}: missing ({path})")

    # Reference data
    ref_full = PLATFORM_DIR / "reference/1000G_full/1000G_full.bed"
    ref_chr22 = PLATFORM_DIR / "reference/1000G/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
    if ref_full.exists():
        ok("1000G reference: genome-wide (full)")
    elif ref_chr22.exists():
        warn("1000G reference: chr22 only — run scripts/setup/download_1000G_full.py")
    else:
        warn("1000G reference: MISSING — run scripts/setup/download_1000G_full.py")


# ═══════════════════════════════════════════════════════════════════════════════
# TEST 4: Variant and sample counts
# ═══════════════════════════════════════════════════════════════════════════════

def test_variant_counts():
    """Verify variant and sample counts across stages."""
    print("\n── Variant & Sample Counts ──")

    # Target sample variants
    bim_initial = PLATFORM_DIR / "plink/cohort.bim"
    bim_qc = PLATFORM_DIR / "qc/qc_filtered.bim"
    bim_pruned = PLATFORM_DIR / "plink/ld_pruned_dataset.bim"
    fam = PLATFORM_DIR / "plink/cohort.fam"

    if fam.exists():
        n_samples = sum(1 for _ in open(fam))
        if n_samples == 1:
            ok(f"Target sample: 1 sample detected")
        elif n_samples >= 2:
            ok(f"Target sample: {n_samples} samples (multi-sample mode)")
        else:
            warn(f"Target sample: 0 samples — no data")

    if bim_initial.exists():
        n = sum(1 for _ in open(bim_initial))
        ok(f"Stage A: {n:,} variants imported from VCF")

    if bim_qc.exists():
        n = sum(1 for _ in open(bim_qc))
        ok(f"Stage B: {n:,} variants after QC")

    if bim_pruned.exists():
        n = sum(1 for _ in open(bim_pruned))
        if n > 20000:
            ok(f"Stage C: {n:,} LD-pruned SNPs (genome-wide)")
        elif n > 300:
            ok(f"Stage C: {n:,} LD-pruned SNPs (chr22)")
        else:
            warn(f"Stage C: {n:,} SNPs — very low, check reference")

    # Reference variants
    ref_bim = PLATFORM_DIR / "reference/1000G_full/1000G_full.bim"
    if ref_bim.exists():
        n = sum(1 for _ in open(ref_bim))
        if n > 80_000_000:
            ok(f"Reference: {n:,} variants (genome-wide)")
        else:
            warn(f"Reference: {n:,} variants — expected 84M")


# ═══════════════════════════════════════════════════════════════════════════════
# TEST 5: Report content
# ═══════════════════════════════════════════════════════════════════════════════

def test_reports():
    """Check HTML reports have expected structure."""
    print("\n── Reports ──")

    for lang, path in [("EN", "reports/comprehensive_report_en.html"),
                       ("ES", "reports/comprehensive_report_es.html")]:
        full = PLATFORM_DIR / path
        if not full.exists():
            warn(f"[{lang}] HTML report: missing")
            continue

        content = full.read_text()
        size_kb = len(content) / 1024
        sections = content.count("section-header")

        if size_kb > 50:
            ok(f"[{lang}] HTML report: {size_kb:.0f} KB, {sections} sections")
        else:
            warn(f"[{lang}] HTML report: only {size_kb:.0f} KB — may be incomplete")


# ═══════════════════════════════════════════════════════════════════════════════
# TEST 6: Code integrity — all key scripts compile without syntax errors
# ═══════════════════════════════════════════════════════════════════════════════

def test_code_integrity():
    """Verify all key Python scripts compile and are importable."""
    print("\n── Code Integrity ──")

    KEY_SCRIPTS = [
        ("prs.py", Path(__file__).parent.parent.parent.parent / "prs.py"),  # root
        ("06_prs_compute.py", PLATFORM_DIR / "scripts" / "prs" / "06_prs_compute.py"),
        ("39_ancestry_model_unified.py", PLATFORM_DIR / "scripts" / "sss" / "39_ancestry_model_unified.py"),
        ("32_global_scientific_validator.py", PLATFORM_DIR / "scripts" / "validation" / "32_global_scientific_validator.py"),
        ("36_prs_core.py", PLATFORM_DIR / "scripts" / "sss" / "36_prs_core.py"),
        ("pca_ancestry_classifier.py", PLATFORM_DIR / "scripts" / "stages" / "pca_ancestry_classifier.py"),
        ("comprehensive_report.py", PLATFORM_DIR / "scripts" / "publication" / "comprehensive_report.py"),
        ("build_reference_distributions.py", PLATFORM_DIR / "scripts" / "utils" / "build_reference_distributions.py"),
        # download_gwas_full.py was removed — URLs expire in 2h, see GWAS_DATASETS.md
        ("download_dbsnp.py", PLATFORM_DIR / "scripts" / "setup" / "download_dbsnp.py"),
        ("27_real_world_calibration.py", PLATFORM_DIR / "scripts" / "benchmarking" / "27_real_world_calibration.py"),
    ]

    for name, path in KEY_SCRIPTS:
        if not path.exists():
            warn(f"Code: {name} — file not found at {path}")
            continue
        try:
            with open(path) as fh:
                source = fh.read()
            compile(source, name, "exec")
            ok(f"Code: {name} compiles")
        except SyntaxError as e:
            bad(f"Code: {name} — {e}")


# ═══════════════════════════════════════════════════════════════════════════════
# TEST 7: PRS output edge cases
# ═══════════════════════════════════════════════════════════════════════════════

def test_prs_edge_cases():
    """Validate PRS output structure and edge cases."""
    print("\n── PRS Edge Cases ──")

    # PRS raw CSV
    prs_csv = PLATFORM_DIR / "prs" / "prs_raw.csv"
    if prs_csv.exists():
        import csv
        try:
            with open(prs_csv) as fh:
                reader = csv.DictReader(fh)
                rows = list(reader)
            n_rows = len(rows)
            cols = set(rows[0].keys()) if rows else set()
            expected_cols = {"individual_id", "trait", "prs_raw", "n_snps", "n_snps_used"}
            missing_cols = expected_cols - cols

            if n_rows == 0:
                warn("PRS CSV: 0 rows — re-run pipeline")
            elif not missing_cols:
                ok(f"PRS CSV: {n_rows} rows, correct columns")
            else:
                bad(f"PRS CSV: missing columns {missing_cols}")

            # Multi-sample detection
            if rows:
                sample_ids = set(r["individual_id"] for r in rows)
                n_samples = len(sample_ids)
                n_traits = len(set(r["trait"] for r in rows))
                if n_samples >= 3:
                    ok(f"PRS CSV: {n_samples} samples × {n_traits} traits (multi-sample)")
                elif n_samples > 1:
                    ok(f"PRS CSV: {n_samples} samples × {n_traits} traits")
                else:
                    warn(f"PRS CSV: only 1 sample, {n_traits} traits")

                # Validate data types
                for r in rows[:5]:
                    try:
                        float(r["prs_raw"])
                        int(r["n_snps"])
                        int(r["n_snps_used"])
                    except (ValueError, TypeError):
                        bad(f"PRS CSV: invalid data type in row {r}")
                        break
                else:
                    ok("PRS CSV: data types valid")

                # Check for NaN
                import math
                nan_rows = [r for r in rows if r["prs_raw"] and math.isnan(float(r["prs_raw"]))]
                if nan_rows:
                    warn(f"PRS CSV: {len(nan_rows)} rows with NaN PRS")
        except Exception as e:
            bad(f"PRS CSV: parse error — {e}")
    else:
        warn("PRS CSV: missing — pipeline not run?")

    # PRS adjusted CSV
    prs_adj = PLATFORM_DIR / "prs" / "prs_adjusted.csv"
    if prs_adj.exists():
        try:
            import csv
            with open(prs_adj) as fh:
                rows = list(csv.DictReader(fh))
            n = len(rows)
            if n > 0:
                ok(f"PRS adjusted: {n} rows")
                if "prs_adjusted" in rows[0]:
                    ok("PRS adjusted: contains prs_adjusted column")
                elif "prs_pca_adjusted" in rows[0]:
                    ok("PRS adjusted: contains prs_pca_adjusted column")
            else:
                warn("PRS adjusted: 0 rows")
        except Exception as e:
            bad(f"PRS adjusted: parse error — {e}")

    # Population calibration
    cal = PLATFORM_DIR / "prs" / "population_calibration_report.json"
    if cal.exists():
        try:
            with open(cal) as fh:
                data = json.load(fh)
            traits = data.get("traits_analyzed", 0)
            method = data.get("methodology", "")
            if traits >= 10:
                ok(f"Calibration: {traits} traits ({method})")
            elif traits > 0:
                warn(f"Calibration: only {traits} traits")
            else:
                warn("Calibration: 0 traits")
        except Exception as e:
            bad(f"Calibration: parse error — {e}")


# ═══════════════════════════════════════════════════════════════════════════════
# TEST 8: ClinVar pathogenic variant annotation
# ═══════════════════════════════════════════════════════════════════════════════

def test_clinvar():
    """Validate ClinVar pathogenic variant output."""
    print("\n── ClinVar ──")

    cv_path = PLATFORM_DIR / "clinvar" / "clinvar_pathogenic_variants.json"
    if not cv_path.exists():
        warn("ClinVar: file missing — run pipeline with --clinvar")
        return

    try:
        with open(cv_path) as fh:
            data = json.load(fh)
    except json.JSONDecodeError as e:
        bad(f"ClinVar: invalid JSON — {e}")
        return

    meta = data.get("metadata", {})
    variants = data.get("pathogenic_variants", [])
    summary = data.get("pathogenic_variant_summary", {})

    if not variants:
        warn("ClinVar: 0 pathogenic variants found — output valid but empty")
        ok("ClinVar: JSON structure valid")
        return

    ok(f"ClinVar: {len(variants)} pathogenic/likely_pathogenic variants")

    # Check required fields in each variant
    required_fields = ["chrom", "pos", "clinical_significance", "disease_name"]
    bad_variants = [v for v in variants if not all(k in v for k in required_fields)]
    if bad_variants:
        bad(f"ClinVar: {len(bad_variants)} variants missing required fields ({required_fields})")
    else:
        ok(f"ClinVar: all {len(variants)} variants have required fields")

    # Check metadata completeness
    n_pos_overlap = meta.get("positional_overlaps", 0)
    n_exact = meta.get("exact_matches", 0)
    n_user = meta.get("user_vcf_total_variants", 0)

    if n_user > 0:
        ok(f"ClinVar: {n_user:,} user variants, {n_pos_overlap:,} ClinVar overlaps, {n_exact:,} exact matches")
    elif n_pos_overlap > 0:
        ok(f"ClinVar: {n_pos_overlap:,} overlaps, {n_exact:,} exact matches")
    else:
        warn("ClinVar: metadata — 0 positional overlaps")

    # Check clinical significance values are valid
    valid_sigs = {"Pathogenic", "Likely_pathogenic", "Pathogenic/Likely_pathogenic", "Risk_allele"}
    sigs_found = set(v.get("clinical_significance", "") for v in variants)
    if sigs_found - valid_sigs:
        warn(f"ClinVar: unexpected CLNSIG values: {sigs_found - valid_sigs}")
    else:
        ok(f"ClinVar: all CLNSIG values valid ({', '.join(sorted(sigs_found))})")

    # Cross-check summary counts match variant list
    expected_total = (
        summary.get("total_pathogenic", 0)
        + summary.get("total_likely_pathogenic", 0)
        + summary.get("total_pathogenic_or_likely_pathogenic", 0)
        + summary.get("total_risk_alleles", 0)
    )
    actual_total = len(variants)
    if expected_total == actual_total:
        ok(f"ClinVar: summary counts consistent ({actual_total})")
    else:
        warn(f"ClinVar: summary counts ({expected_total}) vs actual ({actual_total}) — check classification logic")

    # Check confidence tier fields
    n_tiered = sum(1 for v in variants if "confidence_tier" in v)
    if n_tiered == len(variants):
        ok(f"ClinVar: all {n_tiered} variants have confidence tier")
    else:
        warn(f"ClinVar: {len(variants) - n_tiered} variants missing confidence tier")

    # Check disease description enrichment
    n_desc = sum(1 for v in variants if v.get("disease_description"))
    if n_desc > 0:
        ok(f"ClinVar: {n_desc}/{len(variants)} variants have disease descriptions (MedGen enriched)")


# ═══════════════════════════════════════════════════════════════════════════════
# TEST 9: Pharmacogenomic validation
# ═══════════════════════════════════════════════════════════════════════════════

def test_pharmgkb():
    """Validate pharmacogenomic output."""
    print("\n── PharmGKB ──")

    pgkb_path = PLATFORM_DIR / "pharmgkb" / "pharmgkb_drug_report.json"
    if not pgkb_path.exists():
        warn("PharmGKB: file missing — run pipeline with --clinvar")
        return

    try:
        with open(pgkb_path) as fh:
            data = json.load(fh)
    except json.JSONDecodeError as e:
        bad(f"PharmGKB: invalid JSON — {e}")
        return

    findings = data.get("pharmacogenomic_findings", [])
    summary = data.get("summary", {})

    if not findings:
        warn("PharmGKB: 0 findings — user may not carry known pharmacogenomic variants")
        ok("PharmGKB: JSON structure valid (empty)")
        return

    ok(f"PharmGKB: {len(findings)} drug-response findings")

    # Check required fields
    required = ["gene", "rsid", "drug", "phenotype", "actionability", "recommendation_en"]
    bad_findings = [f for f in findings if not all(k in f for k in required)]
    if bad_findings:
        bad(f"PharmGKB: {len(bad_findings)} findings missing required fields")
    else:
        ok(f"PharmGKB: all findings have required fields ({required})")

    # Check actionability levels are valid
    valid_actions = {"critical", "important", "informative", "normal"}
    actions = set(f.get("actionability", "") for f in findings)
    if actions - valid_actions:
        warn(f"PharmGKB: unexpected actionability values: {actions - valid_actions}")


# ═══════════════════════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════════════════════

def run_all_tests(quick: bool = False):
    print("╔════════════════════════════════════════════╗")
    print("║   PRS PIPELINE TEST SUITE                  ║")
    print("╚════════════════════════════════════════════╝")
    print(f"  Date: {datetime.now().strftime('%Y-%m-%d %H:%M UTC')}")
    print(f"  Root: {PLATFORM_DIR}")

    test_json_integrity()
    test_data_consistency()
    test_pipeline_artifacts()
    test_prs_edge_cases()

    if not quick:
        test_variant_counts()
        test_reports()
        test_code_integrity()
        test_clinvar()
        test_pharmgkb()

    # ── Summary ──
    total = PASS + FAIL + WARN
    print(f"\n{'═'*50}")
    print(f"  RESULTS: {PASS} passed, {WARN} warnings, {FAIL} failed  ({total} total)")
    print(f"{'═'*50}")

    if FAIL == 0 and WARN == 0:
        print(f"\n  ✅ ALL CHECKS PASSED — pipeline healthy")
    elif FAIL == 0:
        print(f"\n  ⚠️  ALL CHECKS PASSED with {WARN} warnings")
    else:
        print(f"\n  ❌ {FAIL} FAILURES — review output above")

    return FAIL


def main():
    import argparse
    p = argparse.ArgumentParser(description="PRS Pipeline Test Suite")
    p.add_argument("--quick", "-q", action="store_true", help="JSON integrity only")
    p.add_argument("--verbose", "-v", action="store_true")
    args = p.parse_args()

    logging.basicConfig(level=logging.DEBUG if args.verbose else logging.INFO,
                       format="%(message)s")

    return run_all_tests(quick=args.quick)


if __name__ == "__main__":
    sys.exit(main())
