#!/usr/bin/env python3
"""
BlueGen — Personal Genomics Platform
Powered by the PRSKit polygenic risk score engine.
Does NOT duplicate scientific logic — router only.

Usage:
  python prs.py run       Full pipeline (VCF → PRS → SSST → reports)
  python prs.py validate  Scientific validation suite
  python prs.py report    Generate bilingual scientific manuscripts
  python prs.py benchmark External benchmarking suite
  python prs.py status    Platform state + output inventory
  python prs.py audit     Export peer review package
"""

import sys, os
# Auto-detect project venv — must happen before non-stdlib imports
_VENV_PYTHON = os.path.join(os.path.dirname(os.path.abspath(__file__)), "venv", "bin", "python3")
if os.path.exists(_VENV_PYTHON):
    _in_venv = any("venv" in p for p in sys.path if "site-packages" in p)
    if not _in_venv:
        _argv0 = sys.argv[0]
        if not os.path.isabs(_argv0):
            _argv0 = os.path.join(os.path.dirname(os.path.abspath(__file__)), os.path.basename(_argv0) or "prs.py")
        os.execv(_VENV_PYTHON, [_VENV_PYTHON, _argv0] + sys.argv[1:])

import json, subprocess, time, shutil
from pathlib import Path
from datetime import datetime, timezone

# ── Configuration ────────────────────────────────────────────────────────────

PROJECT_ROOT = Path(__file__).parent.resolve()
PLATFORM_DIR = PROJECT_ROOT / "prs_research_pipeline"
SCRIPTS = PLATFORM_DIR / "scripts"
PYTHON = sys.executable
START_TIME = time.time()
SNP_DB_PATH = str(PLATFORM_DIR / "data" / "snp_database_annotated.csv")

# Ensure we're in the project root for relative path resolution
os.chdir(str(PROJECT_ROOT))

# ── Terminal ──────────────────────────────────────────────────────────────────

G = '\033[0;32m'; Y = '\033[1;33m'; R = '\033[0;31m'; C = '\033[0;36m'
B = '\033[1m'; D = '\033[2m'; N = '\033[0m'

def ok(s):     print(f"  {G}✓{N} {s}")
def warn(s):   print(f"  {Y}⚠{N}  {s}")
def err(s):    print(f"  {R}✗{N} {s}")
def hdr(s):    print(f"\n{B}{C}── {s} ──{N}")
def info(s):   print(f"  {D}{s}{N}")
def big(s):    print(f"\n{B}{C}{'═'*60}{N}\n{B}{C}  {s}{N}\n{B}{C}{'═'*60}{N}")
def elapsed(): return f"{time.time() - START_TIME:.0f}s"

# ── Command registry ──────────────────────────────────────────────────────────

COMMANDS = {
    "run":      "Full pipeline: VCF → PRS → SSST → reports",
    "validate": "Scientific validation: leakage, validator, adversarial stress",
    "report":   "Generate bilingual scientific manuscripts (EN/ES)",
    "benchmark":"External benchmarking: PGS, GWAS, portability, methods",
    "status":   "Display platform state and output inventory",
    "audit":    "Export audit package for external peer review",
    "test":    "Run test suite — validates JSON integrity and data consistency",
    "pdf":     "Generate PDF report from HTML (WeasyPrint)",
}

ROUTES = {
    # ── Stages (genotype processing) ──
    "stage_a": "stages/01_vcf_to_plink.sh",
    "stage_b": "stages/02_quality_control.sh",
    "stage_c": "stages/03_ld_ancestry_prune.sh",
    "stage_d": "stages/04_pca_1000G.sh",
    "ancestry_classifier": "stages/pca_ancestry_classifier.py",
    # ── PRS computation ──
    "stage_f": "prs/06_prs_compute.py",
    "prs_plink_score": "prs/prs_plink_score.py",
    "stage_g": "prs/pca_adjust_v2.py",
    "stage_h": "prs/population_calibrate_v2.py",
    "stage_i": "validation/13_gwas_ld_consistency_check.py",
    "stage_j": "validation/14_uncertainty_propagation.py",
    "stage_k": "validation/bilingual_interpretation.py",
    # ── Legacy modules (still referenced) ──
    "pca_true": "prs/pca_true_projection.py",
    "ancestry_v2": "stages/ancestry_inference_v2.py",
    "admixture_v2": "stages/admixture_engine_v2.py",
    # ── Validation (Phase 7-8) ──
    "reproducibility": "validation/16_reproducibility_engine.py",
    "scientific_lock": "validation/18_scientific_lock.py",
    "snp_universe": "validation/30_snp_universe_registry.py",
    "leakage_prevention": "validation/31_leakage_prevention.py",
    "global_validator": "validation/32_global_scientific_validator.py",
    "integrity_score": "validation/34_scientific_integrity_score.py",
    "stability_test": "validation/20_stability_repro_test.py",
    "audit_export": "validation/21_audit_exporter.py",
    "methods_gen": "validation/22_methods_generator.py",
    # ── SSST (Phase 9) ──
    "prs_core": "sss/36_prs_core.py",
    "prs_result": "sss/37_prs_result_unified.py",
    "benchmark_reinterp": "sss/38_benchmark_reinterpretation.py",
    "ancestry_unified": "sss/39_ancestry_model_unified.py",
    "leakage_integrated": "sss/40_leakage_integrated.py",
    "report_engine": "sss/41_unified_report_engine.py",
    "consolidation": "sss/42_consolidation_manifest.py",
    # ── Benchmarking ──
    "pgs_benchmark": "benchmarking/23_pgs_catalog_benchmark.py",
    "method_replication": "benchmarking/24_external_prs_replication.py",
    "gwas_consortium": "benchmarking/25_gwas_consortium_validation.py",
    "portability": "benchmarking/26_population_portability_test.py",
    "calibration_validate": "benchmarking/27_real_world_calibration.py",
    "positioning": "benchmarking/28_scientific_positioning.py",
    "quality_delta": "benchmarking/29_quality_delta_analysis.py",
    "pgs_integration": "benchmarking/pgs_catalog_integration.py",
    "gwas_stats": "benchmarking/gwas_summary_stats.py",
    # ── Publication (Phase 10) ──
    "adversarial": "publication/43_adversarial_prs_validation.py",
    "failure_map": "publication/44_failure_mode_map.py",
    "evidence_pack": "publication/45_publication_evidence_pack.py",
    "final_score": "publication/46_final_scientific_score.py",
    "publication_lock": "publication/47_publication_lock.py",
    # ── Reports & Utilities ──
    "comprehensive_report": "publication/comprehensive_report.py",
    "test_suite": "utils/test_suite.py",
    "build_ref_dists": "utils/build_reference_distributions.py",
    "ancestry_normalize": "prs/33_ancestry_aware_normalization.py",
    # ── ClinVar annotation ──
    "download_clinvar": "setup/download_clinvar.py",
    "clinvar_annotate": "clinical/clinvar_annotator.py",
    "medgen_enrich": "clinical/medgen_enrich.py",
    "pharmgkb_annotate": "clinical/pharmgkb_annotator.py",
    "clinpgx_sync": "clinical/clinpgx_sync.py",
    "ancestry_deep": "clinical/ancestry_deep.py",
}

# ═══════════════════════════════════════════════════════════════════════════════
# HELPERS
# ═══════════════════════════════════════════════════════════════════════════════

DRY_RUN = False
DEBUG = False

def run_script(key, *args, shell=False, required=False):
    """Execute a routed script. Returns exit code. Skips if dry-run."""
    script = ROUTES.get(key, key)
    script_path = SCRIPTS / script

    if DRY_RUN:
        info(f"[dry-run] {key} → {script}")
        return 0

    if not script_path.exists():
        if required:
            err(f"Required script not found: {script}")
            return 1
        info(f"Skipping {key} — {script} not found")
        return 0

    cmd = [PYTHON, str(script_path)] + list(args)
    if shell:
        cmd = ["bash", str(script_path)] + list(args)

    try:
        # Debug mode: live output. Normal mode: capture everything, show only stage-level status.
        capture = not DEBUG
        t0 = time.time()

        result = subprocess.run(cmd, cwd=str(PLATFORM_DIR),
                               capture_output=capture, text=True, timeout=7200)

        elapsed = time.time() - t0
        elapsed_str = f" ({elapsed:.0f}s)" if elapsed > 5 else ""

        if result.returncode == 0:
            info(f"  {G}✓{N} {key}{elapsed_str}")
        elif result.returncode == 77:
            warn(f"  {key} — leakage gate blocked (non-fatal)")
        else:
            err(f"  {key} — exit code {result.returncode}")
            # Always write full output to debug log
            log_path = PLATFORM_DIR / "pipeline_debug.log"
            with open(log_path, "a") as lf:
                lf.write(f"\n{'='*60}\n")
                lf.write(f"STAGE: {key} (exit {result.returncode})\n")
                lf.write(f"CMD: {' '.join(cmd)}\n")
                lf.write(f"{'='*60}\n")
                lf.write(result.stdout or "")
                lf.write(result.stderr or "")
                lf.write(f"\n{'='*60}\n")
            output = (result.stdout or "") + (result.stderr or "")
            lines = [l.strip() for l in output.split("\n") if l.strip()]
            if capture:
                info(f"    {D}Full log: pipeline_debug.log{N}")
                for line in lines[-8:]:
                    info(f"    {D}{line[:130]}{N}")
            if required:
                err(f"  Pipeline halted: {key} failed with exit code {result.returncode}")
        return result.returncode
    except subprocess.TimeoutExpired:
        err(f"  {key} — timed out")
        return 1


def run_shell(key, *args, required=False):
    return run_script(key, *args, shell=True, required=required)


def require_output(path, label, stage=""):
    """Validate that a pipeline output file exists. Returns True or halts."""
    if exists(path):
        return True
    err(f"Output not found: {path} ({label})")
    if stage:
        err(f"  Stage '{stage}' did not produce expected output")
    err(f"  Check pipeline_debug.log for details")
    return False


def exists(path):
    return (PLATFORM_DIR / path).exists()


def load_json(path):
    p = PLATFORM_DIR / path
    if p.exists():
        with open(p) as fh: return json.load(fh)
    return {}


def check(path, label):
    """Check and report on a file."""
    if exists(path):
        ok(f"{label}")
        return True
    return False


# ═══════════════════════════════════════════════════════════════════════════════
# COMMAND: run
# ═══════════════════════════════════════════════════════════════════════════════

def cmd_run(args):
    """Full pipeline execution."""
    global DRY_RUN, DEBUG
    DRY_RUN = args.dry_run
    DEBUG = args.debug
    full = args.full
    vcf = args.vcf or "input.vcf.gz"
    sample = args.sample or "SAMPLE_001"
    lang = args.lang or "both"
    clinvar = args.clinvar or full
    update_refs = args.update_references

    # Support multiple VCFs (comma-separated). Merge with bcftools if >1.
    vcf_files = [os.path.abspath(f.strip()) for f in vcf.split(",")]
    if len(vcf_files) > 1:
        info(f"Merging {len(vcf_files)} VCFs...")
        merged_vcf = str(PLATFORM_DIR / "input_merged.vcf.gz")
        # Use bcftools merge for multi-sample VCFs
        cmd = ["bcftools", "merge", "-Oz", "-o", merged_vcf] + vcf_files
        subprocess.run(cmd, check=True, timeout=600)
        subprocess.run(["bcftools", "index", "-t", merged_vcf], check=True, timeout=60)
        vcf = merged_vcf
        n_samples = int(subprocess.check_output(
            ["bcftools", "query", "-l", vcf], text=True).count('\n'))
        info(f"  Merged: {n_samples} samples")
    else:
        vcf = vcf_files[0]
        if not os.path.exists(vcf) and not DRY_RUN:
            err(f"Input VCF required: {vcf}")
            err("Provide path: python prs.py run --vcf your_sample.vcf.gz")
            err("  Or multiple: python prs.py run --vcf s1.vcf.gz,s2.vcf.gz")
            return 1

    big(f"PRS PIPELINE {'(dry-run)' if DRY_RUN else ''}")
    info(f"Sample: {sample} | VCF: {vcf} | Language: {lang}")
    info(f"Mode: {'FULL' if full else 'PIPELINE ONLY'}{' — Stage: '+args.stage if args.stage else ''}")

    # Validate config
    if not DRY_RUN:
        try:
            sys.path.insert(0, str(PLATFORM_DIR / "scripts"))
            from utils.config_validator import validate_config
            validate_config(str(PLATFORM_DIR / "config.yaml"))
            ok("Config validated")
        except SystemExit:
            return 1

    # ── Single stage mode ──
    if args.stage:
        stage = args.stage.lower()
        user_vcf_abs = os.path.abspath(vcf)

        stage_map = {
            "clinvar": lambda: [
                hdr("ClinVar Annotation"),
                run_script("download_clinvar", "--output-dir", "reference/clinvar") if not exists("reference/clinvar/clinvar.vcf.gz") else None,
                run_script("clinvar_annotate", "--vcf", user_vcf_abs, "--clinvar-vcf", str(PLATFORM_DIR / "reference/clinvar/clinvar.vcf.gz"), "--output-dir", "clinvar/"),
            ],
            "medgen": lambda: [
                hdr("MedGen Enrichment"),
                run_script("medgen_enrich", "--download", "--ref-dir", "reference/medgen") if not exists("reference/medgen/NAMES.RRF.gz") else None,
                run_script("medgen_enrich", "--clinvar-json", str(PLATFORM_DIR / "clinvar/clinvar_pathogenic_variants.json"), "--ref-dir", "reference/medgen"),
            ],
            "pharmgkb": lambda: [
                hdr("Pharmacogenomics"),
                run_script("clinpgx_sync", "--sync", "--output-dir", "reference/clinpgx") if not exists("reference/clinpgx/clinicalVariants.zip") else None,
                run_script("pharmgkb_annotate", "--vcf", user_vcf_abs, "--output-dir", "pharmgkb/"),
            ],
            "ancestry": lambda: [
                hdr("Deep Ancestry"),
                run_script("ancestry_deep", "--vcf", user_vcf_abs, "--output-dir", "ancestry/"),
            ],
            "clinpgx": lambda: [
                hdr("ClinPGx Sync"),
                run_script("clinpgx_sync", "--sync", "--output-dir", "reference/clinpgx"),
            ],
        }

        if stage in stage_map:
            stage_map[stage]()
            info(f"\n  ✅ Stage '{stage}' complete")
            return 0
        else:
            err(f"Unknown stage: '{stage}'")
            info(f"Available stages: {', '.join(sorted(stage_map.keys()))}")
            return 1

    # Rotate debug log at start of each run
    log_path = PLATFORM_DIR / "pipeline_debug.log"
    if not DRY_RUN:
        if log_path.exists():
            log_size = log_path.stat().st_size
            if log_size > 1024 * 1024:  # >1 MB — archive old log
                archive_path = PLATFORM_DIR / f"pipeline_debug_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"
                log_path.rename(archive_path)
                info(f"  Debug log archived ({log_size / (1024*1024):.0f} MB)")
            else:
                log_path.write_text("")  # Truncate for fresh run
        # Write run header
        log_path.write_text(f"=== BlueGen Pipeline v1.0.0 — {datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M UTC')} ===\n"
                           f"Sample: {sample} | VCF: {vcf} | Mode: {'FULL' if full else 'PIPELINE'}\n\n")

    # Pre-flight checks
    if not DRY_RUN:
        # Add tools/ and venv/ to PATH so scripts find PLINK + Python deps
        tools_dir = str(Path(__file__).parent / "tools")
        venv_bin = str(Path(__file__).parent / "venv" / "bin")
        os.environ["PATH"] = f"{tools_dir}:{venv_bin}:" + os.environ.get("PATH", "")

        try:
            sys.path.insert(0, str(PLATFORM_DIR / "scripts"))
            from utils.tool_detection import find_plink
            plink_path, plink_ver = find_plink()
            ok(f"PLINK: {plink_path} ({plink_ver})")
        except SystemExit:
            return 1

    snp_db = SNP_DB_PATH

    # Resolve 1000G reference paths once
    g1k_full_bfile = str(PLATFORM_DIR / "reference/1000G_full/1000G_full")
    g1k_vcf = str(PLATFORM_DIR / "reference/1000G/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz")
    pop_panel_full = str(PLATFORM_DIR / "reference/1000G_full/population_panel.txt")
    pop_panel = str(PLATFORM_DIR / "reference/1000G/20130606_g1k_3202_samples_ped_population.txt")

    # Use full genome-wide reference if available, otherwise chr22
    use_full_ref = Path(g1k_full_bfile + ".bed").exists()
    g1k_ref = g1k_full_bfile if use_full_ref else g1k_vcf
    g1k_arg = "--1000g-bfile" if use_full_ref else "--1000g-vcf"
    pop_ref = pop_panel_full if use_full_ref and Path(pop_panel_full).exists() else pop_panel

    if use_full_ref:
        ok(f"1000G Reference: genome-wide ({g1k_full_bfile})")
    else:
        warn(f"1000G Reference: chr22 only — run scripts/setup/download_1000G_full.py for genome-wide")

    # ── Genotype Processing ──
    hdr("Genotype Processing (Stages A–D)")
    run_shell("stage_a", "--input-vcf", vcf, "--out-dir", "plink/", required=True)
    require_output("plink/cohort.bed", "PLINK binary dataset", "stage_a")

    run_shell("stage_b", "--bfile", "plink/cohort", "--out-dir", "qc/", required=True)
    require_output("qc/qc_filtered.bed", "QC-filtered dataset", "stage_b")

    run_shell("stage_c", "--bfile", "qc/qc_filtered",
              g1k_arg, g1k_ref,
              "--population-panel", pop_ref,
              "--out-dir", "plink/", required=True)
    require_output("plink/ld_pruned_dataset.bed", "LD-pruned dataset", "stage_c")

    run_shell("stage_d", g1k_arg, g1k_ref,
              "--target-bfile", "plink/ld_pruned_dataset",
              "--population-panel", pop_ref, "--out-dir", "pca/")
    # Genome-wide PCA ancestry classifier (real ancestry, not allele-frequency hack)
    run_script("ancestry_classifier", "--ref-pcs", "pca/1000G_pcs.eigenvec",
               "--target-pcs", "pca/target_pcs.eigenvec",
               "--pop-panel", pop_ref, "--output-dir", "pca/")

    # Build population reference distributions once (cached after first run)
    ref_dist = "reference/population_distributions/reference_distributions.json"
    if use_full_ref and not exists(ref_dist):
        info(f"  {D}Building population reference distributions (one-time)...{N}")
        run_script("build_ref_dists", "--bfile", g1k_full_bfile,
                   "--snp-db", snp_db, "--pop-panel", pop_ref,
                   "--output-dir", "reference/population_distributions")

    # Deep ancestry (mtDNA + Y-DNA + Neanderthal)
    run_script("ancestry_deep", "--vcf", os.path.abspath(vcf),
               "--output-dir", "ancestry/")

    # ── PRS Computation ──
    hdr("PRS Computation (Stages F–H)")
    if not DRY_RUN:
        plink_bin = str(Path(__file__).parent / "tools" / "plink")
        rc = run_script("prs_plink_score",
                   "--snp-db", snp_db,
                   "--bfile", "qc/qc_filtered",
                   "--output-dir", "prs/",
                   "--plink", plink_bin)
        if rc != 0:
            err("PRS computation failed — check pipeline_debug.log")
            return 1
        require_output("prs/prs_raw.csv", "Raw PRS scores", "prs_plink_score")

    # PCA adjustment — skip if no PCA output
    if exists("prs/prs_raw.csv") and exists("pca/target_pcs.eigenvec"):
        run_script("stage_g", "--prs-data", "prs/prs_raw.csv",
                   "--sample-pcs", "pca/target_pcs.eigenvec",
                   "--output-dir", "prs/", "--sample-id", sample)
        require_output("prs/pca_adjusted_scores.csv", "PCA-adjusted PRS", "stage_g")
    elif exists("prs/prs_raw.csv"):
        info("  Skipping PCA adjustment — no target_pcs.eigenvec (run Stage D first)")

    # Population calibration
    anc_json = "ancestry/classification_report.json"
    if not exists(anc_json):
        anc_json = "pca/ancestry_inference.json"
    cal_done = False
    if exists(anc_json) and exists("prs/pca_adjusted_scores.csv"):
        run_script("stage_h", "--sample-prs", "prs/pca_adjusted_scores.csv",
                   "--ancestry-json", anc_json, "--output-dir", "prs/",
                   "--calibrate-only")
        cal_done = True
    elif exists(anc_json) and exists("prs/prs_raw.csv"):
        run_script("stage_h", "--sample-prs", "prs/prs_raw.csv",
                   "--ancestry-json", anc_json, "--output-dir", "prs/",
                   "--calibrate-only")
        cal_done = True
    if not cal_done and exists("prs/prs_raw.csv"):
        # Fallback: calibrate without ancestry info
        run_script("stage_h", "--sample-prs", "prs/prs_raw.csv",
                   "--output-dir", "prs/", "--calibrate-only")

    # Uncertainty + consistency + interpretation
    prs_cal = "prs/population_calibrated_v2.csv"
    if not exists(prs_cal):
        prs_cal = "prs/population_calibrated.csv"
    if exists(prs_cal) and exists(anc_json):
        run_script("stage_i", "--ancestry", anc_json, "--curated-db", snp_db,
                   "--output-dir", "prs/", "--no-strict")
        run_script("stage_j", "--prs-data", prs_cal, "--ancestry", anc_json,
                   "--snp-db", snp_db, "--vcf", vcf, "--output-dir", "prs/")
        run_script("stage_k", "--prs-calibrated", prs_cal, "--ancestry", anc_json,
                   "--output-dir", "interpretations/", "--lang", lang)
    else:
        if not exists(prs_cal):
            warn("No calibrated PRS — skipping uncertainty and interpretation stages")
        if not exists(anc_json):
            warn("No ancestry data — skipping uncertainty and interpretation stages")

    # ── ClinVar Pathogenic Variant Annotation ──
    if clinvar:
        hdr("ClinVar Pathogenic Variant Annotation")
        clinvar_ref = "reference/clinvar/clinvar.vcf.gz"
        if not exists(clinvar_ref):
            run_script("download_clinvar", "--output-dir", "reference/clinvar")
        user_vcf_abs = os.path.abspath(vcf)
        run_script("clinvar_annotate", "--vcf", user_vcf_abs,
                   "--clinvar-vcf", str(PLATFORM_DIR / clinvar_ref),
                   "--output-dir", "clinvar/")
        # Enrich with disease descriptions from MedGen local database
        medgen_ref = "reference/medgen"
        if update_refs or not exists(medgen_ref + "/NAMES.RRF.gz"):
            medgen_args = ["--download", "--ref-dir", medgen_ref]
            if update_refs:
                medgen_args.append("--force")
            run_script("medgen_enrich", *medgen_args)
        run_script("medgen_enrich", "--clinvar-json",
                   str(PLATFORM_DIR / "clinvar/clinvar_pathogenic_variants.json"),
                   "--ref-dir", medgen_ref)

        # Pharmacogenomic annotation (drug response variants)
        # One-time download of ClinPGx datasets (~1.8 MB total, free, no login)
        if update_refs or not exists("reference/clinpgx/clinicalVariants.zip"):
            run_script("clinpgx_sync", "--sync", "--output-dir", "reference/clinpgx")
        run_script("pharmgkb_annotate", "--vcf", user_vcf_abs,
                   "--output-dir", "pharmgkb/")

    # ── Freeze Layer ──
    hdr("Scientific Freeze (Phase 7)")
    run_script("reproducibility", "--lock", "--fingerprint", "--output-dir", "reproducibility")
    run_script("scientific_lock", "--import-config", "config.yaml", "--output-dir", "science")

    # ── Correction Layers ──
    hdr("Correction Layers (Phase 8)")
    run_script("snp_universe", "--snp-db", snp_db, "--output-dir", "science")
    run_script("leakage_prevention", "--sample-id", sample, "--output-dir", "science")
    run_script("global_validator", "--snp-db", snp_db, "--output-dir", "science")

    # ── SSST Consolidation ──
    hdr("SSST Consolidation (Phase 9)")
    run_script("prs_core", "--snp-db", snp_db, "--output-dir", "science", "--validate-all")
    run_script("ancestry_unified", "--output-dir", "science")
    run_script("prs_result", "--sample-id", sample, "--output-dir", "prs")
    run_script("benchmark_reinterp", "--output-dir", "benchmark")
    run_script("leakage_integrated", "--sample-id", sample, "--output-dir", "science", "--no-hard-stop")
    run_script("consolidation", "--output-dir", "science")

    # ── FULL MODE ONLY ──
    if full:
        hdr("Validation Suite")
        run_script("adversarial", "--output-dir", "science")
        run_script("failure_map", "--output-dir", "science")

        hdr("Benchmarking Suite")
        if exists("prs/pgs_scores.csv"):
            run_script("pgs_benchmark", "--internal-prs", prs_cal,
                       "--pgs-scores", "prs/pgs_scores.csv", "--output-dir", "benchmark")
        run_script("gwas_consortium", "--snp-db", snp_db, "--output-dir", "benchmark")
        run_script("portability", "--output-dir", "benchmark")
        # GWAS summary stats via tabix (no auth needed — remote VCF queries)
        run_script("gwas_stats", "--snp-db", snp_db, "--output-dir", "gwas")
        if exists(prs_cal):
            run_script("calibration_validate", "--calibrated-prs", prs_cal,
                       "--output-dir", "benchmark")
        run_script("positioning", "--output-dir", "benchmark")
        run_script("quality_delta", "--output-dir", "benchmark")
        run_script("pgs_integration", "--bfile", "qc/qc_filtered",
                   "--output-dir", "pgs")

        hdr("Reports")
        run_script("report_engine", "--sample-id", sample,
                   "--prs-core", "science/prs_core_definition.json",
                   "--ancestry", "science/ANCESTRY_MODEL.json",
                   "--prs-result", "prs/PRS_RESULT.json",
                   "--benchmark", "benchmark/VALIDATION_REPORT.json",
                   "--integrity", "science/scientific_integrity_score.json",
                   "--output-dir", "reports/")
        run_script("comprehensive_report", "--sample-id", sample, "--lang", lang,
                   "--output-dir", "reports/")

        hdr("Publication Lock")
        run_script("evidence_pack", "--output-dir", "publication_evidence_pack")
        run_script("final_score")
        run_script("publication_lock")

    # ── Summary ──
    print(f"\n{B}{C}{'═'*50}{N}")
    print(f"{B}  PIPELINE COMPLETE{N}  ({elapsed()})")
    print(f"{B}{'═'*50}{N}\n")

    deliverables = [
        ("PRS_RESULT.json", "prs/PRS_RESULT.json", "Unified PRS output"),
        ("ANCESTRY_MODEL.json", "science/ANCESTRY_MODEL.json", "Ancestry classification"),
        ("PRS_CORE", "science/prs_core_definition.json", "Canonical PRS definition"),
        ("Deep Ancestry", "ancestry/deep_ancestry.json", "mtDNA/Y-DNA haplogroups"),
    ]
    if full:
        deliverables += [
            ("Validation Report", "science/global_validation_report.json", "8-dimension check"),
            ("EN Manuscript", "reports/SCIENTIFIC_MANUSCRIPT_EN.md", "English manuscript"),
            ("ES Manuscript", "reports/SCIENTIFIC_MANUSCRIPT_ES.md", "Spanish manuscript"),
            ("EN Report (HTML)", "reports/comprehensive_report_en.html", "Interactive HTML report"),
            ("ES Report (HTML)", "reports/comprehensive_report_es.html", "Interactive HTML report"),
            ("Failure Map", "science/failure_mode_map.json", "18 failure modes"),
            ("Final Score", "FINAL_SCIENTIFIC_SCORE.json", "Locked integrity score"),
            ("Publication Lock", "PUBLICATION_LOCK.md", "Readiness declaration"),
        ]

    if clinvar:
        deliverables += [
            ("ClinVar Pathogenic", "clinvar/clinvar_pathogenic_variants.json",
             "Pathogenic variant annotation"),
            ("PharmGKB Drugs", "pharmgkb/pharmgkb_drug_report.json",
             "Pharmacogenomic drug response"),
        ]

    for name, path, desc in deliverables:
        icon = f"{G}✓{N}" if exists(path) else f"{D}⬚{N}"
        print(f"  {icon} {name}  {D}→ {desc}{N}")

    if DRY_RUN:
        info("\nThis was a dry-run. No files were created. Remove --dry-run to execute.")
    else:
        info(f"\nAll outputs in: {PLATFORM_DIR}/")
        info(f"Reports:")
        info(f"  HTML: open {PLATFORM_DIR}/reports/comprehensive_report_en.html")
        info(f"  PDF:  open {PLATFORM_DIR}/reports/comprehensive_report_en.pdf")
        info(f"  MD:   open {PLATFORM_DIR}/reports/SCIENTIFIC_MANUSCRIPT_EN.md")
        info(f"Dashboard: cd {Path(__file__).parent} && streamlit run dashboard.py")

    return 0


# ═══════════════════════════════════════════════════════════════════════════════
# COMMAND: validate
# ═══════════════════════════════════════════════════════════════════════════════

def cmd_validate(args):
    big("SCIENTIFIC VALIDATION")
    sample = args.sample or "SAMPLE_001"
    snp_db = SNP_DB_PATH

    stages = [
        ("Leakage Prevention", "leakage_prevention", ["--sample-id", sample, "--output-dir", "science"]),
        ("Integrated Gates", "leakage_integrated", ["--sample-id", sample, "--output-dir", "science", "--no-hard-stop"]),
        ("Global Validator (8 dims)", "global_validator", ["--snp-db", snp_db, "--output-dir", "science"]),
        ("Genome Coverage Check", "snp_universe", ["--snp-db", snp_db, "--output-dir", "science"]),
        ("Adversarial Stress Tests", "adversarial", ["--output-dir", "science"]),
        ("Failure Mode Map", "failure_map", ["--output-dir", "science"]),
    ]

    for label, key, script_args in stages:
        hdr(label)
        run_script(key, *script_args)

    # Stability
    cal = "prs/population_calibrated_v2.csv"
    if not exists(cal):
        cal = "prs/population_calibrated.csv"
    if exists(cal):
        hdr("Stability Analysis")
        run_script("stability_test", "--prs-file", cal, "--output-dir", "validation")

    hdr("Final Score")
    run_script("final_score")

    # Print result
    val = load_json("science/global_validation_report.json")
    print(f"\n{B}{'─'*40}{N}")
    if val:
        score = val.get("overall_score", 0)
        status = val.get("overall_status", "Unknown")
        icon = "✓" if score >= 75 else ("⚠" if score >= 60 else "✗")
        print(f"  {G if score >= 75 else Y if score >= 60 else R}{icon}{N}  Score: {B}{score:.0f}/100{N} — {status}")
        for c in val.get("checks", []):
            ico = "✓" if c.get("passed") else "✗"
            print(f"    {G if c.get('passed') else R}{ico}{N} {c.get('description','?')[:60]}")

    print(f"\n  {elapsed()} — validation complete")
    return 0


# ═══════════════════════════════════════════════════════════════════════════════
# COMMAND: report
# ═══════════════════════════════════════════════════════════════════════════════

def cmd_report(args):
    big("BILINGUAL SCIENTIFIC MANUSCRIPTS")
    sample = args.sample or "SAMPLE_001"

    lang = args.lang or "both"
    run_script("report_engine", "--sample-id", sample,
               "--prs-core", "science/prs_core_definition.json",
               "--ancestry", "science/ANCESTRY_MODEL.json",
               "--prs-result", "prs/PRS_RESULT.json",
               "--benchmark", "benchmark/VALIDATION_REPORT.json",
               "--integrity", "science/scientific_integrity_score.json",
               "--output-dir", "reports/")
    run_script("comprehensive_report", "--sample-id", sample, "--lang", lang,
               "--output-dir", "reports/")
    run_script("methods_gen", "--output-dir", "science")
    run_script("evidence_pack", "--output-dir", "publication_evidence_pack")

    print(f"\n{B}Generated Reports:{N}")
    for lang_code, path in [("EN", "reports/comprehensive_report_en.html"),
                            ("ES", "reports/comprehensive_report_es.html")]:
        if exists(path):
            size = os.path.getsize(str(PLATFORM_DIR / path))
            ok(f"[{lang_code}] HTML: {path} ({size/1024:.0f} KB) → open {PLATFORM_DIR / path}")
    for lang_code, path in [("EN", "reports/comprehensive_report_en.pdf"),
                            ("ES", "reports/comprehensive_report_es.pdf")]:
        if exists(path):
            size = os.path.getsize(str(PLATFORM_DIR / path))
            ok(f"[{lang_code}] PDF:  {path} ({size/1024:.0f} KB)")
    for lang_code, path in [("EN", "reports/SCIENTIFIC_MANUSCRIPT_EN.md"),
                            ("ES", "reports/SCIENTIFIC_MANUSCRIPT_ES.md")]:
        if exists(path):
            lines = sum(1 for _ in open(str(PLATFORM_DIR / path)))
            ok(f"[{lang_code}] MD:   {path} ({lines} lines)")
    check("science/methods_section.md", "Methods section")
    check("science/supplementary_methods.md", "Supplementary methods")
    check("publication_evidence_pack/limitations.md", "Limitations")
    check("publication_evidence_pack/consistency_report.json", "Consistency report")

    return 0


# ═══════════════════════════════════════════════════════════════════════════════
# COMMAND: benchmark
# ═══════════════════════════════════════════════════════════════════════════════

def cmd_benchmark(args):
    big("EXTERNAL BENCHMARKING")
    snp_db = SNP_DB_PATH
    prs_cal = "prs/population_calibrated_v2.csv"
    if not exists(prs_cal):
        prs_cal = "prs/population_calibrated.csv"

    benchmarks = [
        ("PGS Catalog Comparison", "pgs_benchmark",
         ["--internal-prs", prs_cal, "--pgs-scores", "prs/pgs_scores.csv", "--output-dir", "benchmark"],
         exists("prs/pgs_scores.csv")),
        ("GWAS Consortium Validation", "gwas_consortium",
         ["--snp-db", snp_db, "--output-dir", "benchmark"], True),
        ("Method Replication", "method_replication",
         ["--output-dir", "benchmark"], True),
        ("Population Portability", "portability",
         ["--output-dir", "benchmark"], True),
        ("Calibration Validation", "calibration_validate",
         ["--calibrated-prs", prs_cal, "--output-dir", "benchmark"], exists(prs_cal)),
        ("Scientific Positioning", "positioning",
         ["--output-dir", "benchmark"], True),
        ("Quality Delta Analysis", "quality_delta",
         ["--output-dir", "benchmark"], True),
    ]

    for label, key, script_args, should_run in benchmarks:
        hdr(label)
        if should_run:
            run_script(key, *script_args)
        else:
            info(f"Skipping — missing input data")

    run_script("benchmark_reinterp", "--output-dir", "benchmark")

    print(f"\n{B}Benchmark Outputs:{N}")
    for f in ["pgs_comparison.json", "gwas_consortium_validation.json",
              "portability_report.json", "calibration_validation.json",
              "scientific_positioning.md", "quality_delta.json", "VALIDATION_REPORT.json"]:
        check(f"benchmark/{f}", f)

    return 0


# ═══════════════════════════════════════════════════════════════════════════════
# COMMAND: status
# ═══════════════════════════════════════════════════════════════════════════════

def cmd_status(args):
    """Enhanced status dashboard."""
    print(f"\n{B}🧬 BlueGen v1.0.0{N}\n")

    # Environment
    print(f"  {B}Environment{N}")
    print(f"    Python: {sys.version.split()[0]}")
    print(f"    Platform: {PLATFORM_DIR}")

    # Check PLINK
    try:
        sys.path.insert(0, str(PLATFORM_DIR / "scripts"))
        from utils.tool_detection import find_plink
        plink_path, plink_ver = find_plink()
        print(f"    PLINK: {plink_ver}")
    except SystemExit:
        print(f"    PLINK: {R}not found{N}")

    # Pipeline state
    print(f"\n  {B}Pipeline State{N}")

    # Determine state from key outputs
    has_prs = exists("prs/PRS_RESULT.json") or exists("prs/prs_raw.csv")
    has_ancestry = exists("science/ANCESTRY_MODEL.json") or exists("pca/ancestry_inference.json")
    has_validation = exists("science/global_validation_report.json")
    has_benchmark = exists("benchmark/VALIDATION_REPORT.json")
    has_reports = exists("reports/SCIENTIFIC_MANUSCRIPT_EN.md")
    has_lock = exists("PUBLICATION_LOCK.md")

    if has_lock:
        print(f"    {G}⬤{N} PIPELINE COMPLETE — Publication-locked")
    elif has_reports and has_validation:
        print(f"    {G}⬤{N} PIPELINE COMPLETE — Reports generated")
    elif has_prs:
        print(f"    {Y}◐{N} PARTIAL — PRS computed, validation pending")
    else:
        print(f"    {D}○{N} IDLE — Run `python prs.py run --full` to start")

    # Timestamp
    fingerprint = load_json("reproducibility/run_fingerprint.json")
    if fingerprint:
        ts = fingerprint.get("timestamp_utc", "")[:19]
        print(f"    Last run: {ts}")

    # SSST status
    cons = load_json("science/CONSOLIDATION_MANIFEST.json")
    if cons:
        s = cons.get("summary", {})
        active = s.get("active", 0); total = s.get("total", 7)
        complete = s.get("consolidation_complete", False)
        icon = f"{G}✓{N}" if complete else f"{Y}⚠{N}"
        print(f"\n  {B}SSST{N}  {icon} {active}/{total} sources active — "
              f"{'Consolidated' if complete else 'Pending'}")

    # Scores
    fin = load_json("FINAL_SCIENTIFIC_SCORE.json")
    if fin:
        score = fin.get("scientific_integrity_score", 0)
        cat = fin.get("category", "?")
        color = G if score >= 75 else (Y if score >= 60 else R)
        print(f"  {B}Integrity{N} {color}{score:.0f}/100{N} — {cat}")

    pub_lock = load_json("PUBLICATION_LOCK.json")
    if pub_lock:
        ready = pub_lock.get("publication_ready", False)
        icon = f"{G}✓ READY{N}" if ready else f"{R}✗ NOT READY{N}"
        print(f"  {B}Publication{N} {icon}")

    # Output inventory
    print(f"\n  {B}Output Inventory{N}")
    inventory = [
        ("PRS_RESULT.json", "prs/PRS_RESULT.json"),
        ("ANCESTRY_MODEL.json", "science/ANCESTRY_MODEL.json"),
        ("PRS Core Definition", "science/prs_core_definition.json"),
        ("Validation Report", "science/global_validation_report.json"),
        ("EN Manuscript", "reports/SCIENTIFIC_MANUSCRIPT_EN.md"),
        ("ES Manuscript", "reports/SCIENTIFIC_MANUSCRIPT_ES.md"),
        ("Adversarial Report", "science/adversarial_validation_report.json"),
        ("Failure Map", "science/failure_mode_map.json"),
        ("Benchmarks", "benchmark/VALIDATION_REPORT.json"),
        ("Publication Lock", "PUBLICATION_LOCK.md"),
        ("ClinVar Pathogenic", "clinvar/clinvar_pathogenic_variants.json"),
        ("PharmGKB Drugs", "pharmgkb/pharmgkb_drug_report.json"),
        ("Deep Ancestry", "ancestry/deep_ancestry.json"),
    ]
    for name, path in inventory:
        icon = f"{G}✓{N}" if exists(path) else f"{D}⬚{N}"
        print(f"    {icon} {name}")

    # Quick start hint
    if not has_prs:
        print(f"\n  {Y}💡 Run:{N} python prs.py run --full --vcf your_sample.vcf.gz")

    return 0


# ═══════════════════════════════════════════════════════════════════════════════
# COMMAND: audit
# ═══════════════════════════════════════════════════════════════════════════════

def cmd_audit(args):
    big("EXTERNAL AUDIT EXPORT")
    import hashlib
    run_id = hashlib.sha256(str(datetime.now(timezone.utc)).encode()).hexdigest()[:12]
    sample = args.sample or "SAMPLE_001"

    run_script("audit_export", "--run-id", run_id, "--sample-id", sample,
               "--vcf", args.vcf or "input.vcf.gz",
               "--snp-db", str(PLATFORM_DIR / "data/snp_database_annotated.csv"),
               "--output-dir", "audit/")

    audit_dir = PLATFORM_DIR / "audit"
    if audit_dir.exists():
        for f in sorted(audit_dir.glob("*.zip"), reverse=True):
            ok(f"Package: {f.name} ({f.stat().st_size:,} bytes)")
            break

    return 0


# ═══════════════════════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════════════════════

def main():
    import argparse

    parser = argparse.ArgumentParser(
        description="BlueGen v1.0.0 — Personal Genomics Platform",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""Examples:
  python prs.py run --vcf sample.vcf.gz           # Full pipeline
  python prs.py run --full --vcf sample.vcf.gz    # Pipeline + validation + reports
  python prs.py validate                           # Scientific checks
  python prs.py report --lang both                 # Bilingual manuscripts
  python prs.py benchmark                          # External comparisons
  python prs.py status                             # Platform dashboard
  python prs.py audit                              # Peer review package""")

    parser.add_argument("command", nargs="?", default="status",
                       choices=list(COMMANDS.keys()), help="Command")
    parser.add_argument("--vcf", help="Input VCF path (.vcf.gz)")
    parser.add_argument("--sample", help="Sample identifier")
    parser.add_argument("--lang", default="both", choices=["en", "es", "both"])
    parser.add_argument("--full", action="store_true",
                       help="Run pipeline + validation + benchmarks + reports + lock")
    parser.add_argument("--with-1000g", action="store_true",
                       help="Include 1000 Genomes PCA projection")
    parser.add_argument("--dry-run", action="store_true",
                       help="Preview stages without executing")
    parser.add_argument("--debug", "-d", action="store_true",
                       help="Show full PLINK/shell output (default: quiet with stage summary)")
    parser.add_argument("--seed", type=int, default=42,
                       help="Random seed for reproducibility")
    parser.add_argument("--quick", "-q", action="store_true",
                       help="Test: JSON integrity only (fast)")
    parser.add_argument("--verbose", "-v", action="store_true",
                       help="Show detailed output")
    parser.add_argument("--clinvar", action="store_true",
                       help="Include ClinVar pathogenic variant annotation")
    parser.add_argument("--update-references", action="store_true",
                       help="Force update of reference databases (ClinVar, MedGen)")
    parser.add_argument("--stage", help="Run a single stage (e.g., clinvar, pharmgkb, ancestry, medgen, clinpgx)")

    args = parser.parse_args()

    # Global settings
    global DEBUG
    DEBUG = getattr(args, "debug", False)

    # Set seed if specified
    if args.seed != 42:
        import numpy as np; np.random.seed(args.seed)
        os.environ["PYTHONHASHSEED"] = str(args.seed)

    handlers = {
        "run": cmd_run, "validate": cmd_validate, "report": cmd_report,
        "benchmark": cmd_benchmark, "status": cmd_status, "audit": cmd_audit,
        "test": cmd_test, "pdf": cmd_pdf,
    }
    handler = handlers.get(args.command)
    if handler:
        return handler(args) or 0

    print(f"Unknown command: {args.command}")
    return 1


def cmd_pdf(args):
    """Generate PDF from HTML reports."""
    big("PDF REPORT GENERATION")
    lang = args.lang or "both"
    report_dir = PLATFORM_DIR / "reports"

    for lang_code in (["en", "es"] if lang == "both" else [lang]):
        html_path = report_dir / f"comprehensive_report_{lang_code}.html"
        pdf_path = report_dir / f"comprehensive_report_{lang_code}.pdf"
        if not html_path.exists():
            warn(f"HTML report not found: {html_path}")
            info("Run 'python prs.py report' first")
            continue
        try:
            from weasyprint import HTML
            HTML(filename=str(html_path)).write_pdf(str(pdf_path))
            size_kb = pdf_path.stat().st_size / 1024
            ok(f"[{lang_code.upper()}] {pdf_path.name} ({size_kb:.0f} KB)")
        except Exception as e:
            err(f"[{lang_code.upper()}] PDF failed: {e}")
    return 0


def cmd_test(args):
    """Run test suite."""
    return run_script("test_suite", *(["--quick"] if getattr(args, "quick", False) else []),
                      *(["--verbose"] if getattr(args, "verbose", False) else []))


if __name__ == "__main__":
    sys.exit(main())
