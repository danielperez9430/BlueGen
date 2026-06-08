#!/usr/bin/env bash
# ============================================================================
# PRS RESEARCH PLATFORM v6.0.0 — Interactive Orchestrator
# ============================================================================
# Phase 6: Scientific corrections integrated. Full 1000G, true PCA projection,
# empirical calibration, corrected multi-method PRS, orthogonal RQS v2.
#
# Usage:
#   ./run.sh                                    # Interactive mode
#   ./run.sh --input-vcf file.vcf.gz            # Quick start
#   ./run.sh --input-vcf file.vcf.gz --auto     # Non-interactive full pipeline
#   ./run.sh --input-vcf file.vcf.gz --dry-run  # Preview without executing
#   ./run.sh --help                             # Full help
# ============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_DIR="$SCRIPT_DIR"
cd "$PIPELINE_DIR"

# ═══════════════════════════════════════════════════════════════════════════════
# COLOR & DISPLAY
# ═══════════════════════════════════════════════════════════════════════════════
RED='\033[0;31m'; GREEN='\033[0;32m'; YELLOW='\033[1;33m'; BLUE='\033[0;34m'
CYAN='\033[0;36m'; BOLD='\033[1m'; NC='\033[0m'
ok()  { echo -e "  ${GREEN}✅${NC} $1"; }
warn(){ echo -e "  ${YELLOW}⚠️${NC}  $1"; }
err() { echo -e "  ${RED}❌${NC} $1"; }
info(){ echo -e "  ${BLUE}ℹ${NC}  $1"; }
hdr() { echo -e "\n${BOLD}${CYAN}═══ $1 ═══${NC}"; }
section() { echo -e "\n${BOLD}${BLUE}── $1 ──${NC}"; }

# ═══════════════════════════════════════════════════════════════════════════════
# DEFAULTS
# ═══════════════════════════════════════════════════════════════════════════════
INPUT_VCF=""
SAMPLE_ID="SAMPLE_001"
GWAS_CATALOG=""
OUT_DIR="$PIPELINE_DIR"
PLINK_BIN="plink"
PYTHON_BIN="python3"
THREADS=4; MEMORY=8000; NUM_PCS=20
LANG="both"

# Feature flags (all ON by default)
DRY_RUN=false; AUTO_MODE=false; INTERACTIVE=true
SKIP_GWAS=true; SKIP_PGS=true; SKIP_FULL_1000G=true
SKIP_PDF=false; SKIP_PHASE5=false; SKIP_PHASE6=false
VERBOSE=false; REPRODUCIBLE=false; CHECKPOINT=true
SEED=42

# Paths (auto-detected, user-overridable)
G1K_VCF=""
POP_PANEL=""
SNP_DB="$PIPELINE_DIR/data/snp_database_annotated.csv"
G1K_FULL_DIR="$PIPELINE_DIR/reference/1000G_full"
G1K_FULL_BFILE=""

# QC thresholds
GENO=0.1; MIND=0.1; MAF=0.01; HWE=0.000001
LD_WINDOW=50; LD_STEP=5; LD_R2=0.2
RISK_LOW=25; RISK_HIGH=75

# ═══════════════════════════════════════════════════════════════════════════════
# AUTO-DETECTION
# ═══════════════════════════════════════════════════════════════════════════════
auto_detect() {
    # --- PLINK ---
    for candidate in "plink2" "plink" "$HOME/bin/plink" "/usr/local/bin/plink" "/opt/homebrew/bin/plink"; do
        command -v "$candidate" &>/dev/null && { PLINK_BIN="$candidate"; break; }
    done

    # --- Python venv ---
    for VENV_DIR in "$PIPELINE_DIR/venv" "$PIPELINE_DIR/../nutrigenetic_prs_pipeline/venv"; do
        [[ -f "$VENV_DIR/bin/activate" ]] && { source "$VENV_DIR/bin/activate" 2>/dev/null || true; break; }
    done

    # --- 1000 Genomes ---
    if [[ -z "${G1K_VCF:-}" ]]; then
        for candidate in \
            "$G1K_FULL_DIR/1000G_full.vcf.gz" \
            "$G1K_FULL_DIR/1000G_full.bed" \
            "$PIPELINE_DIR/reference/1000G/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"; do
            if [[ -f "$candidate" ]]; then
                if [[ "$candidate" == *"1000G_full"* ]]; then
                    G1K_FULL_BFILE="${candidate%.vcf.gz}"; G1K_FULL_BFILE="${G1K_FULL_BFILE%.bed}"
                    SKIP_FULL_1000G=false
                else
                    G1K_VCF="$candidate"
                fi
                break
            fi
        done
    fi

    # --- Population panel ---
    if [[ -z "${POP_PANEL:-}" ]]; then
        for candidate in \
            "$G1K_FULL_DIR/population_panel.txt" \
            "$PIPELINE_DIR/reference/1000G/20130606_g1k_3202_samples_ped_population.txt"; do
            [[ -f "$candidate" ]] && { POP_PANEL="$candidate"; break; }
        done
    fi

    # --- Phase 6 module availability ---
    HAS_PHASE6_PCA=false; HAS_PHASE6_ANCESTRY=false; HAS_PHASE6_CALIBRATION=false
    [[ -f "$PIPELINE_DIR/scripts/pca_true_projection.py" ]] && HAS_PHASE6_PCA=true
    [[ -f "$PIPELINE_DIR/scripts/ancestry_inference_v2.py" ]] && HAS_PHASE6_ANCESTRY=true
    [[ -f "$PIPELINE_DIR/scripts/population_calibrate_v2.py" ]] && HAS_PHASE6_CALIBRATION=true

    # --- Phase 6 full 1000G available ---
    if [[ "$SKIP_FULL_1000G" == false ]] && [[ -n "${G1K_FULL_BFILE:-}" ]]; then
        HAS_FULL_1000G=true
    else
        HAS_FULL_1000G=false
    fi
}

# ═══════════════════════════════════════════════════════════════════════════════
# ARGUMENT PARSING
# ═══════════════════════════════════════════════════════════════════════════════
parse_args() {
    while [[ $# -gt 0 ]]; do
        case "$1" in
            --input-vcf)       INPUT_VCF="$2"; shift 2 ;;
            --sample-id)       SAMPLE_ID="$2"; shift 2 ;;
            --gwas)            GWAS_CATALOG="$2"; SKIP_GWAS=false; shift 2 ;;
            --output-dir)      OUT_DIR="$2"; shift 2 ;;
            --plink)           PLINK_BIN="$2"; shift 2 ;;
            --threads)         THREADS="$2"; shift 2 ;;
            --memory)          MEMORY="$2"; shift 2 ;;
            --num-pcs)         NUM_PCS="$2"; shift 2 ;;
            --1000g-vcf)       G1K_VCF="$2"; SKIP_FULL_1000G=true; shift 2 ;;
            --1000g-full-dir)  G1K_FULL_DIR="$2"; SKIP_FULL_1000G=false; shift 2 ;;
            --population-panel) POP_PANEL="$2"; shift 2 ;;
            --lang)            LANG="$2"; shift 2 ;;
            --geno)            GENO="$2"; shift 2 ;;
            --mind)            MIND="$2"; shift 2 ;;
            --maf)             MAF="$2"; shift 2 ;;
            --hwe)             HWE="$2"; shift 2 ;;
            --ld-r2)           LD_R2="$2"; shift 2 ;;
            --risk-low)        RISK_LOW="$2"; shift 2 ;;
            --risk-high)       RISK_HIGH="$2"; shift 2 ;;
            --auto)            AUTO_MODE=true; INTERACTIVE=false; shift ;;
            --dry-run)         DRY_RUN=true; shift ;;
            --reproducible)    REPRODUCIBLE=true; shift ;;
            --seed)            SEED="$2"; shift 2 ;;
            --skip-phase5)     SKIP_PHASE5=true; shift ;;
            --skip-phase6)     SKIP_PHASE6=true; shift ;;
            --no-checkpoint)   CHECKPOINT=false; shift ;;
            --with-pgs)        SKIP_PGS=false; shift ;;
            --skip-pdf)        SKIP_PDF=true; shift ;;
            --verbose|-v)      VERBOSE=true; shift ;;
            -h|--help)         show_help; exit 0 ;;
            *) echo -e "${RED}Unknown: $1${NC} (use --help)"; exit 1 ;;
        esac
    done
}

show_help() {
    cat << 'EOF'
╔══════════════════════════════════════════════════════════════════════════════╗
║           PRS RESEARCH PLATFORM v6.0.0 — Interactive Orchestrator            ║
╚══════════════════════════════════════════════════════════════════════════════╝

USAGE:
  ./run.sh                                  Interactive mode (recommended)
  ./run.sh --input-vcf <file> --auto        Full pipeline, non-interactive
  ./run.sh --input-vcf <file> --dry-run     Preview only
  ./run.sh --input-vcf <file> --skip-phase5 Skip validation modules

INTERACTIVE MODE (default):
  No arguments needed. The orchestrator will:
    1. Auto-detect available tools (PLINK, Python, 1000G reference)
    2. Ask for missing required inputs (VCF path)
    3. Present a menu of pipeline stages to run
    4. Show progress with color-coded output
    5. Skip completed stages (checkpointing)
    6. Integrate Phase 6 corrected modules when available

AUTO MODE (--auto):
  Runs the full pipeline without prompts. Use for batch/CI/CD.

KEY OPTIONS:
  --input-vcf FILE          DeepVariant VCF (GRCh37/hg19, bgzip compressed)
  --sample-id ID            Sample identifier (default: SAMPLE_001)
  --gwas FILE               GWAS Catalog summary statistics
  --lang en|es|both         Report language (default: both)
  --auto                    Non-interactive full pipeline
  --dry-run                 Preview stages without executing
  --reproducible            Deterministic mode (seed=42)
  --skip-phase5             Skip Phase 5 validation layer
  --skip-phase6             Skip Phase 6 corrected modules (use v4.0 fallbacks)
  --with-pgs                Enable PGS Catalog scoring
  --no-checkpoint           Re-run all stages (no skipping)

QC PARAMETERS:
  --geno FLOAT    SNP missingness (default: 0.1)
  --mind FLOAT    Individual missingness (default: 0.1)
  --maf FLOAT     MAF threshold (default: 0.01)
  --hwe FLOAT     HWE threshold (default: 1e-6)
  --ld-r2 FLOAT   LD r² threshold (default: 0.2)

PHASE 6 CORRECTIONS (auto-enabled when modules detected):
  ✓ Full 1000 Genomes genome-wide reference
  ✓ True PCA projection (reference-based, not merge)
  ✓ PCA ensemble ancestry inference
  ✓ Active PCA adjustment (not no-op)
  ✓ Empirical population calibration
  ✓ Corrected multi-method PRS (Σ(β×dosage))
  ✓ Orthogonal RQS v2

STAGES:
  A  VCF→PLINK           F  PRS Computation        K  Bilingual Interpretation
  B  Quality Control      G  PCA-Adjusted PRS        L  Bilingual Reports
  C  LD Pruning           H  Population Calibration  M  Admixture Analysis
  D  PCA + Projection     I  GWAS-LD Consistency     N  Multi-Method PRS
  E  Ancestry Inference   J  Uncertainty Propagation O  External Validation
                                                    P  Reproducibility
                                                    Q  Phase 5 Validation
                                                    R  Phase 6 Corrections

OUTPUTS:
  reports/report_en.html      English interactive report
  reports/report_es.html      Spanish interactive report
  reports/dashboard_en.html   Interactive validation dashboard
  prs/population_calibrated.csv  Final PRS with CIs
  pca/ancestry_inference.json   Ancestry classification
  validation/coverage_audit.csv Variant coverage per trait
  research/research_quality_v2.json  RQS v2
EOF
}

# ═══════════════════════════════════════════════════════════════════════════════
# CHECKPOINT / SKIP LOGIC
# ═══════════════════════════════════════════════════════════════════════════════
checkpoint_file() { echo "$OUT_DIR/.checkpoints/$1"; }
stage_done()  { [[ "$CHECKPOINT" == true ]] && [[ -f "$(checkpoint_file "$1")" ]]; }
mark_done()   { mkdir -p "$OUT_DIR/.checkpoints"; touch "$(checkpoint_file "$1")"; }

# ═══════════════════════════════════════════════════════════════════════════════
# INTERACTIVE HELPERS
# ═══════════════════════════════════════════════════════════════════════════════
ask_input_vcf() {
    if [[ -z "$INPUT_VCF" ]]; then
        echo -e "\n${BOLD}${CYAN}Enter path to VCF file:${NC}"
        echo -e "  (drag-and-drop file here, or type path)"
        read -r -p "  VCF path: " INPUT_VCF
    fi
    if [[ ! -f "$INPUT_VCF" ]]; then
        err "VCF not found: $INPUT_VCF"
        return 1
    fi
    ok "VCF: $INPUT_VCF"
}

ask_sample_id() {
    if [[ "$SAMPLE_ID" == "SAMPLE_001" ]]; then
        echo -e "\n${CYAN}Sample ID [SAMPLE_001]:${NC}"
        read -r -p "  > " custom_id
        [[ -n "$custom_id" ]] && SAMPLE_ID="$custom_id"
    fi
    ok "Sample: $SAMPLE_ID"
}

show_menu() {
    echo ""
    echo -e "${BOLD}${CYAN}╔══════════════════════════════════════════════════╗${NC}"
    echo -e "${BOLD}${CYAN}║   PRS Research Platform v6.0.0 — Stage Menu      ║${NC}"
    echo -e "${BOLD}${CYAN}╚══════════════════════════════════════════════════╝${NC}"
    echo ""
    echo -e "  ${BOLD}Core Pipeline:${NC}"
    echo -e "    ${GREEN}1${NC}) Full Pipeline (A→R, recommended)"
    echo -e "    ${GREEN}2${NC}) Genotype Processing (A→C: VCF→PLINK→QC→LD)"
    echo -e "    ${GREEN}3${NC}) Ancestry + PCA (D→E: PCA projection + inference)"
    echo -e "    ${GREEN}4${NC}) PRS Computation (F→H: Compute + Adjust + Calibrate)"
    echo -e "    ${GREEN}5${NC}) Reports (I→L: Consistency + Uncertainty + Interpretation + Reports)"
    echo ""
    echo -e "  ${BOLD}Advanced:${NC}"
    echo -e "    ${GREEN}6${NC}) Multi-Method PRS + Admixture (M→N)"
    echo -e "    ${GREEN}7${NC}) Validation Layer (O→Q: Phase 5 modules)"
    echo -e "    ${GREEN}8${NC}) Phase 6 Corrections (R: scientific fixes)"
    echo -e "    ${GREEN}9${NC}) Specific stage (choose from A→R)"
    echo ""
    echo -e "  ${BOLD}Utilities:${NC}"
    echo -e "    ${GREEN}d${NC}) Dry-run (preview only)"
    echo -e "    ${GREEN}c${NC}) Configure parameters"
    echo -e "    ${GREEN}h${NC}) Help"
    echo -e "    ${GREEN}q${NC}) Quit"
    echo ""
}

select_stage() {
    local choice
    read -r -p "  Select [1]: " choice
    choice="${choice:-1}"
    case "$choice" in
        1)  RUN_STAGES="A B C D E F G H I J K L M N O P Q R" ;;
        2)  RUN_STAGES="A B C" ;;
        3)  RUN_STAGES="D E" ;;
        4)  RUN_STAGES="F G H" ;;
        5)  RUN_STAGES="I J K L" ;;
        6)  RUN_STAGES="M N" ;;
        7)  RUN_STAGES="O P Q" ;;
        8)  RUN_STAGES="R" ;;
        9)  echo -e "  ${CYAN}Enter stage letters (e.g., 'A B F G'):${NC}"
            read -r -p "  > " RUN_STAGES ;;
        d)  DRY_RUN=true; RUN_STAGES="A B C D E F G H I J K L M N O P Q R" ;;
        c)  configure_params; select_stage ;;
        h)  show_help; select_stage ;;
        q)  echo "  Exiting."; exit 0 ;;
        *)  echo "  Invalid choice."; select_stage ;;
    esac
}

configure_params() {
    echo -e "\n  ${CYAN}── Configure Parameters ──${NC}"
    read -r -p "  Threads [$THREADS]: " val; [[ -n "$val" ]] && THREADS="$val"
    read -r -p "  Memory MB [$MEMORY]: " val; [[ -n "$val" ]] && MEMORY="$val"
    read -r -p "  Num PCs [$NUM_PCS]: " val; [[ -n "$val" ]] && NUM_PCS="$val"
    read -r -p "  Language (en/es/both) [$LANG]: " val; [[ -n "$val" ]] && LANG="$val"
    read -r -p "  LD r² threshold [$LD_R2]: " val; [[ -n "$val" ]] && LD_R2="$val"
    echo -e "  ${GREEN}Parameters updated.${NC}"
}

# ═══════════════════════════════════════════════════════════════════════════════
# STAGE HELPERS
# ═══════════════════════════════════════════════════════════════════════════════
run_stage() {
    local stage="$1" description="$2"
    if stage_done "$stage"; then
        info "Stage $stage: $description — already completed (skipping)"
        return 0
    fi
    hdr "STAGE $stage: $description"
    if [[ "$DRY_RUN" == true ]]; then
        info "[DRY-RUN] Would execute stage $stage"
        return 0
    fi
    return 1  # Caller should execute
}

run_python() {
    if [[ "$DRY_RUN" == true ]]; then
        info "[DRY-RUN] $PYTHON_BIN $*"
        return 0
    fi
    "$PYTHON_BIN" "$@"
}

run_plink_cmd() {
    if [[ "$DRY_RUN" == true ]]; then
        info "[DRY-RUN] $PLINK_BIN $*"
        return 0
    fi
    "$PLINK_BIN" "$@"
}

# ═══════════════════════════════════════════════════════════════════════════════
# PIPELINE STAGES
# ═══════════════════════════════════════════════════════════════════════════════

# ── STAGE A: VCF → PLINK ──────────────────────────────────────────────────────
stage_A() {
    run_stage "A" "VCF → PLINK Conversion" || return 0
    bash "$SCRIPT_DIR/scripts/01_vcf_to_plink.sh" \
        --input-vcf "$INPUT_VCF" --out-dir plink/ \
        --plink "$PLINK_BIN" --threads "$THREADS" --memory "$MEMORY"
    mark_done "A"
}

# ── STAGE B: Quality Control ──────────────────────────────────────────────────
stage_B() {
    run_stage "B" "Quality Control" || return 0
    bash "$SCRIPT_DIR/scripts/02_quality_control.sh" \
        --bfile plink/cohort --out-dir qc/ \
        --plink "$PLINK_BIN" --threads "$THREADS" --memory "$MEMORY" \
        --geno "$GENO" --mind "$MIND" --maf "$MAF" --hwe "$HWE"
    mark_done "B"
}

# ── STAGE C: LD Pruning ───────────────────────────────────────────────────────
stage_C() {
    run_stage "C" "Ancestry-Matched LD Pruning" || return 0

    if [[ -f "$G1K_VCF" ]] && [[ -f "$POP_PANEL" ]]; then
        bash "$SCRIPT_DIR/scripts/03_ld_ancestry_prune.sh" \
            --bfile qc/qc_filtered --out-dir plink/ \
            --1000g-vcf "$G1K_VCF" --population-panel "$POP_PANEL" \
            --plink "$PLINK_BIN" --threads "$THREADS" --memory "$MEMORY" \
            --r2 "$LD_R2" --union-mode conservative 2>&1 || {
            warn "Ancestry LD pruning fell back to self-pruning"
            bash "$SCRIPT_DIR/scripts/03_ld_pruning.sh" \
                --bfile qc/qc_filtered --out-dir plink/ \
                --plink "$PLINK_BIN" --threads "$THREADS" --memory "$MEMORY" \
                --window-size "$LD_WINDOW" --step-size "$LD_STEP" --r2 "$LD_R2"
        }
    else
        bash "$SCRIPT_DIR/scripts/03_ld_pruning.sh" \
            --bfile qc/qc_filtered --out-dir plink/ \
            --plink "$PLINK_BIN" --threads "$THREADS" --memory "$MEMORY" \
            --window-size "$LD_WINDOW" --step-size "$LD_STEP" --r2 "$LD_R2"
    fi
    mark_done "C"
}

# ── STAGE D: PCA Projection ──────────────────────────────────────────────────
stage_D() {
    run_stage "D" "PCA Projection" || return 0

    # Phase 6: Use true PCA projection if full 1000G is available
    if [[ "$SKIP_PHASE6" == false ]] && [[ "$HAS_FULL_1000G" == true ]] && [[ -f "$PIPELINE_DIR/scripts/pca_true_projection.py" ]]; then
        info "Phase 6: Using true PCA projection (Price et al. 2006 method)"

        # Check if reference PCA model already trained
        if [[ -f pca/reference_pca_model.pkl ]]; then
            info "Reference PCA model found — projecting sample only"
            run_python "$SCRIPT_DIR/scripts/pca_true_projection.py" \
                --ref-bfile "${G1K_FULL_BFILE:-reference/1000G_full/1000G_full}" \
                --target-bfile plink/ld_pruned_dataset \
                --sample-id "$SAMPLE_ID" \
                --n-components "$NUM_PCS" \
                --output-dir pca/ --plink "$PLINK_BIN" \
                --threads "$THREADS" --memory "$MEMORY" \
                --project-only
        else
            info "Training reference PCA model (this may take a while)..."
            run_python "$SCRIPT_DIR/scripts/pca_true_projection.py" \
                --ref-bfile "${G1K_FULL_BFILE:-reference/1000G_full/1000G_full}" \
                --target-bfile plink/ld_pruned_dataset \
                --sample-id "$SAMPLE_ID" \
                --n-components "$NUM_PCS" \
                --output-dir pca/ --plink "$PLINK_BIN" \
                --threads "$THREADS" --memory "$MEMORY"
        fi
    elif [[ -f "$G1K_VCF" ]]; then
        # Legacy: 1000G chr22 PCA
        warn "Phase 6 full 1000G not available — using chr22 PCA (legacy)"
        bash "$SCRIPT_DIR/scripts/04_pca_1000G.sh" \
            --1000g-vcf "$G1K_VCF" --target-bfile plink/ld_pruned_dataset \
            --population-panel "$POP_PANEL" --out-dir pca/ \
            --plink "$PLINK_BIN" --threads "$THREADS" --memory "$MEMORY" \
            --num-pcs "$NUM_PCS" 2>&1 || {
            warn "1000G PCA projection fell back — using internal PCA"
            bash "$SCRIPT_DIR/scripts/04_pca_analysis.sh" \
                --bfile plink/ld_pruned_dataset --out-dir pca/ \
                --plink "$PLINK_BIN" --num-pcs "$NUM_PCS"
        }
    else
        bash "$SCRIPT_DIR/scripts/04_pca_analysis.sh" \
            --bfile plink/ld_pruned_dataset --out-dir pca/ \
            --plink "$PLINK_BIN" --num-pcs "$NUM_PCS"
    fi
    mark_done "D"
}

# ── STAGE E: Ancestry Inference ───────────────────────────────────────────────
stage_E() {
    run_stage "E" "Ancestry Inference" || return 0

    # Phase 6: PCA ensemble ancestry if available
    if [[ "$SKIP_PHASE6" == false ]] && [[ -f "$PIPELINE_DIR/scripts/ancestry_inference_v2.py" ]]; then
        if [[ -f pca/projected_sample.csv ]] && [[ -f pca/reference_centroids.csv ]]; then
            info "Phase 6: Using PCA ensemble ancestry inference (genome-wide)"
            run_python "$SCRIPT_DIR/scripts/ancestry_inference_v2.py" \
                --projected-pcs pca/projected_sample.csv \
                --ref-pcs pca/1000G_pca.eigenvec \
                --ref-centroids pca/reference_centroids.csv \
                --population-panel "$POP_PANEL" \
                --output-dir ancestry/ --sample-id "$SAMPLE_ID" \
                --n-pcs "$NUM_PCS"
            mark_done "E"
            return 0
        fi
    fi

    # Legacy: 1000G projection-based ancestry
    if [[ -f pca/1000G_pca.eigenvec ]] && [[ -f pca/target_pcs.eigenvec ]]; then
        run_python "$SCRIPT_DIR/scripts/ancestry_inference.py" \
            --reference-pcs pca/1000G_pca.eigenvec \
            --target-pcs pca/target_pcs.eigenvec \
            --population-panel "$POP_PANEL" \
            --output-dir pca/ --num-pcs "$NUM_PCS" 2>&1 || true
    fi
    mark_done "E"
}

# ── STAGE F: PRS Computation ──────────────────────────────────────────────────
stage_F() {
    run_stage "F" "PRS Computation" || return 0

    if [[ -f "$SNP_DB" ]]; then
        run_python -c "
import sys, importlib.util, logging, os
logging.basicConfig(level=logging.INFO, format='%(asctime)s [%(levelname)s] %(message)s')
spec = importlib.util.spec_from_file_location('prs_compute', '$SCRIPT_DIR/scripts/06_prs_compute.py')
m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
vcf = '$INPUT_VCF'
for c in ['$INPUT_VCF', 'input.vcf.gz']:
    if c and os.path.exists(c): vcf = c; break
df = m.compute_prs_from_curated_database(
    bfile='plink/ld_pruned_dataset', snp_database_path='$SNP_DB',
    output_dir='prs/', sample_id='$SAMPLE_ID', vcf_path=vcf)
print(f'  PRS: {len(df[\"trait\"].unique())} traits, {df[\"n_snps_used\"].sum()}/{df[\"n_snps\"].sum()} SNPs')
"
    fi
    mark_done "F"
}

# ── STAGE G: PCA-Adjusted PRS ─────────────────────────────────────────────────
stage_G() {
    run_stage "G" "PCA-Adjusted PRS" || return 0

    local prs_in="prs/prs_raw.csv"
    [[ ! -f "$prs_in" ]] && { warn "No PRS data — skipping"; return 0; }

    # Phase 6: Active PCA adjustment
    if [[ "$SKIP_PHASE6" == false ]] && [[ -f pca/projected_sample.csv ]]; then
        if [[ -f "$PIPELINE_DIR/scripts/pca_adjust_v2.py" ]]; then
            info "Phase 6: Active PCA regression adjustment"
            run_python "$SCRIPT_DIR/scripts/pca_adjust_v2.py" \
                --prs-data "$prs_in" \
                --sample-pcs pca/projected_sample.csv \
                --output-dir prs/ --sample-id "$SAMPLE_ID" \
                --n-pcs "$NUM_PCS"
        else
            # Use original module but with actual regression (not no-op)
            info "Using v4.0 PCA adjustment (full regression)"
            run_python "$SCRIPT_DIR/scripts/07_prs_adjust.py" \
                --prs-raw "$prs_in" \
                --eigenvec pca/pca_results.eigenvec \
                --output-dir prs/ --num-pcs "$NUM_PCS"
        fi
    else
        # Old no-op path (legacy)
        warn "Using legacy PCA adjustment (may be no-op)"
        run_python -c "
import pandas as pd
df = pd.read_csv('$prs_in')
df['prs_adjusted'] = df['prs_raw']
df.to_csv('prs/prs_adjusted.csv', index=False)
print('  prs/prs_adjusted.csv')
"
    fi
    mark_done "G"
}

# ── STAGE H: Population Calibration ───────────────────────────────────────────
stage_H() {
    run_stage "H" "Population Calibration" || return 0

    local prs_in="prs/prs_adjusted.csv"
    [[ ! -f "$prs_in" ]] && prs_in="prs/prs_raw.csv"
    [[ ! -f "$prs_in" ]] && { warn "No PRS data — skipping"; return 0; }

    # Phase 6: Empirical calibration with real distributions
    if [[ "$SKIP_PHASE6" == false ]] && [[ -f "$PIPELINE_DIR/scripts/population_calibrate_v2.py" ]]; then
        local ancestry_json="ancestry/posterior_probabilities.json"
        [[ ! -f "$ancestry_json" ]] && ancestry_json="pca/ancestry_inference.json"

        if [[ -f "$ancestry_json" ]]; then
            info "Phase 6: Empirical population calibration"

            # Check for pre-computed reference distributions
            local ref_dist_dir="references/population_distributions"
            if [[ -f "$ref_dist_dir/reference_distributions.json" ]]; then
                run_python "$SCRIPT_DIR/scripts/population_calibrate_v2.py" \
                    --sample-prs "$prs_in" --ancestry-json "$ancestry_json" \
                    --ref-dist-dir "$ref_dist_dir" \
                    --output-dir prs/ --sample-id "$SAMPLE_ID" \
                    --calibrate-only
            else
                warn "No reference distributions computed yet — using v4.0 fallback"
                run_python "$SCRIPT_DIR/scripts/11_population_calibrate.py" \
                    --prs-adjusted "$prs_in" --ancestry "$ancestry_json" \
                    --output-dir prs/
            fi
        else
            run_python "$SCRIPT_DIR/scripts/08_normalize.py" \
                --prs-adjusted "$prs_in" --output-dir prs/ \
                --method percentile --low-threshold "$RISK_LOW" --high-threshold "$RISK_HIGH"
        fi
    else
        if [[ -f pca/ancestry_inference.json ]]; then
            run_python "$SCRIPT_DIR/scripts/11_population_calibrate.py" \
                --prs-adjusted "$prs_in" --ancestry pca/ancestry_inference.json \
                --output-dir prs/
        else
            run_python "$SCRIPT_DIR/scripts/08_normalize.py" \
                --prs-adjusted "$prs_in" --output-dir prs/ \
                --method percentile --low-threshold "$RISK_LOW" --high-threshold "$RISK_HIGH"
        fi
    fi
    mark_done "H"
}

# ── STAGE I: GWAS-LD Consistency ──────────────────────────────────────────────
stage_I() {
    run_stage "I" "GWAS-LD-Ancestry Consistency" || return 0
    if [[ -f pca/ancestry_inference.json ]] || [[ -f ancestry/classification_report.json ]]; then
        local anc_json="ancestry/classification_report.json"
        [[ ! -f "$anc_json" ]] && anc_json="pca/ancestry_inference.json"
        run_python "$SCRIPT_DIR/scripts/13_gwas_ld_consistency_check.py" \
            --ancestry "$anc_json" --curated-db "$SNP_DB" \
            --output-dir prs/ --no-strict 2>&1 || true
    fi
    mark_done "I"
}

# ── STAGE J: Uncertainty Propagation ──────────────────────────────────────────
stage_J() {
    run_stage "J" "Uncertainty Propagation" || return 0

    local prs_cal="prs/population_calibrated_v2.csv"
    [[ ! -f "$prs_cal" ]] && prs_cal="prs/population_calibrated.csv"
    [[ ! -f "$prs_cal" ]] && prs_cal="prs/prs_normalized.csv"
    [[ ! -f "$prs_cal" ]] && { warn "No calibrated PRS — skipping"; return 0; }

    local anc_json="ancestry/posterior_probabilities.json"
    [[ ! -f "$anc_json" ]] && anc_json="pca/ancestry_inference.json"

    if [[ -f "$prs_cal" ]] && [[ -f "$anc_json" ]]; then
        run_python "$SCRIPT_DIR/scripts/14_uncertainty_propagation.py" \
            --prs-data "$prs_cal" --ancestry "$anc_json" --snp-db "$SNP_DB" \
            --vcf "$INPUT_VCF" --consistency-check prs/consistency_check_report.json \
            --output-dir prs/
    fi
    mark_done "J"
}

# ── STAGE K: Bilingual Interpretation ─────────────────────────────────────────
stage_K() {
    run_stage "K" "Bilingual Interpretation" || return 0

    local prs_cal="prs/population_calibrated_v2.csv"
    [[ ! -f "$prs_cal" ]] && prs_cal="prs/population_calibrated.csv"
    [[ ! -f "$prs_cal" ]] && { warn "No calibrated PRS — skipping"; return 0; }

    local anc_json="ancestry/classification_report.json"
    [[ ! -f "$anc_json" ]] && anc_json="pca/ancestry_inference.json"

    if [[ -f "$prs_cal" ]] && [[ -f "$anc_json" ]]; then
        run_python "$SCRIPT_DIR/scripts/bilingual_interpretation.py" \
            --prs-calibrated "$prs_cal" --ancestry "$anc_json" \
            --output-dir interpretations/ --lang "$LANG"
    fi
    mark_done "K"
}

# ── STAGE L: Bilingual Reports ────────────────────────────────────────────────
stage_L() {
    run_stage "L" "Bilingual Reports" || return 0

    local prs_cal="prs/population_calibrated_v2.csv"
    [[ ! -f "$prs_cal" ]] && prs_cal="prs/population_calibrated.csv"

    if [[ -f interpretations/interpretations.json ]] && [[ -f "$prs_cal" ]]; then
        run_python -c "
import json, pandas as pd; from datetime import datetime; from html import escape as h; from pathlib import Path
Path('reports').mkdir(parents=True, exist_ok=True)
prs = pd.read_csv('$prs_cal')
with open('interpretations/interpretations.json') as fh: interps = json.load(fh)
anc_paths = ['ancestry/classification_report.json','pca/ancestry_inference.json']
ancestry = {}
for ap in anc_paths:
    try:
        with open(ap) as fh: ancestry = json.load(fh); break
    except: pass
try: uncertainty = pd.read_csv('prs/prs_uncertainty.csv')
except: uncertainty = None
try:
    with open('prs/consistency_check_report.json') as fh: consistency = json.load(fh)
except: consistency = None

ICONS = {'Caffeine metabolism':'☕','Lipid metabolism':'❤️','Omega-3 metabolism':'🐟',
         'Folate & methylation':'🧬','Lactose intolerance':'🥛','Obesity predisposition':'⚖️',
         'Dopamine regulation':'🧠','Detoxification':'🛡️','Vitamin D metabolism':'☀️','Glucose metabolism':'🍬'}
COLORS = {'low':{'bar':'#27ae60','bg':'#d5f5e3','text':'#1e8449'},
          'medium':{'bar':'#f39c12','bg':'#fdebd0','text':'#b7950b'},
          'high':{'bar':'#e74c3c','bg':'#fadbd8','text':'#c0392b'}}
POP_NAMES = {'en':{'EUR':'European','AFR':'African','EAS':'East Asian','SAS':'South Asian','AMR':'Admixed American'},
             'es':{'EUR':'Europea','AFR':'Africana','EAS':'Asiática Oriental','SAS':'Sudasiática','AMR':'Americana Mixta'}}
RISK_LABELS = {'en':{'low':'Lower Risk','medium':'Average Risk','high':'Higher Risk'},
               'es':{'low':'Riesgo Bajo','medium':'Riesgo Promedio','high':'Riesgo Elevado'}}

for lang in ['en','es']:
    is_es = lang == 'es'
    interp = interps.get(lang, interps.get('en', {}))
    meta = interp.get('metadata', {})
    pop = meta.get('assigned_population','EUR')
    pop_name = POP_NAMES[lang].get(pop, pop)
    hi, md, lo = meta.get('high_risk_count',0), meta.get('medium_risk_count',0), meta.get('low_risk_count',0)
    date_str = datetime.now().strftime('%Y-%m-%d %H:%M UTC')
    T = lambda en,es: es if is_es else en

    prs_rows = ''
    for _, row in prs.iterrows():
        trait = row['trait']; risk = row.get('risk_category','medium')
        z = float(row.get('z_score_population',0)); pctl = float(row.get('percentile_population',50))
        icon = ICONS.get(trait,'🧬'); c = COLORS.get(risk,COLORS['medium'])
        bar = int(min(max((z+3)/6*100,0),100))
        prs_rows += f'<tr><td>{icon} {h(trait)}</td><td><strong>{z:+.2f}</strong></td><td>{pctl:.1f}%</td><td><span class=\"risk-category-badge {risk}\">{RISK_LABELS[lang].get(risk,risk).upper()}</span></td><td><div style=\"display:flex;align-items:center\"><div style=\"flex:1;height:10px;background:#e9ecef;border-radius:5px;overflow:hidden\"><div style=\"width:{bar}%;height:100%;background:{c['bar']};border-radius:5px\"></div></div></div></td></tr>'

    trait_secs = ''
    ti_map = interp.get('trait_interpretations', {})
    for _, row in prs.iterrows():
        trait = row['trait']; risk = row.get('risk_category','medium')
        z = float(row.get('z_score_population',0)); pctl = float(row.get('percentile_population',50))
        icon = ICONS.get(trait,'🧬'); c = COLORS.get(risk,COLORS['medium'])
        bar = int(min(max((z+3)/6*100,0),100))
        ti = ti_map.get(trait, {})
        dietary = ''.join(f'<li>{h(item)}</li>' for item in ti.get('dietary_context',[]))
        dietary_html = f'<div style=\"background:#eaf2f8;border-radius:8px;padding:1rem;margin-top:1rem\"><h4>{T(\"Dietary & Lifestyle Context\",\"Contexto Dietético\",is_es)}</h4><ul>{dietary}</ul></div>' if dietary else ''
        trait_secs += f'<div style=\"background:#fff;border-radius:8px;box-shadow:0 1px 3px rgba(0,0,0,.08);margin-bottom:1.5rem;overflow:hidden\"><div style=\"padding:1rem 1.25rem;display:flex;justify-content:space-between;align-items:center;background:{c[\"bg\"]}\"><span style=\"font-size:1.1rem;font-weight:600\">{icon} {h(trait)}</span><span><span class=\"risk-category-badge {risk}\">{RISK_LABELS[lang].get(risk,risk).upper()}</span></span></div><div style=\"padding:1.25rem;border-top:1px solid #dee2e6\"><div style=\"display:grid;grid-template-columns:repeat(auto-fit,minmax(110px,1fr));gap:1rem;margin-bottom:1rem\"><div style=\"text-align:center\"><div style=\"font-size:1.3rem;font-weight:700;color:{c['text']}\">{z:+.2f}</div><div style=\"font-size:.75rem;color:#7f8c8d\">Z-Score</div></div><div style=\"text-align:center\"><div style=\"font-size:1.3rem;font-weight:700\">{pctl:.1f}%</div><div style=\"font-size:.75rem;color:#7f8c8d\">{T(\"Percentile\",\"Percentil\",is_es)}</div></div></div><div style=\"display:flex;align-items:center;gap:.5rem;margin:1rem 0\"><span style=\"font-size:.7rem;color:#27ae60\">{T(\"Low\",\"Bajo\",is_es)}</span><div style=\"flex:1;height:10px;background:#e9ecef;border-radius:5px;overflow:hidden\"><div style=\"width:{bar}%;height:100%;background:{c['bar']};border-radius:5px\"></div></div><span style=\"font-size:.7rem;color:#e74c3c\">{T(\"High\",\"Alto\",is_es)}</span></div><div style=\"background:#f8f9fa;border-radius:8px;padding:1rem;margin:1rem 0;font-size:.9rem;line-height:1.7\"><strong>🧬 {h(ti.get(\"gene\",\"\"))}</strong><p style=\"color:#7f8c8d;font-size:.8rem;margin-top:.25rem\">{h(ti.get(\"pathway\",\"\"))}</p><p style=\"margin-top:.75rem;white-space:pre-line\">{h(ti.get(\"description\",\"\"))}</p><p style=\"margin-top:.75rem;white-space:pre-line;font-weight:500\">{h(ti.get(\"narrative\",\"\"))}</p></div>{dietary_html}</div></div>'

    unc_html = ''
    if uncertainty is not None and len(uncertainty) > 0:
        u_rows = ''.join(f'<tr><td>{h(r[\"trait\"])}</td><td>{float(r[\"prs\"]):.3f} ± {float(r[\"prs_se\"]):.3f}</td><td>[{float(r[\"ci_95_lower\"]):.3f}, {float(r[\"ci_95_upper\"]):.3f}]</td><td>{float(r[\"uncertainty_score\"]):.3f}</td></tr>' for _, r in uncertainty.iterrows())
        unc_html = f'<section><h2 style=\"font-size:1.3rem;font-weight:600;margin:2rem 0 1rem;padding-bottom:.5rem;border-bottom:2px solid #dee2e6\">📊 {T(\"Uncertainty Overview\",\"Resumen de Incertidumbre\",is_es)}</h2><table style=\"width:100%;border-collapse:collapse;background:#fff;border-radius:8px;overflow:hidden;box-shadow:0 1px 3px rgba(0,0,0,.08);margin:1rem 0\"><thead><tr style=\"background:#f1f3f5\"><th style=\"padding:.65rem .75rem;text-align:left;font-weight:600;font-size:.8rem;text-transform:uppercase\">{T(\"Trait\",\"Rasgo\",is_es)}</th><th style=\"padding:.65rem .75rem;text-align:left;font-weight:600;font-size:.8rem;text-transform:uppercase\">PRS ± SE</th><th style=\"padding:.65rem .75rem;text-align:left;font-weight:600;font-size:.8rem;text-transform:uppercase\">95% CI</th><th style=\"padding:.65rem .75rem;text-align:left;font-weight:600;font-size:.8rem;text-transform:uppercase\">{T(\"Uncert.\",\"Incert.\",is_es)}</th></tr></thead><tbody>{u_rows}</tbody></table></section>'

    cons_html = ''
    if consistency:
        passed = consistency.get('passed',True)
        sc = '#27ae60' if passed else '#e74c3c'
        st = '✅ '+T('PASSED','COMPATIBLE',is_es) if passed else '❌ '+T('FAILED','INCOMPATIBLE',is_es)
        cons_html = f'<section><h2 style=\"font-size:1.3rem;font-weight:600;margin:2rem 0 1rem;padding-bottom:.5rem;border-bottom:2px solid #dee2e6\">🔬 {T(\"GWAS-Ancestry Compatibility\",\"Compatibilidad GWAS-Ancestral\",is_es)}</h2><div style=\"display:grid;grid-template-columns:repeat(auto-fit,minmax(160px,1fr));gap:1rem;margin:1.5rem 0\"><div style=\"background:#fff;border-radius:8px;padding:1.25rem;text-align:center;box-shadow:0 1px 3px rgba(0,0,0,.08)\"><div style=\"font-size:2rem;font-weight:800;color:{sc}\">{st}</div><div style=\"font-size:.85rem;color:#7f8c8d\">{T(\"Validation Status\",\"Estado de Validación\",is_es)}</div></div></div></section>'

    css = '''*{box-sizing:border-box;margin:0;padding:0}body{font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',Roboto,sans-serif;color:#2c3e50;background:#f8f9fa;line-height:1.6}.report-header{background:linear-gradient(135deg,#1a1a2e,#16213e,#0f3460);color:#fff;padding:3rem 2rem;text-align:center}.container{max-width:960px;margin:0 auto;padding:2rem 1.5rem}.risk-category-badge{display:inline-block;padding:.15rem .6rem;border-radius:4px;font-size:.75rem;font-weight:700;text-transform:uppercase}.risk-category-badge.low{background:#d5f5e3;color:#1e8449}.risk-category-badge.medium{background:#fdebd0;color:#b7950b}.risk-category-badge.high{background:#fadbd8;color:#c0392b}table{width:100%;border-collapse:collapse}td{padding:.65rem .75rem;border-top:1px solid #dee2e6;font-size:.9rem}tr:hover{background:#f8f9fa}.disclaimer-box{background:#fef9e7;border:1px solid #f9e79f;border-radius:8px;padding:1.25rem;margin-top:2rem;font-size:.8rem;color:#7d6608}.report-footer{text-align:center;padding:2rem;font-size:.75rem;color:#7f8c8d;border-top:1px solid #dee2e6;margin-top:2rem}@media print{body{background:#fff}}'''

    html = f'''<!DOCTYPE html><html lang=\"{lang}\">
<head><meta charset=\"UTF-8\"><title>{T(\"PRS Research Report\",\"Informe PRS\",is_es)} — {h(\"$SAMPLE_ID\")}</title><style>{css}</style></head>
<body>
<header class=\"report-header\"><h1 style=\"font-size:2.2rem\">{T(\"PRS Research Report\",\"Informe de Investigación PRS\",is_es)}</h1><div style=\"font-size:1.1rem;opacity:.9\">{T(\"Population-Calibrated PRS Analysis\",\"Análisis de PRS Calibrado por Población\",is_es)}</div><div style=\"margin-top:1rem;font-size:.85rem;opacity:.7\">{T(\"Sample\",\"Muestra\",is_es)}: {h(\"$SAMPLE_ID\")} | {T(\"Date\",\"Fecha\",is_es)}: {date_str} | Pipeline v6.0.0 | GRCh37/hg19</div></header>
<div class=\"container\">
<section><h2 style=\"font-size:1.3rem;font-weight:600;margin:2rem 0 1rem;padding-bottom:.5rem;border-bottom:2px solid #dee2e6\">📊 {T(\"Genetic Overview\",\"Resumen Genético\",is_es)}</h2><div style=\"display:grid;grid-template-columns:repeat(auto-fit,minmax(160px,1fr));gap:1rem;margin-bottom:2rem\"><div style=\"background:#fff;border-radius:8px;padding:1.25rem;text-align:center;box-shadow:0 1px 3px rgba(0,0,0,.08);border-left:4px solid #e74c3c\"><div style=\"font-size:2.2rem;font-weight:800;color:#e74c3c\">{hi}</div><div style=\"font-size:.85rem;color:#7f8c8d\">{T(\"Higher Risk\",\"Riesgo Elevado\",is_es)}</div></div><div style=\"background:#fff;border-radius:8px;padding:1.25rem;text-align:center;box-shadow:0 1px 3px rgba(0,0,0,.08);border-left:4px solid #f39c12\"><div style=\"font-size:2.2rem;font-weight:800;color:#f39c12\">{md}</div><div style=\"font-size:.85rem;color:#7f8c8d\">{T(\"Average Risk\",\"Riesgo Promedio\",is_es)}</div></div><div style=\"background:#fff;border-radius:8px;padding:1.25rem;text-align:center;box-shadow:0 1px 3px rgba(0,0,0,.08);border-left:4px solid #27ae60\"><div style=\"font-size:2.2rem;font-weight:800;color:#27ae60\">{lo}</div><div style=\"font-size:.85rem;color:#7f8c8d\">{T(\"Lower Risk\",\"Riesgo Bajo\",is_es)}</div></div></div><p style=\"white-space:pre-line;background:#fff;padding:1.25rem;border-radius:8px;box-shadow:0 1px 3px rgba(0,0,0,.08)\">{h(interp.get(\"global_summary\",\"\"))}</p></section>
<section><h2 style=\"font-size:1.3rem;font-weight:600;margin:2rem 0 1rem;padding-bottom:.5rem;border-bottom:2px solid #dee2e6\">🌍 {T(\"Ancestry\",\"Ascendencia\",is_es)}</h2><p style=\"font-size:1.1rem\"><strong>{pop_name} ({pop})</strong></p></section>
<section><h2 style=\"font-size:1.3rem;font-weight:600;margin:2rem 0 1rem;padding-bottom:.5rem;border-bottom:2px solid #dee2e6\">📈 {T(\"Population-Calibrated PRS\",\"PRS Calibrado por Población\",is_es)}</h2><table><thead><tr style=\"background:#f1f3f5\"><th style=\"padding:.65rem .75rem;text-align:left;font-weight:600;font-size:.8rem;text-transform:uppercase\">{T(\"Trait\",\"Rasgo\",is_es)}</th><th style=\"padding:.65rem .75rem;text-align:left;font-weight:600;font-size:.8rem;text-transform:uppercase\">Z-Score</th><th style=\"padding:.65rem .75rem;text-align:left;font-weight:600;font-size:.8rem;text-transform:uppercase\">{T(\"Percentile\",\"Percentil\",is_es)}</th><th style=\"padding:.65rem .75rem;text-align:left;font-weight:600;font-size:.8rem;text-transform:uppercase\">{T(\"Risk\",\"Riesgo\",is_es)}</th></tr></thead><tbody>{prs_rows}</tbody></table></section>
<h2 style=\"font-size:1.3rem;font-weight:600;margin:2rem 0 1rem;padding-bottom:.5rem;border-bottom:2px solid #dee2e6\">🔬 {T(\"Nutrigenetic Trait Analysis\",\"Análisis Nutrigenético\",is_es)}</h2>{trait_secs}
{unc_html}{cons_html}
<section><h2 style=\"font-size:1.3rem;font-weight:600;margin:2rem 0 1rem;padding-bottom:.5rem;border-bottom:2px solid #dee2e6\">📐 {T(\"Methodology\",\"Metodología\",is_es)}</h2><p><strong>{T(\"PRS Model:\",\"Modelo PRS:\",is_es)}</strong> PRS = Σ (βⱼ × Gᵢⱼ) {T(\"with\",\"con\",is_es)} LD pruning (r² < {LD_R2}), true PCA projection onto 1000 Genomes, empirical population calibration, 3-layer uncertainty propagation, Phase 6 corrections.</p></section>
<div class=\"disclaimer-box\"><h3 style=\"color:#9a7d0a;margin-bottom:.5rem\">⚠️ {T(\"Important Disclaimer\",\"Aviso Importante\",is_es)}</h3><p style=\"white-space:pre-line\">{h(interp.get(\"disclaimer\",\"\"))}</p></div>
</div>
<footer class=\"report-footer\"><p>{T(\"PRS Research Platform\",\"Plataforma PRS\",is_es)} v6.0.0 — {T(\"Generated\",\"Generado\",is_es)} {date_str}</p><p>{T(\"Sample\",\"Muestra\",is_es)}: {h(\"$SAMPLE_ID\")} | GRCh37/hg19</p></footer>
</body></html>'''

    with open(f'reports/report_{lang}.html', 'w', encoding='utf-8') as fh:
        fh.write(html)

    try:
        from weasyprint import HTML
        HTML(string=html).write_pdf(f'reports/report_{lang}.pdf')
        print(f'✅ [{lang}] reports/report_{lang}.html + reports/report_{lang}.pdf')
    except Exception as e:
        print(f'✅ [{lang}] reports/report_{lang}.html  ⚠️ PDF: {e}')

print('✅ Bilingual reports generated (v6.0.0)')
"
    fi
    mark_done "L"
}

# ── STAGE M: Admixture Analysis ────────────────────────────────────────────────
stage_M() {
    run_stage "M" "Admixture Analysis" || return 0

    # Phase 6: Corrected admixture engine
    if [[ "$SKIP_PHASE6" == false ]] && [[ -f "$PIPELINE_DIR/scripts/admixture_engine_v2.py" ]]; then
        if [[ -f pca/projected_sample.csv ]] && [[ -f pca/reference_centroids.csv ]]; then
            info "Phase 6: PCA softmax admixture (corrected)"
            run_python "$SCRIPT_DIR/scripts/admixture_engine_v2.py" \
                --sample-pcs pca/projected_sample.csv \
                --ref-centroids pca/reference_centroids.csv \
                --output-dir ancestry/ --sample-id "$SAMPLE_ID"
            mark_done "M"
            return 0
        fi
    fi

    # Legacy
    if [[ -f pca/ancestry_inference.json ]]; then
        run_python "$SCRIPT_DIR/scripts/admixture_engine.py" \
            --ancestry pca/ancestry_inference.json --output-dir pca/
    fi
    mark_done "M"
}

# ── STAGE N: Multi-Method PRS ─────────────────────────────────────────────────
stage_N() {
    run_stage "N" "Multi-Method PRS" || return 0

    # Phase 6: Corrected formula
    if [[ "$SKIP_PHASE6" == false ]] && [[ -f "$PIPELINE_DIR/scripts/prs_multi_method_v2.py" ]]; then
        # Check for GWAS score files
        local score_file=""
        for candidate in gwas/scores/*.score; do
            [[ -f "$candidate" ]] && { score_file="$candidate"; break; }
        done
        if [[ -n "$score_file" ]]; then
            info "Phase 6: Corrected multi-method PRS (Σ(β×dosage))"
            run_python "$SCRIPT_DIR/scripts/prs_multi_method_v2.py" \
                --bfile plink/ld_pruned_dataset --score-file "$score_file" \
                --output-dir prs/ --sample-id "$SAMPLE_ID" --plink "$PLINK_BIN"
        else
            warn "No GWAS score files found — using v4.0 multi-method fallback"
            [[ -f "$SNP_DB" ]] && run_python "$SCRIPT_DIR/scripts/prs_multi_method.py" \
                --bfile plink/ld_pruned_dataset --snp-db "$SNP_DB" \
                --output-dir prs/ --sample-id "$SAMPLE_ID" --plink "$PLINK_BIN" 2>&1 || true
        fi
    else
        [[ -f "$SNP_DB" ]] && run_python "$SCRIPT_DIR/scripts/prs_multi_method.py" \
            --bfile plink/ld_pruned_dataset --snp-db "$SNP_DB" \
            --output-dir prs/ --sample-id "$SAMPLE_ID" --plink "$PLINK_BIN" 2>&1 || true
    fi
    mark_done "N"
}

# ── STAGE O: External Validation ──────────────────────────────────────────────
stage_O() {
    run_stage "O" "External Validation" || return 0
    local prs_cal="prs/population_calibrated_v2.csv"
    [[ ! -f "$prs_cal" ]] && prs_cal="prs/population_calibrated.csv"
    [[ -f "$prs_cal" ]] && run_python "$SCRIPT_DIR/scripts/validation_framework.py" \
        --prs-data "$prs_cal" --snp-db "$SNP_DB" --output-dir validation/
    mark_done "O"
}

# ── STAGE P: Reproducibility + Confidence ─────────────────────────────────────
stage_P() {
    run_stage "P" "Reproducibility + Confidence Score" || return 0
    local prs_cal="prs/population_calibrated_v2.csv"
    [[ ! -f "$prs_cal" ]] && prs_cal="prs/population_calibrated.csv"
    if [[ -f "$prs_cal" ]] && [[ -f prs/prs_uncertainty.csv ]]; then
        run_python "$SCRIPT_DIR/scripts/reproducibility.py" \
            --pipeline-dir "$PIPELINE_DIR" --sample-id "$SAMPLE_ID" \
            --prs-data "$prs_cal" --uncertainty-data prs/prs_uncertainty.csv \
            --consistency-report prs/consistency_check_report.json \
            --snp-db "$SNP_DB" --output-dir reports/ \
            $([[ "$REPRODUCIBLE" == true ]] && echo "--reproducible --seed $SEED")
    fi
    mark_done "P"
}

# ── STAGE Q: Phase 5 Validation Layer ─────────────────────────────────────────
stage_Q() {
    run_stage "Q" "Phase 5 Validation Layer" || return 0
    if [[ "$SKIP_PHASE5" == true ]]; then
        info "Phase 5 validation skipped (--skip-phase5)"
        return 0
    fi

    local prs_cal="prs/population_calibrated_v2.csv"
    [[ ! -f "$prs_cal" ]] && prs_cal="prs/population_calibrated.csv"

    # Coverage audit
    if [[ -f "$SNP_DB" ]] && [[ -f "$PIPELINE_DIR/scripts/coverage_audit.py" ]]; then
        section "Coverage Audit"
        run_python "$SCRIPT_DIR/scripts/coverage_audit.py" \
            --snp-db "$SNP_DB" --bfile plink/ld_pruned_dataset --output-dir validation/
    fi

    # Evidence scoring
    if [[ -f "$SNP_DB" ]] && [[ -f "$PIPELINE_DIR/scripts/evidence_scoring.py" ]]; then
        section "Evidence Scoring"
        local anc_json="ancestry/posterior_probabilities.json"
        [[ ! -f "$anc_json" ]] && anc_json="pca/ancestry_inference.json"
        run_python "$SCRIPT_DIR/scripts/evidence_scoring.py" \
            --snp-db "$SNP_DB" --ancestry "$anc_json" \
            --coverage-csv validation/coverage_audit.csv \
            --uncertainty-csv prs/prs_uncertainty.csv --output-dir validation/
    fi

    # Clinical readiness
    if [[ -f "$PIPELINE_DIR/scripts/clinical_readiness.py" ]]; then
        section "Clinical Readiness"
        local anc_json="ancestry/posterior_probabilities.json"
        [[ ! -f "$anc_json" ]] && anc_json="pca/ancestry_inference.json"
        run_python "$SCRIPT_DIR/scripts/clinical_readiness.py" \
            --confidence-json reports/confidence_score.json \
            --evidence-json validation/evidence_scores.json \
            --ancestry-json "$anc_json" \
            --coverage-csv validation/coverage_audit.csv \
            --uncertainty-csv prs/prs_uncertainty.csv --output-dir validation/
    fi

    # Limitations
    if [[ -f "$PIPELINE_DIR/scripts/limitations_engine.py" ]]; then
        section "Limitations Engine"
        local anc_json="ancestry/posterior_probabilities.json"
        [[ ! -f "$anc_json" ]] && anc_json="pca/ancestry_inference.json"
        run_python "$SCRIPT_DIR/scripts/limitations_engine.py" \
            --coverage-csv validation/coverage_audit.csv \
            --evidence-json validation/evidence_scores.json \
            --ancestry-json "$anc_json" \
            --uncertainty-csv prs/prs_uncertainty.csv --output-dir validation/
    fi

    # RQI
    if [[ -f "$PIPELINE_DIR/scripts/research_quality_index.py" ]]; then
        section "Research Quality Index"
        run_python "$SCRIPT_DIR/scripts/research_quality_index.py" \
            --confidence-json reports/confidence_score.json \
            --evidence-json validation/evidence_scores.json \
            --concordance-json validation/concordance_report.json \
            --coverage-csv validation/coverage_audit.csv \
            --readiness-json validation/clinical_readiness.json --output-dir validation/
    fi

    # Validation dashboard
    if [[ -f "$PIPELINE_DIR/scripts/validation_dashboard.py" ]] && [[ -f "$prs_cal" ]]; then
        section "Validation Dashboard"
        local anc_json="ancestry/classification_report.json"
        [[ ! -f "$anc_json" ]] && anc_json="pca/ancestry_inference.json"
        run_python "$SCRIPT_DIR/scripts/validation_dashboard.py" \
            --ancestry-json "$anc_json" --evidence-json validation/evidence_scores.json \
            --readiness-json validation/clinical_readiness.json \
            --rqi-json validation/research_quality_index.json \
            --coverage-csv validation/coverage_audit.csv \
            --uncertainty-csv prs/prs_uncertainty.csv --output-dir reports/
    fi

    mark_done "Q"
}

# ── STAGE R: Phase 6 Scientific Corrections ────────────────────────────────────
stage_R() {
    run_stage "R" "Phase 6 Scientific Corrections" || return 0
    if [[ "$SKIP_PHASE6" == true ]]; then
        info "Phase 6 corrections skipped (--skip-phase6)"
        return 0
    fi

    # Scientific benchmark
    if [[ -f "$PIPELINE_DIR/scripts/scientific_benchmark_v2.py" ]]; then
        local prs_cal="prs/population_calibrated_v2.csv"
        [[ ! -f "$prs_cal" ]] && prs_cal="prs/population_calibrated.csv"
        if [[ -f "$prs_cal" ]]; then
            section "Scientific Benchmarking"
            run_python "$SCRIPT_DIR/scripts/scientific_benchmark_v2.py" \
                --prs-data "$prs_cal" --snp-db "$SNP_DB" \
                --coverage-csv validation/coverage_audit.csv --output-dir validation/
        fi
    fi

    # RQS v2
    if [[ -f "$PIPELINE_DIR/scripts/research_quality_v2.py" ]]; then
        section "Research Quality Score v2"
        local anc_json="ancestry/quality_metrics.json"
        local conc_json="prs/method_concordance_v2.json"
        local cal_csv="prs/population_calibrated_v2.csv"
        [[ ! -f "$cal_csv" ]] && cal_csv="prs/population_calibrated.csv"

        run_python "$SCRIPT_DIR/scripts/research_quality_v2.py" \
            --ancestry-json "$anc_json" --concordance-json "$conc_json" \
            --coverage-csv validation/coverage_audit.csv \
            --benchmark-json validation/scientific_benchmark.json \
            --manifest-json reports/run_manifest.json \
            --calibration-csv "$cal_csv" --output-dir research/
    fi

    mark_done "R"
}

# ═══════════════════════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════════════════════
main() {
    parse_args "$@"
    auto_detect

    echo ""
    echo -e "${BOLD}${CYAN}╔══════════════════════════════════════════════════════════════╗${NC}"
    echo -e "${BOLD}${CYAN}║   PRS RESEARCH PLATFORM v6.0.0 — Phase 6 Corrections Active  ║${NC}"
    echo -e "${BOLD}${CYAN}╚══════════════════════════════════════════════════════════════╝${NC}"
    echo -e "  Started:    $(date)"
    echo -e "  PLINK:      $PLINK_BIN"
    echo -e "  Python:     $($PYTHON_BIN --version 2>&1)"
    echo -e "  Full 1000G: $HAS_FULL_1000G"
    echo -e "  Phase 6:    PCA=$HAS_PHASE6_PCA Ancestry=$HAS_PHASE6_ANCESTRY Cal=$HAS_PHASE6_CALIBRATION"

    # Setup directories
    cd "$OUT_DIR"
    for d in plink qc pca gwas prs interpretations ancestry reports validation research references; do
        mkdir -p "$d"
    done

    # If reproducible mode
    [[ "$REPRODUCIBLE" == true ]] && export PYTHONHASHSEED="$SEED"

    # Interactive mode: ask for inputs and stage selection
    if [[ "$INTERACTIVE" == true ]] && [[ "$AUTO_MODE" == false ]]; then
        ask_input_vcf || exit 1
        ask_sample_id
        show_menu
        select_stage
    else
        # Auto mode: run all stages
        [[ -z "$INPUT_VCF" ]] && { err "--input-vcf required in auto mode"; exit 1; }
        [[ ! -f "$INPUT_VCF" ]] && { err "VCF not found: $INPUT_VCF"; exit 1; }
        RUN_STAGES="A B C D E F G H I J K L M N O P Q R"
        ok "Auto mode: $SAMPLE_ID ($INPUT_VCF)"
    fi

    # Execute selected stages
    for stage in $RUN_STAGES; do
        case "$stage" in
            A) stage_A ;;  B) stage_B ;;  C) stage_C ;;  D) stage_D ;;
            E) stage_E ;;  F) stage_F ;;  G) stage_G ;;  H) stage_H ;;
            I) stage_I ;;  J) stage_J ;;  K) stage_K ;;  L) stage_L ;;
            M) stage_M ;;  N) stage_N ;;  O) stage_O ;;  P) stage_P ;;
            Q) stage_Q ;;  R) stage_R ;;
            *) warn "Unknown stage: $stage" ;;
        esac
    done

    # ── FINAL SUMMARY ─────────────────────────────────────────────────────
    echo ""
    echo -e "${BOLD}${CYAN}╔══════════════════════════════════════════════════════════════╗${NC}"
    echo -e "${BOLD}${CYAN}║   PIPELINE COMPLETE — v6.0.0                                  ║${NC}"
    echo -e "${BOLD}${CYAN}╚══════════════════════════════════════════════════════════════╝${NC}"
    echo -e "  Completed: $(date)"
    echo ""
    echo -e "  ${BOLD}Outputs:${NC}"

    local files_to_check=(
        "reports/report_en.html"
        "reports/report_es.html"
        "reports/dashboard_en.html"
        "reports/run_manifest.json"
        "reports/confidence_score.json"
        "prs/population_calibrated_v2.csv"
        "prs/population_calibrated.csv"
        "prs/prs_uncertainty.csv"
        "prs/prs_all_methods.csv"
        "prs/method_comparison_v2.csv"
        "interpretations/interpretations_en.json"
        "interpretations/interpretations_es.json"
        "ancestry/classification_report.json"
        "ancestry/admixture_results.json"
        "validation/coverage_audit.csv"
        "validation/evidence_scores.json"
        "validation/scientific_benchmark.json"
        "research/research_quality_v2.json"
        "pca/reference_pca_model.pkl"
        "pca/projected_sample.csv"
    )

    for f in "${files_to_check[@]}"; do
        if [[ -f "$f" ]]; then
            ok "$f"
        fi
    done

    echo ""
    echo -e "  ${BOLD}View:${NC}"
    echo -e "    open reports/report_en.html"
    echo -e "    open reports/dashboard_en.html"
    echo -e "    open research/research_quality_v2.md"
    echo -e "    open validation/scientific_benchmark.md"
    echo ""
    echo -e "  ${GREEN}Platform v6.0.0 — Phase 6 Scientific Corrections Active${NC}"
}

# ═══════════════════════════════════════════════════════════════════════════════
main "$@"
