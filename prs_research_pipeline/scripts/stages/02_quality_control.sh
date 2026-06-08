#!/usr/bin/env bash
# ============================================================================
# STAGE B: Quality Control (QC)
# ============================================================================
# Applies standard GWAS QC filters to PLINK dataset:
#   - SNP missingness (--geno)
#   - Individual missingness (--mind)
#   - Minor allele frequency (--maf)
#   - Hardy-Weinberg equilibrium (--hwe)
#
# Usage:
#   ./02_quality_control.sh --bfile plink/cohort --out-dir qc/
#
# Output:
#   qc/qc_filtered.{bed,bim,fam}  — QC-passed PLINK dataset
#   qc/qc_report.txt              — QC statistics
# ============================================================================

set -euo pipefail

# ── Defaults ──────────────────────────────────────────────────────────────
BFILE=""
OUT_DIR="qc"
PLINK_BIN="plink"
THREADS=4
MEMORY=16000

# QC thresholds (standard GWAS defaults)
GENO=0.1        # SNP missingness: exclude SNPs missing in >10% of samples
MIND=0.1        # Individual missingness: exclude individuals missing >10% of SNPs
MAF=0.01        # Minor allele frequency: exclude SNPs with MAF < 1%
HWE=0.000001    # Hardy-Weinberg: exclude SNPs with HWE p < 1e-6

# ── Parse Arguments ───────────────────────────────────────────────────────
while [[ $# -gt 0 ]]; do
    case "$1" in
        --bfile)       BFILE="$2";       shift 2 ;;
        --out-dir)     OUT_DIR="$2";     shift 2 ;;
        --plink)       PLINK_BIN="$2";   shift 2 ;;
        --threads)     THREADS="$2";     shift 2 ;;
        --memory)      MEMORY="$2";      shift 2 ;;
        --geno)        GENO="$2";        shift 2 ;;
        --mind)        MIND="$2";        shift 2 ;;
        --maf)         MAF="$2";         shift 2 ;;
        --hwe)         HWE="$2";         shift 2 ;;
        -h|--help)
            echo "Usage: $0 --bfile <prefix> [--out-dir qc/]"
            echo ""
            echo "STAGE B: Quality Control"
            echo "━━━━━━━━━━━━━━━━━━━━━━━━━━"
            echo "  --bfile      Input PLINK prefix"
            echo "  --out-dir    Output directory (default: qc/)"
            echo "  --geno       SNP missingness threshold (default: 0.1)"
            echo "  --mind       Individual missingness threshold (default: 0.1)"
            echo "  --maf        MAF threshold (default: 0.01)"
            echo "  --hwe        HWE p-value threshold (default: 1e-6)"
            exit 0
            ;;
        *) echo "Unknown argument: $1"; exit 1 ;;
    esac
done

# ── Validation ────────────────────────────────────────────────────────────
if [[ -z "$BFILE" ]]; then
    echo "ERROR: --bfile is required."
    exit 1
fi
for ext in bed bim fam; do
    if [[ ! -f "${BFILE}.${ext}" ]]; then
        echo "ERROR: Missing PLINK file: ${BFILE}.${ext}"
        exit 1
    fi
done

# ── Setup ─────────────────────────────────────────────────────────────────
mkdir -p "$OUT_DIR"
PREFIX="${OUT_DIR}/qc_filtered"
LOG="${OUT_DIR}/02_quality_control.log"

echo "=============================================" | tee "$LOG"
echo "  STAGE B: Quality Control (QC)"             | tee -a "$LOG"
echo "=============================================" | tee -a "$LOG"
echo "  Input:    $BFILE"                          | tee -a "$LOG"
echo "  Output:   $PREFIX"                         | tee -a "$LOG"
echo "  --geno:   $GENO  (SNP missingness)"        | tee -a "$LOG"
echo "  --mind:   $MIND  (Individual missingness)" | tee -a "$LOG"
echo "  --maf:    $MAF  (Minor allele frequency)"  | tee -a "$LOG"
echo "  --hwe:    $HWE  (HWE p-value)"             | tee -a "$LOG"
echo "=============================================" | tee -a "$LOG"

# ── Initial Stats ─────────────────────────────────────────────────────────
echo "" | tee -a "$LOG"
echo "── Initial dataset ──" | tee -a "$LOG"
INIT_VARIANTS=$(wc -l < "${BFILE}.bim" | tr -d ' ')
INIT_SAMPLES=$(wc -l < "${BFILE}.fam" | tr -d ' ')
echo "  Variants: $INIT_VARIANTS" | tee -a "$LOG"
echo "  Samples:  $INIT_SAMPLES"  | tee -a "$LOG"

# ── Step 1: SNP Missingness ──────────────────────────────────────────────
echo "" | tee -a "$LOG"
echo "── Step 1: SNP missingness filter (--geno $GENO) ──" | tee -a "$LOG"

"$PLINK_BIN" \
    --allow-extra-chr \
    --bfile "$BFILE" \
    --geno "$GENO" \
    --make-bed \
    --out "${OUT_DIR}/tmp_geno" \
    --threads "$THREADS" \
    --memory "$MEMORY" \
    2>&1 | tee -a "$LOG"

POST_GENO=$(wc -l < "${OUT_DIR}/tmp_geno.bim" | tr -d ' ')
REMOVED=$(( INIT_VARIANTS - POST_GENO ))
echo "  SNPs removed: $REMOVED ($POST_GENO remaining)" | tee -a "$LOG"

# ── Step 2: Individual Missingness ────────────────────────────────────────
echo "" | tee -a "$LOG"
echo "── Step 2: Individual missingness filter (--mind $MIND) ──" | tee -a "$LOG"

"$PLINK_BIN" \
    --allow-extra-chr \
    --bfile "${OUT_DIR}/tmp_geno" \
    --mind "$MIND" \
    --make-bed \
    --out "${OUT_DIR}/tmp_mind" \
    --threads "$THREADS" \
    --memory "$MEMORY" \
    2>&1 | tee -a "$LOG"

PRE_SAMPLES=$(wc -l < "${OUT_DIR}/tmp_geno.fam" | tr -d ' ')
POST_SAMPLES=$(wc -l < "${OUT_DIR}/tmp_mind.fam" | tr -d ' ')
REMOVED_SAMPLES=$(( PRE_SAMPLES - POST_SAMPLES ))
echo "  Samples removed: $REMOVED_SAMPLES ($POST_SAMPLES remaining)" | tee -a "$LOG"

# ── Step 3: Minor Allele Frequency ────────────────────────────────────────
echo "" | tee -a "$LOG"
echo "── Step 3: MAF filter (--maf $MAF) ──" | tee -a "$LOG"

PRE_MAF=$(wc -l < "${OUT_DIR}/tmp_mind.bim" | tr -d ' ')

"$PLINK_BIN" \
    --allow-extra-chr \
    --bfile "${OUT_DIR}/tmp_mind" \
    --maf "$MAF" \
    --make-bed \
    --out "${OUT_DIR}/tmp_maf" \
    --threads "$THREADS" \
    --memory "$MEMORY" \
    2>&1 | tee -a "$LOG"

POST_MAF=$(wc -l < "${OUT_DIR}/tmp_maf.bim" | tr -d ' ')
REMOVED_MAF=$(( PRE_MAF - POST_MAF ))
echo "  SNPs removed (MAF < $MAF): $REMOVED_MAF ($POST_MAF remaining)" | tee -a "$LOG"

# ── Step 4: Hardy-Weinberg Equilibrium ────────────────────────────────────
echo "" | tee -a "$LOG"
echo "── Step 4: HWE filter (--hwe $HWE) ──" | tee -a "$LOG"

PRE_HWE=$(wc -l < "${OUT_DIR}/tmp_maf.bim" | tr -d ' ')

# HWE filter: use midp variant for case/control or standard for population
"$PLINK_BIN" \
    --allow-extra-chr \
    --bfile "${OUT_DIR}/tmp_maf" \
    --hwe "$HWE" \
    --make-bed \
    --out "$PREFIX" \
    --threads "$THREADS" \
    --memory "$MEMORY" \
    2>&1 | tee -a "$LOG"

POST_HWE=$(wc -l < "${PREFIX}.bim" | tr -d ' ')
REMOVED_HWE=$(( PRE_HWE - POST_HWE ))
echo "  SNPs removed (HWE p < $HWE): $REMOVED_HWE ($POST_HWE remaining)" | tee -a "$LOG"

# ── Step 5: Generate QC Report ────────────────────────────────────────────
echo "" | tee -a "$LOG"
echo "── Step 5: QC report ──" | tee -a "$LOG"

"$PLINK_BIN" \
    --allow-extra-chr \
    --bfile "$PREFIX" \
    --missing \
    --freq \
    --hardy \
    --out "${OUT_DIR}/qc_stats" \
    --threads "$THREADS" \
    --memory "$MEMORY" \
    2>&1 | tee -a "$LOG"

# ── Cleanup Temporary Files ───────────────────────────────────────────────
rm -f "${OUT_DIR}"/tmp_geno.* "${OUT_DIR}"/tmp_mind.* "${OUT_DIR}"/tmp_maf.*

# ── Final Summary ─────────────────────────────────────────────────────────
FINAL_VARIANTS=$(wc -l < "${PREFIX}.bim" | tr -d ' ')
FINAL_SAMPLES=$(wc -l < "${PREFIX}.fam" | tr -d ' ')

echo "" | tee -a "$LOG"
echo "=============================================" | tee -a "$LOG"
echo "  QC Complete"                               | tee -a "$LOG"
echo "=============================================" | tee -a "$LOG"
echo "  Initial:     $INIT_VARIANTS variants, $INIT_SAMPLES samples" | tee -a "$LOG"
echo "  After QC:    $FINAL_VARIANTS variants, $FINAL_SAMPLES samples" | tee -a "$LOG"
echo "  SNPs removed:    $(( INIT_VARIANTS - FINAL_VARIANTS ))"       | tee -a "$LOG"
echo "  Samples removed: $(( INIT_SAMPLES - FINAL_SAMPLES ))"        | tee -a "$LOG"
RETENTION=$(python3 -c "print(f'{($FINAL_VARIANTS/$INIT_VARIANTS)*100:.1f}%')")
echo "  Retention:       $RETENTION" | tee -a "$LOG"
echo "=============================================" | tee -a "$LOG"
echo "✅ STAGE B complete: QC-filtered dataset at $PREFIX" | tee -a "$LOG"
