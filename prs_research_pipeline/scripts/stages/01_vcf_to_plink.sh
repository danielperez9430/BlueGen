#!/usr/bin/env bash
# ============================================================================
# STAGE A: VCF → PLINK Conversion
# ============================================================================
# Converts DeepVariant VCF (GRCh37/hg19) to PLINK binary format (.bed/.bim/.fam)
#
# Usage:
#   ./01_vcf_to_plink.sh --input-vcf input.vcf.gz --out-dir plink/
#
# Output:
#   plink/cohort.bed   — genotype matrix (binary)
#   plink/cohort.bim   — variant map (chr, rsid, cm, pos, ref, alt)
#   plink/cohort.fam   — sample information (FID, IID, father, mother, sex, pheno)
# ============================================================================

set -euo pipefail

# ── Defaults ──────────────────────────────────────────────────────────────
INPUT_VCF=""
OUT_DIR="plink"
PLINK_BIN="plink"
THREADS=4
MEMORY=16000
MIN_GQ=20
MIN_DP=10
REMOVE_REF_CALLS=true

# ── Parse Arguments ───────────────────────────────────────────────────────
while [[ $# -gt 0 ]]; do
    case "$1" in
        --input-vcf)    INPUT_VCF="$2";    shift 2 ;;
        --out-dir)      OUT_DIR="$2";      shift 2 ;;
        --plink)        PLINK_BIN="$2";    shift 2 ;;
        --threads)      THREADS="$2";      shift 2 ;;
        --memory)       MEMORY="$2";       shift 2 ;;
        --min-gq)       MIN_GQ="$2";       shift 2 ;;
        --min-dp)       MIN_DP="$2";       shift 2 ;;
        --keep-ref-calls) REMOVE_REF_CALLS=false; shift ;;
        -h|--help)
            echo "Usage: $0 --input-vcf <file.vcf.gz> [--out-dir plink/]"
            echo ""
            echo "STAGE A: VCF → PLINK Conversion"
            echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
            echo "  --input-vcf      DeepVariant VCF (bgzip compressed)"
            echo "  --out-dir        Output directory for PLINK files (default: plink/)"
            echo "  --plink          Path to PLINK binary (default: plink)"
            echo "  --threads        Threads for PLINK (default: 4)"
            echo "  --memory         Memory in MB (default: 16000)"
            echo "  --min-gq         Min genotype quality (default: 20)"
            echo "  --min-dp         Min read depth (default: 10)"
            echo "  --keep-ref-calls Keep RefCall entries (default: remove them)"
            exit 0
            ;;
        *) echo "Unknown argument: $1"; exit 1 ;;
    esac
done

# ── Validation ────────────────────────────────────────────────────────────
if [[ -z "$INPUT_VCF" ]]; then
    echo "ERROR: --input-vcf is required."
    exit 1
fi
if [[ ! -f "$INPUT_VCF" ]]; then
    echo "ERROR: Input VCF not found: $INPUT_VCF"
    exit 1
fi

# ── Setup ─────────────────────────────────────────────────────────────────
mkdir -p "$OUT_DIR"
PREFIX="${OUT_DIR}/cohort"
LOG="${OUT_DIR}/01_vcf_to_plink.log"

echo "=============================================" | tee "$LOG"
echo "  STAGE A: VCF → PLINK Conversion"            | tee -a "$LOG"
echo "=============================================" | tee -a "$LOG"
echo "  Input:    $INPUT_VCF"                       | tee -a "$LOG"
echo "  Output:   $PREFIX.{bed,bim,fam}"            | tee -a "$LOG"
echo "  PLINK:    $PLINK_BIN"                       | tee -a "$LOG"
echo "  Threads:  $THREADS"                         | tee -a "$LOG"
echo "  Memory:   ${MEMORY}MB"                      | tee -a "$LOG"
echo "  Min GQ:   $MIN_GQ"                          | tee -a "$LOG"
echo "  Min DP:   $MIN_DP"                          | tee -a "$LOG"
echo "=============================================" | tee -a "$LOG"

# ── Step 1: Pre-filter VCF (Quality) ─────────────────────────────────────
# Remove RefCall entries and low-quality variants before PLINK conversion
FILTERED_VCF="${OUT_DIR}/input_filtered.vcf.gz"

echo "" | tee -a "$LOG"
echo "── Step 1: Pre-filtering VCF ──" | tee -a "$LOG"

# Use bcftools for filtering if available, otherwise pass through
if command -v bcftools &> /dev/null; then
    echo "  Using bcftools for VCF pre-filtering..." | tee -a "$LOG"

    FILTER_EXPR="GT!='./.' && GT!='.'"

    # Check if this is a DeepVariant VCF (has RefCall filter + GQ/DP FORMAT)
    HAS_REF_CALL=$(bcftools view -h "$INPUT_VCF" 2>/dev/null | grep -c "RefCall" || true)
    HAS_GQ=$(bcftools view -h "$INPUT_VCF" 2>/dev/null | grep -c "ID=GQ" || true)

    if [[ "$HAS_REF_CALL" -gt 0 ]] && [[ "$REMOVE_REF_CALLS" == "true" ]]; then
        # DeepVariant VCF: remove RefCall entries
        bcftools view -f PASS -e 'FILTER="RefCall"' "$INPUT_VCF" \
            | bcftools view -i "GQ>=${MIN_GQ} && DP>=${MIN_DP}" \
            -Oz -o "$FILTERED_VCF"
    elif [[ "$HAS_GQ" -gt 0 ]]; then
        # Has GQ/DP FORMAT fields (DeepVariant-like)
        bcftools view -i "GQ>=${MIN_GQ} && DP>=${MIN_DP}" \
            -Oz -o "$FILTERED_VCF" "$INPUT_VCF"
    else
        # Standard VCF (e.g., 1000 Genomes): no per-sample quality filters
        # Just remove sites with no-calls in all samples
        bcftools view -c 1 -Oz -o "$FILTERED_VCF" "$INPUT_VCF"
    fi

    # Index
    bcftools index -t "$FILTERED_VCF"

    VAR_BEFORE=$(bcftools view -H "$INPUT_VCF" 2>/dev/null | wc -l | tr -d ' ')
    VAR_AFTER=$(bcftools view -H "$FILTERED_VCF" 2>/dev/null | wc -l | tr -d ' ')
    echo "  Variants before filtering: $VAR_BEFORE" | tee -a "$LOG"
    echo "  Variants after filtering:  $VAR_AFTER"  | tee -a "$LOG"
else
    echo "  bcftools not found — using VCF directly (install: brew install htslib)" | tee -a "$LOG"
    FILTERED_VCF="$INPUT_VCF"
fi

# ── Step 2: VCF → PLINK Conversion ──────────────────────────────────────
echo "" | tee -a "$LOG"
echo "── Step 2: Converting VCF to PLINK binary ──" | tee -a "$LOG"

# Check PLINK version (1.9 vs 2.0 have different flags)
PLINK_VERSION=$("$PLINK_BIN" --version 2>&1 | head -1 || echo "unknown")
echo "  PLINK version: $PLINK_VERSION" | tee -a "$LOG"

# PLINK 2.0 syntax
if echo "$PLINK_VERSION" | grep -q "PLINK v2"; then
    "$PLINK_BIN" \
        --vcf "$FILTERED_VCF" \
        --make-bed \
        --double-id \
        --allow-extra-chr \
        --chr-set 24 \
        --out "$PREFIX" \
        --threads "$THREADS" \
        --memory "$MEMORY" \
        2>&1 | tee -a "$LOG"
else
    # PLINK 1.9 syntax
    "$PLINK_BIN" \
        --vcf "$FILTERED_VCF" \
        --make-bed \
        --double-id \
        --allow-extra-chr \
        --out "$PREFIX" \
        --threads "$THREADS" \
        --memory "$MEMORY" \
        2>&1 | tee -a "$LOG"
fi

# Rename . variant IDs to chr:pos for compatibility with 1000G reference
echo "" | tee -a "$LOG"
echo "-- Renaming variants to chr:pos --" | tee -a "$LOG"
awk -v OFS="\t" '{print $1, $1":"$4, $3, $4, $5, $6}' "${PREFIX}.bim" > "${PREFIX}.bim.tmp"
mv "${PREFIX}.bim.tmp" "${PREFIX}.bim"
echo "  Variants renamed to chr:pos format" | tee -a "$LOG"


# ── Step 3: Verify Output ────────────────────────────────────────────────
echo "" | tee -a "$LOG"
echo "── Step 3: Verifying output ──" | tee -a "$LOG"

for ext in bed bim fam; do
    if [[ -f "${PREFIX}.${ext}" ]]; then
        SIZE=$(wc -c < "${PREFIX}.${ext}" | tr -d ' ')
        echo "  ✅ ${PREFIX}.${ext} (${SIZE} bytes)" | tee -a "$LOG"
    else
        echo "  ❌ ${PREFIX}.${ext} MISSING" | tee -a "$LOG"
        exit 1
    fi
done

# Summary stats
N_VARIANTS=$(wc -l < "${PREFIX}.bim" | tr -d ' ')
N_SAMPLES=$(wc -l < "${PREFIX}.fam" | tr -d ' ')
echo "" | tee -a "$LOG"
echo "  📊 $N_VARIANTS variants, $N_SAMPLES samples" | tee -a "$LOG"
echo "" | tee -a "$LOG"
echo "✅ STAGE A complete: PLINK dataset ready." | tee -a "$LOG"
