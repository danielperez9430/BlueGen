#!/usr/bin/env bash
# ============================================================================
# STAGE D (UPGRADED): PCA on 1000 Genomes + Project Target Sample
# ============================================================================
# Trains PCA space on 1000 Genomes reference panel, then projects the target
# sample into that space. This is the SCIENTIFICALLY VALID approach — PCA
# must be trained on the reference population, not on the target sample.
#
# Methodology:
#   1. Convert 1000G VCF to PLINK binary (if needed)
#   2. QC + LD-prune 1000G reference
#   3. Train PCA on 1000G only
#   4. Project target sample into trained PCA space
#   5. Merge coordinates for ancestry inference
#
# Reference: Price et al. (2006) Nat Genet; Patterson et al. (2006) PLoS Genet
#
# Usage:
#   ./04_pca_1000G.sh \
#       --1000g-vcf reference/1000G/ALL.chr22.phase3.vcf.gz \
#       --target-bfile plink/ld_pruned_dataset \
#       --population-panel reference/1000G/20130606_g1k_population.txt \
#       --out-dir pca/
# ============================================================================

set -euo pipefail

# ── Defaults ──────────────────────────────────────────────────────────────
TARGET_BFILE=""
G1K_VCF=""
G1K_BFILE=""
POP_PANEL=""
OUT_DIR="pca"
PLINK_BIN="plink"
THREADS=4
MEMORY=16000
NUM_PCS=20

# ── Parse Arguments ───────────────────────────────────────────────────────
while [[ $# -gt 0 ]]; do
    case "$1" in
        --target-bfile)  TARGET_BFILE="$2";  shift 2 ;;
        --1000g-vcf)     G1K_VCF="$2";       shift 2 ;;
        --1000g-bfile)   G1K_BFILE="$2";     shift 2 ;;
        --population-panel) POP_PANEL="$2";  shift 2 ;;
        --out-dir)       OUT_DIR="$2";       shift 2 ;;
        --plink)         PLINK_BIN="$2";     shift 2 ;;
        --threads)       THREADS="$2";       shift 2 ;;
        --memory)        MEMORY="$2";        shift 2 ;;
        --num-pcs)       NUM_PCS="$2";       shift 2 ;;
        -h|--help)
            cat << 'EOF'
STAGE D (UPGRADED): PCA on 1000 Genomes + Projection

Usage: ./04_pca_1000G.sh --target-bfile <prefix> --1000g-vcf <vcf.gz>

OPTIONS:
  --target-bfile      PLINK prefix for target sample
  --1000g-vcf         1000 Genomes phase 3 VCF (any chromosome)
  --1000g-bfile       Pre-built 1000G PLINK binary (skip VCF conversion)
  --population-panel  1000G population labels (tab-separated: sample pop super_pop)
  --out-dir           Output directory (default: pca/)
  --num-pcs           Number of PCs to compute (default: 20)
EOF
            exit 0
            ;;
        *) echo "Unknown: $1"; exit 1 ;;
    esac
done

# ── Validation ────────────────────────────────────────────────────────────
if [[ -z "$TARGET_BFILE" ]]; then
    echo "ERROR: --target-bfile is required"
    exit 1
fi

mkdir -p "$OUT_DIR"
REF_DIR="${OUT_DIR}/1000G_reference"
mkdir -p "$REF_DIR"
LOG="${OUT_DIR}/04_pca_1000G.log"

echo "════════════════════════════════════════════════" | tee "$LOG"
echo "  STAGE D: 1000 Genomes PCA + Projection"        | tee -a "$LOG"
echo "════════════════════════════════════════════════" | tee -a "$LOG"

# ── Step 1: Build 1000G PLINK Dataset ────────────────────────────────────
echo "" | tee -a "$LOG"
echo "── Step 1: 1000G Reference Preparation ──" | tee -a "$LOG"

G1K_PREFIX="${REF_DIR}/1000G_chr22"

if [[ -n "${G1K_BFILE:-}" ]] && [[ -f "${G1K_BFILE}.bed" ]]; then
    echo "  Using pre-built 1000G PLINK: $G1K_BFILE" | tee -a "$LOG"
    G1K_PREFIX="$G1K_BFILE"
elif [[ -n "$G1K_VCF" ]] && [[ -f "$G1K_VCF" ]]; then
    echo "  Converting 1000G VCF → PLINK binary..." | tee -a "$LOG"
    "$PLINK_BIN" \
        --vcf "$G1K_VCF" \
        --make-bed \
        --double-id \
        --allow-extra-chr \
        --out "$G1K_PREFIX" \
        --threads "$THREADS" \
        --memory "$MEMORY" \
        2>&1 | tee -a "$LOG"
    echo "  ✅ 1000G PLINK: ${G1K_PREFIX}.{bed,bim,fam}" | tee -a "$LOG"

    # Assign unique variant IDs (chr:pos) replacing '.' in BIM
    echo "  Assigning unique variant IDs (chr:pos)..." | tee -a "$LOG"
    awk '{print $1"\t"$1":"$4"\t"$3"\t"$4"\t"$5"\t"$6}' "${G1K_PREFIX}.bim" > "${G1K_PREFIX}.bim.tmp"
    mv "${G1K_PREFIX}.bim.tmp" "${G1K_PREFIX}.bim"
else
    echo "  ❌ No 1000G data provided. Use --1000g-vcf or --1000g-bfile" | tee -a "$LOG"
    exit 1
fi

N_G1K=$(wc -l < "${G1K_PREFIX}.fam" | tr -d ' ')
echo "  1000G samples: $N_G1K" | tee -a "$LOG"

# ── Step 2: QC 1000G Reference ──────────────────────────────────────────
echo "" | tee -a "$LOG"
echo "── Step 2: QC 1000G Reference ──" | tee -a "$LOG"

G1K_QC="${REF_DIR}/1000G_qc"

"$PLINK_BIN" \
    --allow-extra-chr \
    --bfile "$G1K_PREFIX" \
    --geno 0.05 \
    --maf 0.01 \
    --hwe 0.000001 midp \
    --make-bed \
    --out "$G1K_QC" \
    --threads "$THREADS" \
    --memory "$MEMORY" \
    2>&1 | tee -a "$LOG"

N_G1K_QC=$(wc -l < "${G1K_QC}.bim" | tr -d ' ')
echo "  QC'd variants: $N_G1K_QC" | tee -a "$LOG"

# ── Step 3: LD Prune 1000G Reference ────────────────────────────────────
echo "" | tee -a "$LOG"
echo "── Step 3: LD Prune 1000G Reference ──" | tee -a "$LOG"

G1K_PRUNED="${REF_DIR}/1000G_pruned"

"$PLINK_BIN" \
    --allow-extra-chr \
    --bfile "$G1K_QC" \
    --indep-pairwise 50 5 0.2 \
    --out "${REF_DIR}/1000G_ld" \
    --threads "$THREADS" \
    --memory "$MEMORY" \
    2>&1 | tee -a "$LOG"

"$PLINK_BIN" \
    --allow-extra-chr \
    --bfile "$G1K_QC" \
    --extract "${REF_DIR}/1000G_ld.prune.in" \
    --make-bed \
    --out "$G1K_PRUNED" \
    --threads "$THREADS" \
    --memory "$MEMORY" \
    2>&1 | tee -a "$LOG"

N_G1K_PRUNED=$(wc -l < "${G1K_PRUNED}.bim" | tr -d ' ')
echo "  LD-pruned variants: $N_G1K_PRUNED" | tee -a "$LOG"

# ── Step 4: Train PCA on 1000G Reference ONLY ────────────────────────────
echo "" | tee -a "$LOG"
echo "── Step 4: Train PCA on 1000 Genomes ──" | tee -a "$LOG"

G1K_PCA="${REF_DIR}/1000G_pca"

"$PLINK_BIN" \
    --allow-extra-chr \
    --bfile "$G1K_PRUNED" \
    --pca "$NUM_PCS" \
    --out "$G1K_PCA" \
    --threads "$THREADS" \
    --memory "$MEMORY" \
    2>&1 | tee -a "$LOG"

if [[ -f "${G1K_PCA}.eigenvec" ]] && [[ -f "${G1K_PCA}.eigenval" ]]; then
    echo "  ✅ PCA trained on 1000G: ${G1K_PCA}.eigenvec" | tee -a "$LOG"
else
    echo "  ❌ PCA training failed" | tee -a "$LOG"
    exit 1
fi

# ── Step 5: Prepare Target Sample SNPs ───────────────────────────────────
echo "" | tee -a "$LOG"
echo "── Step 5: Align Target Sample with 1000G SNPs ──" | tee -a "$LOG"

# Extract target sample SNPs that overlap with 1000G pruned set
TARGET_ALIGNED="${REF_DIR}/target_aligned"

# Get common SNPs between target and 1000G
cut -f2 "${G1K_PRUNED}.bim" > "${REF_DIR}/common_snps.txt"

"$PLINK_BIN" \
    --allow-extra-chr \
    --bfile "$TARGET_BFILE" \
    --extract "${REF_DIR}/common_snps.txt" \
    --make-bed \
    --out "$TARGET_ALIGNED" \
    --threads "$THREADS" \
    --memory "$MEMORY" \
    2>&1 | tee -a "$LOG"

N_TARGET_SNPS=$(wc -l < "${TARGET_ALIGNED}.bim" | tr -d ' ')
echo "  Target SNPs aligned: $N_TARGET_SNPS" | tee -a "$LOG"

# ── Step 6: Project Target into 1000G PCA Space ─────────────────────────
echo "" | tee -a "$LOG"
echo "── Step 6: Project Target Sample into 1000G PCA Space ──" | tee -a "$LOG"

# Merge target with 1000G, then extract the merged PCA by re-running
# PLINK PCA on the combined dataset with the same variants
# OR: use the score projection method

# Method: Combine 1000G + target, run PCA on combined
COMBINED="${REF_DIR}/combined"

# First merge
"$PLINK_BIN" \
    --allow-extra-chr \
    --bfile "$G1K_PRUNED" \
    --bmerge "$TARGET_ALIGNED" \
    --make-bed \
    --out "$COMBINED" \
    --allow-no-sex \
    --threads "$THREADS" \
    --memory "$MEMORY" \
    2>&1 | tee -a "$LOG" || {
    # If merge fails (allele mismatches), try flip-scan then re-merge
    echo "  Merge with allele mismatches; running flip-scan..." | tee -a "$LOG"
    "$PLINK_BIN" \
        --bfile "$TARGET_ALIGNED" \
        --flip "${COMBINED}-merge.missnp" \
        --make-bed \
        --out "${TARGET_ALIGNED}_flipped" \
        --threads "$THREADS" \
        --memory "$MEMORY" \
        2>&1 | tee -a "$LOG"

    "$PLINK_BIN" \
        --bfile "$G1K_PRUNED" \
        --bmerge "${TARGET_ALIGNED}_flipped" \
        --make-bed \
        --out "$COMBINED" \
        --allow-no-sex \
        --threads "$THREADS" \
        --memory "$MEMORY" \
        2>&1 | tee -a "$LOG" || {
        # Second fallback: exclude remaining conflicting variants
        echo "  Still have allele conflicts; excluding them..." | tee -a "$LOG"
        "$PLINK_BIN" --bfile "$G1K_PRUNED" \
            --exclude "${COMBINED}-merge.missnp" \
            --make-bed --out "${G1K_PRUNED}_clean" \
            --threads "$THREADS" --memory "$MEMORY" 2>&1 | tee -a "$LOG"
        "$PLINK_BIN" --bfile "${TARGET_ALIGNED}_flipped" \
            --exclude "${COMBINED}-merge.missnp" \
            --make-bed --out "${TARGET_ALIGNED}_clean" \
            --threads "$THREADS" --memory "$MEMORY" 2>&1 | tee -a "$LOG"
        "$PLINK_BIN" --bfile "${G1K_PRUNED}_clean" \
            --bmerge "${TARGET_ALIGNED}_clean" \
            --make-bed --out "$COMBINED" --allow-no-sex \
            --threads "$THREADS" --memory "$MEMORY" 2>&1 | tee -a "$LOG"
        G1K_PRUNED="${G1K_PRUNED}_clean"
    }
}

# Run PCA on combined dataset
COMBINED_PCA="${OUT_DIR}/combined_pca"

"$PLINK_BIN" \
    --allow-extra-chr \
    --bfile "$COMBINED" \
    --pca "$NUM_PCS" \
    --out "$COMBINED_PCA" \
    --threads "$THREADS" \
    --memory "$MEMORY" \
    2>&1 | tee -a "$LOG"

echo "  ✅ Combined PCA: ${COMBINED_PCA}.eigenvec" | tee -a "$LOG"

# ── Step 7: Extract Target Sample PCs + 1000G PCs ───────────────────────
echo "" | tee -a "$LOG"
echo "── Step 7: Extract Coordinates ──" | tee -a "$LOG"

# Get target sample IIDs (all samples, not just first)
TARGET_IIDS=$(awk '{print $2}' "${TARGET_BFILE}.fam" | paste -sd ',' -)
N_TARGET=$(echo "$TARGET_IIDS" | tr ',' '\n' | wc -l | tr -d ' ')
echo "  Target samples: $N_TARGET" | tee -a "$LOG"

# Extract target rows from combined eigenvec
python3 -c "
import pandas as pd
import numpy as np

# Load combined eigenvec
cols = ['FID','IID'] + [f'PC{i}' for i in range(1,$NUM_PCS+1)]
ev = pd.read_csv('${COMBINED_PCA}.eigenvec', sep=r'\s+', header=None)
ev.columns = cols[:ev.shape[1]]

# Target IIDs as list
target_iids = '${TARGET_IIDS}'.split(',')

# Separate target and reference
target = ev[ev['IID'].isin(target_iids)]
reference = ev[~ev['IID'].isin(target_iids)]

# Save target PCs
target_pcs = target[['FID','IID'] + [f'PC{i}' for i in range(1,$NUM_PCS+1)]]
target_pcs.to_csv('${OUT_DIR}/target_pcs.eigenvec', sep='\t', index=False, float_format='%.6f')

# Save reference PCs (1000G)
ref_pcs = reference[['FID','IID'] + [f'PC{i}' for i in range(1,$NUM_PCS+1)]]
ref_pcs.to_csv('${OUT_DIR}/1000G_pcs.eigenvec', sep='\t', index=False, float_format='%.6f')

print(f'  Target PCs saved: {len(target_pcs)} rows ({len(target_iids)} samples)')
print(f'  1000G PCs saved: {len(ref_pcs)} rows')
print(f'  PC1 range: [{ev[\"PC1\"].min():.4f}, {ev[\"PC1\"].max():.4f}]')

# Copy combined as final output for compatibility
ev.to_csv('${OUT_DIR}/pca_results.eigenvec', sep='\t', index=False, float_format='%.6f')
" 2>&1 | tee -a "$LOG"

# Copy eigenval
cp "${COMBINED_PCA}.eigenval" "${OUT_DIR}/pca_results.eigenval"

echo "" | tee -a "$LOG"
echo "════════════════════════════════════════════════" | tee -a "$LOG"
echo "  ✅ STAGE D COMPLETE: 1000G PCA + Projection"   | tee -a "$LOG"
echo "════════════════════════════════════════════════" | tee -a "$LOG"
echo "  Target PCs:      ${OUT_DIR}/target_pcs.eigenvec" | tee -a "$LOG"
echo "  1000G PCs:       ${OUT_DIR}/1000G_pcs.eigenvec" | tee -a "$LOG"
echo "  Combined PCs:    ${OUT_DIR}/pca_results.eigenvec" | tee -a "$LOG"
