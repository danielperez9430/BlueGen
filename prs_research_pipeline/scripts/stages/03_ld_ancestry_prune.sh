#!/usr/bin/env bash
# ============================================================================
# STAGE C (UPGRADED): Ancestry-Matched LD Pruning from 1000 Genomes
# ============================================================================
# Computes LD independently per super-population in 1000 Genomes, then
# extracts the CONSERVATIVE UNION of independent SNPs across all populations.
#
# Rationale:
#   LD patterns differ dramatically across ancestries (e.g., African populations
#   have shorter LD blocks than European). Pruning based on a single population
#   can leave correlated SNPs for individuals of other ancestries. The union
#   approach ensures a SNP is kept ONLY if it is independent in ALL populations.
#
# Method:
#   1. Split 1000 Genomes into 5 super-population strata
#   3. Union: SNP kept if it appears in ALL population prune.in sets
#      (most conservative) or in ANY (most inclusive)
#   4. Apply the union set to target sample
#
# No single-sample fallback is allowed — 1000 Genomes reference is REQUIRED.
#
# Usage:
#   ./03_ld_ancestry_prune.sh \
#       --bfile plink/qc_filtered \
#       --1000g-bfile reference/1000G/1000G_qc \
#       --population-panel reference/1000G/20130606_g1k_population.txt \
#       --out-dir plink/
# ============================================================================

set -euo pipefail

BFILE=""
G1K_BFILE=""
G1K_VCF=""
POP_PANEL=""
OUT_DIR="plink"
PLINK_BIN="plink"
THREADS=4
MEMORY=16000
WINDOW_SIZE=50
STEP_SIZE=5
R2_THRESHOLD=0.2
UNION_MODE="conservative"  # conservative = intersection (keep if independent in ALL)

while [[ $# -gt 0 ]]; do
    case "$1" in
        --bfile) BFILE="$2"; shift 2 ;;
        --1000g-bfile) G1K_BFILE="$2"; shift 2 ;;
        --1000g-vcf) G1K_VCF="$2"; shift 2 ;;
        --population-panel) POP_PANEL="$2"; shift 2 ;;
        --out-dir) OUT_DIR="$2"; shift 2 ;;
        --plink) PLINK_BIN="$2"; shift 2 ;;
        --threads) THREADS="$2"; shift 2 ;;
        --memory) MEMORY="$2"; shift 2 ;;
        --window-size) WINDOW_SIZE="$2"; shift 2 ;;
        --r2) R2_THRESHOLD="$2"; shift 2 ;;
        --union-mode) UNION_MODE="$2"; shift 2 ;;
        -h|--help)
            cat << 'EOF'
STAGE C: Ancestry-Matched LD Pruning (1000 Genomes)

Usage: ./03_ld_ancestry_prune.sh --bfile <target> --1000g-bfile <ref>

REQUIRED:
  --bfile              Target sample PLINK prefix (QC-filtered)
  --1000g-bfile        1000 Genomes PLINK binary prefix (or --1000g-vcf)
  --population-panel   1000G population labels (sample pop super_pop)

OPTIONS:
  --union-mode         conservative (intersection, default) | inclusive (union)
  --r2                 LD threshold (default: 0.2)
EOF
            exit 0 ;;
        *) echo "Unknown: $1"; exit 1 ;;
    esac
done

# ── Validation ────────────────────────────────────────────────────────
[[ -z "$BFILE" ]] && { echo "ERROR: --bfile required"; exit 1; }
mkdir -p "$OUT_DIR"

REF_DIR="${OUT_DIR}/1000G_ld_reference"
mkdir -p "$REF_DIR"
LOG="${OUT_DIR}/03_ld_ancestry_prune.log"

echo "════════════════════════════════════════════" | tee "$LOG"
echo "  STAGE C: Ancestry-Matched LD Pruning"      | tee -a "$LOG"
echo "════════════════════════════════════════════" | tee -a "$LOG"
echo "  Union mode: $UNION_MODE"                   | tee -a "$LOG"
echo "  r² threshold: $R2_THRESHOLD"               | tee -a "$LOG"

# ── Build 1000G PLINK if needed ───────────────────────────────────────
if [[ -n "${G1K_BFILE:-}" ]] && [[ -f "${G1K_BFILE}.bed" ]]; then
    G1K_PREFIX="$G1K_BFILE"
    # Check for cached ancestry-pruned SNP list from this reference
    CACHED_PRUNE="${G1K_BFILE}_ancestry_pruned_snps.txt"
    if [[ -f "$CACHED_PRUNE" ]]; then
        echo "" | tee -a "$LOG"
        echo "✅ Using cached ancestry-pruned SNP list: $CACHED_PRUNE" | tee -a "$LOG"
        FINAL_PRUNE="$CACHED_PRUNE"
        # Skip directly to applying the prune set
        DATASET="${OUT_DIR}/ld_pruned_dataset"
        echo "  Extracting $(wc -l < "$FINAL_PRUNE" | tr -d ' ') ancestry-matched independent SNPs" | tee -a "$LOG"
        "$PLINK_BIN" --bfile "$BFILE" --extract "$FINAL_PRUNE" \
            --make-bed --out "$DATASET" \
            --threads "$THREADS" --memory "$MEMORY" 2>&1 | tee -a "$LOG"
        cp "$FINAL_PRUNE" "${OUT_DIR}/pruned.prune.in" 2>/dev/null || true
        cut -f2 "${BFILE}.bim" | sort > "${OUT_DIR}/target_snps_sorted.txt"
        sort "$FINAL_PRUNE" | comm -23 "${OUT_DIR}/target_snps_sorted.txt" - > "${OUT_DIR}/pruned.prune.out" 2>/dev/null || true
        rm -f "${OUT_DIR}/target_snps_sorted.txt"
        FINAL_VARIANTS=$(wc -l < "${DATASET}.bim" | tr -d ' ')
        INIT_VARIANTS=$(wc -l < "${BFILE}.bim" | tr -d ' ')
        echo "" | tee -a "$LOG"
        echo "════════════════════════════════════════════" | tee -a "$LOG"
        echo "  ✅ Ancestry-Matched LD Pruning Complete (cached)" | tee -a "$LOG"
        echo "════════════════════════════════════════════" | tee -a "$LOG"
        echo "  Initial variants:    $INIT_VARIANTS"       | tee -a "$LOG"
        echo "  Independent SNPs:    $FINAL_VARIANTS"      | tee -a "$LOG"
        echo "  Retention:           $(python3 -c "print(f'{$FINAL_VARIANTS/$INIT_VARIANTS*100:.1f}%')")" | tee -a "$LOG"
        echo "  Dataset:             $DATASET"             | tee -a "$LOG"
        exit 0
    fi
elif [[ -n "${G1K_VCF:-}" ]] && [[ -f "${G1K_VCF:-}" ]]; then
    echo "" | tee -a "$LOG"
    echo "── Converting 1000G VCF → PLINK ──" | tee -a "$LOG"
    G1K_PREFIX="${REF_DIR}/1000G_raw"
    "$PLINK_BIN" --vcf "$G1K_VCF" --make-bed --double-id --allow-extra-chr \
        --out "$G1K_PREFIX" --threads "$THREADS" --memory "$MEMORY" 2>&1 | tee -a "$LOG"

    # Rename '.' variant IDs to unique chr:pos format
    echo "  Assigning unique variant IDs (chr:pos)..." | tee -a "$LOG"
    # Replace column 2 in BIM with chr:pos since all are '.'
    awk '{print $1"\t"$1":"$4"\t"$3"\t"$4"\t"$5"\t"$6}' "${G1K_PREFIX}.bim" > "${G1K_PREFIX}.bim.tmp"
    mv "${G1K_PREFIX}.bim.tmp" "${G1K_PREFIX}.bim"
    G1K_PREFIX="${G1K_PREFIX}"
else
    echo "  ❌ Must provide --1000g-bfile or --1000g-vcf" | tee -a "$LOG"
    exit 1
fi

N_G1K=$(wc -l < "${G1K_PREFIX}.fam" | tr -d ' ')
echo "  1000G samples: $N_G1K" | tee -a "$LOG"

# ── QC 1000G if needed ────────────────────────────────────────────────
G1K_QC="${REF_DIR}/1000G_qc"
if [[ ! -f "${G1K_QC}.bed" ]] || [[ ! -f "${G1K_QC}.bim" ]] || [[ ! -f "${G1K_QC}.fam" ]]; then
    [[ -f "${G1K_QC}.bed" ]] && rm -f "${G1K_QC}".*  # Clean partial QC files
    echo "" | tee -a "$LOG"
    echo "── QC 1000G reference ──" | tee -a "$LOG"
    "$PLINK_BIN" --bfile "$G1K_PREFIX" --geno 0.05 --maf 0.01 \
        --hwe 0.000001 midp --make-bed --out "$G1K_QC" \
        --threads "$THREADS" --memory "$MEMORY" 2>&1 | tee -a "$LOG"
fi

# ── Load population labels ────────────────────────────────────────────
echo "" | tee -a "$LOG"
echo "── Loading population labels ──" | tee -a "$LOG"

# Extract super-population labels for 1000G samples
python3 -c "
import pandas as pd

# Load panel
panel = pd.read_csv('$POP_PANEL', sep='\t', dtype=str)
print(f'  Panel loaded: {len(panel)} entries')

# Load 1000G FAM
fam = pd.read_csv('${G1K_QC}.fam', sep=r'\s+', header=None, dtype=str)
fam.columns = ['FID','IID'] + list(fam.columns[2:])
print(f'  FAM: {len(fam)} samples')

# Map populations
pop_map = dict(zip(panel.iloc[:,0], panel.iloc[:,2]))  # sample → super_pop
fam['super_pop'] = fam['IID'].map(pop_map)
print(f'  Labeled: {fam[\"super_pop\"].notna().sum()}/{len(fam)}')

# Print distribution
for pop in ['EUR','AFR','EAS','SAS','AMR']:
    n = (fam['super_pop'] == pop).sum()
    print(f'    {pop}: {n}')

# Write per-population keep files
for pop in ['EUR','AFR','EAS','SAS','AMR']:
    pop_samples = fam[fam['super_pop'] == pop]
    if len(pop_samples) > 0:
        pop_samples[['FID','IID']].to_csv(
            '${REF_DIR}/keep_{pop}.txt'.format(pop=pop),
            sep='\t', index=False, header=False
        )
        print(f'    Wrote ${REF_DIR}/keep_{pop}.txt: {len(pop_samples)} samples')
" 2>&1 | tee -a "$LOG"

# ── Per-Population LD Pruning ─────────────────────────────────────────
SUPER_POPS=("EUR" "AFR" "EAS" "SAS" "AMR")
PRUNE_FILES=()

echo "" | tee -a "$LOG"
echo "── LD Pruning per Super-Population (parallel) ──" | tee -a "$LOG"

# Run all 5 populations in parallel for ~4x speedup
# PLINK I/O is the bottleneck; each pop gets independent files
PARALLEL_MEMORY=$((MEMORY / 4))  # Split memory across parallel jobs

for POP in "${SUPER_POPS[@]}"; do
    (
        KEEP_FILE="${REF_DIR}/keep_${POP}.txt"
        if [[ ! -f "$KEEP_FILE" ]]; then
            echo "  ⚠️ $POP: no samples found, skipping"
            exit 0
        fi

        POP_PREFIX="${REF_DIR}/1000G_${POP}"

        # Extract population-specific data
        "$PLINK_BIN" --bfile "$G1K_QC" --keep "$KEEP_FILE" \
            --make-bed --out "${POP_PREFIX}" \
            --memory "$PARALLEL_MEMORY" 2>&1 | sed "s/^/  [$POP] /"

        N_POP=$(wc -l < "${POP_PREFIX}.fam" | tr -d ' ')
        if [[ "$N_POP" -lt 2 ]]; then
            echo "  ⚠️ $POP: < 2 samples ($N_POP), skipping"
            exit 0
        fi

        # LD prune this population
        "$PLINK_BIN" --bfile "${POP_PREFIX}" \
            --indep-pairwise "$WINDOW_SIZE" "$STEP_SIZE" "$R2_THRESHOLD" \
            --out "${POP_PREFIX}_ld" \
            --memory "$PARALLEL_MEMORY" 2>&1 | sed "s/^/  [$POP] /"

        N_KEEP=$(wc -l < "${POP_PREFIX}_ld.prune.in" | tr -d ' ')
        echo "  $POP: $N_KEEP independent SNPs (from $N_POP samples)" | tee -a "$LOG"
    ) &
done

# Wait for all background jobs to complete
wait

# Collect results after parallel jobs finish
for POP in "${SUPER_POPS[@]}"; do
    POP_PREFIX="${REF_DIR}/1000G_${POP}"
    PRUNE_FILE="${POP_PREFIX}_ld.prune.in"
    if [[ -f "$PRUNE_FILE" ]]; then
        PRUNE_FILES+=("$PRUNE_FILE")
    fi
done

if [[ ${#PRUNE_FILES[@]} -eq 0 ]]; then
    echo "  ❌ No populations could be LD-pruned. Check 1000G data." | tee -a "$LOG"
    exit 1
fi

# ── Compute Union/Intersection ────────────────────────────────────────
echo "" | tee -a "$LOG"
echo "── Computing ${UNION_MODE} SNP set ──" | tee -a "$LOG"

FINAL_PRUNE="${REF_DIR}/ancestry_pruned_snps.txt"

python3 -c "
import sys

pop_files = [$(printf '"%s",' "${PRUNE_FILES[@]}" | sed 's/,$//')]
print(f'  Population prune files: {len(pop_files)}')

# Load all prune.in sets
snps_per_pop = {}
for pf in pop_files:
    with open(pf) as fh:
        snps_per_pop[pf] = set(line.strip() for line in fh)

# Conservative (intersection): SNP must be independent in ALL populations
# Inclusive (union): SNP is independent in AT LEAST ONE population
mode = '$UNION_MODE'
if mode == 'conservative':
    # Intersection across all populations
    result = snps_per_pop[pop_files[0]]
    for pf in pop_files[1:]:
        result = result & snps_per_pop[pf]
    print(f'  Mode: CONSERVATIVE (intersection) — SNP independent in ALL {len(pop_files)} populations')
else:
    # Union across all populations
    result = set()
    for pf in pop_files:
        result = result | snps_per_pop[pf]
    print(f'  Mode: INCLUSIVE (union) — SNP independent in ANY of {len(pop_files)} populations')

print(f'  Final SNP set: {len(result)}')
for pop_snps, pf in zip(snps_per_pop.values(), pop_files):
    overlap = len(result & pop_snps)
    print(f'    Overlap with {pf.split(\"/\")[-1]}: {overlap}/{len(pop_snps)}')

# Write final prune.in
with open('$FINAL_PRUNE', 'w') as fh:
    for snp in sorted(result):
        fh.write(snp + '\n')
print(f'  Written: $FINAL_PRUNE ({len(result)} SNPs)')
" 2>&1 | tee -a "$LOG"

# Cache the result alongside the reference for future runs
if [[ -n "${G1K_BFILE:-}" ]] && [[ -f "${G1K_BFILE}.bed" ]]; then
    CACHED_PRUNE="${G1K_BFILE}_ancestry_pruned_snps.txt"
    cp "$FINAL_PRUNE" "$CACHED_PRUNE"
    echo "  Cached: $CACHED_PRUNE" | tee -a "$LOG"
fi

# ── Apply to Target Sample ────────────────────────────────────────────
echo "" | tee -a "$LOG"
echo "── Applying ancestry-matched LD pruning to target ──" | tee -a "$LOG"

DATASET="${OUT_DIR}/ld_pruned_dataset"

if [[ -f "$FINAL_PRUNE" ]] && [[ -s "$FINAL_PRUNE" ]]; then
    N_FINAL=$(wc -l < "$FINAL_PRUNE" | tr -d ' ')
    echo "  Extracting $N_FINAL ancestry-matched independent SNPs" | tee -a "$LOG"

    "$PLINK_BIN" --bfile "$BFILE" --extract "$FINAL_PRUNE" \
        --make-bed --out "$DATASET" \
        --threads "$THREADS" --memory "$MEMORY" 2>&1 | tee -a "$LOG"

    # Also save prune.in/out for compatibility
    cp "$FINAL_PRUNE" "${OUT_DIR}/pruned.prune.in"
    # SNPs removed = all SNPs in target minus those kept
    # Use sort+comm instead of grep -f (O(n log n) vs O(n*m))
    cut -f2 "${BFILE}.bim" | sort > "${OUT_DIR}/target_snps_sorted.txt"
    sort "$FINAL_PRUNE" | comm -23 "${OUT_DIR}/target_snps_sorted.txt" - > "${OUT_DIR}/pruned.prune.out" 2>/dev/null || true
    rm -f "${OUT_DIR}/target_snps_sorted.txt"
else
    echo "  ❌ No SNPs in final prune set. Cannot continue." | tee -a "$LOG"
    exit 1
fi

# ── Verify ────────────────────────────────────────────────────────────
FINAL_VARIANTS=$(wc -l < "${DATASET}.bim" | tr -d ' ')
FINAL_SAMPLES=$(wc -l < "${DATASET}.fam" | tr -d ' ')
INIT_VARIANTS=$(wc -l < "${BFILE}.bim" | tr -d ' ')

echo "" | tee -a "$LOG"
echo "════════════════════════════════════════════" | tee -a "$LOG"
echo "  ✅ Ancestry-Matched LD Pruning Complete"   | tee -a "$LOG"
echo "════════════════════════════════════════════" | tee -a "$LOG"
echo "  Initial variants:    $INIT_VARIANTS"       | tee -a "$LOG"
echo "  Independent SNPs:    $FINAL_VARIANTS"      | tee -a "$LOG"
echo "  Retention:           $(python3 -c "print(f'{$FINAL_VARIANTS/$INIT_VARIANTS*100:.1f}%')")" | tee -a "$LOG"
echo "  Samples:             $FINAL_SAMPLES"       | tee -a "$LOG"
echo "  Populations used:    ${#PRUNE_FILES[@]}"   | tee -a "$LOG"
echo "  Dataset:             $DATASET"             | tee -a "$LOG"
