#!/bin/bash
# ╔══════════════════════════════════════════════════════════════════════════════╗
# ║   DOWNLOAD ASHKENAZI JEWISH REFERENCE PANEL                                  ║
# ║   scripts/setup/download_ashkenazi_reference.sh                              ║
# ║                                                                            ║
# ║   Downloads the NearEastPublic dataset (Lazaridis et al. 2016) from the     ║
# ║   Reich Lab at Harvard. This dataset contains the Human Origins panel       ║
# ║   including Ashkenazi Jewish samples plus diverse European populations.     ║
# ║                                                                            ║
# ║   Source: https://reich.hms.harvard.edu/datasets                           ║
# ║   Paper:  Lazaridis et al. 2016, Nature 536, 419-424                       ║
# ║                                                                            ║
# ║   Output (~150 MB):                                                         ║
# ║     reference/human_origins/HumanOriginsPublic2068.{bed,bim,fam}            ║
# ║     reference/human_origins/population_labels.txt                           ║
# ║     reference/human_origins/european_aj_subset.{bed,bim,fam}                ║
# ╚══════════════════════════════════════════════════════════════════════════════╝

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
REF_DIR="$PIPELINE_ROOT/reference/human_origins"
DOWNLOAD_URL="https://reich.hms.harvard.edu/sites/reich.hms.harvard.edu/files/inline-files/NearEastPublic.tar.gz"
TARBALL="$REF_DIR/NearEastPublic.tar.gz"

# ═══ Step 0: Check prerequisites ═══
echo "═══ Downloading Ashkenazi Jewish Reference Panel ═══"
echo ""

mkdir -p "$REF_DIR"
cd "$REF_DIR"

# ═══ Step 1: Download ═══
if [ -f "HumanOriginsPublic2068.geno" ]; then
    echo "✅ HumanOrigins data already downloaded"
else
    echo "📥 Downloading NearEastPublic dataset (~200 MB)..."
    echo "   URL: $DOWNLOAD_URL"
    if command -v wget &>/dev/null; then
        wget -q --show-progress "$DOWNLOAD_URL" -O "$TARBALL" || {
            echo "❌ Download failed. You can download manually from:"
            echo "   $DOWNLOAD_URL"
            echo "   Extract to: $REF_DIR"
            exit 1
        }
    elif command -v curl &>/dev/null; then
        curl -L -o "$TARBALL" "$DOWNLOAD_URL" || {
            echo "❌ Download failed. See manual instructions above."
            exit 1
        }
    else
        echo "❌ Neither wget nor curl found. Install one and retry."
        exit 1
    fi
    echo "✅ Download complete"

    # ═══ Step 2: Extract ═══
    echo "📦 Extracting..."
    tar xzf "$TARBALL"
    rm -f "$TARBALL"
    echo "✅ Extraction complete"
fi

# ═══ Step 3: List available files ═══
echo ""
echo "📁 Files in reference/human_origins/:"
ls -lh "$REF_DIR" | head -20
echo ""

# ═══ Step 4: Identify Ashkenazi Jewish samples ═══
if [ -f "HumanOriginsPublic2068.ind" ]; then
    echo "🔍 Ashkenazi Jewish samples in the dataset:"
    grep -i "ashkenazi\|jew_ashk" HumanOriginsPublic2068.ind || echo "  (Searching...)"

    echo ""
    echo "🔍 European populations available:"
    cut -f3 HumanOriginsPublic2068.ind | sort | uniq -c | sort -rn | head -30
fi

echo ""
echo "═══ Download Complete ═══"
echo ""
echo "📋 Next steps (manual):"
echo "   1. Filter to European + AJ populations using PLINK:"
echo "      plink --bfile HumanOriginsPublic2068 \\"
echo "            --keep european_aj_samples.txt \\"
echo "            --make-bed --out european_aj_subset"
echo ""
echo "   2. Merge with 1000G EUR reference for subcontinental PCA:"
echo "      plink --bfile ../1000G_full/1000G_full \\"
echo "            --bmerge european_aj_subset \\"
echo "            --make-bed --out eur_aj_merged"
echo ""
echo "   3. LD-prune and run PCA:"
echo "      plink --bfile eur_aj_merged --indep-pairwise 50 5 0.2 \\"
echo "            --out eur_aj_pruned"
echo "      plink --bfile eur_aj_merged --extract eur_aj_pruned.prune.in \\"
echo "            --pca 20 --out eur_aj_pca"
echo ""
echo "   Reference: Carmi et al. 2014, Nature Communications 5, 4835"
echo "   Reference: Behar et al. 2010, Nature 466, 238-242"
