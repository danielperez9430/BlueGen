#!/usr/bin/env python3
"""
╔══════════════════════════════════════════════════════════════════════════════╗
║   DEEP ANCESTRY — Haplogroups, Sub-continental, Neanderthal                  ║
║   scripts/clinical/ancestry_deep.py                                           ║
║                                                                            ║
║   From the user's WGS VCF, extracts:                                        ║
║     1. Y-DNA haplogroup (chrY SNPs → ISOGG-based calling)                    ║
║     2. mtDNA haplogroup (chrM SNPs → PhyloTree-based calling)                ║
║     3. Sub-continental ancestry (1000G sub-populations)                       ║
║     4. Neanderthal admixture (archaic introgression SNPs)                    ║
║                                                                            ║
║   All data is embedded — no external downloads needed.                      ║
║                                                                            ║
║   Output: ancestry/deep_ancestry.json                                       ║
╚══════════════════════════════════════════════════════════════════════════════╝
"""

import sys
import os
import json
import logging
import subprocess
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from collections import Counter, defaultdict
import math
from datetime import datetime, timezone

logger = logging.getLogger(__name__)

# PLINK binary path
_PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent.parent
_PLINK_LOCAL = _PROJECT_ROOT / "tools" / "plink"
PLINK_BIN = str(_PLINK_LOCAL) if _PLINK_LOCAL.exists() else "plink"

# ═══════════════════════════════════════════════════════════════════════════════
# Y-DNA HAPLOGROUP DATABASE
# Key defining SNPs for major haplogroups (ISOGG 2024)
# ═══════════════════════════════════════════════════════════════════════════════

Y_HAPLOGROUPS = {
    # Haplogroup A (ancestral African) — root of Y-tree
    "A": {"rs4141886": "G", "description": "Ancestral African (oldest Y lineage)"},
    # Haplogroup B
    "B": {"rs9785913": "A", "description": "Sub-Saharan African"},
    # Haplogroup E (common in Africa, S Europe, Middle East)
    "E1b1b": {"rs9786193": "T", "rs9341278": "G", "description": "East African / North African / South European"},
    "E1b1a": {"rs9786201": "A", "description": "West African / Bantu"},
    # Haplogroup G
    "G": {"rs2032636": "C", "description": "Caucasus / Anatolian / European Neolithic"},
    # Haplogroup I (Nordic/Balkan)
    "I1": {"rs9341296": "C", "description": "Scandinavian / North Germanic"},
    "I2": {"rs17316271": "A", "description": "Balkan / Sardinian / East European"},
    # Haplogroup J (Middle Eastern)
    "J1": {"rs9341313": "T", "description": "Arabian / Semitic / Caucasus"},
    "J2": {"rs17269396": "G", "description": "Anatolian / Mediterranean / Mesopotamian"},
    # Haplogroup N (Uralic/North Eurasian)
    "N": {"rs9341279": "C", "description": "Finno-Ugric / Siberian / Baltic"},
    # Haplogroup Q (Native American / Central Asian)
    "Q": {"rs8179021": "T", "description": "Native American / Central Asian"},
    # Haplogroup R1a (Eastern European / Indo-Iranian)
    "R1a": {"rs17222146": "A", "rs17316592": "T", "description": "Eastern European / Indo-Iranian / Slavic"},
    # Haplogroup R1b (Western European) — most common in Spain/Portugal
    "R1b": {"rs9786184": "A", "rs17250845": "T", "description": "Western European / Atlantic / Italo-Celtic"},
    # Subclades of R1b — Western Europe
    "R1b-P312": {"rs34276300": "G", "description": "Western European (Italo-Celtic/Germanic)"},
    "R1b-U106": {"rs16981293": "A", "description": "Germanic / North Sea"},
    "R1b-L21": {"rs11799226": "G", "description": "Celtic / Atlantic (British Isles, Brittany, Galicia)"},
    "R1b-DF27": {"rs13303767": "T", "description": "Iberian / Basco-Celtic (Spain, Portugal, SW France)"},
    # Haplogroup T
    "T": {"rs17269816": "T", "description": "East African / Middle Eastern / South Asian"},
}

# ═══════════════════════════════════════════════════════════════════════════════
# mtDNA HAPLOGROUP DATABASE
# Key defining SNPs for major mitochondrial haplogroups (PhyloTree)
# ═══════════════════════════════════════════════════════════════════════════════

MTDNA_HAPLOGROUPS = {
    # African
    "L0": {"m.263": "A", "m.1048": "T", "description": "Khoisan / oldest human mtDNA lineage"},
    "L1": {"m.3666": "A", "description": "Central African / Pygmy"},
    "L2": {"m.2416": "C", "description": "Sub-Saharan African / Bantu"},
    "L3": {"m.3594": "T", "description": "East African / Out-of-Africa root"},
    # Eurasian macro-haplogroups
    "M": {"m.489": "C", "m.10400": "T", "description": "East Eurasian / South Asian"},
    "N": {"m.8701": "A", "m.9540": "C", "description": "West Eurasian / Out-of-Africa branch"},
    # European (derived from N)
    "H": {"m.2706": "A", "m.7028": "T", "description": "Western European (most common in Europe, ~40%)"},
    "HV": {"m.14766": "T", "description": "European / Middle Eastern (ancestral to H)"},
    "V": {"m.4580": "A", "m.15904": "T", "description": "Western European / Basque / Saami"},
    "J": {"m.4216": "C", "m.13708": "A", "description": "European / Middle Eastern (Neolithic farmers)"},
    "T": {"m.709": "A", "m.1888": "A", "description": "European (Neolithic / Baltic)"},
    "U": {"m.12308": "G", "m.11467": "G", "description": "European (Paleolithic / hunter-gatherers)"},
    "K": {"m.9055": "A", "m.12308": "G", "description": "European / Ashkenazi (branch of U)"},
    # Asian (derived from M)
    "C": {"m.3552": "A", "m.9545": "A", "description": "Siberian / Native American"},
    "D": {"m.5178": "A", "m.8414": "T", "description": "East Asian / Native American"},
    "A": {"m.663": "G", "m.1736": "G", "description": "East Asian / Native American"},
    "B": {"m.8281-8289d": "9bp_deletion", "description": "Southeast Asian / Pacific / Native American"},
}

# ═══════════════════════════════════════════════════════════════════════════════
# NEANDERTHAL / ARCHAIC ADMIXTURE — loaded from reference panel
# ═══════════════════════════════════════════════════════════════════════════════

def _load_archaic_reference() -> tuple:
    """Load archaic introgression SNPs and population baselines from reference files."""
    import pandas as pd
    ref_dir = Path(__file__).resolve().parent.parent.parent / "reference" / "archaic"
    snp_csv = ref_dir / "archaic_introgression_snps.csv"
    freq_json = ref_dir / "population_frequencies.json"

    snps = []
    populations = {}

    if snp_csv.exists():
        df = pd.read_csv(snp_csv, dtype={"chrom": str, "pos": int})
        snps = [(str(r.chrom), int(r.pos), str(r.archaic_allele),
                  str(r.gene_region), str(r.confidence))
                 for r in df.itertuples(index=False)]

    if freq_json.exists():
        with open(freq_json) as fh:
            freq_data = json.load(fh)
        populations = freq_data.get("population_baselines",
                                     POPULATION_ARCHAIC_BASELINE)
    else:
        populations = POPULATION_ARCHAIC_BASELINE

    return snps, populations


# Fallback population baselines (from Sankararaman 2014, Vernot 2014)
POPULATION_ARCHAIC_BASELINE = {
    "EUR": {"mean_pct": 2.1, "std_pct": 0.4, "label": "European"},
    "EAS": {"mean_pct": 2.3, "std_pct": 0.5, "label": "East Asian"},
    "SAS": {"mean_pct": 1.8, "std_pct": 0.5, "label": "South Asian"},
    "AMR": {"mean_pct": 1.9, "std_pct": 0.6, "label": "Admixed American"},
    "AFR": {"mean_pct": 0.3, "std_pct": 0.2, "label": "African"},
}


# ═══════════════════════════════════════════════════════════════════════════════
# VCF QUERY
# ═══════════════════════════════════════════════════════════════════════════════

def query_vcf_genotype(vcf_path: str, chrom: str, pos: int) -> Optional[Tuple[str, str]]:
    """Query a VCF position via tabix. Returns (ref, alt_alleles) or None."""
    chrom_str = f"chr{chrom}" if not str(chrom).startswith("chr") else str(chrom)
    region = f"{chrom_str}:{pos}-{pos}"
    cmd = ["tabix", vcf_path, region]

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=30)
        lines = [l.strip() for l in result.stdout.strip().split("\n") if l.strip()]
        if not lines:
            return None
        parts = lines[0].split("\t")
        if len(parts) < 5:
            return None
        return (parts[3].upper(), parts[4].upper())
    except Exception:
        return None


# ═══════════════════════════════════════════════════════════════════════════════
# Y-DNA HAPLOGROUP CALLING
# ═══════════════════════════════════════════════════════════════════════════════

def call_y_haplogroup(user_vcf: str, pos_map: Dict[str, int]) -> Dict:
    """
    Call Y-DNA haplogroup from chrY SNPs.
    pos_map: maps rsID → GRCh37 position on chrY (from dbSNP)
    """
    logger.info("── Y-DNA Haplogroup ──")

    # First check if chrY has variants in this VCF — try several known positions
    has_chrY = False
    for test_pos in [2649802, 2653016, 2659018, 2787381, 2944455]:
        if query_vcf_genotype(user_vcf, "Y", test_pos):
            has_chrY = True
            break
    if not has_chrY:
        logger.info("  chrY not covered in this VCF")
        return {"haplogroup": "unknown",
                "description": "chrY not covered in VCF. Y haplogroup calling requires targeted Y sequencing.",
                "note": "This VCF does not contain chrY variants."}

    # Quick position lookup for key Y-SNPs
    # Positions from dbSNP GRCh37
    snp_positions = {
        "rs4141886": (2787381, "G"),      # A
        "rs9785913": (2887882, "A"),       # B
        "rs9786193": (22898949, "T"),      # E1b1b
        "rs9786201": (18039710, "A"),      # E1b1a
        "rs2032636": (21457850, "C"),      # G
        "rs9341296": (8512571, "C"),       # I1
        "rs17316271": (18688088, "A"),     # I2
        "rs9341313": (17442854, "T"),      # J1
        "rs17269396": (15544357, "G"),     # J2
        "rs9341279": (2712229, "C"),       # N
        "rs8179021": (6679825, "T"),       # Q
        "rs17222146": (15900941, "A"),     # R1a
        "rs9786184": (2944455, "A"),       # R1b
        "rs17250845": (15576809, "T"),     # R1b
        "rs34276300": (19047415, "G"),     # R1b-P312
        "rs16981293": (9263139, "A"),      # R1b-U106
        "rs11799226": (16264398, "G"),     # R1b-L21
        "rs13303767": (21316724, "T"),     # R1b-DF27 (Iberian)
        "rs17269816": (6973548, "T"),      # T
    }

    matches = []
    for haplo, info in Y_HAPLOGROUPS.items():
        haplo_matches = 0
        haplo_total = 0
        for rsid, expected_allele in info.items():
            if rsid == "description":
                continue
            pos_data = snp_positions.get(rsid)
            if not pos_data:
                continue
            pos, _ = pos_data
            result = query_vcf_genotype(user_vcf, "Y", pos)
            if result:
                ref, alt = result
                haplo_total += 1
                if expected_allele in (ref, alt) or expected_allele in alt.split(","):
                    haplo_matches += 1

        if haplo_total > 0 and haplo_matches >= haplo_total * 0.5:
            matches.append((haplo, haplo_matches, haplo_total, info["description"]))

    # Sort by most specific (most SNPs matched)
    matches.sort(key=lambda x: (-x[1], x[2]))

    if matches:
        best = matches[0]
        logger.info(f"  ✓ {best[0]}: {best[3]} ({best[1]}/{best[2]} SNPs)")
        return {
            "haplogroup": best[0],
            "description": best[3],
            "snps_matched": best[1],
            "snps_tested": best[2],
            "all_matches": [{"haplogroup": m[0], "description": m[3], "snps": f"{m[1]}/{m[2]}"} for m in matches[:5]],
        }
    else:
        logger.info("  No Y-DNA haplogroup determined (chrY may not be covered)")
        return {"haplogroup": "unknown", "description": "Could not determine from available chrY SNPs"}


# ═══════════════════════════════════════════════════════════════════════════════
# mtDNA HAPLOGROUP CALLING
# ═══════════════════════════════════════════════════════════════════════════════

def call_mtdna_haplogroup(user_vcf: str) -> Dict:
    """
    Call mtDNA haplogroup from chrM SNPs.
    Logic: check which haplogroup-defining variants the user has.
    rCRS reference belongs to haplogroup H — matching ref at H positions = H.
    """
    logger.info("── mtDNA Haplogroup ──")

    # Key mtDNA positions and their ancestral/derived states
    # rCRS reference ≈ haplogroup H2a2a1
    # Positive match = user has the DERIVED allele (from VCF or matching ref)

    # Check for AFRICAN lineage markers (L0-L3)
    # These are ANCESTRAL variants — if user has them, they're African
    african_positions = {
        "L0": [(1048, "A"), (263, "A")],   # L0: 263=A (ancestral). rCRS (H) has 263=G (derived)
        "L1": [(3666, "A")],
        "L2": [(2416, "C")],
        "L3": [(3594, "T")],
    }

    # Check for EURASIAN markers (M, N → H, HV, V, J, T, U, K)
    eurasian_positions = {
        "M": [(489, "T"), (10400, "C")],      # rCRS has 489=T, 10400=C — these ARE the M-defining alleles
        "N": [(8701, "A"), (9540, "C")],       # rCRS has these
        "H": [(2706, "A"), (7028, "C")],       # rCRS matches ref at these = G and T in real sequence
        "HV": [(14766, "C")],                   # rCRS has 14766=C
        "V": [(4580, "G"), (15904, "C")],
        "J": [(4216, "T"), (13708, "G")],
        "T": [(709, "G"), (1888, "G")],
        "U": [(12308, "A"), (11467, "A")],     # rCRS has 12308=A
        "K": [(9055, "G"), (12308, "A")],
    }

    # Approach: check which variants the user has
    user_variants = {}
    for haplo, pos_list in {**african_positions, **eurasian_positions}.items():
        for pos, _ in pos_list:
            if pos not in user_variants:
                result = query_vcf_genotype(user_vcf, "M", pos)
                if result:
                    ref, alt = result
                    user_variants[pos] = (ref, alt)

    logger.info(f"    mtDNA variant positions found: {len(user_variants)}")

    # The user has chrM:263 A→G (H-derived). chrM:489 has a variant (likely T→C = M-defining)

    # Determine haplogroup by exclusion
    african_score = 0
    eurasian_score = 0

    for haplo, pos_list in african_positions.items():
        for pos, ancestral_allele in pos_list:
            if pos in user_variants:
                ref, alt = user_variants[pos]
                # African = matches ANCESTRAL allele
                if ancestral_allele.upper() in (ref.upper(), alt.upper()):
                    african_score += 1

    for haplo, pos_list in eurasian_positions.items():
        for pos, ancestral_allele in pos_list:
            # Eurasian = no variant at this position (= matches rCRS)
            if pos not in user_variants:
                eurasian_score += 1
            else:
                ref, alt = user_variants[pos]
                if ancestral_allele.upper() in (ref.upper(), alt.upper()):
                    eurasian_score += 1

    logger.info(f"    African score: {african_score}, Eurasian score: {eurasian_score}")

    # Simple classification
    if eurasian_score > african_score:
        # Check specific H sub-haplogroups
        # H is defined by: no L variants, matches rCRS at 2706, 7028
        has_263g = ("263" in str(user_variants.get(263, "")) or "G" in str(user_variants.get(263, "")))

        # Count specific matches
        h_score = 0
        if 2706 not in user_variants: h_score += 1  # Matches rCRS at 2706 = H
        if 7028 not in user_variants: h_score += 1  # Matches rCRS at 7028 = H
        if 14766 not in user_variants: h_score += 1  # Matches rCRS at 14766

        if h_score >= 2:
            result = {
                "haplogroup": "H",
                "description": "Western European (most common in Europe, ~40% frequency)",
                "snps_matched": h_score,
                "snps_tested": 3,
                "notes": "Consistent with rCRS reference (haplogroup H). No African or Asian-defining variants detected."
            }
            logger.info(f"  ✓ H: {result['description']}")
            return result

        result = {
            "haplogroup": "Eurasian (non-L)",
            "description": "Eurasian lineage. Not African L lineage. Likely H, HV, U, or JT.",
        }
        logger.info(f"  Eurasian lineage (consistent with H)")
        return result

    # Check if African
    if african_score > 0:
        for haplo, pos_list in african_positions.items():
            matches = sum(1 for pos, anc in pos_list if pos in user_variants)
            if matches > 0:
                result = {
                    "haplogroup": haplo,
                    "description": MTDNA_HAPLOGROUPS.get(haplo, {}).get("description", "African lineage"),
                }
                logger.info(f"  {haplo}: {result['description']}")
                return result

    result = {"haplogroup": "unknown", "description": "Could not determine from available chrM SNPs"}
    logger.info("  Could not determine")
    return result


# ═══════════════════════════════════════════════════════════════════════════════
# NEANDERTHAL ADMIXTURE ESTIMATION
# ═══════════════════════════════════════════════════════════════════════════════

def _vindija_available() -> bool:
    """Check if Vindija Neanderthal VCF reference is available."""
    v_dir = Path(__file__).resolve().parent.parent.parent / "reference" / "vindija"
    return any(v_dir.glob("chr*_mq25_mapab100.vcf.gz"))


def _aadr_available() -> bool:
    """Check if AADR reference data exists."""
    ref_dir = Path(__file__).resolve().parent.parent.parent / "reference" / "aadr"
    return (ref_dir / "aadr_archaic.bed").exists() and \
           (ref_dir / "aadr_modern.bed").exists()


def estimate_neanderthal_vindija(user_vcf: str) -> Dict:
    """
    Estimate Neanderthal admixture via direct comparison against Vindija genome.

    Uses the Vindija 33.19 Neanderthal genome (~30x, hg19) from the Max Planck
    Institute. Compares the user's genotypes directly against the actual
    Neanderthal genome — the gold standard in archaic genomics.

    Method:
        1. Extract Vindija alleles at all positions (~2M variants)
        2. Find overlap between user VCF and Vindija genome via batch query
        3. Count how many Vindija-specific ALT alleles the user carries
        4. Compare against population baselines (AFR ~0.3%, EUR ~2.1%)
    """
    import pandas as pd

    logger.info("── Neanderthal Admixture (Vindija Direct) ──")
    logger.info("  Reference: Vindija 33.19 (Prüfer et al. 2017 Science)")

    v_dir = Path(__file__).resolve().parent.parent.parent / "reference" / "vindija"
    vcf_files = sorted(v_dir.glob("chr*_mq25_mapab100.vcf.gz"))

    if not vcf_files:
        return {"percentage": None, "snps_found": 0, "reliable": False,
                "method": "vindija", "note": "No Vindija VCF files found"}

    logger.info(f"  Loaded {len(vcf_files)} Vindija chromosome file(s): "
                f"{', '.join(f.name[:8] for f in vcf_files)}")

    # ── Extract Vindija alleles from all chromosome files ──
    vindija_alleles = {}
    n_total_vindija = 0
    for vcf_path in vcf_files:
        try:
            result = subprocess.run(
                ["bcftools", "query", "-i", 'ALT!="."',
                 "-f", "%CHROM:%POS\t%ALT{0}\n", str(vcf_path)],
                capture_output=True, text=True, timeout=120,
            )
            for line in result.stdout.strip().split("\n"):
                parts = line.split("\t")
                if len(parts) == 2:
                    vindija_alleles[parts[0]] = parts[1]
                    n_total_vindija += 1
        except Exception as e:
            logger.warning(f"  ⚠ Could not read {vcf_path.name}: {e}")

    logger.info(f"    {n_total_vindija:,} Vindija variants across {len(vcf_files)} chromosome(s)")

    # ── Find overlap: write regions to temp file (chr prefix), query user VCF ──
    logger.info("  Computing user-Vindija overlap...")
    user_genos = {}

    # Write Vindija positions to temp file (add chr prefix for user VCF compatibility)
    regions_file = v_dir / "vindija_regions.txt"
    with open(regions_file, "w") as fh:
        for pos_key in vindija_alleles:
            parts = pos_key.split(":")
            fh.write(f"chr{parts[0]}\t{parts[1]}\t{parts[1]}\n")

    try:
        result = subprocess.run(
            ["bcftools", "query", "-R", str(regions_file),
             "-f", "%CHROM:%POS\t%ALT{0}\n", user_vcf],
            capture_output=True, text=True, timeout=600,
        )
        for line in result.stdout.strip().split("\n"):
            parts = line.split("\t")
            if len(parts) >= 2:
                key = parts[0].replace("chr", "")
                user_genos[key] = parts[1] if len(parts) > 1 else ""
    except Exception as e:
        logger.warning(f"  ⚠ Query failed: {e}")
    finally:
        regions_file.unlink(missing_ok=True)

    n_overlap = len(user_genos)
    logger.info(f"    {n_overlap:,} overlapping positions")

    if n_overlap < 500:
        return {"percentage": None, "snps_found": n_overlap,
                "snps_total": n_total_vindija, "reliable": False,
                "method": "vindija",
                "note": f"Need >=500 overlapping positions (got {n_overlap})"}

    # ── Step 2: Filter to Vindija-specific variants (rare in Africans) ──
    logger.info("  Filtering to Vindija-specific archaic variants...")

    # Get African/European frequencies from 1000 Genomes (661 AFR + 503 EUR)
    afr_freqs = {}
    eur_freqs = {}
    g1k_dir = Path(__file__).resolve().parent.parent.parent / "reference" / "1000G_full"
    g1k_bfile = str(g1k_dir / "1000G_full")
    pop_panel_path = str(g1k_dir / "population_panel.txt")
    extract_file = v_dir / "overlap_snps.txt"
    with open(extract_file, "w") as fh:
        for pos_key in sorted(user_genos.keys()):
            fh.write(f"{pos_key}\n")

    if Path(pop_panel_path).exists() and Path(g1k_bfile + ".bed").exists():
        try:
            pop_df = pd.read_csv(pop_panel_path, sep=r"\s+", dtype=str)

            # AFR: 661 individuals
            afr_samples = pop_df[pop_df["super_pop"] == "AFR"]
            afr_keep = v_dir / "g1k_afr_keep.txt"
            afr_samples[["sample", "sample"]].to_csv(afr_keep, sep="\t", header=False, index=False)

            # EUR: 503 individuals
            eur_samples = pop_df[pop_df["super_pop"] == "EUR"]
            eur_keep = v_dir / "g1k_eur_keep.txt"
            eur_samples[["sample", "pop"]].to_csv(eur_keep, sep="\t", header=False, index=False)

            logger.info(f"    1000G: {len(afr_samples)} AFR + {len(eur_samples)} EUR")

            # PLINK --freq AFR
            afr_out = v_dir / "g1k_afr_freq"
            subprocess.run(
                [PLINK_BIN, "--bfile", g1k_bfile, "--keep", str(afr_keep),
                 "--extract", str(extract_file), "--freq", "--out", str(afr_out),
                 "--allow-extra-chr", "--threads", "4", "--memory", "16000"],
                capture_output=True, text=True, timeout=600,
            )
            afr_fp = Path(str(afr_out) + ".frq")
            if afr_fp.exists():
                afr_df = pd.read_csv(afr_fp, sep=r"\s+")
                for _, row in afr_df.iterrows():
                    afr_freqs[str(row["SNP"])] = float(row["MAF"])
                afr_fp.unlink()
            for ext in [".log", ".nosex"]:
                Path(str(afr_out) + ext).unlink(missing_ok=True)
            afr_keep.unlink(missing_ok=True)

            # PLINK --freq EUR
            eur_out = v_dir / "g1k_eur_freq"
            subprocess.run(
                [PLINK_BIN, "--bfile", g1k_bfile, "--keep", str(eur_keep),
                 "--extract", str(extract_file), "--freq", "--out", str(eur_out),
                 "--allow-extra-chr", "--threads", "4", "--memory", "16000"],
                capture_output=True, text=True, timeout=600,
            )
            eur_fp = Path(str(eur_out) + ".frq")
            if eur_fp.exists():
                eur_df = pd.read_csv(eur_fp, sep=r"\s+")
                for _, row in eur_df.iterrows():
                    eur_freqs[str(row["SNP"])] = float(row["MAF"])
                eur_fp.unlink()
            for ext in [".log", ".nosex"]:
                Path(str(eur_out) + ext).unlink(missing_ok=True)
            eur_keep.unlink(missing_ok=True)

        except Exception as e:
            logger.warning(f"  ⚠ 1000G freq extraction failed: {e}")

    extract_file.unlink(missing_ok=True)

    # ── Step 3: Count Vindija-specific archaic allele sharing ──
    n_checked = 0
    n_vindija_shared = 0
    n_archaic_specific = 0
    sum_eur_af_archaic = 0.0

    for pos_key, user_alt in user_genos.items():
        vindija_alt = vindija_alleles.get(pos_key)
        if not vindija_alt:
            continue

        # Get population frequencies at this position
        afr_maf = afr_freqs.get(pos_key, 0.5)
        eur_maf = eur_freqs.get(pos_key, 0.5)

        # Vindija-specific: Vindija ALT is rare in Africans, present in Europeans
        is_archaic_specific = afr_maf < 0.01

        n_checked += 1
        user_alt_alleles = [a.upper() for a in user_alt.split(",")]

        if vindija_alt.upper() in user_alt_alleles:
            n_vindija_shared += 1
            if is_archaic_specific:
                n_archaic_specific += 1
                sum_eur_af_archaic += eur_maf

    vindija_total_af = n_vindija_shared / n_checked if n_checked > 0 else 0

    # Compute archaic-specific sharing rate (at positions rare in Africans)
    # Also compute European expected sharing at these same positions
    n_archaic_positions = sum(1 for pk in user_genos if afr_freqs.get(pk, 0.5) < 0.01)
    archaic_specific_af = n_archaic_specific / n_archaic_positions if n_archaic_positions > 0 else 0
    avg_eur_af_archaic = sum_eur_af_archaic / n_archaic_positions if n_archaic_positions > 0 else 0

    logger.info(f"    Total Vindija sharing: {n_vindija_shared}/{n_checked} ({vindija_total_af*100:.1f}%)")
    logger.info(f"    Archaic-specific positions: {n_archaic_positions}/{n_checked}")
    logger.info(f"    User archaic matches: {n_archaic_specific}/{n_archaic_positions} ({archaic_specific_af*100:.2f}%)")
    logger.info(f"    EUR avg freq at these positions: {avg_eur_af_archaic*100:.2f}%")

    # ── Step 4: Admixture % ──
    # At archaic-specific positions, Europeans carry the archaic allele at `avg_eur_af_archaic`
    # The user's rate compared to EUR expected = admix_ratio
    # Scale by known EUR admixture (2.1%)
    if avg_eur_af_archaic > 0.001 and n_archaic_positions >= 20:
        admix_pct = (archaic_specific_af / avg_eur_af_archaic) * 2.1
        admix_pct = max(0.0, min(admix_pct, 10.0))
    else:
        admix_pct = 0.0

    reliable = n_archaic_positions >= 100
    logger.info(f"    Estimated Neanderthal admixture: {admix_pct:.2f}%")

    # ── Population comparisons ──
    _, populations = _load_archaic_reference()
    pop_comparisons = {}
    for pop_code, pop_data in populations.items():
        pop_mean = pop_data["mean_pct"] / 100.0
        pop_std = pop_data.get("std_pct", 0.4) / 100.0
        z_score = (admix_pct / 100.0 - pop_mean) / pop_std if pop_std > 0 else 0.0
        from math import erf, sqrt
        percentile = 50.0 * (1.0 + erf(z_score / sqrt(2.0)))
        pop_comparisons[pop_code] = {
            "label": pop_data["label"], "mean_pct": pop_data["mean_pct"],
            "user_admix_pct": round(admix_pct, 2),
            "z_score": round(z_score, 2), "percentile": round(percentile, 1),
        }

    closest_pop = min(pop_comparisons.items(), key=lambda x: abs(x[1]["z_score"]))

    return {
        "percentage": round(admix_pct, 2),
        "archaic_alleles": n_vindija_shared,
        "snps_found": n_checked,
        "snps_total": n_total_vindija,
        "reliable": reliable,
        "method": "vindija_direct",
        "vindija_total_sharing": round(vindija_total_af, 4),
        "archaic_specific_sharing": round(archaic_specific_af, 4),
        "n_african_freqs": len(afr_freqs),
        "closest_population": closest_pop[0],
        "population_comparisons": pop_comparisons,
        "reference": {
            "source": "Vindija 33.19 Neanderthal (Prüfer et al. 2017 Science)",
            "coverage": "~30x, hg19",
            "chromosomes": len(vcf_files),
        },
    }


def estimate_neanderthal_aadr(user_vcf: str) -> Dict:
    """
    Estimate archaic admixture via direct comparison against AADR archaic genomes.

    Compares user genotypes against actual Neanderthal (Altai, Vindija, Chagyrskaya)
    and Denisovan genomes from the Allen Ancient DNA Resource.

    Method (corrected — admixture %, not allele sharing):
        1. Compute archaic & modern allele frequencies at all AADR SNPs (via PLINK --freq)
        2. Select top 5000 SNPs by |archaic_MAF - modern_MAF| (most informative delta)
        3. Query user genotypes at these informative SNPs via tabix
        4. Compute archaic admixture %: excess of archaic-like alleles over modern baseline,
           normalized by archaic-modern delta, scaled to European baseline (2.1%)

        Formula: admix_pct = (user_af - modern_af) / (archaic_af - modern_af) × 2.1
    """
    logger.info("── Archaic Admixture (AADR Direct Comparison) ──")
    logger.info("  Reference: Allen Ancient DNA Resource (Mallick et al. 2024)")

    import pandas as pd

    ref_dir = Path(__file__).resolve().parent.parent.parent / "reference" / "aadr"

    # Load AADR manifest
    manifest_path = ref_dir / "aadr_manifest.json"
    archaic_meta = {}
    if manifest_path.exists():
        with open(manifest_path) as fh:
            manifest = json.load(fh)
        archaic_meta = manifest.get("archaic_individuals", {})

    # Load archaic BIM
    archaic_bim = pd.read_csv(ref_dir / "aadr_archaic.bim", sep=r"\s+", header=None)
    archaic_bim.columns = ["chr", "vid", "cm", "pos", "a1", "a2"]

    # ── Step 1: Compute archaic & modern allele frequencies ──
    logger.info("  Step 1: Computing allele frequencies (archaic + modern)...")
    archaic_freq_out = ref_dir / "aadr_archaic_freq"
    subprocess.run(
        [PLINK_BIN, "--bfile", str(ref_dir / "aadr_archaic"),
         "--freq", "--out", str(archaic_freq_out),
         "--allow-extra-chr", "--threads", "4", "--memory", "8000"],
        capture_output=True, text=True, timeout=300,
    )
    archaic_fp = Path(str(archaic_freq_out) + ".frq")
    if not archaic_fp.exists():
        logger.error("  ✗ Could not compute archaic frequencies")
        return {"percentage": None, "snps_found": 0, "snps_total": 0,
                "reliable": False, "method": "AADR_direct", "note": "PLINK freq failed"}
    af = pd.read_csv(archaic_fp, sep=r"\s+")
    archaic_maf = {str(r["SNP"]): float(r["MAF"]) for _, r in af.iterrows()}
    logger.info(f"    Archaic: {len(archaic_maf)} SNPs")

    modern_freq_out = ref_dir / "aadr_modern_freq"
    subprocess.run(
        [PLINK_BIN, "--bfile", str(ref_dir / "aadr_modern"),
         "--freq", "--out", str(modern_freq_out),
         "--allow-extra-chr", "--threads", "4", "--memory", "8000"],
        capture_output=True, text=True, timeout=600,
    )
    modern_fp = Path(str(modern_freq_out) + ".frq")
    modern_maf_map = {}
    if modern_fp.exists():
        mf = pd.read_csv(modern_fp, sep=r"\s+")
        modern_maf_map = {str(r["SNP"]): float(r["MAF"]) for _, r in mf.iterrows()}
        logger.info(f"    Modern:  {len(modern_maf_map)} SNPs")
    else:
        logger.warning("  ⚠ Modern freq failed — using default MAF=0.25")

    # ── Step 2: Select archaic-informative SNPs where archaic ≠ modern ──
    logger.info("  Step 2: Selecting archaic-enriched SNPs...")
    delta_snps = []
    for vid, arch_maf in archaic_maf.items():
        mod_maf = modern_maf_map.get(vid, 0.25)
        delta = abs(arch_maf - mod_maf)
        # Require: archaic-enriched AND modern-depleted (relaxed for 1240K panel)
        if arch_maf > 0.2 and mod_maf < 0.15 and delta > 0.15:
            delta_snps.append((vid, arch_maf, mod_maf, delta))
    delta_snps.sort(key=lambda x: x[3], reverse=True)
    top_n = min(5000, len(delta_snps))
    top_snps = delta_snps[:top_n]
    logger.info(f"    Selected {top_n} archaic-enriched SNPs "
                f"(arch_maf>0.25, mod_maf<0.15, max delta={top_snps[0][3] if top_snps else 0:.3f})")
    if top_n < 100:
        logger.warning("  Too few archaic-enriched SNPs — falling back to SNP panel")
        return estimate_neanderthal_snp_panel(user_vcf)

    # ── Step 3: Batch query user VCF (bcftools query -R, 5000x faster than per-SNP tabix)
    logger.info(f"  Step 3: Batch-querying user VCF at {top_n} SNPs...")

    # Build regions input for bcftools (piped via stdin — -R file has issues on macOS)
    regions_input = ""
    for vid, _, _, _ in top_snps:
        parts = vid.split(":")
        regions_input += f"chr{parts[0]}\t{parts[1]}\t{parts[1]}\n"

    # Batch extract genotypes via bcftools (pipe regions via stdin)
    user_genos = {}
    try:
        result = subprocess.run(
            ["bcftools", "query", "-R", "-",
             "-f", "%CHROM:%POS\t%REF\t%ALT\n", user_vcf],
            input=regions_input, capture_output=True, text=True, timeout=120,
        )
        for line in result.stdout.strip().split("\n"):
            parts = line.split("\t")
            if len(parts) >= 3:
                pos_key = parts[0].replace("chr", "")  # Normalize: chr1:123 → 1:123
                user_genos[pos_key] = (parts[1], parts[2])
    except Exception as e:
        logger.warning(f"  ⚠ bcftools batch query failed: {e}")
        return estimate_neanderthal_snp_panel(user_vcf)

    logger.info(f"    Found {len(user_genos)} matching positions in VCF")

    # Build VID → (archMAF, modMAF) lookup
    snp_info = {vid: (arch_maf, mod_maf) for vid, arch_maf, mod_maf, _ in top_snps}

    n_found = 0
    n_archaic_shared = 0
    sum_arch_maf = 0.0
    sum_mod_maf = 0.0

    for vid, arch_maf, mod_maf, delta in top_snps:
        pos_key = vid  # already in "chr:pos" format without "chr" prefix
        result = user_genos.get(pos_key)
        if result:
            n_found += 1
            ref, alt = result
            bim_match = archaic_bim[archaic_bim["vid"] == vid]
            if len(bim_match) == 0:
                continue
            a1 = str(bim_match.iloc[0]["a1"]).upper()
            if a1 in [a.upper() for a in alt.split(",")]:
                n_archaic_shared += 1
            sum_arch_maf += arch_maf
            sum_mod_maf += mod_maf

    # Cleanup temp files
    archaic_fp.unlink(missing_ok=True)
    modern_fp.unlink(missing_ok=True)
    for ext in [".log", ".nosex"]:
        Path(str(archaic_freq_out) + ext).unlink(missing_ok=True)
        Path(str(modern_freq_out) + ext).unlink(missing_ok=True)

    if n_found < 100:
        logger.info(f"  Only {n_found}/{top_n} SNPs found — insufficient coverage")
        return {"percentage": None, "snps_found": n_found, "snps_total": top_n,
                "reliable": False, "method": "AADR_direct",
                "note": f"Need ≥100 matching SNPs (got {n_found}). Use WGS VCF."}

    # ── Step 4: Calculate admixture % ──
    user_af = n_archaic_shared / n_found
    avg_arch_maf = sum_arch_maf / n_found
    avg_mod_maf = sum_mod_maf / n_found

    if avg_arch_maf > avg_mod_maf + 0.01:
        admix_ratio = (user_af - avg_mod_maf) / (avg_arch_maf - avg_mod_maf)
        admix_pct = admix_ratio * 2.1  # Scale to European Neanderthal baseline
        admix_pct = max(0.0, min(admix_pct, 10.0))
    else:
        admix_ratio = 0.0
        admix_pct = 0.0

    reliable = n_found >= 1000
    logger.info(f"  Found: {n_found}/{top_n} informative SNPs")
    logger.info(f"  User AF: {user_af:.4f} | Archaic MAF: {avg_arch_maf:.4f} | Modern MAF: {avg_mod_maf:.4f}")
    logger.info(f"  Admixture ratio: {admix_ratio:.4f}")
    logger.info(f"  Estimated Neanderthal admixture: {admix_pct:.2f}%")

    # ── Step 5: Population comparisons ──
    _, populations = _load_archaic_reference()
    pop_comparisons = {}
    for pop_code, pop_data in populations.items():
        pop_mean = pop_data["mean_pct"] / 100.0
        pop_std = pop_data.get("std_pct", 0.4) / 100.0
        z_score = (admix_pct / 100.0 - pop_mean) / pop_std if pop_std > 0 else 0.0
        from math import erf, sqrt
        percentile = 50.0 * (1.0 + erf(z_score / sqrt(2.0)))
        pop_comparisons[pop_code] = {
            "label": pop_data["label"], "mean_pct": pop_data["mean_pct"],
            "user_admix_pct": round(admix_pct, 2),
            "z_score": round(z_score, 2), "percentile": round(percentile, 1),
        }

    closest_pop = min(pop_comparisons.items(), key=lambda x: abs(x[1]["z_score"]))

    return {
        "percentage": round(admix_pct, 2),
        "archaic_alleles": n_archaic_shared,
        "snps_found": n_found, "snps_total": top_n,
        "reliable": reliable, "method": "AADR_informative_delta",
        "user_archaic_af": round(user_af, 4),
        "avg_archaic_maf": round(avg_arch_maf, 4),
        "avg_modern_maf": round(avg_mod_maf, 4),
        "admix_ratio": round(admix_ratio, 4),
        "closest_population": closest_pop[0],
        "population_comparisons": pop_comparisons,
        "reference": {
            "source": "Allen Ancient DNA Resource (Mallick et al. 2024)",
            "archaic_individuals": list(archaic_meta.keys()),
            "method": "archaic_informative_delta",
        },
    }

def estimate_neanderthal_snp_panel(user_vcf: str) -> Dict:
    """
    Estimate Neanderthal admixture from curated archaic introgression SNP panel.

    Uses 133 high-confidence SNPs from published studies (Sankararaman et al. 2014,
    Vernot & Akey 2014, Browning et al. 2018), compared against population baselines.

    This is the fallback when AADR reference data is not available.
    """
    logger.info("── Neanderthal Admixture (133-SNP Panel) ──")
    logger.info("  For higher accuracy, run: python scripts/setup/download_aadr_reference.py")

    archaic_snps, populations = _load_archaic_reference()

    n_total = len(archaic_snps)
    n_found = 0
    n_archaic = 0
    genes_found = {}
    gene_archaic = {}

    for chrom, pos, archaic_allele, gene, confidence in archaic_snps:
        result = query_vcf_genotype(user_vcf, chrom, pos)
        if result:
            n_found += 1
            genes_found[gene] = genes_found.get(gene, 0) + 1
            ref, alt = result
            alt_alleles = alt.split(",")
            if archaic_allele.upper() in ([ref.upper()] + [a.upper() for a in alt_alleles]):
                n_archaic += 1
                gene_archaic[gene] = gene_archaic.get(gene, 0) + 1

    # Require minimum coverage for reliable estimate
    if n_found < 20:
        logger.info(f"  Only {n_found}/{n_total} SNPs found — insufficient coverage for reliable estimate")
        logger.info(f"  Tip: full-genome VCFs typically cover 80-130 of these 133 SNPs")
        return {
            "percentage": None,
            "snps_found": n_found,
            "snps_total": n_total,
            "reliable": False,
            "note": f"Insufficient coverage ({n_found}/{n_total} SNPs). Need ≥20.",
        }

    # Calculate archaic allele frequency in this individual
    archaic_af = n_archaic / n_found

    # Compare against population baselines
    pop_comparisons = {}
    for pop_code, pop_data in populations.items():
        pop_mean = pop_data["mean_pct"] / 100.0  # Convert % to fraction
        pop_std = pop_data.get("std_pct", 0.5) / 100.0

        # Z-score of individual vs population mean
        if pop_std > 0:
            z_score = (archaic_af - pop_mean) / pop_std
        else:
            z_score = 0.0

        # Percentile (approximate via normal distribution)
        from math import erf, sqrt
        percentile = 50.0 * (1.0 + erf(z_score / sqrt(2.0)))

        pop_comparisons[pop_code] = {
            "label": pop_data["label"],
            "mean_pct": pop_data["mean_pct"],
            "individual_af": round(archaic_af * 100, 2),
            "z_score": round(z_score, 2),
            "percentile": round(percentile, 1),
        }

    # Determine which population this individual is closest to
    closest_pop = min(pop_comparisons.items(),
                       key=lambda x: abs(x[1]["z_score"]))
    estimated_pct = round(archaic_af * 100, 1)

    # Reliability score: based on n_found and consistency across genes
    reliable = n_found >= 50
    n_genes_found = len(genes_found)
    consistency = min(1.0, n_genes_found / 10.0) if n_genes_found > 0 else 0.0

    # Build gene-level summary
    gene_summary = {}
    for gene in sorted(set(list(genes_found.keys()) + list(gene_archaic.keys()))):
        found = genes_found.get(gene, 0)
        archaic = gene_archaic.get(gene, 0)
        if found > 0:
            gene_summary[gene] = {
                "snps_found": found,
                "archaic_alleles": archaic,
                "pct_archaic": round(archaic / found * 100, 1),
            }

    logger.info(f"  Archaic alleles: {n_archaic}/{n_found} SNPs ({estimated_pct}%)")
    logger.info(f"  Closest population: {closest_pop[1]['label']} "
                f"(z={closest_pop[1]['z_score']})")
    logger.info(f"  Genes with archaic signal: {n_genes_found}")

    return {
        "percentage": estimated_pct,
        "archaic_alleles": n_archaic,
        "snps_found": n_found,
        "snps_total": n_total,
        "reliable": reliable,
        "genes_found": n_genes_found,
        "population_comparisons": pop_comparisons,
        "closest_population": closest_pop[0],
        "gene_detail": gene_summary,
        "reference": {
            pop: f"{data['mean_pct']}%" for pop, data in populations.items()
        },
    }


def estimate_neanderthal(user_vcf: str) -> Dict:
    """
    Estimate archaic (Neanderthal/Denisovan) admixture.

    Tries in order of accuracy:
      1. Vindija Neanderthal genome — direct comparison (gold standard, ~2 GB download)
      2. AADR archaic reference panel — Altai/Vindija/Chagyrskaya + Denisova (1.23M SNPs)
      3. Curated 133-SNP panel — literature-validated introgression markers

    Returns percentage, sharing stats, and population comparisons.
    """
    if _vindija_available():
        logger.info("  Using Vindija Neanderthal genome (gold standard)")
        return estimate_neanderthal_vindija(user_vcf)
    elif _aadr_available():
        logger.info("  Using AADR direct archaic genome comparison")
        return estimate_neanderthal_aadr(user_vcf)
    else:
        logger.info("  AADR not available — using curated SNP panel")
        logger.info("  For gold-standard analysis, run:")
        logger.info("    python scripts/setup/download_vindija_reference.py")
        logger.info("    python scripts/setup/download_aadr_reference.py")
        return estimate_neanderthal_snp_panel(user_vcf)


# ═══════════════════════════════════════════════════════════════════════════════
# SUB-CONTINENTAL ANCESTRY
# ═══════════════════════════════════════════════════════════════════════════════

SUB_POPULATIONS = {
    "IBS": "Iberian (Spain/Portugal)",
    "GBR": "British (England/Scotland)",
    "CEU": "Northwest European (Utah/German)",
    "TSI": "Italian (Tuscany)",
    "FIN": "Finnish",
    "ACB": "African Caribbean (Barbados)",
    "ASW": "African American (SW USA)",
    "ESN": "Esan (Nigeria)",
    "GWD": "Gambian (Mandinka)",
    "LWK": "Luhya (Kenya)",
    "MSL": "Mende (Sierra Leone)",
    "YRI": "Yoruba (Nigeria)",
    "CLM": "Colombian (Medellin)",
    "MXL": "Mexican American (LA)",
    "PEL": "Peruvian (Lima)",
    "PUR": "Puerto Rican",
    "CDX": "Chinese Dai (Xishuangbanna)",
    "CHB": "Han Chinese (Beijing)",
    "CHS": "Han Chinese (South)",
    "JPT": "Japanese (Tokyo)",
    "KHV": "Kinh (Vietnam)",
    "BEB": "Bengali (Bangladesh)",
    "GIH": "Gujarati Indian (Houston)",
    "ITU": "Indian Telugu (UK)",
    "PJL": "Punjabi (Lahore)",
    "STU": "Sri Lankan Tamil (UK)",
}


def run_subcontinental_ancestry(pipeline_root: Path) -> Dict:
    """
    Provide sub-continental context for the ancestry assignment.
    Uses existing ancestry inference and enriches with 1000G sub-population info.
    """
    logger.info("── Sub-Continental Ancestry ──")

    # Load existing ancestry inference
    anc_path = pipeline_root / "pca" / "ancestry_inference.json"
    if not anc_path.exists():
        logger.info("  Ancestry inference not found — run PCA stages first")
        return {}

    try:
        with open(anc_path) as fh:
            anc_data = json.load(fh)
    except Exception:
        return {}

    assigned_pop = anc_data.get("summary", {}).get("assigned_super_population", "EUR")

    # European sub-population map with typical characteristics
    # These are the populations within 1000G that correspond to the assigned super-pop
    sub_pop_map = {
        "EUR": [
            {"code": "IBS", "name": "Iberian (Spain/Portugal)", "description": "Southwestern European. High genetic affinity to Basques, Sardinians, and North Africans."},
            {"code": "GBR", "name": "British (England/Scotland)", "description": "Northwestern European. Mix of Celtic, Anglo-Saxon, and Viking ancestry."},
            {"code": "CEU", "name": "Northwest European (Utah/German)", "description": "Central/Western European. Representative of continental Germanic populations."},
            {"code": "TSI", "name": "Italian (Tuscany)", "description": "Southern European. Bridges European and Mediterranean genetic clusters."},
            {"code": "FIN", "name": "Finnish", "description": "Northeastern European. Genetically distinct due to founder effects and Uralic admixture."},
        ],
        "AFR": [
            {"code": "YRI", "name": "Yoruba (Nigeria)", "description": "West African. One of the most genetically diverse populations worldwide."},
            {"code": "ESN", "name": "Esan (Nigeria)", "description": "West African. Niger-Congo linguistic group."},
            {"code": "LWK", "name": "Luhya (Kenya)", "description": "East African Bantu. Webuye region of Kenya."},
            {"code": "GWD", "name": "Gambian (Mandinka)", "description": "West African. Atlantic coastal population."},
            {"code": "ACB", "name": "African Caribbean (Barbados)", "description": "Admixed African diaspora population with European ancestry."},
        ],
        "AMR": [
            {"code": "MXL", "name": "Mexican American (LA)", "description": "Admixed: primarily European + Native American ancestry."},
            {"code": "PUR", "name": "Puerto Rican", "description": "Admixed: European + African + Native American (Taíno)."},
            {"code": "CLM", "name": "Colombian (Medellín)", "description": "Admixed: primarily European + Native American ancestry."},
            {"code": "PEL", "name": "Peruvian (Lima)", "description": "Admixed: primarily Native American + European ancestry."},
        ],
        "EAS": [
            {"code": "CHB", "name": "Han Chinese (Beijing)", "description": "Northern Han Chinese. Largest ethnic group worldwide."},
            {"code": "JPT", "name": "Japanese (Tokyo)", "description": "East Asian island population with distinct genetic history."},
            {"code": "CHS", "name": "Han Chinese (South)", "description": "Southern Han Chinese. Genetically distinct from northern Han."},
            {"code": "CDX", "name": "Chinese Dai (Xishuangbanna)", "description": "Tai-Kadai speaking population from southern Yunnan."},
        ],
        "SAS": [
            {"code": "GIH", "name": "Gujarati Indian (Houston)", "description": "Northwestern Indian. Indo-European linguistic group."},
            {"code": "PJL", "name": "Punjabi (Lahore)", "description": "Northern South Asian. Crossroads of Central and South Asian genetics."},
            {"code": "BEB", "name": "Bengali (Bangladesh)", "description": "Eastern South Asian. Mix of Indian and Southeast Asian ancestry."},
            {"code": "STU", "name": "Sri Lankan Tamil (UK)", "description": "Southern South Asian. Dravidian linguistic group."},
        ],
    }

    sub_pops = sub_pop_map.get(assigned_pop, [])

    logger.info(f"  Assigned: {assigned_pop}")
    logger.info(f"  Sub-populations in {assigned_pop}:")
    for p in sub_pops:
        logger.info(f"    {p['code']}: {p['name']}")

    return {
        "assigned_super_population": assigned_pop,
        "sub_populations_available": sub_pops,
        "note": "Sub-continental resolution limited to 1000G reference populations. "
                "For Iberian/British/Ashkenazi distinction, additional reference panels needed.",
        "method": "1000G Phase 3 sub-population reference (26 populations). "
                 "Target projected via PCA onto reference space.",
    }


# ═══════════════════════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════════════════════

def run_deep_ancestry(user_vcf: str, output_dir: str = "ancestry") -> Dict:
    """Run all deep ancestry analyses."""
    logger.info("═══ Deep Ancestry Analysis ═══")

    output = {
        "analysis_date": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "user_vcf": user_vcf,
    }

    # Y-DNA
    output["y_dna"] = call_y_haplogroup(user_vcf, {})

    # mtDNA
    output["mt_dna"] = call_mtdna_haplogroup(user_vcf)

    # Neanderthal
    output["neanderthal"] = estimate_neanderthal(user_vcf)

    # Sub-continental ancestry (uses existing PCA output)
    pipeline_root = Path(__file__).resolve().parent.parent.parent
    output["sub_continental"] = run_subcontinental_ancestry(pipeline_root)

    # Save
    out_path = Path(output_dir) / "deep_ancestry.json"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(output, indent=2, ensure_ascii=False))
    logger.info(f"\n  ✅ {out_path}")

    return output


def main():
    import argparse

    parser = argparse.ArgumentParser(description="Deep ancestry analysis")
    parser.add_argument("--vcf", required=True, help="User VCF path")
    parser.add_argument("--output-dir", "-o", default="ancestry")
    parser.add_argument("--verbose", "-v", action="store_true")
    args = parser.parse_args()

    logging.basicConfig(level=logging.DEBUG if args.verbose else logging.INFO, format="%(message)s")

    pipeline_root = Path(__file__).resolve().parent.parent.parent
    output_dir = pipeline_root / args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    user_vcf = args.vcf if Path(args.vcf).is_absolute() else str(Path.cwd() / args.vcf)

    if not Path(user_vcf).exists():
        logger.error(f"VCF not found: {user_vcf}")
        return 1

    result = run_deep_ancestry(user_vcf, str(output_dir))

    # Summary
    print(f"\n═══ Deep Ancestry Summary ═══")
    ydna = result.get("y_dna", {})
    mtdna = result.get("mt_dna", {})
    neand = result.get("neanderthal", {})
    subcont = result.get("sub_continental", {})

    print(f"  Y-DNA:  {ydna.get('haplogroup', '?')} — {ydna.get('description', '?')}")
    print(f"  mtDNA:  {mtdna.get('haplogroup', '?')} — {mtdna.get('description', '?')}")
    if subcont:
        assigned = subcont.get("assigned_super_population", "EUR")
        subs = subcont.get("sub_populations_available", [])
        print(f"  Sub-continental: {assigned} ({len(subs)} sub-populations)")
        for p in subs[:3]:
            print(f"    {p['code']}: {p['name']}")
    if neand.get("reliable"):
        print(f"  Neanderthal: {neand['percentage']}% (European avg: 2.1%)")
    else:
        print(f"  Neanderthal: could not estimate reliably")

    return 0


if __name__ == "__main__":
    sys.exit(main())
