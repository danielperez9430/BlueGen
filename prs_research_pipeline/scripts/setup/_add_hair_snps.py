#!/usr/bin/env python3
"""Add new hair color traits and SNPs to the SNP database."""
import pandas as pd
from pathlib import Path

csv_path = Path(__file__).parent.parent.parent.parent / "prs_research_pipeline" / "data" / "snp_database_annotated.csv"
db = pd.read_csv(csv_path, dtype=str)
print(f"Current rows: {len(db)}")

new_rows = [
    # ── RED HAIR: add 2 more MC1R SNPs ──
    {"rsid": "rs11547464", "gene": "MC1R", "trait_category": "Hair color (red)",
     "effect_allele": "T", "reference_allele": "C", "risk_genotype": "T/T",
     "effect_direction": "+", "weight": "0.30", "evidence_level": "A",
     "pmid": "23579845", "notes": "MC1R red hair variant R151C; strongest red hair predictor after rs1805007",
     "effect_size": "", "odds_ratio": "", "population": "EUR", "chrom": "chr16", "pos": "89985861"},

    {"rsid": "rs1110400", "gene": "MC1R", "trait_category": "Hair color (red)",
     "effect_allele": "C", "reference_allele": "T", "risk_genotype": "C/C",
     "effect_direction": "+", "weight": "0.20", "evidence_level": "A",
     "pmid": "23579845", "notes": "MC1R red hair variant I155T; contributes to RHC phenotype",
     "effect_size": "", "odds_ratio": "", "population": "EUR", "chrom": "chr16", "pos": "89986154"},

    # ── BLONDE HAIR: add TYR, IRF4, HERC2 ──
    {"rsid": "rs1393350", "gene": "TYR", "trait_category": "Hair color (blonde)",
     "effect_allele": "A", "reference_allele": "G", "risk_genotype": "A/A",
     "effect_direction": "+", "weight": "0.25", "evidence_level": "A",
     "pmid": "17999355", "notes": "TYR variant associated with blonde vs brown hair in Europeans",
     "effect_size": "", "odds_ratio": "", "population": "EUR", "chrom": "chr11", "pos": "89011094"},

    {"rsid": "rs12203592", "gene": "IRF4", "trait_category": "Hair color (blonde)",
     "effect_allele": "T", "reference_allele": "C", "risk_genotype": "T/T",
     "effect_direction": "+", "weight": "0.20", "evidence_level": "A",
     "pmid": "27182965", "notes": "IRF4 variant; associated with lighter hair (blonde) in Europeans",
     "effect_size": "", "odds_ratio": "", "population": "EUR", "chrom": "chr6", "pos": "396321"},

    {"rsid": "rs12913832", "gene": "HERC2", "trait_category": "Hair color (blonde)",
     "effect_allele": "C", "reference_allele": "T", "risk_genotype": "C/C",
     "effect_direction": "+", "weight": "0.25", "evidence_level": "A",
     "pmid": "18483556", "notes": "HERC2/OCA2; C allele associated with lighter hair and blue eyes",
     "effect_size": "", "odds_ratio": "", "population": "EUR", "chrom": "chr15", "pos": "28365618"},

    # ── BROWN HAIR: new trait ──
    {"rsid": "rs12913832", "gene": "HERC2", "trait_category": "Hair color (brown)",
     "effect_allele": "T", "reference_allele": "C", "risk_genotype": "T/T",
     "effect_direction": "+", "weight": "0.30", "evidence_level": "A",
     "pmid": "18483556", "notes": "HERC2/OCA2; T allele associated with brown hair and brown eyes",
     "effect_size": "", "odds_ratio": "", "population": "EUR", "chrom": "chr15", "pos": "28365618"},

    {"rsid": "rs16891982", "gene": "SLC45A2", "trait_category": "Hair color (brown)",
     "effect_allele": "C", "reference_allele": "G", "risk_genotype": "C/C",
     "effect_direction": "+", "weight": "0.25", "evidence_level": "A",
     "pmid": "18483556", "notes": "SLC45A2; C allele associated with darker hair (brown vs blonde) in Europeans",
     "effect_size": "", "odds_ratio": "", "population": "EUR", "chrom": "chr5", "pos": "33952685"},

    {"rsid": "rs12203592", "gene": "IRF4", "trait_category": "Hair color (brown)",
     "effect_allele": "C", "reference_allele": "T", "risk_genotype": "C/C",
     "effect_direction": "+", "weight": "0.20", "evidence_level": "A",
     "pmid": "27182965", "notes": "IRF4; C allele associated with darker hair (brown) in Europeans",
     "effect_size": "", "odds_ratio": "", "population": "EUR", "chrom": "chr6", "pos": "396321"},

    # ── BLACK HAIR: new trait ──
    {"rsid": "rs16891982", "gene": "SLC45A2", "trait_category": "Hair color (black)",
     "effect_allele": "G", "reference_allele": "C", "risk_genotype": "G/G",
     "effect_direction": "+", "weight": "0.30", "evidence_level": "A",
     "pmid": "18483556", "notes": "SLC45A2; G allele near-fixed in Africans; darkest hair/skin pigmentation",
     "effect_size": "", "odds_ratio": "", "population": "Global", "chrom": "chr5", "pos": "33952685"},

    {"rsid": "rs1426654", "gene": "SLC24A5", "trait_category": "Hair color (black)",
     "effect_allele": "G", "reference_allele": "A", "risk_genotype": "G/G",
     "effect_direction": "+", "weight": "0.35", "evidence_level": "A",
     "pmid": "16357253", "notes": "SLC24A5; G allele near-fixed in Africans; strongest pigmentation SNP",
     "effect_size": "", "odds_ratio": "", "population": "Global", "chrom": "chr15", "pos": "48426484"},

    {"rsid": "rs1834640", "gene": "SLC24A5", "trait_category": "Hair color (black)",
     "effect_allele": "G", "reference_allele": "A", "risk_genotype": "G/G",
     "effect_direction": "+", "weight": "0.25", "evidence_level": "A",
     "pmid": "24045820", "notes": "SLC24A5 additional variant; dark pigmentation in African and South Asian populations",
     "effect_size": "", "odds_ratio": "", "population": "Global", "chrom": "chr15", "pos": "48426380"},

    # ── EARLY GRAYING: new trait ──
    {"rsid": "rs12203592", "gene": "IRF4", "trait_category": "Hair graying (early)",
     "effect_allele": "T", "reference_allele": "C", "risk_genotype": "T/T",
     "effect_direction": "+", "weight": "0.30", "evidence_level": "A",
     "pmid": "27182965", "notes": "IRF4; strongest GWAS hit for early hair graying; T allele increases risk ~2x per copy",
     "effect_size": "2.0", "odds_ratio": "2.0", "population": "EUR", "chrom": "chr6", "pos": "396321"},
]

new_db = pd.concat([db, pd.DataFrame(new_rows)], ignore_index=True)
new_db.to_csv(csv_path, index=False)
print(f"New rows: {len(new_db)} (+{len(new_rows)})")

print("\n=== New traits ===")
for trait in ["Hair color (brown)", "Hair color (black)", "Hair graying (early)"]:
    rows = new_db[new_db["trait_category"] == trait]
    snps = rows["rsid"].tolist()
    print(f"  {trait}: {len(rows)} SNPs — {snps}")

print("\n=== Expanded traits ===")
for trait in ["Hair color (red)", "Hair color (blonde)"]:
    rows = new_db[new_db["trait_category"] == trait]
    snps = rows["rsid"].tolist()
    print(f"  {trait}: {len(rows)} SNPs — {snps}")

print(f"\nTotal traits: {new_db['trait_category'].nunique()}")
print(f"Total SNPs: {len(new_db)}")
