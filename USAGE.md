# BlueGen — User Guide

**Pipeline version:** see `prs_research_pipeline/scripts/utils/constants.py::PIPELINE_VERSION` (single source of truth) | **Status:** Research-Grade

⏱️ **Quick links:** [Pipeline Stages](#pipeline-stages) · [ClinVar](#clinvar-pathogenic-variants) · [PharmGKB](#pharmacogenomics) · [Dashboard](#streamlit-dashboard)

---

## Quick Start

```bash
cd bluegen

# 1. Setup (once)
python3 -m venv venv && source venv/bin/activate
pip install -r prs_research_pipeline/requirements.txt
# Optional: download genome-wide 1000G reference (~25 GB)
python3 prs_research_pipeline/scripts/setup/download_1000G_full.py --output-dir prs_research_pipeline/reference/1000G_full

# Optional: download AADR archaic reference for Neanderthal/Denisovan analysis (~7 GB)
python3 prs_research_pipeline/scripts/setup/download_aadr_reference.py

# Optional: download Vindija Neanderthal genome for gold-standard archaic analysis (~44 GB)
python3 prs_research_pipeline/scripts/setup/download_vindija_reference.py --all
# Or use pre-built mirror from archive.org (44 GB):
# https://archive.org/details/bluegen-vindija-reference

# Optional: download European sub-continental reference (Ashkenazi Jewish + European, ~200 MB)
bash prs_research_pipeline/scripts/setup/download_ashkenazi_reference.sh
# Or use pre-built mirror from archive.org (~150 MB):
# https://archive.org/details/bluegen-european-reference

# 2. Run pipeline with DeepVariant VCF
python3 prs.py run --full --vcf aligned/E250090601_L01_91_dv.vcf.gz

# 3. View results
open prs_research_pipeline/reports/comprehensive_report_en.html
```

---

## Installation & Prerequisites

### System Requirements

| Component | Minimum | Recommended |
|-----------|---------|-------------|
| Operating System | macOS 12+ / Linux (x86_64/ARM64) | macOS 14+ / Ubuntu 24.04+ |
| Python | 3.10+ | 3.12+ |
| PLINK | Bundled in `tools/plink` (1.9) | — |
| bcftools + tabix | `brew install htslib` | — |
| Memory | 8 GB | 16 GB |
| Disk | 5 GB | 100+ GB (with full 1000G + BAM + VCFs) |
| Network | For reference downloads | Broadband |

### Required Software

```bash
# macOS
brew install htslib

# Linux (Ubuntu/Debian)
sudo apt install bcftools tabix

# Python dependencies
python3 -m venv venv
source venv/bin/activate
pip install -r prs_research_pipeline/requirements.txt
```

### Full 1000 Genomes Reference (one-time)

```bash
# Download genome-wide 1000G Phase 3 (~25 GB)
python3 prs_research_pipeline/scripts/setup/download_1000G_full.py
```

Without this, the platform uses chr22-only PCA (insufficient for proper ancestry inference).

### Input File Requirements

| Format | Requirement |
|--------|------------|
| File type | VCF (.vcf.gz, bgzip-compressed) |
| Genome build | GRCh37/hg19 |
| Variant caller | DeepVariant, GATK, or bcftools |
| Index | .vcf.gz.tbi required |
| Autosomes | chr1–22 |
| Samples | Single or multi-sample supported |

---

## Command Reference

### `python prs.py run` — Execute Full Pipeline

```bash
# Single sample
python prs.py run --full --vcf sample.vcf.gz

# Multi-sample
python prs.py run --full --vcf sample1.vcf.gz,sample2.vcf.gz

# DeepVariant VCF (current for Daniel)
python prs.py run --full --vcf aligned/E250090601_L01_91_dv.vcf.gz

# Dry-run: preview stages
python prs.py run --dry-run

# English-only reports
python prs.py run --full --lang en
```

**Pipeline Stages:**

| Stage | Description | Output |
|-------|-------------|--------|
| A | VCF → PLINK binary | `plink/cohort.{bed,bim,fam}` |
| B | Quality control | `qc/qc_filtered.{bed,bim,fam}` |
| C | LD pruning (cached, parallel by pop) | `plink/ld_pruned_dataset.{bed,bim,fam}` |
| D | 1000G PCA + projection | `pca/target_pcs.eigenvec` |
| — | Ancestry classifier (PCA ensemble) | `pca/ancestry_classification.json` |
| F | PRS computation (PLINK --score, multi-sample) | `prs/prs_raw.csv` |
| G | PCA adjustment | `prs/pca_adjusted_scores.csv` |
| H | Population calibration | `prs/population_calibration_report.json` |
| I | GWAS-LD consistency check | `prs/consistency_check_report.json` |
| J | Uncertainty propagation | `prs/uncertainty_report.json` |
| K | Bilingual interpretation | `interpretations/` |
| SSST | PRS core + unified result | `prs/PRS_RESULT.json` |
| Reports | Comprehensive HTML (EN/ES) | `reports/comprehensive_report_*.html` |

**`--full` also runs:**

| Category | Modules |
|----------|---------|
| External data | PGS Catalog (30 scores), GWAS summary stats, dbSNP annotations |
| Validation | Leakage prevention, global validator (8 dims), adversarial stress tests |
| Benchmarking | GWAS consortium validation, population portability, quality delta |
| Publication | Scientific manuscripts, evidence pack, final score, publication lock |
| Reports | Comprehensive HTML reports + population calibration |

---

### `python prs.py test` — Test Suite

```bash
python prs.py test          # Full: 72 checks across 7 categories
python prs.py test --quick  # JSON integrity only (fast)
```

Categories: JSON Integrity, Data Consistency, Pipeline Artifacts, PRS Edge Cases, Variant Counts, Reports, Code Integrity.

---

### `python prs.py validate` — Scientific Validation

```bash
python prs.py validate
```

Runs: leakage prevention → global validator → genome coverage → adversarial tests → failure mode map → integrity score.

**Output:**
- `science/leakage_audit.json` — 7 leakage checks
- `science/global_validation_report.json` — 8-dimension validation
- `science/adversarial_validation_report.json` — Stress tests
- `science/failure_mode_map.json` — 18 failure modes
- `FINAL_SCIENTIFIC_SCORE.json` — Locked composite score

---

### `python prs.py benchmark` — External Benchmarking

```bash
python prs.py benchmark
```

Runs: PGS Catalog comparison, GWAS consortium validation, population portability, calibration validation, quality delta.

---

### `python prs.py status` — Platform Status

```bash
python prs.py status
```

Displays: pipeline state, SSST status, last run, validation score, publication readiness, output inventory.

---

### `python prs.py pdf` — Generate PDF Reports
```bash
python3 prs.py pdf              # EN + ES PDFs from HTML reports
python3 prs.py pdf --lang en    # English only
```

### `python prs.py run --stage` — Run Single Module
```bash
python3 prs.py run --stage clinvar --vcf sample.vcf.gz     # ClinVar only
python3 prs.py run --stage pharmgkb --vcf sample.vcf.gz     # PharmGKB only
python3 prs.py run --stage ancestry --vcf sample.vcf.gz     # Deep ancestry + Neanderthal
python3 prs.py run --stage archaic --vcf sample.vcf.gz     # Archaic admixture only
python3 prs.py run --stage medgen --vcf sample.vcf.gz       # MedGen enrichment only
```

### ClinVar Pathogenic Variants
Enabled with `--clinvar` or `--full`. Downloads ClinVar VCF (~187 MB) once.
- Evaluates ALL 6.1M WGS variants against ClinVar (4.4M records)
- Reports Pathogenic, Likely_pathogenic, and risk alleles
- Confidence tiers: Expert Panel → Multi-Lab → Single Lab → No Criteria
- Disease descriptions from MedGen (23K clinical definitions)
- **Output:** `clinvar/clinvar_pathogenic_variants.json`
- **Exported VCF:** `clinvar/clinvar_pathogenic.vcf_filtered.vcf.gz` (IGV-ready)

### Pharmacogenomics (PharmGKB/CPIC)
Enabled with `--clinvar` or `--full`. Downloads ClinPGx data (~2 MB) once.
- 51 variant-drug associations (CPIC Level A/B)
- Reports genotype → phenotype → drug recommendation
- 218 CPIC/DPWG clinical guidelines with full recommendation texts
- **Output:** `pharmgkb/pharmgkb_drug_report.json`

### Deep Ancestry
- **mtDNA haplogroup:** Called from chrM variants
- **Y-DNA haplogroup:** Requires chrY coverage (limited in WGS)
- **Sub-continental ancestry:** Real PCA-based classification into European sub-populations (IBS/GBR/CEU/TSI/FIN) using 1000G reference. For Ashkenazi Jewish classification, download the extended European reference panel:
  ```bash
  bash prs_research_pipeline/scripts/setup/download_ashkenazi_reference.sh
  ```
  After download and processing, `subcontinental_pca.py` automatically detects and uses the extended reference.
- **Neanderthal/Denisovan admixture:** Three-tier analysis
  - **Vindija Gold:** Direct comparison against the Vindija 33.19 Neanderthal genome (~30x, hg19). Requires Vindija VCF download (~44 GB) or archive.org mirror.
  - **AADR Silver:** Comparison against Altai/Vindija/Chagyrskaya Neanderthal + Denisova genomes via 1240K panel (1.23M SNPs). Requires AADR download (~7 GB) or pre-built snapshot (~188 MB).
  - **133-SNP Panel (fallback):** Curated archaic introgression SNPs from published literature. Works without any downloads.
- **Output:** `ancestry/deep_ancestry.json`, `pca/subcontinental_assignment.json`

### Streamlit Dashboard
```bash
cd bluegen
./venv/bin/streamlit run dashboard.py
# Opens at http://localhost:8501
```
- 6 interactive pages (Overview, PRS, ClinVar, PharmGKB, Ancestry, Raw Data)
- Filterable tables, Plotly charts, expandable CPIC guidelines
- Reads all JSON outputs — no recomputation

## Output Directory Contract

### `prs/` — PRS Results

| File | Description |
|------|-------------|
| `PRS_RESULT.json` | **Canonical unified PRS output (SSST)** |
| `PRS_RESULT.csv` | Tabular version |
| `prs_raw.csv` | Raw PLINK scores (multi-sample) |
| `pca_adjusted_scores.csv` | PCA-corrected PRS |
| `population_calibration_report.json` | Population-stratified calibration |
| `uncertainty_report.json` | 3-layer uncertainty propagation |
| **`pgs_scores/`** | **PGS Catalog scores + calibrated results** |
| `pgs_scores/pgs_results.csv` | 30 PGS scores against sample |
| `pgs_scores/pgs_calibrated.csv` | Population-calibrated (z-scores + percentiles) |
| `pgs_scores/pgs_calibration_report.json` | Structured calibration report |

### `reports/` — Human-Readable Reports

| File | Description |
|------|-------------|
| `comprehensive_report_en.html` | **Primary report** — 22 sections, PGS calibrated + clinical context |
| `comprehensive_report_es.html` | Spanish version |
| `SCIENTIFIC_MANUSCRIPT_EN.md` | English manuscript |
| `SCIENTIFIC_MANUSCRIPT_ES.md` | Spanish manuscript |

### `science/` — Scientific Metadata

| File | Description |
|------|-------------|
| `ANCESTRY_MODEL.json` | Unified ancestry model |
| `global_validation_report.json` | 8-dimension validation |
| `adversarial_validation_report.json` | Stress tests |
| `failure_mode_map.json` | 18 failure modes |
| `leakage_audit.json` | Leakage prevention |
| `CONSOLIDATION_MANIFEST.json` | SSST enforcement |
| `assumptions.lock.json` | Frozen scientific parameters |

### `benchmark/` — External Validation

| File | Description |
|------|-------------|
| `VALIDATION_REPORT.json` | Unified benchmark |
| `gwas_consortium_validation.json` | GWAS alignment |
| `portability_report.json` | Cross-population bias |
| `quality_delta.json` | Internal vs external gap |

### Root-Level

| File | Description |
|------|-------------|
| `FINAL_SCIENTIFIC_SCORE.json` | Locked integrity score |
| `PUBLICATION_LOCK.json` | Machine-readable lock state |
| `PUBLICATION_LOCK.md` | Publication readiness declaration |

---

## Interpretation Guide

### Understanding PRS Results (Curated Database)

```
PRS_RESULT.json → prs_entries[] for each trait:
  {
    "trait": "Glucose metabolism",
    "raw_score": 0.48,
    "population_zscore": +0.46,     ← (PRS − μ_pop) / σ_pop
    "population_percentile": 67.1,  ← Rank vs EUR reference
    "risk_category": "medium",      ← low / medium / high
  }
```

**Risk categories:** Low (<25th), Medium (25-75th), High (>75th), based on 1000G EUR distribution.

### Understanding PGS Catalog Results (External Validation)

```
pgs_scores/pgs_calibrated.csv:
  pgs_id: PGS000020 (Type 2 Diabetes)
  z_score: +9.2  →  99.99th percentile  →  HIGH risk
  eur_mean: 0.000116  eur_std: 0.000038  (EUR reference)
  n_snps: 7,502  (filter: only scores <500K SNPs are reliable)
```

**PGS z-score interpretation:**
- Z > 2 (97.7th percentile) = Notable genetic risk → investigate clinically
- Z > 1 (84th percentile) = Elevated → monitor lifestyle factors
- -1 < Z < +1 = Population average
- Z < -1 = Lower genetic burden
- Z < -2 = Protective genotype

**IMPORTANT:** All results are RESEARCH USE ONLY — not clinically validated.

---

## Data Pipeline: FASTQ → VCF → PRS

```
FASTQ (.fq.gz)                         raw_data/fastq/
  → FastQC                             raw_data/qc/
  → BWA-MEM (hg19)                     aligned/*.bam (73 GB)
  → DeepVariant (CNN, 16 shards)       aligned/*_dv.vcf.gz (6.1M variants)
  → bcftools mpileup (backup)          aligned/*_raw.vcf.gz (5.0M variants)
  → BlueGen Pipeline                       prs_research_pipeline/prs/
  → PGS Catalog (30 scores)            prs_research_pipeline/prs/pgs_scores/
  → Population Calibration             1000G Phase 3 (2,504 samples, 5 super-pops)
  → HTML Report                        reports/comprehensive_report_en.html
```

---

## Troubleshooting

| Problem | Solution |
|---------|----------|
| `plink: command not found` | Download PLINK 1.9 from https://www.cog-genomics.org/plink/ and place in `tools/plink` |
| `VCF not found` | File must be bgzip-compressed (.vcf.gz) with .tbi index |
| `ModuleNotFoundError: numpy` | `source venv/bin/activate` from project root |
| `PCA projection failed` | Download 1000G: `python3 scripts/setup/download_1000G_full.py` |
| `PRS: 0 scores` | Check SNP database has chrom/pos data (`data/snp_database_annotated.csv`) |
| `ClinVar download fails` | Check network. Re-run: `python3 scripts/setup/download_clinvar.py --force` |
| `MedGen enrichment fails` | Run `python3 scripts/clinical/medgen_enrich.py --download` first |
| `tabix/bcftools not found` | macOS: `brew install htslib`. Linux: `sudo apt install bcftools tabix` |
| `Invalid chromosome code` | Non-standard contigs (e.g., `chr6_ssto_hap7`) — scripts use `--allow-extra-chr` |
| `Duplicate ID in PLINK` | Auto-handled by dedup logic in PRS computation |
| `Stage C OOM` | Reduce parallel jobs or increase `--memory` in `scripts/stages/03_ld_ancestry_prune.sh` |
| `DeepVariant OOM` | Reduce `--num_shards` or increase Docker memory allocation |
