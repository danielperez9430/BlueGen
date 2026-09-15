# BlueGen

Polygenic Risk Score research platform with 1000 Genomes PCA projection,
population-calibrated scoring, 8-dimension scientific validation,
SSST (Single Source of Scientific Truth) consolidation, and
publication-ready bilingual reporting (EN/ES).

## Quick Start

```bash
# One-time setup
python3 -m venv venv && source venv/bin/activate
pip install -r requirements.txt
python3 scripts/setup/download_1000G_full.py    # ~25 GB, run once

# Full pipeline (PRS + ClinVar + pharmacogenomics + validation + reports)
python3 prs.py run --full --vcf sample.vcf.gz

# Update reference databases before important runs
python3 prs.py run --full --update-references --vcf sample.vcf.gz

# Multi-sample
python3 prs.py run --full --vcf sample1.vcf.gz,sample2.vcf.gz

# Run test suite (90+ checks, 9 categories)
python3 prs.py test

# Platform status
python3 prs.py status
```

## Commands

| Command | Description |
|---|---|
| `run --full --vcf <file.vcf.gz>` | Full pipeline: VCF → PRS → ClinVar → PharmGKB → reports |
| `run --full --update-references --vcf <file>` | Full pipeline + refresh all reference databases |
| `run --clinvar --vcf <file>` | PRS + ClinVar pathogenic variants only |
| `run --full --vcf <a.vcf.gz,b.vcf.gz>` | Multi-sample mode (comma-separated) |
| `validate` | Scientific validation suite |
| `report` | Generate bilingual manuscripts (EN/ES) |
| `benchmark` | External benchmarking |
| `test` | Test suite — 80+ checks across 8 categories |
| `test --quick` | JSON integrity only (fast) |
| `status` | Platform dashboard |
| `audit` | Export peer review package |

## Pipeline Stages

```
FASTQ ──→ [optional] BWA-MEM → DeepVariant → genome-wide VCF
                                              ↓
VCF (GRCh37/38) or array export → [Gate] build detect / liftover / array→VCF
                                              ↓
VCF → [A] PLINK → [B] QC → [C] LD Prune → [D] PCA + Ancestry Classifier
                                    ↓
    [F] PRS joint with 1000G → [G] PCA Adjust → [H] Pop Calibrate (per-run distributions)
                                    ↓
    [ClinVar] Pathogenic Variants → [MedGen] Disease Descriptions
                                    ↓
    [PharmGKB] Drug Response → [ClinPGx] CPIC Guidelines
                                    ↓
    [7] Scientific Freeze → [8] Corrections → [9] SSST → [10] Publication Lock
```

**Input gate (v3.0):** genome build from the VCF header (contig lengths), GRCh38 lifted to GRCh37 with the UCSC chain and hg19 REF check; array exports (23andMe/AncestryDNA/MyHeritage/FTDNA) converted to a VCF with every assayed site.
**QC:** missingness filters always; MAF/HWE only for cohorts (≥50 samples) — on one sample they would delete every homozygous site.
**PRS:** the user is merged into the 1000G reference and every trait scored with one PLINK `--score` for user and reference together (WGS: absent sites = hom-ref; array: genotyped sites only), dosage-weighted Σ(βⱼ × Gᵢⱼ), multi-sample native.
**LD Pruning:** Per-population parallel, conservative intersection, cached by reference + parameters.
**PCA:** 1000G-trained PCA with target projection, 20 PCs genome-wide.
**Ancestry:** PCA ensemble classifier (centroid + k-NN), 5 super-populations.
**Calibration:** 1000G distributions per trait × super-population built in the same run on the same variant set (2,504 samples); raw score z-scored against the inferred super-population; discrete few-SNP distributions get empirical percentiles and a low-confidence flag.
**PGS Catalog:** same joint scoring per published score, calibrated against the inferred super-population, with per-score coverage.

## Project Structure

```
BlueGen/
├── prs.py                          # Main CLI orchestrator
├── .gitignore                      # Excludes sensitive/large data
├── raw_data/                       # Raw FASTQ files (git-ignored)
│   ├── fastq/                      #   .fq.gz paired-end reads
│   └── qc/                         #   FastQC reports
├── aligned/                        # BAM alignment outputs (git-ignored)
├── tools/                          # PLINK binaries (git-ignored)
│   ├── plink                       #   PLINK 1.9
│   └── plink2                      #   PLINK 2.0
├── venv/                           # Python venv (git-ignored)
└── prs_research_pipeline/
    ├── scripts/
    │   ├── stages/                 # Genotype processing (A-D)
    │   │   ├── 01_vcf_to_plink.sh
    │   │   ├── 02_quality_control.sh
    │   │   ├── 03_ld_ancestry_prune.sh
    │   │   ├── 04_pca_1000G.sh
    │   │   └── pca_ancestry_classifier.py
    │   ├── prs/                    # PRS computation (F-H)
    │   │   ├── 06_prs_compute.py
    │   │   ├── pca_adjust_v2.py
    │   │   ├── population_calibrate_v2.py
    │   │   ├── prs_multi_method_v2.py
    │   │   └── 33_ancestry_aware_normalization.py
    │   ├── validation/             # Phase 7-8 corrections
    │   │   ├── 13_gwas_ld_consistency_check.py
    │   │   ├── 14_uncertainty_propagation.py
    │   │   ├── 16_reproducibility_engine.py
    │   │   ├── 18_scientific_lock.py
    │   │   ├── 30_snp_universe_registry.py
    │   │   ├── 31_leakage_prevention.py
    │   │   ├── 32_global_scientific_validator.py
    │   │   ├── 34_scientific_integrity_score.py
    │   │   └── bilingual_interpretation.py
    │   ├── sss/                    # Phase 9 SSST consolidation
    │   │   ├── 36_prs_core.py
    │   │   ├── 37_prs_result_unified.py
    │   │   ├── 38_benchmark_reinterpretation.py
    │   │   ├── 39_ancestry_model_unified.py
    │   │   ├── 40_leakage_integrated.py
    │   │   ├── 41_unified_report_engine.py
    │   │   └── 42_consolidation_manifest.py
    │   ├── publication/            # Phase 10 + reports
    │   │   ├── 43_adversarial_prs_validation.py
    │   │   ├── 44_failure_mode_map.py
    │   │   ├── 46_final_scientific_score.py
    │   │   ├── 47_publication_lock.py
    │   │   └── comprehensive_report.py
    │   ├── benchmarking/           # External benchmarks
    │   │   ├── 25_gwas_consortium_validation.py
    │   │   ├── 26_population_portability_test.py
    │   │   ├── 27_real_world_calibration.py
    │   │   ├── 29_quality_delta_analysis.py
    │   │   ├── gwas_summary_stats.py
    │   │   └── pgs_catalog_integration.py
    │   ├── setup/                  # One-time download scripts
    │   │   ├── download_1000G_full.py
    │   │   ├── download_dbsnp.py
    │   │   ├── download_clinvar.py
    │   │   └── GWAS_DATASETS.md
    │   ├── clinical/               # Clinical annotation
    │   │   └── clinvar_annotator.py
    │   ├── utils/                  # Tools and test suite
    │   │   ├── build_reference_distributions.py
    │   │   └── test_suite.py
    │   └── legacy/                 # Older modules (still usable)
    ├── data/                       # Curated SNP databases
    │   ├── snp_database.csv
    │   ├── snp_database_annotated.csv
    │   ├── snp_database_dbsnp_enriched.csv
    │   └── dbsnp_annotations.json
    ├── reference/                  # Reference data
│   ├── medgen/                 #   MedGen disease definitions (~10 MB, committed)
│   ├── clinvar/                #   ClinVar VCF (GRCh37, ~187 MB, committed)
│   │   ├── clinvar.vcf.gz      #     Germline classifications
│   │   ├── clinvar.vcf.gz.tbi  #     Tabix index
│   │   ├── manifest.json       #     Download metadata
│   │   └── SOURCES.md          #     Attribution & license
│   └── 1000G_full/             #   1000 Genomes (git-ignored, ~25 GB)
    ├── plink/                      # Stage A-C outputs (git-ignored)
    ├── qc/                         # Stage B outputs (git-ignored)
    ├── pca/                        # Stage D outputs (git-ignored)
    ├── prs/                        # PRS computation outputs (git-ignored)
    ├── benchmark/                  # External benchmarking (git-ignored)
    ├── science/                    # SSST validation outputs (git-ignored)
    ├── reports/                    # Generated reports (git-ignored)
    ├── gwas/                       # GWAS extracted data (git-ignored)
    ├── reproducibility/            # Run fingerprints (git-ignored)
    └── pgs/                        # PGS Catalog cache (git-ignored)
```

## Key Outputs

| File | Description |
|---|---|
| `prs/PRS_RESULT.json` | Unified PRS output (all traits) |
| `prs/PRS_RESULT.csv` | PRS results as CSV |
| `prs/prs_raw.csv` | Raw PLINK scores per (sample, trait), scored jointly with 1000G |
| `prs/prs_reference_raw.csv` | The same scores for the 2,504 reference individuals |
| `prs/pca_adjusted_scores.csv` | PCA-adjusted PRS (reported; calibration uses the raw score) |
| `prs/population_calibrated_v2.csv` | Population-stratified calibration (z-score, percentile, low_confidence per sample × trait) |
| `prs/pgs_scores/pgs_calibrated.csv` | PGS Catalog scores calibrated vs the inferred super-population, with coverage |
| `reproducibility/input_build.json` | Genome build / liftover / array conversion record |
| `prs/uncertainty_report.json` | 3-layer variance propagation |
| `pca/ancestry_classification.json` | PCA ensemble ancestry call |
| `science/ANCESTRY_MODEL.json` | Canonical ancestry model |
| `science/global_validation_report.json` | 8-dimension validation |
| `science/CONSOLIDATION_MANIFEST.json` | SSST manifest |
| `science/adversarial_validation_report.json` | Adversarial stress test |
| `science/failure_mode_map.json` | 18 failure modes mapped |
| `FINAL_SCIENTIFIC_SCORE.json` | Locked integrity score |
| `PUBLICATION_LOCK.md` | Publication readiness declaration |
| `clinvar/clinvar_pathogenic_variants.json` | Pathogenic variant annotation (ClinVar) |
| `reports/comprehensive_report_en.html` | Interactive report (EN); `prs.py pdf` adds cover + TOC |
| `reports/comprehensive_report_es.html` | Interactive report (ES) |
| `reports/comparison_report_{en,es}.html` | Samples side by side (multi-sample runs / `--compare`) |
| `reports/SCIENTIFIC_MANUSCRIPT_EN.md` | Publication manuscript (EN) |
| `reports/SCIENTIFIC_MANUSCRIPT_ES.md` | Publication manuscript (ES) |

## Performance

Measured on a genome-wide DeepVariant VCF (4.1 M variants) with the full 1000G reference, Apple Silicon, 2026-09-15 (v3.0):

| Stage | Cold | Cached / typical |
|---|---|---|
| Gate (build detect; liftover only for GRCh38) | 1 s (+9 s per 1 M records lifted) | — |
| A (VCF → PLINK) | 27 s | — |
| B (QC) | 7 s | — |
| C (LD Prune, parallel) | 45 min | 5 s (cache keyed by reference + parameters) |
| D (PCA + projection) | 7–53 min (varies run to run; I/O-bound) | — |
| F–H (joint PRS + calibration) | 1 min | — |
| ClinVar + PharmGKB + deep ancestry | ~8 min | — |
| 7–10 (Validation + SSST) | 2 min | — |
| PGS Catalog joint calibration (56 scores, 7.3 M variants × 2,505 samples) | 37–55 min | — |
| Reports (HTML + PDF) | 1 min | — |
| **Total `--full` (genome-wide VCF, C cached)** | | **78–87 min** |

## Data Requirements

| Resource | Size | Source |
|---|---|---|
| Input VCF | 10 MB – 2 GB | DeepVariant or standard VCF (GRCh37/hg19) |
| 1000 Genomes Phase 3 | ~25 GB | Auto-downloaded by `download_1000G_full.py` |
| Working space | ~30 GB | Pipeline intermediates |
| RAM | 16 GB | Recommended |
| Python 3.10+ | — | With venv |
| PLINK 1.9 | — | Not bundled - download separately, place in `tools/` (see System Tools in the root README) |
| bcftools, tabix | — | `brew install htslib` (macOS) or `sudo apt install bcftools tabix` (Linux) |

## External Data

| Dataset | Status | Source |
|---|---|---|
| **1000 Genomes** | Genome-wide (84M variants) | EBI FTP — `download_1000G_full.py` |
| **dbSNP** | 108/109 rsIDs annotated | NCBI Entrez — `download_dbsnp.py` |
| **GWAS VCFs** | 10 datasets, 525 records | OpenGWAS — manual (see `GWAS_DATASETS.md`) |
| **PGS Catalog** | 57 downloaded, 52 scored + population-calibrated | PGS Catalog REST API — auto-downloaded |
| **PharmGKB** | Detoxification variants | PharmGKB API |
| **ClinVar** | Pathogenic/likely pathogenic variants | NCBI FTP — `download_clinvar.py` |
| **MedGen** | Disease definitions (23K concepts) | NCBI FTP — `medgen_enrich.py --download` |

## Keeping References Updated

Reference databases are updated periodically by their sources. The pipeline
downloads them on first use and caches them locally. To keep them current:

```bash
# Check if references are up to date
python3 scripts/clinical/medgen_enrich.py --check-update

# Update all references (ClinVar + MedGen)
python3 prs.py run --clinvar --update-references --vcf sample.vcf.gz

# Or update individually
python3 scripts/clinical/medgen_enrich.py --download       # MedGen (weekly)
python3 scripts/setup/download_clinvar.py --force           # ClinVar (monthly)
```

| Reference | Update Frequency | Command | Size |
|---|---|---|---|
| ClinVar VCF | Monthly (1st week) | `download_clinvar.py` | 187 MB |
| MedGen RRF | Weekly (Wednesdays) | `medgen_enrich.py --download` | 10 MB |
| 1000 Genomes | Static (Phase 3) | `download_1000G_full.py` | 25 GB |
| dbSNP | On demand | `download_dbsnp.py` | API |

**Pro tip:** Run `--update-references` once a month before generating a final report.

## Scientific Methods

- **PRS:** Σ(βⱼ × Gᵢⱼ) — dosage-weighted scoring via PLINK `--score`
- **LD Pruning:** Ancestry-matched per-population, conservative intersection (`--indep-pairwise 50 5 0.2`)
- **PCA:** 1000G-trained PCA with target projection, 20 PCs, 2M LD-pruned genome-wide variants
- **Ancestry:** PCA ensemble classifier (centroid + logistic + k-NN), 5 super-populations (EUR, AFR, EAS, SAS, AMR)
- **Calibration:** Empirical 1000G population-stratified distributions with z-score normalization
- **Uncertainty:** 3-layer variance propagation (genotype + ancestry + effect size)
- **Validation:** 8 scientific dimensions, adversarial stress testing, leakage audit
- **Reproducibility:** Deterministic seeds, SHA-256 hashing, environment fingerprint

## Test Suite

```bash
python3 prs.py test          # Full: 80+ checks across 8 categories
python3 prs.py test --quick  # JSON integrity only (fast)
```

Categories:
1. **JSON Integrity** — 22 output files, no numpy types
2. **Data Consistency** — PRS results, validation scores, integrity scores
3. **Pipeline Artifacts** — 16 critical files
4. **PRS Edge Cases** — multi-sample detection, data types, NaN detection
5. **Variant Counts** — cross-stage variant/sample tracking
6. **Report Content** — HTML size and structure
7. **Code Integrity** — 10 key scripts compile without syntax errors
8. **ClinVar** — pathogenic variant counts, field validation, summary consistency

Catches regressions: numpy types in JSON, broken indentation, missing outputs, variant mismatches.

## Git Repository

What's committed vs ignored:

| Committed ✅ | Ignored ❌ |
|---|---|
| Source code (`scripts/`) | Raw FASTQ (`raw_data/`) |
| SNP databases (`data/`) | BAM files (`aligned/`) |
| Config, docs, ClinVar ref | Reference data (`reference/1000G_full/`) |
| GWAS dataset index | PGS cache (`pgs/`) |
| Test suite | Pipeline outputs (`plink/`, `qc/`, `prs/`, etc.) |
| | Virtual environment (`venv/`) |
| | Tools (`tools/`) |
| | API tokens (`.opengwas_token`) |

## Dependencies

`pip install -e .` from the repository root installs everything in `requirements.txt` (pandas, numpy, scipy, scikit-learn, matplotlib, cyvcf2, pyliftover, jinja2, weasyprint, pydantic, pyyaml, click, tqdm, …); `pip install -e ".[dashboard]"` adds streamlit + plotly. `requirements.lock` pins the exact validated versions (used by the Docker image).

System: `plink` 1.9 (**not bundled** in the clone — see the root README "System Tools"; it is in the Docker image), `bcftools`, `tabix` (htslib), optionally `plink2`.
All reference databases and the UCSC liftOver chains are downloaded on first use — no API keys required.

## Data Directory Structure

```
BlueGen/
├── raw_data/                   # Your input files
│   ├── fastq/                  #   .fq.gz paired-end reads
│   ├── qc/                     #   FastQC reports
│   ├── vcf/                    #   Input VCF files (.vcf.gz + .tbi)
│   └── deepvariant/            #   DeepVariant output VCFs
├── aligned/                    # BAM alignment files
└── prs_research_pipeline/
    └── reference/              # All reference databases
        ├── clinvar/            #   ~187 MB (committed)
        ├── medgen/             #   ~10 MB (committed)
        ├── clinpgx/            #   ~2 MB (committed)
        ├── population_distributions/ # PRS distributions (committed)
        └── 1000G_full/         #   ~25 GB (git-ignored, download once)
```

## Cross-Platform Notes

- **macOS:** `brew install htslib` for bcftools/tabix
- **Linux:** `sudo apt install bcftools tabix` (Ubuntu/Debian)
- **Windows:** Use WSL2 (Ubuntu recommended)
- **PLINK:** Not bundled (macOS/Linux/Windows alike) - `tools/` is git-ignored. Download the build for your OS and place it at `tools/plink`/`tools/plink2` (see root README's System Tools table).
- **Python:** 3.10+ with venv. All packages in `requirements.txt`.
- All scripts use `sys.executable` — works regardless of venv activation.
- Pipeline auto-detects project root from `__file__` — works from any CWD.

## Version Compatibility

The pipeline version is defined once, in `scripts/utils/constants.py` (`PIPELINE_VERSION`),
and referenced everywhere else (`prs.py`, this README, `config.yaml`, and every script
that stamps a `pipeline_version` field into its output). CI fails if any of those drift
from each other.

| Version | Key Features | Min Python |
|---------|-------------|------------|
| 3.0.1 | Evidence-level legend in the report, docs synced; plus 3.0.0: GRCh38 input (auto-detect + liftover), genotyping-array input (`--raw`), joint user+1000G scoring for the curated PRS and the PGS Catalog, single-sample QC fix, multi-sample comparison report, PDF cover/TOC, automated PubMed citation audit, multi-arch Docker, `bluegen` CLI, 9-page dashboard, ClinVar, PharmGKB, MedGen, Deep Ancestry, bilingual reports | 3.10+ |

### Upgrading
```bash
# Pull latest code, re-download expanded references
git pull
python3 prs.py run --update-references --vcf sample.vcf.gz
```

### Reference data versioning
Each reference directory contains a `manifest.json` with download date and version.
Run `--update-references` to refresh all databases.

## Troubleshooting

| Problem | Solution |
|---|---|
| `NameError: result not defined` | Update `prs.py` to latest version |
| PRS: 0 scores | Check chrom/pos in SNP database; use annotated version |
| `Duplicate ID` in PLINK | Handled automatically by dedup (numbered suffixes) |
| Stage C timeout | Increase timeout in `prs.py` or use cached SNPs |
| VCF not recognized | Stage A auto-detects DeepVariant vs standard VCF |
| venv not detected | `_venv.py` detects via `sys.path`; run with `python3 prs.py` |

## Internal Notes (Daniel)

```bash
# Full pipeline with DeepVariant VCF (procesado en VPS 2026-06-07)
cd bluegen
./venv/bin/python prs.py run --full --vcf aligned/E250090601_L01_91_dv.vcf.gz

# Backup VCF location
ls raw_data/deepvariant/E250090601_L01_91_dv.vcf.gz
```
