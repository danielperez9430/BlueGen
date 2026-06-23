# Changelog

All notable changes to BlueGen.

## [1.1.0] — 2026-06-12

### Added
- **`scripts/prs/prs_plink_score.py`:** Extracted PRS computation from inline code into standalone module
- **`scripts/publication/bilingual_report_generator.py`:** Extracted bilingual HTML/PDF report generation from inline code
- **`scripts/utils/constants.py`:** Centralized paths, traits, thresholds, and version constants
- **`scripts/utils/tool_detection.py`:** Centralized PLINK/bcftools/tabix detection with version validation
- **`scripts/utils/config_validator.py`:** Schema validation for config.yaml (types, ranges, allowed values)
- **`scripts/utils/logging_config.py`:** Unified logging with colored console output, file rotation, and StageProgress context manager
- **`tests/`:** 28 pytest unit tests for constants, tool_detection, config_validator, logging_config
- **`pyproject.toml`:** Ruff linting configuration (Python 3.10+, line length 120)
- **`require_output()`:** Pipeline helper to validate stage outputs exist before proceeding
- **`OPENGWAS_TOKEN` env var:** Environment variable support for OpenGWAS API token (fallback to file)

### Changed
- **`prs.py`:** PLINK detection now uses `utils.tool_detection` instead of hardcoded paths
- **`prs.py`:** SNP_DB path centralized as module-level constant (eliminated 3 duplications)
- **`prs.py`:** Config validated at pipeline startup via `config_validator.py`
- **`prs.py`:** Stages F-H now have `required=True` and output validation between stages
- **`prs.py`:** Stage G now checks for `pca_adjusted_scores.csv` (matching actual script output)
- **CI/CD:** Python version matrix (3.10, 3.11, 3.12), ruff linting, pytest unit tests
- **CI/CD:** Removed `continue-on-error: true` — tests now fail the build
- **All new modules:** Comprehensive docstrings with input/output schemas and usage examples

### Removed
- **`prs_research_pipeline/run.sh`:** Legacy bash orchestrator (replaced by `prs.py` as single entry point)

### Fixed
- **`03_ld_ancestry_prune.sh`:** Added `--allow-extra-chr` to all PLINK commands (fixes `chr6_ssto_hap7` error with DeepVariant VCFs)
- **`population_calibrate_v2.py`:** Added `@dataclass` decorator to `PopulationDistribution` class
- **`population_calibrate_v2.py`:** Normalizes JSON keys (`q25`→`percentile_25`, `q75`→`percentile_75`, `p5`→`percentile_5`, `p95`→`percentile_95`)
- **`population_calibrate_v2.py`:** Computes `iqr` from `percentile_75 - percentile_25` when missing
- **`population_calibrate_v2.py`:** Filters unknown JSON keys (e.g., `method`) before constructing dataclass
- **`gwas_summary_stats.py`:** Added `OPENGWAS_TOKEN` env var support (was file-only)
- **`population_calibrate_v2.py`:** Removed duplicate `@dataclass` decorator on `PopulationDistribution`
- **`README.md`:** Added BlueGen logo (`assets/blugen-logo.png`)
- **`LICENSE`:** Added MIT license file

---

## [1.0.0] — 2026-06-09

### Added
- **Project renamed:** BlueGen — Personal Genomics Platform powered by PRSKit
- **LICENSE:** MIT
- **.gitattributes:** LFS-ready, cross-platform line endings
- **Root README:** badges, quick start, module table for GitHub
- **Cross-platform support:** auto-detect project root, works from any CWD
- **`--stage` flag:** run individual modules (clinvar, pharmgkb, ancestry, medgen)
- **PDF reports:** `python prs.py pdf` via WeasyPrint
- **CI/CD:** GitHub Actions (syntax check, lint)
- **HWE fallback calibration:** traits not in 1000G get Hardy-Weinberg estimates
- **PRS expanded:** 10→56 traits, 107→179 SNPs (wellness, personal traits, disease risk)
- **PGS Catalog expanded:** 30→54 scores
- **Dashboard:** 6-page Streamlit interactive app
- **VCF/TSV export:** ClinVar pathogenic variants as annotated files
- **Debug log rotation:** auto-archive at 1MB, fresh log per run
- **`--update-references` flag:** one-command refresh of all databases
- **Evidence level colors:** A/B/C/D badges with color coding in reports
- **CPIC guideline explainer:** what CPIC/DPWG mean in pharmacogenomics section
- **Per-section limitation notes:** PRS, PGS, ClinVar each explain their caveats
- **raw_data/ structure:** documented directory layout for user data

### Changed
- **README:** version numbers removed from titles — just "BlueGen"
- **USAGE.md:** new sections for ClinVar, PharmGKB, Deep Ancestry, Dashboard, PDF
- **CI tests:** 1000G MISSING → warn instead of fail; `continue-on-error: true`
- **Reports:** Executive Summary grid fixed to 3 columns; evidence levels explained
- **Summary output:** shows HTML + PDF + MD + Dashboard paths
- **Project structure:** references unified, root junk cleaned, VCFs moved to raw_data/
- **Legacy scripts:** 3 referenced scripts moved to proper dirs; legacy/ gitignored

### Removed
- **ClinVar VCF from git:** 187 MB → auto-downloaded on first run
- **PLINK from git:** platform-specific binary → users download for their OS
- **SSH keys:** `clouding`/`clouding.pub` removed from tracking
- **Deprecated scripts:** `archive/` and `legacy/` gitignored

### Fixed
- Trait names with `/` crashing PRS score file creation
- Chromosome prefix mismatch (`chr1` vs `1`) in ClinVar annotator
- Debug log unbounded growth
- CI test failures on clean clone
- Missing `ancestry_pcs_used` attribute in methods generator

## [0.9.0-beta] — 2026-06-08

### Added

#### ClinVar Pathogenic Variant Annotation
- **New module:** `scripts/clinical/clinvar_annotator.py` — genome-wide ClinVar annotation
- Evaluates all 6.1M WGS variants against ClinVar (GRCh37, 4.4M records)
- Two-pass algorithm: `tabix` positional overlap → exact REF/ALT matching
- Reports variants classified as `Pathogenic`, `Likely_pathogenic`, `Pathogenic/Likely_pathogenic`, and risk alleles
- **New module:** `scripts/setup/download_clinvar.py` — downloads ClinVar VCF + index from NCBI FTP (~187 MB)
- **New section:** ClinVar — Pathogenic Variants in comprehensive HTML report
- Bilingual (EN/ES) with confidence tiers, disease descriptions, and review status
- **Flag:** `--clinvar` (auto-enabled with `--full`)
- **Output:** `clinvar/clinvar_pathogenic_variants.json`

#### ClinVar Confidence Tiers
- Each variant classified by evidence quality:
  - 🏅 High (expert panel, practice guideline) → 1 variant found
  - ✓ Moderate (multiple submitters agree) → 9 variants found
  - ⚠️ Low (single submitter) → 14 variants found
  - ❓ Very Low (no evidence criteria) → 83 variants found
- Veracity alert banner showing how many findings have strong evidence
- Separate tables for high/moderate vs low/very-low confidence variants

#### MedGen Disease Descriptions
- **New module:** `scripts/clinical/medgen_enrich.py` — enriches ClinVar variants with disease definitions
- Downloads MedGen RRF files from NCBI FTP (~10 MB, free)
- 23,059 disease names → 23,187 clinical definitions
- Local lookup: instant after first download, no API calls needed
- **80/107 variants enriched** (75%) with clinically curated definitions
- **Flag:** `--download` to fetch/update MedGen, `--check-update` to verify freshness
- **Output:** `reference/medgen/{NAMES,MGDEF,MGSTY}.RRF.gz`

#### Pharmacogenomics (PharmGKB/ClinPGx)
- **New module:** `scripts/clinical/pharmgkb_annotator.py` — annotates user VCF against pharmacogenomic variants
- Curated CSV database: 28 CPIC Level A/B variant-drug associations
- Reports genotype → phenotype → drug recommendation (bilingual EN/ES)
- **New module:** `scripts/clinical/clinpgx_sync.py` — downloads ClinPGx datasets
- Downloads: clinicalVariants.zip (5,190 pairs), guidelineAnnotations (218 CPIC/DPWG), variants.zip (7,615)
- Free, no login required. CC BY-SA 4.0 license
- **New section:** Pharmacogenomics — Drug Response in comprehensive HTML report
- Actionability levels: Critical, Important, Informative
- **Output:** `pharmgkb/pharmgkb_drug_report.json`, `pharmgkb/clinpgx_parsed.json`

#### Reference Management
- **Flag:** `--update-references` — force download of ClinVar, MedGen, and ClinPGx
- Manifest-based caching with staleness checks
- MedGen auto-check: warns if database >30 days old
- Rate limiting: respects ClinPGx 2 req/sec limit

#### Per-Section Limitations
- Added inline limitation notes to PRS, PGS Catalog, and ClinVar sections
- Each section now explains what it CAN and CANNOT tell you

### Changed
- **prs.py:** Added `--clinvar`, `--update-references` flags. 4 new pipeline stages.
- **comprehensive_report.py:** 3 new sections (ClinVar, Pharmacogenomics, enhanced with descriptions)
- **test_suite.py:** Added TEST 8 (ClinVar validation), 22 JSON files tracked
- **README.md:** New sections: Keeping References Updated, External Data (ClinVar, MedGen, ClinPGx)
- **.gitignore:** Exceptions for ClinVar and MedGen reference data, pharmgkb outputs
- **Memory:** Added PRS Research Pipeline memory with project context

### Data Sources Added
| Source | Size | License | Update Frequency |
|--------|------|---------|-----------------|
| ClinVar VCF (GRCh37) | 187 MB | Public Domain (NCBI/NIH) | Monthly |
| MedGen RRF | 10 MB | Public Domain (NCBI/NIH) | Weekly |
| ClinPGx/PharmGKB | 1.8 MB | CC BY-SA 4.0 | Weekly |
| dbSNP (via NCBI Entrez) | API | Public Domain (NCBI/NIH) | On demand |

### Fixed
- `22_methods_generator.py`: Added missing `ancestry_pcs_used` attribute to `MethodsContext`
- `prs.py`: Fixed `SyntaxWarning` for `\s` escape in f-string
- `clinvar_annotator.py`: Fixed chromosome name mismatch (`chr1` vs `1`) between user VCF and ClinVar

---

## [10.0.0] — 2026-06-07

### Added
- Full 1000 Genomes Reference (84M variants, genome-wide)
- dbSNP annotations (108/109 rsIDs via NCBI Entrez)
- GWAS VCF datasets (10 studies, 525 records)
- PCA Ancestry Classifier (ensemble: centroid + k-NN)
- Multi-sample PRS (PLINK --score, N samples native)
- FASTQ → BWA-MEM alignment (73 GB BAM, 99% mapped, 35x)
- bcftools variant calling (5.0M variants)
- DeepVariant CNN variant calling (6.1M variants, 4.7M SNPs)
- PGS Catalog integration (30 scores, 3.87M SNP-trait associations)
- PGS Population Calibration (24 reliable scores, EUR percentiles)
- Comprehensive HTML reports (22 sections, clinical context)
- Code integrity tests (72 checks, 7 categories)
- Duplicate ID handling in PLINK bim files
- `--allow-extra-chr` throughout pipeline stages
