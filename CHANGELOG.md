# Changelog

All notable changes to BlueGen.

## [Unreleased]

Nothing yet.

## [2.0.0] — 2026-08-08

_First formally tagged release. `PIPELINE_VERSION` had already been `2.0.0` since the version-unification fix below (TIER 0.1) - this release just gives that version an actual git tag/GitHub Release for the first time, covering everything accumulated since the previous (untagged) `1.2.1`._

### Added
- **`bluegen/` importable library package**: `calibration.py` (Stage H, `PopulationCalibrationV2`), `ancestry.py` (Stage G, `PCAAdjustmentV2`), `scoring.py` (Stage F, PLINK PRS computation), `schemas.py` (pydantic contracts, see below). The numbered scripts under `scripts/prs/` are now thin argparse wrappers importing from this package instead of standalone shell-out logic — Stages F–H are now unit-testable Python, not just subprocess black boxes. 26 new tests (`test_pca_adjustment.py`, `test_scoring.py`, plus `test_population_calibration.py` updated in place).
- **Jinja2 templating for the HTML report**: `comprehensive_report.py` (2,998 lines of inline f-string HTML) refactored into `report/data_loader.py` (data loading), `report/interpretations.py` (19 pure helpers), `report/render.py` (Jinja2 rendering) + `templates/*.j2` partials — one per report section. File dropped to 1,852 lines; adding a new section is now "one template + one function + one registration line" instead of editing a monolith. 49 new render-regression tests (`test_comprehensive_report_render.py`) verify zero unrendered template artifacts and byte-identical output (whitespace-normalized) at every migration step.
- **`bluegen/schemas.py`**: pydantic shape contracts for the 4 critical inter-stage JSON artifacts (`PRS_RESULT.json`, `ANCESTRY_MODEL.json`, `calibration_validation.json`, `uncertainty_report.json`), enforced at both write time (all 5 producer scripts) and read time (`report/data_loader.py`). A malformed producer output now fails loudly at the source instead of silently degrading every downstream `.get(key, default)` consumer. 20 new tests (`test_schemas.py`).
- **`--research-mode` flag** (`prs.py run --full --research-mode`): separates audit/manuscript-only pipeline stages (scientific lock document, bilingual manuscripts, evidence pack, publication-readiness lock) from the ~50 stages that actually feed the interactive report, based on tracing every stage's output against what `report/data_loader.py` loads. Off by default — same complete report, faster `--full` run; the audit trail remains fully available with the flag or via `prs.py report`.
- **Palindromic-SNP strand-safety regression gate** (`test_allele_strand_consistency.py`): the existing forward-strand consistency check is structurally blind to A/T or C/G SNPs (a flipped palindromic pair looks identical either way). New test requires every palindromic row in the panel to be explicitly triaged (independently strand-verified, or tracked as a known gap) — a new palindromic SNP landing without review now fails CI instead of inheriting an invisible risk.
- **ClinVar body-system grouping** (`disease_taxonomy.py`): heuristic classifier groups pathogenic findings into 12 body systems (Cardiovascular, Oncology, Neurological, etc.) for the report; MedGen CUI persisted per variant with a direct link to the MedGen record.
- **47 of 56 traits now have a curated, evidence-cited actionable recommendation** (`data/trait_recommendations.json`), up from 10 — each verified against a real PubMed/NIH source before being written. The remaining 9 are documented exclusions (fragmentation duplicates or a genuine data limitation), not gaps.
- `prs.py` orchestrator unit tests (`test_prs_orchestrator_halt.py`, 7 tests) pin the new fail-loud behavior below.

### Changed
- **`prs.py run`'s "required stage" halting was completely non-functional — now actually halts.** `run_script(..., required=True)`/`require_output()` existed specifically to stop the pipeline on a critical genotype-processing failure, but only logged an error and returned a value every caller discarded; a failed Stage A/B/C previously let the pipeline keep running against missing/stale files. Fixed inside the two shared helpers (`sys.exit` on failure); verified with real data that a corrupt VCF now halts in ~2s instead of hanging past a 2-minute timeout.
- CI's dependency install list (`.github/workflows/test.yml`) had silently drifted out of sync with `requirements.txt` twice (missing `jinja2` since the 1.6 refactor, then `pydantic` since 2.2) — CI had been failing on every push since Aug 6 without being noticed. Fixed and annotated to prevent recurrence.
- `requirements.txt`/`requirements.lock`: added `pydantic>=2.0`.

### Fixed
- **Badge/color polarity bug**: favorable genetic findings (Morning chronotype, Cognitive function, Longevity & lifespan) were rendering as a red "HIGHER RISK" badge because the report's default badge logic assumes `risk_category=="high"` always means elevated health risk. Fixed via a `_polarity_inverted` trait list read by `risk_badge()`/`risk_color()`/`risk_bar()`, now covering all 3 confirmed cases after a full 53-trait audit of the remaining panel.
- **Non-deterministic `model_hash`** in one of the two `ANCESTRY_MODEL.json` producers (`pca_ancestry_classifier.py`): used Python's salted builtin `hash()`, varying across runs on identical input; switched to the same deterministic SHA-256 scheme the other producer already used.
- **`genotype_variance` silently written as a JSON string instead of a number** in `uncertainty_report.json`, caused by a numpy `float32` scalar (from cyvcf2's GQ arrays) falling through to a `json.dump(default=str)` fallback — fixed as a side effect of the pydantic schema migration, which coerces to the declared type at construction.
- **3 fragmentation-duplicate trait categories retired** ("Dopamine & cognition", "Longevity & aging", "Longevity") — each was an exact or near-exact re-measurement of a signal already captured under another trait_category (COMT/rs4680, and two APOE/FOXO3 cases respectively); retiring them avoids double-counting the same genetic signal under two names. Panel: 208 → 205 rows, 59 → 56 trait categories.
- **Dozens of fabricated or wrong PMIDs found and corrected across the SNP panel** over many curation rounds (citations that resolved to completely unrelated papers — plant morphogenesis, hair loss, starch crystallography, etc. — cited for vitamin D, dopamine, lipid, and inflammation SNPs). Each replacement was verified via a direct PubMed/PMC fetch before being written, never reconstructed from memory.
- **Multiple inverted `effect_allele`/`reference_allele`/`risk_genotype` bugs** found via the same citation-verification process (HFE C282Y, PPARGC1A, APOE, MC4R, TNF, IL1B, CRP, ABO/RHD gene mislabeling, and others) — each caused the affected SNP to score the wrong direction until corrected.
- **`Inflammation (CRP levels)`**: added a second, independent SNP (IL6R, previously depended on a single CRP-gene SNP) closing the last open case in the single-SNP-trait depth audit; also fixed an internally inconsistent `risk_genotype` left over from an earlier PMID correction.
- **3 remaining blank-PMID rows in `Vitamin D metabolism`** (VDR FokI/BsmI/ApaI) closed decisively: confirmed via the largest modern circulating-vitamin-D GWAS (Manousaki et al. 2020, 69 loci) that VDR is not a circulating-25(OH)D locus at all — not a citation-search gap. `evidence_level` corrected from an inconsistent "A" (implies genome-wide significance) to "D" (no supporting citation).

### Note
- **Full 53-trait polarity audit + a real scoring bug fixed (2026-08-08).** Auditing the rest of the panel for the same badge-polarity bug found 2 more cases (`Vitamin B12 status`, `Apolipoprotein A1`, added to `_polarity_inverted`) - and, along the way, a real *scoring* bug distinct from the badge issue: `effect_direction` was never read anywhere in the live PLINK-based scoring path (`bluegen/scoring.py::build_score_rows()`, plus an independent duplicate implementation in `scripts/utils/build_reference_distributions.py` with the same bug), so the 3 rows using `effect_direction='-'` (PCSK9/Apolipoprotein B, SLC17A1/Urate levels, CLOCK/Morning chronotype) were being added instead of subtracted. Fixed in both scoring paths; population reference distributions rebuilt against the full 1000G panel to stay consistent with the corrected convention - this also surfaced that the reference distributions had been stale since well before this session (never regenerated after the many SNP additions from the curation rounds above), so 25 of 56 traits' distributions updated, verified 1:1 against already-curated panel changes rather than a new regression.
- Also found and fixed while auditing: `rs671`/ALDH2 ("Alcohol flush reaction") had `effect_allele` set to the normal-function allele instead of the flush-causing one, plus a fabricated PMID - both corrected, and verified biologically (the rebuilt EAS reference distribution now shows real variance matching known ALDH2 deficiency prevalence in East Asians, vs. near-zero in the EUR reference, where the deficient allele is rare). The `Vitamin B12 status` recommendation text turned out to be **factually backwards** ("tendency toward lower B12, eat more B12-rich foods" when every cited source says the opposite) - rewritten to correctly frame it as a favorable finding, with the same clinical safety-net content preserved.
- 281 tests passing (up from ~90 before this stretch of work).
- Population Portability score dropped 62 → 32 as a *consequence* of the strand fix (below) landing in a benchmarking run for the first time — reflects the corrected panel measuring the real EUR-vs-other-ancestry PRS transfer gap more completely, not a regression.

### Fixed (earlier in this stretch, 2026-07-14 to 2026-08-04)
- **SNP panel curation** (`data/snp_database_annotated.csv`): 108 GRCh37 position errors corrected, 10 corrupted rsIDs resolved by reverse dbSNP lookup, 3 duplicate rsID/trait_category rows removed (each was double-counting a SNP's weight in the PRS sum), 7 new SNPs added (FUT2/B12, HFE H63D, ADORA2A, HNMT, COL1A1, CLOCK, ACE). Panel now 190 rows / 172 unique rsIDs / 59 trait categories.
- **34 strand-orientation mismatches**: `effect_allele`/`reference_allele` recorded in gene-relative (minus-strand) notation instead of GRCh37 forward-strand for minus-strand genes — PLINK `--score` silently drops any SNP that doesn't match its `.bim` alleles. 26 rows mechanically complemented, 6 rows had an independent effect_allele/risk_genotype mismatch fixed against primary literature, 2 rows removed (placeholder/wrong-gene data). User's own genotype scoring went from 23 to 25 scoring traits.
- **Calibration z-score/percentile inconsistency** (`population_calibrate_v2.py`): `percentile_population` was hardcoded to 50.0 in two branches even when a real z-score existed, silently forcing `risk_category="medium"` regardless of the true z. Also fixed a `norm.ppf` domain overflow at the 0/100 percentile boundary that floored the aggregate calibration stability score. Scientific integrity score: 74.1 (NEEDS_REVISION) → 86.7 (RESEARCH_GRADE).
- **GWAS consortium validation false-failures** (`25_gwas_consortium_validation.py`): sub-trait display aliases (e.g. "LDL cholesterol") were being looked up as literal `trait_category` values instead of mapping to the real combined category, reporting false FAILs. 4/17 → 16/17 consortiums passed.
- **PGS percentile bug** (`pgs_population_calibrate.py`): `'scipy' in dir()` with no arguments only checks local scope, so the condition was always False and percentile was hardcoded to 50 for all 30 PGS scores regardless of z-score. Fixed by moving the `scipy` import to module scope.
- **Pipeline version drift**: 9 files independently hardcoded their own `pipeline_version` (v1.0.0 to v10.0.0 across different files). Unified to a single constant (`utils/constants.py::PIPELINE_VERSION`, now `2.0.0`); CI now fails if any file drifts.

### Added (earlier in this stretch)
- `tests/test_allele_strand_consistency.py`, `tests/test_population_calibration.py`, `tests/test_version_consistency.py`, `tests/test_snp_positions.py` (regression coverage for all of the above).

## [1.2.1] — 2026-06-26

### Changed
- **Calibration Stability score recalculated** from actual per-trait data (mean_r²=0.995, was stale 0.887) — score corrected from 55.3 → 99.4
- **Population Portability now computed dynamically** from `reference_distributions.json` (56 traits × 5 populations, 2504 samples) instead of using hardcoded literature estimates — score improved from 54.5 → 71.5
- **Final Scientific Integrity Score: 80.2 → 88.0/100** (RESEARCH_GRADE, approaching PUBLICATION_READY threshold of 90)

### Fixed
- `calibration_validation.json`: aggregate `mean_slope` and `mean_r²` were stale (0.9982/0.8865) — corrected to 1.0015/0.9952 from per-trait values
- `26_population_portability_test.py`: eliminated hardcoded `reference_bias` dict; now computes shift/drift/instability dynamically using Cohen's d with pooled SD

## [1.2.0] — 2026-06-26

### Added
- **Radar chart (SVG)** in Executive Summary showing PRS profile across 9 traits with risk colors, z-scores, and concentric reference rings
- **Confidence indicators for PGS Catalog section**: reliability badges (✓ Reliable / ⚠ Unreliable), SNP coverage bars, risk bars, and summary bar
- **Clinical Actionability Summary** (new section): cross-references ClinVar pathogenic variants + PharmGKB drug findings + high-risk PRS traits into 3 subsections (High-Confidence Findings, PRS-Gene Convergence, Drug-Gene-PRS Intersections)
- **rsID → chr:pos mapper** (`pgs_rsid_mapper.py`): converts rsID-based PGS variants to chr:pos format for PLINK scoring compatibility. Integrated into `pgs_catalog_integration.py`
- **47 new unit tests** (9 radar chart + 1 legend + 2 portability + 8 rsID mapper + 14 uncertainty propagation + 13 existing extensions) — total 109 tests

### Changed
- **Uncertainty saturation fixed** in `14_uncertainty_propagation.py`:
  - Denominator changed from `abs(prs_val) * 0.5` (individual) to `population_sd * 2.0` (population-based)
  - E[G²] now computed from MAF under HWE instead of hardcoded 1.0
  - Evidence SE ratios reduced: A=0.10, B=0.20, C=0.35, D=0.55 (was A=0.20, B=0.33, C=0.50, D=0.75)
- **Calibration validation sign bug fixed** in `27_real_world_calibration.py`: `abs(z_pop)` removed from slope denominator
- **PGS Catalog section** now loads and displays coverage data from `pgs_results.csv`
- **Executive Summary** restructured with radar chart between risk counts and ancestry/integrity cards

### Fixed
- Calibration slope sign bug for below-median traits (Bitter taste, Lactose intolerance)
- 8/9 traits no longer have saturated uncertainty scores

## [1.1.1] — 2026-06-25

### Added
- **Per-trait confidence scores (0-100% with star ratings)** in the comprehensive PRS report, combining SNP coverage, calibration quality, evidence level, and uncertainty
- **Calibration quality flags** (GOOD/FAIR/POOR/INVERTED) per trait based on R² and slope from `calibration_validation.json`
- **Trust tier classification** (T1 High Trust / T2 Moderate / T3 Low Trust) for each PRS trait
- **Multi-layer uncertainty mini-bars** replacing the single "Uncert." column with genotype/ancestry/effect variance decomposition
- **SNP coverage visual bars** replacing plain text ratios
- **Population portability warning banner** at the top of the PRS results section
- **Per-trait limitation badges** (Low SNPs, Inverted cal., High uncert.) in the PRS table
- **Trait-specific limitation notes** in the expanded limitations section
- **Executive summary confidence overview**: Avg Confidence, High Trust (T1), Low Trust (T3) counts, Strongest/Weakest finding
- **Trust tier legend** above the PRS results table
- **`tests/test_comprehensive_report.py`:** 47 new unit tests for confidence helper functions

### Changed
- **`comprehensive_report.py`:** PRS table expanded from 7 to 11 columns with confidence indicators
- **`comprehensive_report.py`:** Now loads `benchmark/calibration_validation.json` for per-trait calibration quality data
- **`comprehensive_report.py`:** Evidence levels loaded from SNP database for confidence scoring
- **`comprehensive_report.py`:** Executive summary restructured with 2-row grid layout
- **`comprehensive_report.py`:** Limitations section now includes per-trait breakdown

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
- **Pipeline version:** Normalized to `v1.1.0` across all output files (was v6/7/9/10 in different JSONs)
- **Tool detection:** `16_reproducibility_engine.py` now uses `tool_detection.detect_all()` — only reports primary PLINK (v1.9), not plink2
- **`README.md`:** Added BlueGen logo (`assets/blugen-logo.png`)
- **`LICENSE`:** Added MIT license file
- **`scripts/setup/download_archaic_reference.py`:** Curated 133-SNP archaic introgression panel with population baselines
- **`scripts/setup/download_aadr_reference.py`:** AADR 1240K panel download + Eigenstrat→PLINK conversion (Altai/Vindija/Chagyrskaya + Denisova)
- **`scripts/setup/download_vindija_reference.py`:** Vindija Neanderthal genome VCF download (chr1-22, hg19, ~44 GB) with parallel support
- **`scripts/clinical/ancestry_deep.py`:** Vindija direct comparison (gold standard) + 1000G African/European frequency filter (661 AFR + 503 EUR)
- **`archive_upload.py`:** Added `bluegen-archaic-reference` item (188 MB PLINK binaries on archive.org)
- **`scripts/clinical/ancestry_deep.py`:** AADR direct archaic genome comparison + 133-SNP panel fallback; batch bcftools query for fast VCF extraction

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
