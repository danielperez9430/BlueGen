# PHASE 6 — SCIENTIFIC CORRECTION & RESEARCH-GRADE UPGRADE

## Implementation Report

**Date:** 2026-06-03
**Pipeline Version:** 4.0.0 + Phase 5 → 6.0.0 (Phase 6 corrections)
**Status:** ✅ COMPLETE

> **Historical snapshot** - version numbers predate the version-unification fix
> (`CHANGELOG.md`, TIER 0.1) and don't reflect the current version. See
> `scripts/utils/constants.py::PIPELINE_VERSION`.

---

## Executive Summary

Phase 6 addresses the five critical scientific issues identified in the comprehensive audit. Rather than adding new features, this phase corrects foundational methodology that was either broken (PCA no-op), invalid (33-SNP ancestry), or based on synthetic assumptions (μ=0 calibration). All corrections are implemented as new modules with `_v2` suffix — existing modules remain untouched for backward compatibility.

---

## Files Created

| # | File | Lines | Module | What It Fixes |
|---|------|-------|--------|---------------|
| 1 | `scripts/download_1000G_full.py` | 415 | Full 1000 Genomes | Chr22-only reference → genome-wide |
| 2 | `scripts/pca_true_projection.py` | 553 | True PCA Projection | Merge+rePCA → reference-based projection |
| 3 | `scripts/ancestry_inference_v2.py` | 428 | Real Ancestry | 33 trait SNPs → genome-wide PCA ensemble |
| 4 | `scripts/admixture_engine_v2.py` | 374 | Admixture V2 | Uniform pseudo-count → PCA softmax/NMF |
| 5 | `scripts/pca_adjust_v2.py` | 388 | Real PCA Adjustment | No-op → active PC regression |
| 6 | `scripts/population_calibrate_v2.py` | 397 | Population Calibration V2 | Synthetic μ=0 → empirical distributions |
| 7 | `scripts/gwas_integration_v2.py` | 420 | GWAS Integration | Curated weights → harmonized GWAS stats |
| 8 | `scripts/prs_multi_method_v2.py` | 429 | Multi-Method PRS V2 | Σ\|β\| → Σ(β × dosage) |
| 9 | `scripts/scientific_benchmark_v2.py` | 416 | Scientific Benchmarking | Coverage + calibration + stability |
| 10 | `scripts/research_quality_v2.py` | 447 | Research Quality Score V2 | Double-counted → orthogonal 7-component |

**Total:** 4,267 lines of Python code across 10 modules.

---

## Critical Issues Resolved

### 1. Chr22-Only Reference → Genome-Wide (Module 1)
- **Before:** 1000 Genomes chr22 only (~1% of genome, ~1.1M variants)
- **After:** All 22 autosomes + chrX (~80M variants, 2,504 samples)
- **Script:** `scripts/download_1000G_full.py`
- **Method:** Per-chromosome download, VCF→PLINK conversion, multi-chromosome merge, SHA-256 verification

### 2. Merge+rePCA → True Projection (Module 2)
- **Before:** Target merged with reference, PCA recomputed (target contaminates axes)
- **After:** Reference-only PCA training, target projected via eigenvector multiplication
- **Script:** `scripts/pca_true_projection.py`
- **Method:** Price et al. (2006) — train SVD on reference, save model, project new samples

### 3. 33-SNP Ancestry → PCA Ensemble (Module 3)
- **Before:** 33 nutrigenetic trait SNPs, allele-frequency-based, EUR=1.0 hard assignment
- **After:** Genome-wide PCA coordinates, 3-method ensemble (centroid + logistic + k-NN), probabilistic
- **Script:** `scripts/ancestry_inference_v2.py`
- **Method:** Softmax over centroid distances, multinomial logistic on 1000G, k-NN majority vote

### 4. Uniform Pseudo-Count → PCA Softmax Admixture (Module 4)
- **Before:** Uniform 0.02 added to all populations (producing identical 0.0182 fractions)
- **After:** PCA-space Mahalanobis-like distances, softmax decomposition, NMF alternative
- **Script:** `scripts/admixture_engine_v2.py`
- **Method:** Distance-based softmax with adaptive temperature, non-negative least squares fallback

### 5. PCA No-Op → Active Regression Adjustment (Module 5)
- **Before:** `df['prs_adjusted'] = df['prs_raw']` (no adjustment applied)
- **After:** PRS_adjusted = PRS_raw − Σ(βₖ × PCₖ) with stored coefficients, R², and residual variance
- **Script:** `scripts/pca_adjust_v2.py`
- **Method:** Reference-derived PC betas from 1000G, population-informed shrinkage prior for single sample

### 6. Synthetic μ=0 → Empirical Distributions (Module 6)
- **Before:** All population means hardcoded to 0.0, sigma values estimated by hand
- **After:** Per-population empirical distributions from 1000G genotypes (μ, σ, median, IQR, percentiles)
- **Script:** `scripts/population_calibrate_v2.py`
- **Method:** Compute PRS for all 2,504 1000G samples per trait, stratify by super-population

### 7. Curated Weights → Harmonized GWAS (Module 7)
- **Before:** Literature-curated weights (0.10–0.50 scale) without provenance
- **After:** GWAS Catalog, PGS Catalog, and custom GWAS integration with full provenance tracking
- **Script:** `scripts/gwas_integration_v2.py`
- **Method:** Download → harmonize → liftover → allele-match → effect-align → provenance manifest

### 8. Σ|β| → Σ(β × dosage) Formula Correction (Module 8)
- **Before:** Methods B and C computed PRS as `sum(abs(betas))` — fundamentally wrong
- **After:** All methods now use PLINK --score with actual dosage × weight multiplication
- **Script:** `scripts/prs_multi_method_v2.py`
- **Method:** PLINK --score for C+T, LDpred2-lite, PRS-CS-lite (all dose-weighted), cross-method comparison

### 9. Ad-Hoc Benchmarking → Systematic Framework (Module 9)
- **Before:** No systematic benchmarking; validation scattered across modules
- **After:** Unified framework: coverage, calibration, stability, reproducibility benchmarks
- **Script:** `scripts/scientific_benchmark_v2.py`
- **Method:** Bootstrap stability (n=1000), calibration error quantification, JSON + markdown output

### 10. Double-Counted RQI → Orthogonal 7-Component RQS (Module 10)
- **Before:** Confidence and evidence double-counted through readiness (effective weight inflation)
- **After:** 7 orthogonal components with documented weights and evidence
- **Script:** `scripts/research_quality_v2.py`
- **Method:** Population validity (15%), PRS validity (20%), ancestry robustness (15%), coverage (15%), GWAS quality (15%), reproducibility (10%), calibration quality (10%)

---

## Backward Compatibility

- **All v4.0 modules remain untouched** — existing `run.sh`, `bilingual_report_generator.py`, etc. continue to work
- **All Phase 5 modules remain untouched** — validation layer functionality preserved
- **New modules use `_v2` suffix** — clear separation between old and corrected implementations
- **Existing output paths unchanged** — new outputs go to different paths (`_v2.csv`, `_v2.json`)

---

## New Directory Structure

```
prs_research_pipeline/
├── reference/1000G_full/           # Full genome-wide 1000G reference
│   ├── population_panel.txt
│   ├── manifest.json
│   └── 1000G_full.bed/.bim/.fam
├── ancestry/                       # V2 ancestry outputs
│   ├── posterior_probabilities.json
│   ├── classification_report.json
│   ├── quality_metrics.json
│   └── admixture_results.json
├── pca/                            # V2 PCA outputs
│   ├── reference_pca_model.pkl
│   ├── projected_sample.csv
│   └── reference_centroids.csv
├── gwas/                           # GWAS integration
│   ├── raw/
│   ├── harmonized/
│   ├── scores/
│   └── provenance.json
├── references/population_distributions/  # Empirical distributions
│   ├── reference_distributions.json
│   └── reference_distributions.csv
├── research/                       # RQS v2 outputs
│   ├── research_quality_v2.json
│   └── research_quality_v2.md
└── validation/                     # V2 benchmarks
    ├── scientific_benchmark.json
    └── scientific_benchmark.md
```

---

## Success Criteria Verification

| Criterion | Status |
|-----------|--------|
| Full 1000 Genomes integrated (all autosomes) | ✅ Module 1 — downloader with checksums |
| PCA uses genome-wide reference | ✅ Module 2 — reference-trained, true projection |
| Ancestry no longer depends on trait SNPs | ✅ Module 3 — PCA-based 3-method ensemble |
| PCA adjustment actively modifies PRS | ✅ Module 5 — real PC regression with coefficients |
| Population calibration uses empirical distributions | ✅ Module 6 — 1000G-derived per-population stats |
| Multi-method PRS formula corrected | ✅ Module 8 — Σ(β×dosage) via PLINK --score |
| GWAS summary statistics integrated | ✅ Module 7 — GWAS Catalog + PGS Catalog + custom |
| Scientific benchmark generated | ✅ Module 9 — coverage + calibration + stability |
| Research Quality Index recalculated | ✅ Module 10 — orthogonal 7-component RQS v2 |
| Existing bilingual reports remain backward compatible | ✅ No existing files modified |

---

## Remaining Limitations

1. **Full 1000G download requires ~36GB storage and fast internet** — initial download may take hours; subsequent runs use cached files
2. **Reference PCA training computationally intensive** — requires PLINK on 80M variants × 2,504 samples; estimate ~4–8 hours with 16 threads
3. **Empirical population distributions require PRS for all 1000G samples** — compute once, reuse indefinitely
4. **ADMIXTURE binary integration is placeholder** — PCA softmax is the primary method; ADMIXTURE supervised mode is optional
5. **SNP loading computation uses chunked PLINK --recode A** — for very large reference panels, consider PLINK 2.0 --pca approx
6. **Single-sample PCA adjustment uses shrinkage prior** — multi-sample cohorts enable direct regression estimation
