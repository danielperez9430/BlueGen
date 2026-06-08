# PRS RESEARCH PLATFORM — COMPREHENSIVE SCIENTIFIC & ENGINEERING AUDIT

**Date:** 2026-06-03
**Platform Version:** 4.0.0 + Phase 5 validation layer
**Auditor:** Automated multi-agent + direct analysis
**Scope:** All 25 scripts, 10 shell pipelines, 11 Python modules, configuration, outputs

---

## EXECUTIVE SUMMARY

The PRS Research Platform is a **moderately mature research tool** with strong conceptual design but significant implementation gaps in several critical areas. The platform correctly implements PLINK-based genotype processing, has a well-structured bilingual reporting system, and adds a thoughtful multi-dimensional validation layer. However, the ancestry inference system operates on insufficient data (33 SNPs, chr22 only), the PCA adjustment step is non-functional, population calibration uses zero-means defeating its purpose, and uncertainty estimates are driven to saturation for most traits. The Phase 5 validation modules add architectural sophistication but inherit the limitations of upstream components.

**Overall Assessment:** The platform is suitable for **exploratory research and methodology development** but has multiple blockers preventing publication-grade or clinical-use readiness.

---

## GLOBAL SCORES

| Dimension | Score | Grade |
|-----------|-------|-------|
| Architecture Design | 7.2/10 | B |
| Code Quality | 6.8/10 | B- |
| Scientific Validity | 5.4/10 | C+ |
| Reproducibility | 7.0/10 | B |
| Computational Performance | 7.5/10 | B |
| External Dependencies | 6.5/10 | B- |
| PRS Methodology | 4.8/10 | C |
| Ancestry Inference | 3.2/10 | D+ |
| Bilingual Reporting | 8.2/10 | A- |
| Benchmark Quality | 5.0/10 | C |

---

## MODULE-BY-MODULE AUDIT

### CORE PIPELINE (Shell Scripts)

#### 01_vcf_to_plink.sh — VCF→PLINK Conversion
| Metric | Score | Notes |
|--------|-------|-------|
| Maturity | 7/10 | Standard PLINK VCF import; handles bgzip correctly |
| Scientific Confidence | 8/10 | Correctly converts VCF to PLINK binary |
| Completeness | 7/10 | Missing explicit handling of multi-allelic sites and sex chromosomes |
| Risk | LOW | |

**Strengths:** Proper bgzip handling, bcftools index verification, preserves sample IDs.
**Weaknesses:** Does not filter out structural variants; no explicit chrX/Y ploidy handling.

#### 02_quality_control.sh — Quality Control
| Metric | Score | Notes |
|--------|-------|-------|
| Maturity | 8/10 | Standard QC filters implemented correctly |
| Scientific Confidence | 8/10 | geno/mind/maf/hwe thresholds are standard |
| Completeness | 7/10 | Missing heterozygosity check, sex check, relatedness |
| Risk | LOW | |

**Strengths:** Configurable thresholds, generates QC statistics files.
**Weaknesses:** No heterozygosity outlier detection; no sample relatedness check (IBD); HWE midp not specified.

#### 03_ld_pruning.sh — Standard LD Pruning
| Metric | Score | Notes |
|--------|-------|-------|
| Maturity | 8/10 | Correct indep-pairwise implementation |
| Scientific Confidence | 7/10 | Standard parameters (50 5 0.2) |
| Completeness | 7/10 | Single-population only |
| Risk | MEDIUM | LD patterns vary by ancestry |

**Strengths:** Clean implementation, logs retention rates.
**Weaknesses:** Single-population pruning may retain correlated SNPs for non-EUR individuals.

#### 03_ld_ancestry_prune.sh — Ancestry-Matched LD Pruning
| Metric | Score | Notes |
|--------|-------|-------|
| Maturity | 7/10 | Well-designed per-population pruning |
| Scientific Confidence | 6/10 | Conservative intersection is methodologically sound |
| Completeness | 5/10 | **Uses chr22-only 1000G reference** |
| Risk | HIGH | Chr22-only LD pruning is not genome-wide |

**Strengths:** Elegant union/intersection logic, per-population keep files, good Python inline processing.
**Weaknesses:** **Critical: chr22-only reference limits LD pruning to a single chromosome (~1% of genome).** LD blocks on chr22 do not represent genome-wide LD patterns.

#### 04_pca_1000G.sh — 1000 Genomes PCA
| Metric | Score | Notes |
|--------|-------|-------|
| Maturity | 6/10 | Correct projection methodology |
| Scientific Confidence | 3/10 | **chr22-only PCA cannot capture genome-wide ancestry** |
| Completeness | 5/10 | Works but on fundamentally insufficient data |
| Risk | **CRITICAL** | PCA on 1 chromosome ≠ genome-wide ancestry |

**Strengths:** Correct projection logic (train on ref, project target), handles strand flips.
**Weaknesses:** **CRITICAL: Chromosome 22 contains ~1% of the genome. PCA computed on chr22 alone cannot reliably distinguish population clusters.** Published methods use genome-wide LD-pruned SNPs (typically 50,000–150,000 SNPs across all autosomes). This is the single most severe scientific limitation in the platform.

#### run.sh — Master Orchestrator
| Metric | Score | Notes |
|--------|-------|-------|
| Maturity | 7/10 | Well-organized 17-stage pipeline |
| Scientific Confidence | 5/10 | Pipeline flow is correct but individual stages have issues |
| Completeness | 7/10 | All stages connected; graceful degradation on failures |
| Risk | MEDIUM | |

**Strengths:** Comprehensive orchestration, clear logging, auto-detection of tools.
**Weaknesses:** **The PCA-adjusted PRS step (Stage G) is a no-op:** `df['prs_adjusted'] = df['prs_raw']`. No actual regression adjustment is performed.

---

### CORE PRS COMPUTATION (Python)

#### 06_prs_compute.py — PRS Computation
| Metric | Score | Notes |
|--------|-------|-------|
| Maturity | 8/10 | Robust PLINK and VCF-based computation |
| Scientific Confidence | 5/10 | Curated weights not from formal GWAS summary stats |
| Completeness | 7/10 | Good fallback logic for VCF without rsIDs |
| Risk | MEDIUM | Weights from literature, not uniform GWAS |

**Strengths:** Position-based VCF extraction is clever for DeepVariant data; correctly handles missing genotypes; dosage computation is accurate.
**Weaknesses:** **Effect weights are curated literature values, not harmonized GWAS summary statistics.** Standard PRS uses genome-wide significant SNPs from large GWAS with uniform QC. Weights of 0.10–0.50 are not on a standardized scale; the PRS is essentially a weighted count of risk alleles.

#### 11_population_calibrate.py — Population Calibration
| Metric | Score | Notes |
|--------|-------|-------|
| Maturity | 6/10 | KDE + bootstrap approach is clever |
| Scientific Confidence | 3/10 | **Reference population μ=0 for all traits** |
| Completeness | 5/10 | Bootstrap CI computation is valid |
| Risk | **CRITICAL** | Calibration is effectively void |

**Strengths:** Bootstrap CIs, KDE-based percentile estimation, handles edge cases.
**Weaknesses:** **CRITICAL: All `population_mu` values are 0.0.** Population calibration requires per-population reference means and standard deviations. Setting μ=0 for all populations means the "calibration" is just dividing by sigma, not actually calibrating against population-specific distributions. **This defeats the entire purpose of population calibration.**

#### 05_harmonize_gwas.py — GWAS Harmonization
| Metric | Score | Notes |
|--------|-------|-------|
| Maturity | 6/10 | Handles strand flips and allele alignment |
| Scientific Confidence | 6/10 | Correct harmonization logic |
| Completeness | 5/10 | Requires external GWAS files not typically available |
| Risk | HIGH | Limited testing with real diverse GWAS formats |

**Strengths:** Palindromic SNP handling, strand flip detection, effect direction alignment.
**Weaknesses:** GWAS file format requirements are rigid; the module is rarely exercised in practice since most runs use the curated database fallback.

#### prs_multi_method.py — Multi-Method PRS
| Metric | Score | Notes |
|--------|-------|-------|
| Maturity | 4/10 | Framework exists but only curated method runs |
| Scientific Confidence | 4/10 | C+T/LDpred2/PRS-CS are approximations, not actual implementations |
| Completeness | 3/10 | **Only curated method produces output in practice** |
| Risk | HIGH | Methods don't run without GWAS files |

**Strengths:** Clean multi-method architecture, good comparison framework.
**Weaknesses:** **C+T, LDpred2-lite, and PRS-CS-lite produce NO output** because they require GWAS summary statistics files that the platform does not provide. The "LDpred2-lite" is a single shrinkage factor (`h²/(h²+M/N)`) applied uniformly, not actual LDpred2 with LD matrices. The "PRS-CS-lite" is soft-thresholding, not the actual PRS-CS continuous shrinkage model.

---

### ANCESTRY & POPULATION GENETICS

#### ancestry_inference.py
| Metric | Score | Notes |
|--------|-------|-------|
| Maturity | 5/10 | Three methods implemented (logistic, centroid, Mahalanobis) |
| Scientific Confidence | 2/10 | **33 nutrigenetic SNPs cannot classify ancestry** |
| Completeness | 4/10 | Output structure is good |
| Risk | **CRITICAL** | Fundamentally insufficient data |

**Strengths:** Multiple classification methods, probability calibration, well-structured output.
**Weaknesses:** **CRITICAL: The current output (EUR=1.0) was derived from 33 trait-associated SNPs, not from PCA coordinates.** The ancestry_inference.json states: "Inferred from 33 nutrigenetic SNPs against 1000G Phase 3 population allele frequencies." This is fundamentally invalid — trait-associated SNPs have allele frequencies that vary by trait selection, not by neutral ancestry. Published ancestry inference uses thousands of LD-pruned, genome-wide SNPs that are approximately neutral. **33 SNPs is orders of magnitude below what's needed for reliable ancestry classification.**

#### admixture_engine.py
| Metric | Score | Notes |
|--------|-------|-------|
| Maturity | 4/10 | PCA softmax approach is reasonable in theory |
| Scientific Confidence | 2/10 | **Non-EUR fractions are uniform (all 0.0182)** |
| Completeness | 3/10 | Produces output but it's not meaningful |
| Risk | **CRITICAL** | Admixture fractions are artifacts |

**Strengths:** PCA softmax is a valid admixture estimation approach.
**Weaknesses:** **CRITICAL: AFR=EAS=SAS=AMR all show exactly 0.0182 admixture fraction.** This is a mathematical artifact — when PCA places a sample in the EUR cluster, the softmax distributes residual probability uniformly across all other populations. These are not genuine admixture signals; they are a byproduct of the softmax over a degenerate PCA space (chr22 only → insufficient separation).

#### 13_gwas_ld_consistency_check.py
| Metric | Score | Notes |
|--------|-------|-------|
| Maturity | 6/10 | Good conceptual framework |
| Scientific Confidence | 4/10 | **Always passes because all traits hardcoded as EUR** |
| Completeness | 5/10 | Detects mismatches but never finds any |
| Risk | MEDIUM | False sense of security |

**Strengths:** Well-structured risk flag system, per-trait checking.
**Weaknesses:** **The curated database hardcodes all GWAS populations as EUR.** Since the sample is also classified as EUR, the consistency check always returns "ok" with 0.0 confidence downgrade. This is not validation — it's self-consistency of a hardcoded assumption.

#### 14_uncertainty_propagation.py
| Metric | Score | Notes |
|--------|-------|-------|
| Maturity | 7/10 | Good three-layer decomposition framework |
| Scientific Confidence | 5/10 | Propagation formula is correct but inputs are problematic |
| Completeness | 6/10 | All three layers implemented |
| Risk | HIGH | Uncertainty scores saturated at 1.0 |

**Strengths:** Correct variance propagation formula (Var(PRS)=Σβ²Var(G)+ΣG²Var(β)), good decomposition, GQ-to-error conversion.
**Weaknesses:** **7/10 traits have uncertainty_score=1.0 (saturated).** This means `PRS_SE > 0.5 × |PRS|` for most traits. The uncertainty is driven by the effect variance component (which dominates in 8/10 traits). This suggests the effect size uncertainty estimates from evidence levels are overly conservative.

---

### BILINGUAL REPORTING

#### bilingual_interpretation.py
| Metric | Score | Notes |
|--------|-------|-------|
| Maturity | 9/10 | Excellent curated bilingual knowledge base |
| Scientific Confidence | 8/10 | Scientifically accurate trait descriptions |
| Completeness | 9/10 | Anti-circularity safeguard is sophisticated |
| Risk | LOW | |

**Strengths:** Curated EN/ES content (not machine-translated), anti-circularity check preventing interpretation→PRS feedback loops, comprehensive trait coverage, evidence citations with PMIDs.
**Weaknesses:** Some trait descriptions are simplified for readability; dietary context is general not personalized.

#### bilingual_report_generator.py / 10_report_generate.py
| Metric | Score | Notes |
|--------|-------|-------|
| Maturity | 8/10 | Well-structured HTML generation |
| Scientific Confidence | 7/10 | Accurate representation of results |
| Completeness | 8/10 | EN/ES HTML + PDF, uncertainty and consistency sections |
| Risk | LOW | |

**Strengths:** Professional CSS, print stylesheets, bilingual parity, inline templates for robustness.
**Weaknesses:** Some duplication between the two report generators; WeasyPrint dependency is fragile on macOS.

---

### PHASE 5 VALIDATION MODULES

#### pgs_catalog_import.py
| Metric | Score | Notes |
|--------|-------|-------|
| Maturity | 6/10 | Good API + FTP download design |
| Scientific Confidence | 7/10 | Correct metadata extraction |
| Completeness | 7/10 | Caching, harmonization, build detection |
| Risk | MEDIUM | Network-dependent; liftOver optional |

**Strengths:** Checksum-based caching, build detection, flexible format parsing, REST API integration.
**Weaknesses:** No automated FTP directory listing; requires PGS IDs to be known in advance. pyliftover build conversion is optional and usually unavailable.

#### pgs_score_runner.py
| Metric | Score | Notes |
|--------|-------|-------|
| Maturity | 7/10 | Clean PLINK score wrapper |
| Scientific Confidence | 7/10 | Correct allele counting and coverage logic |
| Completeness | 7/10 | Coverage warnings, allele mismatch counting |
| Risk | MEDIUM | Requires PGS scores to be downloaded first |

**Strengths:** Strand flip detection, coverage warning at 80%, normalized scores.
**Weaknesses:** Score harmonization to PLINK format is basic (drops non-ACGT alleles); no dosage support beyond hard-called genotypes.

#### concordance_analysis.py
| Metric | Score | Notes |
|--------|-------|-------|
| Maturity | 6/10 | Good multi-metric approach |
| Scientific Confidence | 5/10 | **Single-sample "correlations" are pseudo-metrics** |
| Completeness | 6/10 | HTML report with proper thresholds |
| Risk | HIGH | Single-sample concordance is not statistically meaningful |

**Strengths:** Multi-metric comparison (Pearson, Spearman, Kendall, percentile, direction), proper PASS/WARNING/FAIL thresholds.
**Weaknesses:** **For single samples, Pearson r is computed as a pseudo-correlation** `1 - |z1-z2|/(2*max(|z1|,|z2|,1))`. This is NOT a correlation — it's a distance-based heuristic. True concordance requires multi-sample cohorts where score distributions can be compared.

#### evidence_scoring.py
| Metric | Score | Notes |
|--------|-------|-------|
| Maturity | 6/10 | Reasonable Bayesian-like framework |
| Scientific Confidence | 5/10 | Weights are arbitrary; evidence levels are coarse |
| Completeness | 6/10 | All 6 factors implemented |
| Risk | MEDIUM | Score calibration not empirically validated |

**Strengths:** Multi-factor evidence synthesis, ancestry match scoring, proper evidence level mapping.
**Weaknesses:** Component weights (25/20/20/15/10/10) are not empirically derived. Evidence level → p-value mapping is approximate. Replication count scoring is categorical, not continuous.

#### clinical_readiness.py
| Metric | Score | Notes |
|--------|-------|-------|
| Maturity | 6/10 | Sound dimensionality for the assessment |
| Scientific Confidence | 4/10 | **Readiness score inherits upstream issues** |
| Completeness | 6/10 | Good flags and recommendations |
| Risk | HIGH | Index is only as good as its inputs |

**Strengths:** Five well-chosen dimensions, bilingual recommendations, critical flag detection.
**Weaknesses:** **Readiness scores are dominated by upstream problems** (fake ancestry, saturated uncertainty). The index provides a false sense of precision — a score of 75 vs 73 is not meaningfully different given input uncertainties.

#### RQI (research_quality_index.py) / Limitations (limitations_engine.py) / Dashboard (validation_dashboard.py)
| Metric | Score | Notes |
|--------|-------|-------|
| Maturity | 7/10 | Good composite architecture |
| Scientific Confidence | 5/10 | RQI weights are arbitrary; limitations auto-detection is threshold-based |
| Completeness | 7/10 | Interactive dashboard is well-designed |
| Risk | MEDIUM | |

**Strengths:** Interactive dashboard with JS filtering, curated bilingual limitation text, composable RQI weights.
**Weaknesses:** RQI component weights (15/25/25/15/20) lack empirical justification. Limitations engine uses simple thresholds rather than statistical tests.

---

## CROSS-CUTTING FINDINGS

### CRITICAL ISSUES (Blockers)

1. **🔴 Chr22-only reference for PCA and LD** — All population genetics (PCA, ancestry inference, LD pruning) relies on chromosome 22 only (~1% of genome). Genome-wide ancestry patterns cannot be reconstructed from a single chromosome. This affects: PCA projection, ancestry classification, admixture estimation, population calibration, LD pruning reference.

2. **🔴 Ancestry inference from 33 trait SNPs** — Current output used 33 trait-associated SNPs for ancestry classification. Published methods require thousands of neutral, LD-pruned, genome-wide SNPs. Trait-associated SNPs have distorted allele frequencies that bias ancestry inference.

3. **🔴 Population calibration μ=0** — All population means are set to zero, defeating the purpose of population-specific calibration. The platform computes z = PRS/σ_pop, which is just scaling, not calibration.

4. **🔴 PCA adjustment is a no-op** — Stage G in run.sh copies raw scores directly: `df['prs_adjusted'] = df['prs_raw']`. No regression against PCs is performed despite the config specifying `pca_adjustment.enabled: true`.

5. **🔴 Multi-method PRS doesn't run** — C+T, LDpred2-lite, and PRS-CS-lite engines produce no output without GWAS summary statistics. The prs_all_methods.csv shows only "curated" method results. The LDpred2-lite and PRS-CS-lite implementations are approximations, not faithful to the published methods.

### HIGH-SEVERITY ISSUES

6. **🟠 Uncertainty scores saturated** — 7/10 traits show uncertainty_score=1.0, meaning PRS_SE exceeds 50% of PRS value. Effect variance dominates (50-86% of total variance), suggesting evidence-level-based SE estimation is overly conservative.

7. **🟠 Admixture fractions are artifacts** — AFR, EAS, SAS, and AMR all show exactly 0.0182. This is a softmax artifact, not genuine admixture.

8. **🟠 GWAS-LD consistency always passes** — All traits are hardcoded as EUR→EUR matches. The check provides no real validation.

9. **🟠 Concordance pseudo-correlation** — Single-sample Pearson r is computed as a distance heuristic. Meaningful concordance requires multi-sample comparisons.

10. **🟠 Confidence score misleading** — The 62.5/100 "Moderate Confidence" score is misleading because several sub-components are derived from degraded inputs (e.g., ancestry_gwas_match=20/20 from a check that always passes).

### STRENGTHS

1. ✅ **Excellent bilingual reporting** — Curated EN/ES content, professional HTML/CSS, anti-circularity safeguards, comprehensive trait knowledge base with PMIDs.
2. ✅ **Well-structured architecture** — Clean separation of concerns, modular scripts, consistent patterns (dataclasses, agents, argparse).
3. ✅ **Proper PLINK integration** — Correct VCF→PLINK→QC→LD pipeline with subprocess management.
4. ✅ **Comprehensive validation framework design** — Phase 5 adds thoughtful multi-dimensional quality assessment.
5. ✅ **Reproducibility infrastructure** — Run manifests with SHA-256 hashes, software version tracking, deterministic mode.
6. ✅ **VCF position-based genotype extraction** — Clever solution for DeepVariant VCFs lacking rsID annotations.
7. ✅ **Interactive dashboard** — Well-designed HTML with JS filtering and sorting.
8. ✅ **Scientific limitation auto-detection** — Curated bilingual limitation text with severity classification.

---

## PUBLICATION BLOCKERS

| # | Blocker | Severity | Fix Required |
|---|---------|----------|--------------|
| 1 | Chr22-only PCA → ancestry inference invalid | **Critical** | Download full 1000G (all autosomes); re-run PCA |
| 2 | 33-SNP ancestry classification | **Critical** | Use genome-wide LD-pruned SNPs from PCA projection |
| 3 | PCA adjustment is a no-op | **Critical** | Implement actual PC regression (already coded in 07_prs_adjust.py) |
| 4 | Population calibration μ=0 | **Critical** | Compute per-population reference distributions from 1000G |
| 5 | Multi-method PRS non-functional | **High** | Integrate actual GWAS summary statistics |
| 6 | No multi-sample validation | **High** | Test on at least 100–1000 diverse samples |
| 7 | Concordance uses pseudo-metrics | **High** | Multi-sample cohort needed for real correlations |
| 8 | GWAS effect sizes not harmonized | **High** | Use PGS Catalog or formal GWAS summary stats |

## CLINICAL-USE BLOCKERS

All publication blockers apply, plus:

| # | Blocker | Severity | Note |
|---|---------|----------|------|
| 9 | No prospective clinical validation | **Critical** | No outcome data linked to PRS scores |
| 10 | No regulatory framework compliance | **Critical** | Not FDA/EMA/CLIA compliant |
| 11 | Single-sample design | **Critical** | Clinical PRS requires population reference norms |
| 12 | Nutrigentic, not disease-risk PRS | **High** | Traits are metabolic/intermediary, not clinical endpoints |
| 13 | No clinical decision support integration | **High** | Reports not suitable for clinical workflow |

---

## TOP 10 PRIORITIES

Ordered by: (Scientific Impact × Implementation Feasibility) → Priority Score

| Rank | Priority | Impact | Effort | Module |
|------|----------|--------|--------|--------|
| 1 | Download full 1000G reference (all autosomes) | 10/10 | 6/10 | PCA, ancestry, LD |
| 2 | Implement actual PCA regression adjustment | 9/10 | 3/10 | 07_prs_adjust.py |
| 3 | Compute real population reference distributions | 9/10 | 5/10 | 11_population_calibrate.py |
| 4 | Integrate actual GWAS summary statistics | 8/10 | 7/10 | 05_harmonize_gwas.py |
| 5 | Add multi-sample cohort support | 8/10 | 8/10 | pipeline core |
| 6 | Fix multi-method PRS to run with available data | 7/10 | 4/10 | prs_multi_method.py |
| 7 | Calibrate uncertainty estimates against empirical data | 7/10 | 5/10 | 14_uncertainty_propagation.py |
| 8 | Add genome-wide PGS Catalog scores | 6/10 | 4/10 | pgs_catalog_import.py |
| 9 | Fix admixture engine with full PCA space | 6/10 | 3/10 | admixture_engine.py |
| 10 | Validate evidence scoring weights empirically | 5/10 | 6/10 | evidence_scoring.py |
