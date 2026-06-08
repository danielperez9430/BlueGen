# BlueGen — Future Roadmap

**Date:** 2026-06-09 | **Current:** v1.0.0

---

## Tier 1 — High Impact / Short Term

### 1. Full PGS Calibration
| | |
|---|---|
| **Impact** | 9/10 — Complex disease risk with percentiles |
| **Effort** | 3/10 — 1 day |
| **Plan** | Calibrate all 54 PGS scores against 1000G EUR reference. Query PGS Catalog API for published population means/stds. Currently only 9 position-based scores are computed; 45 rsID-based scores need mapping. |
| **Deliverable** | `pgs/pgs_calibrated.csv` with z-scores + percentiles for all 54 scores |

### 2. gnomAD Allele Frequencies
| | |
|---|---|
| **Impact** | 8/10 — Replace HWE estimates with real population AFs |
| **Effort** | 4/10 — 2 days |
| **Plan** | Download gnomAD v2.1 allele frequencies for the 30 SNPs that failed 1000G calibration. Build a local AF lookup. More accurate than default 0.25 MAF. |
| **Deliverable** | `reference/gnomad_af.json` + improved calibration for HWE-fallback traits |

### 3. rsID-based PGS Support
| | |
|---|---|
| **Impact** | 7/10 — Unlock 45/54 PGS scores currently unusable |
| **Effort** | 5/10 — 2-3 days |
| **Plan** | Build rsID → chr:pos mapping using the 1000G BIM file or dbSNP API. Convert rsID-based PGS scoring files to position-based format for PLINK. |
| **Deliverable** | 54/54 PGS scores computable |

### 4. Y-DNA Haplogroup from Genotyping Arrays
| | |
|---|---|
| **Impact** | 6/10 — Complete ancestry picture |
| **Effort** | 3/10 — 1 day |
| **Plan** | The genotyping array VCF has 47K chrY variants but uses different positions than ISOGG. Build a mapping of array probe IDs → Y-SNP positions. |
| **Deliverable** | Y-DNA haplogroup in deep ancestry output |

---

## Tier 2 — Medium Impact

### 5. Docker Container
| | |
|---|---|
| **Impact** | 8/10 — One-command setup on any machine |
| **Effort** | 4/10 — 1-2 days |
| **Plan** | `Dockerfile` with Python 3.12, PLINK, bcftools, tabix, and all Python deps. Mount data as volumes. |
| **Deliverable** | `docker run bluegen run --full --vcf /data/sample.vcf.gz` |

### 6. PDF Polish
| | |
|---|---|
| **Impact** | 6/10 — Professional-looking PDF |
| **Effort** | 2/10 — Few hours |
| **Plan** | Print CSS: page breaks, headers/footers, cover page, table of contents. |
| **Deliverable** | Polished `comprehensive_report.pdf` |

### 7. Multi-Sample Comparison
| | |
|---|---|
| **Impact** | 7/10 — Family/trio analysis |
| **Effort** | 5/10 — 2-3 days |
| **Plan** | Accept 2+ VCFs, compute PRS for each, generate comparison table. |
| **Deliverable** | `prs.py run --vcf a.vcf.gz,b.vcf.gz --compare` |

### 8. SNP Position Auto-Fixer
| | |
|---|---|
| **Impact** | 6/10 — More traits pass empirical calibration |
| **Effort** | 3/10 — 1 day |
| **Plan** | When a SNP fails 1000G matching, query dbSNP API for correct GRCh37 position. Auto-update CSV. |
| **Deliverable** | Fewer HWE fallbacks, more empirical calibrations |

---

## Tier 3 — Long Term

### 9. Cloud Pipeline (AWS/OCI)
| | |
|---|---|
| **Impact** | 8/10 — Scale to hundreds of samples |
| **Effort** | 8/10 — 2-3 weeks |
| **Plan** | Terraform + Batch/Compute. S3 for VCFs. Step Functions for orchestration. |

### 10. All of Us / UK Biobank Integration
| | |
|---|---|
| **Impact** | 10/10 — Massive population references |
| **Effort** | 10/10 — Requires data access (IRB, DUA) |

### 11. Plugin System
| | |
|---|---|
| **Impact** | 7/10 — Community contributions |
| **Effort** | 6/10 — 1-2 weeks |
| **Plan** | Plugin API: `plugins/<name>/manifest.json` + `annotate.py`. |

### 12. Mobile Dashboard (PWA)
| | |
|---|---|
| **Impact** | 5/10 — View reports on mobile |
| **Effort** | 5/10 — 1 week |
| **Plan** | Convert Streamlit to static site. Host on GitHub Pages. |

---

## Completed ✅

| Item | Version |
|------|---------|
| ClinVar pathogenic variants + confidence tiers | v1.0.0 |
| MedGen disease descriptions (local DB) | v1.0.0 |
| Pharmacogenomics (PharmGKB + ClinPGx + CPIC guidelines) | v1.0.0 |
| Deep Ancestry (mtDNA, sub-continental) | v1.0.0 |
| PRS expanded: 10→56 traits, 107→179 SNPs | v1.0.0 |
| PGS Catalog: 30→54 scores | v1.0.0 |
| Streamlit Dashboard (6 pages) | v1.0.0 |
| PDF reports (WeasyPrint) | v1.0.0 |
| CI/CD (GitHub Actions) | v1.0.0 |
| HWE fallback calibration | v1.0.0 |
| Cross-platform CWD auto-detect | v1.0.0 |
| `--stage` individual module execution | v1.0.0 |
| `--update-references` flag | v1.0.0 |
| Evidence level colors + CPIC explainer | v1.0.0 |
| Per-section limitation notes | v1.0.0 |
