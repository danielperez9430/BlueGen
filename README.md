<p align="center">
  <img src="assets/blugen-logo.png" alt="BlueGen Logo" width="300">
</p>

# 🧬 BlueGen

**Personal Genomics Platform** — Powered by the PRSKit polygenic risk score engine.

[![Python](https://img.shields.io/badge/python-3.10+-blue)](https://www.python.org/)
[![License](https://img.shields.io/badge/license-MIT-green)](LICENSE)
[![Version](https://img.shields.io/badge/version-2.0.1-orange)]()
[![PRSKit](https://img.shields.io/badge/engine-PRSKit-purple)]()

Turn a WGS VCF into a comprehensive personal genomics report: polygenic risk scores, pathogenic variants, pharmacogenomics, ancestry, and wellness traits — all offline, all free.

## 🚀 Quick Start

```bash
git clone https://github.com/danielperez9430/BlueGen.git
cd BlueGen
python3 -m venv venv && source venv/bin/activate
pip install -r prs_research_pipeline/requirements.txt

# Recommended: download the genome-wide 1000G reference first (~25 GB, one-time).
# Without it, the pipeline silently falls back to a chr22-only reference for
# PCA/ancestry/LD pruning (see IMPROVEMENT_PLAN.md TIER 0.2).
python3 prs_research_pipeline/scripts/setup/download_1000G_full.py

# Full analysis
python3 prs.py run --full --vcf your_sample.vcf.gz

# Interactive dashboard
venv/bin/streamlit run dashboard.py
```

## 📊 What It Does

| Module | Description | Data Source |
|--------|-------------|-------------|
| 🧬 **PRS Engine** | ~56 traits, 206 SNP rows (187 unique rsIDs), population-calibrated z-scores | Curated GWAS + 1000 Genomes |
| 🔬 **ClinVar** | Pathogenic/likely pathogenic variants with confidence tiers | NCBI ClinVar (4.4M records) |
| 💊 **PharmGKB** | Drug-gene interactions, CPIC guideline recommendations | CPIC/DPWG (218 guidelines) |
| 🌍 **Ancestry** | PCA + mtDNA/Y-DNA haplogroups + sub-continental | 1000 Genomes (26 populations) |
| 🦴 **Archaic DNA** | Neanderthal/Denisovan admixture via AADR direct comparison | Allen Ancient DNA Resource (1.23M SNPs) |
| 🩺 **PGS Catalog** | 52/57 scored and population-calibrated (z-score + percentile vs 1000G EUR); a handful excluded for being genome-wide-scale (>500K variants, impractical for single-sample scoring) | PGS Catalog (EBI) |
| 📖 **MedGen** | Disease definitions for ClinVar findings | NCBI MedGen (23K concepts) |
| 📊 **Dashboard** | 6-page interactive Streamlit app | All JSON outputs |

## 🏗️ Architecture

```mermaid
flowchart LR
    VCF["WGS VCF<br/>(.vcf.gz)"]
    QC["Genotype processing<br/>QC filter → LD-prune → PCA vs 1000 Genomes"]

    PRS["PRS + PGS scoring<br/>PLINK --score → PCA-adjust → population-calibrate<br/>56 curated traits + 52 PGS Catalog scores"]
    CLIN["ClinVar + MedGen<br/>pathogenic/likely-pathogenic screen"]
    PHARM["PharmGKB<br/>CPIC drug-gene guidelines"]
    ANC["Ancestry + Archaic DNA<br/>PCA population, mtDNA/Y-DNA, Neanderthal/Denisovan"]

    SSST["SSST consolidation<br/>schema-validated canonical JSON"]
    VALID["Validation suite<br/>adversarial + calibration checks"]
    REPORT["Comprehensive report<br/>bilingual HTML + Streamlit dashboard"]

    VCF --> QC
    QC --> PRS
    QC --> CLIN
    QC --> PHARM
    QC --> ANC
    PRS --> SSST
    CLIN --> SSST
    PHARM --> SSST
    ANC --> SSST
    SSST --> VALID
    VALID --> REPORT
```

Everything above runs locally — no genetic data ever leaves the machine.

## 📁 Documentation

Full docs: [`prs_research_pipeline/README.md`](prs_research_pipeline/README.md)

Quick reference: [`USAGE.md`](USAGE.md)

Changelog: [`CHANGELOG.md`](CHANGELOG.md)

Data sources: [`prs_research_pipeline/reference/SOURCES.md`](prs_research_pipeline/reference/SOURCES.md)

## 🔧 Requirements

- **Python 3.10+** — `pip install -r prs_research_pipeline/requirements.txt`
- **PLINK v1.90b7.2+** — [Download](https://www.cog-genomics.org/plink/) (free, GPLv3)
- **PLINK 2.0** — [Download](https://www.cog-genomics.org/plink/2.0/) (optional, for advanced QC)
- **bcftools + tabix** — `brew install bcftools tabix` (macOS) or `apt install bcftools tabix` (Linux)
- **BWA** — `brew install bwa` (macOS) or `apt install bwa` (Linux) (optional, for FASTQ → BAM alignment)
- macOS/Linux (Windows via WSL2)
- ~200 MB reference data (auto-downloaded on first run)
- ~65 GB reference data bundle (optional — from [archive.org](https://archive.org/details/bluegen-reference-data) snapshot or `prs_research_pipeline/scripts/setup/` to fetch latest from public sources)
  - [`bluegen-reference-data`](https://archive.org/details/bluegen-reference-data) — 1000 Genomes, hg19, ClinVar, MedGen, ClinPGx
  - [`bluegen-pgs-cache`](https://archive.org/details/bluegen-pgs-cache) — 56 PGS Catalog scoring files
  - [`bluegen-archaic-reference`](https://archive.org/details/bluegen-archaic-reference) — AADR 1240K archaic panel (Neanderthal + Denisovan), pre-converted to PLINK
  - [`bluegen-vindija-reference`](https://archive.org/details/bluegen-vindija-reference) — Vindija Neanderthal genome VCFs (chr1–22, hg19, ~44 GB)
  - > **Maintainer:** `python archive_upload.py -j 8` to refresh snapshots

### System Tools

PLINK is **not bundled** — download the correct build for your OS:

| Tool | Version | macOS (Apple Silicon) | macOS (Intel) | Linux |
|------|---------|----------------------|---------------|-------|
| PLINK 1.9 | [v1.90b7.2](https://www.cog-genomics.org/plink/) | `plink` (Rosetta) | `plink_mac` | `plink_linux` |
| PLINK 2.0 | [v2.0.0-a.7.1](https://www.cog-genomics.org/plink/2.0/) | `plink2_mac_arm64` | `plink2_mac` | `plink2_linux` |

Place the binaries in your `PATH` or symlink them into the project root.

## ⚠️ Disclaimer

**RESEARCH USE ONLY — NOT FOR CLINICAL DIAGNOSIS.** This platform is for research and educational purposes. Clinical decisions should not be based on its output without confirmation by a certified laboratory.
