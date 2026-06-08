# 🧬 BlueGen

**Personal Genomics Platform** — Powered by the PRSKit polygenic risk score engine.

[![Python](https://img.shields.io/badge/python-3.10+-blue)](https://www.python.org/)
[![License](https://img.shields.io/badge/license-MIT-green)](LICENSE)
[![Version](https://img.shields.io/badge/version-1.0.0-orange)]()
[![PRSKit](https://img.shields.io/badge/engine-PRSKit-purple)]()

Turn a WGS VCF into a comprehensive personal genomics report: polygenic risk scores, pathogenic variants, pharmacogenomics, ancestry, and wellness traits — all offline, all free.

## 🚀 Quick Start

```bash
git clone https://github.com/danielperez9430/BlueGen.git
cd BlueGen
python3 -m venv venv && source venv/bin/activate
pip install -r prs_research_pipeline/requirements.txt

# Full analysis
python3 prs.py run --full --vcf your_sample.vcf.gz

# Interactive dashboard
venv/bin/streamlit run dashboard.py
```

## 📊 What It Does

| Module | Description | Data Source |
|--------|-------------|-------------|
| 🧬 **PRS Engine** | 56 traits, 179 SNPs, population-calibrated z-scores | Curated GWAS + 1000 Genomes |
| 🔬 **ClinVar** | Pathogenic/likely pathogenic variants with confidence tiers | NCBI ClinVar (4.4M records) |
| 💊 **PharmGKB** | Drug-gene interactions, CPIC guideline recommendations | CPIC/DPWG (218 guidelines) |
| 🌍 **Ancestry** | PCA + mtDNA haplogroup + sub-continental | 1000 Genomes (26 populations) |
| 🩺 **PGS Catalog** | 54 published polygenic scores for complex diseases | PGS Catalog (EBI) |
| 📖 **MedGen** | Disease definitions for ClinVar findings | NCBI MedGen (23K concepts) |
| 📊 **Dashboard** | 6-page interactive Streamlit app | All JSON outputs |

## 📁 Documentation

Full docs: [`prs_research_pipeline/README.md`](prs_research_pipeline/README.md)

Quick reference: [`USAGE.md`](USAGE.md)

Changelog: [`CHANGELOG.md`](CHANGELOG.md)

Data sources: [`prs_research_pipeline/reference/SOURCES.md`](prs_research_pipeline/reference/SOURCES.md)

## 🔧 Requirements

- Python 3.10+, PLINK 1.9 (bundled), bcftools + tabix
- macOS/Linux (Windows via WSL2)
- ~200 MB reference data (auto-downloaded on first run)
- ~25 GB optional: 1000 Genomes reference (one-time download)

## ⚠️ Disclaimer

**RESEARCH USE ONLY — NOT FOR CLINICAL DIAGNOSIS.** This platform is for research and educational purposes. Clinical decisions should not be based on its output without confirmation by a certified laboratory.
