# GWAS VCF Datasets — Manual Download Guide

Data already extracted: **525 SNP records** across 10 studies (2026-06-05).
Stored at: `prs_research_pipeline/gwas/extracted/gwas_extracted.json`

## How to refresh (manual)

The OpenGWAS API has **no VCF download endpoint**. Pre-signed OCI URLs expire
in ~2 hours. To re-download, visit each dataset page, click **Download**,
select **VCF files**, and pipe the fresh URL to tabix:

```bash
# Per dataset:
curl -L -o temp.vcf.gz "FRESH_URL_HERE"
tabix -p vcf temp.vcf.gz
tabix temp.vcf.gz chr22 > temp_chr22.vcf
python3 -c "
import json, subprocess
# Extract records matching our SNP database
...
"
rm temp.vcf.gz temp.vcf.gz.tbi
```

## Dataset index

| Dataset ID | Trait | Samples | OpenGWAS Page |
|---|---|---|---|
| `ebi-a-GCST006867` | Type 2 diabetes | 62,892 | [link](https://opengwas.io/datasets/ebi-a-GCST006867) |
| `ieu-b-40` | BMI | 681,275 | [link](https://opengwas.io/datasets/ieu-b-40) |
| `ebi-a-GCST90025986` | Hypertension | 484,598 | [link](https://opengwas.io/datasets/ebi-a-GCST90025986) |
| `ieu-b-109` | HDL cholesterol | 403,943 | [link](https://opengwas.io/datasets/ieu-b-109) |
| `ieu-b-110` | LDL cholesterol | 403,943 | [link](https://opengwas.io/datasets/ieu-b-110) |
| `ieu-b-111` | Triglycerides | 403,943 | [link](https://opengwas.io/datasets/ieu-b-111) |
| `ebi-a-GCST90000618` | Coronary artery disease | 547,261 | [link](https://opengwas.io/datasets/ebi-a-GCST90000618) |
| `ukb-b-5237` | Glucose | 460,749 | [link](https://opengwas.io/datasets/ukb-b-5237) |
| `ebi-a-GCST90013997` | Alzheimer's disease | 472,868 | [link](https://opengwas.io/datasets/ebi-a-GCST90013997) |
| `ukb-b-11349` | C-reactive protein | 462,897 | [link](https://opengwas.io/datasets/ukb-b-11349) |

The extracted data at `gwas/extracted/gwas_extracted.json` does not expire.
Only re-download if you need fresh results from an updated GWAS release.
