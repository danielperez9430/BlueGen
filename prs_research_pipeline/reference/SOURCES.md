# Reference Data Sources & Attribution

> **Note:** The canonical, complete attribution document (including software tools,
> mirrors, and citation guidance) lives at the repository root: [`SOURCES.md`](../../SOURCES.md).
> This file is a quick reference for the `reference/` data directory.

This project uses the following public databases. All are free and publicly accessible.

## ClinVar — Clinical Variant Classifications

- **Source:** NCBI ClinVar, National Institutes of Health (NIH)
- **URL:** https://www.ncbi.nlm.nih.gov/clinvar/
- **Download:** `ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar.vcf.gz`
- **License:** Public Domain (U.S. Government work)
- **Citation:** Landrum MJ et al. *Nucleic Acids Research*. 2018. doi:10.1093/nar/gkx1153
- **Update frequency:** Monthly
- **Local path:** `reference/clinvar/`
- **Used by:** `scripts/clinical/clinvar_annotator.py`

## MedGen — Disease Definitions

- **Source:** NCBI MedGen, National Library of Medicine (NLM)
- **URL:** https://www.ncbi.nlm.nih.gov/medgen/
- **Download:** `ftp://ftp.ncbi.nlm.nih.gov/pub/medgen/{NAMES,MGDEF,MGSTY}.RRF.gz`
- **License:** Public Domain (U.S. Government work)
- **Update frequency:** Weekly (Wednesdays)
- **Local path:** `reference/medgen/`
- **Used by:** `scripts/clinical/medgen_enrich.py`

## ClinPGx / PharmGKB — Pharmacogenomics

- **Source:** ClinPGx (formerly PharmGKB), Stanford University
- **URL:** https://www.clinpgx.org/
- **Download:** `https://api.clinpgx.org/v1/download/file/data/`
- **License:** Creative Commons Attribution-ShareAlike 4.0 (CC BY-SA 4.0)
- **Citation:** Whirl-Carrillo M et al. *Clinical Pharmacology & Therapeutics*. 2012. doi:10.1038/clpt.2012.96
- **Update frequency:** Weekly
- **Local path:** `reference/clinpgx/`
- **Used by:** `scripts/clinical/clinpgx_sync.py`, `scripts/clinical/pharmgkb_annotator.py`

## 1000 Genomes Phase 3 — Population Reference Panel

- **Source:** International Genome Sample Resource (IGSR)
- **URL:** https://www.internationalgenome.org/
- **Download:** EBI FTP (via `scripts/setup/download_1000G_full.py`)
- **License:** Public Domain (EMBL-EBI)
- **Citation:** 1000 Genomes Project Consortium. *Nature*. 2015. doi:10.1038/nature15393
- **Update frequency:** Static (Phase 3 final release)
- **Local path:** `reference/1000G_full/`

## dbSNP — Reference SNP Identifiers

- **Source:** NCBI dbSNP, National Institutes of Health (NIH)
- **URL:** https://www.ncbi.nlm.nih.gov/snp/
- **Access:** NCBI Entrez E-utilities API (free, no key required)
- **License:** Public Domain (U.S. Government work)
- **Citation:** Sherry ST et al. *Nucleic Acids Research*. 2001. doi:10.1093/nar/29.1.308
- **Used by:** `scripts/setup/download_dbsnp.py`

## PGS Catalog — Polygenic Scores

- **Source:** PGS Catalog, University of Cambridge
- **URL:** https://www.pgscatalog.org/
- **License:** CC0 / Open Access
- **Citation:** Lambert SA et al. *Nature Genetics*. 2021. doi:10.1038/s41588-021-00783-5
- **Used by:** `scripts/benchmarking/pgs_catalog_integration.py`

## GWAS Catalog — Genome-Wide Association Studies

- **Source:** NHGRI-EBI GWAS Catalog
- **URL:** https://www.ebi.ac.uk/gwas/
- **License:** Public Domain (EMBL-EBI)
- **Citation:** Buniello A et al. *Nucleic Acids Research*. 2019. doi:10.1093/nar/gky1120
- **Used by:** `scripts/benchmarking/gwas_summary_stats.py`

## AADR — Allen Ancient DNA Resource (Archaic Reference Panel)

- **Source:** Harvard Dataverse / Reich Lab, Harvard Medical School
- **URL:** https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/FFIDCW
- **Download:** `python scripts/setup/download_aadr_reference.py` (builds from 1240K panel, ~7 GB download)
- **Pre-built:** https://archive.org/details/bluegen-archaic-reference (~188 MB PLINK binaries)
- **License:** CC0 1.0 Universal (Public Domain Dedication)
- **Citation:** Mallick S et al. *Scientific Data*. 2024. doi:10.1038/s41597-024-03031-7
- **Contents:** Altai, Vindija, Chagyrskaya Neanderthals + Denisova Denisovan + 341 modern reference individuals (27 populations)
- **Local path:** `reference/aadr/`
- **Used by:** `scripts/clinical/ancestry_deep.py` (AADR direct comparison mode)

## Archaic Introgression SNP Panel (Curated)

- **Source:** Published scientific literature
- **References:**
  - Sankararaman S et al. *Nature*. 2014. doi:10.1038/nature12961
  - Vernot B & Akey JM. *Science*. 2014. doi:10.1126/science.1245938
  - Browning SR et al. *Cell*. 2018. doi:10.1016/j.cell.2018.02.031
- **Local path:** `reference/archaic/archaic_introgression_snps.csv` (133 high-confidence SNPs)
- **Used by:** `scripts/clinical/ancestry_deep.py` (133-SNP panel fallback)

## Vindija Neanderthal Genome — Direct Archaic Comparison

- **Source:** Max Planck Institute for Evolutionary Anthropology
- **URL:** http://ftp.eva.mpg.de/neandertal/Vindija/VCF/Vindija33.19/
- **Download:** `python scripts/setup/download_vindija_reference.py --all` (~44 GB)
- **Pre-built mirror:** https://archive.org/details/bluegen-vindija-reference
- **License:** Public Domain (scientific data)
- **Citation:** Prüfer K et al. *Science*. 2017. doi:10.1126/science.aao1887
- **Specimen:** Vindija 33.19, ~30x coverage, hg19/GRCh37, MQ≥25 filter
- **Local path:** `reference/vindija/`
- **Used by:** `scripts/clinical/ancestry_deep.py` (Vindija direct comparison — gold standard)

---

## Disclaimer

> ⚠️ **RESEARCH USE ONLY — NOT FOR CLINICAL DIAGNOSIS**
>
> All data from these sources is used for research purposes. Clinical decisions should not be based on the output of this pipeline without confirmation by a certified clinical laboratory. This project is not affiliated with or endorsed by any of the listed data providers.
