# Data Sources, Licenses & Attribution

This document lists every external data source and third-party tool used by BlueGen,
with its provider, license, and the citation required (or recommended) when publishing
work derived from this pipeline. BlueGen's **code** is released under the [MIT License](LICENSE);
the MIT license does **not** extend to bundled data files, which remain under the
licenses listed below (notably `reference/clinpgx/`, which is CC BY-SA 4.0).

> All data sources listed here are free and publicly accessible. BlueGen is not
> affiliated with, endorsed by, or sponsored by any of the listed providers.

## Overview

| Source | Type | License | Used for |
|---|---|---|---|
| [1000 Genomes Phase 3](#1-1000-genomes-project-phase-3) | Population panel | Open access (IGSR) | PCA, ancestry, PRS calibration |
| [ClinVar](#2-clinvar) | Clinical variants | Public domain | Pathogenic variant screen |
| [MedGen](#3-medgen) | Disease concepts | Public domain | Disease definitions |
| [ClinPGx / PharmGKB](#4-clinpgx--pharmgkb-cpic) | Pharmacogenomics | CC BY-SA 4.0 | Drug–gene guidelines |
| [dbSNP](#5-dbsnp) | rsID reference | Public domain | SNP identifiers |
| [PGS Catalog](#6-pgs-catalog) | Polygenic scores | EBI open / CC0 metadata | PGS scoring (52 scores) |
| [GWAS Catalog](#7-nhgri-ebi-gwas-catalog) | GWAS associations | EMBL-EBI open terms | Curated PRS panel |
| [IEU OpenGWAS](#8-mrc-ieu-opengwas) | GWAS summary stats | Per-dataset (open) | Effect-size validation |
| [AADR](#9-aadr--allen-ancient-dna-resource) | Ancient DNA panel | CC0 1.0 | Archaic admixture |
| [Vindija Neanderthal](#10-vindija-neanderthal-genome) | Archaic genome | Open (MPI-EVA) | Direct archaic comparison |
| [Archaic SNP panel](#11-curated-archaic-introgression-snp-panel) | Curated literature | Published data | Archaic fallback mode |
| [Human Origins](#12-human-origins--neareastpublic-panel) | Population panel | Open (Reich Lab) | Sub-continental ancestry |
| [hg19 / GRCh37](#13-hg19--grch37-reference-genome) | Reference genome | Public domain | Alignment, coordinates |
| [ISOGG Y-tree](#14-isogg-y-dna-haplogroup-tree) | Y-DNA phylogeny | Free with attribution | Y-DNA haplogroups |
| [PhyloTree](#15-phylotree-mtdna) | mtDNA phylogeny | Free with attribution | mtDNA haplogroups |

---

## Data Sources

### 1. 1000 Genomes Project (Phase 3)

- **Provider:** International Genome Sample Resource (IGSR), EMBL-EBI
- **URL:** https://www.internationalgenome.org/
- **Download:** EBI FTP, via `prs_research_pipeline/scripts/setup/download_1000G_full.py`
- **Version:** Phase 3 final release (2,504 samples, 26 populations, 5 super-populations)
- **License:** Open access — data are fully public per IGSR data reuse policy
- **Citation:** 1000 Genomes Project Consortium; Auton A, et al. A global reference for human genetic variation. *Nature*. 2015;526(7571):68–74. doi:[10.1038/nature15393](https://doi.org/10.1038/nature15393)
- **Local path:** `prs_research_pipeline/reference/1000G/`, `reference/1000G_full/`
- **Used for:** PCA reference space, population assignment, LD pruning, PRS/PGS population calibration (`reference/population_distributions/` is derived from this panel)

### 2. ClinVar

- **Provider:** NCBI, U.S. National Library of Medicine, NIH
- **URL:** https://www.ncbi.nlm.nih.gov/clinvar/
- **Download:** `https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar.vcf.gz` via `scripts/setup/download_clinvar.py`
- **Version used:** `clinvar_20260530.vcf.gz` (~4.4M variant records; updated monthly)
- **License:** Public domain (U.S. Government work) — "available for unrestricted use"
- **Citation:** Landrum MJ, et al. ClinVar: improving access to variant interpretations and supporting evidence. *Nucleic Acids Res*. 2018;46(D1):D1062–D1067. doi:[10.1093/nar/gkx1153](https://doi.org/10.1093/nar/gkx1153)
- **Local path:** `prs_research_pipeline/reference/clinvar/` (see its own [SOURCES.md](prs_research_pipeline/reference/clinvar/SOURCES.md) for field-level details)
- **Used by:** `scripts/clinical/clinvar_annotator.py`

### 3. MedGen

- **Provider:** NCBI, U.S. National Library of Medicine, NIH
- **URL:** https://www.ncbi.nlm.nih.gov/medgen/
- **Download:** `https://ftp.ncbi.nlm.nih.gov/pub/medgen/{NAMES,MGDEF,MGSTY}.RRF.gz`
- **License:** Public domain (U.S. Government work). Some definitions aggregated by MedGen originate from third parties (e.g. GeneReviews, Orphanet) and retain their own attribution requirements.
- **Citation:** NCBI Resource Coordinators. Database resources of the National Center for Biotechnology Information. *Nucleic Acids Res*. 2018;46(D1):D8–D13. doi:[10.1093/nar/gkx1095](https://doi.org/10.1093/nar/gkx1095)
- **Local path:** `prs_research_pipeline/reference/medgen/` (~23K concepts used)
- **Used by:** `scripts/clinical/medgen_enrich.py`

### 4. ClinPGx / PharmGKB (CPIC)

- **Provider:** ClinPGx (formerly PharmGKB), Stanford University; guidelines from CPIC and DPWG
- **URL:** https://www.clinpgx.org/
- **Download:** `https://api.clinpgx.org/v1/download/file/data/` via `scripts/clinical/clinpgx_sync.py`
- **License:** Creative Commons Attribution-ShareAlike 4.0 ([CC BY-SA 4.0](https://creativecommons.org/licenses/by-sa/4.0/)) — attribution required; derived data must be shared under the same license
- **Citations:**
  - Whirl-Carrillo M, et al. Pharmacogenomics knowledge for personalized medicine. *Clin Pharmacol Ther*. 2012;92(4):414–417. doi:[10.1038/clpt.2012.96](https://doi.org/10.1038/clpt.2012.96)
  - Whirl-Carrillo M, et al. An evidence-based framework for evaluating pharmacogenomics knowledge for personalized medicine. *Clin Pharmacol Ther*. 2021;110(3):563–572. doi:[10.1002/cpt.2350](https://doi.org/10.1002/cpt.2350)
  - Relling MV, Klein TE. CPIC: Clinical Pharmacogenetics Implementation Consortium of the Pharmacogenomics Research Network. *Clin Pharmacol Ther*. 2011;89(3):464–467. doi:[10.1038/clpt.2010.279](https://doi.org/10.1038/clpt.2010.279)
- **Local path:** `prs_research_pipeline/reference/clinpgx/` (218 guidelines used)
- **Used by:** `scripts/clinical/pharmgkb_annotator.py`

### 5. dbSNP

- **Provider:** NCBI, NIH
- **URL:** https://www.ncbi.nlm.nih.gov/snp/
- **Access:** NCBI Entrez E-utilities API (free, no key required) via `scripts/setup/download_dbsnp.py`
- **License:** Public domain (U.S. Government work)
- **Citation:** Sherry ST, et al. dbSNP: the NCBI database of genetic variation. *Nucleic Acids Res*. 2001;29(1):308–311. doi:[10.1093/nar/29.1.308](https://doi.org/10.1093/nar/29.1.308)
- **Used for:** rsID resolution and position validation of the curated SNP panel

### 6. PGS Catalog

- **Provider:** PGS Catalog (EMBL-EBI / University of Cambridge / HDR-UK)
- **URL:** https://www.pgscatalog.org/
- **License:** Open access under EMBL-EBI terms of use; catalog metadata distributed as CC0. Individual scoring files inherit the terms of their source publications.
- **Citation:** Lambert SA, et al. The Polygenic Score Catalog as an open database for reproducibility and systematic evaluation. *Nat Genet*. 2021;53(4):420–425. doi:[10.1038/s41588-021-00783-5](https://doi.org/10.1038/s41588-021-00783-5)
- **Used by:** `scripts/benchmarking/pgs_catalog_integration.py` (52 scores computed and population-calibrated). Per-score attribution (PGS IDs and source publications) is embedded in the pipeline's JSON output.

### 7. NHGRI-EBI GWAS Catalog

- **Provider:** NHGRI and EMBL-EBI
- **URL:** https://www.ebi.ac.uk/gwas/
- **License:** Open access under EMBL-EBI terms of use
- **Citation:** Buniello A, et al. The NHGRI-EBI GWAS Catalog of published genome-wide association studies, targeted arrays and summary statistics 2019. *Nucleic Acids Res*. 2019;47(D1):D1005–D1012. doi:[10.1093/nar/gky1120](https://doi.org/10.1093/nar/gky1120)
- **Used for:** Source of effect sizes and risk alleles for the curated PRS panel (~56 traits, 206 SNP rows, 187 unique rsIDs). The underlying primary GWAS publications for each SNP are recorded in `prs_research_pipeline/interpretations/` and the science documentation.

### 8. MRC IEU OpenGWAS

- **Provider:** MRC Integrative Epidemiology Unit, University of Bristol
- **URL:** https://opengwas.io/
- **License:** Per-dataset; the 10 datasets used (UK Biobank and EBI GWAS derivatives — see `prs_research_pipeline/scripts/setup/GWAS_DATASETS.md`) are openly redistributable summary statistics
- **Citation:** Elsworth B, et al. The MRC IEU OpenGWAS data infrastructure. *bioRxiv*. 2020. doi:[10.1101/2020.08.10.244293](https://doi.org/10.1101/2020.08.10.244293)
- **Local path:** `prs_research_pipeline/gwas/extracted/gwas_extracted.json` (525 SNP records, 10 studies)
- **Used for:** Independent effect-size validation of the curated panel

### 9. AADR — Allen Ancient DNA Resource

- **Provider:** Reich Lab, Harvard Medical School / Harvard Dataverse
- **URL:** https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/FFIDCW
- **Download:** `scripts/setup/download_aadr_reference.py` (1240K panel, ~7 GB); pre-built PLINK subset at https://archive.org/details/bluegen-archaic-reference (~188 MB)
- **License:** CC0 1.0 Universal (Public Domain Dedication)
- **Citations** (the AADR requests both):
  - Mallick S, Micco A, Mah M, et al. The Allen Ancient DNA Resource (AADR) a curated compendium of ancient human genomes. *Sci Data*. 2024;11:182. doi:[10.1038/s41597-024-03031-7](https://doi.org/10.1038/s41597-024-03031-7)
  - Mallick S, Reich D. The Allen Ancient DNA Resource (AADR): A curated compendium of ancient human genomes. Harvard Dataverse, V54.1. doi:[10.7910/DVN/FFIDCW](https://doi.org/10.7910/DVN/FFIDCW)
- **Contents used:** Altai, Vindija and Chagyrskaya Neanderthals; Denisova; 341 modern reference individuals (27 populations); 1.23M SNPs
- **Local path:** `prs_research_pipeline/reference/aadr/`
- **Used by:** `scripts/clinical/ancestry_deep.py` (AADR direct-comparison mode)

### 10. Vindija Neanderthal Genome

- **Provider:** Max Planck Institute for Evolutionary Anthropology (MPI-EVA)
- **URL:** http://ftp.eva.mpg.de/neandertal/Vindija/VCF/Vindija33.19/
- **Download:** `scripts/setup/download_vindija_reference.py --all` (~44 GB); mirror at https://archive.org/details/bluegen-vindija-reference
- **License:** Freely available scientific data (MPI-EVA public FTP); cite the source publication
- **Citation:** Prüfer K, et al. A high-coverage Neandertal genome from Vindija Cave in Croatia. *Science*. 2017;358(6363):655–658. doi:[10.1126/science.aao1887](https://doi.org/10.1126/science.aao1887)
- **Specimen:** Vindija 33.19, ~30× coverage, hg19/GRCh37 coordinates, MQ≥25 filter
- **Local path:** `prs_research_pipeline/reference/vindija/`
- **Used by:** `scripts/clinical/ancestry_deep.py` (direct archaic comparison — gold standard mode)

### 11. Curated Archaic Introgression SNP Panel

- **Provider:** Compiled by this project from published literature (133 high-confidence SNPs)
- **Citations:**
  - Sankararaman S, et al. The genomic landscape of Neanderthal ancestry in present-day humans. *Nature*. 2014;507(7492):354–357. doi:[10.1038/nature12961](https://doi.org/10.1038/nature12961)
  - Vernot B, Akey JM. Resurrecting surviving Neandertal lineages from modern human genomes. *Science*. 2014;343(6174):1017–1021. doi:[10.1126/science.1245938](https://doi.org/10.1126/science.1245938)
  - Browning SR, et al. Analysis of human sequence data reveals two pulses of archaic Denisovan admixture. *Cell*. 2018;173(1):53–61.e9. doi:[10.1016/j.cell.2018.02.031](https://doi.org/10.1016/j.cell.2018.02.031)
- **Local path:** `prs_research_pipeline/reference/archaic/archaic_introgression_snps.csv`
- **Used by:** `scripts/clinical/ancestry_deep.py` (fallback mode when AADR/Vindija are absent)

### 12. Human Origins / NearEastPublic Panel

- **Provider:** Reich Lab, Harvard Medical School
- **URL:** https://reich.hms.harvard.edu/datasets
- **Download:** `scripts/setup/download_ashkenazi_reference.sh` (NearEastPublic dataset, ~150 MB)
- **License:** Freely downloadable for research use per Reich Lab dataset terms
- **Citation:** Lazaridis I, et al. Genomic insights into the origin of farming in the ancient Near East. *Nature*. 2016;536(7617):419–424. doi:[10.1038/nature19310](https://doi.org/10.1038/nature19310)
- **Contents used:** Human Origins genotype panel (2,068 individuals) including Ashkenazi Jewish and diverse European populations
- **Local path:** `prs_research_pipeline/reference/human_origins/`
- **Used for:** Sub-continental / Ashkenazi ancestry refinement

### 13. hg19 / GRCh37 Reference Genome

- **Provider:** Genome Reference Consortium (GRC); distributed via UCSC Genome Browser
- **URL:** https://hgdownload.soe.ucsc.edu/goldenPath/hg19/
- **License:** Freely available; sequence data are public domain
- **Citation:** Church DM, et al. Modernizing reference genome assemblies. *PLoS Biol*. 2011;9(7):e1001091. doi:[10.1371/journal.pbio.1001091](https://doi.org/10.1371/journal.pbio.1001091)
- **Local path:** `prs_research_pipeline/reference/hg19/` (FASTA + BWA index)
- **Used for:** FASTQ→BAM alignment (optional) and as the coordinate system for the whole pipeline (GRCh37)

### 14. ISOGG Y-DNA Haplogroup Tree

- **Provider:** International Society of Genetic Genealogy (ISOGG)
- **URL:** https://isogg.org/tree/
- **Version:** ISOGG Y-DNA Haplogroup Tree 2024 (defining SNPs for major haplogroups)
- **License:** Freely accessible; attribution to ISOGG requested when referencing the tree
- **Used by:** `scripts/clinical/ancestry_deep.py` (Y-DNA haplogroup calling; SNP table embedded in code)

### 15. PhyloTree (mtDNA)

- **Provider:** van Oven M — PhyloTree.org
- **URL:** https://www.phylotree.org/
- **License:** Freely accessible for research with citation
- **Citation:** van Oven M, Kayser M. Updated comprehensive phylogenetic tree of global human mitochondrial DNA variation. *Hum Mutat*. 2009;30(2):E386–E394. doi:[10.1002/humu.20921](https://doi.org/10.1002/humu.20921)
- **Used by:** `scripts/clinical/ancestry_deep.py` (mtDNA haplogroup calling; defining-SNP table embedded in code)

---

## Third-Party Software

| Tool | License | Citation / URL |
|---|---|---|
| **PLINK 1.9** | GPLv3 | Chang CC, et al. Second-generation PLINK. *GigaScience*. 2015;4:7. doi:[10.1186/s13742-015-0047-8](https://doi.org/10.1186/s13742-015-0047-8) — https://www.cog-genomics.org/plink/ |
| **PLINK 2.0** | GPLv3 | Same citation — https://www.cog-genomics.org/plink/2.0/ |
| **bcftools / HTSlib** | MIT/Expat + BSD | Danecek P, et al. Twelve years of SAMtools and BCFtools. *GigaScience*. 2021;10(2):giab008. doi:[10.1093/gigascience/giab008](https://doi.org/10.1093/gigascience/giab008) |
| **tabix** | MIT/Expat | Li H. Tabix: fast retrieval of sequence features from generic TAB-delimited files. *Bioinformatics*. 2011;27(5):718–719. doi:[10.1093/bioinformatics/btq671](https://doi.org/10.1093/bioinformatics/btq671) |
| **BWA** | GPLv3 | Li H, Durbin R. Fast and accurate short read alignment with Burrows-Wheeler transform. *Bioinformatics*. 2009;25(14):1754–1760. doi:[10.1093/bioinformatics/btp324](https://doi.org/10.1093/bioinformatics/btp324) |
| **Python scientific stack** (NumPy, SciPy, pandas, scikit-learn, Streamlit) | BSD / Apache 2.0 | See `prs_research_pipeline/requirements.txt` for exact versions |

---

## Mirrors

Reference-data snapshots are mirrored by this project on the Internet Archive for
reproducibility. The mirrors redistribute only data whose license permits it
(public domain / CC0 / open access), unmodified except for format conversion:

- https://archive.org/details/bluegen-reference-data — 1000 Genomes, hg19, ClinVar, MedGen, ClinPGx bundle (~65 GB)
- https://archive.org/details/bluegen-archaic-reference — AADR-derived PLINK subset (~188 MB)
- https://archive.org/details/bluegen-vindija-reference — Vindija VCF mirror

---

## How to Cite

If you use BlueGen in academic work, cite the individual data sources above
according to their requirements (ClinVar, PharmGKB/CPIC, PGS Catalog and AADR
explicitly require citation in publications), plus the specific source publications
of any PGS Catalog scores or GWAS effect sizes you report.

---

## Disclaimer

> ⚠️ **RESEARCH USE ONLY — NOT FOR CLINICAL DIAGNOSIS**
>
> All data from these sources is used for research purposes. Clinical decisions
> must not be based on the output of this pipeline without confirmation by a
> certified clinical laboratory and consultation with a qualified healthcare
> professional. This project is not affiliated with or endorsed by any of the
> listed data providers.
