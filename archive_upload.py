#!/usr/bin/env python3
"""
BlueGen Archive.org Upload Script
==================================
Uploads public reference datasets to archive.org under the `bluegen` collection.

Items created:
  - bluegen-reference-data  (~61GB)  1000 Genomes, hg19, ClinVar, MedGen, ClinPGx
  - bluegen-pgs-cache        (~2GB)   PGS Catalog scoring files (56 scores)

Requirements:
  pip install internetarchive
  ia configure   # or set IA_ACCESS_KEY / IA_SECRET_KEY env vars

Usage:
  python archive_upload.py --dry-run        # validate without uploading
  python archive_upload.py                  # upload (skips duplicates via checksum)
  python archive_upload.py --only ref        # only reference data
  python archive_upload.py --only pgs        # only PGS cache

Archive.org S3 credentials:
  Get your keys at https://archive.org/account/s3.php
  Then either:
    export IA_ACCESS_KEY=xxx
    export IA_SECRET_KEY=yyy
  Or run `ia configure` once.
"""

import hashlib
import json
import os
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from threading import Lock

ROOT = Path(__file__).resolve().parent
REFERENCE_DIR = ROOT / "prs_research_pipeline" / "reference"
PGS_DIR = ROOT / "prs_research_pipeline" / "pgs"

# Items appear under https://archive.org/details/@danielperez9430

# ── archive.org items ──────────────────────────────────────────────────────

ITEMS = {
    "bluegen-vindija-reference": {
        "title": "BlueGen — Vindija Neanderthal Reference (hg19, MQ≥25)",
        "description": (
            "Vindija 33.19 Neanderthal genome VCF files (per-chromosome) for direct "
            "archaic admixture analysis in the BlueGen Personal Genomics Platform.<br><br>"
            "<b>Contents:</b><br>"
            "<ul>"
            "<li>22 VCF files (chr1–chr22) with MQ≥25 and MAPQ≥100 filters</li>"
            "<li>Tabix indexes (.tbi) for all chromosomes</li>"
            "<li>~2.6M variant positions across the genome</li>"
            "</ul>"
            "<b>Source:</b> Max Planck Institute for Evolutionary Anthropology<br>"
            "URL: http://ftp.eva.mpg.de/neandertal/Vindija/VCF/Vindija33.19/<br>"
            "<b>Specimen:</b> Vindija 33.19, ~30x coverage, hg19/GRCh37<br>"
            "<b>Citation:</b> Prüfer et al. (2017) <i>Science</i>. "
            "doi:<a href='https://doi.org/10.1126/science.aao1887'>10.1126/science.aao1887</a><br><br>"
            "<b>License:</b> Public Domain (scientific data)<br><br>"
            "<b>Why this exists:</b> The original MPI FTP server can be slow or unavailable. "
            "This archive.org mirror ensures reliable access. Users can also download directly "
            "from MPI via <code>scripts/setup/download_vindija_reference.py</code>."
        ),
        "subject": [
            "genomics", "Neanderthal", "Vindija", "ancient DNA", "archaic hominin",
            "hg19", "bioinformatics", "paleogenomics",
        ],
        "creator": "BlueGen Pipeline / Max Planck Institute EVA",
        "date": "2026-06-24",
        "licenseurl": "https://creativecommons.org/publicdomain/zero/1.0/",
        "mediatype": "data",
        "collection": "opensource",
        "source_paths": [
            ("vindija", "Vindija Neanderthal VCFs + indexes + manifest (44 GB)"),
        ],
    },
    "bluegen-archaic-reference": {
        "title": "BlueGen — Archaic Reference Panel (AADR 1240K — Neanderthal/Denisovan)",
        "description": (
            "Archaic hominin reference genotypes extracted from the Allen Ancient DNA "
            "Resource (AADR v66.p1) for Neanderthal and Denisovan admixture analysis "
            "in the BlueGen Personal Genomics Platform.<br><br>"
            "<b>Contents:</b><br>"
            "<ul>"
            "<li><b>aadr_archaic</b> — PLINK binaries (bed/bim/fam) for Altai, Vindija, "
            "Chagyrskaya Neanderthals and Denisova Denisovan (~1.24M SNPs, hg19)</li>"
            "<li><b>aadr_modern</b> — PLINK binaries for ~500 modern reference individuals "
            "from global populations</li>"
            "<li><b>aadr_manifest.json</b> — provenance, version, individual metadata</li>"
            "</ul>"
            "<b>Source:</b> Mallick et al. (2024) The Allen Ancient DNA Resource: "
            "A curated compendium of ancient human genomes. <i>Scientific Data</i>.<br>"
            "DOI: <a href='https://doi.org/10.7910/DVN/FFIDCW'>10.7910/DVN/FFIDCW</a><br><br>"
            "<b>License:</b> <a href='https://creativecommons.org/publicdomain/zero/1.0/'>"
            "CC0 1.0 Universal</a> (Public Domain Dedication)<br><br>"
            "<b>Why this exists:</b> The full AADR dataset is ~12 GB and requires "
            "Eigenstrat→PLINK conversion. This pre-built PLINK snapshot enables direct "
            "archaic admixture analysis in BlueGen without the multi-hour download and "
            "conversion step. Users who need the complete AADR with all samples can "
            "re-run <code>scripts/setup/download_aadr_reference.py</code>."
        ),
        "subject": [
            "genomics", "ancient DNA", "Neanderthal", "Denisovan", "archaic hominin",
            "AADR", "Allen Ancient DNA Resource", "PLINK", "bioinformatics",
        ],
        "creator": "BlueGen Pipeline",
        "date": "2026-06-23",
        "licenseurl": "https://creativecommons.org/publicdomain/zero/1.0/",
        "mediatype": "data",
        "collection": "opensource",
        "source_paths": [
            ("aadr/aadr_archaic.bed", "Archaic genotypes — PLINK binary (CC0)"),
            ("aadr/aadr_archaic.bim", "Archaic variant info — PLINK (CC0)"),
            ("aadr/aadr_archaic.fam", "Archaic individual info — PLINK (CC0)"),
            ("aadr/aadr_modern.bed", "Modern reference genotypes — PLINK (CC0)"),
            ("aadr/aadr_modern.bim", "Modern reference variants — PLINK (CC0)"),
            ("aadr/aadr_modern.fam", "Modern reference individuals — PLINK (CC0)"),
            ("aadr/aadr_manifest.json", "Provenance & metadata (CC0)"),
        ],
    },
    "bluegen-reference-data": {
        "title": "BlueGen — Public Reference Data (1000 Genomes, hg19, ClinVar, MedGen, ClinPGx)",
        "description": (
            "Frozen snapshot of public reference datasets used by BlueGen v1.0.0 "
            "(Personal Genomics Platform).<br><br>"
            "<b>Contents &amp; Licenses (per dataset):</b><br>"
            "<ul>"
            "<li><b>1000 Genomes Phase 3</b> — merged PLINK binaries + ancestry panel<br>"
            "&nbsp;&nbsp;License: Public Domain (EMBL-EBI / IGSR)<br>"
            "&nbsp;&nbsp;Citation: 1000 Genomes Project Consortium. <i>Nature</i> 526, 68-74 (2015). "
            "doi:<a href='https://doi.org/10.1038/nature15393'>10.1038/nature15393</a></li>"
            "<li><b>hg19/GRCh37</b> — BWA index (fa, bwt, pac, sa, amb, ann, fai)<br>"
            "&nbsp;&nbsp;License: Public Domain (UCSC Genome Browser)<br>"
            "&nbsp;&nbsp;URL: <a href='https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/'>UCSC hg19</a></li>"
            "<li><b>ClinVar</b> — pathogenic variant classifications (VCF + annotations)<br>"
            "&nbsp;&nbsp;License: Public Domain (U.S. Government / NIH / NCBI)<br>"
            "&nbsp;&nbsp;Citation: Landrum MJ et al. <i>Nucleic Acids Research</i> 46(D1), D1062-D1067 (2018). "
            "doi:<a href='https://doi.org/10.1093/nar/gkx1153'>10.1093/nar/gkx1153</a></li>"
            "<li><b>MedGen</b> — disease definitions (RRF)<br>"
            "&nbsp;&nbsp;License: Public Domain (U.S. Government / NIH / NLM)</li>"
            "<li><b>ClinPGx / PharmGKB</b> — pharmacogenomic annotations<br>"
            "&nbsp;&nbsp;License: <a href='https://creativecommons.org/licenses/by-sa/4.0/'>CC BY-SA 4.0</a> "
            "(Stanford University)<br>"
            "&nbsp;&nbsp;Citation: Whirl-Carrillo M et al. <i>Clinical Pharmacology &amp; Therapeutics</i> "
            "92(4), 414-417 (2012). doi:<a href='https://doi.org/10.1038/clpt.2012.96'>10.1038/clpt.2012.96</a></li>"
            "<li><b>Population distributions</b> — 26-population allele frequency panels<br>"
            "&nbsp;&nbsp;License: Public Domain (derived from 1000 Genomes)</li>"
            "</ul>"
            "<b>Overall item license:</b> This item aggregates datasets with different licenses. "
            "The most restrictive applicable license is "
            "<a href='https://creativecommons.org/licenses/by-sa/4.0/'>CC BY-SA 4.0</a> "
            "(covering ClinPGx/PharmGKB content). All other datasets are Public Domain.<br><br>"
            "<b>Why this exists:</b> These public databases may change or go offline. "
            "This frozen copy ensures exact reproducibility of BlueGen pipeline results. "
            "Users who want the freshest data can re-run the download scripts in "
            "<code>prs_research_pipeline/scripts/setup/</code>.<br><br>"
            "<b>Full provenance:</b> See <code>SOURCES.md</code> included in this upload."
        ),
        "subject": [
            "genomics", "reference genome", "1000 Genomes", "GRCh37", "hg19",
            "ClinVar", "MedGen", "ClinPGx", "PharmGKB", "PLINK", "bioinformatics",
        ],
        "creator": "BlueGen Pipeline",
        "date": "2026-06-04",  # manifest download_date
        "licenseurl": "https://creativecommons.org/licenses/by-sa/4.0/",
        "mediatype": "data",
        "collection": "opensource",
        "source_paths": [
            ("1000G", "1000 Genomes chr22 VCF — test/validation (Public Domain)"),
            ("1000G_full", "1000 Genomes Phase 3 — merged PLINK binaries (Public Domain)"),
            ("hg19", "hg19/GRCh37 reference genome — BWA index (Public Domain)"),
            ("clinvar", "ClinVar pathogenic variant classifications (Public Domain)"),
            ("medgen", "MedGen disease definitions (Public Domain)"),
            ("clinpgx", "ClinPGx pharmacogenomic annotations (CC BY-SA 4.0)"),
            ("population_distributions", "Population allele frequency panels — 26 pops (Public Domain)"),
            ("SOURCES.md", "Data provenance & attribution"),
        ],
    },
    "bluegen-european-reference": {
        "title": "BlueGen — European Sub-Continental Reference Panel (Human Origins + 1000G)",
        "description": (
            "Extended European sub-continental ancestry reference panel for the BlueGen "
            "Personal Genomics Platform. Combines 1000 Genomes European populations with "
            "the Human Origins dataset (Lazaridis et al. 2016) to enable fine-scale "
            "classification including Ashkenazi Jewish ancestry.<br><br>"
            "<b>Contents:</b><br>"
            "<ul>"
            "<li><b>european_aj_subset</b> — PLINK binaries (bed/bim/fam) with European + "
            "Ashkenazi Jewish samples from Human Origins (~600K SNPs, hg19)</li>"
            "<li><b>population_labels.txt</b> — Population assignments for all samples "
            "(IBS, GBR, CEU, TSI, FIN, AJ)</li>"
            "<li><b>download_ashkenazi_reference.sh</b> — Script to reproduce the download "
            "and processing from source</li>"
            "</ul>"
            "<b>Populations:</b> Iberian (IBS), British (GBR), NW European (CEU), "
            "Tuscan Italian (TSI), Finnish (FIN), Ashkenazi Jewish (AJ)<br>"
            "<b>Source (Human Origins):</b> Lazaridis et al. 2016, Nature 536, 419-424<br>"
            "URL: <a href='https://reich.hms.harvard.edu/datasets'>https://reich.hms.harvard.edu/datasets</a><br>"
            "<b>Source (1000G):</b> 1000 Genomes Phase 3, The International Genome Sample Resource<br>"
            "<b>License:</b> Public Domain (scientific data)<br><br>"
            "<b>Usage:</b> Place in <code>reference/human_origins/</code>. The subcontinental "
            "PCA classifier (<code>subcontinental_pca.py</code>) automatically detects and "
            "uses this extended reference when available."
        ),
        "subject": [
            "genomics", "ancestry", "Ashkenazi Jewish", "European", "PCA",
            "population genetics", "1000 Genomes", "Human Origins",
        ],
        "creator": "BlueGen Pipeline / Reich Lab (Harvard Medical School)",
        "date": "2026-06-26",
        "licenseurl": "https://creativecommons.org/publicdomain/zero/1.0/",
        "mediatype": "data",
        "collection": "opensource",
        "source_paths": [
            ("aj_centroid.json", "Ashkenazi Jewish PCA centroid (literature-based)"),
            (".gitkeep", "Directory structure marker"),
        ],
    },
    "bluegen-pgs-cache": {
        "title": "BlueGen — PGS Catalog Scoring Files (56 polygenic scores)",
        "description": (
            "Frozen snapshot of 56 PGS Catalog scoring files used by BlueGen v1.0.0 "
            "(Personal Genomics Platform).<br><br>"
            "Each subdirectory contains a PGS scoring file (PGS######.txt.gz) downloaded "
            "from <a href='https://www.pgscatalog.org/'>PGS Catalog</a> (University of "
            "Cambridge).<br><br>"
            "<b>License:</b> <a href='https://creativecommons.org/publicdomain/zero/1.0/'>CC0 1.0 Universal</a> "
            "(Public Domain Dedication). All PGS Catalog data is openly available for any use.<br>"
            "<b>Citation:</b> Lambert SA et al. <i>Nature Genetics</i> 53, 420-425 (2021). "
            "doi:<a href='https://doi.org/10.1038/s41588-021-00783-5'>10.1038/s41588-021-00783-5</a><br><br>"
            "<b>Why this exists:</b> PGS weights are updated periodically by the PGS Catalog consortium. "
            "This frozen copy ensures exact reproducibility of polygenic risk scores computed by "
            "the BlueGen pipeline. Users who want the latest scores can re-run "
            "<code>prs_research_pipeline/scripts/benchmarking/pgs_catalog_integration.py</code>."
        ),
        "subject": [
            "genomics", "polygenic risk scores", "PRS", "PGS Catalog",
            "GWAS", "bioinformatics",
        ],
        "creator": "BlueGen Pipeline",
        "date": "2026-06-08",
        "licenseurl": "https://creativecommons.org/publicdomain/zero/1.0/",
        "mediatype": "data",
        "collection": "opensource",
        "source_paths": [
            (".", "PGS Catalog scoring files — 56 scores (CC0)"),
        ],
    },
}


# ── helpers ─────────────────────────────────────────────────────────────────

# Files to EXCLUDE from upload — these contain personal DNA-derived data
EXCLUDE_PATTERNS = {
    # PGS test scores run against personal genotype
    "test_score.log", "test_score.nopred", "test_score.nosex",
    "test_score.txt", "test_score.profile",
    "*.nopred", "*.nosex", "*.profile",  # per-score test results
    # Aggregate results derived from personal data
    "concordance.json", "pgs_results.csv",
}


def human(n: int) -> str:
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if n < 1024:
            return f"{n:.1f} {unit}"
        n /= 1024
    return f"{n:.1f} PB"


def _is_excluded(remote_name: str) -> bool:
    """Check if a file should be excluded based on its name."""
    fname = Path(remote_name).name
    if fname in EXCLUDE_PATTERNS:
        return True
    for pat in EXCLUDE_PATTERNS:
        if "*" in pat:
            from fnmatch import fnmatch
            if fnmatch(fname, pat):
                return True
    return False


def collect_files(item_cfg: dict, base_dir: Path) -> list[tuple[Path, str]]:
    """
    Walk source_paths and return [(absolute_path, remote_name), ...].
    remote_name preserves directory structure for subdirectories.
    """
    files: list[tuple[Path, str]] = []
    excluded_count = 0
    for rel, _desc in item_cfg["source_paths"]:
        src = base_dir / rel
        if not src.exists():
            print(f"  ⚠  SKIP: {rel} (not found at {src})")
            continue
        if src.is_file():
            if not _is_excluded(rel):
                files.append((src, rel))
        else:
            for f in sorted(src.rglob("*")):
                if f.is_file() and not f.name.startswith("."):
                    remote = str(f.relative_to(base_dir))
                    if _is_excluded(remote):
                        excluded_count += 1
                        continue
                    files.append((f, remote))
    if excluded_count:
        print(f"  ⚠  Excluded {excluded_count} personal-DNA-derived file(s)")
    return files


def build_metadata(item_cfg: dict) -> dict:
    """Return archive.org metadata dict."""
    return {
        "title": item_cfg["title"],
        "description": item_cfg["description"],
        "subject": item_cfg["subject"],
        "creator": item_cfg["creator"],
        "date": item_cfg["date"],
        "licenseurl": item_cfg["licenseurl"],
        "mediatype": item_cfg["mediatype"],
        "collection": item_cfg["collection"],
    }


# ── upload logic ────────────────────────────────────────────────────────────


def _get_s3_keys() -> tuple[str, str]:
    """Get archive.org S3 access/secret keys from ia config or env."""
    import internetarchive.config
    cfg = internetarchive.config.get_config()
    s3 = cfg.get("s3", {})
    access = s3.get("access") or os.environ.get("IA_ACCESS_KEY")
    secret = s3.get("secret") or os.environ.get("IA_SECRET_KEY")
    if not access or not secret:
        raise RuntimeError(
            "S3 keys not found. Run `ia configure` or set IA_ACCESS_KEY / IA_SECRET_KEY."
        )
    return access, secret


def _remote_md5(item, remote: str) -> str | None:
    """Get the MD5 hash of a remote file on archive.org, or None if not found."""
    try:
        for f in item.files:
            if f.get("name") == remote:
                return f.get("md5")
    except Exception:
        pass
    return None


def _local_md5(fp: Path) -> str:
    """Compute MD5 of a local file."""
    h = hashlib.md5()
    with open(fp, "rb") as f:
        for chunk in iter(lambda: f.read(8_388_608), b""):
            h.update(chunk)
    return h.hexdigest()


def _upload_file(
    item, remote: str, fp: Path, metadata: dict,
    large_threshold: int = 100 * 1024 * 1024,  # 100 MB
) -> tuple[str, bool]:
    """
    Upload a single file. Uses ia CLI for large files (>100 MB) with
    progress bars; uses Python library for small files. Skips if remote
    MD5 matches (verified on archive.org).
    """
    # Check if already uploaded with same content
    remote_md5 = _remote_md5(item, remote)
    if remote_md5:
        local_md5 = _local_md5(fp)
        if local_md5 == remote_md5:
            return (remote, True)

    file_size = fp.stat().st_size
    identifier = item.identifier

    # ── Large files: use ia CLI (robust multipart upload, progress bar) ──
    if file_size > large_threshold:
        try:
            result = subprocess.run(
                ["ia", "upload", "--retries", "5", identifier, str(fp)],
                timeout=86400,  # 24h timeout for very large files
                capture_output=True, text=True,
            )
            ok = result.returncode == 0
            if not ok:
                print(f"    ia CLI stderr: {result.stderr[:200]}")
            return (remote, ok)
        except FileNotFoundError:
            print(f"    ⚠ 'ia' CLI not found. Install: pip install internetarchive")
            print(f"    Falling back to Python library (may fail for large files)...")
        except Exception as e:
            print(f"    ⚠ ia CLI failed: {e}")
            print(f"    Falling back to Python library...")

    # ── Small files (or fallback): use Python library ──
    try:
        access, secret = _get_s3_keys()
        result = item.upload(
            {remote: str(fp)},
            metadata=metadata,
            access_key=access,
            secret_key=secret,
            verbose=False,
            retries=5,
            retries_sleep=30,
            checksum=True,
        )
        ok = result and result[0].status_code in (200, 201)
        return (remote, ok)
    except Exception:
        return (remote, False)


def upload_item(
    identifier: str,
    item_cfg: dict,
    base_dir: Path,
    dry_run: bool = False,
    parallel: int = 4,
):
    """Upload one archive.org item with all its files.

    Strategy:
      1. Small files (≤500 MB) → Python lib, parallel, clean one-line ✓ per file
      2. Large files (>500 MB) → ia CLI, serial, real-time progress bar
    """
    print(f"\n{'='*70}")
    print(f"  Item: {identifier}")
    print(f"  Item:      https://archive.org/details/{identifier}")
    print(f"{'='*70}")

    all_files = collect_files(item_cfg, base_dir)
    # Sort smallest first (quick wins, builds momentum)
    all_files.sort(key=lambda x: x[0].stat().st_size)
    total_bytes = sum(p.stat().st_size for p, _ in all_files)
    print(f"  Files to upload: {len(all_files)}  ({human(total_bytes)})")

    if dry_run:
        print("\n  [DRY RUN] Would upload:")
        for fp, remote in all_files[:10]:
            print(f"    {human(fp.stat().st_size):>10}  {remote}")
        if len(all_files) > 10:
            print(f"    ... and {len(all_files) - 10} more")
        skipped = _find_skippable()
        if skipped:
            print(f"\n  [DRY RUN] Directories intentionally EXCLUDED (personal DNA):")
            for d, size in sorted(skipped.items(), key=lambda x: -x[1]):
                print(f"    {human(size):>10}  {d}/")
        return

    import internetarchive

    # ── Connect ──────────────────────────────────────────────────────────
    print(f"\n  Connecting to archive.org ...")
    try:
        item = internetarchive.get_item(identifier)
    except Exception:
        item = None
    if item is None or item.metadata.get("identifier") != identifier:
        print(f"  Creating new item: {identifier}")
    else:
        print(f"  Item exists: https://archive.org/details/{identifier}")

    metadata = {"title": item_cfg["title"]}
    print_lock = Lock()
    failed: list[str] = []
    uploaded = 0

    print(f"  Uploading with {parallel} workers (checksum skip enabled) ...")
    with ThreadPoolExecutor(max_workers=parallel) as pool:
        futures = {}
        for fp, remote in all_files:
            fut = pool.submit(_upload_file, item, remote, fp, metadata)
            futures[fut] = (remote, fp)

        for fut in as_completed(futures):
            remote, fp = futures[fut]
            size = fp.stat().st_size
            rname, ok = fut.result()
            if ok:
                with print_lock:
                    print(f"  ✓ {human(size):>10}  {rname}")
                uploaded += 1
            else:
                with print_lock:
                    print(f"  ✗ {human(size):>10}  {rname}")
                failed.append(rname)

    # ── Summary ──────────────────────────────────────────────────────────
    print(f"{'='*70}")
    print(f"  ✓ Uploaded: {uploaded}/{len(all_files)}")
    if failed:
        print(f"  ✗ Failed ({len(failed)}):")
        for f in failed:
            print(f"    - {f}")
    print(f"{'='*70}")



# ── exclusion scanner ───────────────────────────────────────────────────────

def _find_skippable() -> dict[str, int]:
    """Report directories intentionally NOT uploaded (personal DNA)."""
    skippable_dirs = [
        ROOT / "raw_data",
        ROOT / "aligned",
        ROOT / "prs_research_pipeline" / "plink",
        ROOT / "prs_research_pipeline" / "pca",
        ROOT / "prs_research_pipeline" / "qc",
        ROOT / "prs_research_pipeline" / "prs",
        ROOT / "prs_research_pipeline" / "raw_data",
        ROOT / "prs_research_pipeline" / "ancestry",
    ]
    result = {}
    for d in skippable_dirs:
        if d.exists():
            total = sum(
                f.stat().st_size
                for f in d.rglob("*")
                if f.is_file() and not f.name.startswith(".")
            )
            if total > 0:
                result[str(d.relative_to(ROOT))] = total
    return result


# ── CLI ─────────────────────────────────────────────────────────────────────

def main():
    import argparse

    parser = argparse.ArgumentParser(
        description="BlueGen — upload public reference data to archive.org",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        "--dry-run", action="store_true",
        help="Validate files and metadata without uploading anything."
    )
    parser.add_argument(
        "--only", choices=("ref", "pgs", "archaic", "vindija", "european"), default=None,
        help="Upload only reference data, PGS cache, archaic panel, Vindija, or European reference."
    )
    parser.add_argument(
        "--parallel", "-j", type=int, default=4,
        help="Number of parallel upload workers (default: 4). Use 1 for serial."
    )
    args = parser.parse_args()

    print("BlueGen Archive.org Upload")
    print(f"  Root: {ROOT}")
    print(f"  Dry run: {args.dry_run}")
    print(f"  Parallel: {args.parallel}")

    if not args.dry_run:
        try:
            import internetarchive  # noqa: F401
        except ImportError:
            print("\n❌  'internetarchive' not installed. Run:")
            print("      pip install internetarchive")
            print("   Then configure credentials:")
            print("      ia configure")
            print("   Or set env vars: IA_ACCESS_KEY + IA_SECRET_KEY")
            sys.exit(1)

    for label, dir_path in [
        ("Reference", REFERENCE_DIR),
        ("PGS", PGS_DIR),
    ]:
        if args.only and label.lower()[:3] != args.only[:3]:
            continue
        if not dir_path.exists():
            print(f"\n❌  {label} directory not found: {dir_path}")
            sys.exit(1)

    start = time.time()

    if args.only is None or args.only == "ref":
        upload_item(
            "bluegen-reference-data",
            ITEMS["bluegen-reference-data"],
            REFERENCE_DIR,
            dry_run=args.dry_run,
            parallel=args.parallel,
        )

    if args.only is None or args.only == "european":
        upload_item(
            "bluegen-european-reference",
            ITEMS["bluegen-european-reference"],
            REFERENCE_DIR / "human_origins",
            dry_run=args.dry_run,
            parallel=args.parallel,
        )

    if args.only is None or args.only == "pgs":
        upload_item(
            "bluegen-pgs-cache",
            ITEMS["bluegen-pgs-cache"],
            PGS_DIR,
            dry_run=args.dry_run,
            parallel=args.parallel,
        )

    if args.only is None or args.only == "archaic":
        upload_item(
            "bluegen-archaic-reference",
            ITEMS["bluegen-archaic-reference"],
            REFERENCE_DIR,
            dry_run=args.dry_run,
            parallel=args.parallel,
        )

    if args.only is None or args.only == "vindija":
        upload_item(
            "bluegen-vindija-reference",
            ITEMS["bluegen-vindija-reference"],
            REFERENCE_DIR,
            dry_run=args.dry_run,
            parallel=args.parallel,
        )

    elapsed = time.time() - start
    print(f"\n{'='*70}")
    print(f"  Done in {elapsed:.0f}s.")
    if not args.dry_run:
        print(f"  Reference:  https://archive.org/details/bluegen-reference-data")
        print(f"  PGS cache:  https://archive.org/details/bluegen-pgs-cache")
        print(f"  Archaic:    https://archive.org/details/bluegen-archaic-reference")
        print(f"  Vindija:    https://archive.org/details/bluegen-vindija-reference")
    print(f"{'='*70}\n")


if __name__ == "__main__":
    main()
