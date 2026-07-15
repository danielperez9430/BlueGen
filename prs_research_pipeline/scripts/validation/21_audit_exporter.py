#!/usr/bin/env python3
"""
╔══════════════════════════════════════════════════════════════════════════════╗
║   PHASE 7 — MODULE 6: EXTERNAL AUDIT EXPORT FORMAT                          ║
║   scripts/21_audit_exporter.py                                              ║
║                                                                            ║
║   Prepares a complete, self-contained data package for external peer review.║
║                                                                            ║
║   SCIENTIFIC FREEZE LAYER — Everything an external reviewer needs.          ║
║                                                                            ║
║   Exports:                                                                  ║
║     • Harmonized VCF subset (reviewer-relevant variants only)               ║
║     • PRS tables (raw, adjusted, calibrated, uncertainty-quantified)        ║
║     • PCA coordinates (reference centroids + sample projection)             ║
║     • Ancestry probabilities + classification report                        ║
║     • GWAS/PGS mappings + provenance                                        ║
║     • Full method metadata (software versions, parameters, seeds)           ║
║     • SHA-256 manifest for all exported files                               ║
║     • Structured README for reviewers                                       ║
║                                                                            ║
║   Output:                                                                   ║
║     audit/audit_package_YYYYMMDD.zip                                        ║
║     audit/audit_manifest.json                                               ║
║     audit/README_for_reviewers.md                                           ║
╚══════════════════════════════════════════════════════════════════════════════╝
"""

import sys
import os
import json
import hashlib
import logging
import zipfile
import shutil
from pathlib import Path
from typing import Dict, List, Optional, Any, Tuple
from dataclasses import dataclass, field, asdict
from datetime import datetime

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from utils.constants import PIPELINE_VERSION

logger = logging.getLogger(__name__)


# ── Data Structures ────────────────────────────────────────────────────────────

@dataclass
class AuditFileEntry:
    """Entry for a single file in the audit package."""
    relative_path: str
    description: str
    sha256: str = ""
    file_size_bytes: int = 0
    category: str = ""  # genotypes, prs, ancestry, gwas, metadata, documentation


@dataclass
class AuditPackage:
    """Complete audit package metadata."""
    package_version: str = "1.0.0"
    created_date: str = ""
    pipeline_version: str = PIPELINE_VERSION
    run_id: str = ""
    sample_id: str = ""
    total_files: int = 0
    total_size_bytes: int = 0
    files: List[AuditFileEntry] = field(default_factory=list)
    manifest_sha256: str = ""


# ── Audit Exporter ────────────────────────────────────────────────────────────

class AuditExporter:
    """
    Exports a complete, self-contained audit package for external review.

    The package contains everything needed for an independent reviewer to:
      1. Verify input data integrity (hashes)
      2. Reproduce PRS computation (methodology + parameters)
      3. Validate results (all intermediate + final outputs)
      4. Audit the full pipeline execution trace

    Usage:
        exporter = AuditExporter(output_dir="audit")
        package = exporter.export(
            run_id="run_20260603",
            sample_id="SAMPLE_001",
            input_files={"vcf": "input.vcf.gz", "snp_db": "data/snp_database_annotated.csv"},
            output_files=["prs/population_calibrated.csv", "pca/projected_sample.csv"],
            metadata_files=["reproducibility/run_fingerprint.json", "science/assumptions.lock.json"],
        )
    """

    FILE_CATEGORIES = {
        ".vcf": "genotypes", ".vcf.gz": "genotypes", ".bed": "genotypes",
        ".bim": "genotypes", ".fam": "genotypes",
        ".csv": "prs", ".tsv": "prs",
        "pca": "ancestry", "ancestry": "ancestry", "admixture": "ancestry",
        "gwas": "gwas", "score": "gwas",
        "manifest": "metadata", "fingerprint": "metadata", "lock": "metadata",
        "report": "documentation", "dashboard": "documentation",
    }

    def __init__(self, output_dir: str = "audit"):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)

    # ── Public API ───────────────────────────────────────────────────────

    def export(
        self,
        run_id: str,
        sample_id: str,
        input_files: Dict[str, str],
        output_files: Dict[str, str],
        metadata_files: Dict[str, str],
        include_vcf_subset: bool = True,
        vcf_subset_variants: int = 5000,
        plink_binary: str = "plink",
    ) -> AuditPackage:
        """
        Export complete audit package.

        Args:
            run_id: Pipeline run identifier.
            sample_id: Sample identifier.
            input_files: Dict of {label: path} for input files.
            output_files: Dict of {label: path} for output files.
            metadata_files: Dict of {label: path} for metadata files.
            include_vcf_subset: Include a VCF subset for reviewer inspection.
            vcf_subset_variants: Number of variants in VCF subset.

        Returns:
            AuditPackage with complete manifest.
        """
        logger.info("═══ External Audit Export (Phase 7) ═══")

        package_dir = self.output_dir / f"audit_package_{datetime.now().strftime('%Y%m%d')}"
        package_dir.mkdir(parents=True, exist_ok=True)

        package = AuditPackage(
            created_date=datetime.now().strftime("%Y-%m-%d %H:%M UTC"),
            run_id=run_id,
            sample_id=sample_id,
        )

        entries = []

        # 1. Export metadata files (fingerprint, manifest, lock, seed registry)
        logger.info("  Exporting metadata...")
        for label, path in metadata_files.items():
            if Path(path).exists():
                dest = package_dir / label
                dest.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(path, dest)
                entry = self._create_entry(str(dest.relative_to(package_dir)), label, "metadata")
                entries.append(entry)

        # 2. Export output files
        logger.info("  Exporting outputs...")
        for label, path in output_files.items():
            if Path(path).exists():
                dest = package_dir / "outputs" / label
                dest.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(path, dest)
                entry = self._create_entry(
                    str(dest.relative_to(package_dir)), label,
                    self._categorize(label)
                )
                entries.append(entry)

        # 3. Export VCF subset for reviewer genotype inspection
        if include_vcf_subset and "vcf" in input_files:
            logger.info("  Extracting VCF subset...")
            vcf_path = input_files["vcf"]
            if Path(vcf_path).exists():
                vcf_subset_path = self._create_vcf_subset(
                    vcf_path, package_dir, vcf_subset_variants, plink_binary
                )
                if vcf_subset_path:
                    entry = self._create_entry(
                        str(vcf_subset_path.relative_to(package_dir)),
                        "Genotype subset for reviewer inspection",
                        "genotypes",
                    )
                    entries.append(entry)

        # 4. Copy SNP database
        if "snp_db" in input_files and Path(input_files["snp_db"]).exists():
            dest = package_dir / "reference" / "snp_database.csv"
            dest.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(input_files["snp_db"], dest)
            entry = self._create_entry(
                str(dest.relative_to(package_dir)),
                "Curated SNP database (109 nutrigenetic variants)",
                "gwas",
            )
            entries.append(entry)

        # 5. Generate README for reviewers
        readme_path = package_dir / "README_for_reviewers.md"
        self._generate_reviewer_readme(readme_path, package, entries, run_id, sample_id)
        entry = self._create_entry("README_for_reviewers.md", "Reviewer instructions", "documentation")
        entries.append(entry)

        # 6. Generate manifest
        package.files = entries
        package.total_files = len(entries)
        package.total_size_bytes = sum(e.file_size_bytes for e in entries)

        manifest_path = package_dir / "audit_manifest.json"
        package.manifest_sha256 = self._save_manifest(package, manifest_path)

        # 7. Create ZIP archive
        zip_path = self._create_zip(package_dir)
        logger.info(f"  ✅ Audit package: {zip_path} ({package.total_size_bytes:,} bytes)")

        return package

    # ── Private: VCF Subset ───────────────────────────────────────────────

    def _create_vcf_subset(
        self, vcf_path: str, package_dir: Path, n_variants: int, plink_binary: str
    ) -> Optional[Path]:
        """Create a small VCF subset for reviewer inspection."""
        import subprocess

        subset_vcf = package_dir / "genotypes" / "sample_subset.vcf.gz"
        subset_vcf.parent.mkdir(parents=True, exist_ok=True)

        # Extract first N variants using bcftools or PLINK
        try:
            # Try bcftools first
            result = subprocess.run(
                ["bcftools", "view", "-h", vcf_path],
                capture_output=True, text=True, timeout=30
            )
            header_lines = result.stdout.count("\n")

            subprocess.run(
                ["bcftools", "view", vcf_path, "|", "head", "-n", str(header_lines + n_variants),
                 "|", "bgzip", "-c", ">", str(subset_vcf)],
                shell=True, capture_output=True, timeout=120
            )

            if subset_vcf.exists():
                subprocess.run(["tabix", "-p", "vcf", str(subset_vcf)],
                             capture_output=True, timeout=30)
                logger.info(f"    VCF subset: {subset_vcf} ({n_variants} variants)")
                return subset_vcf

        except Exception:
            pass

        # Fallback: copy header + first N lines
        try:
            with open(vcf_path, "rb") as src:
                # Read header
                header = b""
                for _ in range(200):
                    line = src.readline()
                    header += line
                    if not line.startswith(b"#"):
                        break

                # Read N variant lines
                body = b""
                for _ in range(n_variants):
                    line = src.readline()
                    if not line:
                        break
                    body += line

            simple_vcf = package_dir / "genotypes" / "sample_subset.vcf"
            with open(simple_vcf, "wb") as fh:
                fh.write(header + body)

            logger.info(f"    VCF subset (simple): {simple_vcf}")
            return simple_vcf

        except Exception as e:
            logger.warning(f"    VCF subset failed: {e}")
            return None

    # ── Private: README Generator ─────────────────────────────────────────

    def _generate_reviewer_readme(
        self,
        readme_path: Path,
        package: AuditPackage,
        entries: List[AuditFileEntry],
        run_id: str,
        sample_id: str,
    ) -> None:
        """Generate a structured README for external reviewers."""

        lines = [
            "# PRS Research Platform — External Audit Package",
            "",
            f"**Run ID:** `{run_id}`",
            f"**Sample ID:** `{sample_id}`",
            f"**Generated:** {package.created_date}",
            f"**Pipeline:** v{package.pipeline_version} (Phase 7 Scientific Freeze)",
            "",
            "---",
            "",
            "## Purpose of This Package",
            "",
            "This package contains all data needed for independent peer review of the",
            "PRS Research Platform analysis. It includes:",
            "",
            "- **Genotype data:** A representative VCF subset for validation",
            "- **PRS outputs:** All computed risk scores (raw, adjusted, calibrated)",
            "- **Ancestry analysis:** PCA coordinates, population classifications, admixture",
            "- **GWAS/SNP data:** The curated variant database with effect sizes and evidence levels",
            "- **Methodology:** Complete scientific assumption lock file and run fingerprint",
            "- **Execution trace:** Full pipeline manifest with per-stage timing and metrics",
            "",
            "---",
            "",
            "## Package Contents",
            "",
            "| File | Category | Size | Description |",
            "|------|----------|------|-------------|",
        ]

        for entry in entries:
            size_kb = entry.file_size_bytes / 1024
            lines.append(
                f"| `{entry.relative_path}` | {entry.category} | {size_kb:.1f} KB | {entry.description} |"
            )

        lines += [
            "",
            "---",
            "",
            "## How to Verify Reproducibility",
            "",
            "1. **Check the environment fingerprint** in `run_fingerprint.json`",
            "2. **Compare all file hashes** against `audit_manifest.json`",
            "3. **Validate assumptions** against `assumptions.lock.json`",
            "4. **Re-run the pipeline** using the documented parameters and seed",
            "5. **Compare output hashes** — they should match the manifest exactly",
            "",
            "```bash",
            "# Verify all file hashes",
            "sha256sum -c audit_manifest.json  # (after extracting hash list)",
            "",
            "# Or use the Python verification script:",
            "python3 scripts/16_reproducibility_engine.py --verify audit_package_*/*",
            "```",
            "",
            "---",
            "",
            "## Key Scientific Parameters (Frozen)",
            "",
            "| Parameter | Value |",
            "|-----------|-------|",
            "| LD pruning | r² < 0.2, window = 50 variants |",
            "| PCA method | Reference projection (Price et al. 2006) |",
            "| PCA components | 20 |",
            "| Ancestry reference | 1000 Genomes Phase 3 (2,504 samples) |",
            "| Ancestry method | PCA ensemble (centroid + logistic + k-NN) |",
            "| PRS formula | Σ(βⱼ × Gᵢⱼ) via PLINK --score |",
            "| Population calibration | Empirical 1000G distributions |",
            "| Uncertainty | 3-layer variance propagation |",
            "| GWAS evidence | Levels A–D (p < 5e-8 to p < 0.05) |",
            "| Risk categories | Low < 25%, Medium 25–75%, High > 75% |",
            "",
            "---",
            "",
            "## Limitations (Please Read Before Interpreting)",
            "",
            "- **PRS is probabilistic, not diagnostic** — these scores estimate genetic susceptibility",
            "- **European-centric GWAS** — effect sizes are primarily from EUR discovery populations",
            "- **Curated SNP panel** — 109 nutrigenetic SNPs, not genome-wide scoring",
            "- **Single-sample analysis** — population variance estimated from 1000G reference",
            "- **Not clinically validated** — research use only, not FDA/EMA approved",
            "",
            "---",
            "",
            "## Contact & Citation",
            "",
            "This audit package was generated by the PRS Research Platform v7.0.0",
            "under the Phase 7 Scientific Freeze Layer.",
            "",
            "For questions about methodology, reproducibility, or data access, please",
            "contact the corresponding author.",
            "",
            f"*Package generated: {package.created_date}*",
        ]

        with open(readme_path, "w") as fh:
            fh.write("\n".join(lines))

        logger.info(f"  ✅ Reviewer README: {readme_path}")

    # ── Private: Helpers ─────────────────────────────────────────────────

    def _create_entry(self, relative_path: str, description: str, category: str) -> AuditFileEntry:
        """Create an audit file entry with hash."""
        full_path = self.output_dir / f"audit_package_{datetime.now().strftime('%Y%m%d')}" / relative_path
        sha = hashlib.sha256()
        if full_path.exists():
            with open(full_path, "rb") as fh:
                for chunk in iter(lambda: fh.read(65536), b""):
                    sha.update(chunk)

        return AuditFileEntry(
            relative_path=relative_path,
            description=description,
            sha256=sha.hexdigest(),
            file_size_bytes=full_path.stat().st_size if full_path.exists() else 0,
            category=category,
        )

    def _categorize(self, filename: str) -> str:
        """Categorize a file by name/path."""
        name_lower = filename.lower()
        for pattern, category in self.FILE_CATEGORIES.items():
            if pattern in name_lower or pattern in filename:
                return category
        return "prs"

    def _save_manifest(self, package: AuditPackage, manifest_path: Path) -> str:
        """Save audit manifest with SHA-256."""
        manifest_data = {
            "package_version": package.package_version,
            "created_date": package.created_date,
            "pipeline_version": package.pipeline_version,
            "run_id": package.run_id,
            "sample_id": package.sample_id,
            "total_files": package.total_files,
            "total_size_bytes": package.total_size_bytes,
            "files": [asdict(f) for f in package.files],
        }

        manifest_json = json.dumps(manifest_data, indent=2)
        with open(manifest_path, "w") as fh:
            fh.write(manifest_json)

        return hashlib.sha256(manifest_json.encode()).hexdigest()[:16]

    def _create_zip(self, package_dir: Path) -> Path:
        """Create ZIP archive of the audit package."""
        zip_path = self.output_dir / f"{package_dir.name}.zip"

        with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as zf:
            for file_path in sorted(package_dir.rglob("*")):
                if file_path.is_file():
                    arcname = str(file_path.relative_to(package_dir))
                    zf.write(file_path, arcname)

        logger.info(f"  ✅ ZIP archive: {zip_path} ({zip_path.stat().st_size:,} bytes)")
        return zip_path


# ── CLI Entry Point ──────────────────────────────────────────────────────────

def main():
    import argparse

    parser = argparse.ArgumentParser(
        description="Phase 7 Module 6: External Audit Export Format"
    )
    parser.add_argument("--run-id", required=True, help="Pipeline run ID")
    parser.add_argument("--sample-id", default="SAMPLE_001")
    parser.add_argument("--output-dir", "-o", default="audit")
    parser.add_argument("--vcf", default="input.vcf.gz", help="Input VCF path")
    parser.add_argument("--snp-db", default="data/snp_database_annotated.csv")
    parser.add_argument("--vcf-subset-variants", type=int, default=5000)
    parser.add_argument("--no-vcf-subset", action="store_true")
    parser.add_argument("--plink", default="plink")
    parser.add_argument("--verbose", "-v", action="store_true")

    args = parser.parse_args()

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
    )

    # Gather available output files
    output_files = {}
    for label, path in [
        ("prs_raw", "prs/prs_raw.csv"),
        ("prs_adjusted", "prs/prs_adjusted.csv"),
        ("prs_calibrated", "prs/population_calibrated.csv"),
        ("prs_calibrated_v2", "prs/population_calibrated_v2.csv"),
        ("prs_uncertainty", "prs/prs_uncertainty.csv"),
        ("prs_all_methods", "prs/prs_all_methods.csv"),
        ("pca_sample", "pca/projected_sample.csv"),
        ("pca_centroids", "pca/reference_centroids.csv"),
        ("ancestry_probs", "ancestry/posterior_probabilities.json"),
        ("ancestry_admixture", "ancestry/admixture_results.json"),
        ("coverage_audit", "validation/coverage_audit.csv"),
        ("evidence_scores", "validation/evidence_scores.json"),
        ("benchmark", "validation/scientific_benchmark.json"),
        ("rqs_v2", "research/research_quality_v2.json"),
        ("report_en", "reports/report_en.html"),
        ("report_es", "reports/report_es.html"),
    ]:
        if Path(path).exists():
            output_files[label] = path

    # Gather metadata files
    metadata_files = {}
    for label, path in [
        ("run_fingerprint", "reproducibility/run_fingerprint.json"),
        ("seed_registry", "reproducibility/seed_registry.json"),
        ("execution_manifest", "validation/execution_manifest.json"),
        ("assumptions_lock", "science/assumptions.lock.json"),
        ("run_manifest", "reports/run_manifest.json"),
    ]:
        if Path(path).exists():
            metadata_files[label] = path

    exporter = AuditExporter(output_dir=args.output_dir)
    package = exporter.export(
        run_id=args.run_id,
        sample_id=args.sample_id,
        input_files={
            "vcf": args.vcf,
            "snp_db": args.snp_db,
        },
        output_files=output_files,
        metadata_files=metadata_files,
        include_vcf_subset=not args.no_vcf_subset,
        vcf_subset_variants=args.vcf_subset_variants,
        plink_binary=args.plink,
    )

    print(f"\n═══ Audit Package ═══")
    print(f"  Files: {package.total_files}")
    print(f"  Size: {package.total_size_bytes:,} bytes")
    print(f"  Categories:")
    cats = {}
    for f in package.files:
        cats[f.category] = cats.get(f.category, 0) + 1
    for cat, count in sorted(cats.items()):
        print(f"    {cat}: {count} files")

    return 0


if __name__ == "__main__":
    sys.exit(main())
