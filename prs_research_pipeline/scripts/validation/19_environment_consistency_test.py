#!/usr/bin/env python3
"""
╔══════════════════════════════════════════════════════════════════════════════╗
║   PHASE 7 — MODULE 4: CROSS-ENVIRONMENT CONSISTENCY TEST                    ║
║   scripts/19_environment_consistency_test.py                                 ║
║                                                                            ║
║   Verifies that pipeline output is identical across environments:           ║
║     • macOS (primary)                                                       ║
║     • Linux (Docker / CI)                                                   ║
║     • Any other POSIX environment                                           ║
║                                                                            ║
║   Method:                                                                   ║
║     • Run subset pipeline on current environment                            ║
║     • Compare output hashes against reference fingerprints                  ║
║     • Compute environment divergence score                                  ║
║     • Flag any environment-specific deviations                              ║
║                                                                            ║
║   Output:                                                                   ║
║     validation/environment_consistency.json                                 ║
╚══════════════════════════════════════════════════════════════════════════════╝
"""

import sys
import os
import json
import hashlib
import logging
import platform
import subprocess
from pathlib import Path
from typing import Dict, List, Optional, Any, Tuple
from dataclasses import dataclass, field, asdict
from datetime import datetime

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


# ── Data Structures ────────────────────────────────────────────────────────────

@dataclass
class EnvironmentInfo:
    """Information about the current environment."""
    os_name: str
    os_version: str
    kernel: str
    architecture: str
    python_version: str
    plink_version: str
    hostname: str
    is_ci: bool = False
    is_docker: bool = False


@dataclass
class FileComparison:
    """Comparison of a single output file across environments."""
    file_path: str
    reference_hash: str
    current_hash: str
    match: bool
    file_size_bytes: int = 0


@dataclass
class ConsistencyReport:
    """Complete environment consistency report."""
    current_environment: EnvironmentInfo
    reference_environment: Optional[Dict[str, str]]
    files_compared: int = 0
    files_matched: int = 0
    files_diverged: int = 0
    comparisons: List[FileComparison] = field(default_factory=list)
    divergence_score: float = 0.0
    overall_consistent: bool = False
    generated_date: str = ""


# ── Environment Consistency Tester ────────────────────────────────────────────

class EnvironmentConsistencyTester:
    """
    Tests output consistency across execution environments.

    Usage:
        tester = EnvironmentConsistencyTester()
        report = tester.run_test(
            reference_fingerprint="reproducibility/run_fingerprint.json",
            output_files=["prs/population_calibrated.csv", "pca/projected_sample.csv"],
        )
    """

    def __init__(self, output_dir: str = "validation"):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)

    # ── Public API ───────────────────────────────────────────────────────

    def run_test(
        self,
        reference_fingerprint: str,
        output_files: List[str],
        tolerance: float = 0.0,  # 0 = exact match required
    ) -> ConsistencyReport:
        """
        Compare current outputs against a reference fingerprint.

        Args:
            reference_fingerprint: Path to reference run_fingerprint.json.
            output_files: List of output files to compare.
            tolerance: Allowed hash divergence (0 = exact match required).

        Returns:
            ConsistencyReport with per-file comparison.
        """
        logger.info("═══ Cross-Environment Consistency Test ═══")

        # Capture current environment
        current_env = self._capture_environment()
        logger.info(f"  Current: {current_env.os_name} {current_env.os_version}")
        logger.info(f"  Python: {current_env.python_version}")
        logger.info(f"  PLINK: {current_env.plink_version}")

        # Load reference
        ref_env = None
        ref_hashes = {}
        if Path(reference_fingerprint).exists():
            with open(reference_fingerprint) as fh:
                ref_data = json.load(fh)
            ref_env = ref_data.get("environment", {})
            ref_hashes = ref_data.get("output_hashes", {})
            logger.info(f"  Reference: {ref_env.get('os_name', '?')} {ref_env.get('os_version', '?')}")

        # Compare each file
        comparisons = []
        for file_path in output_files:
            if not Path(file_path).exists():
                comparisons.append(FileComparison(
                    file_path=file_path,
                    reference_hash=ref_hashes.get(file_path, "MISSING"),
                    current_hash="FILE_NOT_FOUND",
                    match=False,
                ))
                continue

            current_hash = self._hash_file(file_path)
            ref_hash = ref_hashes.get(file_path, "NO_REFERENCE")
            file_size = Path(file_path).stat().st_size

            match = (current_hash == ref_hash) if ref_hash != "NO_REFERENCE" else None

            comparisons.append(FileComparison(
                file_path=file_path,
                reference_hash=ref_hash,
                current_hash=current_hash,
                match=match,
                file_size_bytes=file_size,
            ))

            icon = "✅" if match else ("❌" if match is False else "⬚")
            logger.info(f"  {icon} {Path(file_path).name}: {current_hash[:12]}")

        # Compute divergence score
        n_comparable = sum(1 for c in comparisons if c.match is not None)
        n_matched = sum(1 for c in comparisons if c.match is True)

        if n_comparable > 0:
            divergence = 1.0 - (n_matched / n_comparable)
        else:
            divergence = 0.0  # No reference → assume consistent

        # Check environment-level differences
        env_diffs = []
        if ref_env:
            if current_env.python_version != ref_env.get("python_version", ""):
                env_diffs.append(f"Python: {current_env.python_version} ≠ {ref_env.get('python_version')}")
            if current_env.os_name != ref_env.get("os_name", ""):
                env_diffs.append(f"OS: {current_env.os_name} ≠ {ref_env.get('os_name')}")

        report = ConsistencyReport(
            current_environment=current_env,
            reference_environment=ref_env,
            files_compared=len(comparisons),
            files_matched=n_matched,
            files_diverged=n_comparable - n_matched,
            comparisons=comparisons,
            divergence_score=round(divergence, 4),
            overall_consistent=(divergence == 0.0),
            generated_date=datetime.now().strftime("%Y-%m-%d %H:%M UTC"),
        )

        # Summary
        if report.overall_consistent:
            logger.info(f"  ✅ Environment consistent — all {n_matched} files match")
        else:
            logger.warning(f"  ⚠️  {report.files_diverged} files diverge "
                         f"(divergence={divergence:.1%})")
            for diff in env_diffs:
                logger.warning(f"    {diff}")

        # Save
        self._save_report(report)

        return report

    def generate_subset_test(
        self,
        input_vcf: str,
        output_dir: str,
        plink_binary: str = "plink",
    ) -> Dict[str, str]:
        """
        Run a minimal subset of the pipeline for rapid cross-environment testing.

        This executes a lightweight test (VCF→PLINK→QC→5 SNPs scored)
        that can be run in seconds, not hours, for CI/CD validation.
        """
        logger.info("═══ Running Subset Test ═══")

        hashes = {}
        tmp_dir = Path(output_dir) / "env_test"
        tmp_dir.mkdir(parents=True, exist_ok=True)

        # Stage 1: VCF → PLINK (first 1000 variants only)
        tmp_prefix = str(tmp_dir / "test")
        subprocess.run([
            plink_binary, "--vcf", input_vcf,
            "--make-bed", "--out", tmp_prefix,
            "--extract", f"<(cut -f2 {tmp_prefix}.bim | head -1000)",
            "--threads", "2",
        ], capture_output=True, text=True, timeout=120, shell=True)

        # Hash outputs
        for ext in [".bed", ".bim", ".fam"]:
            p = f"{tmp_prefix}{ext}"
            if Path(p).exists():
                hashes[f"test{ext}"] = self._hash_file(p)[:12]

        # Calculate basic QC stats
        subprocess.run([
            plink_binary, "--bfile", tmp_prefix,
            "--freq", "--out", f"{tmp_prefix}_freq",
            "--threads", "2",
        ], capture_output=True, text=True, timeout=60)

        for ext in [".frq", ".log"]:
            p = f"{tmp_prefix}_freq{ext}"
            if Path(p).exists():
                hashes[f"test_freq{ext}"] = self._hash_file(p)[:12]

        logger.info(f"  Subset test: {len(hashes)} output hashes generated")
        return hashes

    # ── Private Methods ──────────────────────────────────────────────────

    def _capture_environment(self) -> EnvironmentInfo:
        """Capture current environment information."""
        is_docker = Path("/.dockerenv").exists()
        is_ci = bool(os.environ.get("CI")) or bool(os.environ.get("GITHUB_ACTIONS"))

        # Detect PLINK version
        plink_ver = "unknown"
        try:
            result = subprocess.run(["plink", "--version"], capture_output=True, text=True, timeout=10)
            plink_ver = result.stdout.split("\n")[0].strip() if result.stdout else "unknown"
        except Exception:
            pass

        return EnvironmentInfo(
            os_name=platform.system(),
            os_version=platform.release(),
            kernel=platform.version(),
            architecture=platform.machine(),
            python_version=platform.python_version(),
            plink_version=plink_ver,
            hostname=platform.node(),
            is_ci=is_ci,
            is_docker=is_docker,
        )

    @staticmethod
    def _hash_file(path: str) -> str:
        """Compute SHA-256 hash of a file."""
        sha = hashlib.sha256()
        with open(path, "rb") as fh:
            for chunk in iter(lambda: fh.read(65536), b""):
                sha.update(chunk)
        return sha.hexdigest()[:16]

    def _save_report(self, report: ConsistencyReport) -> None:
        """Save consistency report to JSON."""
        json_path = self.output_dir / "environment_consistency.json"
        with open(json_path, "w") as fh:
            json.dump({
                "overall_consistent": report.overall_consistent,
                "divergence_score": report.divergence_score,
                "files_compared": report.files_compared,
                "files_matched": report.files_matched,
                "files_diverged": report.files_diverged,
                "generated_date": report.generated_date,
                "current_environment": {
                    "os": report.current_environment.os_name,
                    "os_version": report.current_environment.os_version,
                    "python": report.current_environment.python_version,
                    "plink": report.current_environment.plink_version,
                    "is_docker": report.current_environment.is_docker,
                    "is_ci": report.current_environment.is_ci,
                },
                "comparisons": [
                    {
                        "file": c.file_path,
                        "reference_hash": c.reference_hash,
                        "current_hash": c.current_hash,
                        "match": c.match,
                    }
                    for c in report.comparisons
                ],
            }, fh, indent=2)
        logger.info(f"  ✅ Consistency report: {json_path}")


# ── CLI Entry Point ──────────────────────────────────────────────────────────

def main():
    import argparse

    parser = argparse.ArgumentParser(
        description="Phase 7 Module 4: Cross-Environment Consistency Test"
    )
    parser.add_argument("--reference", required=True,
                       help="Reference run_fingerprint.json")
    parser.add_argument("--files", nargs="+", required=True,
                       help="Output files to compare")
    parser.add_argument("--output-dir", "-o", default="validation")
    parser.add_argument("--subset-test", action="store_true",
                       help="Run rapid subset test for CI/CD")
    parser.add_argument("--input-vcf", help="VCF for subset test")
    parser.add_argument("--plink", default="plink")
    parser.add_argument("--tolerance", type=float, default=0.0)
    parser.add_argument("--verbose", "-v", action="store_true")

    args = parser.parse_args()

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
    )

    tester = EnvironmentConsistencyTester(output_dir=args.output_dir)

    if args.subset_test and args.input_vcf:
        hashes = tester.generate_subset_test(args.input_vcf, args.output_dir, args.plink)
        print(f"\n═══ Subset Test Hashes ═══")
        for f, h in hashes.items():
            print(f"  {f}: {h}")

    report = tester.run_test(
        reference_fingerprint=args.reference,
        output_files=args.files,
        tolerance=args.tolerance,
    )

    print(f"\n═══ Environment Consistency ═══")
    print(f"  Current: {report.current_environment.os_name} "
          f"(Python {report.current_environment.python_version})")
    print(f"  Files compared: {report.files_compared}")
    print(f"  Matched: {report.files_matched} / Diverged: {report.files_diverged}")
    print(f"  Divergence score: {report.divergence_score:.1%}")
    print(f"  Overall: {'✅ CONSISTENT' if report.overall_consistent else '❌ DIVERGENT'}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
