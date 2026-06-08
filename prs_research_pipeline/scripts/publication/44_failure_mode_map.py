#!/usr/bin/env python3
"""
╔══════════════════════════════════════════════════════════════════════════════╗
║   PHASE 10 — SCIENTIFIC FAILURE MODE MAP                                     ║
║   scripts/44_failure_mode_map.py                                             ║
║                                                                            ║
║   Maps every failure mode of the platform:                                  ║
║     • Where PRS becomes invalid                                              ║
║     • Where ancestry inference fails                                         ║
║     • Where calibration breaks                                               ║
║     • Where benchmarking becomes unreliable                                  ║
║                                                                            ║
║   Each failure includes: threshold, severity, biological/statistical reason ║
║                                                                            ║
║   Output:                                                                    ║
║     science/failure_mode_map.json                                            ║
║     science/failure_mode_map.md                                              ║
╚══════════════════════════════════════════════════════════════════════════════╝
"""

import sys, os, json, logging
from pathlib import Path
from typing import Dict, List
from dataclasses import dataclass, field, asdict
from datetime import datetime

logger = logging.getLogger(__name__)

FAILURE_MODES = [
    # --- PRS FAILURES ---
    {"id": "FM-PRS-001", "component": "PRS computation", "failure": "All SNPs missing from sample",
     "threshold": "coverage < 5%", "severity": "CRITICAL",
     "cause": "No variants in common between genotype data and PRS model.",
     "effect": "PRS = 0 for all traits. Complete failure of scoring."},
    {"id": "FM-PRS-002", "component": "PRS computation", "failure": "Effect allele strand ambiguity unresolved",
     "threshold": "> 30% A/T or C/G SNPs", "severity": "HIGH",
     "cause": "Palindromic SNPs (A/T, C/G) cannot be resolved without LD reference.",
     "effect": "Randomized effect directions for ambiguous SNPs — score variance inflated."},
    {"id": "FM-PRS-003", "component": "PRS computation", "failure": "Absolute beta summation error",
     "threshold": "PRS = Σ|β| detected", "severity": "CRITICAL",
     "cause": "Incorrect formula: sum of absolute betas instead of dosage-weighted sum.",
     "effect": "PRS becomes a constant — no individual variation captured."},
    {"id": "FM-PRS-004", "component": "PRS computation", "failure": "LD pruning removes all trait SNPs",
     "threshold": "0 SNPs retained after pruning", "severity": "HIGH",
     "cause": "r² threshold too aggressive for the trait's genomic architecture.",
     "effect": "Zero SNPs available for scoring — PRS undefined."},
    {"id": "FM-PRS-005", "component": "PRS normalization", "failure": "Synthetic μ=0 calibration",
     "threshold": "population_mu = 0 for all populations", "severity": "CRITICAL",
     "cause": "Reference distributions not computed from actual data.",
     "effect": "Population calibration is meaningless — cross-ancestry comparisons invalid."},

    # --- ANCESTRY FAILURES ---
    {"id": "FM-ANC-001", "component": "Ancestry inference", "failure": "Insufficient PCA dimensions",
     "threshold": "< 5 PCs used", "severity": "HIGH",
     "cause": "PCA does not capture enough variance to separate populations.",
     "effect": "Population clusters overlap — classification confidence degraded."},
    {"id": "FM-ANC-002", "component": "Ancestry inference", "failure": "Chr22-only reference bias",
     "threshold": "> 95% of variants from chr22", "severity": "CRITICAL",
     "cause": "Chromosome 22 is ~1% of genome — cannot capture genome-wide ancestry.",
     "effect": "Ancestry misclassification for non-EUR populations. Admixture undetectable."},
    {"id": "FM-ANC-003", "component": "Ancestry inference", "failure": "Trait-SNP-based inference",
     "threshold": "< 100 SNPs used", "severity": "CRITICAL",
     "cause": "Trait-associated SNPs have distorted allele frequencies.",
     "effect": "Ancestry call reflects trait biology, not population genetics."},
    {"id": "FM-ANC-004", "component": "Ancestry inference", "failure": "Reference panel mismatch",
     "threshold": "Sample ancestry not in reference populations", "severity": "HIGH",
     "cause": "Reference panel does not contain the sample's ancestral population.",
     "effect": "Sample forced to nearest population — admixture fractions unreliable."},

    # --- CALIBRATION FAILURES ---
    {"id": "FM-CAL-001", "component": "Population calibration", "failure": "Reference sample size too small",
     "threshold": "< 30 samples per population", "severity": "HIGH",
     "cause": "Insufficient reference data to estimate population parameters.",
     "effect": "μ_pop and σ_pop estimates have large uncertainty — z-scores unreliable."},
    {"id": "FM-CAL-002", "component": "Population calibration", "failure": "Target population different from all references",
     "threshold": "Max ancestry probability < 0.50", "severity": "CRITICAL",
     "cause": "No single population adequately represents the sample.",
     "effect": "Population-specific calibration not applicable — use global percentiles only."},
    {"id": "FM-CAL-003", "component": "Population calibration", "failure": "Calibration distribution non-normal",
     "threshold": "Shapiro-Wilk p < 0.001", "severity": "MODERATE",
     "cause": "PRS distribution is skewed or heavy-tailed in reference population.",
     "effect": "Normal-theory percentiles inaccurate in distribution tails."},

    # --- BENCHMARK FAILURES ---
    {"id": "FM-BMK-001", "component": "Benchmarking", "failure": "Circular validation",
     "threshold": "Benchmark reference = internal platform output", "severity": "CRITICAL",
     "cause": "Training data used as validation reference.",
     "effect": "Benchmarking provides no independent evidence of validity."},
    {"id": "FM-BMK-002", "component": "Benchmarking", "failure": "Single-sample pseudo-correlation",
     "threshold": "N < 30 for concordance analysis", "severity": "HIGH",
     "cause": "Correlation computed from single observation — statistically meaningless.",
     "effect": "Reported r values are distance heuristics, not correlations."},
    {"id": "FM-BMK-003", "component": "Benchmarking", "failure": "PGS Catalog score mismatch",
     "threshold": "Trait mapping confidence < 0.50", "severity": "MODERATE",
     "cause": "Platform trait categories don't map cleanly to PGS Catalog traits.",
     "effect": "Comparisons are between different phenotypes — concordance lost."},

    # --- REPRODUCIBILITY FAILURES ---
    {"id": "FM-REP-001", "component": "Reproducibility", "failure": "Unseeded randomness",
     "threshold": "PYTHONHASHSEED not set", "severity": "HIGH",
     "cause": "Python hash randomization affects dict ordering, set iteration.",
     "effect": "Non-deterministic outputs across runs — not bit-reproducible."},
    {"id": "FM-REP-002", "component": "Reproducibility", "failure": "PLINK version mismatch",
     "threshold": "PLINK 1.9 ≠ PLINK 2.0 output format differences", "severity": "MODERATE",
     "cause": "PLINK versions produce slightly different numeric outputs.",
     "effect": "Cross-version hashes don't match — environment-dependent results."},
]

@dataclass
class FailureModeMap:
    failures: List[Dict]; n_critical: int; n_high: int; n_moderate: int
    most_vulnerable_component: str; generated_date: str

class FailureModeMapper:
    """Maps every known failure mode with thresholds and severities."""

    def __init__(self, output_dir: str = "science"):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)

    def map(self, adversarial_json: Optional[str] = None,
            leakage_json: Optional[str] = None) -> FailureModeMap:
        logger.info("═══ Scientific Failure Mode Map ═══")

        # Cross-reference with adversarial results
        adversarial_findings = set()
        if adversarial_json and Path(adversarial_json).exists():
            try:
                with open(adversarial_json) as fh:
                    adv = json.load(fh)
                adversarial_findings = set(adv.get("critical_findings", []))
            except Exception: pass

        failures = FAILURE_MODES.copy()
        # Annotate with adversarial validation
        for f in failures:
            f["adversarial_validated"] = any(
                tag in f["id"] for tag in adversarial_findings)

        n_crit = sum(1 for f in failures if f["severity"] == "CRITICAL")
        n_high = sum(1 for f in failures if f["severity"] == "HIGH")
        n_mod = sum(1 for f in failures if f["severity"] == "MODERATE")

        # Count by component
        comps = {}
        for f in failures:
            comp = f["component"]
            comps[comp] = comps.get(comp, 0) + (3 if f["severity"] == "CRITICAL" else (2 if f["severity"] == "HIGH" else 1))
        most_vulnerable = max(comps, key=comps.get) if comps else "Unknown"

        map_result = FailureModeMap(
            failures=failures, n_critical=n_crit, n_high=n_high, n_moderate=n_mod,
            most_vulnerable_component=most_vulnerable,
            generated_date=datetime.now().strftime("%Y-%m-%d %H:%M UTC"))

        self._save(map_result)
        return map_result

    def _save(self, fm: FailureModeMap) -> None:
        with open(self.output_dir / "failure_mode_map.json", "w") as fh:
            json.dump({
                "n_failures": len(fm.failures),
                "n_critical": fm.n_critical, "n_high": fm.n_high,
                "n_moderate": fm.n_moderate,
                "most_vulnerable_component": fm.most_vulnerable_component,
                "generated_date": fm.generated_date,
                "failures": fm.failures,
            }, fh, indent=2)

        # Markdown
        lines = [
            "# Scientific Failure Mode Map",
            f"\n**Total failure modes:** {len(fm.failures)}",
            f"**Critical:** {fm.n_critical} | **High:** {fm.n_high} | **Moderate:** {fm.n_moderate}",
            f"**Most vulnerable:** {fm.most_vulnerable_component}",
            "\n## All Failure Modes",
            "\n| ID | Component | Failure | Threshold | Severity |",
            "|----|-----------|---------|-----------|----------|",
        ]
        for f in fm.failures:
            lines.append(f"| {f['id']} | {f['component']} | {f['failure']} | {f['threshold']} | {f['severity']} |")

        lines += [
            "\n## Critical Failures (Publication Blockers)",
        ]
        for f in fm.failures:
            if f["severity"] == "CRITICAL":
                lines.append(f"\n### {f['id']}: {f['failure']}")
                lines.append(f"\n**Cause:** {f['cause']}")
                lines.append(f"\n**Effect:** {f['effect']}")

        lines += ["\n---", "\n*Phase 10 Failure Mode Map — Peer Review Resilience*"]
        with open(self.output_dir / "failure_mode_map.md", "w") as fh:
            fh.write("\n".join(lines))

        logger.info(f"  ✅ Failure map: {len(fm.failures)} modes, {fm.n_critical} critical")

def main():
    import argparse
    parser = argparse.ArgumentParser(description="Phase 10: Scientific Failure Mode Map")
    parser.add_argument("--adversarial-json", default="science/adversarial_validation_report.json")
    parser.add_argument("--leakage-json", default="science/leakage_audit.json")
    parser.add_argument("--output-dir", "-o", default="science")
    parser.add_argument("--verbose", "-v", action="store_true")
    args = parser.parse_args()
    logging.basicConfig(level=logging.DEBUG if args.verbose else logging.INFO,
                        format="%(asctime)s [%(levelname)s] %(message)s")
    mapper = FailureModeMapper(args.output_dir)
    fm = mapper.map(args.adversarial_json, args.leakage_json)
    print(f"\n═══ Failure Mode Map ═══")
    print(f"  Total: {len(fm.failures)} | Critical: {fm.n_critical} | High: {fm.n_high} | Moderate: {fm.n_moderate}")
    print(f"  Most vulnerable: {fm.most_vulnerable_component}")
    return 0

if __name__ == "__main__":
    sys.exit(main())
