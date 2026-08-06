"""Unit tests for bluegen.ancestry (Stage G, PCA adjustment) - IMPROVEMENT_PLAN.md 2.1.

Before this file, pca_adjust_v2.py had zero direct unit tests - only
indirect coverage via the full pipeline. These exercise the real formulas:
adjusted = raw - sum(beta_k * PC_k), the shrinkage-prior fallback used when
no reference betas exist for a trait, and the OLS regression that estimates
reference betas from a cohort.
"""

import csv
import json
import math
import sys
import tempfile
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "prs_research_pipeline"))

from bluegen.ancestry import PCAAdjustmentV2


def _write_prs_csv(tmpdir, rows):
    path = Path(tmpdir) / "prs_raw.csv"
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["individual_id", "trait", "prs_raw"])
        for trait, prs_raw in rows:
            w.writerow(["SAMPLE_001", trait, prs_raw])
    return str(path)


def _write_sample_pcs(tmpdir, pcs):
    path = Path(tmpdir) / "target_pcs.eigenvec"
    header = "FID\tIID\t" + "\t".join(f"PC{i+1}" for i in range(len(pcs)))
    values = "SAMPLE_001\tSAMPLE_001\t" + "\t".join(str(v) for v in pcs)
    path.write_text(header + "\n" + values + "\n")
    return str(path)


def _write_ref_betas(tmpdir, betas):
    path = Path(tmpdir) / "reference_pc_betas.json"
    path.write_text(json.dumps(betas))
    return str(path)


class TestAdjustSubtractsKnownBetas:
    def test_adjusted_equals_raw_minus_pc_contribution(self):
        pcs = [1.0, 2.0]
        betas = {"Test trait": {"PC1": 0.5, "PC2": -0.25}}
        with tempfile.TemporaryDirectory() as tmpdir:
            prs_data = _write_prs_csv(tmpdir, [("Test trait", 10.0)])
            sample_pcs = _write_sample_pcs(tmpdir, pcs)
            ref_betas = _write_ref_betas(tmpdir, betas)

            adjuster = PCAAdjustmentV2(n_pcs=2)
            report = adjuster.adjust(
                prs_data=prs_data, sample_pcs=sample_pcs, output_dir=tmpdir,
                ref_beta_path=ref_betas, now="2020-01-01 00:00 UTC",
            )

        # expected: 10.0 - (0.5*1.0 + -0.25*2.0) = 10.0 - 0.0 = 10.0
        result = report.results[0]
        assert result.raw_prs == 10.0
        assert result.adjusted_prs == 10.0
        assert result.delta == 0.0

    def test_nonzero_pc_contribution_changes_adjusted_value(self):
        pcs = [2.0]
        betas = {"Test trait": {"PC1": 1.5}}
        with tempfile.TemporaryDirectory() as tmpdir:
            prs_data = _write_prs_csv(tmpdir, [("Test trait", 10.0)])
            sample_pcs = _write_sample_pcs(tmpdir, pcs)
            ref_betas = _write_ref_betas(tmpdir, betas)

            adjuster = PCAAdjustmentV2(n_pcs=1)
            report = adjuster.adjust(
                prs_data=prs_data, sample_pcs=sample_pcs, output_dir=tmpdir,
                ref_beta_path=ref_betas, now="2020-01-01 00:00 UTC",
            )

        # expected: 10.0 - (1.5*2.0) = 10.0 - 3.0 = 7.0
        result = report.results[0]
        assert result.adjusted_prs == 7.0
        assert result.delta == 3.0
        assert result.is_significant is True


class TestInjectableGeneratedDate:
    def test_now_parameter_is_used_verbatim(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            prs_data = _write_prs_csv(tmpdir, [("Test trait", 1.0)])
            sample_pcs = _write_sample_pcs(tmpdir, [0.0])

            adjuster = PCAAdjustmentV2(n_pcs=1)
            report = adjuster.adjust(
                prs_data=prs_data, sample_pcs=sample_pcs, output_dir=tmpdir,
                now="2020-01-01 00:00 UTC",
            )
        assert report.generated_date == "2020-01-01 00:00 UTC"


class TestEstimateBetasFromPcsShrinkagePrior:
    def test_matches_exact_shrinkage_formula(self):
        adjuster = PCAAdjustmentV2(n_pcs=3)
        prs_value = 4.0
        pcs = np.array([2.0, -1.0, 0.0])

        betas = adjuster._estimate_betas_from_pcs(prs_value, pcs, "Test trait")

        tau = 0.02 * abs(prs_value) / max(np.sqrt(3), 1.0)
        assert math.isclose(betas["PC1"], tau * (1.0 / math.sqrt(1)) * 1.0, rel_tol=1e-9)
        assert math.isclose(betas["PC2"], tau * (1.0 / math.sqrt(2)) * -1.0, rel_tol=1e-9)
        assert betas["PC3"] == 0.0  # sign(0) special-cased to 0.0, not NaN

    def test_shrinkage_scales_with_prs_magnitude(self):
        adjuster = PCAAdjustmentV2(n_pcs=1)
        small = adjuster._estimate_betas_from_pcs(1.0, np.array([1.0]), "t")
        large = adjuster._estimate_betas_from_pcs(10.0, np.array([1.0]), "t")
        assert abs(large["PC1"]) > abs(small["PC1"])


class TestComputeRefBetasRecoversKnownSlope:
    def test_ols_recovers_true_slope_from_noiseless_cohort(self):
        # PRS = 3.0 * PC1 + 0.0 (no noise, no other PCs) for 20 synthetic samples.
        true_beta = 3.0
        n_samples = 20
        pc1_values = np.linspace(-2, 2, n_samples)

        with tempfile.TemporaryDirectory() as tmpdir:
            eigenvec_path = Path(tmpdir) / "ref.eigenvec"
            with open(eigenvec_path, "w") as fh:
                for i, pc1 in enumerate(pc1_values):
                    fh.write(f"FAM{i}\tSAMPLE{i}\t{pc1}\t0.0\n")  # PC1, PC2(=0)

            prs_path = Path(tmpdir) / "ref_prs.csv"
            with open(prs_path, "w", newline="") as fh:
                w = csv.writer(fh)
                w.writerow(["individual_id", "trait", "prs_raw"])
                for i, pc1 in enumerate(pc1_values):
                    w.writerow([f"SAMPLE{i}", "Test trait", true_beta * pc1])

            adjuster = PCAAdjustmentV2(n_pcs=2)
            ref_betas = adjuster.compute_ref_betas(
                ref_prs_data=str(prs_path), ref_pcs_path=str(eigenvec_path),
                population_panel="unused", output_dir=tmpdir,
            )

        assert math.isclose(ref_betas["Test trait"]["PC1"], true_beta, abs_tol=1e-6)
        assert math.isclose(ref_betas["Test trait"]["PC2"], 0.0, abs_tol=1e-6)
        assert ref_betas["Test trait"]["r_squared"] >= 0.999  # noiseless -> perfect fit
