"""Unit tests for bluegen.scoring (Stage F, PLINK PRS computation) - IMPROVEMENT_PLAN.md 2.1.

Before this file, prs_plink_score.py had zero direct unit tests - only
indirect coverage (tests/test_snp_positions.py, test_allele_strand_consistency.py
test the upstream SNP panel data, not this module's matching logic itself).

build_score_rows() is the pure part: no I/O, no PLINK, fully testable with a
fake bim_ids set. write_score_file()/parse_plink_profile() are tested against
real temp files but never invoke PLINK itself.
"""

import sys
import tempfile
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "prs_research_pipeline"))

from bluegen.scoring import build_score_rows, write_score_file, parse_plink_profile


def _panel(rows):
    """rows: list of dicts with chrom/pos/effect_allele/weight (all str, matching
    the real pd.read_csv(snp_db, dtype=str) upstream)."""
    return pd.DataFrame(rows, dtype=str)


class TestBuildScoreRowsMatching:
    def test_matches_by_chrom_colon_pos(self):
        panel = _panel([
            {"chrom": "1", "pos": "1000", "effect_allele": "A", "weight": "0.5"},
        ])
        rows = build_score_rows(panel, bim_ids={"1:1000"})
        assert rows == [("1:1000", "A", 0.5)]

    def test_position_not_in_bim_is_silently_dropped(self):
        panel = _panel([
            {"chrom": "1", "pos": "1000", "effect_allele": "A", "weight": "0.5"},
            {"chrom": "2", "pos": "2000", "effect_allele": "T", "weight": "0.3"},
        ])
        rows = build_score_rows(panel, bim_ids={"1:1000"})  # "2:2000" not genotyped
        assert rows == [("1:1000", "A", 0.5)]

    def test_chr_prefix_is_stripped_before_matching(self):
        panel = _panel([
            {"chrom": "chr1", "pos": "1000", "effect_allele": "A", "weight": "0.5"},
        ])
        rows = build_score_rows(panel, bim_ids={"1:1000"})
        assert rows == [("1:1000", "A", 0.5)]

    def test_invalid_weight_is_silently_skipped(self):
        panel = _panel([
            {"chrom": "1", "pos": "1000", "effect_allele": "A", "weight": "not_a_number"},
        ])
        rows = build_score_rows(panel, bim_ids={"1:1000"})
        assert rows == []

    def test_missing_effect_allele_is_skipped(self):
        panel = _panel([
            {"chrom": "1", "pos": "1000", "effect_allele": "", "weight": "0.5"},
        ])
        rows = build_score_rows(panel, bim_ids={"1:1000"})
        assert rows == []

    def test_missing_weight_defaults_to_1_0(self):
        panel = pd.DataFrame([{"chrom": "1", "pos": "1000", "effect_allele": "A"}])
        rows = build_score_rows(panel, bim_ids={"1:1000"})
        assert rows == [("1:1000", "A", 1.0)]

    def test_multiple_matching_rows_preserve_order(self):
        panel = _panel([
            {"chrom": "1", "pos": "1000", "effect_allele": "A", "weight": "0.1"},
            {"chrom": "1", "pos": "2000", "effect_allele": "G", "weight": "0.2"},
        ])
        rows = build_score_rows(panel, bim_ids={"1:1000", "1:2000"})
        assert rows == [("1:1000", "A", 0.1), ("1:2000", "G", 0.2)]


class TestWriteScoreFile:
    def test_tab_separated_format(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "test.score"
            write_score_file(path, [("1:1000", "A", 0.5), ("2:2000", "T", 0.3)])
            content = path.read_text()
        assert content == "1:1000\tA\t0.5\n2:2000\tT\t0.3\n"


class TestParsePlinkProfile:
    def test_reads_score_column(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "test.profile"
            path.write_text(" FID   IID   PHENO    CNT   CNT2   SCORE\n"
                             " FAM1  SAMPLE_001  -9    10    8    1.2345\n")
            results = parse_plink_profile(path, "Test trait")
        assert len(results) == 1
        assert results[0] == {
            "individual_id": "SAMPLE_001", "trait": "Test trait",
            "prs_raw": 1.2345, "n_snps": 10, "n_snps_used": 8,
        }

    def test_falls_back_to_scoresum_when_score_absent(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "test.profile"
            path.write_text(" FID   IID   PHENO    CNT   CNT2   SCORESUM\n"
                             " FAM1  SAMPLE_001  -9    10    8    9.87\n")
            results = parse_plink_profile(path, "Test trait")
        assert results[0]["prs_raw"] == 9.87
