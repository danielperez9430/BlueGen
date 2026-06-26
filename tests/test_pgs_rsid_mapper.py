"""Unit tests for rsID → chr:pos mapper."""

import sys
import os
import json
import tempfile
from pathlib import Path
import pytest

# Add scripts directories to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..",
    "prs_research_pipeline", "scripts", "benchmarking"))

from pgs_rsid_mapper import RsidMapper


class TestRsidMapper:
    def test_lookup_from_in_memory(self):
        mapper = RsidMapper()
        mapper._mapping["rs123"] = "10:114758349"
        mapper._mapping["rs456"] = "16:89985861"
        assert mapper.lookup("rs123") == "10:114758349"
        assert mapper.lookup("rs456") == "16:89985861"

    def test_lookup_not_found(self):
        mapper = RsidMapper()
        assert mapper.lookup("rs999999") is None

    def test_lookup_without_rs_prefix(self):
        mapper = RsidMapper()
        mapper._mapping["rs123"] = "10:114758349"
        assert mapper.lookup("123") == "10:114758349"

    def test_batch_lookup(self):
        mapper = RsidMapper()
        mapper._mapping["rs1"] = "1:100"
        mapper._mapping["rs2"] = "2:200"
        result = mapper.batch_lookup(["rs1", "rs2", "rs3"])
        assert result["rs1"] == "1:100"
        assert result["rs2"] == "2:200"
        assert result["rs3"] is None

    def test_cache_save_load_roundtrip(self):
        mapper = RsidMapper()
        mapper._mapping["rsA"] = "1:100"
        mapper._mapping["rsB"] = "2:200"

        with tempfile.NamedTemporaryFile(mode="w", suffix=".json", delete=False) as f:
            cache_path = f.name
        try:
            mapper.save_cache(cache_path)
            # Load into a new mapper
            mapper2 = RsidMapper()
            mapper2.load_from_cache(cache_path)
            assert mapper2.lookup("rsA") == "1:100"
            assert mapper2.lookup("rsB") == "2:200"
        finally:
            os.unlink(cache_path)

    def test_load_from_cache_missing_file(self):
        mapper = RsidMapper()
        n = mapper.load_from_cache("/nonexistent/path/cache.json")
        assert n == 0

    def test_n_entries(self):
        mapper = RsidMapper()
        assert mapper.n_entries == 0
        mapper._mapping["rsX"] = "1:1"
        assert mapper.n_entries == 1

    def test_stats(self):
        mapper = RsidMapper()
        assert mapper.stats["miss"] == 0
        mapper.lookup("rs_nonexistent")
        assert mapper.stats["miss"] == 1
