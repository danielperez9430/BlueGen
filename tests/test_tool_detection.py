"""Tests for utils.tool_detection module."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "prs_research_pipeline" / "scripts"))

from utils.tool_detection import (
    _find_plink_raw,
    find_bcftools,
    find_tabix,
    detect_all,
)


def test_find_plink_raw_returns_tuple():
    result = _find_plink_raw()
    assert isinstance(result, tuple)
    assert len(result) == 2


def test_find_plink_raw_none_or_strings():
    path, ver = _find_plink_raw()
    if path is not None:
        assert isinstance(path, str)
        assert len(path) > 0
    if ver is not None:
        assert isinstance(ver, str)


def test_find_plink_raw_handles_missing():
    path, ver = _find_plink_raw()
    # Should return (None, None) when PLINK not found, or (str, str) when found
    assert path is None or isinstance(path, str)
    assert ver is None or isinstance(ver, str)


def test_find_bcftools_returns_string_or_none():
    result = find_bcftools()
    assert result is None or isinstance(result, str)


def test_find_tabix_returns_string_or_none():
    result = find_tabix()
    assert result is None or isinstance(result, str)


def test_detect_all_returns_dict():
    result = detect_all()
    assert isinstance(result, dict)
    assert "plink" in result
    assert "bcftools" in result
    assert "tabix" in result


def test_detect_all_plink_has_fields():
    result = detect_all()
    plink = result["plink"]
    assert "path" in plink
    assert "version" in plink
