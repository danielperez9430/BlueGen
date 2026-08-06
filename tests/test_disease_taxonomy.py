"""Unit tests for the ClinVar body-system classifier (IMPROVEMENT_PLAN.md 1.3)."""

import os
import sys

sys.path.insert(0, os.path.join(
    os.path.dirname(__file__),
    "..", "prs_research_pipeline", "scripts", "clinical"
))

import disease_taxonomy as dt  # noqa: E402


class TestClassifyBySymbol:
    def test_known_gene_list_classifies_correctly(self):
        assert dt.classify_body_system(["BRCA1"]) == "oncology"
        assert dt.classify_body_system(["CFTR"]) == "respiratory"
        assert dt.classify_body_system(["LDLR"]) == "cardiovascular"

    def test_raw_gene_info_string_is_parsed(self):
        # clinvar rows sometimes store gene_info as "SYMBOL:entrez_id|OTHER:id"
        assert dt.classify_body_system("CFTR:1080|CFTR-AS1:111082987") == "respiratory"

    def test_gene_list_checks_all_entries_not_just_first(self):
        assert dt.classify_body_system(["LOC12345", "BRCA2"]) == "oncology"


class TestClassifyByDiseaseNameFallback:
    def test_unknown_gene_falls_back_to_disease_keyword(self):
        assert dt.classify_body_system(["ZZZNOTAGENE"], "Renal cell carcinoma") == "oncology"
        assert dt.classify_body_system(["ZZZNOTAGENE"], "Cystic_fibrosis") == "respiratory"

    def test_underscore_disease_names_are_normalized(self):
        assert dt.classify_body_system([], "Holt-Oram_syndrome") == "cardiovascular"


class TestClassifyFallbackToOther:
    def test_unknown_gene_and_disease_returns_other(self):
        assert dt.classify_body_system(["ZZZNOTAGENE"], "not_provided") == "other"

    def test_empty_input_returns_other(self):
        assert dt.classify_body_system([], "") == "other"


class TestSystemLabel:
    def test_known_system_has_bilingual_label(self):
        en = dt.system_label("oncology", "en")
        es = dt.system_label("oncology", "es")
        assert "Oncology" in en
        assert "Oncología" in es

    def test_unknown_system_key_falls_back_to_other_label(self):
        assert dt.system_label("not_a_real_system", "en") == dt.system_label("other", "en")
