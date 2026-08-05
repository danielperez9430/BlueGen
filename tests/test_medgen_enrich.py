"""Unit tests for MedGen CUI persistence (IMPROVEMENT_PLAN.md 1.3).

find_disease_definition() already computed a medgen_cui for every matched
disease, but enrich_clinvar_output() discarded it after using it to look up
the definition text - so the report had no way to link out to the source
MedGen record. These tests guard the fix: the CUI must survive onto the
variant record.
"""

import json
import os
import sys

sys.path.insert(0, os.path.join(
    os.path.dirname(__file__),
    "..", "prs_research_pipeline", "scripts", "clinical"
))

import medgen_enrich  # noqa: E402


FAKE_DB = {
    "name_to_cui": {"hereditary hemochromatosis": "C0018995"},
    "cui_to_def": {"C0018995": "A disorder of iron metabolism causing iron overload."},
    "cui_to_name": {"C0018995": "Hereditary hemochromatosis"},
}


class TestFindDiseaseDefinitionReturnsCui:
    def test_exact_match_includes_cui(self):
        result = medgen_enrich.find_disease_definition("Hereditary hemochromatosis", FAKE_DB)
        assert result is not None
        assert result["medgen_cui"] == "C0018995"


class TestEnrichClinvarOutputPersistsCui:
    def test_matched_variant_gets_medgen_cui_field(self, tmp_path, monkeypatch):
        clinvar_json = tmp_path / "clinvar_pathogenic_variants.json"
        clinvar_json.write_text(json.dumps({
            "pathogenic_variants": [
                {"rsid": "rs1800562", "disease_name": "Hereditary hemochromatosis"},
            ]
        }))

        monkeypatch.setattr(medgen_enrich, "check_medgen_update", lambda ref_dir: {"status": "ok"})
        monkeypatch.setattr(medgen_enrich, "load_medgen_database", lambda ref_dir: FAKE_DB)

        result = medgen_enrich.enrich_clinvar_output(str(clinvar_json), ref_dir="unused")
        assert result["enriched"] == 1

        saved = json.loads(clinvar_json.read_text())
        variant = saved["pathogenic_variants"][0]
        assert variant["medgen_cui"] == "C0018995"
        assert "iron overload" in variant["disease_description"]

    def test_unmatched_variant_gets_empty_cui_not_missing_key(self, tmp_path, monkeypatch):
        clinvar_json = tmp_path / "clinvar_pathogenic_variants.json"
        clinvar_json.write_text(json.dumps({
            "pathogenic_variants": [
                {"rsid": "rs0000000", "disease_name": "not_provided"},
            ]
        }))

        monkeypatch.setattr(medgen_enrich, "check_medgen_update", lambda ref_dir: {"status": "ok"})
        monkeypatch.setattr(medgen_enrich, "load_medgen_database", lambda ref_dir: FAKE_DB)

        medgen_enrich.enrich_clinvar_output(str(clinvar_json), ref_dir="unused")

        saved = json.loads(clinvar_json.read_text())
        variant = saved["pathogenic_variants"][0]
        assert variant["medgen_cui"] == ""
