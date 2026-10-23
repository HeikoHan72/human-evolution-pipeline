"""Small unit tests for the pure helper functions (no Spark, no databases, no raw data needed)."""
import importlib
import sys
from pathlib import Path

import pytest

SRC = Path(__file__).resolve().parents[1] / "src"
sys.path.insert(0, str(SRC))

cfg = importlib.import_module("00_config")
gf = pytest.importorskip("pysam") and importlib.import_module("genom_funktionen")


def test_ancestral_is_reference():
    assert gf.ahnenzustand_status("C", "T", "C") == "Referenz=ancestral"


def test_ancestral_is_alternative_and_case_insensitive():
    assert gf.ahnenzustand_status("C", "T", "t") == "Alt=ancestral"


def test_ancestral_matches_neither_allele():
    assert gf.ahnenzustand_status("C", "T", "G") == "weder_noch"


def test_ancestral_unknown_character():
    assert gf.ahnenzustand_status("C", "T", "N") == "kein_Ahnenzustand"


def test_genotype_text_never_turns_missing_into_homozygous_reference():
    assert gf.gt_text((0, 1)) == "0/1"
    assert gf.gt_text((None, 1)) == "./."
    assert gf.gt_text(None) == "./."


def test_first_value_of_info_tuple():
    assert gf._erster_wert((0.25, 0.5)) == 0.25
    assert gf._erster_wert(()) is None
    assert gf._erster_wert(None) is None


def test_config_regions_cover_their_genes():
    assert set(cfg.REGIONEN) == {"TLR", "STAT2", "OAS"}
    assert len(cfg.GENE) == 7
    for name, region in cfg.REGIONEN.items():
        assert region["start"] < region["ende"], name


def test_missing_secret_stops_with_clear_message(monkeypatch):
    monkeypatch.delenv("POSTGRES_PASSWORD", raising=False)
    with pytest.raises(SystemExit):
        cfg.env("POSTGRES_PASSWORD")
