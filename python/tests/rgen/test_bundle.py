"""The whole bundle: the demo case, the documents, the refusals and the round trip."""

from __future__ import annotations

from pathlib import Path

import pytest

from azoth.rgen import Bundle, BundleRefused, build, canonical, inputs, payloads
from azoth.rgen.__main__ import main

ANTOINE: dict[str, float | str] = {
    "A": 5.23243,
    "B": 891.0098,
    "C": 332.0975,
    "D": 0.0,
    "E": 0.0,
    "form": "pow10",
    "Tc": 190.56,
    "Pc": 4_599_000.0,
    "T": 150.0,
}


@pytest.fixture(scope="module")
def bundle() -> Bundle:
    return build(
        "eos.antoine_vapor_pressure",
        ANTOINE,
        project="azoth-demo",
        jurisdiction="US-TX",
    )


def test_the_demo_has_five_validated_cases(bundle: Bundle) -> None:
    assert bundle.manifest["quantity"] == {"unit": "validated-cases", "value": 5}
    assert len(bundle.documents["cases"]["declared"]["cases"]) == 5
    assert all(report["passed"] for report in bundle.documents["cases"]["replay"])


def test_the_demo_is_valid_and_scores_well(bundle: Bundle) -> None:
    assert bundle.manifest["verdict"] == "valid"
    quality = bundle.manifest["quality"]
    assert set(quality) == {"completeness", "correctness", "traceability", "fitness"}
    assert quality["completeness"] == 5
    assert quality["correctness"] == 5
    assert quality["fitness"] == 5


def test_the_statement_carries_the_result_provenance_block(bundle: Bundle) -> None:
    provenance = bundle.documents["statement"]["result"]["provenance"]
    assert provenance["calc_id"] == "eos.antoine_vapor_pressure"
    assert len(provenance["spec_sha256"]) == 64
    assert len(provenance["python_sha256"]) == 64
    assert len(provenance["rust_sha256"]) == 64


def test_seven_documents_each_hash_to_their_record(bundle: Bundle) -> None:
    assert set(bundle.documents) == {
        "definition",
        "build",
        "statement",
        "cases",
        "agreement",
        "quality",
        "validation-evidence",
    }
    for name, document in bundle.documents.items():
        assert canonical.hash_value(document) == bundle.hashes[name]


def test_the_payloads_are_integer_only_and_the_claim_counts_five(bundle: Bundle) -> None:
    payloads.integers_only(bundle.claim)
    payloads.integers_only(bundle.validation)
    assert bundle.claim["quantity"]["value"] == 5
    assert set(bundle.validation["quality"]) == {
        "completeness",
        "correctness",
        "traceability",
        "fitness",
    }


def test_the_write_round_trip_verifies(tmp_path: Path, bundle: Bundle) -> None:
    bundle.write(tmp_path)
    for name in bundle.documents:
        assert (tmp_path / "documents" / f"{name}.json").is_file()
    assert (tmp_path / "payloads" / "claim.json").is_file()
    assert (tmp_path / "manifest.json").is_file()
    assert main(["verify", "--dir", str(tmp_path)]) == 0


def test_the_bundle_is_byte_reproducible(bundle: Bundle) -> None:
    second = build(
        "eos.antoine_vapor_pressure",
        ANTOINE,
        project="azoth-demo",
        jurisdiction="US-TX",
    )
    for name, document in bundle.documents.items():
        assert canonical.dumps(document) == canonical.dumps(second.documents[name])


def test_a_validation_needs_a_different_implementation() -> None:
    with pytest.raises(BundleRefused, match="different implementation"):
        build(
            "eos.antoine_vapor_pressure",
            ANTOINE,
            project="azoth-demo",
            jurisdiction="US-TX",
            claim_backend="python",
            validation_backend="python",
        )


def test_a_record_list_input_is_refused_by_name() -> None:
    with pytest.raises(inputs.UnsupportedCalc, match="zones"):
        build(
            "thermal.zone_rating",
            {"zones": []},
            project="azoth-demo",
            jurisdiction="US-TX",
        )


def test_an_exploratory_run_is_recorded_and_scores_at_most_three() -> None:
    exploratory = build(
        "eos.antoine_vapor_pressure",
        {**ANTOINE, "T": 151.0},
        project="azoth-demo",
        jurisdiction="US-TX",
    )
    assert exploratory.documents["agreement"]["run_matched"] is None
    assert exploratory.manifest["verdict"] == "partial"
    assert exploratory.manifest["quality"]["correctness"] <= 3
    assert exploratory.manifest["quantity"]["value"] == 5
