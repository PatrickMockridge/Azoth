"""The input codec: the declared kinds, SI canonicalisation and the named refusals."""

from __future__ import annotations

import pytest

from azoth.rgen import inputs, runner

ANTOINE = {
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


def test_a_unitful_input_is_stored_in_si_and_decoded_to_a_quantity() -> None:
    encoded = inputs.encode("eos.antoine_vapor_pressure", ANTOINE)
    assert encoded["Tc"] == {"value": 190.56, "unit": "K"}
    assert encoded["Pc"] == {"value": 4_599_000.0, "unit": "Pa"}
    assert encoded["A"] == 5.23243
    assert encoded["form"] == "pow10"

    decoded = inputs.decode("eos.antoine_vapor_pressure", encoded)
    result = runner.run("eos.antoine_vapor_pressure", decoded, "python")
    assert result["result"]["p_sat"]["magnitude"] == pytest.approx(929295.9523003073, rel=1e-12)


def test_an_optional_input_is_omitted_rather_than_sent_as_null() -> None:
    encoded = inputs.encode(
        "hydraulics.darcy_weisbach",
        {"f": 0.02, "L": 100.0, "D": 0.1, "rho": 998.0, "v": 1.5},
    )
    assert "mu" not in encoded
    result = runner.run(
        "hydraulics.darcy_weisbach",
        inputs.decode("hydraulics.darcy_weisbach", encoded),
        "python",
    )
    assert result["result"]["re"] is None


def test_a_supplied_optional_input_survives_the_round_trip() -> None:
    encoded = inputs.encode(
        "hydraulics.darcy_weisbach",
        {"f": 0.02, "L": 100.0, "D": 0.1, "rho": 998.0, "v": 1.5, "mu": 0.001002},
    )
    assert encoded["mu"] == {"value": 0.001002, "unit": "Pa*s"}
    result = runner.run(
        "hydraulics.darcy_weisbach",
        inputs.decode("hydraulics.darcy_weisbach", encoded),
        "python",
    )
    assert result["result"]["regime"] == "turbulent"


def test_a_fitting_list_is_a_list_of_names() -> None:
    encoded = inputs.encode("hydraulics.crane_k_factors", {"fittings": ["90_elbow"], "f_t": 0.019})
    assert encoded["fittings"] == ["90_elbow"]
    decoded = inputs.decode("hydraulics.crane_k_factors", encoded)
    result = runner.run("hydraulics.crane_k_factors", decoded, "python")
    assert result["result"]["k_total"] > 0


def test_a_record_list_input_is_refused_by_name() -> None:
    with pytest.raises(inputs.UnsupportedCalc, match="zones"):
        inputs.encode("thermal.zone_rating", {"zones": []})


def test_a_model_id_is_not_a_calculation() -> None:
    with pytest.raises(inputs.UnknownCalc):
        inputs.encode("eos.rachford_rice", {})


def test_an_undeclared_keyword_is_refused() -> None:
    with pytest.raises(inputs.InputError, match="not_a_field"):
        inputs.encode("eos.antoine_vapor_pressure", {**ANTOINE, "not_a_field": 1})


def test_a_missing_required_input_is_named() -> None:
    incomplete = dict(ANTOINE)
    incomplete.pop("T")
    with pytest.raises(inputs.InputError, match="T"):
        inputs.encode("eos.antoine_vapor_pressure", incomplete)


def test_an_enum_outside_its_vocabulary_is_refused() -> None:
    with pytest.raises(inputs.InputError, match="not_a_form"):
        inputs.encode("eos.antoine_vapor_pressure", {**ANTOINE, "form": "not_a_form"})
