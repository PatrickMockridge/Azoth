"""Spec-driven tests for ``eos.ideal_gas_cp``."""

from __future__ import annotations

from typing import Any

import pytest

import _helpers as h
from azoth import ureg, use_backend
from azoth._registry_gen import spec as spec_for
from azoth.core.errors import OutOfRangeError
from azoth.core.result import IdealGasCpResult
from azoth.eos import ideal_gas_cp
from azoth.eos.reference.pr_molar_volume import MOLAR_GAS_CONSTANT

CALC_ID = "eos.ideal_gas_cp"
Q = ureg.Quantity

SPEC = spec_for(CALC_ID)


def call(case: dict[str, Any]) -> IdealGasCpResult:
    inputs = case["inputs"]
    return ideal_gas_cp(inputs["a"], inputs["b"], inputs["c"], inputs["d"], Q(inputs["T"], "K"))


def _cases() -> list[dict[str, Any]]:
    """The worked example plus every other test that pins values.

    The `tests` list opens with a *pointer* to the worked example - an entry with an
    id, a type and nothing else - because that is the shape every calc's spec uses.
    Reading it as a case gives a case with no inputs, so it is dropped rather than
    skipped, and the worked example itself is taken from the block it points at.
    """
    block = SPEC["worked_example"]
    others = [t for t in SPEC["tests"] if "inputs" in t]
    return [{"id": "worked_example", **block, "type": "worked_example"}, *others]


ACTIVE = [
    c
    for c in _cases()
    if c.get("status") == "active" and c.get("type") in ("worked_example", "reference")
]


@pytest.mark.parametrize("case", ACTIVE, ids=lambda c: c["id"])
def test_spec_case(case: dict[str, Any]) -> None:
    result = call(case)
    h.assert_close(
        result.cp_over_r,
        case["expected"]["cp_over_r"],
        case["tolerance"],
        f"{case['id']} (cp_over_r)",
    )
    h.assert_close(
        result.cp.to("J/(mol*K)").magnitude,
        case["expected"]["cp"],
        case["tolerance"],
        f"{case['id']} (cp)",
    )
    h.assert_consistent(result, case["id"])


def test_every_case_ran() -> None:
    assert len(ACTIVE) >= 3, f"expected several cases, found {len(ACTIVE)}"


def test_the_polynomial_is_linear_in_its_coefficients() -> None:
    """Linear in ``(a, b, c, d)``, which holds for every coefficient set.

    Catches an implementation which multiplies ``d`` by ``theta`` twice, or drops a
    term - both of which pass at the values one worked example happens to use.
    """
    base = ideal_gas_cp(4.0, 1.0, -0.5, 0.1, Q(500.0, "K"))
    doubled = ideal_gas_cp(8.0, 2.0, -1.0, 0.2, Q(500.0, "K"))
    negated = ideal_gas_cp(-4.0, -1.0, 0.5, -0.1, Q(500.0, "K"))
    h.assert_close(doubled.cp_over_r, 2.0 * base.cp_over_r, 1e-15, "doubled")
    h.assert_close(negated.cp_over_r, -base.cp_over_r, 1e-15, "negated")


def test_the_reference_temperature_is_a_thousand_kelvin() -> None:
    """At 1000 K the polynomial collapses to ``a + b + c + d``.

    The one temperature at which the four coefficients can be checked against the
    answer by addition alone, and the check that catches an implementation using
    ``T/100`` or ``T`` - all of which agree elsewhere after rescaling coefficients.
    """
    result = ideal_gas_cp(4.0, -20.0, 10.0, 1.0, Q(1000.0, "K"))
    h.assert_close(result.cp_over_r, -5.0, 1e-15, "a + b + c + d")
    h.assert_close(result.cp.to("J/(mol*K)").magnitude, -5.0 * MOLAR_GAS_CONSTANT, 1e-12, "cp")


def test_a_non_positive_heat_capacity_warns() -> None:
    """The failure this catches is the one that matters.

    A negative heat capacity fed into ``integral Cp dT`` produces an enthalpy wrong by
    an amount nobody can see, so it must not pass in silence.
    """
    result = ideal_gas_cp(4.0, -20.0, 10.0, 1.0, Q(1000.0, "K"))
    assert result.cp.to("J/(mol*K)").magnitude < 0.0
    assert result.has_warning(__import__("azoth").WarningCode.OUT_OF_VALID_RANGE), (
        f"a negative heat capacity should say so: {result.warnings}"
    )
    assert ideal_gas_cp(4.0, 1.0, -0.5, 0.1, Q(500.0, "K")).is_clean


@pytest.mark.parametrize("t", [0.0, -1.0])
def test_a_non_positive_temperature_is_refused(t: float) -> None:
    with pytest.raises(OutOfRangeError) as excinfo:
        ideal_gas_cp(4.0, 1.0, -0.5, 0.1, Q(t, "K"))
    assert excinfo.value.field() == "T"


@pytest.mark.parametrize(
    ("b", "c", "d"),
    [(1.5, -0.5, 0.1), (1.0, -0.9, 0.1), (1.0, -0.5, 0.5)],
    ids=["b", "c", "d"],
)
def test_every_coefficient_changes_the_answer(b: float, c: float, d: float) -> None:
    """A dropped term is a small correction at most coefficient sets, so it needs a test.

    Each of the four is perturbed on its own; an implementation that silently ignored
    one would be indistinguishable from a correct one at the values a single worked
    example uses.
    """
    base = ideal_gas_cp(4.0, 1.0, -0.5, 0.1, Q(700.0, "K"))
    changed = ideal_gas_cp(4.0, b, c, d, Q(700.0, "K"))
    assert abs(changed.cp_over_r - base.cp_over_r) > 1e-6


def test_the_constant_coefficient_changes_the_answer() -> None:
    """``a`` gets its own test because the parametrisation above holds it fixed."""
    base = ideal_gas_cp(4.0, 1.0, -0.5, 0.1, Q(700.0, "K"))
    changed = ideal_gas_cp(4.5, 1.0, -0.5, 0.1, Q(700.0, "K"))
    assert abs(changed.cp_over_r - base.cp_over_r) > 1e-6


def test_the_two_backends_agree() -> None:
    for case in ACTIVE:
        with use_backend("python"):
            py = call(case)
        with use_backend("rust"):
            rs = call(case)
        h.assert_close(py.cp_over_r, rs.cp_over_r, 1e-15, f"{case['id']} (cp_over_r)")
        h.assert_close(
            py.cp.to("J/(mol*K)").magnitude,
            rs.cp.to("J/(mol*K)").magnitude,
            1e-15,
            f"{case['id']} (cp)",
        )
        assert [w.code for w in py.warnings] == [w.code for w in rs.warnings]


def test_the_unit_round_trips_through_a_non_si_unit() -> None:
    """`T` crosses in kelvin, and a caller may supply it in another absolute unit.

    Unlike the rest of this namespace, this calc is dimensional, so the round trip is
    real rather than vacuous. A degree-Celsius input must give the same answer as the
    equivalent kelvin one, because `T` here is an absolute temperature and not an
    interval - which is the distinction `azoth.core.units` refuses to guess at.
    """
    kelvin = ideal_gas_cp(4.0, 1.0, -0.5, 0.1, Q(500.0, "K"))
    celsius = ideal_gas_cp(4.0, 1.0, -0.5, 0.1, Q(226.85, "degC"))
    h.assert_close(kelvin.cp_over_r, celsius.cp_over_r, 1e-12, "K against degC")
