"""Spec-driven tests for ``eos.ideal_gas_cp``.

The calc is a port of NeqSim's ``Component.getCp0``, and it is dimensional: each
coefficient carries a power of temperature, which is what lets the ``CPA``-``CPE``
columns of NeqSim's ``COMP.csv`` be read as they are stored. Every case below comes
from that file rather than from a table somebody typed.
"""

from __future__ import annotations

from typing import Any

import pytest

import _helpers as h
from azoth import ureg, use_backend
from azoth._registry_gen import spec as spec_for
from azoth.core.errors import OutOfRangeError
from azoth.core.result import IdealGasCpResult
from azoth.eos import ideal_gas_cp

CALC_ID = "eos.ideal_gas_cp"
Q = ureg.Quantity

SPEC = spec_for(CALC_ID)

#: The unit of each coefficient, taken from the spec rather than restated, so a spec
#: that changes one changes what these tests supply.
UNITS = {name: declaration["unit"] for name, declaration in SPEC["inputs"].items()}


def q(value: float, name: str) -> Any:
    """``value`` in the unit the spec declares for input ``name``."""
    return Q(value, UNITS[name])


def call(case: dict[str, Any]) -> IdealGasCpResult:
    inputs = case["inputs"]
    return ideal_gas_cp(
        q(inputs["cp_a"], "cp_a"),
        q(inputs["cp_b"], "cp_b"),
        q(inputs["cp_c"], "cp_c"),
        q(inputs["cp_d"], "cp_d"),
        q(inputs["cp_e"], "cp_e"),
        q(inputs["T"], "T"),
    )


def cp_of(cp_a: float, cp_b: float, cp_c: float, cp_d: float, cp_e: float, t: float) -> float:
    """The result of one call, as a magnitude in J/(mol*K)."""
    return (
        ideal_gas_cp(
            q(cp_a, "cp_a"),
            q(cp_b, "cp_b"),
            q(cp_c, "cp_c"),
            q(cp_d, "cp_d"),
            q(cp_e, "cp_e"),
            q(t, "T"),
        )
        .cp.to("J/(mol*K)")
        .magnitude
    )


def _cases() -> list[dict[str, Any]]:
    """The worked example plus every other test that pins values.

    The `tests` list opens with a *pointer* to the worked example - an entry with an
    id, a type and nothing else - because that is the shape every calc's spec uses.
    Reading it as a case gives a case with no inputs, so it is dropped rather than
    skipped, and the worked example itself is taken from the block it points at.
    """
    block = SPEC["worked_example"]
    others = [t for t in SPEC["tests"] if "inputs" in t]
    return [
        {"id": "worked_example", "status": "active", **block, "type": "worked_example"},
        *others,
    ]


ACTIVE = [
    c
    for c in _cases()
    if c.get("status") == "active" and c.get("type") in ("worked_example", "reference")
]

#: Methane's coefficients as NeqSim ships them, and the pair most of these tests use.
METHANE = (37.978352, -0.07461815, 0.000301881, -2.83e-07, 9.070574e-11)


@pytest.mark.parametrize("case", ACTIVE, ids=lambda c: c["id"])
def test_spec_case(case: dict[str, Any]) -> None:
    result = call(case)
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
    """Linear in all five coefficients, which holds for every coefficient set.

    Catches an implementation which multiplies a coefficient by the temperature twice,
    or drops a term - both of which pass at the values one worked example happens to
    use.
    """
    base = cp_of(METHANE[0], METHANE[1], METHANE[2], METHANE[3], METHANE[4], 500.0)
    doubled = cp_of(
        2.0 * METHANE[0],
        2.0 * METHANE[1],
        2.0 * METHANE[2],
        2.0 * METHANE[3],
        2.0 * METHANE[4],
        500.0,
    )
    negated = cp_of(-METHANE[0], -METHANE[1], -METHANE[2], -METHANE[3], -METHANE[4], 500.0)
    h.assert_close(doubled, 2.0 * base, 1e-12, "doubled")
    h.assert_close(negated, -base, 1e-12, "negated")


def test_at_one_kelvin_the_polynomial_collapses_to_the_sum_of_its_coefficients() -> None:
    """At 1 K every power of the temperature is one.

    The one temperature at which the five coefficients can be checked against the
    answer by addition alone. An implementation using a reduced temperature - a
    different scale, or a divisor - agrees elsewhere after rescaling its
    coefficients, and fails here.
    """
    coefficients = (4.0, -20.0, 10.0, 1.0, 0.5)
    h.assert_close(cp_of(*coefficients, 1.0), sum(coefficients), 1e-12, "the sum")


@pytest.mark.parametrize(
    ("index", "value"),
    [(0, 4.5), (1, 1.5), (2, -0.9), (3, 0.5), (4, 1e-11)],
    ids=["cp_a", "cp_b", "cp_c", "cp_d", "cp_e"],
)
def test_every_coefficient_changes_the_answer(index: int, value: float) -> None:
    """A dropped term is a small correction at most coefficient sets, so it needs a test.

    Each of the five is perturbed on its own, because an implementation that silently
    ignored one would be indistinguishable from a correct one at the values a single
    worked example uses - and `cp_e` is the one a four-term implementation drops.
    """
    perturbed = list(METHANE)
    perturbed[index] = value
    assert (
        abs(
            cp_of(perturbed[0], perturbed[1], perturbed[2], perturbed[3], perturbed[4], 700.0)
            - cp_of(*METHANE, 700.0)
        )
        > 1e-9
    )


def test_a_non_positive_heat_capacity_warns() -> None:
    """The failure this catches is the one that matters.

    A negative heat capacity fed into ``integral Cp dT`` gives an enthalpy wrong by an
    amount nobody can see, so it must not pass in silence.
    """
    result = ideal_gas_cp(
        q(4.0, "cp_a"),
        q(-20.0, "cp_b"),
        q(10.0, "cp_c"),
        q(1.0, "cp_d"),
        q(0.0, "cp_e"),
        q(1.0, "T"),
    )
    assert result.cp.to("J/(mol*K)").magnitude < 0.0
    assert result.has_warning(__import__("azoth").WarningCode.OUT_OF_VALID_RANGE), (
        f"a negative heat capacity should say so: {result.warnings}"
    )
    assert ideal_gas_cp(
        q(METHANE[0], "cp_a"),
        q(METHANE[1], "cp_b"),
        q(METHANE[2], "cp_c"),
        q(METHANE[3], "cp_d"),
        q(METHANE[4], "cp_e"),
        q(300.0, "T"),
    ).is_clean


@pytest.mark.parametrize("t", [0.0, -1.0])
def test_a_non_positive_temperature_is_refused(t: float) -> None:
    with pytest.raises(OutOfRangeError) as excinfo:
        ideal_gas_cp(
            q(METHANE[0], "cp_a"),
            q(METHANE[1], "cp_b"),
            q(METHANE[2], "cp_c"),
            q(METHANE[3], "cp_d"),
            q(METHANE[4], "cp_e"),
            q(t, "T"),
        )
    assert excinfo.value.field() == "T"


def test_the_two_backends_agree() -> None:
    for case in ACTIVE:
        with use_backend("python"):
            py = call(case)
        with use_backend("rust"):
            rs = call(case)
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
    kelvin = cp_of(*METHANE, 500.0)
    celsius = (
        ideal_gas_cp(
            q(METHANE[0], "cp_a"),
            q(METHANE[1], "cp_b"),
            q(METHANE[2], "cp_c"),
            q(METHANE[3], "cp_d"),
            q(METHANE[4], "cp_e"),
            Q(226.85, "degC"),
        )
        .cp.to("J/(mol*K)")
        .magnitude
    )
    h.assert_close(kelvin, celsius, 1e-9, "K against degC")
