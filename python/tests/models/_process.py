"""Shared fixtures for the unit-operation tests.

Every unit operation in :mod:`azoth.process` is called with the same things - a mixture,
an ideal-gas datum, a temperature, a pressure, a molar flow and a composition - so the
eight test modules that exercise them would otherwise each build the same objects and
each run the same ``test_spec_case`` loop.

**A helper module rather than a ``conftest`` fixture, on purpose.** These tests do not
want a fixture's lifetime or its indirection; they want a function they can read, and one
that shows exactly which numbers a model is being given.

The two components are the ones every other model spec in this tree uses -
methane/n-butane with the same ``kij`` - so a reader comparing a unit operation's case
against ``eos.ph_flash``'s is comparing the same mixture.
"""

from __future__ import annotations

from typing import Any

import _helpers as h
from azoth import ureg
from azoth.core.result import Phase
from azoth.eos import Component, IdealGasModel, Mixture, mixture

Q = ureg.Quantity

METHANE = Component(Q(190.56, "K"), Q(4_599_000.0, "Pa"), 0.0115)
BUTANE = Component(Q(425.12, "K"), Q(3_796_000.0, "Pa"), 0.2002)

#: The units a unit operation's *stream* arguments are declared in, by argument name.
#: Only the ones that are not `dimensionless` need an entry; a spec that renamed one of
#: these would fail the contract tests rather than silently converting with the wrong
#: unit.
#:
#: `T_ref` and `P_ref` are deliberately **absent**. They are the ideal-gas datum's
#: reference state rather than arguments to the unit operation - they go to
#: `IdealGasModel`, which is why `an_ideal_gas` builds them - and an entry here would
#: pass them to a function that does not take them.
UNITS = {
    "T": "K",
    "P": "Pa",
    "n": "mol/s",
    "pressure_drop": "Pa",
    "heat_duty": "W",
    "outlet_pressure": "Pa",
}

#: The unit each result field is a quantity in, for the ones that are quantities at all.
#: Everything else a result carries - a flow, a power, a fraction, a composition - is a
#: plain number in its SI base unit, and comparing it needs no conversion.
RESULT_UNITS = {
    "T": "K",
    "P": "Pa",
    "isentropic_temperature": "K",
}


def a_mixture() -> Mixture:
    """The mixture the spec cases describe."""
    return mixture([METHANE, BUTANE], kij={(0, 1): 0.01289789})


def an_ideal_gas() -> IdealGasModel:
    """A deliberately trivial ideal-gas model: five zero coefficients.

    Nothing about a unit operation depends on the coefficients being physical - its
    arithmetic is the same whatever they are - so they are chosen to make a failure
    readable rather than to describe a substance.
    """
    return IdealGasModel(
        cp_a=(3.0, 5.0),
        cp_b=(0.0, 0.0),
        cp_c=(0.0, 0.0),
        cp_d=(0.0, 0.0),
        cp_e=(0.0, 0.0),
    )


def a_quantity(name: str, value: Any) -> Any:
    """One state argument as a quantity, or as a list of them if the spec declares one.

    The vector case is what makes a multi-inlet unit testable from its spec: `T`, `P` and
    `n` are one entry per inlet in `process.mixer`, and the same helper has to produce
    either shape without the caller saying which.
    """
    unit = UNITS[name]
    if isinstance(value, list):
        return [Q(item, unit) for item in value]
    return Q(value, unit)


def call(model: dict[str, Any], case: dict[str, Any]) -> Any:
    """Run one case declared in a unit-operation spec, through the public API.

    Builds the arguments from the case's own ``inputs`` rather than from a hand-written
    call, so a case that gains an input cannot leave the test calling an older shape.
    """
    import azoth

    inputs = case["inputs"]
    arguments: dict[str, Any] = {"mixture": a_mixture()}
    if "cp_a" in inputs:
        arguments["ideal_gas"] = an_ideal_gas()
    for name, value in inputs.items():
        if name in UNITS:
            arguments[name] = a_quantity(name, value)
        elif name in ("z", "fractions", "efficiency"):
            arguments[name] = value
    return getattr(azoth.process, model["id"].split(".")[1])(**arguments)


def assert_case(model: dict[str, Any], case: dict[str, Any]) -> None:
    """Run one case and compare every recorded value against the result.

    The loop every unit-operation test module runs, in one place. A `T` or a `P` is
    compared as a magnitude in its own unit; anything else is compared as it stands,
    including a list.
    """
    result = call(model, case)
    for name, expected in case["expected"].items():
        got = getattr(result, name)
        if name in RESULT_UNITS:
            got = got.to(RESULT_UNITS[name]).magnitude
        if isinstance(expected, list):
            assert len(got) == len(expected), f"{case['id']} ({name}): length differs"
            for index, (a, b) in enumerate(zip(got, expected, strict=True)):
                h.assert_close(a, b, case["tolerance"], f"{case['id']} ({name}[{index}])")
        else:
            h.assert_close(got, expected, case["tolerance"], f"{case['id']} ({name})")


__all__ = [
    "BUTANE",
    "METHANE",
    "UNITS",
    "Phase",
    "a_mixture",
    "a_quantity",
    "an_ideal_gas",
    "assert_case",
    "call",
    "h",
]
