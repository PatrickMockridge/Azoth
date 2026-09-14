"""Spec-driven tests for the ``eos.molar_enthalpy_entropy`` model."""

from __future__ import annotations

from typing import Any

import pytest

import _helpers as h
from azoth import _models_gen, ureg, use_backend
from azoth.core.errors import InvalidInputError, OutOfRangeError
from azoth.core.result import MolarEnthalpyEntropyResult
from azoth.eos import (
    IdealGasModel,
    Mixture,
    component,
    from_names,
    mixture,
    molar_enthalpy_entropy,
    pr_departure,
    pr_kappa,
)
from azoth.eos.reference.pr_molar_volume import MOLAR_GAS_CONSTANT

MODEL_ID = "eos.molar_enthalpy_entropy"
Q = ureg.Quantity

SPEC = _models_gen.model(MODEL_ID)
CASES = SPEC["cases"]

# The substances the cases and the identities below use, resolved through the
# databank rather than typed here. A `Component` written out longhand is a second
# copy of NeqSim's table, and this one had drifted: its methane was 0.01142 and
# 4 599 200 Pa, where COMP.csv says 0.0115 and 4 599 000.
METHANE = component("methane")
BUTANE = component("n-butane")


def methane_butane() -> Mixture:
    """The methane/n-butane pair, with the interaction parameter NeqSim fits for it.

    Resolved by name rather than assembled here, so the pair a sweep runs and the pair
    a case runs are the same fluid. The `kij` used to be an "illustrative" 0.05 - a
    number in a test file that described no fluid, and that disagreed with `INTER.csv`.
    """
    return from_names(["methane", "n-butane"])


def ideal_gas(**overrides: Any) -> IdealGasModel:
    fields: dict[str, Any] = {
        "cp_a": (4.0, 4.0),
        "cp_b": (1.0, 1.0),
        "cp_c": (-0.5, -0.5),
        "cp_d": (0.1, 0.1),
        "cp_e": (0.0, 0.0),
    }
    fields.update(overrides)
    return IdealGasModel(**fields)


def call(case: dict[str, Any]) -> MolarEnthalpyEntropyResult:
    """Run one case declared in the model spec.

    Through :func:`_helpers.model_kwargs`, which is the one place a case's
    declared inputs become arguments: it resolves `components` against the
    databank and hands over the mixture and the ideal-gas model the function
    takes. A hand-built mixture here would be a second fluid, described by the
    case file rather than by NeqSim's tables.
    """
    return molar_enthalpy_entropy(**h.model_kwargs(SPEC, case["inputs"]))


@pytest.mark.parametrize("case", CASES, ids=lambda c: c["id"])
def test_spec_case(case: dict[str, Any]) -> None:
    result = call(case)
    for name in ("h", "s", "h_ideal", "s_ideal", "h_departure", "s_departure"):
        h.assert_close(
            getattr(result, name).to("J/mol" if name.startswith("h") else "J/(mol*K)").magnitude,
            case["expected"][name],
            case["tolerance"],
            f"{case['id']} ({name})",
        )
    h.assert_close(result.psi_bar, case["expected"]["psi_bar"], 1e-12, f"{case['id']} (psi_bar)")
    h.assert_consistent(result, case["id"])


def test_every_case_ran() -> None:
    assert len(CASES) >= 2, f"expected several cases, found {len(CASES)}"


def test_the_enthalpy_is_the_departure_when_the_ideal_gas_terms_are_off() -> None:
    """The identity that ties this model to ``eos.pr_departure``: ``h = R*T*h_dep_rt``.

    Exact rather than approximate, because both sides are the same product - so a
    defect in the departure arithmetic or an extra term in the enthalpy shows up as a
    disagreement rather than as a tolerance.
    """
    fluid = methane_butane()
    off = ideal_gas(cp_a=(0.0, 0.0), cp_b=(0.0, 0.0), cp_c=(0.0, 0.0), cp_d=(0.0, 0.0))
    result = molar_enthalpy_entropy(
        fluid,
        off,
        T=Q(330.0, "K"),
        P=Q(2.5e6, "Pa"),
        z=[0.6, 0.4],
        compressibility=0.8274482588400789,
    )
    assert result.h_ideal.to("J/mol").magnitude == 0.0
    h.assert_close(
        result.h.to("J/mol").magnitude,
        result.h_departure.to("J/mol").magnitude,
        1e-12,
        "h should be the departure",
    )

    # And the entropy is the departure *plus* the two ideal-gas terms that no
    # coefficient switches off - the pair a reader might mistake for a defect.
    import math

    pressure_and_mixing = -MOLAR_GAS_CONSTANT * math.log(2.5e6 / 101325.0) - MOLAR_GAS_CONSTANT * (
        0.6 * math.log(0.6) + 0.4 * math.log(0.4)
    )
    h.assert_close(result.s_ideal.to("J/(mol*K)").magnitude, pressure_and_mixing, 1e-9, "s_ideal")
    h.assert_close(
        result.s.to("J/(mol*K)").magnitude,
        pressure_and_mixing + result.s_departure.to("J/(mol*K)").magnitude,
        1e-9,
        "s",
    )


def test_the_departures_reduce_to_pr_departure_at_one_component() -> None:
    """At ``N = 1`` the departures are the registered pure-component calc's.

    `h_departure / (R*T)` must equal ``pr_departure``'s ``h_dep_rt`` exactly - the
    mixture form is the pure form with ``psi`` replaced by an average that at one
    component is that component's own ``psi``.
    """
    from azoth.eos.reference._mixture_state import phase_state, reduced_parameters

    for substance, t_c, p_pa in ((METHANE, 200.0, 3.0e6), (BUTANE, 350.0, 1.0e6)):
        fluid = mixture([substance])
        reduced = reduced_parameters(fluid, t_c, p_pa)
        state = phase_state(reduced, fluid.kij, [1.0], liquid=False)

        result = molar_enthalpy_entropy(
            fluid,
            ideal_gas(
                cp_a=(4.0,),
                cp_b=(1.0,),
                cp_c=(-0.5,),
                cp_d=(0.1,),
                cp_e=(0.0,),
            ),
            T=Q(t_c, "K"),
            P=Q(p_pa, "Pa"),
            z=[1.0],
            compressibility=state.z,
        )
        pure = pr_departure(
            reduced.a[0],
            reduced.b[0],
            state.z,
            pr_kappa(substance.omega).kappa,
            t_c / substance.Tc.to_base_units().magnitude,
        )

        assert result.psi_bar == reduced.psi[0], "psi_bar should be the substance's own psi"
        assert (
            result.h_departure.to("J/mol").magnitude / (MOLAR_GAS_CONSTANT * t_c) == pure.h_dep_rt
        ), "the departure enthalpy should be bit-identical"
        h.assert_close(
            result.s_departure.to("J/(mol*K)").magnitude / MOLAR_GAS_CONSTANT,
            pure.s_dep_r,
            1e-14,
            "the entropy departures",
        )


def test_the_ideal_gas_enthalpy_is_zero_at_neqs_im_reference_temperature() -> None:
    """``integral Cp dT`` from `T_ref` to `T_ref` is zero, whatever the coefficients.

    ``h_ideal`` is measured from NeqSim's fixed ``referenceTemperature`` of 273.15 K, so
    at that temperature the ideal-gas enthalpy is exactly nothing - and the entropy is
    then only the pressure and mixing terms. This is the check that would fail if a
    reference were mixed into the integral rather than being its lower limit, and it
    pins the constant: 298.15 K, a rounder number, gives a non-zero answer.
    """
    fluid = methane_butane()
    args: dict[str, Any] = {
        "P": Q(2.5e6, "Pa"),
        "z": [0.6, 0.4],
        "compressibility": 0.8274482588400789,
    }
    reference = molar_enthalpy_entropy(fluid, ideal_gas(), T=Q(273.15, "K"), **args)
    assert reference.h_ideal.to("J/mol").magnitude == 0.0

    other = molar_enthalpy_entropy(fluid, ideal_gas(), T=Q(298.15, "K"), **args)
    assert abs(other.h_ideal.to("J/mol").magnitude) > 1.0
    assert other.h_departure.to("J/mol").magnitude != reference.h_departure.to("J/mol").magnitude


def test_a_malformed_input_is_refused() -> None:
    fluid = methane_butane()
    args: dict[str, Any] = {
        "T": Q(330.0, "K"),
        "P": Q(2.5e6, "Pa"),
        "z": [0.6, 0.4],
        "compressibility": 0.8274482588400789,
    }

    with pytest.raises(InvalidInputError):
        molar_enthalpy_entropy(fluid, ideal_gas(cp_b=(1.0,)), **args)  # wrong length
    for z in ([0.6, 0.5], [0.6, -0.6], [0.6, 0.4, 0.0]):
        with pytest.raises(InvalidInputError):
            molar_enthalpy_entropy(fluid, ideal_gas(), **{**args, "z": z})
    with pytest.raises(OutOfRangeError) as excinfo:
        molar_enthalpy_entropy(fluid, ideal_gas(), **{**args, "compressibility": 0.01})
    assert excinfo.value.field() == "z"


def test_the_two_backends_agree_on_every_spec_case() -> None:
    for case in CASES:
        with use_backend("python"):
            py = call(case)
        with use_backend("rust"):
            rs = call(case)
        for name in ("h", "s", "h_ideal", "s_ideal", "h_departure", "s_departure"):
            h.assert_close(
                getattr(py, name).to("J/mol" if name.startswith("h") else "J/(mol*K)").magnitude,
                getattr(rs, name).to("J/mol" if name.startswith("h") else "J/(mol*K)").magnitude,
                1e-12,
                f"{case['id']} ({name})",
            )
        h.assert_close(py.psi_bar, rs.psi_bar, 1e-15, f"{case['id']} (psi_bar)")


def test_the_cases_are_announced_in_the_model_docs() -> None:
    from pathlib import Path

    root = Path(__file__).resolve().parents[3]
    page = root / "docs" / "src" / "eos" / "molar_enthalpy_entropy.md"
    assert page.is_file(), f"{page} was not generated"
    summary = (root / "docs" / "src" / "SUMMARY.md").read_text(encoding="utf-8")
    assert "eos/molar_enthalpy_entropy.md" in summary
