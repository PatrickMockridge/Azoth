"""Spec-driven tests for the ``eos.molar_enthalpy_entropy`` model."""

from __future__ import annotations

from typing import Any

import pytest

import _helpers as h
from azoth import _models_gen, ureg, use_backend
from azoth.core.errors import InvalidInputError, OutOfRangeError
from azoth.core.result import MolarEnthalpyEntropyResult
from azoth.eos import (
    Component,
    IdealGasModel,
    Mixture,
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

METHANE = Component(Q(190.56, "K"), Q(4_599_200.0, "Pa"), 0.01142)
BUTANE = Component(Q(425.12, "K"), Q(3_796_000.0, "Pa"), 0.2002)


def methane_butane() -> Mixture:
    return mixture([METHANE, BUTANE], kij={(0, 1): 0.05})


def ideal_gas(**overrides: Any) -> IdealGasModel:
    fields: dict[str, Any] = {
        "cp_a": (4.0, 4.0),
        "cp_b": (1.0, 1.0),
        "cp_c": (-0.5, -0.5),
        "cp_d": (0.1, 0.1),
        "h_ref": (0.0, 0.0),
        "s_ref": (0.0, 0.0),
        "T_ref": Q(298.15, "K"),
        "P_ref": Q(101325.0, "Pa"),
    }
    fields.update(overrides)
    return IdealGasModel(**fields)


def call(case: dict[str, Any]) -> MolarEnthalpyEntropyResult:
    inputs = case["inputs"]
    fluid = Mixture(
        components=tuple(
            Component(Q(tc, "K"), Q(pc, "Pa"), omega)
            for tc, pc, omega in zip(inputs["Tc"], inputs["Pc"], inputs["omega"], strict=True)
        ),
        kij=tuple(tuple(row) for row in inputs["kij"]),
    )
    return molar_enthalpy_entropy(
        fluid,
        ideal_gas(
            cp_a=tuple(inputs["cp_a"]),
            cp_b=tuple(inputs["cp_b"]),
            cp_c=tuple(inputs["cp_c"]),
            cp_d=tuple(inputs["cp_d"]),
            h_ref=tuple(inputs["h_ref"]),
            s_ref=tuple(inputs["s_ref"]),
            T_ref=Q(inputs["T_ref"], "K"),
            P_ref=Q(inputs["P_ref"], "Pa"),
        ),
        T=Q(inputs["T"], "K"),
        P=Q(inputs["P"], "Pa"),
        z=list(inputs["z"]),
        compressibility=inputs["compressibility"],
    )


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

    for component, t_c, p_pa in ((METHANE, 200.0, 3.0e6), (BUTANE, 350.0, 1.0e6)):
        fluid = mixture([component])
        reduced = reduced_parameters(fluid, t_c, p_pa)
        state = phase_state(reduced, fluid.kij, [1.0], liquid=False)

        result = molar_enthalpy_entropy(
            fluid,
            ideal_gas(
                cp_a=(4.0,),
                cp_b=(1.0,),
                cp_c=(-0.5,),
                cp_d=(0.1,),
                h_ref=(0.0,),
                s_ref=(0.0,),
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
            pr_kappa(component.omega).kappa,
            t_c / component.Tc.to_base_units().magnitude,
        )

        assert result.psi_bar == reduced.psi[0], "psi_bar should be the component's own psi"
        assert (
            result.h_departure.to("J/mol").magnitude / (MOLAR_GAS_CONSTANT * t_c) == pure.h_dep_rt
        ), "the departure enthalpy should be bit-identical"
        h.assert_close(
            result.s_departure.to("J/(mol*K)").magnitude / MOLAR_GAS_CONSTANT,
            pure.s_dep_r,
            1e-14,
            "the entropy departures",
        )


def test_the_datum_shifts_the_enthalpy_and_leaves_the_entropy_alone() -> None:
    """The reference is an enthalpy datum, and nothing else.

    Raising every ``h_ref`` by a constant must raise ``h`` by exactly that constant and
    leave ``s`` and both departures untouched. It is what would fail if the reference
    were mixed into the integral rather than added to it - invisible at a zero datum,
    which is the only datum a spec case can state.
    """
    fluid = methane_butane()
    args: dict[str, Any] = {
        "T": Q(330.0, "K"),
        "P": Q(2.5e6, "Pa"),
        "z": [0.6, 0.4],
        "compressibility": 0.8274482588400789,
    }
    plain = molar_enthalpy_entropy(fluid, ideal_gas(), **args)
    shifted = molar_enthalpy_entropy(fluid, ideal_gas(h_ref=(-285_830.0, -285_830.0)), **args)

    h.assert_close(
        shifted.h.to("J/mol").magnitude,
        plain.h.to("J/mol").magnitude - 285_830.0,
        1e-9,
        "the datum should shift h by exactly what was put in",
    )
    assert shifted.h_departure.to("J/mol").magnitude == plain.h_departure.to("J/mol").magnitude
    assert shifted.s.to("J/(mol*K)").magnitude == plain.s.to("J/(mol*K)").magnitude


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
    with pytest.raises(OutOfRangeError):
        molar_enthalpy_entropy(fluid, ideal_gas(T_ref=Q(0.0, "K")), **args)


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
