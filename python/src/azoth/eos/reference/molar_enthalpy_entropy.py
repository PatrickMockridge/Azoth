"""``eos.molar_enthalpy_entropy`` - the absolute enthalpy and entropy of a mixture.

Spec: ``specs/models/eos/molar_enthalpy_entropy.yaml``, which carries the assembly - which
term belongs on which side - and what the caller is responsible for: the datum, the
coefficients and the compressibility factor.

.. code-block:: text

    H = H_ig(T_ref) + integral Cp dT      + H_dep
    S = S_ig(T_ref) + integral Cp/T dT    - R ln(P/P_ref) - R sum z_i ln z_i + S_dep

The departure functions are ``_mixture_state``'s and the integrals are exact integrals of
the polynomial :func:`azoth.eos.ideal_gas_cp` registers.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from azoth.core.errors import InvalidInputError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import MolarEnthalpyEntropyResult
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning
from azoth.eos.mixture import Mixture
from azoth.eos.reference._mixture_state import phase_state_at, reduced_parameters
from azoth.eos.reference.ideal_gas_cp import REFERENCE_TEMPERATURE
from azoth.eos.reference.pr_molar_volume import MOLAR_GAS_CONSTANT

MODEL_ID = "eos.molar_enthalpy_entropy"


@dataclass(frozen=True, slots=True)
class IdealGasModel:
    """The caller's ideal-gas model and the datum it is referenced to.

    One object rather than eight arguments because the eight belong together: a
    coefficient set without a reference state is not a thermodynamic model, and
    passing them separately invites a call that supplies four of the six vectors.

    Attributes:
        cp_a: the constant term of each component's ``Cp/R`` polynomial.
        cp_b: the coefficient of ``theta`` in each component's polynomial.
        cp_c: the coefficient of ``theta**2``.
        cp_d: the coefficient of ``theta**3``.
        h_ref: each component's ideal-gas molar enthalpy at ``T_ref``, in J/mol.
        s_ref: each component's ideal-gas molar entropy at ``T_ref`` and ``P_ref``,
            in J/(mol*K).
        T_ref: the temperature the reference values are given at.
        P_ref: the pressure ``s_ref`` is given at. It does not enter the enthalpy.
    """

    cp_a: tuple[float, ...] | list[float]
    cp_b: tuple[float, ...] | list[float]
    cp_c: tuple[float, ...] | list[float]
    cp_d: tuple[float, ...] | list[float]
    h_ref: tuple[float, ...] | list[float]
    s_ref: tuple[float, ...] | list[float]
    T_ref: Q
    P_ref: Q

    def validate(self, n: int) -> None:
        """Check the six vectors agree in length.

        Raises:
            InvalidInputError: if any vector's length differs from ``n``.
        """
        for name in ("cp_a", "cp_b", "cp_c", "cp_d", "h_ref", "s_ref"):
            vector = getattr(self, name)
            if len(vector) != n:
                raise InvalidInputError(
                    name,
                    f"a mixture of {n} components needs one value per component, but "
                    f"{name} has {len(vector)}",
                )


def molar_enthalpy_entropy(
    mixture: Mixture,
    ideal_gas: IdealGasModel,
    T: Q,
    P: Q,
    z: list[float],
    compressibility: float,
) -> MolarEnthalpyEntropyResult:
    """The absolute molar enthalpy and entropy of a mixture at a state.

    Args:
        mixture: the components and their interaction parameters.
        ideal_gas: the heat-capacity coefficients, the reference values and the
            reference state. **This is the datum**, and two results computed from
            different ones are not comparable.
        T: absolute temperature of the state.
        P: absolute pressure of the state.
        z: the mixture's mole fractions, checked rather than renormalised.
        compressibility: the cubic's root for the phase wanted. Not solved for here:
            which root describes the phase is a choice, and a caller holding a ``Z``
            has already made it.

    Returns:
        The enthalpy and entropy, split into their ideal-gas and departure parts.

    Raises:
        InvalidInputError: if the composition or any of the ideal-gas vectors is the
            wrong length, or if ``z`` is not a composition.
        OutOfRangeError: if ``T``, ``P`` or either reference coordinate is not
            positive, or if ``compressibility`` is not admissible.
    """
    from azoth import _models_gen

    spec = _models_gen.model(MODEL_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    t_si = input_to_si(spec, "T", T)
    p_si = input_to_si(spec, "P", P)
    t_ref = input_to_si(spec, "T_ref", ideal_gas.T_ref)
    p_ref = input_to_si(spec, "P_ref", ideal_gas.P_ref)

    apply_checks(
        checks.on_input,
        {
            "T": t_si,
            "P": p_si,
            "T_ref": t_ref,
            "P_ref": p_ref,
            "compressibility": compressibility,
        }.get,
        warnings,
    )

    n = len(mixture)
    ideal_gas.validate(n)
    _check_composition(z, n)

    reduced = reduced_parameters(mixture, t_si, p_si)
    warnings.extend(reduced.warnings)
    state = phase_state_at(reduced, mixture.kij, list(z), compressibility)

    # The ideal-gas part: each component's reference value plus the exact integral of
    # its polynomial. Closed forms rather than quadrature, because it is a polynomial.
    theta = t_si / REFERENCE_TEMPERATURE
    theta_ref = t_ref / REFERENCE_TEMPERATURE
    h_ideal = 0.0
    s_ideal = 0.0
    for i in range(n):
        a, b, c, d = (
            ideal_gas.cp_a[i],
            ideal_gas.cp_b[i],
            ideal_gas.cp_c[i],
            ideal_gas.cp_d[i],
        )
        # `T = REFERENCE_TEMPERATURE * theta`, so `dT = REFERENCE_TEMPERATURE dtheta`.
        dh = (
            MOLAR_GAS_CONSTANT
            * REFERENCE_TEMPERATURE
            * (
                a * (theta - theta_ref)
                + b * (theta**2 - theta_ref**2) / 2.0
                + c * (theta**3 - theta_ref**3) / 3.0
                + d * (theta**4 - theta_ref**4) / 4.0
            )
        )
        # `integral Cp/T dT` is `R * integral (a + b theta + ...)/theta dtheta`.
        ds = MOLAR_GAS_CONSTANT * (
            a * math.log(theta / theta_ref)
            + b * (theta - theta_ref)
            + c * (theta**2 - theta_ref**2) / 2.0
            + d * (theta**3 - theta_ref**3) / 3.0
        )
        h_ideal += z[i] * (ideal_gas.h_ref[i] + dh)
        s_ideal += z[i] * (ideal_gas.s_ref[i] + ds)

    # The two ideal-gas terms that no coefficient switches off.
    s_ideal -= MOLAR_GAS_CONSTANT * math.log(p_si / p_ref)
    s_ideal -= MOLAR_GAS_CONSTANT * sum(zi * math.log(zi) for zi in z)

    h_departure = MOLAR_GAS_CONSTANT * t_si * state.h_dep_rt
    s_departure = MOLAR_GAS_CONSTANT * state.s_dep_r

    return MolarEnthalpyEntropyResult(
        h=from_si(h_ideal + h_departure, "J/mol"),
        s=from_si(s_ideal + s_departure, "J/(mol*K)"),
        h_ideal=from_si(h_ideal, "J/mol"),
        s_ideal=from_si(s_ideal, "J/(mol*K)"),
        h_departure=from_si(h_departure, "J/mol"),
        s_departure=from_si(s_departure, "J/(mol*K)"),
        psi_bar=state.psi_bar,
        warnings=tuple(warnings),
    )


def _check_composition(values: list[float], n: int) -> None:
    if len(values) != n:
        raise InvalidInputError(
            "z", f"a mixture of {n} components needs {n} mole fractions, but z has {len(values)}"
        )
    for i, value in enumerate(values):
        if value < 0.0:
            raise InvalidInputError(
                "z", f"z[{i}] is {value} but a mole fraction cannot be negative"
            )
    if abs(sum(values) - 1.0) > 1.0e-09:
        raise InvalidInputError(
            "z",
            f"the composition sums to {sum(values)}, not to one. Renormalising it here "
            f"would make a caller's error invisible in every number downstream, so it "
            f"is refused instead",
        )
