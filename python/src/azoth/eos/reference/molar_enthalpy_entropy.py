"""``eos.molar_enthalpy_entropy`` - the absolute enthalpy and entropy of a mixture.

Spec: ``specs/models/eos/molar_enthalpy_entropy.toml``, which carries the assembly - which
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
from azoth.eos.reference.pr_molar_volume import MOLAR_GAS_CONSTANT

MODEL_ID = "eos.molar_enthalpy_entropy"

#: The temperature NeqSim's ideal-gas integrals are measured from, in kelvin.
#:
#: ``ThermodynamicConstantsInterface.referenceTemperature``. Fixed rather than a
#: caller's choice: it is where ``getHID`` and ``getIdEntropy`` both start, so it is
#: part of the correlation rather than a datum a caller supplies.
REFERENCE_TEMPERATURE = 273.15

#: The pressure NeqSim's ideal-gas entropy is measured from, in pascals.
#:
#: ``ThermodynamicConstantsInterface.referencePressure``, which is 1.01325 bar.
REFERENCE_PRESSURE = 1.01325e5


@dataclass(frozen=True, slots=True)
class IdealGasModel:
    """The ideal-gas heat-capacity coefficients of every component in the mixture.

    Five vectors rather than eight arguments because they belong together: a
    coefficient set without the other four is not a polynomial, and passing them
    separately invites a call that supplies three of the five.

    **There is no datum here, and that is NeqSim's design rather than an omission.**
    Its ``getHID`` is ``integral Cp dT`` from a fixed ``referenceTemperature`` of
    273.15 K, and it multiplies the formation enthalpy by zero; its entropy is the same
    integral of ``Cp/T`` from the same temperature, less
    ``R ln(P / referencePressure)``. So the ideal-gas state is fully determined by the
    five coefficients, and ``h_ref``, ``s_ref``, ``T_ref`` and ``P_ref`` were azoth's
    inventions.

    Attributes:
        cp_a: the constant term of each component's ``Cp``, in J/(mol*K).
        cp_b: the coefficient of ``T``, in J/(mol*K**2).
        cp_c: the coefficient of ``T**2``, in J/(mol*K**3).
        cp_d: the coefficient of ``T**3``, in J/(mol*K**4).
        cp_e: the coefficient of ``T**4``, in J/(mol*K**5).
    """

    cp_a: tuple[float, ...] | list[float]
    cp_b: tuple[float, ...] | list[float]
    cp_c: tuple[float, ...] | list[float]
    cp_d: tuple[float, ...] | list[float]
    cp_e: tuple[float, ...] | list[float]

    def validate(self, n: int) -> None:
        """Check the five vectors agree in length.

        Raises:
            InvalidInputError: if any vector's length differs from ``n``.
        """
        for name in ("cp_a", "cp_b", "cp_c", "cp_d", "cp_e"):
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

    apply_checks(
        checks.on_input,
        {
            "T": t_si,
            "P": p_si,
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

    # The ideal-gas part, which is NeqSim's `getHID` and `getIdEntropy`: each
    # component's polynomial integrated exactly from the fixed reference temperature to
    # the state. Closed forms rather than quadrature, because it is a polynomial - and a
    # five-term one, so the enthalpy gains a fifth-power term and the entropy a fourth.
    t_ref, p_ref = REFERENCE_TEMPERATURE, REFERENCE_PRESSURE
    t2, t3, t4, t5 = t_si**2, t_si**3, t_si**4, t_si**5
    r2, r3, r4, r5 = t_ref**2, t_ref**3, t_ref**4, t_ref**5
    h_ideal = 0.0
    s_ideal = 0.0
    cp_ideal = 0.0
    for i in range(n):
        a, b, c, d, e = (
            ideal_gas.cp_a[i],
            ideal_gas.cp_b[i],
            ideal_gas.cp_c[i],
            ideal_gas.cp_d[i],
            ideal_gas.cp_e[i],
        )
        # `integral Cp dT` from `T_ref` to `T`.
        h_ideal += z[i] * (
            a * (t_si - t_ref)
            + b * (t2 - r2) / 2.0
            + c * (t3 - r3) / 3.0
            + d * (t4 - r4) / 4.0
            + e * (t5 - r5) / 5.0
        )
        # `integral Cp / T dT`, the same limits.
        s_ideal += z[i] * (
            a * math.log(t_si / t_ref)
            + b * (t_si - t_ref)
            + c * (t2 - r2) / 2.0
            + d * (t3 - r3) / 3.0
            + e * (t4 - r4) / 4.0
        )
        # The polynomial itself at the state. It is the integrand of the enthalpy, so
        # reporting it costs one evaluation and no new assumption.
        cp_ideal += z[i] * (a + b * t_si + c * t2 + d * t3 + e * t4)

    # The two ideal-gas terms that no coefficient switches off.
    s_ideal -= MOLAR_GAS_CONSTANT * math.log(p_si / p_ref)
    s_ideal -= MOLAR_GAS_CONSTANT * sum(zi * math.log(zi) for zi in z)

    h_departure = MOLAR_GAS_CONSTANT * t_si * state.h_dep_rt
    s_departure = MOLAR_GAS_CONSTANT * state.s_dep_r
    cp_departure = MOLAR_GAS_CONSTANT * state.cp_dep_r

    return MolarEnthalpyEntropyResult(
        h=from_si(h_ideal + h_departure, "J/mol"),
        s=from_si(s_ideal + s_departure, "J/(mol*K)"),
        h_ideal=from_si(h_ideal, "J/mol"),
        s_ideal=from_si(s_ideal, "J/(mol*K)"),
        h_departure=from_si(h_departure, "J/mol"),
        s_departure=from_si(s_departure, "J/(mol*K)"),
        psi_bar=state.psi_bar,
        cp=from_si(cp_ideal + cp_departure, "J/(mol*K)"),
        cp_ideal=from_si(cp_ideal, "J/(mol*K)"),
        cp_departure=from_si(cp_departure, "J/(mol*K)"),
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
