"""``eos.wax_solid_fugacity`` - a wax cut's solid fugacity coefficient.

```text
phi_wax = phi_liq(T, P) exp( -dH_fus/(R T) (1 - T/T_tp)
                           + dCp_SL/R (T_tp/T - 1 - ln(T_tp/T))
                           - (v_liq - v_sol)(P - P_ref)/(R T) )
```

Spec: ``specs/calcs/eos/wax_solid_fugacity.toml``, which carries the fits, the one-bar
reference pressure, and the measurement that pins all of it.

NeqSim's ``ComponentWax.fugcoef2``, the default ``PhaseWax`` component model. **The
coefficient is a pure-component quantity**: ``SolidFug = x f_liq exp(...)`` and the reported
coefficient is ``SolidFug/(P x)``, so the mole fraction cancels.

**The units are where this goes wrong.** With ``v`` in m³/mol against a pressure in bar the
volume term comes out ``1e5`` too small and the coefficient is out by a per cent that grows
with the cut's molar mass; NeqSim's ``refPressure`` is ``1.0`` in its bar-valued code, so in
pascals the term is ``(P - 1.0e5)``.
"""

from __future__ import annotations

import math

from azoth._registry_gen import spec as _spec_for
from azoth.core.errors import InvalidInputError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import WaxSolidFugacityResult
from azoth.core.units import Q, input_to_si
from azoth.core.warnings import Warning
from azoth.eos.cubic import PR, SRK
from azoth.eos.mixture import Component, Mixture
from azoth.eos.reference._mixture_state import phase_state, reduced_parameters

CALC_ID = "eos.wax_solid_fugacity"

#: NeqSim's ``R``, which the class's compiled form states as a literal.
R = 8.3144621

#: One bar, the pressure the volume term is referred to.
REFERENCE_PRESSURE = 1.0e5

#: The solid's molar volume over the liquid's, which is NeqSim's own shortcut.
SOLID_VOLUME_RATIO = 0.9


def wax_solid_fugacity(
    molar_mass: Q,
    tc: Q,
    pc: Q,
    omega: float,
    heat_of_fusion: Q,
    triple_point_temperature: Q,
    T: Q,
    P: Q,
    eos: str = "srk",
) -> WaxSolidFugacityResult:
    """A wax cut's solid fugacity coefficient at a state.

    Args:
        molar_mass: the cut's molar mass.
        tc: the cut's critical temperature, which the reference liquid is a one-component
            phase at.
        pc: the cut's critical pressure.
        omega: the cut's acentric factor, which sets the reference liquid's alpha.
        heat_of_fusion: the cut's heat of fusion.
        triple_point_temperature: the cut's triple-point temperature.
        T: absolute temperature.
        P: absolute pressure.
        eos: the cubic the fluid runs, which is also the one the reference liquid is built
            from - NeqSim clones the host phase's class, so this is a real choice.

    Returns:
        The coefficient, dimensionless, and any caveats.

    Raises:
        InvalidInputError: if ``eos`` is not a cubic this reaches.
        OutOfRangeError: if any of the cut's constants or the state is not positive.
    """
    spec = _spec_for(CALC_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []
    mass_si = input_to_si(spec, "molar_mass", molar_mass)
    fusion_si = input_to_si(spec, "heat_of_fusion", heat_of_fusion)
    tc_si = input_to_si(spec, "tc", tc)
    pc_si = input_to_si(spec, "pc", pc)
    ttp_si = input_to_si(spec, "triple_point_temperature", triple_point_temperature)
    t_si = input_to_si(spec, "T", T)
    p_si = input_to_si(spec, "P", P)
    apply_checks(
        checks.on_input,
        {
            "molar_mass": mass_si,
            "tc": tc_si,
            "pc": pc_si,
            "heat_of_fusion": fusion_si,
            "triple_point_temperature": ttp_si,
            "T": t_si,
            "P": p_si,
        }.get,
        warnings,
    )

    if eos not in ("pr", "srk"):
        raise InvalidInputError(
            "eos",
            f"`{eos}` is not a cubic this reaches: NeqSim builds the reference liquid from "
            "the *host phase's* class, and only `pr` and `srk` have a wax route",
        )

    # **The reference liquid: this component alone, on the fluid's own cubic, on the liquid
    # root.** NeqSim's `setSolidRefFluidPhase` clones the host phase's class and adds one
    # component, so the alpha is that cubic's own Soave form at the cut's acentric factor.
    reference = Mixture(
        components=(Component(Tc=tc, Pc=pc, omega=omega, molar_mass=molar_mass),),
        cubic=PR if eos == "pr" else SRK,
        alpha="pr" if eos == "pr" else "srk",
    )
    reduced = reduced_parameters(reference, t_si, p_si)
    state = phase_state(reduced, reference.kij, [1.0], liquid=True)
    phi_liquid = math.exp(state.ln_phi[0])
    v_liquid = state.z * R * t_si / p_si

    v_solid = SOLID_VOLUME_RATIO * v_liquid
    pressure_term = -(v_liquid - v_solid) * (p_si - REFERENCE_PRESSURE) / R / t_si
    molar_mass_grams = mass_si * 1000.0
    delta_cp_sl = (0.3033 * molar_mass_grams - 4.635e-4 * molar_mass_grams * t_si) * 4.184
    triple_point_ratio = ttp_si / t_si
    heat_capacity_term = delta_cp_sl / R * (triple_point_ratio - 1.0 - math.log(triple_point_ratio))
    fusion_term = -fusion_si / (R * t_si) * (1.0 - t_si / ttp_si)

    return WaxSolidFugacityResult(
        fugacity_coefficient=phi_liquid
        * math.exp(fusion_term + heat_capacity_term + pressure_term),
        warnings=tuple(warnings),
    )


__all__ = ["wax_solid_fugacity"]
