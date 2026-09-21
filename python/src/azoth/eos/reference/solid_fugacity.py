"""``eos.solid_fugacity`` - a pure solid's fugacity coefficient, from tabulated properties.

```text
phi_solid = phi_liq(T, P) exp( -dH_fus/(R T) (1 - T/T_tp)
                             + dCp_SL/(R T) (T_tp - T)
                             - dCp_SL/R ln(T_tp/T)
                             - dV_sl (P_bar - 1)/(R T) )
```

Spec: ``specs/calcs/eos/solid_fugacity.toml``, which carries the heat-capacity difference's
two routes and the measurement that pins the whole expression.

NeqSim's ``ComponentSolid.fugcoef2``. **The reference is a *liquid*** - ``fugcoef2``
initialises it ``PhaseType.LIQUID``, where the class's other entry (``fugcoef``) uses a gas at
the component's solid vapour pressure and is a different model - and it is the host's own
class, so ``eos`` is an input.
"""

from __future__ import annotations

import math

from azoth._registry_gen import spec as _spec_for
from azoth.core.errors import InvalidInputError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import SolidFugacityResult
from azoth.core.units import Q, input_to_si
from azoth.core.warnings import Warning
from azoth.eos.cubic import PR, SRK
from azoth.eos.mixture import Component, Mixture
from azoth.eos.reference._mixture_state import phase_state, reduced_parameters

CALC_ID = "eos.solid_fugacity"

#: NeqSim's ``R``, which ``Component`` states as a literal.
R = 8.3144621

#: One atmosphere in bar: the pressure the volume term is referred to in NeqSim's own code.
REFERENCE_PRESSURE_BAR = 1.0


def solid_fugacity(
    heat_of_fusion: Q,
    triple_point_temperature: Q,
    delta_cp_sl: Q,
    delta_solid_volume: Q,
    tc: Q,
    pc: Q,
    omega: float,
    T: Q,
    P: Q,
    eos: str = "srk",
) -> SolidFugacityResult:
    """A pure solid's fugacity coefficient at a state.

    Args:
        heat_of_fusion: the component's heat of fusion.
        triple_point_temperature: its triple-point temperature.
        delta_cp_sl: the solid-liquid heat-capacity difference **at the triple point**, which
            the caller states: NeqSim takes it from the table and overrides it to ``37.12``
            for water.
        delta_solid_volume: the solid's molar volume less the liquid's, at the triple point.
        tc: the component's critical temperature.
        pc: its critical pressure.
        omega: its acentric factor.
        T: absolute temperature.
        P: absolute pressure.
        eos: the cubic the fluid runs, which is also the one the reference liquid is built
            from.

    Raises:
        InvalidInputError: if ``eos`` is not a cubic this reaches.
        OutOfRangeError: if any of the constants or the state is not positive.
    """
    spec = _spec_for(CALC_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []
    fusion_si = input_to_si(spec, "heat_of_fusion", heat_of_fusion)
    ttp_si = input_to_si(spec, "triple_point_temperature", triple_point_temperature)
    dcp_si = input_to_si(spec, "delta_cp_sl", delta_cp_sl)
    dvol_si = input_to_si(spec, "delta_solid_volume", delta_solid_volume)
    t_si = input_to_si(spec, "T", T)
    p_si = input_to_si(spec, "P", P)
    apply_checks(
        checks.on_input,
        {
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
            f"`{eos}` is not a cubic this reaches: NeqSim builds the reference liquid from the "
            "*host phase's* class, and only `pr` and `srk` are ported",
        )

    reference = Mixture(
        components=(Component(Tc=tc, Pc=pc, omega=omega),),
        cubic=PR if eos == "pr" else SRK,
        alpha="pr" if eos == "pr" else "srk",
    )
    reduced = reduced_parameters(reference, t_si, p_si)
    state = phase_state(reduced, reference.kij, [1.0], liquid=True)
    phi_liquid = math.exp(state.ln_phi[0])

    pressure_bar = p_si / 1.0e5
    fusion = -fusion_si / (R * t_si) * (1.0 - t_si / ttp_si)
    heat_capacity = dcp_si / (R * t_si) * (ttp_si - t_si) - dcp_si / R * math.log(ttp_si / t_si)
    volume = -dvol_si * (pressure_bar - REFERENCE_PRESSURE_BAR) / (R * t_si)

    return SolidFugacityResult(
        fugacity_coefficient=phi_liquid * math.exp(fusion + heat_capacity + volume),
        warnings=tuple(warnings),
    )


__all__ = ["solid_fugacity"]
