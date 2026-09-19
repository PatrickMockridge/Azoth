"""The UMR-CPA phase state, the pure-Python reference.

The second, independent expression of the physics the Rust ``azoth_eos::umr_cpa_phase``
computes: a Peng-Robinson cubic whose attraction is mixed by the UMR universal rule over
the UNIFAC-UMR-PRU activity coefficients, with the Wertheim association contribution
added to the residual Helmholtz energy.

**There is no flash for this model here.** ``eos.pt_flash``'s Jacobian is ``d ln phi /
d n``, which :func:`azoth.eos.reference._mixture_state.phase_derivatives` refuses for
every excess-Gibbs rule because that derivative is the excess Gibbs energy's own
Hessian. The model is a phase state and the states it is checked at are the ones NeqSim's
own ``SystemUMRCPAEoS`` reports.

The state, the parameters and the mixing are :func:`azoth.eos.components.umr_cpa_mixture_of`'s;
the arithmetic is ``_cpa_phase``'s, the same kernel ``eos.srk_cpa_phase`` and
``eos.pr_cpa_phase`` use.
"""

from __future__ import annotations

from collections.abc import Sequence

from azoth import _models_gen
from azoth.core.errors import InvalidInputError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import UmrCpaPhaseResult
from azoth.core.units import Q, input_to_si, ureg
from azoth.core.warnings import Warning
from azoth.eos.components import umr_cpa_mixture_of
from azoth.eos.reference._mixture_state import phase_state, reduced_parameters

MODEL_ID = "eos.umr_cpa_phase"

#: The molar gas constant, in J/(mol*K) - NeqSim's, which the departure needs to reach
#: joules from its reduced form.
R = 8.3144621

#: The roots this model has, as the databank's cubic root convention names them.
_PHASES = {"liquid": True, "vapour": False, "vapor": False}


def umr_cpa_phase(
    components: list[str],
    T: Q,
    P: Q,
    z: Sequence[float],
    compressed_phase: str,
) -> UmrCpaPhaseResult:
    """One UMR-CPA phase's state at a temperature, pressure and composition.

    Args:
        components: the substance names, resolved against the databank with their
            UMR-CPA parameter set, their ``UMRCPA_MC1..5`` coefficients and their
            ``UNIFACcompUMRPRU`` group decomposition.
        T: absolute temperature.
        P: absolute pressure.
        z: the mole fractions, checked rather than renormalised.
        compressed_phase: ``"liquid"`` or ``"vapour"``.

    Returns:
        The compressibility factor, the fugacity coefficients as logarithms, and the
        residual enthalpy and entropy.

    Raises:
        InvalidInputError: if ``z`` is not a composition, ``compressed_phase`` is neither
            spelling, or a name is not in the databank.
        PropertyUnavailableError: if a component carries no ``UMRCPA_MC1..5`` set.
        OutOfRangeError: if ``T`` or ``P`` is not positive, or no volume root exists.

    See :func:`azoth.eos.umr_cpa_phase`.
    """
    spec = _models_gen.model(MODEL_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    t_si = input_to_si(spec, "T", T)
    p_si = input_to_si(spec, "P", P)
    apply_checks(checks.on_input, {"T": t_si, "P": p_si}.get, warnings)

    if len(components) != len(z):
        raise InvalidInputError(
            "z", f"{len(components)} components but {len(z)} fractions; the two must match"
        )
    total = sum(z)
    if abs(total - 1.0) > 1.0e-9:
        raise InvalidInputError(
            "z",
            f"the mole fractions sum to {total}, not to one. Renormalising them here "
            "would make a composition error invisible in every number downstream, so "
            "it is refused instead",
        )
    if compressed_phase not in _PHASES:
        raise InvalidInputError(
            "compressed_phase",
            f"is {compressed_phase!r}; the roots this model has are {sorted(_PHASES)}",
        )

    fluid, _ = umr_cpa_mixture_of(list(components))
    reduced = reduced_parameters(fluid, t_si, p_si)
    warnings.extend(reduced.warnings)
    state = phase_state(reduced, fluid.kij, list(z), liquid=_PHASES[compressed_phase])

    r_t = R * t_si
    return UmrCpaPhaseResult(
        z_factor=state.z,
        ln_phi=tuple(state.ln_phi),
        h_res=ureg.Quantity(state.h_dep_rt * r_t, "J/mol"),
        s_res=ureg.Quantity(state.s_dep_r * R, "J/(mol*K)"),
        warnings=tuple(warnings),
    )


__all__ = ["umr_cpa_phase"]
