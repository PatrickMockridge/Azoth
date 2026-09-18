"""The Soave-Redlich-Kwong CPA phase state, the pure-Python reference.

The second, independent expression of the physics the Rust `azoth_eos::srk_cpa_phase`
computes: the cubic's attraction and covolume replaced by each component's fitted
`aCPA`/`bCPA`, mixed with the `cpakij_SRK` column, and the Wertheim association
contribution added to the residual Helmholtz energy.

**The fluid is resolved from names here, and again in Rust**, which is the
`eos.eos_cg_phase` precedent. For this model it matters more than for any other in the
library: an associating mixture mixes with `cpakij_SRK` and a classical one with `KIJPR`,
and on water/methanol those differ by a factor of two - so a comparison that passed a
resolved fluid from one side to the other would compare the two kernels while assuming
away the likeliest way for them to disagree.
"""

from __future__ import annotations

from azoth import _models_gen
from azoth.core.errors import InvalidInputError, OutOfRangeError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import SrkCpaPhaseResult
from azoth.core.units import Q, input_to_si, ureg
from azoth.core.warnings import Warning
from azoth.eos import components as databank
from azoth.eos.reference._association import R
from azoth.eos.reference._mixture_state import phase_state, reduced_parameters

MODEL_ID = "eos.srk_cpa_phase"

#: The two roots this model has, and the `liquid` flag each selects.
_PHASES: dict[str, bool] = {"liquid": True, "vapour": False}


def srk_cpa_phase(
    components: list[str],
    T: Q,
    P: Q,
    z: list[float],
    compressed_phase: str,
) -> SrkCpaPhaseResult:
    """One CPA phase's state at a temperature, pressure and composition.

    Args:
        components: the substance names, one per entry, resolved against the databank
            with the fitted association set and the `cpakij_SRK` column.
        T: absolute temperature.
        P: absolute pressure.
        z: the mole fractions, checked rather than renormalised.
        compressed_phase: ``"liquid"`` or ``"vapour"``. Not a compressibility factor: the
            association moves the root, so a caller could not supply one without solving
            this model first.

    Raises:
        OutOfRangeError: if ``T`` or ``P`` is not positive, or no volume root exists above
            the mixture's covolume.
        InvalidInputError: if ``z`` is not a composition, a name is not in the databank, or
            ``compressed_phase`` is neither spelling.

    See :func:`azoth.eos.srk_cpa_phase`.
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
        raise OutOfRangeError(
            "z",
            total,
            f"the mole fractions sum to {total}, not to one. Renormalising them here would "
            f"make a composition error invisible in every number downstream, so it is "
            f"refused instead",
        )
    if compressed_phase not in _PHASES:
        raise InvalidInputError(
            "compressed_phase",
            f"is {compressed_phase!r}; the roots this model has are {sorted(_PHASES)}",
        )

    mixture = databank.from_names(list(components), eos="srk", associating=True)
    reduced = reduced_parameters(mixture, t_si, p_si)
    warnings.extend(reduced.warnings)
    state = phase_state(reduced, mixture.kij, list(z), liquid=_PHASES[compressed_phase])

    r_t = R * t_si
    return SrkCpaPhaseResult(
        z_factor=state.z,
        ln_phi=tuple(state.ln_phi),
        h_res=ureg.Quantity(state.h_dep_rt * r_t, "J/mol"),
        s_res=ureg.Quantity(state.s_dep_r * R, "J/(mol*K)"),
        warnings=tuple(warnings),
    )


__all__ = ["srk_cpa_phase"]
