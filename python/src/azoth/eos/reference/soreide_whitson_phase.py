"""``eos.soreide_whitson_phase`` - the phase state of a Soreide-Whitson fluid.

Spec: ``specs/models/eos/soreide_whitson_phase.toml``. A *direct* model: the root is the
cubic's, chosen by ordering, so no iteration and no algorithm block.

This is the pure-Python reference: a second, independent expression of the same physics as
the Rust kernel, held to it by ``test_cross_impl.py``. What is here is the assembly - the
fluid is ``soreide_whitson_mixture_of``'s, the reduced parameters and the phase state are
``_mixture_state``'s, and the one line this model adds is that the interaction matrix is
resolved **at the phase's composition** before either is evaluated.
"""

from __future__ import annotations

from collections.abc import Sequence

from azoth import _models_gen
from azoth.core.errors import InvalidInputError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import SoreideWhitsonPhaseResult
from azoth.core.units import Q, input_to_si
from azoth.core.warnings import Warning
from azoth.eos.components import soreide_whitson_mixture_of
from azoth.eos.reference._mixture_state import phase_state, reduced_parameters

MODEL_ID = "eos.soreide_whitson_phase"

#: The roots this model has, as the databank's cubic root convention names them.
_PHASES = {"liquid": True, "vapour": False, "vapor": False}


def soreide_whitson_phase(
    components: Sequence[str],
    T: Q,
    P: Q,
    x: Sequence[float],
    salinity: Q,
    compressed_phase: str,
) -> SoreideWhitsonPhaseResult:
    """One Soreide-Whitson phase's state at a temperature, pressure and composition.

    Args:
        components: the fluid's substances, by name. Each name decides the role it takes
            in the aqueous correlation.
        T: absolute temperature.
        P: absolute pressure.
        x: the phase's mole fractions, checked rather than renormalised.
        salinity: the equivalent-NaCl molality of the brine, in mol per kg of water.
        compressed_phase: ``"liquid"`` or ``"vapour"``.

    Returns:
        The compressibility factor at the chosen root and the fugacity coefficients.

    Raises:
        InvalidInputError: if ``x`` is not a composition, ``compressed_phase`` is neither
            spelling, a name is not in the databank, an ion is named, or ``salinity`` is
            negative.
        PropertyUnavailableError: if a component carries no heat-capacity coefficients.
        OutOfRangeError: if ``T`` or ``P`` is not positive, or no volume root exists.

    See :func:`azoth.eos.soreide_whitson_phase`.
    """
    spec = _models_gen.model(MODEL_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    t_si = input_to_si(spec, "T", T)
    p_si = input_to_si(spec, "P", P)
    salinity_si = input_to_si(spec, "salinity", salinity)
    apply_checks(
        checks.on_input, {"T": t_si, "P": p_si, "salinity": salinity_si}.get, warnings
    )

    if len(components) != len(x):
        raise InvalidInputError(
            "x", f"{len(components)} components but {len(x)} fractions; the two must match"
        )
    total = sum(x)
    if abs(total - 1.0) > 1.0e-9:
        raise InvalidInputError(
            "x",
            f"the mole fractions sum to {total}, not to one. Renormalising them here "
            "would make a composition error invisible in every number downstream, so "
            "it is refused instead",
        )
    if compressed_phase not in _PHASES:
        raise InvalidInputError(
            "compressed_phase",
            f"is {compressed_phase!r}; the roots this model has are {sorted(_PHASES)}",
        )

    fluid, _ = soreide_whitson_mixture_of(list(components), salinity_si)
    reduced = reduced_parameters(fluid, t_si, p_si)
    warnings.extend(reduced.warnings)
    # `phase_kij` and not `fluid.kij`: the rule replaces the water-gas entries of a
    # water-rich phase, so the matrix is a function of the composition this state is at.
    state = phase_state(
        reduced, fluid.phase_kij(reduced, list(x)), list(x), liquid=_PHASES[compressed_phase]
    )

    return SoreideWhitsonPhaseResult(
        z_factor=state.z,
        ln_phi=tuple(state.ln_phi),
        warnings=tuple(warnings),
    )


__all__ = ["soreide_whitson_phase"]
