"""``eos.furst_electrolyte_phase`` - the phase state of a Fürst electrolyte fluid.

Spec: ``specs/models/eos/furst_electrolyte_phase.toml``. A *direct* model: the root is
solved by the shared cubic-plus-pressure iteration and there is no outer loop, so no
algorithm block.

This is the pure-Python reference: a second, independent expression of the same physics as
the Rust kernel, held to it by ``test_cross_impl.py``. What is here is the assembly - the
fluid is ``furst_mixture_of``'s, the reduced parameters and the phase state are
``_mixture_state``'s, and the one line this model adds is that the fluid's mixing rule and
electrolyte term are the model's own rather than a caller's.
"""

from __future__ import annotations

from collections.abc import Sequence

from azoth import _models_gen
from azoth.core.errors import InvalidInputError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import FurstElectrolytePhaseResult
from azoth.core.units import Q, input_to_si
from azoth.core.warnings import Warning
from azoth.eos.components import furst_mixture_of
from azoth.eos.reference._mixture_state import phase_state, reduced_parameters

MODEL_ID = "eos.furst_electrolyte_phase"

#: The roots this model has, as the databank's cubic root convention names them.
_PHASES = {"liquid": True, "vapour": False, "vapor": False}


def furst_electrolyte_phase(
    components: Sequence[str],
    T: Q,
    P: Q,
    x: Sequence[float],
    compressed_phase: str,
) -> FurstElectrolytePhaseResult:
    """One Fürst electrolyte phase's state at a temperature, pressure and composition.

    Args:
        components: the substances, by name, **ions included** - the salt is a component
            here and not a scalar, which is the opposite of the Soreide-Whitson model.
        T: absolute temperature.
        P: absolute pressure.
        x: the phase's mole fractions, checked rather than renormalised.
        compressed_phase: ``"liquid"`` or ``"vapour"``.

    Returns:
        The compressibility factor at the chosen root and the fugacity coefficients.

    Raises:
        InvalidInputError: if ``x`` is not a composition, ``compressed_phase`` is neither
            spelling, or a name is not in the databank.
        PropertyUnavailableError: if a component carries no heat-capacity coefficients.
        OutOfRangeError: if ``T`` or ``P`` is not positive, or no volume root exists.

    See :func:`azoth.eos.furst_electrolyte_phase`.
    """
    spec = _models_gen.model(MODEL_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    t_si = input_to_si(spec, "T", T)
    p_si = input_to_si(spec, "P", P)
    apply_checks(checks.on_input, {"T": t_si, "P": p_si}.get, warnings)

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

    fluid, _ = furst_mixture_of(list(components))
    reduced = reduced_parameters(fluid, t_si, p_si)
    warnings.extend(reduced.warnings)
    state = phase_state(reduced, fluid.kij, list(x), liquid=_PHASES[compressed_phase])

    return FurstElectrolytePhaseResult(
        z_factor=state.z,
        ln_phi=tuple(state.ln_phi),
        warnings=tuple(warnings),
    )


__all__ = ["furst_electrolyte_phase"]
