"""The CPA phase state, cubic family aside.

``eos.srk_cpa_phase`` and ``eos.pr_cpa_phase`` are the same model with a different fitted
parameter set: the cubic's attraction and covolume replaced by each component's fitted
``aCPA``/``bCPA``, mixed with that family's ``cpakij`` column, and the Wertheim association
contribution added to the residual Helmholtz energy. Everything else - the validation, the
root choice, the departure functions - is one thing, so it lives here once and the two
models are the two lines that name the family and the result type.

**The fluid is resolved from names here, and again in Rust.** ``eos.eos_cg_phase`` sets
the precedent: the names cross the boundary unresolved and each implementation looks them
up in its own databank. That is what makes the two-kernel comparison cover the
*resolution* as well as the arithmetic, which matters more for this model than for any
other in the library - an associating mixture mixes with the ``cpa`` column and a
classical one with ``KIJPR``, and on water/methanol those differ by a factor of two.

**The family is not a mixing rule.** ``AssociationCubic::of`` reads it off the cubic, and
the two fitted sets are separate fits rather than one converted into the other - water's
``kappa_AB`` is 0.0692 for SRK against 0.046473789 for PR - so choosing is a selection and
reading the wrong family is a different fluid.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

from azoth.core.errors import InvalidInputError, OutOfRangeError
from azoth.core.range import apply_checks, checks_for
from azoth.core.units import Q, input_to_si
from azoth.core.warnings import Warning
from azoth.eos import components as databank
from azoth.eos.reference._association import R
from azoth.eos.reference._mixture_state import phase_state, reduced_parameters

#: The two roots this model has, and the ``liquid`` flag each selects.
_PHASES: dict[str, bool] = {"liquid": True, "vapour": False}


@dataclass(frozen=True)
class CpaPhaseState:
    """One CPA phase's state, before it is wrapped in the family's own result type."""

    z_factor: float
    ln_phi: tuple[float, ...]
    h_res: float
    s_res: float
    warnings: tuple[Warning, ...]


def cpa_phase(
    model_id: str,
    eos: str,
    components: list[str],
    t: Q,
    p: Q,
    z: list[float],
    compressed_phase: str,
) -> CpaPhaseState:
    """One CPA phase's state at a temperature, pressure and composition.

    Args:
        model_id: the model this is being computed for, for its spec's range checks.
        eos: the cubic whose family is read - ``"srk"`` or ``"pr"`` - and with it the
            fitted set and the interaction column.
        components: the substance names, resolved against the databank.
        t: absolute temperature.
        p: absolute pressure.
        z: the mole fractions, checked rather than renormalised.
        compressed_phase: ``"liquid"`` or ``"vapour"``.

    Raises:
        OutOfRangeError: if ``t`` or ``p`` is not positive, or no volume root exists above
            the mixture's covolume.
        InvalidInputError: if ``z`` is not a composition, a name is not in the databank, or
            ``compressed_phase`` is neither spelling.
    """
    spec = _spec(model_id)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    t_si = input_to_si(spec, "T", t)
    p_si = input_to_si(spec, "P", p)
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

    mixture = databank.from_names(list(components), eos=eos, associating=True)
    reduced = reduced_parameters(mixture, t_si, p_si)
    warnings.extend(reduced.warnings)
    state = phase_state(reduced, mixture.kij, list(z), liquid=_PHASES[compressed_phase])

    r_t = R * t_si
    return CpaPhaseState(
        z_factor=state.z,
        ln_phi=tuple(state.ln_phi),
        h_res=state.h_dep_rt * r_t,
        s_res=state.s_dep_r * R,
        warnings=tuple(warnings),
    )


def phase_state_of(
    mixture: Any,
    t_si: float,
    p_si: float,
    z: list[float],
    liquid: bool,
) -> CpaPhaseState:
    """The same arithmetic with the resolution lifted out, for a caller holding a fluid.

    This is what the transport calls when a caller has built the mixture themselves, and
    what a test uses to pin the two paths apart.
    """
    reduced = reduced_parameters(mixture, t_si, p_si)
    state = phase_state(reduced, mixture.kij, list(z), liquid=liquid)
    r_t = R * t_si
    return CpaPhaseState(
        z_factor=state.z,
        ln_phi=tuple(state.ln_phi),
        h_res=state.h_dep_rt * r_t,
        s_res=state.s_dep_r * R,
        warnings=tuple(reduced.warnings),
    )


def _spec(model_id: str) -> dict[str, Any]:
    """The model's spec, imported here so the two family modules stay two lines."""
    from azoth import _models_gen

    return _models_gen.model(model_id)


__all__ = ["CpaPhaseState", "cpa_phase", "phase_state_of"]
