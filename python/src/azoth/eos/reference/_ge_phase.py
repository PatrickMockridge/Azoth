"""The arithmetic every activity-coefficient phase shares.

NeqSim's ``ComponentGE.fugcoef`` sets ``phi_i = gamma_i P0_i / P`` for a component whose
``REFERENCESTATETYPE`` is ``solvent``, and *every* GE phase inherits that method - none of
``PhaseGENRTL``, ``PhaseGEUnifac``, ``PhaseGEUniquac``, ``PhaseGEWilson`` or
``PhaseGEVanLaarAcid`` overrides ``fugcoef``. What differs between the phases is which
activity model supplies ``gamma_i`` and which correlation supplies ``P0_i``; the
composition is one line and is here rather than written once per phase.

What is *not* here is the reference-state branch. A component tagged otherwise takes a
Henry's-law coefficient in NeqSim, which this library does not implement, so the resolvers
that build the parameters refuse such a component before a phase is evaluated.
"""

from __future__ import annotations

from typing import Any, NamedTuple

from azoth.core.errors import InvalidInputError
from azoth.core.warnings import Warning
from azoth.eos.components import form_from_type
from azoth.eos.reference.antoine_vapor_pressure import antoine_vapor_pressure


class GeFugacities(NamedTuple):
    """``ln_phi`` and the saturation pressures it was built from, in Pa."""

    #: ``ln(gamma_i P0_i / P)`` per component.
    ln_phi: list[float]
    #: The pure-component saturation pressure at the state's temperature, in Pa.
    p_sat: list[float]
    #: Caveats from the correlations.
    warnings: list[Warning]


class AntoineSubset(NamedTuple):
    """A per-component Antoine table, as :func:`saturation` reads one."""

    antoine_type: tuple[str, ...]
    antoine_coefficients: tuple[float, ...]
    antoine_tc: tuple[float, ...]
    antoine_pc: tuple[float, ...]


def uncovered(params: Any) -> AntoineSubset:
    """The components of ``params`` that take the Antoine fallback.

    **A phase that overrides some components' ``P0`` from another correlation must not ask
    :func:`saturation` about those.** Upstream's acid phase takes the three acids' vapour
    pressures from
    :mod:`azoth.eos.reference.nitric_sulfuric_acid_vapor_pressure` and reaches Antoine only
    for anything else, and the acids carry upstream's ``none`` marker in the databank - so
    asking refuses the phase for its own components. The resolver agrees: it requires an
    Antoine row only of a component the acid correlation does not cover.
    """
    keep = [i for i, index in enumerate(params.acid_index) if index == 0]
    return AntoineSubset(
        antoine_type=tuple(params.antoine_type[i] for i in keep),
        antoine_coefficients=tuple(
            params.antoine_coefficients[i * 5 + k] for i in keep for k in range(5)
        ),
        antoine_tc=tuple(params.antoine_tc[i] for i in keep),
        antoine_pc=tuple(params.antoine_pc[i] for i in keep),
    )


def saturation(antoine: Any, t_k: float) -> tuple[list[float], list[Warning]]:
    """The pure-component saturation pressures, one evaluation per component.

    Split out of :func:`ge_fugacities` because not every GE phase takes its ``P0`` from
    Antoine. ``eos.ge_van_laar_acid_phase`` takes it from the Taleb correlation for the
    three acids it models and from Antoine only for the species it does not - so it shares
    the *composition* and not this, and the composition is the part worth sharing.
    """
    from azoth.core.units import from_si

    warnings: list[Warning] = []
    p_sat: list[float] = []
    for i in range(len(antoine.antoine_type)):
        # **A row with no correlation is refused rather than evaluated.** Upstream's marker
        # for an *unavailable* correlation - `83b64e5`, PR #3775, which now covers 331 of its
        # 389 rows - and an activity-coefficient phase's standard state **is** `P0`, so a
        # component without one cannot be in this phase. Sent to Wagner the row's five zeros
        # give `exp(0) * Pc = Pc`: a plausible number four orders of magnitude wrong, which is
        # the failure this library exists to make impossible. The rust kernel refuses at the
        # same place.
        #
        # **The label is the raw one and the mapping happens here**, which is what the rust
        # kernel does and what `form_from_type` is for: the parameters carry the databank's
        # own `AntoineVapPresLiqType` label, and `""` from the mapping is the refusal.
        start = i * 5
        [a, b, c, d, e] = antoine.antoine_coefficients[start : start + 5]
        form = form_from_type(antoine.antoine_type[i], e)
        if not form:
            raise InvalidInputError(
                "components",
                "the component at index "
                f"{i} has no vapour-pressure correlation: its `AntoineVapPresLiqType` is "
                "`none`, upstream's marker for unavailable data. An activity-coefficient "
                "phase's standard state is the pure liquid's saturation pressure, so a "
                "component without one cannot be in this phase",
            )
        saturated = antoine_vapor_pressure(
            a,
            b,
            c,
            d,
            e,
            form,
            from_si(antoine.antoine_tc[i], "K"),
            from_si(antoine.antoine_pc[i], "Pa"),
            from_si(t_k, "K"),
        )
        warnings.extend(saturated.warnings)
        p_sat.append(saturated.p_sat.to_base_units().magnitude)
    return p_sat, warnings


def combine(gamma: list[float], p_sat: list[float], p_pa: float) -> list[float]:
    """``ln phi_i = ln gamma_i + ln(P0_i / P)``, the composition itself.

    One line, and the reason this module exists: every activity-coefficient phase in this
    library is this expression, and the phases differ only in where ``gamma`` and ``P0``
    come from.
    """
    import math

    return [math.log(g) + math.log(p0 / p_pa) for g, p0 in zip(gamma, p_sat, strict=True)]


def ge_fugacities(gamma: list[float], antoine: Any, t_k: float, p_pa: float) -> GeFugacities:
    """``phi_i = gamma_i P0_i / P``, at a state and composition.

    ``gamma`` is the activity coefficients the phase's own model produced, in component
    order; ``antoine`` is the resolved phase record the per-component vapour-pressure
    columns are read from, in the same order. The saturation pressures come back beside
    the coefficients rather than folded into them, so a caller can check the correlation
    and the arithmetic separately: a wrong ``P0`` and a wrong ``gamma`` produce the same
    kind of wrong ``phi``, and only the parts tell them apart.
    """
    p_sat, warnings = saturation(antoine, t_k)
    return GeFugacities(ln_phi=combine(gamma, p_sat, p_pa), p_sat=p_sat, warnings=warnings)
