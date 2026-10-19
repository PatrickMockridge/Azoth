"""The arithmetic every activity-coefficient phase shares.

NeqSim's ``ComponentGE.fugcoef`` sets ``phi_i = gamma_i P0_i / P`` for a component whose
``REFERENCESTATETYPE`` is ``solvent``, and *every* GE phase inherits that method - none of
``PhaseGENRTL``, ``PhaseGEUnifac``, ``PhaseGEUniquac``, ``PhaseGEWilson`` or
``PhaseGEVanLaarAcid`` overrides ``fugcoef``. What differs between the phases is which
activity model supplies ``gamma_i`` and which correlation supplies ``P0_i``; the
composition is one line and is here rather than written once per phase.

What is here is the reference-state branch: ``ReferenceState`` carries each component's
arm and :func:`ge_fugacities_with_references` composes the two. A phase whose resolver has
not moved to it yet still refuses a solute before it is evaluated.
"""

from __future__ import annotations

import math
from typing import Any, NamedTuple

from azoth.core.errors import InvalidInputError
from azoth.core.warnings import Warning
from azoth.eos.components import form_from_type
from azoth.eos.reference.antoine_vapor_pressure import antoine_vapor_pressure


class GeFugacities(NamedTuple):
    """The composition, and the reference quantities it was built from."""

    #: ``ln(gamma_i P0_i / P)`` on the solvent arm and
    #: ``ln(gamma_i / gamma_i^inf) + ln(H_i / P)`` on the Henry arm.
    ln_phi: list[float]
    #: The pure-component saturation pressure at the state's temperature, in Pa, and
    #: **zero where the component's arm reads a Henry coefficient instead**.
    p_sat: list[float]
    #: The Henry coefficient each component's arm read, in Pa, and **zero where the arm
    #: reads none**.
    henry: list[float]
    #: The infinite-dilution activity coefficient each Henry arm divided by, and **one
    #: where the arm does not divide by one**.
    gamma_inf: list[float]
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
    return GeFugacities(
        ln_phi=combine(gamma, p_sat, p_pa),
        p_sat=p_sat,
        henry=[0.0] * len(gamma),
        gamma_inf=[1.0] * len(gamma),
        warnings=warnings,
    )


#: One bar in pascals: the databank's Henry correlation is in bar, the phase's pressure
#: in pascals, and the ratio ``H / P`` is dimensionless in either.
BAR_TO_PA = 1.0e5


class ReferenceComposition(NamedTuple):
    """The composition a Henry arm's ``gamma^inf`` is evaluated at.

    The first index is always the solute, and ``x`` runs over ``indices`` only.
    """

    indices: tuple[int, ...]
    x: tuple[float, ...]


def reference_composition(k: int, x: list[float], water: int | None) -> ReferenceComposition:
    """``PhaseGE.getActivityCoefficientInfDilWater`` and ``getActivityCoefficientInfDil``.

    With water, the reference is the two-component phase ``Phase.initRefPhases`` builds:
    the solute at ``1e-10`` mol and water at ``10.0``, so ``x = [1e-11, 1 - 1e-11]``.
    Without water, the phase's own composition is cloned with the solute at
    ``x_k = 1e-10`` and the rest renormalised: ``x'_k = 1e-10 / denom`` and
    ``x'_j = x_j / denom``, where ``denom = sum_{j != k} x_j + 1e-10``.

    The no-water form drops components at zero: they contribute nothing to any of the
    activity models, and UNIFAC divides by ``x_i``.
    """
    if water is not None:
        return ReferenceComposition((k, water), (1.0e-11, 1.0 - 1.0e-11))
    others = sum(value for j, value in enumerate(x) if j != k)
    denominator = others + 1.0e-10
    indices = [k]
    values = [1.0e-10 / denominator]
    for j, value in enumerate(x):
        if j != k and value > 0.0:
            indices.append(j)
            values.append(value / denominator)
    return ReferenceComposition(tuple(indices), tuple(values))


def henry_coefficient(
    kind: str,
    name: str,
    h0: float,
    h1: float,
    h2: float,
    h3: float,
    t_k: float,
) -> float:
    """``H(T)`` for one component, in pascals, with the cap applied.

    ``henry_iapws`` fails closed to the insoluble limit outside its row's fitted window,
    which is ``IapwsHenryLaw.isUsable``'s false; the database arms cap a non-finite,
    non-positive or over-limit value and every ion.
    """
    from azoth.core.errors import PropertyUnavailableError
    from azoth.core.units import from_si
    from azoth.eos.reference._henry import (
        HenryRecord,
        Insoluble,
        coefficient,
        is_capped_value,
    )

    if kind == "henry_iapws":
        from azoth.eos.reference.iapws_henry_law import (
            TABLE,
            gas_from_name,
            iapws_henry_law,
        )

        gas = gas_from_name(name)
        if gas is None:
            raise ValueError(f"`{name}` is not an IAPWS Henry's-law species")
        row = TABLE[gas]
        if t_k < row[3] or t_k > row[4]:
            return Insoluble.HENRY_COEFFICIENT.value * BAR_TO_PA
        return iapws_henry_law(gas, from_si(t_k, "K")).henry.to_base_units().magnitude
    record = HenryRecord(h0, h1, h2, h3)
    raw = coefficient(record, t_k)
    if is_capped_value(raw, kind == "henry_database_ion"):
        return Insoluble.HENRY_COEFFICIENT.value * BAR_TO_PA
    if not record.is_fitted():
        raise PropertyUnavailableError(
            name,
            "a Henry coefficient",
            "the compiled table carries `HenryCoef1..4` as zeros for it, which is the "
            "table's marker for a substance it fits no correlation for",
        )
    return raw * BAR_TO_PA


def ge_fugacities_with_references(
    gamma: list[float],
    params: Any,
    x: list[float],
    water: int | None,
    t_k: float,
    p_pa: float,
    reference_gamma: Any,
) -> GeFugacities:
    """The fugacity coefficients of a GE phase, from each component's reference state.

    ``params`` carries the flattened reference columns (``names``, ``reference_kind``,
    the Antoine columns and ``henry_h0``-``henry_h3``); ``reference_gamma(solute,
    composition)`` evaluates the phase's own activity model on the reference composition
    and returns the solute's ``gamma^inf``. A model that divides by ``x_i`` (UNIFAC,
    UNIQUAC) must evaluate on the pair or subset, not the full list with the other
    components zeroed.

    ``p_sat`` and ``henry`` are zero where the component's arm reads none, and
    ``gamma_inf`` is one where it does not divide by one - the convention
    ``eos.pitzer_phase`` already uses for the two Henry quantities.
    """
    n = len(params.names)
    ln_phi: list[float] = []
    p_sat: list[float] = []
    henry: list[float] = []
    gamma_inf: list[float] = []
    warnings: list[Warning] = []
    for i in range(n):
        kind = params.reference_kind[i]
        if kind == "solvent":
            subset = AntoineSubset(
                antoine_type=(params.antoine_type[i],),
                antoine_coefficients=tuple(params.antoine_coefficients[i * 5 : i * 5 + 5]),
                antoine_tc=(params.antoine_tc[i],),
                antoine_pc=(params.antoine_pc[i],),
            )
            value, caveats = saturation(subset, t_k)
            warnings.extend(caveats)
            ln_phi.append(math.log(gamma[i]) + math.log(value[0] / p_pa))
            p_sat.append(value[0])
            henry.append(0.0)
            gamma_inf.append(1.0)
        else:
            coefficient = henry_coefficient(
                kind,
                params.names[i],
                params.henry_h0[i],
                params.henry_h1[i],
                params.henry_h2[i],
                params.henry_h3[i],
                t_k,
            )
            composition = reference_composition(i, x, water)
            infinite = reference_gamma(i, composition)
            ln_phi.append(math.log(gamma[i] / infinite) + math.log(coefficient / p_pa))
            p_sat.append(0.0)
            henry.append(coefficient)
            gamma_inf.append(infinite)
    return GeFugacities(
        ln_phi=ln_phi,
        p_sat=p_sat,
        henry=henry,
        gamma_inf=gamma_inf,
        warnings=warnings,
    )
