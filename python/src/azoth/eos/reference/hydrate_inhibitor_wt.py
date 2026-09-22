"""``eos.hydrate_inhibitor_wt``: the Python reference.

The mirror of ``crates/azoth-eos/src/hydrate_inhibitor_wt.rs``. NeqSim's
``HydrateInhibitorwtFlash``: a **secant on the inhibitor's moles** whose residual is the aqueous
phase's own mass fraction against the target.

**It is the one model here that reads a phase's label.** The composition comes from
``getPhase(PhaseType.AQUEOUS)``, and a cubic phase carries that label from ``PhaseEos.init``'s
three branches - reproduced by :func:`phase_label`, including the decade its first branch
carries.
"""

from __future__ import annotations

from azoth import _models_gen
from azoth.core.errors import InvalidInputError, SolverNotConvergedError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import HydrateInhibitorWtResult, Phase, PtFlashResult
from azoth.core.units import Q, input_to_si, to_si
from azoth.core.warnings import Warning
from azoth.eos.mixture import Mixture
from azoth.eos.reference._mixture_state import (
    ReducedParameters,
    mixture_parameters,
    reduced_parameters,
)
from azoth.eos.reference.pt_flash import pt_flash

MODEL_ID = "eos.hydrate_inhibitor_wt"

#: The residual the secant stops on, in mass fraction. NeqSim's own ``1e-5``.
TOLERANCE = 1.0e-5

#: The step cap and the floor, NeqSim's own ``iter < 100`` and ``|| iter < 3``.
MAXIMUM_STEPS = 100
MINIMUM_STEPS = 3

#: The bar NeqSim compares ``getVolume() / getB()`` against.
#:
#: **It is the ratio the name says, with no unit slip in it.** Both halves are in NeqSim's own
#: ``1e-5 m3/mol``: methane at 300 K and 50 bar prints ``46.0917080848140`` and
#: ``2.98485123264092``, whose quotient is ``15.4418778332323``. That is ``Z / b_r``, where
#: ``b_r`` is the reduced covolume the reduced parameters carry.
GAS_VOLUME_OVER_B = 1.75

#: NeqSim's three phase labels.
GAS = "gas"
OIL = "oil"
AQUEOUS = "aqueous"


def phase_label(
    mixture: Mixture, reduced: ReducedParameters, composition: list[float], z_factor: float
) -> str:
    """NeqSim's ``PhaseEos.init`` rule, on one phase.

    Raises:
        InvalidInputError: if a component carries no ``COMPTYPE``, which a card's substance
            does not - and a label that silently assumed one would decide between an oil and an
            aqueous phase on nothing.
    """
    _a_mix, b_mix = mixture_parameters(reduced.a, reduced.b, mixture.kij, composition)
    if b_mix <= 0.0:
        raise InvalidInputError(
            "composition", f"a phase with no covolume ({b_mix}) has no volume ratio to compare"
        )
    # `Z / b_r` is `V_m/b`: the reduced covolume already carries `P/(R T)`, so no unit enters.
    if z_factor / b_mix > GAS_VOLUME_OVER_B:
        return GAS

    hydrocarbons = 0.0
    aqueous = 0.0
    for component, fraction in zip(mixture.components, composition, strict=True):
        if not component.component_class:
            raise InvalidInputError(
                "components",
                "a component carries no `COMPTYPE`, and a phase's label is decided by whether "
                "its hydrocarbons outweigh its aqueous components",
            )
        # "Hydrocarbon" is wider than the name: `inert` is counted on the same side, which
        # is what puts a pure CO2 liquid on the OIL branch.
        if component.component_class in ("hc", "inert") and fraction > 0.0:
            hydrocarbons += fraction
        else:
            aqueous += fraction
    return OIL if hydrocarbons > aqueous else AQUEOUS


def _phases_of(flash: PtFlashResult) -> list[tuple[float, list[float], float]]:
    """The phases a two-phase flash reported, as ``(amount, composition, z factor)``.

    **Vapour first, which is NeqSim's own order** and the order the probe prints.
    """
    if flash.phase is Phase.TWO_PHASE:
        beta = flash.beta or 0.0
        return [
            (beta, list(flash.y), flash.z_vapour),
            (1.0 - beta, list(flash.x), flash.z_liquid),
        ]
    if flash.phase is Phase.ALL_LIQUID:
        return [(1.0, list(flash.x), flash.z_liquid)]
    return [(1.0, list(flash.y), flash.z_vapour)]


def _weight_fraction(
    mixture: Mixture, composition: list[float], inhibitor: int, water: int
) -> float:
    """The inhibitor's mass fraction of the inhibitor-and-water pair of one phase."""
    components = mixture.components

    def grams(index: int) -> float:
        mass = components[index].molar_mass
        if mass is None:
            raise InvalidInputError(
                "components",
                "a dosing fraction is a mass fraction of the inhibitor and the water, and one "
                "of them carries no molar mass",
            )
        return composition[index] * float(mass.magnitude)

    inhibitor_mass = grams(inhibitor)
    water_mass = grams(water)
    total = inhibitor_mass + water_mass
    if total <= 0.0:
        raise InvalidInputError(
            "components", "this phase holds no inhibitor and no water, so it has no fraction"
        )
    return inhibitor_mass / total


def hydrate_inhibitor_wt(
    components: list[str],
    moles: list[Q],
    inhibitor: str,
    wt_target: float,
    T: Q,
    P: Q,
    eos: str = "srk",
) -> HydrateInhibitorWtResult:
    """The moles of inhibitor that put the aqueous phase at a target mass fraction.

    Raises:
        InvalidInputError: if ``moles`` does not match ``components``, if ``inhibitor`` is not
            one of them, if the feed has no water, or if ``P`` is not positive.
        SolverNotConvergedError: if no trial produces an aqueous phase, or the secant reaches
            its step cap without landing.
    """
    from azoth.eos import components as databank

    spec = _models_gen.model(MODEL_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []
    t_si = input_to_si(spec, "T", T)
    p_si = input_to_si(spec, "P", P)
    apply_checks(checks.on_input, {"T": t_si, "P": p_si, "wt_target": wt_target}.get, warnings)

    if len(moles) != len(components):
        raise InvalidInputError(
            "moles",
            f"{len(components)} component(s) and {len(moles)} mole number(s), and the secant "
            f"walks them side by side",
        )
    mixture, _ = databank.mixture_of(components, eos=eos)
    lowered = [name.strip().lower() for name in components]
    wanted = inhibitor.strip().lower()
    if wanted not in lowered:
        raise InvalidInputError(
            "inhibitor",
            f"{inhibitor!r} is not one of the feed's components, so there is nothing to add to it",
        )
    index = lowered.index(wanted)
    if "water" not in lowered:
        raise InvalidInputError(
            "components",
            "no water: this reads an aqueous phase's own inhibitor fraction, so a feed without "
            "one has nothing to hold",
        )
    water = lowered.index("water")

    from azoth.core.units import quantity

    walk = [to_si(value, "mol", "moles") for value in moles]
    error = 1.0
    old_error = 1.0
    old_c = walk[index]
    iterations = 0

    while True:
        iterations += 1
        c = walk[index]
        # NeqSim's own secant, degenerate on its first step for the same reason its sibling's is.
        derivative = (error - old_error) / (c - old_c) if c != old_c else float("nan")
        old_error = error
        old_c = c

        step = error * 0.01 if iterations < MINIMUM_STEPS + 1 else -_ratio(error, derivative) * 0.5
        walk[index] += step

        total = sum(walk)
        fractions = [value / total for value in walk]
        flash = pt_flash(mixture, quantity(t_si, "K"), quantity(p_si, "Pa"), fractions)
        reduced = reduced_parameters(mixture, t_si, p_si)

        found = None
        for _amount, composition, z_factor in _phases_of(flash):
            if phase_label(mixture, reduced, composition, z_factor) == AQUEOUS:
                found = _weight_fraction(mixture, composition, index, water)
                break
        if found is None:
            raise SolverNotConvergedError(iterations, error, TOLERANCE)
        reached = found
        error = -(reached - wt_target)

        if not (
            (abs(error) > TOLERANCE and iterations < MAXIMUM_STEPS) or iterations < MINIMUM_STEPS
        ):
            break

    if abs(error) > TOLERANCE:
        raise SolverNotConvergedError(iterations, error, TOLERANCE)

    return HydrateInhibitorWtResult(
        inhibitor_moles=walk[index],
        weight_fraction=reached,
        phases=len(_phases_of(flash)),
        iterations=iterations,
        residual=error,
        warnings=tuple(warnings),
    )


def _ratio(numerator: float, denominator: float) -> float:
    """``numerator / denominator``, giving ``NaN`` on a zero denominator as Rust and Java do."""
    if denominator == 0.0:
        return float("nan")
    return numerator / denominator
