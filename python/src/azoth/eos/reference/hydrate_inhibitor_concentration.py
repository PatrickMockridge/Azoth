"""``eos.hydrate_inhibitor_concentration``: the Python reference.

The mirror of ``crates/azoth-eos/src/hydrate_inhibitor_concentration.rs``. NeqSim's
``HydrateInhibitorConcentrationFlash``: a **secant on the inhibitor's moles** whose residual is
``T_hydrate - T_target``.

**The feed is in moles, not fractions.** The secant walks an absolute amount - NeqSim's first
three steps are ``error * 0.01`` moles and the rest are ``-error / dError/dC * 0.5`` moles - so
a normalised feed would reproduce the equation and not the path, and would stop somewhere else.

**Nothing here reads a phase type**, unlike ``HydrateInhibitorwtFlash``, which reads the
aqueous phase's own composition: this only asks for a hydrate temperature.
"""

from __future__ import annotations

from azoth import _models_gen
from azoth.core.errors import InvalidInputError, SolverNotConvergedError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import HydrateInhibitorConcentrationResult
from azoth.core.units import Q, from_si, input_to_si, to_si
from azoth.core.warnings import Warning
from azoth.eos.reference.hydrate_formation_temperature import hydrate_formation_temperature

MODEL_ID = "eos.hydrate_inhibitor_concentration"

#: The residual the secant stops on, in kelvin. NeqSim's own ``1e-3``.
TOLERANCE = 1.0e-3

#: The step cap and the floor, NeqSim's own ``iter < 100`` and ``|| iter < 3``.
MAXIMUM_STEPS = 100
MINIMUM_STEPS = 3


def hydrate_inhibitor_concentration(
    components: list[str],
    moles: list[float],
    inhibitor: str,
    T_target: Q,
    P: Q,
    eos: str = "srk",
    hydrate_model: str = "pvtsim",
) -> HydrateInhibitorConcentrationResult:
    """The moles of inhibitor that hold a hydrate temperature down to a target.

    Raises:
        InvalidInputError: if ``moles`` does not match ``components``, if ``inhibitor`` is not
            one of them, if the feed has no water, or if ``P`` is not positive.
        SolverNotConvergedError: if the secant reaches its step cap without landing.
    """
    from azoth.eos import components as databank

    spec = _models_gen.model(MODEL_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []
    t_target_si = input_to_si(spec, "T_target", T_target)
    p_si = input_to_si(spec, "P", P)
    apply_checks(checks.on_input, {"T_target": t_target_si, "P": p_si}.get, warnings)

    if len(moles) != len(components):
        raise InvalidInputError(
            "moles",
            f"{len(components)} component(s) and {len(moles)} mole number(s), and the secant "
            f"walks them side by side",
        )
    # **The names cross unresolved**, as for every hydrate model: the guest tables are keyed by
    # name. The index is resolved against the caller's own list and not the mixture's, which
    # carries no names here - the same move `tp_solid_flash` makes.
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
            "no water: an inhibitor holds a hydrate temperature down by moving water's "
            "fugacity, so a feed without it has nothing to inhibit",
        )
    water = lowered.index("water")

    # **A unit-bearing vector arrives as quantities**, one per entry, and the walk is in SI
    # moles: this is the library's first declared vector with a unit, so the conversion is
    # stated here rather than assumed.
    walk = [to_si(value, "mol", "moles") for value in moles]
    error = 1.0
    old_error = 1.0
    old_c = walk[index]
    iterations = 0
    hydrate_temperature = float("nan")

    while True:
        iterations += 1
        c = walk[index]
        # NeqSim's own secant. **Its first step's denominator is zero** - both amounts are the
        # feed's - so its derivative is `0/0`, which Rust and Java make `NaN` and Python
        # raises on. The three steps that never read it are why that is harmless rather than a
        # defect, and `_ratio` is what keeps the two kernels at the same arithmetic.
        derivative = _ratio(error - old_error, c - old_c)
        old_error = error
        old_c = c

        step = error * 0.01 if iterations < MINIMUM_STEPS + 1 else -_ratio(error, derivative) * 0.5
        walk[index] += step

        total = sum(walk)
        fractions = [value / total for value in walk]
        solved = hydrate_formation_temperature(
            components, P=from_si(p_si, "Pa"), z=fractions, eos=eos, hydrate_model=hydrate_model
        )
        hydrate_temperature = float(solved.temperature.to("K").magnitude)
        error = hydrate_temperature - t_target_si

        if not (
            (abs(error) > TOLERANCE and iterations < MAXIMUM_STEPS) or iterations < MINIMUM_STEPS
        ):
            break

    if abs(error) > TOLERANCE:
        # `SolverNotConvergedError` is positional: iterations, residual, tolerance.
        raise SolverNotConvergedError(iterations, error, TOLERANCE)

    return HydrateInhibitorConcentrationResult(
        inhibitor_moles=walk[index],
        weight_fraction=_weight_fraction(mixture, walk, index, water),
        hydrate_temperature=from_si(hydrate_temperature, "K"),
        iterations=iterations,
        residual=error,
        warnings=tuple(warnings),
    )


def _ratio(numerator: float, denominator: float) -> float:
    """``numerator / denominator``, giving ``NaN`` on a zero denominator as Rust and Java do.

    ``0.0 / 0.0`` is ``NaN`` in IEEE arithmetic and both of those, and a ``ZeroDivisionError``
    here - so a secant whose first step is the degenerate one would abort the twin at a state
    the other kernel carries past.
    """
    if denominator == 0.0:
        return float("nan")
    return numerator / denominator


def _weight_fraction(mixture: object, moles: list[float], inhibitor: int, water: int) -> float:
    """The inhibitor's mass fraction of the inhibitor-and-water pair, NeqSim's own denominator."""
    components = mixture.components  # type: ignore[attr-defined]
    inhibitor_mass = moles[inhibitor] * float(components[inhibitor].molar_mass.magnitude)
    water_mass = moles[water] * float(components[water].molar_mass.magnitude)
    return inhibitor_mass / (inhibitor_mass + water_mass)
