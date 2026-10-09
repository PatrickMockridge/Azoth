"""``characterization.tbp_closure`` - a cut's molar mass from its boiling point and gravity.

Spec: ``specs/models/characterization/tbp_closure.toml``. Oracle:
``validation/neqsim/captures/tbp_closure_probe.tsv``.

Mirrors ``crates/azoth-characterization/src/tbp_closure.rs`` line for line.

# The bisection is NeqSim's

Bracket ``[0.010, 0.800]`` kg/mol, a ``1e-9`` tolerance on the *bracket width*, at most 200
halvings, and the **final bracket's midpoint** returned - not the last midpoint evaluated.
"""

from __future__ import annotations

import math

from azoth import _models_gen
from azoth.core.errors import InvalidInputError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import TbpClosureResult
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning

MODEL_ID = "characterization.tbp_closure"

SEARCH_LOWER = 0.010
SEARCH_UPPER = 0.800
TOLERANCE = 1.0e-9
MAX_ITERATIONS = 200

CLOSURES = ("riazi_daubert_1980", "riazi_daubert_1987", "soreide", "tbp_model")


def soreide_boiling_point(molar_mass: float, density: float) -> float:
    """``TbpClosure.calcBoilingPointSoreide``, in K. ``molar_mass`` is kg/mol."""
    molar_mass_gmol = molar_mass * 1000.0
    rankine = 1928.3 - 1.695e5 * molar_mass_gmol**-0.03522 * density**3.266 * math.exp(
        -4.922e-3 * molar_mass_gmol - 4.7685 * density + 3.462e-3 * molar_mass_gmol * density
    )
    return rankine / 1.8


def _forward_boiling_point(
    closure: str, molar_mass: float, density: float, model: str | None
) -> float:
    """The boiling point a member's bisection is inverting, in K."""
    if closure == "tbp_model":
        if model is None:
            raise InvalidInputError(
                "model",
                "must be supplied for the `tbp_model` closure, which bisects that model's own "
                "boiling point",
            )
        # The same correlated boiling point `characterization.tbp_cut_properties` publishes, so
        # the two ids cannot disagree about the switch at 540 g/mol.
        from azoth.characterization.reference.tbp_cut_properties import _evaluate

        return float(_evaluate(model, molar_mass * 1000.0, density, None)["tb"])
    return soreide_boiling_point(molar_mass, density)


def _solve_molar_mass(
    closure: str, boiling_point: float, density: float, model: str | None
) -> float:
    """``TbpClosure.solveMolarMass``, reproducing its stopping rules and its return."""
    lower = SEARCH_LOWER
    upper = SEARCH_UPPER
    f_lower = _forward_boiling_point(closure, lower, density, model) - boiling_point
    f_upper = _forward_boiling_point(closure, upper, density, model) - boiling_point
    if f_lower * f_upper > 0.0:
        raise InvalidInputError(
            "boiling_point",
            f"{boiling_point} K is not attainable with this closure at a specific gravity of "
            f"{density}. Its range over the molar-mass bracket is {f_lower + boiling_point} to "
            f"{f_upper + boiling_point} K, for {SEARCH_LOWER * 1000.0} to "
            f"{SEARCH_UPPER * 1000.0} g/mol",
        )
    iterations = 0
    while iterations < MAX_ITERATIONS and (upper - lower) > TOLERANCE:
        molar_mass = 0.5 * (lower + upper)
        f_mid = _forward_boiling_point(closure, molar_mass, density, model) - boiling_point
        if f_mid == 0.0:
            return molar_mass
        if f_lower * f_mid < 0.0:
            upper = molar_mass
        else:
            lower = molar_mass
            f_lower = f_mid
        iterations += 1
    # The **bracket's** midpoint, not the last one evaluated.
    return 0.5 * (lower + upper)


def _riazi_daubert_1980(boiling_point: float, density: float) -> float:
    """The 1980 pair, in kg/mol. Specific gravity in g/cm3."""
    return 4.5673e-5 * (boiling_point * 1.8) ** 2.1962 * density**-1.0164 / 1000.0


def _riazi_daubert_1987(boiling_point: float, density: float) -> float:
    """The 1987 pair, in kg/mol."""
    molar_mass_gmol = (
        42.965
        * math.exp(
            2.097e-4 * boiling_point - 7.78712 * density + 2.08476e-3 * boiling_point * density
        )
        * boiling_point**1.26007
        * density**4.98308
    )
    return molar_mass_gmol / 1000.0


def tbp_closure(
    closure: str,
    boiling_point: Q,
    density: Q,
    model: str | None = None,
) -> TbpClosureResult:
    """A cut's molar mass, from its normal boiling point and its specific gravity.

    ``closure`` is one of ``riazi_daubert_1980``, ``riazi_daubert_1987``, ``soreide`` or
    ``tbp_model``; the last needs ``model``, from
    :func:`azoth.characterization.tbp_cut_properties`'s own vocabulary.

    Raises:
        InvalidInputError: if a `tbp_model` closure has no `model`, or the boiling point is not
            attainable over the search bracket.
        OutOfRangeError: if `boiling_point` or `density` is not positive.

    See :func:`azoth.characterization.reference.tbp_cut_properties` for the model vocabulary.
    """
    spec = _models_gen.model(MODEL_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    values = {
        "boiling_point": input_to_si(spec, "boiling_point", boiling_point),
        "density": input_to_si(spec, "density", density),
    }
    apply_checks(checks.on_input, values.get, warnings)

    # The boundary: the correlations take g/cm3, the same number a thousand smaller.
    specific_gravity = values["density"] / 1000.0
    bp = values["boiling_point"]
    if closure == "riazi_daubert_1980":
        molar_mass = _riazi_daubert_1980(bp, specific_gravity)
    elif closure == "riazi_daubert_1987":
        molar_mass = _riazi_daubert_1987(bp, specific_gravity)
    else:
        molar_mass = _solve_molar_mass(closure, bp, specific_gravity, model)

    apply_checks(checks.derived, lambda _name: None, warnings)

    return TbpClosureResult(
        molar_mass=from_si(molar_mass, "kg/mol"),
        warnings=tuple(warnings),
    )
