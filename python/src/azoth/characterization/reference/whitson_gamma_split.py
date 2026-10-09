"""``characterization.whitson_gamma_split`` - a plus fraction split by Whitson's gamma.

Spec: ``specs/models/characterization/whitson_gamma_split.toml``. Oracle:
``validation/neqsim/captures/plus_fraction_probe.tsv``.

Mirrors ``crates/azoth-characterization/src/whitson_gamma_split.rs`` line for line.

# A window in molar mass, not a carbon number

The cuts are windows that begin at ``eta`` and step by the SCN increment of 14 g/mol, with the
last widened to 10 000 to catch the tail. The two carbon numbers set **how many** windows there
are and nothing else - ``WhitsonGammaModel`` never reads the plus component's name, so the first
one stays at the outer object's own value - and the abundance is the gamma density integrated
across each window.

# Two correlations, one of which falls back to the other

``densityUOP`` is ``6.0108 * M^0.17947 * Kw^-1.18241`` on the plus fraction's Watson factor;
``densitySoreide`` is ``0.2855 + C_f*(M - 66)^0.13`` clamped to ``[0.6, 1.2]``, and calls the
Watson form when its own ``C_f`` is non-positive.

# The gamma function is an eight-term fit

``gamma`` reduces its argument to ``[0, 1)``, evaluates the classical polynomial for
``Gamma(1+z)``, and reflects below one by dividing by ``X``. Reproduced rather than replaced with
a library call, because the polynomial's error is part of every ``P0`` and ``P1``.
"""

from __future__ import annotations

import math

from azoth import _models_gen
from azoth.core.errors import InvalidInputError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import WhitsonGammaSplitResult
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning

MODEL_ID = "characterization.whitson_gamma_split"

#: The SCN increment the windows step by, g/mol.
SCN_INCREMENT = 14.0

#: The upper edge the last window is widened to, g/mol.
LAST_WINDOW_UPPER = 10_000.0

#: The floor the integrated window abundance is held above, before renormalisation.
WINDOW_FLOOR = 1.0e-15

#: The denominator below which a window's first moment is replaced by its midpoint.
MOMENT_FLOOR = 1.0e-15

#: The term the incomplete-gamma series stops at, and its own cap.
SERIES_TOLERANCE = 1.0e-8
SERIES_MAX_TERMS = 10_000

#: ``WhitsonGammaModel``'s constructor values for the two shape parameters.
DEFAULT_SHAPE = 1.0
DEFAULT_ETA_GMOL = 90.0

#: The Søreide correlation's floor, intercept and clamp.
SOREIDE_MOLAR_MASS_FLOOR = 66.0
SOREIDE_INTERCEPT = 0.2855
SOREIDE_LOWER = 0.6
SOREIDE_UPPER = 1.2

#: The eight-term fit for ``Gamma(1+z)`` on ``[0, 1)``.
GAMMA_COEFFICIENTS = (
    -0.577191652,
    0.988205891,
    -0.897056937,
    0.918206857,
    -0.756704078,
    0.482199394,
    -0.193527818,
    0.035868343,
)

DENSITY_MODELS = ("uop", "soreide")


def gamma(x: float) -> float:
    """``WhitsonGammaModel.gamma``, the fit with its reflection."""
    constant = 1.0
    reduced = x + 1.0 if x < 1.0 else x
    while reduced >= 2.0:
        constant *= reduced - 1.0
        reduced -= 1.0
    reduced -= 1.0
    polynomial = 1.0
    power = reduced
    for coefficient in GAMMA_COEFFICIENTS:
        polynomial += coefficient * power
        power *= reduced
    value = constant * polynomial
    return value / x if x < 1.0 else value


def p0_p1(molar_mass: float, eta: float, shape: float, scale: float) -> list[float]:
    """``WhitsonGammaModel.P0P1``: the incomplete-gamma pair at one molar mass, g/mol."""
    if molar_mass == eta:
        return [0.0, 0.0]
    y = (molar_mass - eta) / scale
    q = math.exp(-y) * y**shape / gamma(shape)
    term = 1.0 / shape
    total = term
    for j in range(1, SERIES_MAX_TERMS + 1):
        term *= y / (shape + j)
        total += term
        if abs(term) <= SERIES_TOLERANCE:
            return [q * total, q * (total - 1.0 / shape)]
    # The class leaves P0 and P1 at zero when the series has not broken by the cap.
    return [0.0, 0.0]


def watson_factor(m_plus: float, dens_plus: float) -> float:
    """``getWatsonKFactor``. ``dens_plus`` is a specific gravity in g/cm3."""
    return 4.5579 * (m_plus * 1000.0) ** 0.15178 * dens_plus**-1.18241


def estimate_shape(m_plus: float, dens_plus: float) -> float:
    """``WhitsonGammaModel.estimateAlpha``, the four bands on the Watson factor."""
    kw = watson_factor(m_plus, dens_plus)
    if kw >= 12.5:
        return 0.5 + 0.1 * (kw - 12.5)
    if kw >= 11.5:
        return 1.0 + 0.5 * (kw - 11.5)
    if kw >= 10.5:
        return 1.5 + 0.5 * (kw - 10.5)
    return 2.0 + 0.5 * (10.5 - kw)


def _density_uop(molar_masses: list[float], kw: float) -> list[float]:
    """``densityUOP``: every cut's gravity from the plus fraction's Watson factor, in g/cm3."""
    return [
        6.0108 * (m * 1000.0) ** 0.17947 * kw**-1.18241 if m > 0 else 0.0 for m in molar_masses
    ]


def _density_soreide(
    molar_masses: list[float], m_plus: float, dens_plus: float, kw: float
) -> list[float]:
    """``densitySoreide``, which falls back to the Watson form on a non-positive factor."""
    exponent = max(m_plus * 1000.0 - SOREIDE_MOLAR_MASS_FLOOR, 1.0)
    factor = (dens_plus - SOREIDE_INTERCEPT) / exponent**0.13
    if factor <= 0.0 or math.isnan(factor):
        return _density_uop(molar_masses, kw)
    out = []
    for m in molar_masses:
        if m <= 0:
            out.append(0.0)
            continue
        argument = max(m * 1000.0 - SOREIDE_MOLAR_MASS_FLOOR, 1.0)
        value = SOREIDE_INTERCEPT + factor * argument**0.13
        out.append(min(max(value, SOREIDE_LOWER), SOREIDE_UPPER))
    return out


def _dimensionless(spec: dict[str, object], name: str, value: float | Q) -> float:
    """A declared dimensionless input as its SI magnitude, quantity or not."""
    if isinstance(value, int | float):
        return float(value)
    return input_to_si(spec, name, value)


def whitson_gamma_split(
    molar_mass: Q,
    density: Q,
    mole_fraction: float,
    first_carbon_number: int,
    last_carbon_number: int,
    alpha: float | None = None,
    eta: Q | None = None,
    density_model: str | None = None,
    auto_estimate_shape: bool | None = None,
) -> WhitsonGammaSplitResult:
    """A plus fraction split into cuts by Whitson's three-parameter gamma distribution.

    ``first_carbon_number`` labels the row and ``last_carbon_number`` is the only bound the
    windowing respects, because the windows step from ``eta`` by 14 g/mol of molar mass rather
    than by carbon number. ``alpha`` is the shape and ``eta`` the minimum molar mass; absent, the
    class's own ``1.0`` and ``0.090`` kg/mol.

    Raises:
        InvalidInputError: if the plus fraction's molar mass is not above ``eta``, so the derived
            gamma scale is not positive.
        OutOfRangeError: if an input is outside its declared range.
    """
    spec = _models_gen.model(MODEL_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    values = {
        "molar_mass": input_to_si(spec, "molar_mass", molar_mass),
        "density": input_to_si(spec, "density", density),
        "mole_fraction": _dimensionless(spec, "mole_fraction", mole_fraction),
        "first_carbon_number": _dimensionless(
            spec, "first_carbon_number", first_carbon_number
        ),
        "last_carbon_number": _dimensionless(spec, "last_carbon_number", last_carbon_number),
        "alpha": None if alpha is None else _dimensionless(spec, "alpha", alpha),
        "eta": None if eta is None else input_to_si(spec, "eta", eta),
    }
    apply_checks(checks.on_input, values.get, warnings)

    first = int(values["first_carbon_number"])
    last = int(values["last_carbon_number"])
    if last <= first:
        raise InvalidInputError(
            "last_carbon_number",
            f"{last} leaves no cut above {first}: the range is half-open, so at least one cut has "
            f"to lie between them",
        )

    # The class works in g/mol and g/cm3 from here on, the conversion made once.
    m_plus = values["molar_mass"] * 1000.0
    dens_plus = values["density"] / 1000.0
    given_eta = values["eta"]
    eta_gmol = DEFAULT_ETA_GMOL if given_eta is None else given_eta * 1000.0
    # The estimate runs before the scale is derived, so it moves both parameters.
    shape = (
        estimate_shape(values["molar_mass"], dens_plus)
        if auto_estimate_shape
        else (DEFAULT_SHAPE if values["alpha"] is None else values["alpha"])
    )
    scale = (m_plus - eta_gmol) / shape
    if scale <= 0.0:
        raise InvalidInputError(
            "eta",
            f"the derived gamma scale is {scale}, because the plus fraction's {m_plus} g/mol is "
            f"not above eta's {eta_gmol} g/mol",
        )

    cut_z: list[float] = []
    cut_molar_mass: list[float] = []
    window_lower = eta_gmol
    for i in range(first, last):
        window_upper = window_lower + SCN_INCREMENT
        if i == last - 1:
            window_upper = LAST_WINDOW_UPPER

        low = p0_p1(window_lower, eta_gmol, shape, scale)
        high = p0_p1(window_upper, eta_gmol, shape, scale)
        window = high[0] - low[0]
        if window < WINDOW_FLOOR:
            window = WINDOW_FLOOR
        cut_z.append(window * mole_fraction)

        # The moment's denominator is the **unfloored** difference.
        denominator = high[0] - low[0]
        if abs(denominator) < MOMENT_FLOOR:
            mean = 0.5 * (window_lower + window_upper)
        else:
            mean = eta_gmol + shape * scale * (high[1] - low[1]) / denominator
        cut_molar_mass.append(mean / 1000.0)
        window_lower = window_upper

    total = sum(cut_z)
    if total > 0.0:
        cut_z = [value * mole_fraction / total for value in cut_z]

    kw = watson_factor(values["molar_mass"], dens_plus)
    if (density_model or "uop") == "soreide":
        densities = _density_soreide(cut_molar_mass, values["molar_mass"], dens_plus, kw)
    else:
        densities = _density_uop(cut_molar_mass, kw)

    apply_checks(checks.derived, lambda _name: None, warnings)

    return WhitsonGammaSplitResult(
        cut_z=tuple(cut_z),
        cut_molar_mass=tuple(from_si(m, "kg/mol") for m in cut_molar_mass),
        cut_density=tuple(from_si(d * 1000.0, "kg/m**3") for d in densities),
        shape=shape,
        minimum_molar_mass=from_si(eta_gmol / 1000.0, "kg/mol"),
        scale=from_si(scale / 1000.0, "kg/mol"),
        warnings=tuple(warnings),
    )
