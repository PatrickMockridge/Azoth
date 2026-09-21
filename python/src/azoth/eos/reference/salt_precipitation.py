"""``eos.salt_precipitation`` - how much of one mineral precipitates from a brine.

```text
e  <-  the largest extent at which SR(e) = 1, up to min_i n_i / stoc_i
```

Spec: ``specs/models/eos/salt_precipitation.toml``, which carries the bracket, the refusal and
the captured state.

NeqSim's ``CalcSaltSatauration.precipitate()``, the per-mineral half of
``MultiSaltPrecipitation``'s complementarity loop: the loop holds one amount per mineral and
takes the largest ``|SR - 1|`` each step, so porting the loop is looping over this.

**The brine is an electrolyte and not a cubic mixture**, which is why this takes names:
``mixture_of`` refuses a cubic over an ion and rightly, so the coefficients come from
:func:`azoth.eos.reference.pitzer_phase` over the same composition rather than from a caller.
"""

from __future__ import annotations

import math

from azoth import _models_gen
from azoth.core.errors import InvalidInputError, SolverNotConvergedError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import SaltPrecipitationResult
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning
from azoth.eos.reference.pitzer_phase import pitzer_phase
from azoth.eos.reference.scale_saturation_ratio import scale_saturation_ratio

MODEL_ID = "eos.salt_precipitation"

#: The bisection's stopping rule on `|log10 SR|`, NeqSim's own.
LOG10_TOLERANCE = 1.0e-8

#: How far inside the maximum extent the bracket's upper end sits.
BRACKET_INSET = 1.0e-12


def salt_precipitation(
    components: list[str], salt: str, T: Q, P: Q, z: list[float]
) -> SaltPrecipitationResult:
    """The solid one mineral takes from a brine.

    Args:
        components: the brine's substances, by name, including its ions.
        salt: the mineral's ``compsalt`` name.
        T: absolute temperature.
        P: absolute pressure.
        z: the brine's mole fractions, taken as the aqueous phase's own.

    Returns:
        The solid taken in moles per mole of feed, the ratios that bound it, and the ledger.

    Raises:
        InvalidInputError: if the salt is not in the table, or the brine does not name both of
            its ions and its water.
        OutOfRangeError: if ``T`` or ``P`` is not positive.
        SolverNotConvergedError: if the ratio is not under one at the maximum extent.
    """
    from azoth.eos import components as databank

    spec = _models_gen.model(MODEL_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []
    t_si = input_to_si(spec, "T", T)
    p_si = input_to_si(spec, "P", P)
    apply_checks(checks.on_input, {"T": t_si, "P": p_si}.get, warnings)

    n = len(components)
    if len(z) != n:
        raise InvalidInputError(
            "z", f"a brine of {n} components needs {n} mole fractions and got {len(z)}"
        )

    record = databank.salt(salt)
    if record is None:
        raise InvalidInputError(
            "salt", f"`{salt}` is not a row of `compsalt`, so it has no solubility product"
        )

    def index(name: str) -> int:
        key = name.strip().lower()
        for position, candidate in enumerate(components):
            if candidate.strip().lower() == key:
                return position
        raise InvalidInputError(
            "components",
            f"`{salt}` is built from `{name}`, and this brine does not carry it. A mineral the "
            "fluid cannot form is refused rather than reported as a zero, which would be an "
            "answer about a mineral that is not there",
        )

    first, second, water = index(record.cation), index(record.anion), index("water")
    hydrogen: float | None = None
    for position, candidate in enumerate(components):
        if candidate.strip().lower() == "h3o+":
            water_mass = databank.entry("water").molar_mass
            assert water_mass is not None  # the databank carries it for every row
            hydrogen = float(z[position] / (z[water] * float(water_mass.to("kg/mol").magnitude)))
            break

    algorithm = spec["algorithm"]
    tolerance = float(algorithm["tolerance"])

    def ratio(extent: float) -> float:
        moles = list(z)
        moles[first] -= extent * record.cation_stoichiometry
        moles[second] -= extent * record.anion_stoichiometry
        total = sum(moles)
        x = [value / total for value in moles]
        # **The coefficients are the brine's own**, re-solved at the composition the extent
        # leaves: `eos.pitzer_phase` over the same names and the same state.
        coefficients = pitzer_phase(components, T, x)
        return scale_saturation_ratio(
            salt,
            x[first],
            x[second],
            x[water],
            coefficients.gamma[first],
            coefficients.gamma[second],
            x[water] * coefficients.gamma[water],
            T,
            P,
            h3o_molality=None if hydrogen is None else from_si(hydrogen, "mol/kg"),
        ).saturation_ratio

    initial = ratio(0.0)
    if initial <= 1.0:
        return SaltPrecipitationResult(
            precipitated_moles=0.0,
            initial_saturation_ratio=initial,
            final_saturation_ratio=initial,
            iterations=0,
            extent_of_maximum=0.0,
            warnings=tuple(warnings),
        )

    maximum = min(
        z[first] / record.cation_stoichiometry,
        z[second] / record.anion_stoichiometry,
    )
    if record.water_stoichiometry > 0.0:
        maximum = min(maximum, z[water] / record.water_stoichiometry)

    upper = maximum * (1.0 - BRACKET_INSET)
    at_upper = ratio(upper)
    if not at_upper < 1.0:
        raise SolverNotConvergedError(1, at_upper, LOG10_TOLERANCE)

    low, high, extent, iterations = 0.0, upper, math.nan, 0
    for step in range(1, int(algorithm["max_iterations"]) + 1):
        iterations = step
        extent = 0.5 * (low + high)
        value = ratio(extent)
        if abs(math.log10(value)) <= LOG10_TOLERANCE:
            break
        if value > 1.0:
            low = extent
        else:
            high = extent

    final = ratio(extent)
    # **How much of the bracket the answer used.** One means an ion ran out and less means the
    # ratio crossed one first.
    extent_of_maximum = extent / maximum
    _ = tolerance
    return SaltPrecipitationResult(
        precipitated_moles=extent,
        initial_saturation_ratio=initial,
        final_saturation_ratio=final,
        iterations=iterations,
        extent_of_maximum=extent_of_maximum,
        warnings=tuple(warnings),
    )


__all__ = ["salt_precipitation"]
