"""``characterization.lumping`` - a cut table grouped into equal-mass lumps.

Spec: ``specs/models/characterization/lumping.toml``. Oracle:
``validation/neqsim/captures/plus_fraction_probe.tsv``.

Mirrors ``crates/azoth-characterization/src/lumping.rs`` line for line.

# Equal mass, not equal cut counts

``StandardLumpingModel.generateLumpedComposition`` walks the cuts accumulating ``z*M`` and
closes a lump when the running sum reaches a target. The target is ``W/N`` to begin with and
``(W - accumulated)/(N - k - 1)`` after each lump, so every lump but the last aims at the same
mass and the last takes whatever is left.

# The total mass is the fluid's, not the table's

``W`` and the mole-fraction total are accumulated over the *system's* TBP and plus rows, which
for a fluid with one plus end is ``z_plus`` and ``z_plus * M_plus``. The cut table's own
``sum(z*M)`` carries the abundance solve's residual instead and differs by about ``1e-11``, and
on the heavy rows that moves a partition boundary - the last lump reads ``0.041863`` from the
fluid's numbers against ``0.048901`` from the table's. That is why ``molar_mass`` and
``mole_fraction`` are inputs and why they are the *plus fraction's*.

# What is not here

The class ends by removing the plus component and adding one ``addTBPfraction`` per lump, a
mutation of a ``SystemInterface``. This id stops at the table.
"""

from __future__ import annotations

from azoth import _models_gen
from azoth.core.errors import InvalidInputError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import LumpingResult
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning

MODEL_ID = "characterization.lumping"

#: ``LumpingModel.numberOfPseudocomponents``' own default.
DEFAULT_LUMPS = 7

#: The denominator floor the class adds to the lump count.
COUNT_FLOOR = 1.0e-10


def _dimensionless(spec: dict[str, object], name: str, value: float | Q) -> float:
    """A declared dimensionless input as its SI magnitude, quantity or not."""
    if isinstance(value, int | float):
        return float(value)
    return input_to_si(spec, name, value)


def lumping(
    molar_mass: Q,
    mole_fraction: float,
    cut_z: list[float],
    cut_molar_mass: list[Q],
    cut_density: list[Q],
    number_of_lumps: int | None = None,
) -> LumpingResult:
    """A cut table grouped into equal-mass lumps.

    ``molar_mass`` and ``mole_fraction`` are the **plus fraction's** own two numbers rather than
    the table's sums - the class reads them off the fluid and the difference is enough to move a
    partition boundary on the heavy rows. ``number_of_lumps`` is the class's
    ``numberOfPseudocomponents``; absent takes its own 7.

    Raises:
        InvalidInputError: if the three vectors are not one table, or if more lumps are asked for
            than there are cuts.
        OutOfRangeError: if an input is outside its declared range.
    """
    spec = _models_gen.model(MODEL_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    values = {
        "molar_mass": input_to_si(spec, "molar_mass", molar_mass),
        "mole_fraction": _dimensionless(spec, "mole_fraction", mole_fraction),
        "number_of_lumps": (
            None if number_of_lumps is None else _dimensionless(spec, "number_of_lumps", number_of_lumps)
        ),
    }
    apply_checks(checks.on_input, values.get, warnings)

    z = [float(value) for value in cut_z]
    masses = [input_to_si(spec, "cut_molar_mass", value) for value in cut_molar_mass]
    densities = [input_to_si(spec, "cut_density", value) for value in cut_density]

    count = len(z)
    if count == 0:
        raise InvalidInputError(
            "cut_z", "a cut table to group is at least one cut; there is nothing to lump here"
        )
    if len(masses) != count or len(densities) != count:
        raise InvalidInputError(
            "cut_molar_mass",
            f"the three vectors are one table and must be the same length: {count} cuts against "
            f"{len(masses)} molar masses and {len(densities)} densities",
        )
    lumps = DEFAULT_LUMPS if values["number_of_lumps"] is None else int(values["number_of_lumps"])
    if lumps > count:
        raise InvalidInputError(
            "number_of_lumps",
            f"{lumps} lumps from {count} cuts: the partition closes a lump only on a cut, so "
            f"asking for more groups than there are cuts leaves a count the loop cannot reach",
        )

    weight_total = values["mole_fraction"] * values["molar_mass"]
    mole_fraction_total = values["mole_fraction"]
    target = weight_total / (lumps + COUNT_FLOOR)

    fraction_of_heavy_end: list[float] = []
    lump_mole_fraction: list[float] = []
    lump_molar_mass: list[float] = []
    lump_density: list[float] = []
    accumulated = 0.0
    running = 0.0
    lump_z = 0.0
    lump_weight = 0.0
    denominator = 0.0
    opened = 1
    for index in range(count):
        w = z[index] * masses[index]
        running += w
        accumulated += w
        lump_z += z[index]
        lump_weight += w
        denominator += w / densities[index]

        if not ((running >= target and lumps != opened) or index == count - 1):
            continue
        remaining = lumps - opened
        # A division by zero here is the Java one: an infinity the class never reads, because the
        # loop ends on the same cut.
        target = (weight_total - accumulated) / remaining if remaining else float("inf")
        fraction_of_heavy_end.append(lump_z / mole_fraction_total)
        lump_mole_fraction.append(lump_z)
        lump_molar_mass.append(lump_weight / lump_z)
        # The gravity is the mass-weighted harmonic mean in the same scale the table uses: a pure
        # rescaling leaves the ratio unchanged, so no unit boundary is crossed here.
        lump_density.append(lump_weight / denominator)
        running = 0.0
        lump_z = 0.0
        lump_weight = 0.0
        denominator = 0.0
        opened += 1

    apply_checks(checks.derived, lambda _name: None, warnings)

    return LumpingResult(
        fraction_of_heavy_end=tuple(fraction_of_heavy_end),
        lump_mole_fraction=tuple(lump_mole_fraction),
        lump_molar_mass=tuple(from_si(value, "kg/mol") for value in lump_molar_mass),
        lump_density=tuple(from_si(value, "kg/m**3") for value in lump_density),
        warnings=tuple(warnings),
    )
