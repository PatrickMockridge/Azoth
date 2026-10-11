"""``characterization.assay_mass_fractions`` - an oil assay's declared fractions in mass.

Spec: ``specs/models/characterization/assay_mass_fractions.toml``. Oracle:
``validation/neqsim/captures/oil_assay_probe.tsv``.

Mirrors ``crates/azoth-characterization/src/assay_mass_fractions.rs`` line for line.

# One basis, one closure, one answer

An assay declares each cut's yield either as a mass fraction or as a liquid-volume fraction,
never as both and never as neither, and the whole assay shares one basis. The row is normalised,
and the closure - ``|sum - 1| <= 1e-3`` - is a **hard error**: a tenth of a percent of closure is
a broken assay and not a rounding artefact.

A mass basis then resolves to itself. A volume basis multiplies each normalised volume fraction
by its cut's specific gravity and renormalises.

# The cut is a record

NeqSim's ``AssayCut`` holds twelve nullable Doubles, and the one that matters here is the
density: it is optional *per cut*, which a parallel vector cannot say - a vector can be absent
per vector and not per entry. A volume basis needs every cut's, so the first that states none is
refused by position; a mass basis reads them for the bulk density alone, so a row may state some
and not others.

# What is not here

``OilAssayCharacterisation.apply`` resolves a molar mass per cut, pulls standard components out
of the databank and adds a pseudo-component for each. The row this id takes is the two fields it
reads; the API gravity, the boiling point or Watson factor, the molar mass, the sulfur and the
nitrogen are still owed, and ``apply`` is the id that will want them.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any

from azoth import _models_gen
from azoth.core.errors import InvalidInputError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import AssayMassFractionsResult
from azoth.core.units import from_si, input_records_to_si
from azoth.core.warnings import Warning

MODEL_ID = "characterization.assay_mass_fractions"

#: ``OilAssayCharacterisation.ASSAY_CLOSURE_TOLERANCE``, on the declared fractions' sum.
ASSAY_CLOSURE_TOLERANCE = 1.0e-3

BASES = ("mass", "volume")


def assay_mass_fractions(
    basis: str,
    cuts: Sequence[Mapping[str, Any]],
) -> AssayMassFractionsResult:
    """An oil assay's declared yields, resolved to a mass basis.

    ``basis`` is ``mass`` or ``volume``. ``cuts`` is one mapping per cut, keyed by the field names
    the spec declares: ``declared_fraction`` on every cut, and ``density`` where the assay states
    one, each in the unit that field's declaration fixes - a density is a quantity in kg/m3.

    Raises:
        InvalidInputError: if the assay has no cuts, if a volume basis reaches a cut that states
            no density, if the declared fractions do not close on one within ``1e-3``, or if the
            row normalises to nothing.
        OutOfRangeError: if the declared fractions' sum falls outside its declared range.
    """
    spec = _models_gen.model(MODEL_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    # **The one conversion site, and the bridge calls the same one.** A record field's unit lives
    # in the spec's `fields` block and not in `inputs`, which is why the reader here is
    # `input_records_to_si` rather than `input_to_si` - but it is one function, so the two backends
    # cannot disagree about what a field is in.
    rows: list[dict[str, Any]] = input_records_to_si(spec, "cuts", cuts)
    if not rows:
        raise InvalidInputError("cuts", "an assay of no cuts has no fractions to resolve")

    declared = [float(row["declared_fraction"]) for row in rows]
    # A field the record omits crosses as `None`, which is a cut that states no density rather
    # than a zero - and the two must not be conflated, because zero is not a density.
    densities = [None if row["density"] is None else float(row["density"]) for row in rows]

    total = sum(declared)
    # The spec's one bound is on this sum, which is not an input - so it is a derived check and
    # the closure resolves it from what was computed rather than from what was handed over.
    apply_checks(checks.derived, {"total_declared_fraction": total}.get, warnings)
    if total != total or total in (float("inf"), float("-inf")):
        raise InvalidInputError("cuts", f"the declared fractions sum to {total}, which is no assay")
    if abs(total - 1.0) > ASSAY_CLOSURE_TOLERANCE:
        raise InvalidInputError(
            "cuts",
            f"assay fractions must sum to 1.0 within {ASSAY_CLOSURE_TOLERANCE}; "
            f"supplied sum={total}",
        )

    normalised = [value / total for value in declared]
    if basis == "mass":
        mass_fraction = normalised
    else:
        # **Every cut needs a density, and the refusal names which one does not.** A volume
        # fraction has no mass to weight by without one, so a cut that states none is a hole in
        # the conversion rather than a cut that merely misses the bulk density.
        relative = []
        for index, (fraction, density) in enumerate(zip(normalised, densities, strict=True)):
            if density is None:
                raise InvalidInputError(
                    f"cuts[{index}].density",
                    "a volume-basis assay needs every cut's density to reach a mass fraction",
                )
            relative.append(fraction * density)
        total_relative = sum(relative)
        if total_relative != total_relative or total_relative <= 0.0:
            raise InvalidInputError(
                "cuts",
                "the volume-basis row has no mass to normalise: every cut's density is zero",
            )
        mass_fraction = [mass / total_relative for mass in relative]

    # The bulk density is the mass-weighted harmonic mean, and it is answerable only when **every**
    # cut carries a density - so a mass-basis assay that states one for some cuts and not others
    # loses it rather than averaging over the ones it has.
    bulk = None
    if all(density is not None for density in densities):
        reciprocal = sum(
            mass / density
            for mass, density in zip(mass_fraction, densities, strict=True)
            if density is not None
        )
        bulk = 1.0 / reciprocal
        if bulk != bulk or bulk in (float("inf"), float("-inf")) or bulk <= 0.0:
            raise InvalidInputError(
                "bulk_density", f"the bulk density of these fractions and densities is {bulk}"
            )

    return AssayMassFractionsResult(
        mass_fraction=tuple(mass_fraction),
        total_declared_fraction=total,
        bulk_density=None if bulk is None else from_si(bulk, "kg/m**3"),
        warnings=tuple(warnings),
    )
