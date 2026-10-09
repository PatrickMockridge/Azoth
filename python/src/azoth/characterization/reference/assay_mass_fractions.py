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

# What is not here

``OilAssayCharacterisation.apply`` resolves a molar mass per cut, pulls standard components out
of the databank and adds a pseudo-component for each. ``AssayCut`` holds six optional fields per
cut, and a parallel vector cannot say "absent" per entry - only per vector - so the assay's
per-cut properties cannot be stated at this boundary at all.
"""

from __future__ import annotations

from azoth import _models_gen
from azoth.core.errors import InvalidInputError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import AssayMassFractionsResult
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning

MODEL_ID = "characterization.assay_mass_fractions"

#: ``OilAssayCharacterisation.ASSAY_CLOSURE_TOLERANCE``, on the declared fractions' sum.
ASSAY_CLOSURE_TOLERANCE = 1.0e-3

BASES = ("mass", "volume")


def assay_mass_fractions(
    basis: str,
    declared_fraction: list[float],
    density: list[Q] | None = None,
) -> AssayMassFractionsResult:
    """An oil assay's declared yields, resolved to a mass basis.

    ``basis`` is ``mass`` or ``volume``. A volume basis needs ``density``, one entry per cut, in
    the same scale the class carries in g/cm3 - crossed here as kg/m3. A mass basis needs it only
    if the bulk density is wanted.

    Raises:
        InvalidInputError: if the densities do not match the cuts, if a volume basis has none, if
            the declared fractions do not close on one within ``1e-3``, or if the row normalises
            to nothing.
        OutOfRangeError: if the declared fractions' sum falls outside its declared range.
    """
    spec = _models_gen.model(MODEL_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    if not declared_fraction:
        raise InvalidInputError(
            "declared_fraction", "an assay of no cuts has no fractions to resolve"
        )
    densities = (
        None if density is None else [input_to_si(spec, "density", value) for value in density]
    )
    if densities is not None and len(densities) != len(declared_fraction):
        raise InvalidInputError(
            "density",
            f"one density per cut is what the volume conversion needs: "
            f"{len(declared_fraction)} cuts against {len(densities)} densities",
        )
    if basis == "volume" and densities is None:
        raise InvalidInputError(
            "density",
            "a volume-basis assay needs each cut's density to reach a mass fraction",
        )

    total = sum(float(value) for value in declared_fraction)
    # The spec's one bound is on this sum, which is not an input - so it is a derived check and
    # the closure resolves it from what was computed rather than from what was handed over.
    apply_checks(checks.derived, {"total_declared_fraction": total}.get, warnings)
    if total != total or total in (float("inf"), float("-inf")):
        raise InvalidInputError(
            "declared_fraction", f"the declared fractions sum to {total}, which is no assay"
        )
    if abs(total - 1.0) > ASSAY_CLOSURE_TOLERANCE:
        raise InvalidInputError(
            "declared_fraction",
            f"assay fractions must sum to 1.0 within {ASSAY_CLOSURE_TOLERANCE}; "
            f"supplied sum={total}",
        )

    normalised = [float(value) / total for value in declared_fraction]
    if basis == "mass":
        mass_fraction = normalised
    else:
        # The densities are present: the guard above refused their absence.
        assert densities is not None  # noqa: S101 - a narrowing the guard above has established
        relative = [fraction * rho for fraction, rho in zip(normalised, densities, strict=True)]
        total_relative = sum(relative)
        if total_relative != total_relative or total_relative <= 0.0:
            raise InvalidInputError(
                "density",
                "the volume-basis row has no mass to normalise: every cut's density is zero",
            )
        mass_fraction = [mass / total_relative for mass in relative]

    # The bulk density is the mass-weighted harmonic mean, and it is only answerable when every
    # cut carries a density - which a mass-basis assay need not state.
    bulk = None
    if densities is not None:
        reciprocal = sum(
            mass / rho for mass, rho in zip(mass_fraction, densities, strict=True)
        )
        bulk = 1.0 / reciprocal
        if bulk != bulk or bulk in (float("inf"), float("-inf")) or bulk <= 0.0:
            raise InvalidInputError(
                "density", f"the bulk density of these fractions and densities is {bulk}"
            )

    return AssayMassFractionsResult(
        mass_fraction=tuple(mass_fraction),
        total_declared_fraction=total,
        bulk_density=None if bulk is None else from_si(bulk, "kg/m**3"),
        warnings=tuple(warnings),
    )
