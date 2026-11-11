"""Thermal calculations, over arrays.

Same arrangement as :mod:`azoth.batch.hydraulics`: plain numbers in the spec's canonical
units in, SI base magnitudes out, and a loop over the scalar reference rather than a
second implementation of it.
"""

from __future__ import annotations

from array import array
from collections.abc import Sequence
from dataclasses import dataclass

from azoth.batch._core import run, sequence
from azoth.batch._result import BatchResult
from azoth.core.warnings import Warning

__all__ = [
    "ConductionPlaneWallBatch",
    "LmtdCorrectionFactorBatch",
    "ShahCondensationBatch",
    "TwoPhasePressureDropBatch",
    "conduction_plane_wall",
    "lmtd_correction_factor",
    "shah_condensation",
    "two_phase_pressure_drop",
]

_CONDUCTION_PLANE_WALL = "thermal.conduction_plane_wall"
_LMTD_CORRECTION_FACTOR = "thermal.lmtd_correction_factor"
_SHAH_CONDENSATION = "thermal.shah_condensation"
_TWO_PHASE_PRESSURE_DROP = "thermal.two_phase_pressure_drop"


@dataclass(frozen=True, slots=True, eq=False, repr=False)
class ConductionPlaneWallBatch(BatchResult):
    """Result of a batch :func:`azoth.thermal.conduction_plane_wall`."""

    #: Heat flow rate per element, in watts. Signed, following the sign of `dT`.
    q: array[float]


def _build(
    columns: dict[str, object],
    units: dict[str, str],
    warnings: tuple[tuple[Warning, ...], ...],
) -> ConductionPlaneWallBatch:
    return ConductionPlaneWallBatch(warnings=warnings, units=units, q=columns["q"])  # type: ignore[arg-type]


def conduction_plane_wall(
    *,
    k: Sequence[float],
    A: Sequence[float],
    dT: Sequence[float],
    L: Sequence[float],
) -> ConductionPlaneWallBatch:
    """Steady conduction through a plane wall, over arrays.

    Inputs in W/(m*K), m**2, K and m. ``dT`` is a temperature **difference**: a bare
    number here means kelvin of difference and no 273.15 is added.

    The hazard this avoids one layer up is worth naming. The scalar API takes a `pint`
    quantity for `dT`, and a difference and an absolute temperature share a dimension,
    so ``Q(30, "degC")`` is the tempting mistake - it converts to 303.15 K and returns a
    plausible answer ten times too large. The spec marks that input ``interval: true``
    and the scalar boundary refuses the offset unit. A batch call cannot make the
    mistake at all, because it takes no quantities: there is no unit to be wrong about.
    Inputs are plain numbers in the spec's canonical unit.

    See :func:`azoth.thermal.conduction_plane_wall` for the calculation itself.
    """
    result: ConductionPlaneWallBatch = run(
        _CONDUCTION_PLANE_WALL,
        {
            "k": sequence(k, "k"),
            "A": sequence(A, "A"),
            "dT": sequence(dT, "dT"),
            "L": sequence(L, "L"),
        },
        _build,
    )
    return result


@dataclass(frozen=True, slots=True, eq=False, repr=False)
class LmtdCorrectionFactorBatch(BatchResult):
    """Result of a batch :func:`azoth.thermal.lmtd_correction_factor`."""

    #: The capacity ratio per element.
    r: array[float]
    #: The thermal effectiveness per element.
    p: array[float]
    #: The correction factor for the stated shell passes.
    ft: array[float]
    #: The correction factor for one shell pass.
    ft_1_shell_pass: array[float]
    #: The correction factor for two shell passes.
    ft_2_shell_passes: array[float]
    #: The first pass count reaching 0.75, or -1.
    required_shell_passes: array[float]


def _build_lmtd(
    columns: dict[str, object],
    units: dict[str, str],
    warnings: tuple[tuple[Warning, ...], ...],
) -> LmtdCorrectionFactorBatch:
    return LmtdCorrectionFactorBatch(
        warnings=warnings,
        units=units,
        r=columns["r"],  # type: ignore[arg-type]
        p=columns["p"],  # type: ignore[arg-type]
        ft=columns["ft"],  # type: ignore[arg-type]
        ft_1_shell_pass=columns["ft_1_shell_pass"],  # type: ignore[arg-type]
        ft_2_shell_passes=columns["ft_2_shell_passes"],  # type: ignore[arg-type]
        required_shell_passes=columns["required_shell_passes"],  # type: ignore[arg-type]
    )


def lmtd_correction_factor(
    *,
    hot_inlet_temperature: Sequence[float],
    hot_outlet_temperature: Sequence[float],
    cold_inlet_temperature: Sequence[float],
    cold_outlet_temperature: Sequence[float],
    shell_passes: Sequence[float] | None = None,
) -> LmtdCorrectionFactorBatch:
    """The LMTD correction factor, over arrays. Inputs in K.

    ``shell_passes`` is optional **per batch, not per element**: present for every
    element or for none. Omit it and every element is one shell pass.

    See :func:`azoth.thermal.lmtd_correction_factor` for the calculation itself.
    """
    arrays: dict[str, object] = {
        "hot_inlet_temperature": sequence(hot_inlet_temperature, "hot_inlet_temperature"),
        "hot_outlet_temperature": sequence(hot_outlet_temperature, "hot_outlet_temperature"),
        "cold_inlet_temperature": sequence(cold_inlet_temperature, "cold_inlet_temperature"),
        "cold_outlet_temperature": sequence(cold_outlet_temperature, "cold_outlet_temperature"),
    }
    if shell_passes is not None:
        arrays["shell_passes"] = sequence(shell_passes, "shell_passes")
    result: LmtdCorrectionFactorBatch = run(_LMTD_CORRECTION_FACTOR, arrays, _build_lmtd)
    return result


@dataclass(frozen=True, slots=True, eq=False, repr=False)
class ShahCondensationBatch(BatchResult):
    """Result of a batch :func:`azoth.thermal.shah_condensation`."""

    #: The all-liquid Dittus-Boelter coefficient per element, W/(m**2*K).
    liquid_only_htc: array[float]
    #: `G*D/mu_l` per element.
    liquid_only_reynolds: array[float]
    #: The local Shah coefficient per element, W/(m**2*K).
    local_htc: array[float]
    #: The Simpson's-rule average per element, W/(m**2*K).
    average_htc: array[float]
    #: The vertical-tube coefficient per element, W/(m**2*K).
    vertical_htc: array[float]
    #: Whether each element is inside the validated bands, as 1.0 or 0.0.
    in_valid_range: array[float]


def _build_shah(
    columns: dict[str, object],
    units: dict[str, str],
    warnings: tuple[tuple[Warning, ...], ...],
) -> ShahCondensationBatch:
    return ShahCondensationBatch(
        warnings=warnings,
        units=units,
        liquid_only_htc=columns["liquid_only_htc"],  # type: ignore[arg-type]
        liquid_only_reynolds=columns["liquid_only_reynolds"],  # type: ignore[arg-type]
        local_htc=columns["local_htc"],  # type: ignore[arg-type]
        average_htc=columns["average_htc"],  # type: ignore[arg-type]
        vertical_htc=columns["vertical_htc"],  # type: ignore[arg-type]
        in_valid_range=columns["in_valid_range"],  # type: ignore[arg-type]
    )


def shah_condensation(
    *,
    mass_flux: Sequence[float],
    tube_inner_diameter: Sequence[float],
    liquid_density: Sequence[float],
    liquid_viscosity: Sequence[float],
    liquid_cp: Sequence[float],
    liquid_conductivity: Sequence[float],
    vapor_density: Sequence[float],
    reduced_pressure: Sequence[float],
    vapor_quality: Sequence[float],
    quality_in: Sequence[float],
    quality_out: Sequence[float],
    intervals: Sequence[float] | None = None,
) -> ShahCondensationBatch:
    """Shah's condensation correlation, over arrays.

    Inputs in kg/(m**2*s), m, kg/m**3, Pa*s, J/(kg*K), W/(m*K), kg/m**3 and the three
    dimensionless qualities and pressure. ``intervals`` is optional per batch.

    See :func:`azoth.thermal.shah_condensation` for the calculation itself.
    """
    arrays: dict[str, object] = {
        "mass_flux": sequence(mass_flux, "mass_flux"),
        "tube_inner_diameter": sequence(tube_inner_diameter, "tube_inner_diameter"),
        "liquid_density": sequence(liquid_density, "liquid_density"),
        "liquid_viscosity": sequence(liquid_viscosity, "liquid_viscosity"),
        "liquid_cp": sequence(liquid_cp, "liquid_cp"),
        "liquid_conductivity": sequence(liquid_conductivity, "liquid_conductivity"),
        "vapor_density": sequence(vapor_density, "vapor_density"),
        "reduced_pressure": sequence(reduced_pressure, "reduced_pressure"),
        "vapor_quality": sequence(vapor_quality, "vapor_quality"),
        "quality_in": sequence(quality_in, "quality_in"),
        "quality_out": sequence(quality_out, "quality_out"),
    }
    if intervals is not None:
        arrays["intervals"] = sequence(intervals, "intervals")
    result: ShahCondensationBatch = run(_SHAH_CONDENSATION, arrays, _build_shah)
    return result


@dataclass(frozen=True, slots=True, eq=False, repr=False)
class TwoPhasePressureDropBatch(BatchResult):
    """Result of a batch :func:`azoth.thermal.two_phase_pressure_drop`."""

    #: The Friedel frictional gradient per element, Pa/m.
    friedel_gradient: array[float]
    #: The Friedel drop at the stated quality per element, Pa.
    friedel_pressure_drop: array[float]
    #: The averaged Friedel drop over the quality range per element, Pa.
    friedel_average_pressure_drop: array[float]
    #: The Muller-Steinhagen-Heck frictional gradient per element, Pa/m.
    muller_steinhagen_heck_gradient: array[float]
    #: The gravitational gradient per element, Pa/m.
    gravitational_gradient: array[float]
    #: The acceleration term per element, Pa; negative when the quality falls.
    acceleration_pressure_drop: array[float]


def _build_two_phase_pressure_drop(
    columns: dict[str, object],
    units: dict[str, str],
    warnings: tuple[tuple[Warning, ...], ...],
) -> TwoPhasePressureDropBatch:
    return TwoPhasePressureDropBatch(
        warnings=warnings,
        units=units,
        friedel_gradient=columns["friedel_gradient"],  # type: ignore[arg-type]
        friedel_pressure_drop=columns["friedel_pressure_drop"],  # type: ignore[arg-type]
        friedel_average_pressure_drop=columns["friedel_average_pressure_drop"],  # type: ignore[arg-type]
        muller_steinhagen_heck_gradient=columns["muller_steinhagen_heck_gradient"],  # type: ignore[arg-type]
        gravitational_gradient=columns["gravitational_gradient"],  # type: ignore[arg-type]
        acceleration_pressure_drop=columns["acceleration_pressure_drop"],  # type: ignore[arg-type]
    )


def two_phase_pressure_drop(
    *,
    mass_flux: Sequence[float],
    vapor_quality: Sequence[float],
    tube_inner_diameter: Sequence[float],
    tube_length: Sequence[float],
    liquid_density: Sequence[float],
    vapor_density: Sequence[float],
    liquid_viscosity: Sequence[float],
    vapor_viscosity: Sequence[float],
    surface_tension: Sequence[float],
    quality_in: Sequence[float],
    quality_out: Sequence[float],
    intervals: Sequence[float] | None = None,
) -> TwoPhasePressureDropBatch:
    """Friedel and Muller-Steinhagen-Heck two-phase pressure drop, over arrays.

    Inputs in kg/(m**2*s), m, kg/m**3, Pa*s, N/m and the three dimensionless qualities.
    ``intervals`` is optional per batch.

    See :func:`azoth.thermal.two_phase_pressure_drop` for the calculation itself.
    """
    arrays: dict[str, object] = {
        "mass_flux": sequence(mass_flux, "mass_flux"),
        "vapor_quality": sequence(vapor_quality, "vapor_quality"),
        "tube_inner_diameter": sequence(tube_inner_diameter, "tube_inner_diameter"),
        "tube_length": sequence(tube_length, "tube_length"),
        "liquid_density": sequence(liquid_density, "liquid_density"),
        "vapor_density": sequence(vapor_density, "vapor_density"),
        "liquid_viscosity": sequence(liquid_viscosity, "liquid_viscosity"),
        "vapor_viscosity": sequence(vapor_viscosity, "vapor_viscosity"),
        "surface_tension": sequence(surface_tension, "surface_tension"),
        "quality_in": sequence(quality_in, "quality_in"),
        "quality_out": sequence(quality_out, "quality_out"),
    }
    if intervals is not None:
        arrays["intervals"] = sequence(intervals, "intervals")
    result: TwoPhasePressureDropBatch = run(
        _TWO_PHASE_PRESSURE_DROP, arrays, _build_two_phase_pressure_drop
    )
    return result
