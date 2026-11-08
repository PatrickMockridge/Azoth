"""Heat transfer calculations.

The first namespace beside :mod:`azoth.hydraulics`. Its existence is the point:
everything here - the spec, the generated range checks, the test cases, the
documentation - flows from the same ``specs/calcs/`` tree and the same tooling,
with no hydraulics-shaped assumption left in the pipeline.

* :func:`conduction_plane_wall` - steady conduction through a slab

# Which implementation answers

As in the hydraulics namespace, each function dispatches to the Rust extension
when it is built and to :mod:`azoth.thermal.reference` otherwise. Both are always
reachable - see :func:`azoth.backends` and :func:`azoth.use_backend`.
"""

from __future__ import annotations

from azoth._dispatch import resolve
from azoth.core.result import ConductionPlaneWallResult, ShellAndTubeRatingResult
from azoth.core.units import Q

__all__ = [
    "conduction_plane_wall",
    "shell_and_tube_rating",
]

_CONDUCTION_PLANE_WALL = "thermal.conduction_plane_wall"
_SHELL_AND_TUBE_RATING = "thermal.shell_and_tube_rating"


def conduction_plane_wall(k: Q, A: Q, dT: Q, L: Q) -> ConductionPlaneWallResult:
    """Steady heat flow through a plane wall.

    ``dT`` is a temperature *difference*, not an absolute temperature: pass
    ``Q(30, "delta_degC")`` or ``Q(30, "K")``. ``q`` follows the sign of ``dT``.

    Raises:
        OutOfRangeError: if ``k``, ``A`` or ``L`` is not positive.

    See :func:`azoth.thermal.reference.conduction_plane_wall`.
    """
    return resolve(_CONDUCTION_PLANE_WALL)(k=k, A=A, dT=dT, L=L)  # type: ignore[no-any-return]


def shell_and_tube_rating(
    tube_outside_diameter: Q,
    tube_inside_diameter: Q,
    tube_length: Q,
    tube_count: float,
    tube_passes: float,
    tube_pitch: Q,
    triangular_pitch: bool,
    shell_inside_diameter: Q,
    baffle_spacing: Q,
    baffle_count: float,
    baffle_cut: float,
    tube_wall_conductivity: Q,
    tube_density: Q,
    tube_viscosity: Q,
    tube_cp: Q,
    tube_conductivity: Q,
    tube_mass_flow: Q,
    shell_density: Q,
    shell_viscosity: Q,
    shell_cp: Q,
    shell_conductivity: Q,
    shell_mass_flow: Q,
    fouling_tube: float,
    fouling_shell: float,
    shell_side_method: str,
    tube_fouling_layer_thickness: Q | None = None,
    tube_fouling_layer_conductivity: Q | None = None,
    shell_viscosity_wall: Q | None = None,
    has_sealing: bool | None = None,
    sealing_pairs: float | None = None,
    tube_to_baffle_clearance: Q | None = None,
    shell_to_baffle_clearance: Q | None = None,
    bypass_area: Q | None = None,
) -> ShellAndTubeRatingResult:
    """A shell-and-tube exchanger's thermal-hydraulic rating.

    ``fouling_tube`` and ``fouling_shell`` are in m**2*K/W, declared dimensionless
    because no vocabulary unit expresses the area-specific resistance.

    Raises:
        InvalidInputError: for a state the class's ``calculateStrict`` refuses.
        OutOfRangeError: for an input outside the spec's declared bounds.

    See :func:`azoth.thermal.reference.shell_and_tube_rating`.
    """
    return resolve(_SHELL_AND_TUBE_RATING)(  # type: ignore[no-any-return]
        tube_outside_diameter=tube_outside_diameter,
        tube_inside_diameter=tube_inside_diameter,
        tube_length=tube_length,
        tube_count=tube_count,
        tube_passes=tube_passes,
        tube_pitch=tube_pitch,
        triangular_pitch=triangular_pitch,
        shell_inside_diameter=shell_inside_diameter,
        baffle_spacing=baffle_spacing,
        baffle_count=baffle_count,
        baffle_cut=baffle_cut,
        tube_wall_conductivity=tube_wall_conductivity,
        tube_density=tube_density,
        tube_viscosity=tube_viscosity,
        tube_cp=tube_cp,
        tube_conductivity=tube_conductivity,
        tube_mass_flow=tube_mass_flow,
        shell_density=shell_density,
        shell_viscosity=shell_viscosity,
        shell_cp=shell_cp,
        shell_conductivity=shell_conductivity,
        shell_mass_flow=shell_mass_flow,
        fouling_tube=fouling_tube,
        fouling_shell=fouling_shell,
        shell_side_method=shell_side_method,
        tube_fouling_layer_thickness=tube_fouling_layer_thickness,
        tube_fouling_layer_conductivity=tube_fouling_layer_conductivity,
        shell_viscosity_wall=shell_viscosity_wall,
        has_sealing=has_sealing,
        sealing_pairs=sealing_pairs,
        tube_to_baffle_clearance=tube_to_baffle_clearance,
        shell_to_baffle_clearance=shell_to_baffle_clearance,
        bypass_area=bypass_area,
    )
