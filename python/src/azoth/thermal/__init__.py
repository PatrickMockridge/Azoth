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
from azoth.core.result import ConductionPlaneWallResult
from azoth.core.units import Q

__all__ = [
    "conduction_plane_wall",
]

_CONDUCTION_PLANE_WALL = "thermal.conduction_plane_wall"


def conduction_plane_wall(k: Q, A: Q, dT: Q, L: Q) -> ConductionPlaneWallResult:
    """Steady heat flow through a plane wall.

    ``dT`` is a temperature *difference*, not an absolute temperature: pass
    ``Q(30, "delta_degC")`` or ``Q(30, "K")``. ``q`` follows the sign of ``dT``.

    Raises:
        OutOfRangeError: if ``k``, ``A`` or ``L`` is not positive.

    See :func:`azoth.thermal.reference.conduction_plane_wall`.
    """
    return resolve(_CONDUCTION_PLANE_WALL)(k=k, A=A, dT=dT, L=L)  # type: ignore[no-any-return]
