"""``thermal.conduction_plane_wall`` - steady conduction through a plane wall.

```text
q = k * A * dT / L
```

Fourier, J. (1822). "Théorie analytique de la chaleur." Paris: Firmin Didot.

Spec: ``specs/calcs/thermal/conduction_plane_wall.yaml``

# The sign convention

``dT`` may be negative and ``q`` follows its sign. This calc models a temperature
difference across a slab, not a named hot face and cold face, so it cannot say
which side is which - and inventing a convention the inputs do not carry would be
worse than returning the signed answer and letting the caller interpret it. The
spec records the alternative and why it was not taken.

# What this calc has no bounds for

Every hydraulics calc carries at least one warning-severity bound, marking a band
where an empirical correlation stops being covered by its accuracy claim.
Fourier's law is exact within its assumptions, so there is no such band, and the
assumptions are listed as assumptions rather than dressed up as checked bounds.
"""

from __future__ import annotations

from azoth._registry_gen import spec as _spec_for
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import ConductionPlaneWallResult
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning

CALC_ID = "thermal.conduction_plane_wall"


def conduction_plane_wall(k: Q, A: Q, dT: Q, L: Q) -> ConductionPlaneWallResult:
    """Steady heat flow through a plane wall.

    Args:
        k: thermal conductivity, taken as constant across the temperature range.
        A: area of the wall face the heat flows through.
        dT: temperature difference across the wall. A *difference*, not an absolute
            temperature: pass ``Q(30, "delta_degC")``, ``Q(30, "delta_degF")`` or
            ``Q(30, "K")``. An absolute ``Q(30, "degC")`` is **refused** with
            `UnitMismatchError`, because the spec marks this input ``interval: true``
            - converting it would silently give 303.15 K where 30 was meant.
        L: wall thickness in the direction of heat flow.

    Raises:
        OutOfRangeError: if ``k``, ``A`` or ``L`` is not positive. ``dT`` is
            deliberately unbounded; a negative value reverses the flow.

    Example:
        >>> import azoth
        >>> q = azoth.ureg.Quantity
        >>> r = conduction_plane_wall(
        ...     q(45.0, "W/(m*K)"), q(2.0, "m**2"), q(30.0, "K"), q(0.05, "m")
        ... )
        >>> round(r.q.magnitude, 1)
        54000.0
    """
    spec = _spec_for(CALC_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    values = {
        "k": input_to_si(spec, "k", k),
        "A": input_to_si(spec, "A", A),
        # Reads the spec's `interval: true`, which is what makes an absolute
        # `Q(30, "degC")` a refusal rather than a silent 303.15 K.
        "dT": input_to_si(spec, "dT", dT),
        "L": input_to_si(spec, "L", L),
    }

    apply_checks(checks.on_input, values.get, warnings)

    # Guarded by the checks above, so L is non-zero here.
    q = values["k"] * values["A"] * values["dT"] / values["L"]

    apply_checks(checks.derived, lambda name: q if name == "q" else None, warnings)

    # `q` is an SI base magnitude, so it is rebuilt with `from_si` - the direction
    # that means "this number is SI base", not "this number is already in watts".
    # The two coincide for watts, which is exactly why the distinction is easy to
    # lose; see azoth.core.units.
    return ConductionPlaneWallResult(q=from_si(q, "W"), warnings=tuple(warnings))
