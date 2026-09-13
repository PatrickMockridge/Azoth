"""``hydraulics.control_valve_cv`` - liquid flow through a control valve.

```text
q = CV_TO_SI * Cv * sqrt(dP / SG)
```

Spec: ``specs/calcs/hydraulics/control_valve_cv.yaml``

# The coefficient is an input, and its units are the hard part

``Cv`` is supplied by the caller: IEC 60534's published coefficient values are a
table, and this library does not reproduce tables. That much is the same as ``Cd``
in :mod:`~azoth.hydraulics.reference.orifice_flow`.

What makes this calc different is that ``Cv`` is *not* dimensionless. It is defined
as the flow of water in US gallons per minute at a pressure drop of one pound per
square inch, so it carries the units ``gpm / sqrt(psi)``. The spec has to declare it
``dimensionless`` because the unit vocabulary cannot express a fractional power of a
non-SI unit - so undoing that declaration correctly is this calculation's whole job,
and it is why the conversion constant below is named, derived and tested rather than
folded into an input description.

# Which convention

This takes the US ``Cv``. The metric ``Kv`` is defined against ``m**3/h`` and
``bar`` and they differ by a factor of about 1.156, so a caller holding ``Kv``
should convert once and deliberately: the error is large, silent, and produces a
perfectly ordinary-looking flow.
"""

from __future__ import annotations

import math
from typing import Final

from azoth._registry_gen import spec as _spec_for
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import ControlValveCvResult
from azoth.core.units import Q, from_si, to_si
from azoth.core.warnings import Warning

CALC_ID = "hydraulics.control_valve_cv"

#: One US gallon in cubic metres. Exact, by definition.
GALLON_M3: Final[float] = 3.785411784e-3

#: One pound-force in newtons. Exact, by definition.
POUND_FORCE_N: Final[float] = 4.4482216152605

#: One inch in metres. Exact, by definition.
INCH_M: Final[float] = 0.0254

#: One pound per square inch in pascals, from the two definitions above.
PSI_PA: Final[float] = POUND_FORCE_N / INCH_M**2

#: The conversion from the US `Cv` convention to SI.
#:
#: ``(1 gpm in m**3/s) / sqrt(1 psi in Pa)``, so that
#: ``q[m**3/s] = CV_TO_SI * Cv * sqrt(dP[Pa] / SG)`` is the conventional relation
#: ``Q[gpm] = Cv * sqrt(dP[psi] / SG)`` written for SI inputs. Every factor in it
#: is a definition rather than a measurement, so the constant is arithmetic.
#:
#: Derived here from the named definitions rather than written as a literal, because
#: Python can. The Rust side cannot - `f64::sqrt` is not available in a `const`
#: context - so it carries the literal and a test that re-derives it; the two
#: implementations agreeing is what a cross-language test then checks.
CV_TO_SI: Final[float] = (GALLON_M3 / 60.0) / math.sqrt(PSI_PA)


def control_valve_cv(Cv: float, dP: Q, SG: float) -> ControlValveCvResult:
    """Volumetric flow through a control valve.

    Args:
        Cv: valve flow coefficient in the US convention. Dimensionless only in the
            sense that no unit expresses it; see the module docstring.
        dP: pressure drop across the valve. A magnitude; a negative value is
            rejected.
        SG: specific gravity of the liquid, relative to water at 15.6 C.

    Raises:
        OutOfRangeError: if ``Cv`` or ``SG`` is not positive, or if ``dP`` is
            negative.

    Example:
        >>> import azoth
        >>> q = azoth.ureg.Quantity
        >>> r = control_valve_cv(
        ...     10.0, q(6894.757293168361, "Pa"), 1.0
        ... )
        >>> round(r.q.magnitude / 6.30901964e-5, 9)   # in US gallons per minute
        10.0
    """
    spec = _spec_for(CALC_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    values = {
        # Dimensionless by declaration, not by nature: `Cv` carries gpm/sqrt(psi)
        # and `SG` is a ratio, so neither has a unit to convert.
        "Cv": Cv,
        "dP": to_si(dP, "Pa", "dP"),
        "SG": SG,
    }

    apply_checks(checks.on_input, values.get, warnings)

    # Guarded by the checks above, so SG is positive and dP is not negative.
    q = CV_TO_SI * values["Cv"] * math.sqrt(values["dP"] / values["SG"])

    apply_checks(checks.derived, lambda name: q if name == "q" else None, warnings)

    return ControlValveCvResult(q=from_si(q, "m**3/s"), warnings=tuple(warnings))
