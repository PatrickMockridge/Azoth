"""``eos.pr_mass_density`` - mass density from a molar volume.

```text
rho = M/v
```

Spec: ``specs/calcs/eos/pr_mass_density.toml``, which carries why one division earns
a calc and the ``g/mol``-versus-``kg/mol`` trap that the units layer cannot catch.
"""

from __future__ import annotations

from azoth._registry_gen import spec as _spec_for
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import PrMassDensityResult
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning

CALC_ID = "eos.pr_mass_density"


def pr_mass_density(M: Q, v: Q) -> PrMassDensityResult:
    """Mass density from a molar mass and a molar volume.

    Args:
        M: molar mass. **Supplied by the caller** - this library ships no component
            data. Note the unit: ``kg/mol``, not the ``g/mol`` tables usually quote.
        v: molar volume, from :func:`azoth.eos.pr_molar_volume`.

    Returns:
        The mass density, in ``kg/m**3``.

    Raises:
        OutOfRangeError: if ``M <= 0`` or ``v <= 0``. Both are strictly positive for
            any real substance and state, and ``v`` is a divisor.

    Example:
        >>> import azoth
        >>> q = azoth.ureg.Quantity
        >>> r = pr_mass_density(q(0.0440956, "kg/mol"), q(0.0018317107825229842, "m**3/mol"))
        >>> round(r.rho.magnitude, 15)
        24.073451125981286
    """
    spec = _spec_for(CALC_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    values = {"M": input_to_si(spec, "M", M), "v": input_to_si(spec, "v", v)}

    apply_checks(checks.on_input, values.get, warnings)

    rho = values["M"] / values["v"]

    apply_checks(checks.derived, lambda name: rho if name == "rho" else None, warnings)

    return PrMassDensityResult(rho=from_si(rho, "kg/m**3"), warnings=tuple(warnings))
