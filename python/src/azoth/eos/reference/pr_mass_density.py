"""``eos.pr_mass_density`` - mass density from a molar volume.

```text
rho = M/v
```

Spec: ``specs/calcs/eos/pr_mass_density.yaml``

# Why one division earns a calc

It is the last step of the path a caller actually wants - ``pr_z_factor`` gives a
compressibility factor, ``pr_molar_volume`` turns it into a volume, and this turns
that into the density somebody asked for. It is also the step where the units are
load-bearing in both directions, and the one that gives ``kg/mol`` its first
consumer: that unit has been in the vocabulary since the vocabulary was made
checkable, with a correct conversion path and no calculation using it. A unit whose
conversion is tested is not a defect, but it is a unit nothing has ever proved works
end to end. This is that proof.

# The unit that catches people

Molar masses are conventionally quoted in ``g/mol`` and this takes ``kg/mol``, a
factor of 1000 apart. The units layer cannot catch it - both are plausible
magnitudes of the same dimension - so the spec says so in the input's description,
in the derivation and in the assumptions.
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
        OutOfRangeError: if ``M <= 0`` or ``v <= 0``.

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
