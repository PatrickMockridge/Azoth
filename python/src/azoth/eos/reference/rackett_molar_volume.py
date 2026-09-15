"""``eos.rackett_molar_volume`` - the saturated liquid molar volume from the Rackett
equation.

```text
v = (R*Tc/Pc)*Z_RA**(1 + (1 - Tr)**(2/7)),  Z_RA = 0.29056 - 0.08775*omega
```

Spec: ``specs/calcs/eos/rackett_molar_volume.toml``.
"""

from __future__ import annotations

from azoth._registry_gen import spec as _spec_for
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import RackettMolarVolumeResult
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning
from azoth.eos.reference.pr_molar_volume import MOLAR_GAS_CONSTANT

CALC_ID = "eos.rackett_molar_volume"


def rackett_molar_volume(omega: float, Tc: Q, Pc: Q, T: Q) -> RackettMolarVolumeResult:
    """The saturated liquid molar volume of a pure component, from the Rackett equation.

    Args:
        omega: the Pitzer acentric factor. Dimensionless.
        Tc: critical temperature.
        Pc: critical pressure.
        T: absolute temperature, below the critical temperature.

    Returns:
        The saturated liquid molar volume, in ``m**3/mol``.

    Raises:
        OutOfRangeError: if ``Tc``, ``Pc`` or ``T`` is not positive.

    Example:
        >>> import azoth
        >>> q = azoth.ureg.Quantity
        >>> r = rackett_molar_volume(0.152, q(369.83, "K"), q(4_248_000.0, "Pa"), q(298.15, "K"))
        >>> round(r.v.magnitude, 18)
        8.991467403942754e-05
    """
    spec = _spec_for(CALC_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    values = {
        "omega": omega,
        "Tc": input_to_si(spec, "Tc", Tc),
        "Pc": input_to_si(spec, "Pc", Pc),
        "T": input_to_si(spec, "T", T),
    }

    apply_checks(checks.on_input, values.get, warnings)

    z_ra = 0.29056 - 0.08775 * omega
    tr = values["T"] / values["Tc"]
    exponent = 1.0 + (1.0 - tr) ** (2.0 / 7.0)
    v = MOLAR_GAS_CONSTANT * values["Tc"] / values["Pc"] * z_ra**exponent

    apply_checks(checks.derived, lambda name: v if name == "v" else None, warnings)

    return RackettMolarVolumeResult(v=from_si(v, "m**3/mol"), warnings=tuple(warnings))
