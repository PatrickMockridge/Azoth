"""``eos.srk_peneloux_shift`` - the Soave-Redlich-Kwong Peneloux volume-translation
parameter.

```text
c = 0.40768*(0.29441 - Z_RA)*R*Tc/Pc,  Z_RA = 0.29056 - 0.08775*omega
```

Spec: ``specs/calcs/eos/srk_peneloux_shift.toml``.
"""

from __future__ import annotations

from azoth._registry_gen import spec as _spec_for
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import SrkPenelouxShiftResult
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning
from azoth.eos.reference.pr_molar_volume import MOLAR_GAS_CONSTANT

CALC_ID = "eos.srk_peneloux_shift"


def srk_peneloux_shift(omega: float, Tc: Q, Pc: Q) -> SrkPenelouxShiftResult:
    """The Soave-Redlich-Kwong Peneloux volume-translation parameter for a pure component.

    The shift is subtracted from the untranslated molar volume: ``v_corr = v - c``.

    Args:
        omega: the Pitzer acentric factor. Dimensionless.
        Tc: critical temperature.
        Pc: critical pressure.

    Returns:
        The volume-translation parameter, in ``m**3/mol``.

    Raises:
        OutOfRangeError: if ``Tc`` or ``Pc`` is not positive.

    Example:
        >>> import azoth
        >>> q = azoth.ureg.Quantity
        >>> r = srk_peneloux_shift(0.152, q(369.83, "K"), q(4_248_000.0, "Pa"))
        >>> round(r.c.magnitude, 20)
        5.072202290436595e-06
    """
    spec = _spec_for(CALC_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    values = {
        "omega": omega,
        "Tc": input_to_si(spec, "Tc", Tc),
        "Pc": input_to_si(spec, "Pc", Pc),
    }

    apply_checks(checks.on_input, values.get, warnings)

    z_ra = 0.29056 - 0.08775 * omega
    c = 0.40768 * (0.29441 - z_ra) * MOLAR_GAS_CONSTANT * values["Tc"] / values["Pc"]

    apply_checks(checks.derived, lambda name: c if name == "c" else None, warnings)

    return SrkPenelouxShiftResult(c=from_si(c, "m**3/mol"), warnings=tuple(warnings))
