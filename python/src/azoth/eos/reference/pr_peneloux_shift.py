"""``eos.pr_peneloux_shift`` - the Peng-Robinson Peneloux volume-translation parameter.

```text
c = 0.50033*(0.25969 - Z_RA)*R*Tc/Pc,  Z_RA = 0.29056 - 0.08775*omega
```

Spec: ``specs/calcs/eos/pr_peneloux_shift.toml``, which records the dead NeqSim
expression this is *not* and why the Rackett compressibility is the correlation
rather than a databank value.
"""

from __future__ import annotations

from azoth._registry_gen import spec as _spec_for
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import PrPenelouxShiftResult
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning
from azoth.eos.reference.pr_molar_volume import MOLAR_GAS_CONSTANT

CALC_ID = "eos.pr_peneloux_shift"


def pr_peneloux_shift(omega: float, Tc: Q, Pc: Q) -> PrPenelouxShiftResult:
    """The Peng-Robinson Peneloux volume-translation parameter for a pure component.

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
        >>> r = pr_peneloux_shift(0.152, q(369.83, "K"), q(4_248_000.0, "Pa"))
        >>> round(r.c.magnitude, 20)
        -6.349504285100194e-06
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
    c = 0.50033 * (0.25969 - z_ra) * MOLAR_GAS_CONSTANT * values["Tc"] / values["Pc"]

    apply_checks(checks.derived, lambda name: c if name == "c" else None, warnings)

    return PrPenelouxShiftResult(c=from_si(c, "m**3/mol"), warnings=tuple(warnings))
