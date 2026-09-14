"""``eos.pr_molar_volume`` - molar volume from a compressibility factor.

```text
v = z*R*T/P
```

Spec: ``specs/calcs/eos/pr_molar_volume.toml``, which carries why this is the one
dimensional calc in the namespace, why the gas constant is computed from its two
defining constants, and why the evaluation order matters.
"""

from __future__ import annotations

from azoth._registry_gen import spec as _spec_for
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import PrMolarVolumeResult
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning

CALC_ID = "eos.pr_molar_volume"

#: The Avogadro constant, ``6.02214076e23 /mol``. Exact by the 2019 SI definition.
AVOGADRO_PER_MOL = 6.02214076e23

#: The Boltzmann constant, ``1.380649e-23 J/K``. Exact by the 2019 SI definition.
BOLTZMANN_J_PER_K = 1.380649e-23

#: The molar gas constant, ``N_A * k_B``. Exact, because both factors are.
MOLAR_GAS_CONSTANT = AVOGADRO_PER_MOL * BOLTZMANN_J_PER_K


def pr_molar_volume(z: float, T: Q, P: Q) -> PrMolarVolumeResult:
    """Molar volume at a state, from its compressibility factor.

    ``z`` comes from :func:`azoth.eos.pr_z_factor`; either admissible root may be
    passed, the vapour one giving the vapour volume and the liquid one the liquid
    volume. That is why this takes a ``z`` rather than a phase - the caller already
    chose what the answer means.

    Args:
        z: the compressibility factor. Dimensionless.
        T: absolute temperature.
        P: absolute pressure.

    Returns:
        The molar volume, in ``m**3/mol``.

    Raises:
        OutOfRangeError: if ``z <= 0``, or if ``T`` or ``P`` is not positive. A
            molar volume is positive, so a non-positive input can only be a caller
            error, and ``P`` is a divisor besides.

    Example:
        >>> import azoth
        >>> q = azoth.ureg.Quantity
        >>> r = pr_molar_volume(0.7907789662973796, q(295.864, "K"), q(1_062_000.0, "Pa"))
        >>> round(r.v.magnitude, 19)
        0.0018317107825229842
    """
    spec = _spec_for(CALC_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    values = {
        "z": z,
        "T": input_to_si(spec, "T", T),
        "P": input_to_si(spec, "P", P),
    }

    apply_checks(checks.on_input, values.get, warnings)

    # Written in the equation's own order, left to right; the spec's worked example
    # records why the order is stated rather than incidental.
    v = z * MOLAR_GAS_CONSTANT * values["T"] / values["P"]

    apply_checks(checks.derived, lambda name: v if name == "v" else None, warnings)

    return PrMolarVolumeResult(v=from_si(v, "m**3/mol"), warnings=tuple(warnings))
