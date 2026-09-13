"""``eos.vdw1f_mix_binary`` - van der Waals one-fluid mixing, for a binary.

```text
a_mix = z1**2*a1 + 2*z1*z2*(1 - k12)*sqrt(a1*a2) + z2**2*a2
b_mix = z1*b1 + z2*b2
```

Spec: ``specs/calcs/eos/vdw1f_mix_binary.yaml``

# Why two components and not N

The quadratic form for ``a_mix`` is a double sum, so a general version takes a
composition vector and a matrix of ``kij``. The calc registry is scalar - its inputs
are named quantities with units, and there is no vector or matrix type - so two
components is the largest it can express. The general case belongs to the model
layer, and that layer's spec requires it to reduce to this calc at N = 2.

# The two terms a reader should check

Written longhand rather than as a loop, so the double sum's three distinct parts are
visible: the two pure terms, and the cross term carrying ``k12``. An implementation
that got the cross term's factor of two wrong, or that transposed ``z1`` and ``z2``
in it, would still look like a weighted average at a mid-range composition - which
is why the pure-component limits are a spec test.
"""

from __future__ import annotations

from azoth._registry_gen import spec as _spec_for
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import Vdw1fMixBinaryResult
from azoth.core.warnings import Warning

CALC_ID = "eos.vdw1f_mix_binary"


def vdw1f_mix_binary(
    z1: float, a1: float, a2: float, b1: float, b2: float, k12: float
) -> Vdw1fMixBinaryResult:
    """The van der Waals one-fluid mixture parameters for a binary.

    Args:
        z1: mole fraction of component 1. Component 2 is ``1 - z1``; there is no
            ``z2`` input, because passing both would allow a pair that does not sum
            to one.
        a1: the cubic's attraction parameter ``A`` for component 1, from
            :func:`azoth.eos.pr_alpha_ab`.
        a2: the same for component 2.
        b1: the cubic's repulsion parameter ``B`` for component 1.
        b2: the same for component 2.
        k12: the binary interaction parameter. **Supplied by the caller - this
            library ships no values for it**, because a table of fitted binary
            parameters is the databank it deliberately does not have.

    Returns:
        ``a_mix`` and ``b_mix``, the cubic's parameters for the mixture.

    Raises:
        OutOfRangeError: if ``z1`` is outside ``[0, 1]``, or if the resulting
            ``a_mix`` is negative - which happens when ``k12`` is outside ``[0, 2]``
            and the composition is unfavourable.

    Example:
        >>> m = vdw1f_mix_binary(0.6, 0.20206500174625697, 0.08448417260831159,
        ...                      0.02431127309496514, 0.025932024634629486, 0.05)
        >>> round(m.a_mix, 17)
        0.14584053499701302
        >>> round(m.b_mix, 17)
        0.02495957371083088
    """
    spec = _spec_for(CALC_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    values = {"z1": z1, "a1": a1, "a2": a2, "b1": b1, "b2": b2, "k12": k12}

    apply_checks(checks.on_input, values.get, warnings)

    z2 = 1.0 - z1
    a_mix = z1 * z1 * a1 + 2.0 * z1 * z2 * (1.0 - k12) * (a1 * a2) ** 0.5 + z2 * z2 * a2
    b_mix = z1 * b1 + z2 * b2

    apply_checks(checks.derived, lambda name: a_mix if name == "a_mix" else None, warnings)

    return Vdw1fMixBinaryResult(a_mix=a_mix, b_mix=b_mix, warnings=tuple(warnings))
