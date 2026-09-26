"""``eos.effective_diffusion`` - the effective diffusion coefficients of a phase.

Spec: ``specs/models/eos/effective_diffusion.toml``

The Python twin of ``crates/azoth-eos/src/effective_diffusion.rs``, written to mirror it
line for line.

# The assembly

```text
D_eff_i = (1 - x_i) / sum_{j != i} x_j / D_ij
```

``PhysicalProperties.calcEffectiveDiffusionCoefficients`` delegates to the phase's
diffusivity model, and every model in the family writes the same eight lines.

# The input that was recorded as unreachable

``reactions.kinetics`` takes the effective vector as an input, on the ground that the
binary matrix behind it could not be read from outside its class - a probe called
``getFickDiffusionCoefficient``, or the *liquid* ``getFickBinaryDiffusionCoefficient``,
and got a diagonal array.

**That was the wrong method.** ``DiffusivityInterface.calcDiffusionCoefficients(int, int)``
*returns* the matrix and ``PhysicalProperties.diffusivityCalc`` is a public field, so no
cast is needed either, and the diagonal came from the liquid override's Maxwell-Stefan
correction with its derivative uninitialised.
``validation/neqsim/EffectiveDiffusionProbe.java`` prints the matrix, the vector the class
computes from it and the vector recomputed from it, in one run.

# The matrix is not symmetric

``D_ij`` is component ``i``'s coefficient at infinite dilution in ``j``, so the two
directions differ - ``1.432e-7`` against ``1.428e-7`` on the captured oil - and only row
``i`` decides ``D_eff_i``.
"""

from __future__ import annotations

from collections.abc import Sequence
from typing import Any

from azoth.core.errors import InvalidInputError, OutOfRangeError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import EffectiveDiffusionResult
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning


def effective_diffusion(
    binary_diffusion: Sequence[Sequence[Q]],
    x: Sequence[float],
) -> EffectiveDiffusionResult:
    """The effective diffusion coefficients of a phase.

    Args:
        binary_diffusion: the pair coefficients, **row-major and not symmetric**. ``D_ij``
            is component ``i``'s coefficient at infinite dilution in ``j``, and only row
            ``i`` is read for ``D_eff_i``.
        x: the phase's mole fractions.

    Returns:
        One effective coefficient per component, in ``x``'s order, in m**2/s.

    Raises:
        InvalidInputError: for fewer than two components, a matrix that is not square or
            too small for ``x``, or a negative mole fraction.
        OutOfRangeError: for a zero pair coefficient a sum would divide by, or a component
            every other one is absent from.

    Example:
        >>> import azoth
        >>> q = azoth.ureg.Quantity
        >>> r = effective_diffusion(
        ...     [[q(1.0e-9, "m**2/s"), q(2.0e-9, "m**2/s")],
        ...      [q(3.0e-9, "m**2/s"), q(4.0e-9, "m**2/s")]],
        ...     [0.4, 0.6],
        ... )
        >>> round(r.effective_diffusion[0].magnitude, 12)
        2.2e-09
    """
    spec = _spec()
    checks = checks_for(spec)
    warnings: list[Warning] = []

    count = len(x)
    if count < 2:
        raise InvalidInputError(
            "x",
            f"an effective coefficient needs a second component to diffuse into and there "
            f"are {count}",
        )
    if len(binary_diffusion) != count:
        raise InvalidInputError(
            "binary_diffusion",
            f"{len(binary_diffusion)} row(s) against {count} component(s)",
        )
    for row, values in enumerate(binary_diffusion):
        if len(values) != count:
            raise InvalidInputError(
                "binary_diffusion",
                f"row {row} has {len(values)} column(s) against {count} component(s)",
            )

    # The declared unit is m**2/s, so every entry is a quantity until it is converted -
    # a matrix is the one shape whose boundary is nested.
    matrix = [
        [input_to_si(spec, "binary_diffusion", value) for value in row] for row in binary_diffusion
    ]
    fractions = [float(value) for value in x]

    apply_checks(
        checks.on_input,
        # The class does not check a negative fraction and a negative one silently
        # subtracts from a sum, so the spec's bound is what catches it.
        {"x": min(fractions)}.get,
        warnings,
    )

    effective_diffusion: list[Q] = []
    for i in range(count):
        total = 0.0
        for j in range(count):
            if i == j:
                continue
            pair = matrix[i][j]
            if pair == 0.0:
                raise OutOfRangeError(
                    "binary_diffusion",
                    pair,
                    f"D[{i}][{j}] is zero and the sum divides by it; a zero pair "
                    "coefficient is not a state this can take a quotient of",
                )
            total += fractions[j] / pair
        if total == 0.0:
            raise OutOfRangeError(
                "x",
                fractions[i],
                f"every component other than {i} carries a zero mole fraction, so there is "
                "nothing for it to diffuse into",
            )
        effective_diffusion.append(from_si((1.0 - fractions[i]) / total, "m**2/s"))

    return EffectiveDiffusionResult(
        effective_diffusion=tuple(effective_diffusion),
        warnings=tuple(warnings),
    )


def _spec() -> dict[str, Any]:
    from azoth._models_gen import model

    return model("eos.effective_diffusion")
