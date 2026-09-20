"""The `E_theta` integral of Pitzer's same-sign electrostatic mixing.

The Python mirror of ``crates/azoth-eos/src/pitzer_electrostatic.rs``.

``PitzerElectrostaticMixing``, and it is the one piece of the Pitzer kernel that is
**numerical rather than algebraic**: a Clenshaw recurrence over a 42-term Chebyshev
expansion, evaluating

.. code-block:: text

    E_theta(z_j, z_k, I) = z_j z_k ( J(x_jk) - J(x_jj)/2 - J(x_kk)/2 ) / (4 I)

with ``x_ab = 6 A_phi sqrt(I) z_a z_b`` and ``J`` the integral the recurrence computes.
The second returned entry is ``dE_theta/dI``, which NeqSim carries through the same loop.

**`J` has two branches.** ``x`` in ``(0, 1]`` is evaluated through ``x^0.2`` and ``x > 1``
through ``x^-0.1``, each against its own half of the coefficient table, which is what makes
the expansion converge over the whole range.

The fractional exponents are total: ``x`` is positive by construction, since the function
returns before reaching them unless ``A_phi > 0``, ``I > 0`` and the two charges are
same-sign and non-zero.
"""

from __future__ import annotations

#: The Chebyshev coefficients, in NeqSim's order: ``[0..21]`` for ``x <= 1``,
#: ``[21..42]`` for ``x > 1``.
COEFFICIENTS: tuple[float, ...] = (
    1.925154014814667,
    -0.060076477753119,
    -0.029779077456514,
    -0.007299499690937,
    0.000388260636404,
    0.000636874599598,
    0.000036583601823,
    -0.000045036975204,
    -0.000004537895710,
    0.000002937706971,
    0.000000396566462,
    -0.000000202099617,
    -0.000000025267769,
    0.000000013522610,
    0.000000001229405,
    -0.000000000821969,
    -0.000000000050847,
    0.000000000046333,
    0.000000000001943,
    -0.000000000002563,
    -0.000000000010991,
    0.628023320520852,
    0.462762985338493,
    0.150044637187895,
    -0.028796057604906,
    -0.036552745910311,
    -0.001668087945272,
    0.006519840398744,
    0.001130378079086,
    -0.000887171310131,
    -0.000242107641309,
    0.000087294451594,
    0.000034682122751,
    -0.000004583768938,
    -0.000003548684306,
    -0.000000250453880,
    0.000000216991779,
    0.000000080779570,
    0.000000004558555,
    -0.000000006944757,
    -0.000000002849257,
    0.000000000237816,
)

#: Below this the two charges are the same ion, and ``E_theta`` is zero.
CHARGE_TOLERANCE = 1.0e-12

#: Where the recurrence's ``b`` and ``d`` arrays sit in the workspace, and how long each
#: is. NeqSim's own offsets, kept because the recurrence writes two slots past both ends
#: of each array and the offsets are what keep them from colliding.
_B_OFFSET = 6
_D_OFFSET = 28
_WORKSPACE = 50


def calculate(
    charge_j: float, charge_k: float, ionic_strength: float, a_phi: float
) -> tuple[float, float]:
    """``E_theta`` and ``dE_theta/dI`` for a same-sign pair.

    NeqSim throws where this returns ``(0.0, 0.0)`` for a non-finite charge, an
    opposite-sign pair or a negative ionic strength. A caller that reached here with those
    has a topology error the dataset selection already refuses upstream.
    """
    if abs(charge_j - charge_k) < CHARGE_TOLERANCE or ionic_strength == 0.0 or a_phi == 0.0:
        return (0.0, 0.0)

    x_constant = 6.0 * a_phi * (ionic_strength**0.5)
    product = charge_j * charge_k

    workspace = [0.0] * _WORKSPACE
    _evaluate_integral(x_constant * product, workspace, 0)
    _evaluate_integral(x_constant * charge_j * charge_j, workspace, 2)
    _evaluate_integral(x_constant * charge_k * charge_k, workspace, 4)

    first = (
        product * (workspace[0] - 0.5 * workspace[2] - 0.5 * workspace[4]) / (4.0 * ionic_strength)
    )
    second = (
        product
        * (workspace[1] - 0.5 * workspace[3] - 0.5 * workspace[5])
        / (8.0 * ionic_strength * ionic_strength)
        - first / ionic_strength
    )
    return (first, second)


def _evaluate_integral(x: float, workspace: list[float], offset: int) -> None:
    """One ``J`` evaluation, filling ``workspace[offset]`` and ``[offset + 1]``."""
    if x <= 1.0:
        power = x**0.2
        transformed, derivative_scale, coefficients = 4.0 * power - 2.0, 0.4 * power, 0
    else:
        power = x**-0.1
        transformed, derivative_scale, coefficients = (
            (40.0 * power - 22.0) / 9.0,
            -2.0 * power / 9.0,
            21,
        )

    b, d = _B_OFFSET, _D_OFFSET
    workspace[b + 21] = 0.0
    workspace[d + 20] = 0.0
    workspace[d + 21] = 0.0
    workspace[b + 20] = COEFFICIENTS[coefficients + 20]
    workspace[b + 19] = (
        transformed * COEFFICIENTS[coefficients + 20] + COEFFICIENTS[coefficients + 19]
    )
    workspace[d + 19] = COEFFICIENTS[coefficients + 20]
    for i in range(18, -1, -1):
        workspace[b + i] = (
            transformed * workspace[b + i + 1]
            - workspace[b + i + 2]
            + COEFFICIENTS[coefficients + i]
        )
        workspace[d + i] = (
            workspace[b + i + 1] + transformed * workspace[d + i + 1] - workspace[d + i + 2]
        )

    workspace[offset] = x / 4.0 - 1.0 + 0.5 * (workspace[b] - workspace[b + 2])
    workspace[offset + 1] = x / 4.0 + derivative_scale * (workspace[d] - workspace[d + 2])
