"""The cubic equation of state's shape: the constants that fix one cubic.

The general cubic

```text
P = R T/(v - b) - a/((v + delta1 b)(v + delta2 b))
```

is fixed by four numbers - ``omega_a``, ``omega_b``, ``delta1``, ``delta2`` - and the
three models this library ports are three settings of them. The temperature
dependence of ``a`` is a separate axis, in :mod:`azoth.eos.alpha_term`.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

_SQRT_2 = math.sqrt(2.0)


@dataclass(frozen=True)
class Cubic:
    """A cubic equation of state's shape.

    The model layer reads its geometry through one of these rather than inlining
    ``sqrt(2)`` the way the Peng-Robinson-only code once did, so a second cubic is a
    new instance rather than a new branch in every fugacity and Helmholtz expression.
    """

    #: The short name that crosses the Python boundary: ``"pr"`` or ``"srk"``.
    name: str
    omega_a: float
    omega_b: float
    delta1: float
    delta2: float
    delta_sum: float
    delta_prod: float
    delta_diff: float

    @property
    def half_delta_diff(self) -> float:
        """``(delta1 - delta2)/2``, the ``sqrt(2)`` the Helmholtz terms carry bare."""
        return self.delta_diff / 2.0

    def i_term(self, z: float, b: float) -> float:
        """``ln((z + delta1 B)/(z + delta2 B))``, the "I" term of the fugacity."""
        return math.log((z + self.delta1 * b) / (z + self.delta2 * b))

    def coefficient(self, a: float, b: float) -> float:
        """``A/((delta1 - delta2) B)``, the coefficient the fugacity term is scaled by."""
        return a / (self.delta_diff * b)

    def z_coefficients(self, a: float, b: float) -> tuple[float, float, float]:
        """The monic cubic ``z**3 + c2 z**2 + c1 z + c0`` in the reduced parameters."""
        ds = self.delta_sum
        dp = self.delta_prod
        c2 = (ds - 1.0) * b - 1.0
        c1 = dp * b * b - ds * b * b - ds * b + a
        c0 = -dp * b * b * b - dp * b * b - a * b
        return c2, c1, c0

    def df_dz(self, z: float, a: float, b: float) -> float:
        """``dF/dz`` of the z-cubic, at a root ``z``."""
        c2, c1, _ = self.z_coefficients(a, b)
        return 3.0 * z * z + 2.0 * c2 * z + c1

    def df_da(self, z: float, b: float) -> float:
        """``dF/dA`` of the z-cubic, at a root ``z``.

        ``A`` appears in ``c1`` with coefficient one and in ``c0`` as ``-A B``, so this
        is ``z - B`` for every cubic.
        """
        return z - b

    def df_db(self, z: float, a: float, b: float) -> float:
        """``dF/dB`` of the z-cubic, at a root ``z``."""
        ds = self.delta_sum
        dp = self.delta_prod
        # From `c2 = (delta1 + delta2 - 1) B - 1`, `c1` and `c0` in `z_coefficients`.
        dc2_db = ds - 1.0
        dc1_db = 2.0 * (dp - ds) * b - ds
        dc0_db = -3.0 * dp * b * b - 2.0 * dp * b - a
        return dc2_db * z * z + dc1_db * z + dc0_db

    def t_dfdt(self, z: float, a: float, b: float, t_da: float, t_db: float) -> float:
        """``T * dF/dT`` of the z-cubic, with ``t_da = T da/dT`` and ``t_db = T db/dT``."""
        ds = self.delta_sum
        dp = self.delta_prod
        # `T * d(c2)/dT = (delta1 + delta2 - 1) t_db`.
        t_dc2 = (ds - 1.0) * t_db
        # `T * d(c1)/dT = (2 B (delta1 delta2 - delta1 - delta2) - (delta1 + delta2)) t_db + t_da`.
        t_dc1 = (2.0 * b * (dp - ds) - ds) * t_db + t_da
        # `T * d(c0)/dT = (-3 delta1 delta2 B**2 - 2 delta1 delta2 B - A) t_db - t_da B`.
        t_dc0 = (-3.0 * dp * b * b - 2.0 * dp * b - a) * t_db - t_da * b
        return t_dc2 * z * z + t_dc1 * z + t_dc0

    def helmholtz_g(self, b: float) -> float:
        """``G(b) = ln((1 + delta1 b)/(1 + delta2 b))``, the Helmholtz geometry term."""
        return math.log((1.0 + self.delta1 * b) / (1.0 + self.delta2 * b))

    def helmholtz_q(self, b: float) -> float:
        """``q(b) = (1 + delta1 b)(1 + delta2 b)``, the denominator ``G'`` and ``G''`` share."""
        return 1.0 + self.delta_sum * b + self.delta_prod * b * b

    def helmholtz_g_prime(self, b: float) -> float:
        """``G'(b) = (delta1 - delta2)/q``."""
        return self.delta_diff / self.helmholtz_q(b)

    def helmholtz_g_second(self, b: float) -> float:
        """``G''(b) = -(delta1 - delta2)(delta1 + delta2 + 2 delta1 delta2 b)/q**2``."""
        q = self.helmholtz_q(b)
        return -self.delta_diff * (self.delta_sum + 2.0 * self.delta_prod * b) / (q * q)


#: Peng-Robinson: ``delta = (1 + sqrt(2), 1 - sqrt(2))``. The ``omega`` pair is NeqSim
#: 3.20.0's, not the paper's - see ``azoth.eos.reference.pr_alpha_ab.OMEGA_A``.
PR = Cubic(
    name="pr",
    omega_a=0.45724333333,
    omega_b=0.077803333,
    delta1=1.0 + _SQRT_2,
    delta2=1.0 - _SQRT_2,
    delta_sum=2.0,
    delta_prod=-1.0,
    delta_diff=2.0 * _SQRT_2,
)

#: Soave-Redlich-Kwong: ``delta = (1, 0)``. NeqSim's ``ComponentSrk`` computes the
#: ``omega`` pair from ``Math.pow(2.0, 1.0/3.0)`` at construction; the two decimals are
#: that expression to full double precision.
SRK = Cubic(
    name="srk",
    omega_a=0.4274802335403413,
    omega_b=0.08664034996495773,
    delta1=1.0,
    delta2=0.0,
    delta_sum=1.0,
    delta_prod=0.0,
    delta_diff=1.0,
)

#: Redlich-Kwong: the original, with the same ``omega`` and ``delta`` as Soave's but
#: ``alpha = 1/sqrt(Tr)`` instead of Soave's correlation.
RK = Cubic(
    name="rk",
    omega_a=0.4274802335403413,
    omega_b=0.08664034996495773,
    delta1=1.0,
    delta2=0.0,
    delta_sum=1.0,
    delta_prod=0.0,
    delta_diff=1.0,
)

#: The cubics this library runs, keyed by the short name that crosses the boundary.
CUBICS: dict[str, Cubic] = {"pr": PR, "srk": SRK, "rk": RK}
