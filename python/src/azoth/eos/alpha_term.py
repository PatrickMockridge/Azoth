"""The temperature dependence of the attraction term, as a strategy object.

NeqSim's ``AttractiveTerm*`` classes own ``alpha(T)`` and its derivatives; this is the
port of that seam. A cubic's shape is in :mod:`azoth.eos.cubic` and its
temperature dependence is here, so a new alpha correlation is a new type rather than a
new branch in the model layer.
"""

from __future__ import annotations

import math
from enum import Enum


class Soave:
    """Soave's correlation: ``alpha = (1 + m(1 - sqrt(Tr)))**2``.

    Peng-Robinson and Soave-Redlich-Kwong share this form and differ only in the
    coefficient ``m`` - PR's ``0.37464 + 1.54226 w - 0.26992 w**2``, SRK's
    ``0.48 + 1.574 w - 0.176 w**2`` - so the term carries ``m`` as a plain number and
    the correlation that produced it lives in the ``kappa`` calc.
    """

    def __init__(self, kappa: float) -> None:
        self.kappa = kappa

    def alpha(self, tr: float) -> float:
        """The attraction scale at reduced temperature ``Tr``."""
        attraction = 1.0 + self.kappa * (1.0 - math.sqrt(tr))
        return attraction**2

    def psi(self, tr: float) -> float:
        """``d ln alpha / d ln T``, the logarithmic derivative."""
        sqrt_tr = math.sqrt(tr)
        return -self.kappa * sqrt_tr / (1.0 + self.kappa * (1.0 - sqrt_tr))

    def psi_t(self, tr: float) -> float:
        """``T * d(psi)/dT``, already multiplied by ``T``."""
        sqrt_tr = math.sqrt(tr)
        return (
            -self.kappa
            * (1.0 + self.kappa)
            * tr
            / (2.0 * sqrt_tr * (1.0 + self.kappa * (1.0 - sqrt_tr)) ** 2)
        )


class Danesh:
    """Danesh's correlation: Soave's form with ``m`` scaled by 1.21 above the critical
    temperature. The base ``m`` is ``eos.pr78_kappa``.
    """

    def __init__(self, kappa: float) -> None:
        self.kappa = kappa

    def _m_mod(self, tr: float) -> float:
        return self.kappa * 1.21 if tr > 1.0 else self.kappa

    def alpha(self, tr: float) -> float:
        m_mod = self._m_mod(tr)
        return (1.0 + m_mod * (1.0 - math.sqrt(tr))) ** 2

    def psi(self, tr: float) -> float:
        m_mod = self._m_mod(tr)
        sqrt_tr = math.sqrt(tr)
        return -m_mod * sqrt_tr / (1.0 + m_mod * (1.0 - sqrt_tr))

    def psi_t(self, tr: float) -> float:
        m_mod = self._m_mod(tr)
        sqrt_tr = math.sqrt(tr)
        return -m_mod * (1.0 + m_mod) * tr / (2.0 * sqrt_tr * (1.0 + m_mod * (1.0 - sqrt_tr)) ** 2)


class RkAlpha:
    """Redlich-Kwong's original correlation: ``alpha = 1/sqrt(Tr)``.

    Kappa-free, which is the whole difference from :class:`Soave`: the logarithmic
    derivative is the constant ``-1/2`` and its temperature derivative is zero.
    """

    def alpha(self, tr: float) -> float:
        return 1.0 / math.sqrt(tr)

    def psi(self, _tr: float) -> float:
        return -0.5

    def psi_t(self, _tr: float) -> float:
        return 0.0


class TwuCoon:
    """Twu-Coon's correlation, a non-Soave ``alpha`` with the acentric factor as ``m``."""

    def __init__(self, omega: float) -> None:
        self.omega = omega

    def _values(self, tr: float) -> tuple[float, float, float]:
        tr_c = tr**2.29528
        tr_f = tr**2.63165
        a_term = tr**-0.201158 * math.exp(0.141599 * (1.0 - tr_c))
        d_term = tr**-0.660145 * math.exp(0.500315 * (1.0 - tr_f))

        la = -0.201158 / tr - 0.141599 * 2.29528 * tr_c / tr
        ld = -0.660145 / tr - 0.500315 * 2.63165 * tr_f / tr

        alpha = a_term + self.omega * (d_term - a_term)
        d_alpha = a_term * la + self.omega * (d_term * ld - a_term * la)

        la2 = 0.201158 / (tr * tr) - 0.141599 * 2.29528 * 1.29528 * tr_c / (tr * tr)
        ld2 = 0.660145 / (tr * tr) - 0.500315 * 2.63165 * 1.63165 * tr_f / (tr * tr)
        d2_alpha = a_term * (la * la + la2) + self.omega * (
            d_term * (ld * ld + ld2) - a_term * (la * la + la2)
        )
        return alpha, d_alpha, d2_alpha

    def alpha(self, tr: float) -> float:
        return self._values(tr)[0]

    def psi(self, tr: float) -> float:
        alpha, d_alpha, _ = self._values(tr)
        return tr * d_alpha / alpha

    def psi_t(self, tr: float) -> float:
        alpha, d_alpha, d2_alpha = self._values(tr)
        return tr * d_alpha / alpha + tr * tr * (d2_alpha * alpha - d_alpha * d_alpha) / (
            alpha * alpha
        )


class Gassem2001:
    """Gassem et al. (2001)'s correlation, a non-Soave ``alpha`` with the acentric factor
    in the exponent: ``alpha = exp((A + B Tr)(1 - Tr**g))`` with ``g = C + D w + E w**2``.
    """

    A = 2.0
    B = 0.836
    C = 0.134
    D = 0.508
    E = -0.0467

    def __init__(self, omega: float) -> None:
        self.omega = omega

    def _g(self) -> float:
        return self.C + self.D * self.omega + self.E * self.omega**2

    def alpha(self, tr: float) -> float:
        return math.exp((self.A + self.B * tr) * (1.0 - math.pow(tr, self._g())))

    def psi(self, tr: float) -> float:
        g = self._g()
        return tr * (
            self.B * (1.0 - math.pow(tr, g)) - g * math.pow(tr, g - 1.0) * (self.A + self.B * tr)
        )

    def psi_t(self, tr: float) -> float:
        g = self._g()
        tr_g = math.pow(tr, g)
        return self.B * tr - g * g * self.A * tr_g - self.B * (1.0 + g) ** 2 * tr_g * tr


def _squared_psi_t(s: float, ds: float, dds: float, tr: float) -> float:
    """`T * d(psi)/dT` for a squared form `alpha = S**2`, from `S` and its derivatives."""
    return 2.0 * tr * ds / s + 2.0 * tr * tr * (s * dds - ds * ds) / (s * s)


def _matcop_polynomial(params: tuple[float, ...], tr: float) -> tuple[float, float, float]:
    """The Mathias-Copeman polynomial `S = 1 + sum_k c_k (1 - sqrt(Tr))**k` and its
    first two `Tr` derivatives."""
    sqrt_tr = math.sqrt(tr)
    u = 1.0 - sqrt_tr
    du = -0.5 / sqrt_tr
    ddu = 0.25 / (tr * sqrt_tr)
    s = 1.0
    ds = 0.0
    dds = 0.0
    u_pow = 1.0
    u_pow_minus = 0.0
    for i, c in enumerate(params):
        k = i + 1.0
        u_k = u_pow * u
        s += c * u_k
        ds += c * k * u_pow * du
        dds += c * (k * (k - 1.0) * u_pow_minus * du * du + k * u_pow * ddu)
        u_pow_minus = u_pow
        u_pow = u_k
    return s, ds, dds


def matcop_kappa(omega: float) -> float:
    """Mathias-Copeman's base ``m`` for the SRK variant."""
    return 0.48 + 1.574 * omega - 0.175 * omega * omega


def umr_kappa(omega: float) -> float:
    """The UMR-PRU quartic ``m``."""
    return (
        0.384401 + 1.52276 * omega - 0.213808 * omega**2 + 0.034616 * omega**3 - 0.001976 * omega**4
    )


class Schwartzentruber:
    """Schwartzentruber's correlation, a squared form with three fitted parameters."""

    def __init__(self, omega: float, params: tuple[float, ...]) -> None:
        self.omega = omega
        self.params = params

    def _kappa(self) -> float:
        return 0.48508 + 1.55191 * self.omega - 0.15613 * self.omega**2

    def _polynomial(self, tr: float) -> tuple[float, float, float]:
        m = self._kappa()
        p0 = self.params[0] if self.params else 0.0
        p1 = self.params[1] if len(self.params) > 1 else 0.0
        p2 = self.params[2] if len(self.params) > 2 else 0.0
        sqrt_tr = math.sqrt(tr)
        q = 1.0 + (p1 - 1.0) * tr + (p2 - p1) * tr * tr - p2 * tr * tr * tr
        dq = p1 - 1.0 + 2.0 * (p2 - p1) * tr - 3.0 * p2 * tr * tr
        ddq = 2.0 * (p2 - p1) - 6.0 * p2 * tr
        s = 1.0 + m * (1.0 - sqrt_tr) - p0 * q
        ds = -m / (2.0 * sqrt_tr) - p0 * dq
        dds = m / (4.0 * tr * sqrt_tr) - p0 * ddq
        return s, ds, dds

    def alpha(self, tr: float) -> float:
        s, _, _ = self._polynomial(tr)
        return s * s

    def psi(self, tr: float) -> float:
        s, ds, _ = self._polynomial(tr)
        return 2.0 * tr * ds / s

    def psi_t(self, tr: float) -> float:
        s, ds, dds = self._polynomial(tr)
        return _squared_psi_t(s, ds, dds, tr)


class Mollerup:
    """Mollerup's correlation, a three-parameter non-Soave form."""

    def __init__(self, params: tuple[float, ...]) -> None:
        self.params = params

    def _values(self, tr: float) -> tuple[float, float, float]:
        p0 = self.params[0] if self.params else 0.0
        p1 = self.params[1] if len(self.params) > 1 else 0.0
        p2 = self.params[2] if len(self.params) > 2 else 0.0
        ln_tr = math.log(tr)
        alpha = 1.0 + p0 * (1.0 / tr - 1.0) + p1 * tr * ln_tr + p2 * (tr - 1.0)
        d_alpha = -p0 / (tr * tr) + p1 * (ln_tr + 1.0) + p2
        d2_alpha = 2.0 * p0 / (tr * tr * tr) + p1 / tr
        return alpha, d_alpha, d2_alpha

    def alpha(self, tr: float) -> float:
        return self._values(tr)[0]

    def psi(self, tr: float) -> float:
        alpha, d_alpha, _ = self._values(tr)
        return tr * d_alpha / alpha

    def psi_t(self, tr: float) -> float:
        alpha, d_alpha, d2_alpha = self._values(tr)
        return tr * d_alpha / alpha + tr * tr * (d2_alpha * alpha - d_alpha * d_alpha) / (
            alpha * alpha
        )


class MatCopFallback(Enum):
    """When Mathias-Copeman's polynomial is replaced by the base Soave alpha."""

    NONE = "none"
    SUPERCRITICAL = "supercritical"
    UNSET = "unset"
    ALL_UNSET = "all_unset"


class MatCop:
    """Mathias-Copeman's correlation, a squared polynomial with three or five fitted
    parameters."""

    def __init__(self, kappa: float, params: tuple[float, ...], fallback: MatCopFallback) -> None:
        self.kappa = kappa
        self.params = params
        self.fallback = fallback

    def _coefficients(self) -> list[float]:
        c = list(self.params)
        if c and abs(c[0]) < 1e-12:
            c[0] = self.kappa
        return c

    def _use_soave(self, c: list[float], tr: float) -> bool:
        if self.fallback is MatCopFallback.NONE:
            return False
        if self.fallback is MatCopFallback.SUPERCRITICAL:
            return tr > 1.0 or (c[0] < 1e-20 if c else True)
        if self.fallback is MatCopFallback.UNSET:
            return c[0] < 1e-20 if c else True
        return all(abs(p) < 1e-20 for p in c)

    def alpha(self, tr: float) -> float:
        c = self._coefficients()
        if self._use_soave(c, tr):
            return Soave(self.kappa).alpha(tr)
        s, _, _ = _matcop_polynomial(tuple(c), tr)
        return s * s

    def psi(self, tr: float) -> float:
        c = self._coefficients()
        if self._use_soave(c, tr):
            return Soave(self.kappa).psi(tr)
        s, ds, _ = _matcop_polynomial(tuple(c), tr)
        return 2.0 * tr * ds / s

    def psi_t(self, tr: float) -> float:
        c = self._coefficients()
        if self._use_soave(c, tr):
            return Soave(self.kappa).psi_t(tr)
        s, ds, dds = _matcop_polynomial(tuple(c), tr)
        return _squared_psi_t(s, ds, dds, tr)


class Delft1998:
    """Delft (1998)'s correlation: a fitted cubic for methane, the 1978 Peng-Robinson
    Soave form for everything else."""

    def __init__(self, kappa: float, is_methane: bool) -> None:
        self.kappa = kappa
        self.is_methane = is_methane

    def alpha(self, tr: float) -> float:
        if self.is_methane:
            return 0.969617 + 0.20089 * tr - 0.3256987 * tr * tr + 0.06653 * tr * tr * tr
        return Soave(self.kappa).alpha(tr)

    def psi(self, tr: float) -> float:
        if self.is_methane:
            alpha = 0.969617 + 0.20089 * tr - 0.3256987 * tr * tr + 0.06653 * tr * tr * tr
            d_alpha = 0.20089 - 2.0 * 0.3256987 * tr + 3.0 * 0.06653 * tr * tr
            return tr * d_alpha / alpha
        return Soave(self.kappa).psi(tr)

    def psi_t(self, tr: float) -> float:
        if self.is_methane:
            alpha = 0.969617 + 0.20089 * tr - 0.3256987 * tr * tr + 0.06653 * tr * tr * tr
            d_alpha = 0.20089 - 2.0 * 0.3256987 * tr + 3.0 * 0.06653 * tr * tr
            d2_alpha = -2.0 * 0.3256987 + 6.0 * 0.06653 * tr
            return tr * d_alpha / alpha + tr * tr * (d2_alpha * alpha - d_alpha * d_alpha) / (
                alpha * alpha
            )
        return Soave(self.kappa).psi_t(tr)
