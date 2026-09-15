"""The temperature dependence of the attraction term, as a strategy object.

NeqSim's ``AttractiveTerm*`` classes own ``alpha(T)`` and its derivatives; this is the
port of that seam. A cubic's shape is in :mod:`azoth.eos.cubic` and its
temperature dependence is here, so a new alpha correlation is a new type rather than a
new branch in the model layer.
"""

from __future__ import annotations

import math


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
