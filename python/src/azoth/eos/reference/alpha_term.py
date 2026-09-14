"""The temperature dependence of the attraction term, as a strategy object.

NeqSim's ``AttractiveTerm*`` classes own ``alpha(T)`` and its derivatives; this is the
port of that seam. A cubic's shape is in :mod:`azoth.eos.reference.cubic` and its
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
