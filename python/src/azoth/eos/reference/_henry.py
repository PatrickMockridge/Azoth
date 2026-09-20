"""The Henry reference state, as NeqSim's ``ComponentGE`` carries it.

The Python mirror of ``crates/azoth-eos/src/henry.rs``. **Held to the same constants
rather than to each other**: ``test_cross_impl.py`` compares registered models, and this
is a shared function, so both sides' tests assert ``H(298.15) = 51.6545 bar`` for CO2 and
the same finite cap - the numbers, not the agreement.

**There is no single branch to port.** Three of the tranche's models override
``ComponentGE.fugcoef`` and disagree about the reference state, so where a model is a
*solute* is that model's own arithmetic:

===========================  ==============================  ====================
NeqSim component             a neutral solute's ``phi``       an ion's ``phi``
===========================  ==============================  ====================
``ComponentGE``              ``(gamma / gamma_inf) H / P``   ``1e12 / P``
``ComponentGePitzer``        ``gamma * H * (m / x) / P``     ``1e12 / P``
``ComponentKentEisenberg``   ``H / P``, with ``gamma = 1``   ``1e8``
``ComponentDesmukhMather``   ``(gamma / gamma_inf) H / P``   ``1e-15``
===========================  ==============================  ====================

What they share is ``H(T)`` and the cap, and that is what this module is. The
constants an insoluble ion gets are three *different* numbers, and each belongs to the
model that chose it rather than to a rule here.

The data behind it is thin, and that is a finding: of the 348 compiled rows, 296 carry
no correlation and the other 52 carry one of four distinct sets - a ``900`` sentinel
whose ``exp`` overflows into the cap (40 rows), a constant (5), ``0.059565`` shared by
``h2s``, ``sf6``, ``r12`` and ``r134a`` (4), and CO2's real correlation, which is
shared with ``co`` and ``cos`` (3).
"""

from __future__ import annotations

import math
from typing import TYPE_CHECKING, NamedTuple

from azoth.core.errors import PropertyUnavailableError

if TYPE_CHECKING:
    from azoth.eos.components import DatabankEntry

#: The Henry coefficient a model uses for a substance with no usable correlation.
#:
#: ``ComponentGE.INSOLUBLE_HENRY_COEFFICIENT``, in bar. Effectively insoluble: the
#: coefficient enters as ``H / P``, so ``1e12 bar`` at a process pressure gives a
#: fugacity coefficient of order ``1e7`` and a mole fraction of order ``1e-7``.
INSOLUBLE_HENRY_COEFFICIENT = 1.0e12

#: NeqSim's own two factors, kept unmultiplied: water's molar mass in kg/mol and the
#: `100` beside it (`Component.java:2214`). Folded to `1.802` they would be one number
#: nobody could check against the source.
WATER_MOLAR_MASS_KG_PER_MOL = 0.01802
NEQSIM_FACTOR = 100.0


class HenryRecord(NamedTuple):
    """The four coefficients of NeqSim's Henry correlation, dimensionless, in bar.

    ``HenryCoef1``..``HenryCoef4`` of ``COMP.csv``, which ``Component.java:2210`` reads as
    ``henryCoefParameter[0..3]``. **All four are zero on 296 of the 348 rows**, which is
    the table's marker for a substance it fits no correlation for - and a zero set
    evaluated as a polynomial gives ``1.802 bar`` for every one of them, so a caller must
    refuse rather than evaluate. :func:`effective_coefficient` is that refusal.
    """

    h0: float
    h1: float
    h2: float
    h3: float

    def is_fitted(self) -> bool:
        """Whether the table fits this substance at all."""
        return self.h0 != 0.0 or self.h1 != 0.0 or self.h2 != 0.0 or self.h3 != 0.0


def coefficient(record: HenryRecord, temperature_k: float) -> float:
    """``H(T)``, the correlation alone, in bar.

    NeqSim's own expression (``Component.java:2214``):

    .. code-block:: text

        H(T) = exp(h0 + h1/T + h2 ln T + h3 T) * 0.01802 * 100

    The value is **not** capped - see :func:`effective_coefficient`, which is what a
    model takes. A ``900`` constant overflows ``exp`` to infinity, which is the table's
    sentinel for a substance it lists and does not fit; that is not an error here,
    because the cap turns it into :data:`INSOLUBLE_HENRY_COEFFICIENT`.
    """
    exponent = (
        record.h0
        + record.h1 / temperature_k
        + record.h2 * math.log(temperature_k)
        + record.h3 * temperature_k
    )
    try:
        return math.exp(exponent) * WATER_MOLAR_MASS_KG_PER_MOL * NEQSIM_FACTOR
    except OverflowError:
        # `math.exp` raises where Rust's `.exp()` returns infinity. The cap decides both,
        # so the two languages agree; raising here would refuse 40 rows NeqSim evaluates.
        return math.inf


def coefficient_dt(record: HenryRecord, temperature_k: float) -> float:
    """``dH/dT``, in bar per kelvin.

    NeqSim's ``getHenryCoefdT``, which is ``H(T)`` times the correlation's own logarithmic
    derivative - so a caller wanting ``d(ln H)/dT`` divides by :func:`coefficient`.
    """
    return coefficient(record, temperature_k) * (
        -record.h1 / (temperature_k * temperature_k) + record.h2 / temperature_k + record.h3
    )


def is_capped(entry: DatabankEntry, value: float) -> bool:
    """Whether a correlation must fail closed to the insoluble limit.

    ``ComponentGE.isHenryCoefficientCapped``: a value that is not finite, not positive, or
    above the limit, **or a substance the table classes an ion**. The ion clause is what
    makes an ion insoluble whatever its row says, and it is why the cap is not merely a
    numerical guard - 40 of the 52 fitted rows reach it through the ``900`` sentinel, and
    62 more through their class.
    """
    from azoth.eos.components import ION

    return (
        not math.isfinite(value)
        or value <= 0.0
        or value > INSOLUBLE_HENRY_COEFFICIENT
        or entry.component_type == ION
    )


def effective_coefficient(entry: DatabankEntry, temperature_k: float) -> float:
    """The Henry coefficient a model takes: ``ComponentGE.getEffectiveHenryCoefficient``.

    The database arm of NeqSim's two, and the only one this library ports. The other
    selects the IAPWS pure-water table for a supported neutral solute in a water-bearing
    phase, which is a second correlation over a second table and is refused by name rather
    than approximated here.

    Raises:
        PropertyUnavailableError: if the table fits the substance no correlation - the
            four zeros on 296 rows. Refused rather than evaluated, because the zeros
            would give ``1.802 bar`` for every such substance.
    """
    if not entry.henry.is_fitted():
        raise PropertyUnavailableError(
            entry.name,
            "a Henry coefficient",
            "the compiled table carries `HenryCoef1..4` as zeros for it, which is the "
            "table's marker for a substance it fits no correlation for. Evaluating the "
            "zeros as a polynomial gives 1.802 bar for every such substance, which is a "
            "number rather than an absence",
        )
    value = coefficient(entry.henry, temperature_k)
    return INSOLUBLE_HENRY_COEFFICIENT if is_capped(entry, value) else value
