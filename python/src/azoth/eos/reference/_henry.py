"""The Henry reference state, as NeqSim's ``ComponentGE`` carries it.

The Python mirror of ``crates/azoth-eos/src/henry.rs``. **Held to the same constants
rather than to each other**: ``test_cross_impl.py`` compares registered models, and this
is a shared function, so both sides' tests assert ``H(298.15) = 28.95 bar`` for CO2 and
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

The data behind it is sourced, and its history is worth keeping: of the 348 compiled
rows, 269 carry no correlation and the other 79 carry 75 distinct two-term sets.
**Upstream PR #4045 replaced the four placeholder sets this module used to describe** -
forty rows of ``900``, the constant ``2.92``, the ``0.059565`` quartet and the one CO2 set
copied onto ``co`` and ``cos`` - with sourced coefficients. Methane is now ``723.75 bar``
at 298.15 K rather than ``exp(900)`` through the cap, and CO2 is ``28.95 bar`` where the
copied set gave ``51.65``. The cap remains the model's rule for an ion.
"""

from __future__ import annotations

import enum
import math
from typing import TYPE_CHECKING, NamedTuple

from azoth.core.errors import PropertyUnavailableError

if TYPE_CHECKING:
    from azoth.eos.components import DatabankEntry


class Insoluble(enum.Enum):
    """What a model gives an ion it has no data for, by model.

    **Three numbers, and none of them a rule.** ``ComponentGE`` and ``ComponentGePitzer``
    give an insoluble ion the capped Henry coefficient - ``1e12`` bar, which is also the
    number the old ``exp(900)`` sentinel overflowed to, so the cap and the sentinel were
    one constant in two roles - and the other two models each carry a fugacity coefficient
    of their own.
    ``ComponentKentEisenberg``'s ``1e8`` says "sparingly present" and
    ``ComponentDesmukhMather``'s ``1e-15`` says "absent", which are different claims about
    the same ion.

    They were three loose floats under two names, both spelled ``INSOLUBLE_ION``, in a
    package where a ``from ... import`` decides which one a caller gets: ``1e8`` against
    ``1e-15`` is a factor of ``1e23``, it changes a vapour composition rather than raising,
    and nothing compared the two modules that each defined the name. A model now names the
    one it takes, and this module is where the three are stated.
    """

    #: The Henry coefficient an ion is capped to, **in bar** - and the cap itself, which
    #: :func:`is_capped` compares against. ``ComponentGE.INSOLUBLE_HENRY_COEFFICIENT``.
    HENRY_COEFFICIENT = 1.0e12
    #: ``ComponentKentEisenberg.fugcoef``'s constant, **dimensionless**.
    KENT_EISENBERG = 1.0e8
    #: ``ComponentDesmukhMather.fugcoef``'s constant, **dimensionless**, the smallest of
    #: the three.
    DESMUKH_MATHER = 1.0e-15


#: NeqSim's own two factors, kept unmultiplied: water's molar mass in kg/mol and the
#: `100` beside it (`Component.java:2214`). Folded to `1.802` they would be one number
#: nobody could check against the source.
WATER_MOLAR_MASS_KG_PER_MOL = 0.01802
NEQSIM_FACTOR = 100.0


class HenryRecord(NamedTuple):
    """The four coefficients of NeqSim's Henry correlation, dimensionless, in bar.

    ``HenryCoef1``..``HenryCoef4`` of ``COMP.csv``, which ``Component.java:2210`` reads as
    ``henryCoefParameter[0..3]``. **All four are zero on 269 of the 348 rows**, which is
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
    because the cap turns it into :attr:`Insoluble.HENRY_COEFFICIENT`.
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

    return is_capped_value(value, entry.component_type == ION)


def is_capped_value(value: float, ion: bool) -> bool:
    """:func:`is_capped` from the value and the one component fact it reads.

    A GE phase's resolved parameters carry the ion flag rather than the entry.
    """
    return (
        not math.isfinite(value) or value <= 0.0 or value > Insoluble.HENRY_COEFFICIENT.value or ion
    )


def effective_coefficient(entry: DatabankEntry, temperature_k: float) -> float:
    """The Henry coefficient a model takes: ``ComponentGE.getEffectiveHenryCoefficient``.

    The database arm of NeqSim's two. The other is the IAPWS pure-water table, which
    :mod:`azoth.eos.reference.iapws_henry_law` carries; the GE phases select it through
    ``_ge_phase.henry_coefficient``, and this function is the database arm alone.

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
    return Insoluble.HENRY_COEFFICIENT.value if is_capped(entry, value) else value
