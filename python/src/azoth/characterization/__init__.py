"""Petroleum-fraction characterisation.

The front door for a fluid that has no databank row. A crude is an *assay* - a list of cuts,
each a molar mass and a density - and everything a cubic needs about a cut follows from those
two numbers by correlation. NeqSim's ``thermo/characterization/`` is that subsystem.

It is a namespace of its own rather than a family under :mod:`azoth.eos`, mirroring NeqSim's
own package split: what is here is not an equation of state and is not reached from one.

* :func:`tbp_cut_properties` - a cut's critical properties, by any of NeqSim's ten models

# Which implementation answers

As in the other namespaces, each function dispatches to the Rust extension when it is built and
to :mod:`azoth.characterization.reference` otherwise. Both are always reachable - see
:func:`azoth.backends` and :func:`azoth.use_backend`.
"""

from __future__ import annotations

from azoth._dispatch import resolve
from azoth.core.result import TbpCutPropertiesResult
from azoth.core.units import Q

__all__ = [
    "tbp_cut_properties",
]

_TBP_CUT_PROPERTIES = "characterization.tbp_cut_properties"


def tbp_cut_properties(
    molar_mass: Q,
    density: Q,
    model: str | None = None,
    boiling_point: Q | None = None,
) -> TbpCutPropertiesResult:
    """A TBP cut's critical properties, by any of NeqSim's ten models.

    ``model`` is one of ``pedersen_srk``, ``pedersen_srk_heavy_oil``, ``pedersen_pr``,
    ``pedersen_pr2``, ``pedersen_pr_heavy_oil``, ``riazi_daubert``, ``lee_kesler``, ``twu``,
    ``cavett`` or ``standing``; ``None`` takes NeqSim's own default, ``pedersen_srk``.

    ``attraction_exponent`` is ``None`` for the five models that have no exponent to offer.

    Raises:
        OutOfRangeError: if ``molar_mass`` or ``density`` is not positive, or if a supplied
            ``boiling_point`` is not positive.

    See :func:`azoth.characterization.reference.tbp_cut_properties`.
    """
    return resolve(_TBP_CUT_PROPERTIES)(  # type: ignore[no-any-return]
        model=model, molar_mass=molar_mass, density=density, boiling_point=boiling_point
    )
