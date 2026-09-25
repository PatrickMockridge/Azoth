"""Standard and regulatory calculations.

The namespace NeqSim files under ``standards/``, whose first member is ISO 6976: the
calorific values, density and relative density of a natural gas from its composition. It is
here rather than in :mod:`azoth.eos` because it is not a property the equation of state knows
- it is a *standard's* tabulated constants, summed.

See :func:`azoth.standards.iso6976`.
"""

from __future__ import annotations

from azoth.core.result import Iso6976Result

__all__ = ["Iso6976Result", "iso6976"]


def iso6976(
    components: list[str],
    z: list[float],
    volumetric_reference_temperature: object,
    energy_reference_temperature: object,
) -> Iso6976Result:
    """The calorific values and density of a natural gas, by ISO 6976.

    ``components`` and ``z`` are the gas. The two reference temperatures are **in kelvin**:
    the volumetric one is the temperature the densities and the relative density are formed
    at, and the energy one the temperature the calorific values are stated at. The standard
    states both in degrees Celsius and the offset is exact.

    See :func:`azoth.standards.reference.iso6976`.
    """
    from azoth._dispatch import resolve

    return resolve("standards.iso6976")(  # type: ignore[no-any-return]
        components=components,
        z=z,
        volumetric_reference_temperature=volumetric_reference_temperature,
        energy_reference_temperature=energy_reference_temperature,
    )
