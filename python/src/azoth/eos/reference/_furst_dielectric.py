"""The dielectric surface a Furst electrolyte phase is built on: the Python reference.

A private module, the mirror of ``crates/azoth-eos/src/furst_dielectric.rs``. The
arithmetic is the same in both languages and the cross-implementation tests are what hold
them together; this docstring does not restate what the Rust module's does, except for the
two things a reader of *this* file needs:

* **NeqSim's own constants**, not the standard ones. ``avagadroNumber`` is ``6.023e23``
  where CODATA says ``6.02214076e23``, and ``pi`` is ``3.14159265`` where ``math.pi``
  carries sixteen digits. Both land in the ninth significant digit of an ion's radius, and
  the truncated ``pi`` is the one a port would be least likely to notice.
* **The other two mixing rules**: ``VOLUME_AVERAGE`` and ``LOOYENGA`` are implemented and
  reachable, and NeqSim's own comment says neither is thermodynamically consistent with
  complete composition derivatives - which is why ``MOLAR_AVERAGE`` is the default and the
  only one whose derivatives this library assembles.
"""

from __future__ import annotations

import math
from enum import Enum

from azoth.core.errors import InvalidInputError, OutOfRangeError

#: NeqSim's ``avagadroNumber``, from ``ThermodynamicConstantsInterface``.
NEQSIM_AVOGADRO = 6.023e23

#: NeqSim's ``pi``, from the same interface. Not ``math.pi``.
NEQSIM_PI = 3.14159265

#: The reference temperature ``Wij(T)`` is written against.
T_REFERENCE = 298.15


class MixingRule(Enum):
    """How the solvents' dielectric constants are combined."""

    molar_average = "MOLAR_AVERAGE"
    volume_average = "VOLUME_AVERAGE"
    looyenga = "LOOYENGA"


def default_for_the_model() -> MixingRule:
    """The rule ``PhaseModifiedFurstElectrolyteEos`` is constructed with."""
    return MixingRule.molar_average


def component_dielectric(coefficients: tuple[float, ...], temperature: float) -> float:
    """``eps_i(T)``, one component's dielectric constant from its five coefficients."""
    d0, d1, d2, d3, d4 = coefficients
    return (
        d0
        + d1 / temperature
        + d2 * temperature
        + d3 * temperature * temperature
        + d4 * temperature**3
    )


def component_dielectric_dt(coefficients: tuple[float, ...], temperature: float) -> float:
    """``d eps_i / dT``."""
    _, d1, d2, d3, d4 = coefficients
    return -d1 / temperature**2 + d2 + 2.0 * d3 * temperature + 3.0 * d4 * temperature**2


def component_dielectric_dtdt(coefficients: tuple[float, ...], temperature: float) -> float:
    """``d^2 eps_i / dT^2``."""
    _, d1, _, d3, d4 = coefficients
    return 2.0 * d1 / temperature**3 + 2.0 * d3 + 6.0 * d4 * temperature


def _checked_lengths(mole_numbers: list[float], eps_i: list[float], is_ion: list[bool]) -> None:
    if len(mole_numbers) == len(eps_i) == len(is_ion):
        return
    raise InvalidInputError(
        "components",
        f"{len(mole_numbers)} mole numbers, {len(eps_i)} dielectric constants and "
        f"{len(is_ion)} ion flags; the three are one per component",
    )


def solvent_dielectric(
    rule: MixingRule,
    mole_numbers: list[float],
    eps_i: list[float],
    is_ion: list[bool],
    critical_volumes: list[float],
) -> float:
    """The solvent mixture's dielectric constant under the rule named.

    ``critical_volumes`` is read only by the two volume-fraction rules.

    Raises:
        InvalidInputError: if the three slices are not one length.
        OutOfRangeError: if no solvent component carries any moles, which is the divisor
            and where NeqSim returns a number of order ``1e50`` instead.
    """
    _checked_lengths(mole_numbers, eps_i, is_ion)
    solvent_total = sum(n for n, ion in zip(mole_numbers, is_ion, strict=True) if not ion)
    if not math.isfinite(solvent_total) or solvent_total <= 0.0:
        raise OutOfRangeError(
            "mole_numbers",
            solvent_total,
            "a dielectric constant is a mole-number average over the solvent components, "
            "and none of them carries any moles. NeqSim divides by its `1e-50` floor here, "
            "which returns a number of order `1e50` that is not a dielectric constant",
        )
    if rule is MixingRule.molar_average:
        return (
            sum(n * eps for n, eps, ion in zip(mole_numbers, eps_i, is_ion, strict=True) if not ion)
            / solvent_total
        )

    if len(critical_volumes) != len(mole_numbers):
        raise InvalidInputError(
            "critical_volumes",
            f"the rule {rule.value} weights by volume fraction and needs one critical "
            f"volume per component, but got {len(critical_volumes)} for "
            f"{len(mole_numbers)} components",
        )
    total_volume = sum(
        n * vc for n, vc, ion in zip(mole_numbers, critical_volumes, is_ion, strict=True) if not ion
    )
    if not math.isfinite(total_volume) or total_volume <= 0.0:
        raise OutOfRangeError(
            "critical_volumes",
            total_volume,
            "the volume-weighted dielectric constant divides by the solvent's total "
            "critical volume and it came out zero",
        )
    total = 0.0
    for n, eps, vc, ion in zip(mole_numbers, eps_i, critical_volumes, is_ion, strict=True):
        if ion:
            continue
        fraction = n * vc / total_volume
        total += fraction * (math.cbrt(eps) if rule is MixingRule.looyenga else eps)
    return total**3 if rule is MixingRule.looyenga else total


def packing_fraction(
    diameters_m: list[float],
    mole_numbers: list[float],
    is_ion: list[bool],
    total_moles: float,
    molar_volume: float,
    *,
    ions_only: bool,
) -> float:
    """The packing fraction ``(N_A pi/6) sum_i n_i sigma_i^3 / (n V)``.

    ``diameters`` are in metres, ``molar_volume`` in m3/mol.

    Raises:
        InvalidInputError: if the slices are not one length.
        OutOfRangeError: if ``total_moles * molar_volume`` is not positive.
    """
    if len(diameters_m) != len(mole_numbers) or len(is_ion) != len(mole_numbers):
        raise InvalidInputError(
            "diameters",
            f"{len(diameters_m)} diameters, {len(mole_numbers)} mole numbers and "
            f"{len(is_ion)} ion flags; the three are one per component",
        )
    volume = total_moles * molar_volume
    if not math.isfinite(volume) or volume <= 0.0:
        raise OutOfRangeError(
            "molar_volume",
            volume,
            "the packing fraction is a volume ratio and its denominator is the phase's "
            "total volume, which must be positive",
        )
    total = sum(
        n * diameter**3
        for diameter, n, ion in zip(diameters_m, mole_numbers, is_ion, strict=True)
        if (ion or not ions_only)
    )
    return NEQSIM_AVOGADRO * NEQSIM_PI / 6.0 * total / volume


def packing_fraction_dv(packing: float, total_moles: float, molar_volume: float) -> float:
    """``d eps / dV`` for a packing fraction, at constant composition."""
    return -packing / (total_moles * molar_volume)


def packing_fraction_dvdv(packing: float, total_moles: float, molar_volume: float) -> float:
    """``d^2 eps / dV^2``, which for ``C / V`` is ``2 eps / V^2``."""
    return 2.0 * packing / (total_moles * molar_volume) ** 2


def phase_dielectric(solvent: float, ionic_packing: float) -> float:
    """The phase's dielectric constant, with the ions' own volume taken out of it."""
    return 1.0 + (solvent - 1.0) * (1.0 - ionic_packing) / (1.0 + ionic_packing / 2.0)
