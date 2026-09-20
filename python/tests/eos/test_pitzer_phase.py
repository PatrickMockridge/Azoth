"""Pitzer's activity-coefficient arithmetic, as `_pitzer_phase` carries it.

The Python mirror of the Rust tests in `crates/azoth-eos/src/pitzer_phase.rs`, asserting
the same constants - the rule a *shared* function follows rather than a registered model.

The end-to-end numbers come from `validation/neqsim/PitzerArithmetic.java` on a live
`SystemPitzer` at 298.15 K, so they are NeqSim's rather than a restatement of its formulas.
"""

from __future__ import annotations

import pytest

from azoth.eos.reference import _pitzer_catalog as catalog
from azoth.eos.reference import _pitzer_phase as phase

WATER_MOLAR_MASS = 0.018015


def molalities(x: list[float]) -> list[float]:
    """The molalities of a composition, per kilogram of water."""
    mass_of_water = x[0] * WATER_MOLAR_MASS
    return [xi / mass_of_water for xi in x]


class CatalogParameters:
    """The catalogue-backed pair parameters, which is what a selected brine reads.

    A zero for an absent row is NeqSim's own behaviour - its parameter arrays are
    zero-initialised and only the rows it read are written.
    """

    def __init__(self, names: list[str]) -> None:
        self.species = [catalog.canonical_species(name) for name in names]

    def _form(self, family: str, names: list[str]) -> tuple[float, ...]:
        return catalog.find(family, names) or (0.0,) * 6

    def beta0(self, first: int, second: int, temperature: float) -> float:
        return phase.catalog_value(
            self._form("B0", [self.species[first], self.species[second]]), temperature
        )

    def beta1(self, first: int, second: int, temperature: float) -> float:
        return phase.catalog_value(
            self._form("B1", [self.species[first], self.species[second]]), temperature
        )

    def cphi(self, first: int, second: int, temperature: float) -> float:
        return phase.catalog_value(
            self._form("C0", [self.species[first], self.species[second]]), temperature
        )

    def beta2(self, first: int, second: int, temperature: float) -> float:
        return phase.catalog_value(
            self._form("B2", [self.species[first], self.species[second]]), temperature
        )

    def theta(self, first: int, second: int, temperature: float) -> float:
        return phase.catalog_value(
            self._form("THETA", [self.species[first], self.species[second]]), temperature
        )

    def psi(self, first: int, second: int, third: int, temperature: float) -> float:
        return phase.catalog_value(
            self._form("PSI", [self.species[first], self.species[second], self.species[third]]),
            temperature,
        )


def test_the_debye_huckel_parameter_matches_neqsim() -> None:
    """The oracle's column, and the literature's `0.3915` against NeqSim's `0.392034`."""
    for temperature, expected in (
        (273.15, 0.377466966949275),
        (298.15, 0.392034451863750),
        (373.15, 0.456800085250619),
    ):
        assert phase.debye_huckel_a_phi(temperature) == pytest.approx(expected, abs=1e-14)


def test_the_ion_activity_coefficient_matches_neqsim() -> None:
    """water + Na+ + Cl-, the plain path: no same-sign pair, no non-2:2 `beta2`."""
    names = ["water", "Na+", "Cl-"]
    molality = molalities([0.88, 0.06, 0.06])
    charge = [0.0, 1.0, -1.0]
    context = {
        "common_ion_terms": True,
        "non_two_two_beta2": False,
        "unequal_charge_same_sign": False,
    }
    parameters = CatalogParameters(names)
    for index, expected in ((1, -0.267512432075052), (2, -0.267512432075052)):
        got = phase.ln_gamma(
            molality,
            charge,
            298.15,
            phase.debye_huckel_a_phi(298.15),
            parameters,
            index,
            **context,
        )
        assert got == pytest.approx(expected, abs=1e-10), names[index]


def test_the_same_sign_and_beta2_terms_match_neqsim() -> None:
    """The brine `CaCl2` creates: a same-sign pair of unequal charge, and a non-2:2 `beta2`."""
    names = ["water", "Na+", "Ca++", "Cl-"]
    molality = molalities([0.88, 0.03, 0.03, 0.06])
    charge = [0.0, 1.0, 2.0, -1.0]
    context = {
        "common_ion_terms": True,
        "non_two_two_beta2": True,
        "unequal_charge_same_sign": True,
    }
    parameters = CatalogParameters(names)
    for index, expected in (
        (1, -0.504259131023400),
        (2, -1.84132192677047),
        (3, 0.726599630955942),
    ):
        got = phase.ln_gamma(
            molality,
            charge,
            298.15,
            phase.debye_huckel_a_phi(298.15),
            parameters,
            index,
            **context,
        )
        assert got == pytest.approx(expected, abs=1e-10), names[index]
