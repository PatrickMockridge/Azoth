"""Pitzer's activity-coefficient arithmetic, as `_pitzer_phase` carries it.

The Python mirror of the Rust tests in `crates/azoth-eos/src/pitzer_phase.rs`, asserting
the same constants - the rule a *shared* function follows rather than a registered model.

The end-to-end numbers come from `validation/neqsim/PitzerArithmetic.java` on a live
`SystemPitzer` at 298.15 K, so they are NeqSim's rather than a restatement of its formulas.
"""

from __future__ import annotations

import math

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


def test_the_water_activity_coefficient_matches_neqsim() -> None:
    """**The solvent route is not the ion one.**

    `getWaterGamma` builds the solvent's coefficient from the Pitzer *osmotic* coefficient,
    whose binary function is `beta0 + beta1 exp(-alpha sqrt(I))` - where the ion branch's is
    `beta0 + beta1 g(alpha sqrt(I))`. A port that shared the two would agree on neither
    number, which is what these two cases pin.
    """
    for names, x, charges, expected in (
        (["water", "Na+", "Cl-"], [0.88, 0.06, 0.06], [0.0, 1.0, -1.0], -0.0220339385649971),
        (
            ["water", "Na+", "Ca++", "Cl-"],
            [0.88, 0.03, 0.03, 0.06],
            [0.0, 1.0, 2.0, -1.0],
            -0.0563151245495846,
        ),
    ):
        four = len(charges) == 4
        got = phase.ln_gamma_water(
            molalities(x),
            charges,
            298.15,
            phase.debye_huckel_a_phi(298.15),
            CatalogParameters(names),
            0,
            0.88,
            0.0,
            non_two_two_beta2=four,
            unequal_charge_same_sign=four,
            neutral_interactions_active=False,
        )
        assert got == pytest.approx(expected, abs=1e-10), names


def test_the_osmotic_coefficient_matches_neqsim() -> None:
    """Inverted from `ln a_w = -phi M_w sum m`, so the number checked is the model's own."""
    for names, x, charges, expected in (
        (["water", "Na+", "Cl-"], [0.88, 0.06, 0.06], [0.0, 1.0, -1.0], 1.09902694054914),
        (
            ["water", "Na+", "Ca++", "Cl-"],
            [0.88, 0.03, 0.03, 0.06],
            [0.0, 1.0, 2.0, -1.0],
            1.35042230443611,
        ),
    ):
        four = len(charges) == 4
        molality = molalities(x)
        ln_a_w = phase.ln_gamma_water(
            molality,
            charges,
            298.15,
            phase.debye_huckel_a_phi(298.15),
            CatalogParameters(names),
            0,
            0.88,
            0.0,
            non_two_two_beta2=four,
            unequal_charge_same_sign=four,
            neutral_interactions_active=False,
        ) + math.log(0.88)
        sum_m = sum(m for m, z in zip(molality, charges, strict=False) if z != 0.0)
        phi = -ln_a_w / (phase.WATER_MOLAR_MASS * sum_m)
        assert phi == pytest.approx(expected, abs=1e-9), names


def test_the_neutral_layer_matches_neqsim() -> None:
    """**The λ/ζ/μ/η layer, on the one topology the catalogue covers for it.**

    The oracle is `PitzerArithmetic.java`'s CO2-brine section. **The chloride case falls
    back to the legacy dataset**, because the catalogue has no `ZETA(CO2, Na+, Cl-)` row -
    it pairs CO2 and H2S with *sulphate* - so the layer is reachable only through a
    sulphate-bearing brine, and a CO2/NaCl brine gets no neutral physics at all.
    """
    names = ["water", "Na+", "SO4--", "CO2"]
    charge = [0.0, 1.0, -2.0, 0.0]
    molality = molalities([0.86, 0.06, 0.03, 0.05])

    interactions = phase.catalogue_interactions(names, charge, [1, 2], [3])
    assert interactions is not None, "the catalogue covers CO2 with Na+ and SO4--"
    assert len(interactions) == 4, "one LAMBDA(CO2,CO2), two LAMBDA(CO2,ion), one ZETA"

    assert phase.osmotic_neutral(interactions, molality, 298.15) == pytest.approx(
        1.09825174371, abs=1e-9
    )
    for component, expected in (
        (1, 0.454900106480425),
        (2, 0.296616109274644),
        (3, 0.749844524371078),
    ):
        got = phase.ln_gamma_neutral(interactions, molality, component, 298.15)
        assert got == pytest.approx(expected, abs=1e-9), names[component]


def test_the_repetition_structure_decides_the_coefficients() -> None:
    """The repeated-species case is the one that is not `[2, 2]`.

    PHREEQC differentiates both slots before accumulating them, so the same component twice
    pairs with one each and a **half** - not two each and one.
    """
    form = (1.0, 0.0, 0.0, 0.0, 0.0, 0.0)

    distinct = phase.NeutralInteraction(phase.LAMBDA, [0, 1], form)
    assert distinct.osmotic_contribution([2.0, 3.0], 298.15) == 6.0

    repeated = phase.NeutralInteraction(phase.LAMBDA, [0, 0], form)
    assert repeated.osmotic_contribution([2.0, 3.0], 298.15) == 2.0
    assert repeated.log_gamma_contribution([2.0, 3.0], 0, 298.15) == 4.0

    # `Mu`'s multiplicity counts the tuple's distinct permutations: 1, 3 and 6.
    assert (
        phase.NeutralInteraction(phase.MU, [0, 0, 0], form).osmotic_contribution(
            [2.0, 0.0, 0.0], 298.15
        )
        == 8.0
    )
    assert (
        phase.NeutralInteraction(phase.MU, [0, 0, 1], form).osmotic_contribution(
            [2.0, 3.0, 0.0], 298.15
        )
        == 36.0
    )
    assert (
        phase.NeutralInteraction(phase.MU, [0, 1, 2], form).osmotic_contribution(
            [2.0, 3.0, 4.0], 298.15
        )
        == 144.0
    )


def test_an_uncovered_zeta_row_abandons_the_layer() -> None:
    """The coverage rule rather than a defect, and the reason a CO2/NaCl brine loses the
    whole family."""
    assert catalog.find("ZETA", ["CO2", "Na+", "Cl-"]) is None
    assert (
        phase.catalogue_interactions(
            ["water", "Na+", "Cl-", "CO2"], [0.0, 1.0, -1.0, 0.0], [1, 2], [3]
        )
        is None
    )
    assert (
        phase.catalogue_interactions(
            ["water", "Na+", "SO4--", "CO2"], [0.0, 1.0, -2.0, 0.0], [1, 2], [3]
        )
        is not None
    )
