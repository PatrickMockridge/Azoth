"""The molality-scale surface, as `azoth.eos.reference._electrolyte` carries it.

The Python mirror of the Rust tests in `crates/azoth-eos/src/electrolyte.rs`, asserting
the same constants - the rule a *shared* function follows rather than a registered model.
"""

from __future__ import annotations

import math

import pytest

from azoth.core.errors import InvalidInputError
from azoth.eos.reference import _electrolyte as electrolyte

WATER = 0.018015
SODIUM = 0.02299
CHLORIDE = 0.03545


def test_the_composition_is_per_kilogram_of_water() -> None:
    """**The molality is `n_i / m_solvent`, and the solvent's own entry is not zero.**

    The solvent's molality is its reciprocal molar mass, which is the part of this a
    reader is most likely to assume away - a sum over "the molalities" that included it
    would carry `55.5 mol/kg` of water into an ionic strength it does not belong to.
    """
    c = electrolyte.composition(
        ["water", "na+", "cl-"], [0.9, 0.05, 0.05], [WATER, SODIUM, CHLORIDE], [0.0, 1.0, -1.0]
    )

    assert c.solvent_mass == pytest.approx(0.9 * WATER)
    assert c.molality[1] == pytest.approx(0.05 / (0.9 * WATER))
    assert c.molality[2] == pytest.approx(c.molality[1]), "1:1"
    assert c.molality[0] == pytest.approx(1.0 / WATER), "the solvent's own molality"
    # I = 1/2 (m_na + m_cl) = the solute molality for a 1:1 salt.
    assert c.ionic_strength == pytest.approx(c.molality[1])


def test_the_ionic_strength_weights_by_charge_squared() -> None:
    """A 2:1 salt is not its molality, which is what makes the weighting worth a test."""
    c = electrolyte.composition(
        ["water", "ca++", "cl-", "cl-"],
        [0.9, 0.033, 0.033, 0.033],
        [WATER, 0.04008, CHLORIDE, CHLORIDE],
        [0.0, 2.0, -1.0, -1.0],
    )
    m_calcium = 0.033 / (0.9 * WATER)
    expected = 0.5 * (m_calcium * 4.0 + m_calcium + m_calcium)
    assert c.ionic_strength == pytest.approx(expected)
    assert c.ionic_strength > 1.5 * m_calcium


def test_a_neutral_solute_does_not_reach_the_ionic_strength() -> None:
    """It contributes zero because its charge is zero, with no test of which are ions."""
    c = electrolyte.composition(["water", "co2"], [0.99, 0.01], [WATER, 0.04401], [0.0, 0.0])
    assert c.molality[1] > 0.0
    assert c.ionic_strength == 0.0


def test_a_phase_without_a_solvent_is_refused() -> None:
    """Refused rather than evaluated as a solution of nothing, which is where NeqSim
    returns `I = 0` and a finite answer with no symptom."""
    with pytest.raises(InvalidInputError):
        electrolyte.composition(["methane", "co2"], [0.6, 0.4], [0.016043, 0.04401], [0.0, 0.0])
    with pytest.raises(InvalidInputError):
        electrolyte.composition(["water", "co2"], [0.6, 0.4], [WATER], [0.0, 0.0])


def test_the_water_activity_is_the_osmotic_relation() -> None:
    """`ln a_w = -phi M_w sum m`, against NeqSim's own constants.

    `phi = 1` over pure water's own molality gives `ln a_w = -1` exactly, which is a test
    of the arithmetic rather than of a physical value - the expression is applied to the
    sum the model built, and this pins the factors.
    """
    ln_a_w = electrolyte.ln_water_activity(1.0, 1.0 / WATER, WATER)
    assert ln_a_w == pytest.approx(-1.0, abs=1.0e-12)
    assert math.exp(ln_a_w) == pytest.approx(math.exp(-1.0), abs=1.0e-12)

    gamma = electrolyte.water_activity_coefficient(ln_a_w, 0.9)
    assert gamma == pytest.approx(math.exp(ln_a_w) / 0.9)


def test_a_vanished_water_mole_fraction_is_refused() -> None:
    """Where NeqSim returns the ideal `gamma = 1` rather than the value this state has."""
    with pytest.raises(InvalidInputError):
        electrolyte.water_activity_coefficient(-1.0, 0.0)
