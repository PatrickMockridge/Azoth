"""`oil_assay_to_mixture`: an assay's cuts derived, carded and built into a fluid.

The derivation and the card are each tested on their own; this file holds the composition - that
a pseudo cut arrives with its card's constants, a standard substance as a base component, and
that the mixture is one a flash takes.
"""

from __future__ import annotations

import pytest

from azoth import characterization, eos, ureg

q = ureg.Quantity


def pseudo(name: str, fraction: float, molar_mass: float, density: float) -> dict[str, object]:
    return {
        "name": name,
        "declared_fraction": fraction,
        "density": q(density, "kg/m**3"),
        "molar_mass": q(molar_mass, "kg/mol"),
    }


def standard(name: str, substance: str, fraction: float) -> dict[str, object]:
    return {
        "name": name,
        "declared_fraction": fraction,
        "standard_component_name": substance,
    }


def assay() -> list[dict[str, object]]:
    return [
        pseudo("C7", 0.5, 0.095, 730.0),
        pseudo("C10", 0.3, 0.134, 810.0),
        standard("methane", "methane", 0.2),
    ]


def test_an_assay_becomes_the_mixture_its_cuts_describe() -> None:
    mixture, ideal_gas = characterization.oil_assay_to_mixture("mass", assay(), q(100.0, "kg"))

    assert len(mixture.components) == 3
    assert len(ideal_gas.cp_a) == 3
    card = characterization.tbp_fraction_component(q(0.095, "kg/mol"), q(730.0, "kg/m**3"))
    assert mixture.components[0].Tc.to("K").magnitude == pytest.approx(
        card.critical_temperature.to("K").magnitude, rel=1e-12
    )
    assert mixture.components[0].omega == pytest.approx(card.acentric_factor, rel=1e-12)


def test_the_mixture_flashes() -> None:
    mixture, _ = characterization.oil_assay_to_mixture("mass", assay(), q(100.0, "kg"))

    flash = eos.pt_flash(mixture, q(300.0, "K"), q(5.0e6, "Pa"), [0.25, 0.15, 0.6])
    assert flash.z_liquid > 0.0


def test_a_cut_with_no_gravity_is_refused_before_a_card_is_built() -> None:
    cuts = [{"name": "C7", "declared_fraction": 1.0, "molar_mass": q(0.095, "kg/mol")}]
    with pytest.raises(Exception, match="density"):
        characterization.oil_assay_to_mixture("mass", cuts, q(100.0, "kg"))
