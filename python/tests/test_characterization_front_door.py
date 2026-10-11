"""A characterisation cut becomes a fluid, end to end.

This is the claim the card width exists for. Nine ids landed that compute a crude's cut
table - molar masses, normal liquid densities, critical constants, the cubic's attraction
exponent - and before the card could state those numbers, every one of them was a value
with nowhere to go: the only route for a substance the databank does not carry is a
per-substance card, and a card could state `Tc`, `Pc` and `omega`, which a cut cannot be
built from.

So the test is not "the loader stores a number". It is: **the namespace publishes the
numbers, the card carries them, and a flash runs on them** - with the two properties that
make it an arithmetic test rather than a plumbing one:

* the same fluid built by hand from the same published numbers gives the same answer, so
  the card route is the direct route and not a second one; and
* an exponent the card states **changes the answer** when it is read, and the same call
  with Soave's `m` gives a different one - so a reader that dropped `alpha_param_1` and
  fell back cannot satisfy both.

The expected critical constants are the spec's own worked example, from
``validation/neqsim/captures/characterization_probe.tsv`` cut[4] - a `C19`-scale cut at
`M = 500.0` g/mol and `d = 0.88`, which is the same cut ``eos.tbp_fraction_properties``
works. They are asserted here rather than restated, so a drift in the namespace fails this
file rather than being carried into a card.
"""

from __future__ import annotations

import pytest

import azoth
from azoth import eos, keycard
from azoth.core.errors import PropertyUnavailableError
from azoth.core.units import ureg
from azoth.eos.cubic import SRK
from azoth.eos.mixture import Component, mixture

q = ureg.Quantity

#: The cut: a `C19`-scale pseudo-component, which is what the split publishes.
CUT_MOLAR_MASS = 0.5
CUT_DENSITY = 880.0


#: The cut as a card states it: `(value, unit)` per key, each in the unit the schema declares.
def the_cut() -> dict[str, tuple[float, str]]:
    """The card entry for one cut, from the namespace's own answer."""
    cut = azoth.characterization.tbp_cut_properties(
        CUT_MOLAR_MASS * q("kg/mol"), CUT_DENSITY * q("kg/m**3")
    )
    exponent = cut.attraction_exponent
    assert exponent is not None, "the cut states an attraction exponent"
    return {
        "Tc": (cut.tc.to("K").magnitude, "K"),
        "Pc": (cut.pc.to("Pa").magnitude, "Pa"),
        "omega": (cut.acentric_factor, "dimensionless"),
        "molar_mass": (CUT_MOLAR_MASS, "kg/mol"),
        # **The one this whole change exists for.** A cut's attraction exponent is fitted,
        # and before the card could state it a fluid built from a cut ran Soave's `m`.
        "alpha_param_1": (exponent, "dimensionless"),
    }


def a_card(entry: dict[str, tuple[float, str]]) -> keycard.Keycard:
    stated = {name: {"value": value, "unit": unit} for name, (value, unit) in entry.items()}
    return keycard.use({"schema_version": 2, "components": {"c19": stated}})


def test_the_namespace_publishes_the_cut_the_card_needs() -> None:
    """The five numbers, against the capture the namespace is held to.

    Not a restatement of the kernel: these are the spec's worked example, so a change in
    the correlation fails here rather than being carried into a card that still loads.
    """
    cut = azoth.characterization.tbp_cut_properties(
        CUT_MOLAR_MASS * q("kg/mol"), CUT_DENSITY * q("kg/m**3")
    )
    assert cut.tc.to("K").magnitude == pytest.approx(891.945256085426, rel=1e-12)
    assert cut.pc.to("Pa").magnitude == pytest.approx(1196293.59599742, rel=1e-12)
    assert cut.acentric_factor == pytest.approx(1.36586663262266, rel=1e-12)
    # **The one the card could not state before this change**, and the reason for it: a
    # cut's attraction exponent is fitted, not a function of `omega`, and `Alpha::SrkFitted`
    # is the correlation that reads it.
    assert cut.attraction_exponent == pytest.approx(2.228110216, rel=1e-12)


def test_a_cut_becomes_a_component_and_its_mass_survives() -> None:
    """The front door: the card carries the cut's identity into the fluid."""
    fluid = eos.from_names(["c19"], card=a_card(the_cut()), eos="srk", alpha="srk_fitted")
    carried = fluid.components[0].molar_mass
    assert carried is not None, "the card stated a molar mass and the fluid lost it"
    assert carried.to("kg/mol").magnitude == pytest.approx(CUT_MOLAR_MASS)


def test_the_card_route_is_the_direct_route() -> None:
    """**The non-tautological half.** The same numbers, built without a card at all.

    A card that resolved to something other than what it states would agree with itself
    and disagree here: the two fluids are the same substance, one reached through the
    loader and one constructed directly, and a flash has to answer the same for both.
    """
    through_the_card = eos.from_names(
        ["c19"], card=a_card(the_cut()), eos="srk", alpha="srk_fitted"
    )
    by_hand = mixture(
        [
            Component(
                Tc=the_cut()["Tc"][0] * q("K"),
                Pc=the_cut()["Pc"][0] * q("Pa"),
                omega=the_cut()["omega"][0],
                molar_mass=the_cut()["molar_mass"][0] * q("kg/mol"),
                alpha_params=(the_cut()["alpha_param_1"][0],),
            )
        ],
        cubic=SRK,
        alpha="srk_fitted",
    )
    T, P, z = q(400.0, "K"), q(2_000_000.0, "Pa"), [1.0]
    loaded = eos.pt_flash(through_the_card, T, P, z)
    direct = eos.pt_flash(by_hand, T, P, z)
    assert loaded.z_vapour == pytest.approx(direct.z_vapour, rel=1e-12)
    assert loaded.z_liquid == pytest.approx(direct.z_liquid, rel=1e-12)
    assert loaded.vapour_fraction == pytest.approx(direct.vapour_fraction, rel=1e-12)


def test_the_exponent_is_read_rather_than_falling_back_to_soave() -> None:
    """The cut's fitted exponent moves the answer, so a reader cannot be dropping it.

    `srk_fitted` reads `alpha_params[0]`; `srk` is Soave's `m` from `omega` alone. A card
    whose `alpha_param_1` reached the fluid and a fluid whose `Alpha` remained Soave would
    give the *same* answer - so differing answers are what says the number is used, and
    this cannot be satisfied by both paths agreeing.
    """
    T, P, z = q(400.0, "K"), q(2_000_000.0, "Pa"), [1.0]
    fitted = eos.pt_flash(
        eos.from_names(["c19"], card=a_card(the_cut()), eos="srk", alpha="srk_fitted"),
        T,
        P,
        z,
    )
    soave = eos.pt_flash(
        eos.from_names(["c19"], card=a_card(the_cut()), eos="srk", alpha="srk"),
        T,
        P,
        z,
    )
    assert fitted.z_vapour != pytest.approx(soave.z_vapour, rel=1e-9), (
        "a fitted exponent and Soave's give the same Z, so the card's value is not "
        "reaching the alpha function"
    )


def test_a_card_that_adds_a_substance_without_a_mass_is_refused_by_name() -> None:
    """**Both readers, one message.** The refusal the width is organised around.

    A cubic reads only `Tc`, `Pc` and `omega`, so a substance without a molar mass builds
    and a `pt_flash` answers correctly - while every mass-bearing correlation refuses far
    from where the omission was made. The card states the substance, so the card is where
    the omission is caught.
    """
    entry = the_cut()
    entry.pop("molar_mass")
    card = a_card(entry)
    with pytest.raises(PropertyUnavailableError) as raised:
        eos.components.entry("c19", card=card)
    assert "molar_mass" in str(raised.value)

    # And the other reader, through the extension's own resolution of the same card.
    import azoth._core as core

    overlay = azoth._rust_bridge.overlay_from(card)
    with pytest.raises(Exception) as rust_raised:
        core.overlay_entry_row("c19", overlay)
    assert "molar_mass" in str(rust_raised.value)


def test_both_readers_carry_every_new_parameter() -> None:
    """The width, held across the boundary rather than sampled.

    `test_card_agreement.py` compares the two readers' resolution field by field; this is
    the one field per group the front door actually depends on, checked here so a
    regression in the crossing fails a test whose subject is the characterisation chain.
    """
    import azoth._core as core

    entry = the_cut()
    entry["normal_liquid_density"] = (CUT_DENSITY, "kg/m**3")
    entry["parachor"] = (400.0, "dimensionless")
    entry["rackett_z"] = (0.27, "dimensionless")
    card = a_card(entry)
    python_side = eos.components.entry("c19", card=card)

    rust_side = core.overlay_entry_row("c19", azoth._rust_bridge.overlay_from(card))

    assert python_side.molar_mass is not None
    assert python_side.liquid_density is not None
    assert python_side.molar_mass.to("kg/mol").magnitude == pytest.approx(CUT_MOLAR_MASS)
    assert rust_side.molar_mass == pytest.approx(CUT_MOLAR_MASS)
    assert python_side.liquid_density.to("kg/m**3").magnitude == pytest.approx(CUT_DENSITY)
    assert rust_side.normal_liquid_density_kg_per_m3 == pytest.approx(CUT_DENSITY)
    assert python_side.parachor == pytest.approx(400.0)
    assert rust_side.parachor == pytest.approx(400.0)
    assert python_side.rackett_z == pytest.approx(0.27)
    assert rust_side.rackett_z == pytest.approx(0.27)
    # The exponent, which the characterisation chain publishes and the card now carries.
    assert rust_side.alpha_params[0] == pytest.approx(the_cut()["alpha_param_1"][0])
