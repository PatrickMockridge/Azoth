"""The component databank's formation columns, as the reactions tier reads them.

The three columns are `GIBBSENERGYOFFORMATION`, `ENTHALPYOFFORMATION` and
`ABSOLUTEENTROPY`, and the reference potentials a reactive flash builds need them. They
are read here rather than through :mod:`azoth.eos` because the two tables spell their
names differently - the element table keeps NeqSim's `CO2` and the component databank is
lowercased by its generator - which is the same difference `ionic_charge` settles.

**A zero is a value and not an absence.** The Rust reader in
`crates/azoth-reactions/src/databank.rs` asserts the same two cases, so the twins are
pinned to the same facts rather than to each other; they cannot call each other, and the
id that consumes them is what `test_cross_impl` will compare when it lands.
"""

from __future__ import annotations

from azoth.reactions.reference import _tables


def test_the_formation_properties_are_read_by_name_in_the_databanks_own_spelling() -> None:
    # The element table spells it `CO2`, the component databank spells it `co2`.
    co2 = _tables.formation_properties("CO2")
    assert co2 == _tables.FormationProperties(
        gibbs_energy_of_formation=-394_359.0,
        enthalpy_of_formation=-393_509.0,
        absolute_entropy=213.8,
    )
    assert _tables.formation_properties("no-such-substance") is None


def test_a_zero_formation_property_is_returned_as_a_value() -> None:
    # Oxygen's zeros are the standard state's definition, and `H+` is zero in all three
    # columns by the aqueous convention. Both are answers.
    oxygen = _tables.formation_properties("oxygen")
    assert oxygen is not None
    assert (oxygen.gibbs_energy_of_formation, oxygen.enthalpy_of_formation) == (0.0, 0.0)
    assert oxygen.absolute_entropy == 205.1

    # The same column carries a zero that is *not* an answer: formic acid's fitted
    # enthalpy of formation is a number while its Gibbs energy of formation is zero, and
    # nothing in the value says which kind of zero it is.
    formic_acid = _tables.formation_properties("formic acid")
    assert formic_acid is not None
    assert formic_acid.gibbs_energy_of_formation == 0.0
    assert formic_acid.enthalpy_of_formation == -378_700.0
