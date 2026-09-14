"""The keycard loader: what it overrides, and what it refuses.

A keycard is the extension surface, and the failure it has to be organised against is
the one this repository names everywhere else - **a declaration that is stored and
read by nothing**. A model that names a cubic the code does not implement, a component
parameter no implementation reads, a unit the checkpoint accepts and the loader reads
as something else: each of those is a keycard that loads, looks like it is doing
something, and changes no answer. So most of this file asserts *refusals*, and the
positive cases are there to prove the refusals are about the rule under test rather
than about a loader that rejects everything.

Two tests hold the loader's vocabularies to the schema's. The loader cannot read the
schema at runtime - it is a dev-time artefact and is not shipped in the wheel - so the
vocabularies are a copy, and a copy without a test is the drift this project answers
with two implementations compared against each other. The units are the sharpest case:
`pint` parses `kelvin`, so a loader that trusted `pint` would accept a keycard that
`tools/check_user_data.py` then rejects. `UNIT_VOCABULARY` is generated from
`specs/vocabulary/vocabulary.toml` by `tools/gen_vocabulary.py` and ships with the
package, which is why it can be read at runtime at all; the model vocabularies below
are still hand-written beside this module, and are the ones the comparison protects.
"""

from __future__ import annotations

import inspect
import json
import tomllib
from pathlib import Path
from typing import Any

import pytest

from azoth import eos, hydraulics, keycard
from azoth.core.errors import InvalidInputError, KeycardError, PropertyUnavailableError
from azoth.core.units import ureg
from azoth.eos.mixture import Component

REPO_ROOT = Path(__file__).resolve().parents[2]
SCHEMA = REPO_ROOT / "specs" / "schema" / "keycard.schema.json"
CALC_SCHEMA = REPO_ROOT / "specs" / "schema" / "calc.schema.json"
#: The unit enum lives in its own schema, which `calc.schema.json` `$ref`s. It is
#: generated from `specs/vocabulary/vocabulary.toml` by `tools/gen_vocabulary.py`.
UNIT_SCHEMA = REPO_ROOT / "specs" / "schema" / "unit.schema.json"
TEMPLATE = REPO_ROOT / "keycard.example.toml"

q = ureg.Quantity


def minimal(**sections: Any) -> dict[str, Any]:
    """The smallest keycard carrying the given sections."""
    return {"schema_version": 2, **sections}


def a_card(**sections: Any) -> keycard.Keycard:
    """A keycard for a call to read, built from the given sections.

    A helper rather than a fixture, and there is no autouse fixture here at all: the
    card is a value, so a test that wants one holds it and passes it. There is
    nothing a test could leave behind for the next test to inherit, which is what the
    fixture this replaced existed to clean up.
    """
    return keycard.use(minimal(**sections))


# ---------------------------------------------------------------------------
# The vocabularies, held to the schema
# ---------------------------------------------------------------------------


def test_the_unit_vocabulary_is_the_schema_s_enum() -> None:
    """Compared in both directions, because either alone is satisfiable by accident.

    A loader offering a unit the schema forbids accepts a file the checker rejects; a
    loader missing a unit the schema allows rejects a file the checker accepts. Each
    passes a one-way test.
    """
    schema = json.loads(UNIT_SCHEMA.read_text(encoding="utf-8"))
    declared = tuple(schema["enum"])
    assert declared == keycard.UNIT_VOCABULARY, (
        "the loader's unit vocabulary and the unit schema's enum have drifted. One of "
        "them accepts what the other refuses."
    )


@pytest.mark.parametrize(
    ("name", "schema_path"),
    [
        ("MODEL_KINDS", "kind"),
        ("MODEL_SHAPES", "shape"),
        ("MODEL_ALPHAS", "alpha"),
        ("MODEL_MIXING_RULES", "mixing_rule"),
    ],
)
def test_the_model_vocabularies_are_the_schema_s_enums(name: str, schema_path: str) -> None:
    """The same both-ways comparison, for the model vocabularies.

    This is the test that would have caught the state this file was written in: the
    schema listed five cubic shapes and six alpha functions, the code implemented one,
    and a keycard declaring `soave_redlich_kwong` would have loaded, been stored, and
    been evaluated as Peng-Robinson.
    """
    schema = json.loads(SCHEMA.read_text(encoding="utf-8"))
    declared = schema["properties"]["models"]["additionalProperties"]["properties"]
    assert tuple(getattr(keycard, name)) == tuple(declared[schema_path]["enum"])


def test_the_model_keys_are_the_schema_s_properties() -> None:
    """Every key a definition may carry, compared in both directions.

    `additionalProperties: false` is what makes this bite. A key the schema admits and
    the loader refuses rejects a file the checker accepts; a key the loader admits and
    the schema refuses is stored nowhere, which is the failure `critical_rule`,
    `volume_translation` and `root_selection` were: accepted, never read, and gone from
    the format now.
    """
    schema = json.loads(SCHEMA.read_text(encoding="utf-8"))
    model = schema["properties"]["models"]["additionalProperties"]
    declared = set(model["properties"]) | set(model["required"])
    assert declared == keycard.MODEL_KEYS, (
        "the loader's accepted model keys and keycard.schema.json's properties have "
        "drifted. One of them accepts what the other refuses."
    )


def test_the_component_parameters_are_what_the_implementation_reads() -> None:
    """Every parameter the keycard accepts must reach a `Component`.

    A parameter that the loader accepts and `Component` has no field for is a value
    silently dropped - the exact failure the closed list exists to prevent.
    """
    assert set(keycard.COMPONENT_PARAMETERS) == {"Tc", "Pc", "omega"}
    fields = set(Component.__dataclass_fields__)
    assert set(keycard.COMPONENT_PARAMETERS) <= fields


# ---------------------------------------------------------------------------
# Components
# ---------------------------------------------------------------------------


def test_a_component_override_wins_over_the_databank() -> None:
    """The point of the section: a user's value replaces the shipped one."""
    shipped = eos.components.entry("methane").Pc
    card = a_card(components={"methane": {"Pc": {"value": 4_600_000.0, "unit": "Pa"}}})
    overridden = eos.components.entry("methane", card=card)
    assert overridden.Pc.to("Pa").magnitude == pytest.approx(4_600_000.0)
    assert overridden.Pc != shipped
    assert overridden.source == "keycard"


def test_an_override_keeps_the_parameters_it_does_not_name() -> None:
    """Parameter by parameter, not record by record.

    A user correcting one value should not have to restate the others - and, more to
    the point, should not silently lose them if they do not.
    """
    shipped = eos.components.entry("methane")
    card = a_card(components={"methane": {"omega": {"value": 0.5, "unit": "dimensionless"}}})
    now = eos.components.entry("methane", card=card)
    assert now.omega == pytest.approx(0.5)
    assert now.Tc == shipped.Tc
    assert now.Pc == shipped.Pc


def test_a_component_is_converted_not_assumed() -> None:
    """A `Tc` in the wrong unit is a shift with no symptom, so it is converted.

    This is the whole reason a keycard declares units rather than bare numbers.
    """
    card = a_card(
        components={
            "methane": {
                "Tc": {"value": 190.0, "unit": "K"},
                "Pc": {"value": 4.599e6, "unit": "Pa"},
                "omega": {"value": 0.0115, "unit": "dimensionless"},
            }
        }
    )
    resolved = eos.components.entry("methane", card=card)
    assert resolved.Tc.to("K").magnitude == pytest.approx(190.0)
    assert resolved.Pc.to("Pa").magnitude == pytest.approx(4.599e6)


def test_a_component_the_databank_lacks_is_added() -> None:
    """An added substance has no CAS number, and `None` says so.

    Filling those fields from a similar substance would be inventing data, which is
    what this module is arranged against.
    """
    card = a_card(
        components={
            "unobtainium": {
                "Tc": {"value": 500.0, "unit": "K"},
                "Pc": {"value": 2.0e6, "unit": "Pa"},
                "omega": {"value": 0.3, "unit": "dimensionless"},
            }
        }
    )
    assert "unobtainium" in eos.components.available(card=card)
    added = eos.components.entry("unobtainium", card=card)
    assert added.cas is None
    assert added.molar_mass is None
    assert added.source == "keycard"


def test_a_partial_new_component_is_refused() -> None:
    """Missing a parameter a cubic reads is an error, not a default."""
    card = a_card(components={"unobtainium": {"Tc": {"value": 500.0, "unit": "K"}}})
    with pytest.raises(PropertyUnavailableError, match="is missing"):
        eos.components.entry("unobtainium", card=card)


def test_a_parameter_nothing_reads_is_refused() -> None:
    """A stored-but-unread value is the failure this whole file is about."""
    with pytest.raises(KeycardError, match="not a parameter this build reads"):
        keycard.use(
            minimal(components={"methane": {"critical_pressure": {"value": 1.0, "unit": "Pa"}}})
        )


def test_an_unknown_unit_is_refused() -> None:
    """`kelvin` is the spelling people reach for, and the vocabulary says `K`."""
    with pytest.raises(KeycardError, match="not in the vocabulary"):
        keycard.use(minimal(components={"methane": {"Tc": {"value": 1.0, "unit": "kelvin"}}}))


# ---------------------------------------------------------------------------
# Interaction parameters
# ---------------------------------------------------------------------------


def test_a_keycard_kij_wins_over_the_databank() -> None:
    """The parameters a mixing rule actually reads."""
    databank_value = eos.components.kij_for(("methane", "n-butane"))[(0, 1)]
    card = a_card(kij=[{"component_a": "methane", "component_b": "n-butane", "value": 0.5}])
    assert eos.components.kij_for(("methane", "n-butane"), card=card)[(0, 1)] == pytest.approx(0.5)
    assert databank_value != pytest.approx(0.5)


def test_a_kij_pair_is_order_independent() -> None:
    """`a`/`b` and `b`/`a` name one pair, so they are stored as one."""
    card = keycard.use(
        minimal(kij=[{"component_a": "n-butane", "component_b": "methane", "value": 0.02}])
    )
    assert card.kij_for("methane", "n-butane") == pytest.approx(0.02)
    assert card.kij_for("n-butane", "methane") == pytest.approx(0.02)


def test_a_component_paired_with_itself_is_refused() -> None:
    """A self-pair would silently rescale that component's attraction."""
    with pytest.raises(KeycardError, match="does not interact with itself"):
        keycard.use(
            minimal(kij=[{"component_a": "methane", "component_b": "methane", "value": 0.1}])
        )


# ---------------------------------------------------------------------------
# Coefficients
# ---------------------------------------------------------------------------


def test_an_explicit_coefficient_wins_and_the_keycard_is_not_consulted() -> None:
    """The precedence the schema states, which is the direction that matters.

    A caller who writes a value in their own call is stating it for that call.
    Overriding it from a file they loaded an hour ago would be the worst kind of
    surprise, so the keycard is not merely outranked here - it is not read.
    """
    card = a_card(
        coefficients={"hydraulics.orifice_flow": {"Cd": {"value": 0.1, "unit": "dimensionless"}}}
    )
    explicit = hydraulics.orifice_flow(
        q(50, "mm"), q(20, "kPa"), q(998, "kg/m**3"), 0.61, card=card
    )
    from_card = hydraulics.orifice_flow(q(50, "mm"), q(20, "kPa"), q(998, "kg/m**3"), card=card)
    assert explicit.q != from_card.q


def test_a_coefficient_the_keycard_supplies_is_used_when_omitted() -> None:
    """The positive half: omitting it picks up the keycard's value."""
    card = a_card(
        coefficients={"hydraulics.orifice_flow": {"Cd": {"value": 0.61, "unit": "dimensionless"}}}
    )
    omitted = hydraulics.orifice_flow(q(50, "mm"), q(20, "kPa"), q(998, "kg/m**3"), card=card)
    explicit = hydraulics.orifice_flow(q(50, "mm"), q(20, "kPa"), q(998, "kg/m**3"), 0.61)
    assert omitted.q == explicit.q


def test_a_coefficient_that_is_omitted_with_no_keycard_is_an_error() -> None:
    """A plausible default nobody chose is a wrong answer with no symptom."""
    with pytest.raises(InvalidInputError, match="none was passed"):
        hydraulics.orifice_flow(q(50, "mm"), q(20, "kPa"), q(998, "kg/m**3"))


def test_a_convention_is_carried_with_the_coefficient() -> None:
    """Where a quantity is not a true ratio, the convention has to travel with it."""
    card = keycard.use(
        minimal(
            coefficients={
                "hydraulics.orifice_flow": {
                    "Cd": {
                        "value": 0.61,
                        "unit": "dimensionless",
                        "convention": "iso_5167_corner",
                    }
                }
            }
        )
    )
    assert card.conventions["hydraulics.orifice_flow"]["Cd"] == "iso_5167_corner"


# ---------------------------------------------------------------------------
# Models
# ---------------------------------------------------------------------------


def a_model(**overrides: Any) -> dict[str, Any]:
    """A model definition that is valid unless a test says otherwise."""
    body: dict[str, Any] = {
        "kind": "cubic_eos",
        "shape": "peng_robinson",
        "alpha": "peng_robinson",
        "mixing_rule": "classical_kij",
        "components": ["methane", "n-butane"],
    }
    body.update(overrides)
    return body


def test_a_model_resolves_to_a_mixture() -> None:
    """The positive half, and the reason the section exists: name a mixture once."""
    card = a_card(models={"vendor_gas": a_model()})
    fluid = eos.from_model("vendor_gas", card=card)
    assert len(fluid) == 2
    assert fluid.kij[0][1] == pytest.approx(eos.components.kij_for(("methane", "n-butane"))[(0, 1)])


def test_a_model_sees_the_keycard_s_component_overrides() -> None:
    """The sections compose, which is the point of resolving components by name."""
    card = a_card(
        components={"methane": {"Tc": {"value": 123.0, "unit": "K"}}},
        models={"vendor_gas": a_model()},
    )
    resolved = eos.from_model("vendor_gas", card=card)
    assert resolved.components[0].Tc.to("K").magnitude == pytest.approx(123.0)


@pytest.mark.parametrize("field", ["kind", "shape", "alpha", "mixing_rule"])
def test_a_model_name_this_build_does_not_implement_is_refused(field: str) -> None:
    """Refused when the keycard loads, not when the model is used.

    The alternative - accepting the declaration and evaluating something else - is a
    wrong answer with no symptom, which is the failure mode this project treats as
    its worst.
    """
    with pytest.raises(KeycardError, match="does not implement"):
        keycard.use(minimal(models={"vendor_gas": a_model(**{field: "not_a_thing"})}))


def test_a_fitted_alpha_parameter_is_refused() -> None:
    """The alpha this build runs takes no parameters, so supplying one is a mismatch."""
    with pytest.raises(KeycardError, match="takes none"):
        keycard.use(minimal(models={"vendor_gas": a_model(alpha_parameters={"eta": 1.0})}))


def test_an_undeclared_model_names_what_is_declared() -> None:
    """The error has to say what exists, or it is a puzzle rather than a message."""
    card = a_card(models={"vendor_gas": a_model()})
    with pytest.raises(PropertyUnavailableError, match="vendor_gas"):
        eos.from_model("something_else", card=card)


# ---------------------------------------------------------------------------
# The document itself
# ---------------------------------------------------------------------------


def test_there_is_no_card_in_force() -> None:
    """The module holds no card, and offers no way to set one.

    This is the property, asserted rather than assumed. `load` and `use` return what
    they read and store nothing, so there is no global for a second call to inherit -
    and the `clear()` that used to undo one is gone with it, because a function whose
    only job is to unset a global exists only because the global does.

    Checked by *building a card and looking* rather than by asserting that `current`
    raises. Two weaker tests are passed by defects this one catches:
    asserting only for the missing accessor passes with a `_current` still being
    written and never read; and looking for a bound card *without* building one first
    passes with a `_current` declared and left `None`, which is the same dead binding
    one step earlier. So a card is built here, and then the module is inspected.
    """
    keycard.use(minimal(components={"methane": {"omega": {"value": 0.5, "unit": "dimensionless"}}}))

    live = [name for name, value in vars(keycard).items() if isinstance(value, keycard.Keycard)]
    assert not live, (
        f"azoth.keycard holds {live} at module scope. A card in force is a card whose "
        f"answer depends on what somebody loaded before the call."
    )
    assert not hasattr(keycard, "current"), "azoth.keycard still offers `current`"
    assert not hasattr(keycard, "clear"), "azoth.keycard still offers `clear`"


def test_a_call_with_no_card_reads_what_the_library_ships() -> None:
    """The baseline is the databank, and it is what a cardless call reads.

    Both directions, again: a card changes the answer, and without one the answer is
    the shipped one. A test that only checked the first would pass if a card leaked
    into a later call.
    """
    shipped = eos.components.entry("methane").Pc
    card = a_card(components={"methane": {"Pc": {"value": 1.0, "unit": "Pa"}}})
    assert eos.components.entry("methane", card=card).Pc != shipped
    assert eos.components.entry("methane").Pc == shipped
    assert eos.components.entry("methane").source == "databank"


def test_an_unknown_section_is_refused() -> None:
    """A misspelled section is data that looks in use and is read by nothing."""
    with pytest.raises(KeycardError, match="does not define"):
        keycard.use(minimal(component={"methane": {}}))


@pytest.mark.parametrize("version", [1, 3, None, "2"])
def test_a_version_this_build_does_not_read_is_refused(version: Any) -> None:
    """Version 1 is refused rather than upgraded: it had no `components`."""
    document = {"schema_version": version, "fittings": []}
    with pytest.raises(KeycardError, match="schema_version"):
        keycard.use(document)


def test_a_document_that_is_not_a_mapping_is_refused() -> None:
    """A list of rows is a plausible mistake and an unhelpful one to guess at."""
    with pytest.raises(KeycardError, match="must be a mapping"):
        keycard.use(["methane"])


def test_the_template_is_a_keycard_this_build_loads() -> None:
    """`keycard.example.toml` must pass its own loader, not just its checker.

    Same argument as the template passing `check_user_data.py`: a template that does
    not load teaches the wrong shape, and it is the file a user copies.
    """
    card = keycard.load(TEMPLATE)
    assert card.keyholder == "Example Engineering Ltd"
    assert card.component("methane") is not None
    assert card.kij_for("methane", "n-butane") is not None
    assert card.coefficient("hydraulics.orifice_flow", "Cd") is not None


def test_the_template_parses_as_toml_with_a_version() -> None:
    """A guard against the file being replaced by prose during an edit.

    Parsed rather than read as text, so a template that stopped being a TOML document
    while still containing the right words fails here.
    """
    document = tomllib.loads(TEMPLATE.read_text(encoding="utf-8"))
    assert document["schema_version"] == keycard.SCHEMA_VERSION


def test_a_file_that_is_not_toml_is_refused_as_such(tmp_path: Path) -> None:
    """The parse failure names the format it tried.

    The message is the whole of what a user gets when their file does not parse, so it
    has to name what was expected. The fixture is a keycard written in the format this
    repository read before TOML, which is the document most likely to be lying around.
    """
    stale = tmp_path / "keycard.toml"
    stale.write_text("schema_version: 2\n", encoding="utf-8")

    with pytest.raises(KeycardError, match="does not parse as TOML"):
        keycard.load(stale)


# ---------------------------------------------------------------------------
# Which card a call reads
# ---------------------------------------------------------------------------


def test_two_cards_in_one_process_are_two_answers() -> None:
    """A card is a value, so two datasets are two calls rather than two processes.

    This is the property `spec.md` states as a rule and that a process-wide card
    cannot have: *"a library whose answers depend on call order is a library that
    returns two results for one calculation."* With the card as an argument the two
    results are asked for explicitly and in one process, and the second call does not
    depend on the first having happened.

    Both directions are asserted, and a test that checked only one would pass with the
    card ignored: the two cards must give *different* answers, and each must give its
    own answer whatever order they are asked in.
    """
    lean = a_card(components={"methane": {"omega": {"value": 0.5, "unit": "dimensionless"}}})
    rich = a_card(components={"methane": {"omega": {"value": 0.25, "unit": "dimensionless"}}})

    assert eos.components.entry("methane", card=lean).omega == pytest.approx(0.5)
    assert eos.components.entry("methane", card=rich).omega == pytest.approx(0.25)
    # Reversed, to show the first call did not leave anything behind for the second.
    assert eos.components.entry("methane", card=rich).omega == pytest.approx(0.25)
    assert eos.components.entry("methane", card=lean).omega == pytest.approx(0.5)


def test_the_card_reaches_the_coefficient_path_too() -> None:
    """A coefficient, not just a component, comes from the card the call passed.

    The two read different sections of the same object, so this is a check that the
    parameter is threaded to both rather than that either works. `orifice_flow` is the
    calc with a coefficient a card can supply.
    """
    from azoth.core.units import quantity

    def a_cd(value: float) -> keycard.Keycard:
        return a_card(
            coefficients={
                "hydraulics.orifice_flow": {"Cd": {"value": value, "unit": "dimensionless"}}
            }
        )

    at_60 = a_cd(0.60)
    at_61 = a_cd(0.61)

    # Through the public calc, so the card is threaded the whole way. `Cd` scales the
    # flow linearly, so the ratio between two cards is the ratio of their coefficients
    # and nothing else - which makes this a statement about *which* card answered
    # rather than about the arithmetic.
    bore = quantity(50.0, "mm")
    drop = quantity(10_000.0, "Pa")
    density = quantity(998.0, "kg/m**3")
    sixty = hydraulics.orifice_flow(bore, drop, density, card=at_60).q.to("m**3/s").magnitude
    sixty_one = hydraulics.orifice_flow(bore, drop, density, card=at_61).q.to("m**3/s").magnitude
    assert sixty_one == pytest.approx(sixty * 0.61 / 0.60)

    # And with no card, no coefficient is supplied: an error rather than a default.
    with pytest.raises(InvalidInputError, match="needs a value"):
        hydraulics.orifice_flow(bore, drop, density)


def test_a_kernel_does_not_take_a_card() -> None:
    """The reference implementations take spec inputs only, and no capability.

    A kernel that took a card would be reading the caller's authority itself rather
    than the values it was handed. The card is resolved at the boundary, and
    `test_registry_contract.py` asserts the same asymmetry across every calc.
    """
    from azoth.hydraulics.reference import orifice_flow as reference

    assert "card" not in inspect.signature(reference).parameters


def test_a_card_reaches_the_physics() -> None:
    """Not just a lookup: a card changes what a calculation answers.

    Everything else in this file is about the loader and the data layer. This is the
    claim the mechanism exists for - **that the values a card supplies are the values
    a calculation runs on** - and it is checkable end to end, because a calculation
    names its components and resolves them through `mixture_of`.

    A `Tc` shifted by a hundred kelvin moves a flash's vapour fraction by a lot, so
    the assertion below is about arithmetic rather than about plumbing. Both
    directions: the card's answer differs from the shipped one, and the shipped
    answer is unchanged by a card the call was not given.
    """
    names = ["methane", "n-butane"]
    feed = [0.6, 0.4]
    T = q(300.0, "K")
    P = q(2_000_000.0, "Pa")

    shipped_mixture, shipped_gas = eos.components.mixture_of(names)
    baseline = eos.pt_flash(shipped_mixture, T, P, feed)

    card = a_card(components={"methane": {"Tc": {"value": 300.0, "unit": "K"}}})
    shifted_mixture, shifted_gas = eos.components.mixture_of(names, card=card)
    shifted = eos.pt_flash(shifted_mixture, T, P, feed)

    assert shifted.beta != pytest.approx(baseline.beta), (
        "a card that shifts methane's critical temperature by 110 K left the vapour "
        "fraction where it was, so the card is not reaching the calculation"
    )

    # The card is not sticky: the same call without it answers as it did before.
    again = eos.pt_flash(eos.components.mixture_of(names)[0], T, P, feed)
    assert again.beta == pytest.approx(baseline.beta)

    # And the ideal-gas model came from the same card, so an enthalpy would move too.
    assert shifted_gas.cp_a == shipped_gas.cp_a
