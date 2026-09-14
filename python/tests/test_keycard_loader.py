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
`specs/vocabulary/vocabulary.yaml` by `tools/gen_vocabulary.py` and ships with the
package, which is why it can be read at runtime at all; the model vocabularies below
are still hand-written beside this module, and are the ones the comparison protects.
"""

from __future__ import annotations

import inspect
import json
from pathlib import Path
from typing import Any

import pytest
import yaml

from azoth import eos, hydraulics, keycard
from azoth.core.errors import InvalidInputError, KeycardError, PropertyUnavailableError
from azoth.core.units import ureg
from azoth.eos.mixture import Component

REPO_ROOT = Path(__file__).resolve().parents[2]
SCHEMA = REPO_ROOT / "specs" / "schema" / "keycard.schema.json"
CALC_SCHEMA = REPO_ROOT / "specs" / "schema" / "calc.schema.json"
#: The unit enum lives in its own schema, which `calc.schema.json` `$ref`s. It is
#: generated from `specs/vocabulary/vocabulary.yaml` by `tools/gen_vocabulary.py`.
UNIT_SCHEMA = REPO_ROOT / "specs" / "schema" / "unit.schema.json"
TEMPLATE = REPO_ROOT / "keycard.example.yaml"

q = ureg.Quantity


@pytest.fixture(autouse=True)
def _no_keycard_left_loaded() -> Any:
    """Every test starts and ends with the library on the data it ships.

    Not a convenience: `load()` sets a process-wide keycard, so a test that leaves one
    behind changes the answer every later test gets, and the failure appears in
    whichever test happens to run next.
    """
    keycard.clear()
    yield
    keycard.clear()


def minimal(**sections: Any) -> dict[str, Any]:
    """The smallest keycard carrying the given sections."""
    return {"schema_version": 2, **sections}


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
    keycard.use(minimal(components={"methane": {"Pc": {"value": 4_600_000.0, "unit": "Pa"}}}))
    overridden = eos.components.entry("methane")
    assert overridden.Pc.to("Pa").magnitude == pytest.approx(4_600_000.0)
    assert overridden.Pc != shipped
    assert overridden.source == "keycard"


def test_an_override_keeps_the_parameters_it_does_not_name() -> None:
    """Parameter by parameter, not record by record.

    A user correcting one value should not have to restate the others - and, more to
    the point, should not silently lose them if they do not.
    """
    shipped = eos.components.entry("methane")
    keycard.use(minimal(components={"methane": {"omega": {"value": 0.5, "unit": "dimensionless"}}}))
    now = eos.components.entry("methane")
    assert now.omega == pytest.approx(0.5)
    assert now.Tc == shipped.Tc
    assert now.Pc == shipped.Pc


def test_a_component_is_converted_not_assumed() -> None:
    """A `Tc` in the wrong unit is a shift with no symptom, so it is converted.

    This is the whole reason a keycard declares units rather than bare numbers.
    """
    keycard.use(
        minimal(
            components={
                "methane": {
                    "Tc": {"value": 190.0, "unit": "K"},
                    "Pc": {"value": 4.599e6, "unit": "Pa"},
                    "omega": {"value": 0.0115, "unit": "dimensionless"},
                }
            }
        )
    )
    resolved = eos.components.entry("methane")
    assert resolved.Tc.to("K").magnitude == pytest.approx(190.0)
    assert resolved.Pc.to("Pa").magnitude == pytest.approx(4.599e6)


def test_a_component_the_databank_lacks_is_added() -> None:
    """An added substance has no CAS number, and `None` says so.

    Filling those fields from a similar substance would be inventing data, which is
    what this module is arranged against.
    """
    keycard.use(
        minimal(
            components={
                "unobtainium": {
                    "Tc": {"value": 500.0, "unit": "K"},
                    "Pc": {"value": 2.0e6, "unit": "Pa"},
                    "omega": {"value": 0.3, "unit": "dimensionless"},
                }
            }
        )
    )
    assert "unobtainium" in eos.components.available()
    added = eos.components.entry("unobtainium")
    assert added.cas is None
    assert added.molar_mass is None
    assert added.source == "keycard"


def test_a_partial_new_component_is_refused() -> None:
    """Missing a parameter a cubic reads is an error, not a default."""
    keycard.use(minimal(components={"unobtainium": {"Tc": {"value": 500.0, "unit": "K"}}}))
    with pytest.raises(PropertyUnavailableError, match="is missing"):
        eos.components.entry("unobtainium")


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
    keycard.use(minimal(kij=[{"component_a": "methane", "component_b": "n-butane", "value": 0.5}]))
    assert eos.components.kij_for(("methane", "n-butane"))[(0, 1)] == pytest.approx(0.5)
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
    keycard.use(
        minimal(
            coefficients={
                "hydraulics.orifice_flow": {"Cd": {"value": 0.1, "unit": "dimensionless"}}
            }
        )
    )
    explicit = hydraulics.orifice_flow(q(50, "mm"), q(20, "kPa"), q(998, "kg/m**3"), 0.61)
    from_keycard = hydraulics.orifice_flow(q(50, "mm"), q(20, "kPa"), q(998, "kg/m**3"))
    assert explicit.q != from_keycard.q


def test_a_coefficient_the_keycard_supplies_is_used_when_omitted() -> None:
    """The positive half: omitting it picks up the keycard's value."""
    keycard.use(
        minimal(
            coefficients={
                "hydraulics.orifice_flow": {"Cd": {"value": 0.61, "unit": "dimensionless"}}
            }
        )
    )
    omitted = hydraulics.orifice_flow(q(50, "mm"), q(20, "kPa"), q(998, "kg/m**3"))
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
    keycard.use(minimal(models={"vendor_gas": a_model()}))
    fluid = eos.from_model("vendor_gas")
    assert len(fluid) == 2
    assert fluid.kij[0][1] == pytest.approx(eos.components.kij_for(("methane", "n-butane"))[(0, 1)])


def test_a_model_sees_the_keycard_s_component_overrides() -> None:
    """The sections compose, which is the point of resolving components by name."""
    keycard.use(
        minimal(
            components={"methane": {"Tc": {"value": 123.0, "unit": "K"}}},
            models={"vendor_gas": a_model()},
        )
    )
    assert eos.from_model("vendor_gas").components[0].Tc.to("K").magnitude == pytest.approx(123.0)


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
    keycard.use(minimal(models={"vendor_gas": a_model()}))
    with pytest.raises(PropertyUnavailableError, match="vendor_gas"):
        eos.from_model("something_else")


# ---------------------------------------------------------------------------
# The document itself
# ---------------------------------------------------------------------------


def test_a_component_override_is_undone_by_clear() -> None:
    """`clear()` returns the library to the data it ships, with no residue."""
    shipped = eos.components.entry("methane").Pc
    keycard.use(minimal(components={"methane": {"Pc": {"value": 1.0, "unit": "Pa"}}}))
    assert eos.components.entry("methane").Pc != shipped
    keycard.clear()
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
    """`keycard.example.yaml` must pass its own loader, not just its checker.

    Same argument as the template passing `check_user_data.py`: a template that does
    not load teaches the wrong shape, and it is the file a user copies.
    """
    card = keycard.load(TEMPLATE)
    assert card.keyholder == "Example Engineering Ltd"
    assert card.component("methane") is not None
    assert card.kij_for("methane", "n-butane") is not None
    assert card.coefficient("hydraulics.orifice_flow", "Cd") is not None


def test_the_template_is_valid_yaml_with_a_version() -> None:
    """A guard against the file being replaced by prose during an edit."""
    document = yaml.safe_load(TEMPLATE.read_text(encoding="utf-8"))
    assert document["schema_version"] == keycard.SCHEMA_VERSION


# ---------------------------------------------------------------------------
# Which card a call reads
# ---------------------------------------------------------------------------


def a_card(omega: float) -> keycard.Keycard:
    """A keycard overriding one parameter, so two of them differ observably."""
    return keycard.use(
        minimal(components={"methane": {"omega": {"value": omega, "unit": "dimensionless"}}})
    )


def test_an_explicit_card_wins_over_the_loaded_one() -> None:
    """The precedence `spec.md` states for a coefficient, applied to the card itself.

    Both directions are asserted, and they are different. A call that passes no card
    reads the loaded one - which is what "one keycard in force" means, and what every
    existing caller relies on. A call that *does* pass one reads that card and not the
    loaded one, which is the capability being a value: a caller who hands a card to a
    call is stating it for that call, and reading a file loaded an hour ago instead
    would be the surprise the keycard rules exist to prevent.

    A test that only checked the second would pass with the parameter ignored and the
    global read anyway, so both are here.
    """
    loaded = a_card(0.5)
    explicit = a_card(0.25)

    # No card passed: the loaded one answers. `a_card` sets the global, so the second
    # call left 0.25 in force and the first card is the one being held explicitly.
    assert eos.components.entry("methane").omega == pytest.approx(0.25)

    # A card passed: it answers, and the loaded one is not consulted.
    assert eos.components.entry("methane", card=loaded).omega == pytest.approx(0.5)
    assert eos.components.entry("methane", card=explicit).omega == pytest.approx(0.25)

    # And the global is unchanged by a call that passed one, which is what makes the
    # precedence a precedence rather than a mutation.
    assert eos.components.entry("methane").omega == pytest.approx(0.25)


def test_the_card_reaches_the_coefficient_path_too() -> None:
    """A coefficient, not just a component, comes from the card the call passed.

    The two read different sections of the same object, and the resolver is the same
    function, so this is a check that the parameter was threaded rather than that the
    resolver works. `orifice_flow` is the calc with a coefficient a card can supply.
    """
    from azoth.core.units import quantity
    from azoth.hydraulics import orifice_flow
    from azoth.keycard import coefficient_value

    def a_cd(value: float) -> keycard.Keycard:
        return keycard.use(
            minimal(
                coefficients={
                    "hydraulics.orifice_flow": {"Cd": {"value": value, "unit": "dimensionless"}}
                }
            )
        )

    # The second `use` leaves 0.61 in force, so `explicit` is also the loaded card and
    # `loaded` is the one only being held explicitly.
    loaded = a_cd(0.60)
    explicit = a_cd(0.61)

    assert coefficient_value("hydraulics.orifice_flow", "Cd", None, card=loaded) == pytest.approx(
        0.60
    )
    assert coefficient_value("hydraulics.orifice_flow", "Cd", None, card=explicit) == pytest.approx(
        0.61
    )
    assert coefficient_value("hydraulics.orifice_flow", "Cd", None) == pytest.approx(0.61)

    # Through the public calc, so the argument is threaded the whole way rather than
    # only as far as the helper above. `Cd` scales the flow linearly, so the ratio
    # between two cards is the ratio of their coefficients and nothing else - which is
    # what makes this a statement about *which* card answered rather than about the
    # arithmetic.
    bore = quantity(50.0, "mm")
    drop = quantity(10_000.0, "Pa")
    density = quantity(998.0, "kg/m**3")
    at_60 = orifice_flow(bore, drop, density, card=loaded).q.to("m**3/s").magnitude
    at_61 = orifice_flow(bore, drop, density, card=explicit).q.to("m**3/s").magnitude
    assert at_61 == pytest.approx(at_60 * 0.61 / 0.60)

    # And with no card passed, the loaded one answers - `explicit`, since `use` was
    # called on it last.
    assert orifice_flow(bore, drop, density).q.to("m**3/s").magnitude == pytest.approx(at_61)


def test_a_kernel_does_not_take_a_card() -> None:
    """The reference implementations take spec inputs only, and no capability.

    A kernel that took a card would be reading the caller's authority itself rather
    than the values it was handed. The card is resolved at the boundary, and
    `test_registry_contract.py` asserts the same asymmetry across every calc.
    """
    from azoth.hydraulics.reference import orifice_flow as reference

    assert "card" not in inspect.signature(reference).parameters
