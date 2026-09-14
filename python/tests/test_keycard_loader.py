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
vocabularies are a deliberate copy, and a copy without a test is the drift this project
answers with two implementations compared against each other. `UNIT_VOCABULARY` is the
sharpest case: `pint` parses `kelvin`, so a loader that trusted `pint` would accept a
keycard that `tools/check_user_data.py` then rejects.
"""

from __future__ import annotations

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
    schema = json.loads(CALC_SCHEMA.read_text(encoding="utf-8"))
    declared = tuple(schema["$defs"]["unit"]["enum"])
    assert declared == keycard.UNIT_VOCABULARY, (
        "the loader's unit vocabulary and calc.schema.json's enum have drifted. One of "
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
