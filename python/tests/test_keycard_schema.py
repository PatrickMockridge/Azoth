"""The keycard schema, and the one definition it shares with the calc schema.

A schema is a claim about which files are valid, and a claim nothing tests is a
claim in a comment. So this file does two jobs.

**It proves every constraint fires.** A table of malformed documents, one per rule,
each asserted to be *rejected* - the discipline `check_input_flags` applies to
`interval: true`. A constraint nobody has seen reject anything is
indistinguishable from a constraint with a typo in it, and the failure mode is the
one this repository is organised against: a rule that looks like validation and
does nothing.

**It holds the two schemas to one citation definition.** `keycard.schema.json`
`$ref`s `calc.schema.json#/$defs/citation` rather than restating what a citation is,
because a user's data row and the repository's record the same kind of thing. The test
asserts the reference resolves to *that* definition - so a later refactor that inlines a
copy fails here rather than quietly creating a second answer to what a citation is.

**And it asserts the absence of the old regime.** A keycard used to require a
`verify_status` on every row, one of `verified` / `unverified` / `estimated_dummy`,
with a rule that the citation had to agree with it. That was provenance policing: a
field that could not be checked by a tool, required anyway, which taught people to fill
it in rather than to know the answer. It is gone, and a test says so, because a field
deleted from a schema and left in the tests comes back.
"""

from __future__ import annotations

import copy
import importlib
import json
import sys
from pathlib import Path
from types import ModuleType
from typing import Any, cast

import pytest
from jsonschema import Draft202012Validator

REPO_ROOT = Path(__file__).resolve().parents[2]

CALC_SCHEMA = REPO_ROOT / "specs" / "schema" / "calc.schema.json"
KEYCARD_SCHEMA = REPO_ROOT / "specs" / "schema" / "keycard.schema.json"

CITATION_REF = "calc.schema.json#/$defs/citation"


def load(path: Path) -> dict[str, Any]:
    return cast("dict[str, Any]", json.loads(path.read_text(encoding="utf-8")))


def schema_registry() -> ModuleType:
    """The tools' registry builder, loaded the way the generator tests load theirs.

    That it *is* the tools' builder matters: it reads every schema in
    `specs/schema/` rather than being told about two of them, and the keycard
    schema reaches the calc schema, which reaches the unit schema in turn.
    """
    sys.path.insert(0, str(REPO_ROOT / "tools"))
    try:
        return importlib.import_module("schema_registry")
    finally:
        sys.path.pop(0)


@pytest.fixture(scope="module")
def keycard_validator() -> Draft202012Validator:
    """A validator for the keycard schema with the calc schema resolvable.

    The schemas are separate documents and the references between them are by
    `$id`, so a registry has to carry all of them. Built here rather than in the
    schema because a `$ref` that only resolves at validation time is exactly what
    this fixture exists to exercise.
    """
    keycard = load(KEYCARD_SCHEMA)
    Draft202012Validator.check_schema(keycard)
    registry = schema_registry().registry()
    return Draft202012Validator(keycard, registry=registry)


def a_valid_keycard() -> dict[str, Any]:
    """Every section populated, so a test can break one field at a time."""
    return {
        "schema_version": 2,
        "keyholder": {
            "name": "Acme Process Engineering Ltd",
            "licence": "https://acme.example/licences/iso-5167",
        },
        "fittings": [
            {
                "id": "gate_valve_open",
                "family": "valve",
                "name": "gate valve, fully open",
                "n_ld": 8.0,
                "f_t_basis": "f_t",
                "citation": "Read from the source named below.",
            }
        ],
        "fluids": {
            "seawater": [
                {
                    "temperature_c": 0.0,
                    "density_kg_m3": 1028.0,
                    "dynamic_viscosity_pa_s": 0.00188,
                    "citation": "Read from the source named below.",
                },
                {
                    "temperature_c": 20.0,
                    "density_kg_m3": 1024.0,
                    "dynamic_viscosity_pa_s": 0.00107,
                    "citation": "Read from the source named below.",
                },
            ]
        },
        "components": {
            "methane": {
                "Tc": {
                    "value": 190.564,
                    "unit": "K",
                    "citation": "NeqSim v3.20.0 COMP.csv (Equinor/NTNU), Apache-2.0",
                }
            }
        },
        "coefficients": {
            "hydraulics.orifice_flow": {
                "Cd": {
                    "value": 0.61,
                    "unit": "dimensionless",
                    "convention": "iso_5167_corner",
                    "citation": "Read from the source named below.",
                }
            }
        },
        "models": {
            "vendor_pr": {
                "kind": "cubic_eos",
                "shape": "peng_robinson",
                "alpha": "peng_robinson",
                "mixing_rule": "classical_kij",
                "components": ["methane"],
            }
        },
    }


def errors_for(validator: Draft202012Validator, document: Any) -> list[str]:
    return [error.message for error in validator.iter_errors(document)]


def test_the_schema_is_itself_valid(keycard_validator: Draft202012Validator) -> None:
    """A malformed schema is the worst kind of broken: it validates nothing."""
    assert keycard_validator is not None


def test_a_well_formed_keycard_is_accepted(keycard_validator: Draft202012Validator) -> None:
    """The non-vacuous half. Every case below asserts a rejection, so without this
    one a schema that rejected *everything* would pass the whole file."""
    assert errors_for(keycard_validator, a_valid_keycard()) == []


def test_the_two_schemas_share_one_citation_definition() -> None:
    """The keycard must `$ref` the calc schema's definition, not carry a copy.

    Asserted on the reference *and* on the target's presence, because either alone is
    satisfiable by accident: a `$ref` to a definition that does not exist, or a copy
    that happens to sit under the same name, would each pass one of them.
    """
    keycard = load(KEYCARD_SCHEMA)
    calc = load(CALC_SCHEMA)
    assert "citation" in calc["$defs"], (
        "calc.schema.json no longer defines citation; the keycard schema references "
        "it and the reference would dangle"
    )
    references = json.dumps(keycard)
    assert CITATION_REF in references, (
        f"keycard.schema.json does not reference {CITATION_REF!r}. If the definition "
        f"were inlined, there are now two answers to what a citation is."
    )


#: A citation is a string and nothing else. The rules that used to live beside it -
#: a machine-fetchable `source_ref`, a `source_locator`, a `verify_status` and the
#: agreement between the last two - are gone, so what is left to reject is a citation
#: of the wrong shape and a field that no longer exists.
@pytest.mark.parametrize(
    ("label", "mutate"),
    [
        ("a citation that is not a string", lambda: _row(citation={"source": "tables"})),
        ("a citation that is a number", lambda: _row(citation=42)),
        ("a verify_status, which no longer exists", lambda: _row(verify_status="verified")),
        ("a source_ref, which no longer exists", lambda: _row(source_ref="doi:10.1/x")),
    ],
)
def test_the_removed_provenance_fields_are_gone(
    keycard_validator: Draft202012Validator, label: str, mutate: Any
) -> None:
    """The old regime must be *rejected*, not merely permitted.

    A field deleted from the schema but still accepted is a field that comes back.
    `additionalProperties: false` on the row is what makes this a test rather than a
    hope, and these four cases are what prove the setting is live.
    """
    document = a_valid_keycard()
    document["fittings"] = [mutate()]
    assert errors_for(keycard_validator, document), f"{label} was accepted"


def test_a_row_needs_no_citation_at_all(keycard_validator: Draft202012Validator) -> None:
    """The positive half, and the point of the change.

    A keycard row with no citation, no status and no reference to anything is valid.
    The engineer supplying the data owns its provenance; this library does not ask.
    """
    row = _row()
    del row["citation"]
    document = a_valid_keycard()
    document["fittings"] = [row]
    assert errors_for(keycard_validator, document) == []


def _row(**overrides: Any) -> dict[str, Any]:
    """A fitting row, valid unless a test says otherwise."""
    row = {
        "id": "gate_valve_open",
        "family": "valve",
        "name": "gate valve, fully open",
        "n_ld": 8.0,
        "f_t_basis": "f_t",
        "citation": "Read from the source named below.",
    }
    row.update(overrides)
    return row


@pytest.mark.parametrize(
    ("label", "mutate"),
    [
        ("an unknown top-level section", lambda d: d.update(fitting=[])),
        ("a section misspelt", lambda d: d.update(component={})),
        ("version 1, before the new sections existed", lambda d: d.update(schema_version=1)),
        ("a keyholder with an empty name", lambda d: d["keyholder"].update(name="")),
        (
            "a coefficient with no unit",
            lambda d: d["coefficients"]["hydraulics.orifice_flow"]["Cd"].pop("unit"),
        ),
        (
            "a coefficient in a unit the vocabulary has no conversion for",
            lambda d: d["coefficients"]["hydraulics.orifice_flow"]["Cd"].update(unit="furlong"),
        ),
        (
            "a component parameter with no value",
            lambda d: d["components"]["methane"]["Tc"].pop("value"),
        ),
        (
            "a model kind outside the closed list",
            lambda d: d["models"]["vendor_pr"].update(kind="gerg_2008"),
        ),
        (
            "an alpha function outside the closed list",
            lambda d: d["models"]["vendor_pr"].update(alpha="whatever_works"),
        ),
        (
            "a mathias_copeman alpha with no eta",
            lambda d: d["models"]["vendor_pr"].update(alpha="mathias_copeman"),
        ),
        (
            "a twu alpha missing one of its three parameters",
            lambda d: d["models"]["vendor_pr"].update(
                alpha="twu", alpha_parameters={"L": 1.0, "M": 2.0}
            ),
        ),
        ("a model with no components", lambda d: d["models"]["vendor_pr"].pop("components")),
    ],
)
def test_the_keycard_rules_reject_what_they_should(
    keycard_validator: Draft202012Validator, label: str, mutate: Any
) -> None:
    """One malformed keycard per constraint, each proven to fire.

    The `if/then` branches are the ones that matter here. A schema that accepted
    `alpha: mathias_copeman` without `eta` would let a caller select one equation
    and have the implementation silently evaluate another - which is the failure
    `check_input_flags` exists to catch for `interval: true`, in a different file.
    """
    document = copy.deepcopy(a_valid_keycard())
    mutate(document)
    assert errors_for(keycard_validator, document), f"{label} was accepted"
