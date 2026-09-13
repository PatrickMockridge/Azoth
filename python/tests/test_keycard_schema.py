"""The keycard schema, and the one definition it shares with the calc schema.

A schema is a claim about which files are valid, and a claim nothing tests is a
claim in a comment. So this file does two jobs.

**It proves every constraint fires.** A table of malformed documents, one per rule,
each asserted to be *rejected* - the discipline `check_input_flags` applies to
`interval: true`. A constraint nobody has seen reject anything is
indistinguishable from a constraint with a typo in it, and the failure mode is the
one this repository is organised against: a rule that looks like validation and
does nothing.

**It holds the two schemas to one provenance definition.** `keycard.schema.json`
`$ref`s `calc.schema.json#/$defs/provenance_row` rather than restating the
citation rules, because a user's data row and the repository's are the same kind of
claim. The test asserts the reference resolves to *that* definition - so a later
refactor that inlines a copy fails here rather than quietly creating a second
answer to what a valid citation is. That is the same drift that had already
happened once between `check_source` and its duplicate in the user-data checker.
"""

from __future__ import annotations

import copy
import json
from pathlib import Path
from typing import Any, cast

import pytest
from jsonschema import Draft202012Validator
from referencing import Registry, Resource

REPO_ROOT = Path(__file__).resolve().parents[2]
CALC_SCHEMA = REPO_ROOT / "specs" / "schema" / "calc.schema.json"
KEYCARD_SCHEMA = REPO_ROOT / "specs" / "schema" / "keycard.schema.json"

PROVENANCE_REF = "calc.schema.json#/$defs/provenance_row"


def load(path: Path) -> dict[str, Any]:
    return cast("dict[str, Any]", json.loads(path.read_text(encoding="utf-8")))


@pytest.fixture(scope="module")
def keycard_validator() -> Draft202012Validator:
    """A validator for the keycard schema with the calc schema resolvable.

    The two schemas are separate documents and the reference between them is by
    `$id`, so a registry has to carry both. Built here rather than in the schema
    because a `$ref` that only resolves at validation time is exactly what this
    fixture exists to exercise.
    """
    calc = load(CALC_SCHEMA)
    keycard = load(KEYCARD_SCHEMA)
    Draft202012Validator.check_schema(keycard)
    registry = Registry()
    for document in (calc, keycard):
        registry = registry.with_resource(document["$id"], Resource.from_contents(document))
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
                "verify_status": "verified",
                "source_ref": "arweave:AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
                "source_locator": "Table 2, gate valve, fully open",
            }
        ],
        "fluids": {
            "seawater": [
                {
                    "temperature_c": 0.0,
                    "density_kg_m3": 1028.0,
                    "dynamic_viscosity_pa_s": 0.00188,
                    "citation": "Read from the source named below.",
                    "verify_status": "unverified",
                    "source_ref": "https://example.test/seawater",
                    "source_locator": "Table 1",
                },
                {
                    "temperature_c": 20.0,
                    "density_kg_m3": 1024.0,
                    "dynamic_viscosity_pa_s": 0.00107,
                    "citation": "Read from the source named below.",
                    "verify_status": "unverified",
                    "source_ref": "https://example.test/seawater",
                    "source_locator": "Table 1",
                },
            ]
        },
        "components": {
            "methane": {
                "Tc": {
                    "value": 190.564,
                    "unit": "K",
                    "citation": "NeqSim v3.20.0 COMP.csv (Equinor/NTNU), Apache-2.0",
                    "verify_status": "unverified",
                    "source_ref": "https://github.com/equinor/neqsim",
                    "source_locator": "COMP.csv, methane",
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
                    "verify_status": "verified",
                    "source_ref": "https://iso.example/5167-2",
                    "source_locator": "Table 4, corner tappings, beta = 0.5",
                }
            }
        },
        "models": {
            "vendor_pr": {
                "kind": "cubic_eos",
                "shape": "peng_robinson",
                "alpha": "peng_robinson_1978",
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


def test_the_two_schemas_share_one_provenance_definition() -> None:
    """The keycard must `$ref` the calc schema's rules, not carry a copy.

    Asserted on the reference *and* on the target's presence, because either alone
    is satisfiable by accident: a `$ref` to a definition that does not exist, or a
    copy that happens to sit under the same name, would each pass one of them.
    """
    keycard = load(KEYCARD_SCHEMA)
    calc = load(CALC_SCHEMA)
    assert "provenance_row" in calc["$defs"], (
        "calc.schema.json no longer defines provenance_row; the keycard schema "
        "references it and the reference would dangle"
    )
    references = json.dumps(keycard)
    assert PROVENANCE_REF in references, (
        f"keycard.schema.json does not reference {PROVENANCE_REF!r}. If the rules "
        f"were inlined, there are now two definitions of what a valid citation is."
    )


#: Every rule in `provenance_row`, one malformed row each. Mirrors the branches
#: `tools/spec_lint.py::check_source` enforces on the CSV files, because the two
#: must agree and only one of them is machine-checked here.
@pytest.mark.parametrize(
    ("label", "mutate"),
    [
        (
            "placeholder whose citation does not say DUMMY",
            lambda: _row(citation="a plausible number", verify_status="estimated_dummy"),
        ),
        (
            "placeholder that also claims a source",
            lambda: _row(verify_status="estimated_dummy", source_ref="https://x.test/y"),
        ),
        (
            "verified with no source_ref",
            lambda: _row(source_ref=""),
        ),
        (
            "source_ref in prose rather than a fetchable form",
            lambda: _row(source_ref="the standard, section 4"),
        ),
        (
            "a source_ref that is not a whole reference",
            lambda: _row(source_ref="arweave:tooshort"),
        ),
        (
            "a source_ref with no locator",
            lambda: _row(source_locator=""),
        ),
        (
            "promoted to verified without editing the citation",
            lambda: _row(citation="DUMMY value, not from any source"),
        ),
        (
            "an unrecognised verify_status",
            lambda: _row(verify_status="probably_fine"),
        ),
    ],
)
def test_the_provenance_rules_reject_what_they_should(
    keycard_validator: Draft202012Validator, label: str, mutate: Any
) -> None:
    document = a_valid_keycard()
    document["fittings"] = [mutate()]
    assert errors_for(keycard_validator, document), f"{label} was accepted"


def _row(**overrides: str) -> dict[str, Any]:
    """A fitting row, valid unless a test says otherwise."""
    row = {
        "id": "gate_valve_open",
        "family": "valve",
        "name": "gate valve, fully open",
        "n_ld": 8.0,
        "f_t_basis": "f_t",
        "citation": "Read from the source named below.",
        "verify_status": "verified",
        "source_ref": "https://example.test/table",
        "source_locator": "Table 2",
    }
    row.update(overrides)
    return row


@pytest.mark.parametrize(
    ("label", "mutate"),
    [
        ("an unknown top-level section", lambda d: d.update(fitting=[])),
        ("a section misspelt", lambda d: d.update(component={})),
        ("version 1, before the new sections existed", lambda d: d.update(schema_version=1)),
        ("no keyholder", lambda d: d.pop("keyholder")),
        ("a keyholder with an empty name", lambda d: d["keyholder"].update(name="")),
        (
            "a licence that is not a fetchable reference",
            lambda d: d["keyholder"].update(licence="ask us"),
        ),
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
