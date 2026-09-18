#!/usr/bin/env python3
"""Compile the unit vocabulary table into Rust, Python and JSON Schema.

Reads `specs/vocabulary/vocabulary.toml` and emits:

  specs/schema/unit.schema.json              the closed enum of canonical units
  crates/azoth-core/src/unit_vocab_gen.rs    UNIT_NAMES, dimensions, conversions
  python/src/azoth/core/_units_gen.py        the same, in Python
  lean/Azoth/Vocabulary.lean                 the same, as a dimension per unit
  lean/Azoth/Gate.lean                       one `#print axioms` line per unit

The table is the one hand-written source of that list. Before it existed the same
24 strings were maintained by hand in four places - the calc schema's
`$defs.unit.enum`, Rust's `UNIT_NAMES`, Python's `CANONICAL_UNITS` and the
keycard loader's `UNIT_VOCABULARY` - plus a fifth place that was worse than the
other four: `CONVERSION_PATHS` in `crates/azoth-core/src/units.rs` paired each
name with a hand-typed SI magnitude. Two of those lists had already drifted into
disagreeing with each other twice.

**What this generator deliberately does not emit: a conversion factor.** A factor
is a number `uom` and `pint` each already know. The generated Rust carries the
conversion a calculation runs and no expected value for it; the check that the
conversion is right is a test that compares it against `pint`'s own answer for
the same name, in the language where `pint` is available. Writing the number down
here would be the same defect in a new file.

**The uom quantity type is derived, not declared.** The table gives a unit's
exponents and its path within uom; the quantity type comes from `UOM_TYPES`
below, keyed by the exponents. That ordering is what makes the generated
assertion bite: a unit declared with the wrong dimension produces an assignment
whose two sides are different types, so it fails to compile rather than being
waved through by a name that happens to match.

Usage:
    python tools/gen_vocabulary.py            # write the files
    python tools/gen_vocabulary.py --check    # fail if they are out of date
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
import tomllib
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent.parent
VOCAB_PATH = ROOT / "specs" / "vocabulary" / "vocabulary.toml"
SCHEMA_PATH = ROOT / "specs" / "schema" / "vocabulary.schema.json"
UNIT_SCHEMA_OUT = ROOT / "specs" / "schema" / "unit.schema.json"
RUST_OUT = ROOT / "crates" / "azoth-core" / "src" / "unit_vocab_gen.rs"
PY_OUT = ROOT / "python" / "src" / "azoth" / "core" / "_units_gen.py"
LEAN_OUT = ROOT / "lean" / "Azoth" / "Vocabulary.lean"
LEAN_GATE_OUT = ROOT / "lean" / "Azoth" / "Gate.lean"

GENERATED_BANNER = "GENERATED FILE - DO NOT EDIT BY HAND."

#: The uom quantity type for each dimension, in `slots` order.
#:
#: This is the map that makes the table's exponents load-bearing. A dimension
#: absent from it is refused rather than skipped, so adding a dimension to the
#: vocabulary means saying whether uom carries it - and `None` is an answer, not
#: an omission. `None` means uom has no such quantity, which is a structural fact
#: about uom rather than a gap here: `quantity!` has to run inside `uom::si`, and
#: the `Units` trait that makes `.new::<u>()` compile is assembled by `system!`
#: from a closed list, so no crate outside uom can add a quantity to `uom::si`.
#:
#: uom's module for a quantity is its name in snake case, which is why the table
#: records a path rather than a module and a unit: the module is checked against
#: this map rather than restated beside it.
UOM_TYPES: dict[tuple[int, ...], str | None] = {
    (0, 0, 0, 0, 0, 0, 0): None,  # dimensionless: a ratio, carried as a bare f64
    (1, 0, 0, 0, 0, 0, 0): "Length",
    (2, 0, 0, 0, 0, 0, 0): "Area",
    (3, 0, -1, 0, 0, 0, 0): "VolumeRate",
    (0, 1, -1, 0, 0, 0, 0): "MassRate",
    (0, 0, -1, 0, 0, 1, 0): None,  # MolarFlow: uom has no molar-flow quantity
    (-3, 1, 0, 0, 0, 0, 0): "MassDensity",
    (1, 0, -1, 0, 0, 0, 0): "Velocity",
    (2, 0, -1, 0, 0, 0, 0): "DiffusionCoefficient",
    (-1, 1, -2, 0, 0, 0, 0): "Pressure",
    (-1, 1, -1, 0, 0, 0, 0): "DynamicViscosity",
    (0, 0, 0, 0, 1, 0, 0): "ThermodynamicTemperature",
    (2, 1, -3, 0, 0, 0, 0): "Power",
    (2, 0, -2, 0, -1, 0, 0): "SpecificHeatCapacity",
    (1, 1, -3, 0, -1, 0, 0): "ThermalConductivity",
    (0, 1, -3, 0, -1, 0, 0): "HeatTransfer",
    (0, 1, -2, 0, 0, 0, 0): "SurfaceTension",
    (0, 1, 0, 0, 0, -1, 0): "MolarMass",
    (3, 0, 0, 0, 0, -1, 0): "MolarVolume",
    (2, 1, -2, 0, 0, -1, 0): "MolarEnergy",
    # The equation of state's attraction parameter: uom carries a pressure and a
    # molar volume but no quantity for their product with a length to the sixth,
    # so this is the third dimension azoth names itself.
    (5, 1, -2, 0, 0, -2, 0): None,
    (2, 1, -2, 0, -1, -1, 0): "MolarHeatCapacity",
    # The Cp polynomial coefficients: uom carries no quantity for a fractional
    # power of the temperature unit either, so these are the second group azoth
    # names itself. They are already their own SI base unit, so the generated
    # conversion is the identity.
    (2, 1, -2, 0, -2, -1, 0): None,
    (2, 1, -2, 0, -3, -1, 0): None,
    (2, 1, -2, 0, -4, -1, 0): None,
    (2, 1, -2, 0, -5, -1, 0): None,
}


#: The `lean-units` dimension for each dimension in the table, as a Lean term.
#:
#: This is what makes the emitted theorem a *check* rather than a tautology, and it
#: is keyed by **unit** rather than by dimension on purpose.
#:
#: Keyed by dimension it would prove nothing: a table entry whose exponents were
#: wrong would select the expression for those wrong exponents, and the theorem
#: would hold. Keyed by unit, this is the hand-written claim about what each unit
#: *is* - the counterpart of `rust_ctor`, which is what makes the Rust assertion
#: bite - and the theorem forces the table's exponents to agree with it.
#:
#: So a unit declared an area when this map says it is a length fails to prove, and
#: so does a slot order that moves the exponents under it. A unit with no entry here
#: is refused rather than skipped, which is the property that keeps the map total as
#: the table grows.
#:
#: Every unit in the vocabulary is expressible this way, which is why there is no
#: `None` here as there is in `UOM_TYPES`: `uom` names about twenty quantities, and
#: `lean-units` names the base dimensions and enough derived ones to build the rest
#: by division.
LEAN_DIMENSIONS: dict[str, str] = {
    "dimensionless": "0",
    "m": "Dimension.Length",
    "mm": "Dimension.Length",
    "m**2": "Dimension.Area",
    "m**3/s": "Dimension.Volume / Dimension.Time",
    "kg/s": "Dimension.Mass / Dimension.Time",
    "mol/s": "Dimension.AmountOfSubstance / Dimension.Time",
    "kg/m**3": "Dimension.Mass / Dimension.Volume",
    "m/s": "Dimension.Speed",
    "m**2/s": "Dimension.Area / Dimension.Time",
    "Pa": "Dimension.Pressure",
    "Pa*s": "Dimension.Pressure * Dimension.Time",
    "K": "Dimension.Temperature",
    "W": "Dimension.Power",
    "J/(kg*K)": "Dimension.Energy / (Dimension.Mass * Dimension.Temperature)",
    "W/(m*K)": "Dimension.Power / (Dimension.Length * Dimension.Temperature)",
    "W/(m**2*K)": "Dimension.Power / (Dimension.Area * Dimension.Temperature)",
    "N/m": "Dimension.Force / Dimension.Length",
    "kg/mol": "Dimension.Mass / Dimension.AmountOfSubstance",
    "m**3/mol": "Dimension.Volume / Dimension.AmountOfSubstance",
    "J/mol": "Dimension.Energy / Dimension.AmountOfSubstance",
    "Pa*m**6/mol**2": "Dimension.Pressure * Dimension.Length ^ 6 / Dimension.AmountOfSubstance ^ 2",
    "J/(mol*K)": "Dimension.Energy / (Dimension.AmountOfSubstance * Dimension.Temperature)",
    "J/(mol*K**2)": "Dimension.Energy / (Dimension.AmountOfSubstance * Dimension.Temperature ^ 2)",
    "J/(mol*K**3)": "Dimension.Energy / (Dimension.AmountOfSubstance * Dimension.Temperature ^ 3)",
    "J/(mol*K**4)": "Dimension.Energy / (Dimension.AmountOfSubstance * Dimension.Temperature ^ 4)",
    "J/(mol*K**5)": "Dimension.Energy / (Dimension.AmountOfSubstance * Dimension.Temperature ^ 5)",
}


def snake(name: str) -> str:
    """`VolumeRate` -> `volume_rate`, which is uom's module name for it."""
    return "".join(f"_{c.lower()}" if c.isupper() else c for c in name).lstrip("_")


def load_table() -> dict[str, Any]:
    """Read and validate the table, and refuse a shape this generator cannot compile."""
    try:
        from jsonschema import Draft202012Validator
    except ImportError:  # pragma: no cover
        sys.exit("gen_vocabulary requires jsonschema: pip install jsonschema")

    table: dict[str, Any] = tomllib.loads(VOCAB_PATH.read_text(encoding="utf-8"))
    schema = json.loads(SCHEMA_PATH.read_text(encoding="utf-8"))
    Draft202012Validator.check_schema(schema)
    errors = sorted(Draft202012Validator(schema).iter_errors(table), key=lambda e: list(e.path))
    if errors:
        for error in errors:
            print(
                f"gen_vocabulary: {VOCAB_PATH.name}: {list(error.path)}: {error.message}",
                file=sys.stderr,
            )
        sys.exit("gen_vocabulary: the vocabulary table does not validate")

    slots = table["slots"]
    width = len(slots)

    # Duplicates are refused rather than resolved. Every generated map keys on the
    # name, so a second entry would silently drop the first - and the entry that
    # disappeared would be the one a reviewer had just read.
    dim_ids = [d["id"] for d in table["dimensions"]]
    dupes = {i for i in dim_ids if dim_ids.count(i) > 1}
    if dupes:
        sys.exit(f"gen_vocabulary: duplicate dimension id(s): {sorted(dupes)}")
    unit_ids = [u["id"] for u in table["units"]]
    dupes = {i for i in unit_ids if unit_ids.count(i) > 1}
    if dupes:
        sys.exit(f"gen_vocabulary: duplicate unit id(s): {sorted(dupes)}")

    # Every unit needs a hand-written lean-units expression, or the theorem the
    # generator emits for it has no independent right-hand side to check against.
    unnamed = sorted(set(unit_ids) - set(LEAN_DIMENSIONS))
    if unnamed:
        sys.exit(
            f"gen_vocabulary: unit(s) {unnamed} have no entry in LEAN_DIMENSIONS, so "
            f"there is nothing for their dimension to be checked against. Add one: "
            f"that map is the hand-written claim about what each unit is, and the "
            f"generated theorem is what makes the table agree with it."
        )

    by_id = {d["id"]: d for d in table["dimensions"]}
    for dimension in table["dimensions"]:
        if len(dimension["exponents"]) != width:
            sys.exit(
                f"gen_vocabulary: dimension {dimension['id']!r} has "
                f"{len(dimension['exponents'])} exponents for {width} slot(s)"
            )
        key = tuple(dimension["exponents"])
        if key not in UOM_TYPES:
            sys.exit(
                f"gen_vocabulary: dimension {dimension['id']!r} {list(key)} is not in "
                f"UOM_TYPES, so this generator does not know whether uom carries it. "
                f"Add it with a quantity type, or with None if uom has none."
            )

    for unit in table["units"]:
        dimension = by_id.get(unit["dimension"])
        if dimension is None:
            sys.exit(
                f"gen_vocabulary: unit {unit['id']!r} names dimension "
                f"{unit['dimension']!r}, which is not declared"
            )
        exponents = tuple(dimension["exponents"])
        quantity = UOM_TYPES[exponents]
        uom = unit.get("uom")
        if (uom is None) != (quantity is None):
            sys.exit(
                f"gen_vocabulary: unit {unit['id']!r} names uom path {uom!r} but "
                f"dimension {unit['dimension']!r} maps to {quantity!r}. A unit's uom "
                f"path and its dimension's quantity type have to agree in both "
                f"directions, or the generated conversion and the generated assertion "
                f"would be about different things."
            )
        if uom is not None:
            module = uom.split("::")[0]
            if module != snake(quantity):  # type: ignore[arg-type]
                sys.exit(
                    f"gen_vocabulary: unit {unit['id']!r} gives the uom path "
                    f"{uom!r}, whose module is {module!r}, but dimension "
                    f"{unit['dimension']!r} is uom's {quantity!r}, whose module is "
                    f"{snake(quantity)!r}"
                )

    return table


def emit_rust(table: dict[str, Any]) -> str:
    slots = table["slots"]
    by_id = {d["id"]: d for d in table["dimensions"]}
    units = table["units"]

    out = [
        f"//! {GENERATED_BANNER}",
        "//!",
        "//! Generated by `tools/gen_vocabulary.py` from",
        "//! `specs/vocabulary/vocabulary.toml`.",
        "",
        "/// The base dimensions, in the order every exponent tuple here is written in.",
        "///",
        "/// uom::si::ISQ's order, and not cosmetic: the tuples below and the uom",
        "/// quantity types the table's assertions are written against are both built",
        "/// from it, so a different order permutes all of them at once.",
        f"pub const SLOTS: [&str; {len(slots)}] = [{', '.join(json.dumps(s) for s in slots)}];",
        "",
        "/// Every named dimension a channel field may declare, in the vocabulary's own",
        "/// words. A dimension id is the canonical name for one exponent tuple; two",
        "/// fields are compatible exactly when their dimensions' tuples are equal.",
        "pub const DIMENSION_IDS: &[&str] = &[",
    ]
    for dimension in table["dimensions"]:
        out.append(f"    {json.dumps(dimension['id'])},")
    out += [
        "];",
        "",
        "/// Each dimension id, as exponents in [`SLOTS`] order.",
        f"pub const DIMENSION_EXPONENTS: &[(&str, [i8; {len(slots)}])] = &[",
    ]
    for dimension in table["dimensions"]:
        out.append(f"    ({json.dumps(dimension['id'])}, {rust_i8_array(dimension['exponents'])}),")
    out += [
        "];",
        "",
        "/// The exponents of one named dimension, or `None` if the id is not in the",
        "/// vocabulary.",
        "#[must_use]",
        f"pub fn dimension_exponents(id: &str) -> Option<[i8; {len(slots)}]> {{",
        "    DIMENSION_EXPONENTS",
        "        .iter()",
        "        .find(|(n, _)| *n == id)",
        "        .map(|(_, d)| *d)",
        "}",
        "",
        "/// Every canonical unit string a spec may declare.",
        "///",
        "/// A name here is a claim that this crate has a correct conversion path for",
        "/// it - see `conversion`, and the tests in `crate::units` that hold the two",
        "/// lists to each other.",
        "pub const UNIT_NAMES: &[&str] = &[",
    ]
    for unit in units:
        out.append(f"    {json.dumps(unit['id'])},")
    out += [
        "];",
        "",
        "/// Each unit's dimension, as exponents in [`SLOTS`] order.",
        f"pub const UNIT_DIMENSIONS: &[(&str, [i8; {len(slots)}])] = &[",
    ]
    for unit in units:
        exps = by_id[unit["dimension"]]["exponents"]
        out.append(f"    ({json.dumps(unit['id'])}, {rust_i8_array(exps)}),")
    out += [
        "];",
        "",
        "/// One row of [`CONVERSION_PATHS`]: a canonical unit string, and the",
        "/// conversion this crate performs for it.",
        "pub type ConversionPath = (&'static str, fn(f64) -> f64);",
        "",
        "/// Each unit's conversion from the declared unit to the SI base magnitude.",
        "///",
        "/// A unit whose dimension uom carries goes through the hand-written",
        "/// constructor in [`crate::units`]; one whose dimension it does not is",
        "/// already its own SI base unit, so the conversion is the identity.",
        "///",
        "/// There is deliberately no expected magnitude beside these. A factor is a",
        "/// number `uom` and `pint` each already know, and the check that these are",
        "/// right is `python/tests/test_units_cross_library.py`, which compares the",
        "/// result of calling one of these against `pint`'s answer for the same name.",
        "/// Writing the number down here would be the defect this generated file",
        "/// exists to remove.",
        "pub const CONVERSION_PATHS: &[ConversionPath] = &[",
    ]
    for unit in units:
        ctor = unit.get("rust_ctor")
        body = "|v| v" if ctor is None else f"|v| crate::units::{ctor}(v).value"
        out.append(f"    ({json.dumps(unit['id'])}, {body}),")
    out += [
        "];",
        "",
        "/// The conversion for one canonical unit, or `None` if the name is not in",
        "/// the vocabulary.",
        "#[must_use]",
        "pub fn conversion(name: &str) -> Option<fn(f64) -> f64> {",
        "    CONVERSION_PATHS",
        "        .iter()",
        "        .find(|(n, _)| *n == name)",
        "        .map(|(_, f)| *f)",
        "}",
        "",
        "/// The SI base magnitude one of `name` is worth, by running the same",
        "/// conversion a calculation runs.",
        "///",
        "/// This is what makes the conversion checkable from outside: the caller that",
        "/// reads it is `azoth._core.unit_si_factor`, and what it is compared against",
        "/// is `pint`.",
        "#[must_use]",
        "pub fn si_factor(name: &str) -> Option<f64> {",
        "    conversion(name).map(|f| f(1.0))",
        "}",
        "",
        "/// The dimension of one canonical unit, or `None` if the name is not in the",
        "/// vocabulary.",
        "#[must_use]",
        f"pub fn dimension(name: &str) -> Option<[i8; {len(slots)}]> {{",
        "    UNIT_DIMENSIONS",
        "        .iter()",
        "        .find(|(n, _)| *n == name)",
        "        .map(|(_, d)| *d)",
        "}",
        "",
        "#[cfg(test)]",
        "mod dimension_assertions {",
        "    //! One assertion per unit, and the whole point of deriving the uom type",
        "    //! from the exponents rather than declaring it beside them.",
        "    //!",
        "    //! Each line says: the constructor the conversion calls produces the uom",
        "    //! quantity that THIS unit's exponents name. Change an exponent in the",
        "    //! table and the two sides stop being the same type, so this fails to",
        "    //! compile rather than passing on a name that happens to match.",
        "",
        "    #[test]",
        "    fn the_table_agrees_with_uom() {",
    ]
    for unit in units:
        exps = by_id[unit["dimension"]]["exponents"]
        quantity = UOM_TYPES[tuple(exps)]
        if quantity is None:
            continue
        out.append(
            f"        let _: uom::si::f64::{quantity} = crate::units::{unit['rust_ctor']}(1.0);"
        )
    out += [
        "    }",
        "}",
        "",
    ]
    return "\n".join(out)


def emit_python(table: dict[str, Any]) -> str:
    slots = table["slots"]
    units = table["units"]

    out = [
        f'"""{GENERATED_BANNER}',
        "",
        "Generated by `tools/gen_vocabulary.py` from",
        "`specs/vocabulary/vocabulary.toml`.",
        "",
        "The same data the Rust side reads from `azoth_core::unit_vocab_gen`, and",
        "generated from the same table so the two cannot name different sets.",
        '"""',
        "",
        "from __future__ import annotations",
        "",
        "from typing import Final",
        "",
        "#: The base dimensions, in the order every exponent tuple here is written in.",
        "#: uom::si::ISQ's order; see the table.",
        "SLOTS: Final[tuple[str, ...]] = (" + ", ".join(json.dumps(s) for s in slots) + ")",
        "",
        "#: Canonical unit strings, keyed by the strings the spec schema allows, mapped",
        "#: to the unit's name in `pint`'s registry. This is the unit a spec and its",
        "#: generated docs speak in - what a worked example's numbers mean, and what a",
        "#: caller gets back. It is deliberately not the unit anything is calculated",
        "#: in: calculations work in SI base magnitudes, and `to_si`/`from_si` are the",
        "#: only places the two representations meet.",
        "CANONICAL_UNITS: Final[dict[str, str]] = {",
    ]
    for unit in units:
        out.append(f"    {json.dumps(unit['id'])}: {json.dumps(unit['pint'])},")
    out += [
        "}",
        "",
        "#: The canonical unit strings, in table order. The keycard loader validates",
        "#: against this rather than reading the schema, which is a dev-time artefact",
        "#: and does not ship in the wheel.",
        "UNIT_VOCABULARY: Final[tuple[str, ...]] = tuple(CANONICAL_UNITS)",
        "",
        "",
        "#: `SLOTS` is the table's slot order, and it is here because something reads",
        "#: it: `python/tests/test_units_cross_library.py` maps `pint`'s base dimensions",
        "#: onto these names, and compares them with the Rust side's through",
        "#: `azoth._core.unit_slots`. The table's *exponents* are deliberately not",
        "#: carried here - Python's runtime need is a name and its `pint` spelling, and a",
        "#: dimension is something only a caller asking for one needs. They arrive when",
        "#: such a caller does, not before.",
    ]
    return "\n".join(out) + "\n"


def emit_lean(table: dict[str, Any]) -> str:
    """The vocabulary as Lean data, one dimension per canonical unit.

    Data only: the table, and a lookup. Theorems *about* it are hand-written,
    because "everything that names them is generated" and a proven statement is
    not a naming - so the generator emits what a unit is and a proof says it means
    what the calculus says it means.
    """
    by_id = {d["id"]: d for d in table["dimensions"]}
    units = table["units"]

    width = len(table["slots"])
    out = [
        f"-- {GENERATED_BANNER}",
        "--",
        "-- Generated by `tools/gen_vocabulary.py` from",
        "-- `specs/vocabulary/vocabulary.toml`.",
        "--",
        "-- One entry per canonical unit a spec may declare, with the dimension it",
        "-- carries as exponents in `Azoth.Dim.slots` order. The theorems that say",
        "-- these mean what the calculus says they mean are hand-written, in",
        "-- `Azoth/Units.lean`, because a proven statement is not a naming.",
        "",
        "import Azoth.Dim",
        "",
        "namespace Azoth.Vocabulary",
        "",
        "open Units (Dimension)",
        "",
        "/-- Every canonical unit string, with the dimension it carries.",
        "",
        "Each exponent vector is written in `Azoth.Dim.slots` order, and has one",
        f"entry per slot - {width} of them.",
        "-/",
        # `\u00d7` and not the character itself: it is Lean's product type, and
        # writing it literally makes RUF001 read it as an ambiguous multiplication
        # sign in a string that is not prose.
        "def units : List (String \u00d7 Units.Dimension) :=",
        "  [",
    ]
    for index, unit in enumerate(units):
        exponents = by_id[unit["dimension"]]["exponents"]
        rendered = "[" + ", ".join(str(e) for e in exponents) + "]"
        comma = "," if index + 1 < len(units) else ""
        out.append(f"    ({json.dumps(unit['id'])}, Dim.ofExponents {rendered}){comma}")
    out += [
        "  ]",
        "",
        "/-- The dimension of a canonical unit, or `none` if the name is not in the",
        "vocabulary.",
        "",
        "A name that is not there is `none` rather than an error: asking about one is a",
        "question about the vocabulary, and the answer is that it has no dimension",
        "because it is not in it.",
        "-/",
        "def dimOf (name : String) : Option Units.Dimension :=",
        "  (units.find? (fun entry => entry.1 == name)).map (fun entry => entry.2)",
        "",
        emit_lean_theorems(table),
        "end Azoth.Vocabulary",
        "",
    ]
    return "\n".join(out)


#: What the per-unit proof unfolds, and nothing else. The named dimensions of
#: `lean-units` come first: unfolding them is what turns the right-hand side into
#: base-dimension arithmetic the tactics below can normalise.
_DIMENSION_LEMMAS = " ".join(
    [
        "dimOf,",
        "units,",
        "Azoth.Dim.ofExponents,",
        "Azoth.Dim.ofExponentsOn,",
        "Azoth.slots,",
        "Units.Dimension.Acceleration,",
        "Units.Dimension.AmountOfSubstance,",
        "Units.Dimension.Area,",
        "Units.Dimension.Energy,",
        "Units.Dimension.Force,",
        "Units.Dimension.Length,",
        "Units.Dimension.Mass,",
        "Units.Dimension.Power,",
        "Units.Dimension.Pressure,",
        "Units.Dimension.Speed,",
        "Units.Dimension.Temperature,",
        "Units.Dimension.Time,",
        "Units.Dimension.Volume,",
        "Units.Dimension.ofString,",
        "Units.Dimension.div_eq_sub,",
        "Units.Dimension.mul_eq_add,",
        "Units.Dimension.npow_eq_nsmul,",
        "sub_eq_add_neg,",
        "List.zip_cons_cons,",
        "List.zip_nil_right,",
        "List.map_cons,",
        "List.map_nil,",
        "List.sum_cons,",
        "List.sum_nil",
    ]
)


def lean_name(unit_id: str) -> str:
    """A canonical unit string as a Lean identifier.

    `J/(mol*K**2)` is not one, so the punctuation a unit may contain is spelled
    out rather than dropped - dropping it would map `m**2` and `m**3` onto the
    same name, and the second would shadow the first in the namespace.
    """
    spelled = (
        unit_id.replace("**", "_pow_")
        .replace("*", "_times_")
        .replace("/", "_per_")
        .replace("(", "_")
        .replace(")", "_")
    )
    cleaned = "".join(c if (c.isalnum() or c == "_") else "_" for c in spelled)
    return "u_" + cleaned.strip("_") + "_dimension"


def emit_lean_theorems(table: dict[str, Any]) -> str:
    """One theorem per canonical unit, checking its dimension against `lean-units`.

    The counterpart of the compile-time assertion on the Rust side, and a *check*
    rather than a restatement: the right-hand side is a dimension with a name
    `lean-units` wrote, in a system that was in this tree before this table
    existed. A wrong slot order, or an exponent in the wrong coordinate, fails to
    prove - `ofExponents [1, 0, 0, 0, 0, 0, 0]` is `Length` only because the slot
    order says the first slot is length.

    The proof is `simp` to push `_impl` through the module operations and reduce
    `\u211a`-scaled sums, then `module` to normalise the result. `abel` does not
    close these: it fails on a goal containing a negation, which every unit whose
    dimension has a negative exponent produces.
    """
    out: list[str] = []
    for unit in table["units"]:
        expression = LEAN_DIMENSIONS[unit["id"]]
        out += [
            f"/-- {json.dumps(unit['id'])} carries the dimension `lean-units` calls",
            f"    `{expression}`. -/",
            f"theorem {lean_name(unit['id'])} :",
            # Parenthesised whether or not the expression needs it, so the shape
            # does not depend on which units happen to have compound dimensions.
            f"    dimOf {json.dumps(unit['id'])} = some ({expression}) := by",
            # Chained on one line and each step wrapped in `try`, because `simp
            # only` closes some of these outright and a tactic applied where there
            # is nothing left to do is an error rather than a no-op. `try` cannot
            # hide a failure: a step that does not close the goal leaves it open,
            # and the theorem then does not prove.
            f"  simp only [{_DIMENSION_LEMMAS}] <;> simp <;> module",
            "",
        ]
    return "\n".join(out)


def emit_lean_gate(table: dict[str, Any]) -> str:
    """One `#print axioms` line per unit theorem, so the gate cannot miss one.

    Generated for the reason the rest of this file is: a hand-maintained list of
    twenty-four names goes stale the first time a unit is added, and it goes stale
    *silently* - the theorem is proved and nothing gates it. Here the gate grows
    with the table.
    """
    out = [
        f"-- {GENERATED_BANNER}",
        "--",
        "-- Generated by `tools/gen_vocabulary.py` from",
        "-- `specs/vocabulary/vocabulary.toml`.",
        "--",
        "-- One `#print axioms` per unit theorem in `Azoth/Vocabulary.lean`, so the",
        "-- gate covers every unit rather than the ones somebody remembered. Run by",
        "-- `tools/check_lean_axioms.py`, which refuses any axiom set outside",
        "-- `propext`, `Classical.choice` and `Quot.sound`.",
        "",
        "import Azoth.Vocabulary",
        "",
    ]
    out += [f"#print axioms Azoth.Vocabulary.{lean_name(u['id'])}" for u in table["units"]]
    out.append("")
    return "\n".join(out)


def emit_unit_schema(table: dict[str, Any]) -> str:
    units = [unit["id"] for unit in table["units"]]
    schema = {
        "$schema": "https://json-schema.org/draft/2020-12/schema",
        "$id": (
            "https://raw.githubusercontent.com/PatrickMockridge/Azoth/main/"
            "specs/schema/unit.schema.json"
        ),
        "title": "Azoth canonical unit",
        "description": (
            f"{GENERATED_BANNER} Generated by `tools/gen_vocabulary.py` from "
            "`specs/vocabulary/vocabulary.toml`, which is the one hand-written "
            "source of this list. A spec's `unit:` must name one of these, so a spec "
            "cannot name a unit no implementation has a conversion path for. What "
            "dimension each one carries is in the table and in "
            "`azoth_core::unit_vocab_gen`."
        ),
        "enum": units,
    }
    return json.dumps(schema, indent=2) + "\n"


def rust_i8_array(exponents: list[int]) -> str:
    return "[" + ", ".join(str(e) for e in exponents) + "]"


def rustfmt(source: str) -> str:
    """Format generated Rust with rustfmt.

    The output is committed and CI runs `cargo fmt --check` over the whole
    workspace, so generated code has to be formatted like hand-written code.
    Formatting here rather than emitting carefully aligned strings means the
    generator never has to know rustfmt's rules - and if rustfmt is unavailable
    the source is still valid, just unformatted, which `cargo fmt --check` will
    report rather than something silently drifting.
    """
    try:
        proc = subprocess.run(
            ["rustfmt", "--edition", "2024", "--emit", "stdout"],
            input=source,
            capture_output=True,
            text=True,
            check=False,
        )
    except FileNotFoundError:
        print("gen_vocabulary: rustfmt not found; emitting unformatted Rust", file=sys.stderr)
        return source
    if proc.returncode != 0:
        print(f"gen_vocabulary: rustfmt failed:\n{proc.stderr}", file=sys.stderr)
        return source
    return proc.stdout


def main() -> int:
    global VOCAB_PATH

    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "--check",
        action="store_true",
        help="do not write; fail if the generated files are out of date",
    )
    parser.add_argument(
        "--vocabulary",
        type=Path,
        default=VOCAB_PATH,
        help=(
            "the table to compile (default: specs/vocabulary/vocabulary.toml). "
            "Exists so the test suite can compile a mutated copy and check that "
            "the invariant it breaks is one that fails."
        ),
    )
    args = parser.parse_args()

    if args.vocabulary != VOCAB_PATH:
        VOCAB_PATH = args.vocabulary

    table = load_table()
    outputs: dict[Path, str] = {
        RUST_OUT: rustfmt(emit_rust(table)),
        PY_OUT: emit_python(table),
        UNIT_SCHEMA_OUT: emit_unit_schema(table),
        LEAN_OUT: emit_lean(table),
        LEAN_GATE_OUT: emit_lean_gate(table),
    }

    stale: list[Path] = []
    for path, content in outputs.items():
        path.parent.mkdir(parents=True, exist_ok=True)
        if args.check:
            current = path.read_text(encoding="utf-8") if path.exists() else ""
            if current != content:
                stale.append(path)
        else:
            path.write_text(content, encoding="utf-8")

    if args.check:
        if stale:
            for path in stale:
                print(
                    f"gen_vocabulary: {path.relative_to(ROOT)} is out of date",
                    file=sys.stderr,
                )
            print("Run `python tools/gen_vocabulary.py` to regenerate.", file=sys.stderr)
            return 1
        print(f"gen_vocabulary: {len(table['units'])} unit(s), generated files up to date")
        return 0

    for path in outputs:
        print(f"gen_vocabulary: wrote {path.relative_to(ROOT)}")
    print(f"gen_vocabulary: {len(table['units'])} unit(s)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
