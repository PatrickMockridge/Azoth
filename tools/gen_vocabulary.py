#!/usr/bin/env python3
"""Compile the unit vocabulary table into Rust, Python and JSON Schema.

Reads `specs/vocabulary/vocabulary.yaml` and emits:

  specs/schema/unit.schema.json              the closed enum of canonical units
  crates/azoth-core/src/unit_vocab_gen.rs    UNIT_NAMES, dimensions, conversions
  python/src/azoth/core/_units_gen.py        the same, in Python
  lean/Azoth/Vocabulary.lean                 the same, as a dimension per unit

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
from pathlib import Path
from typing import Any

try:
    import yaml
except ImportError:  # pragma: no cover
    sys.exit("gen_vocabulary requires PyYAML: pip install pyyaml")

ROOT = Path(__file__).resolve().parent.parent
VOCAB_PATH = ROOT / "specs" / "vocabulary" / "vocabulary.yaml"
SCHEMA_PATH = ROOT / "specs" / "schema" / "vocabulary.schema.json"
UNIT_SCHEMA_OUT = ROOT / "specs" / "schema" / "unit.schema.json"
RUST_OUT = ROOT / "crates" / "azoth-core" / "src" / "unit_vocab_gen.rs"
PY_OUT = ROOT / "python" / "src" / "azoth" / "core" / "_units_gen.py"
LEAN_OUT = ROOT / "lean" / "Azoth" / "Vocabulary.lean"

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
    (-1, 1, -2, 0, 0, 0, 0): "Pressure",
    (-1, 1, -1, 0, 0, 0, 0): "DynamicViscosity",
    (0, 0, 0, 0, 1, 0, 0): "ThermodynamicTemperature",
    (2, 1, -3, 0, 0, 0, 0): "Power",
    (2, 0, -2, 0, -1, 0, 0): "SpecificHeatCapacity",
    (1, 1, -3, 0, -1, 0, 0): "ThermalConductivity",
    (0, 1, -3, 0, -1, 0, 0): "HeatTransfer",
    (0, 1, 0, 0, 0, -1, 0): "MolarMass",
    (3, 0, 0, 0, 0, -1, 0): "MolarVolume",
    (2, 1, -2, 0, 0, -1, 0): "MolarEnergy",
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


def snake(name: str) -> str:
    """`VolumeRate` -> `volume_rate`, which is uom's module name for it."""
    return "".join(f"_{c.lower()}" if c.isupper() else c for c in name).lstrip("_")


def load_table() -> dict[str, Any]:
    """Read and validate the table, and refuse a shape this generator cannot compile."""
    try:
        from jsonschema import Draft202012Validator
    except ImportError:  # pragma: no cover
        sys.exit("gen_vocabulary requires jsonschema: pip install jsonschema")

    table: dict[str, Any] = yaml.safe_load(VOCAB_PATH.read_text(encoding="utf-8"))
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
        if (unit["uom"] is None) != (quantity is None):
            sys.exit(
                f"gen_vocabulary: unit {unit['id']!r} declares uom {unit['uom']!r} but "
                f"dimension {unit['dimension']!r} maps to {quantity!r}. A unit's uom "
                f"path and its dimension's quantity type have to agree in both "
                f"directions, or the generated conversion and the generated assertion "
                f"would be about different things."
            )
        if unit["uom"] is not None:
            module = unit["uom"].split("::")[0]
            if module != snake(quantity):  # type: ignore[arg-type]
                sys.exit(
                    f"gen_vocabulary: unit {unit['id']!r} gives the uom path "
                    f"{unit['uom']!r}, whose module is {module!r}, but dimension "
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
        "//! `specs/vocabulary/vocabulary.yaml`.",
        "",
        "/// The base dimensions, in the order every exponent tuple here is written in.",
        "///",
        "/// uom::si::ISQ's order, and not cosmetic: the tuples below and the uom",
        "/// quantity types the table's assertions are written against are both built",
        "/// from it, so a different order permutes all of them at once.",
        f"pub const SLOTS: [&str; {len(slots)}] = [{', '.join(json.dumps(s) for s in slots)}];",
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
        "`specs/vocabulary/vocabulary.yaml`.",
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
        "-- `specs/vocabulary/vocabulary.yaml`.",
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
    ]
    out += [
        "end Azoth.Vocabulary",
        "",
    ]
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
            "`specs/vocabulary/vocabulary.yaml`, which is the one hand-written "
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
            "the table to compile (default: specs/vocabulary/vocabulary.yaml). "
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
