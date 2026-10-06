#!/usr/bin/env python3
"""Generate the Python result dataclasses.

`python/src/azoth/core/result.py` carried one frozen dataclass per registered id, and each is the
same three statements: a `CALC_ID`, one field per output in `CalcResult::FIELDS` order, and the
warnings. Nothing in one is a decision, so every line of it is a function of declarations that
already exist:

* **the fields are `FIELDS`**, in its order, which is what `test_registry_contract` holds the class
  to and what `dataclasses.fields` then gives `serialise` to walk.
* **each annotation is the transport's Rust type and the unit the spec declares** - `Q` where the
  unit is dimensioned and the bare Python type where it is not, `tuple[...]` once per vector depth,
  and `X | None` for an optional. The mapping is measured, not chosen: over the 194 classes it has
  one exception, `PitzerDataset`, which is a string at this boundary and is tabled with its reason.
* **each doc is the spec's own `description`**, because a hand paraphrase is a second copy of prose
  that already has one home. The transport generator made the same call for the Rust structs, and
  `warnings` - the one field no spec declares - gets the module's sentence.

**What a class docstring loses is commentary, not physics.** Each one opened with
``Result of ``<id>``.``, which is kept, and the paragraphs after it were notes about the Python
class. The physics is in the spec and is rendered from it; a note that is true of every result is
stated once in this module's header instead of 194 times, and one that is true of a single class
belongs in that class's spec.

    python tools/gen_python_result.py            # write
    python tools/gen_python_result.py --check    # fail if the tree is not what this emits
    python tools/gen_python_result.py --survey   # what it covers, and what it leaves
"""

from __future__ import annotations

import re
import sys
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parent))
import rust_index
from gen_python_bridge import core_signatures
from gen_python_transport import dimensioned, element_kind, parse_type, specs

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "python" / "src" / "azoth" / "core" / "result_gen.py"

#: The element type an integer width is, on this side.
INTEGERS = ("u8", "u16", "u32", "u64", "usize", "i8", "i16", "i32", "i64", "isize")

#: The element type each non-numeric kind is, named rather than derived because it is a name.
PY_ELEMENT = {"bool": "bool", "string": "str", "warning": "Warning", "kcomponent": "KComponent"}

#: The one enum that crosses as a string rather than as its own class, and the reason.
#:
#: The spec declares `dataset` an enum and `PitzerPhaseResult.dataset` is annotated `str`, so a
#: caller reads the name. It is a Python-boundary fact with no signal in the Rust type, which is
#: why it is a table of one and not a rule - and `gen_python_bridge` reads it from here rather
#: than keeping the second copy it used to.
ENUMS_THAT_STAY_STRINGS = {"PitzerDataset"}

#: The classes that are comparable with `==` where the other 193 are not.
#:
#: The module's own rule is that a result is not comparable, because equality would compare `pint`
#: quantities; `reactions.reactive_ph_flash`'s class has been the exception since the commit that
#: added the id, and **no reason for it is recorded anywhere in the tree**. Kept as it is: a
#: generated file may not quietly withdraw a behaviour a caller could be relying on.
COMPARABLE = {"ReactivePhFlashResult"}


def field_annotation(kind: str, depth: int, optional: bool, unit: object, element: str) -> str:
    """One field's Python annotation, from its transport type and the spec's unit."""
    if kind in ("number", "quantity"):
        leaf = "Q" if dimensioned(unit) else "float"
    elif kind in INTEGERS:
        leaf = "int"
    elif kind == "enum":
        leaf = PY_ELEMENT["string"] if element in ENUMS_THAT_STAY_STRINGS else element
    else:
        leaf = PY_ELEMENT[kind]
    annotation = leaf
    for _ in range(depth):
        annotation = f"tuple[{annotation}, ...]"
    return f"{annotation} | None" if optional else annotation


def field_doc(name: str, declaration: dict[str, Any]) -> str:
    """One field's `#:` line: the spec's description, or the module's sentence for `warnings`."""
    if name == "warnings":
        return "    #: Caveats."
    text = declaration.get("description") or name
    # A description is prose and may be a paragraph; a `#:` comment is one line per line.
    return "\n".join(f"    #: {line}".rstrip() for line in text.splitlines())


def python_name(calc_id: str) -> str:
    """The class name this result carries in Python, which is not always the Rust struct's.

    `eos.tp_flash_saft` answers a `SaftFlashResult` in Rust and a `TpFlashSaftResult` here, so
    the name is read from where the Python side is written down: `_core.pyi`'s return annotation
    for the id, which is the same reader `gen_python_bridge` uses to find it.
    """
    return core_signatures()[calc_id.rpartition(".")[2]][1]


def emit_class(result: rust_index.ResultType, spec: dict[str, Any]) -> tuple[str, set[str]]:
    """One dataclass, and the names it needs imported."""
    fields = list(zip([name for name, _ in result.fields], result.rust_field_types, strict=True))
    outputs = spec.get("outputs") or {}
    name = python_name(result.calc_id)
    eq = "True" if name in COMPARABLE else "False"
    lines = [
        f"@dataclass(frozen=True, slots=True, eq={eq})",
        f"class {name}(_HasWarnings):",
        f'    """Result of ``{result.calc_id}``."""',
        "",
        f'    CALC_ID: ClassVar[str] = "{result.calc_id}"',
    ]
    used: set[str] = set()
    for field, ty in fields:
        declaration = outputs.get(field) or {}
        optional, depth, element = parse_type(ty)
        kind = element_kind(element, declaration.get("type", "scalar"))
        annotation = field_annotation(kind, depth, optional, declaration.get("unit"), element)
        lines += ["", field_doc(field, declaration), f"    {field}: {annotation}"]
        # Every capitalised name in the annotation is one this module has to import: `Q`, the
        # warnings tuple's `Warning`, the fitting record, and whichever enum the field is. The
        # bare Python types are lower-case and need nothing.
        used.update(set(re.findall(r"\b[A-Z]\w*\b", annotation)) - {"None"})
    return "\n".join(lines) + "\n", used


def covered() -> tuple[list[tuple[str, str, set[str]]], list[str]]:
    """`(covered, refused)`. A result the spec does not describe is refused rather than guessed."""
    all_specs = specs()
    out: list[tuple[str, str, set[str]]] = []
    left: list[str] = []
    for result in rust_index.result_types():
        spec = all_specs.get(result.calc_id)
        if spec is None:
            left.append(f"{result.calc_id}: no spec")
            continue
        try:
            source, used = emit_class(result, spec)
        except KeyError as missing:
            left.append(f"{result.calc_id}: the transport type {missing} has no rule")
            continue
        out.append((result.calc_id, source, used))
    return out, left


def emit() -> str:
    """`result_gen.py`: every result this file can derive, in registry order."""
    out, _ = covered()
    bodies = "\n\n".join(source for _, source, _ in out)
    used: set[str] = set()
    for _, _, names in out:
        used |= names

    from_base = sorted(used - {"Q", "Warning"})
    imports = [
        "from __future__ import annotations",
        "",
        "from dataclasses import dataclass",
        "from typing import ClassVar",
        "",
        "from azoth.core.result_base import (",
    ]
    imports += [f"    {name}," for name in from_base]
    # isort puts a leading underscore after the capitals, so the mixin is named last.
    imports += ["    _HasWarnings,", ")"]
    imports += ["from azoth.core.units import Q", "from azoth.core.warnings import Warning"]

    header = (
        '"""The result dataclasses, emitted from the extension\'s own declarations.\n'
        "\n"
        "Generated by `tools/gen_python_result.py` - see its module documentation for the three\n"
        "declarations each one is a function of. **Do not edit by hand.** Several notes that used\n"
        "to sit on one class and are true of all of them:\n"
        "\n"
        "* **why `eq=False`.** `__eq__` would compare `pint` quantities, whose equality is not\n"
        "  what a numerical test wants; results are compared within a tolerance, by\n"
        "  :func:`azoth.testing.assert_results_equal`, which makes the tolerance visible.\n"
        "* **why warnings are tuples.** Immutable, and they match `Vec<Warning>`; a result whose\n"
        "  caveats can be edited after the fact is a result whose caveats can be removed.\n"
        "* **why one class per calc id.** The spec-to-result contract is asserted per id, so a\n"
        "  shared class would make it impossible to tell which calc a result came from.\n"
        "\n"
        "`python/src/azoth/core/result.py` re-exports these beside the hand-written enums, the\n"
        "`_HasWarnings` mixin and the one record that is not a result.\n"
        '"""\n'
        "\n" + "\n".join(imports) + "\n\n\n"
    )
    return header + bodies


def main() -> None:
    if "--survey" in sys.argv:
        out, left = covered()
        print(f"covered: {len(out)} of {len(out) + len(left)}")
        for why in left:
            print(f"   left hand-written: {why}")
        return
    text = emit()
    if "--check" in sys.argv:
        current = OUT.read_text(encoding="utf-8") if OUT.exists() else ""
        if current != text:
            raise SystemExit("gen_python_result: result_gen.py is out of date")
        print("gen_python_result: result_gen.py is current")
        return
    OUT.write_text(text, encoding="utf-8")
    print(f"gen_python_result: wrote {OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
