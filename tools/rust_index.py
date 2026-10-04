#!/usr/bin/env python3
"""The Rust result types, as the generators that emit against them read them.

Every site the registration generators touch needs the same pair: the **public** field name - the
one in the spec's `outputs`, in `CalcResult::FIELDS` and in the Python dataclass - and the **Rust
struct's** field name, which is internal and may differ. On the `eos` flashes it does: `FIELDS`
carries a symbol where the struct spells the word out (`T` against `temperature` - or the shorter
`t`, on three of them - `P` against `pressure`, `V` against `v`), and three more are behind a
macro. A generator that emitted a
transport struct or a `From` impl from the spec alone would name a field the struct does not have.

Everything else those generators need is a function of spec data that already exists - the id, the
outputs' kinds, units and descriptions, the unit-to-constructor map. **This is the one index that
has to be read from the source**, and it is read positionally: `FIELDS` and the struct's `pub`
fields are the same list in two spellings, so pairing them by index recovers the mapping with no
exception table. A name rule would need one entry per flash, and a *spec* field naming each Rust
struct would be a second registration of what the type already declares - the thing
`docs/src/architecture/spec-files.md` forbids.

**A Python test already reads Rust source this way** (`test_registration_completeness.py`'s
`_rust_pyclasses`), and `tools/check_doc_claims.py` parses `kernels/mod.rs` and `dispatch.rs` the
same way. This module is the shared reader for the result types, so the generators do not each grow
one.
"""

from __future__ import annotations

import importlib.util
import re
import sys
import tomllib
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
CRATES = ROOT / "crates"
VOCABULARY = ROOT / "specs" / "vocabulary" / "vocabulary.toml"

_CHAR = re.compile(r"'(?:\\.|[^'\\])'")
_IMPL = re.compile(r"impl\s+CalcResult\s+for\s+([A-Za-z_][A-Za-z0-9_]*)\s*\{")
_STRUCT = re.compile(r"pub\s+struct\s+([A-Za-z_][A-Za-z0-9_]*)\s*(?:<[^{;]*>)?\s*\{")
_MACRO = re.compile(r"macro_rules!\s*([A-Za-z_][A-Za-z0-9_]*)\s*\{")
_STRUCT_TOKEN = re.compile(r"pub\s+struct\s+(\$[A-Za-z_][A-Za-z0-9_]*)\s*(?:<[^{;]*>)?\s*\{")
_CALC_ID = re.compile(r"const\s+CALC_ID\s*:\s*&'static\s+str\s*=\s*([^;]+);")
_FIELDS = re.compile(
    r"const\s+FIELDS\s*:\s*&(?:'[A-Za-z_][A-Za-z0-9_]*\s*)?\[&(?:'[A-Za-z_][A-Za-z0-9_]*\s*)?str\]"
    r"\s*=\s*&\[(.*?)\];",
    re.DOTALL,
)
_STRING = re.compile(r'"([^"]*)"')
_FIELD = re.compile(r"^\s*pub\s+([a-z_][A-Za-z0-9_]*)\s*:", re.MULTILINE)


@dataclass(frozen=True)
class ResultType:
    """One registered result: the id, the Rust type, where it lives, and its field pairing."""

    calc_id: str
    rust_name: str
    item_path: str
    module_file: Path
    fields: tuple[tuple[str, str], ...]
    rust_field_types: tuple[str, ...] = ()

    @property
    def public_fields(self) -> tuple[str, ...]:
        """The declared names, in `FIELDS` order."""
        return tuple(public for public, _ in self.fields)

    @property
    def rust_fields(self) -> tuple[str, ...]:
        """The Rust struct's own names, in declaration order."""
        return tuple(rust for _, rust in self.fields)


def _blank_comments(source: str) -> str:
    """Blank `//` and `/* */` comments at unchanged offsets, leaving literals whole.

    Offsets are kept so a later reader can point at a line number, and a doc comment is blanked
    like any other: a struct's prose names fields, and a scan that read it would pair a word from
    an explanation with a word from a declaration.
    """
    out = list(source)
    i, n = 0, len(source)
    while i < n:
        ch = source[i]
        if ch == '"':
            i += 1
            while i < n and source[i] != '"':
                i += 2 if source[i] == "\\" else 1
            i += 1
        elif ch == "'" and (char := _CHAR.match(source, i)):
            i = char.end()
        elif source.startswith("//", i):
            end = source.find("\n", i)
            end = n if end < 0 else end
            out[i:end] = " " * (end - i)
            i = end
        elif source.startswith("/*", i):
            end = source.find("*/", i + 2)
            end = n if end < 0 else end + 2
            out[i:end] = " " * (end - i)
            i = end
        else:
            i += 1
    return "".join(out)


def _block(source: str, brace: int, where: str) -> str:
    """The text between the `{` at `brace` and its matching `}`, exclusive.

    Braces inside a string literal are skipped, so a `format!` in an impl's other methods cannot
    unbalance the scan for the `const`s above it.
    """
    depth = 0
    i, n = brace, len(source)
    while i < n:
        ch = source[i]
        if ch == '"':
            i += 1
            while i < n and source[i] != '"':
                i += 2 if source[i] == "\\" else 1
            i += 1
            continue
        if ch == "'" and (char := _CHAR.match(source, i)):
            i = char.end()
            continue
        if ch == "{":
            depth += 1
        elif ch == "}":
            depth -= 1
            if depth == 0:
                return source[brace + 1 : i]
        i += 1
    sys.exit(f"rust_index: unbalanced braces from offset {brace} in {where}")


def _literal(raw: str, where: str) -> str:
    """A `CALC_ID`'s value, refused where it is not a literal this reader can resolve."""
    match = _STRING.fullmatch(raw.strip())
    if match is None:
        sys.exit(f"rust_index: {where}: CALC_ID is {raw.strip()!r}, which is not a string literal")
    return match.group(1)


def _strings(body: str) -> list[str]:
    """The string literals in a `FIELDS` array, in order."""
    return _STRING.findall(body)


def _public_fields(body: str) -> list[str]:
    """A struct body's `pub` field names, in declaration order."""
    return [match.group(1) for match in _FIELD.finditer(body)]


def _field_decls(body: str) -> list[tuple[str, str]]:
    """A struct body's `pub` fields as `(name, type)`, in declaration order.

    The type runs to the first `,` at bracket depth zero, so `Option<Pressure>` and
    `Vec<(String, f64)>` come back whole - a line-wise split would stop at the comma inside
    the brackets. Most of the tree's types are a bare `f64` or a `uom` alias; the container
    and the alias together are what a transport struct is emitted from.
    """
    out: list[tuple[str, str]] = []
    for match in _FIELD.finditer(body):
        start = match.end()
        depth = 0
        i = start
        while i < len(body):
            ch = body[i]
            if ch in "<([":
                depth += 1
            elif ch in ">)]":
                if depth == 0:
                    break
                depth -= 1
            elif ch == "," and depth == 0:
                break
            i += 1
        out.append((match.group(1), body[start:i].strip()))
    return out


def _item_path(module_file: Path, rust_name: str) -> str:
    """The path the type is reachable at, derived from the file the type is declared in."""
    relative = module_file.relative_to(CRATES)
    crate = relative.parts[0].replace("-", "_")
    inside = Path(*relative.parts[1:]).relative_to("src")
    parts = list(inside.with_suffix("").parts)
    if parts and parts[-1] == "mod":
        parts.pop()
    return "::".join([crate, *parts, rust_name])


def _parse_file(
    path: Path,
) -> tuple[dict[str, tuple[str, list[str]]], dict[str, list[tuple[str, str]]]]:
    """One file's `(impls, structs)`, with the result macros it defines expanded in place.

    `impls` maps a Rust type name to `(calc_id, FIELDS)` and `structs` to its `pub` fields as
    `(name, type)` - every `pub struct` in the file, since the two are resolved against each
    other by the caller.
    """
    source = _blank_comments(path.read_text(encoding="utf-8"))
    where = str(path.relative_to(ROOT))

    impls: dict[str, tuple[str, list[str]]] = {}
    structs: dict[str, list[tuple[str, str]]] = {}
    for match in _IMPL.finditer(source):
        name = match.group(1)
        body = _block(source, match.end() - 1, where)
        id_match = _CALC_ID.search(body)
        fields_match = _FIELDS.search(body)
        if id_match is None or fields_match is None:
            sys.exit(f"rust_index: {where}: `impl CalcResult for {name}` has no CALC_ID/FIELDS")
        impls[name] = (_literal(id_match.group(1), where), _strings(fields_match.group(1)))

    for match in _MACRO.finditer(source):
        body = _block(source, match.end() - 1, where)
        if "impl CalcResult for" not in body:
            continue
        struct_match = _STRUCT_TOKEN.search(body)
        id_match = _CALC_ID.search(body)
        fields_match = _FIELDS.search(body)
        if struct_match is None or id_match is None or fields_match is None:
            sys.exit(
                f"rust_index: {where}: macro `{match.group(1)}` declares a CalcResult without a "
                "`pub struct $name`, a `$id` and a `FIELDS` list this reader can expand"
            )
        if not id_match.group(1).strip().startswith("$"):
            sys.exit(f"rust_index: {where}: macro `{match.group(1)}`'s CALC_ID is not its `$id`")
        struct_fields = _field_decls(_block(body, struct_match.end() - 1, where))
        fields = _strings(fields_match.group(1))
        invoke = re.compile(
            rf"\b{re.escape(match.group(1))}!\s*\(\s*([A-Za-z_][A-Za-z0-9_]*)\s*,\s*\"([^\"]+)\"\s*\)\s*;"
        )
        expansions = list(invoke.finditer(source))
        if not expansions:
            sys.exit(f"rust_index: {where}: macro `{match.group(1)}` is never invoked")
        for call in expansions:
            impls[call.group(1)] = (call.group(2), list(fields))
            structs.setdefault(call.group(1), list(struct_fields))

    for match in _STRUCT.finditer(source):
        structs.setdefault(match.group(1), _field_decls(_block(source, match.end() - 1, where)))

    return impls, structs


@lru_cache(maxsize=1)
def result_types() -> tuple[ResultType, ...]:
    """Every registered Rust result type, by calc id.

    **Fails loudly rather than skipping.** A type whose `FIELDS` and struct fields do not pair
    one-to-one, or two types claiming one id, is a tree this index cannot describe - and a
    generator that shrugged would emit a transport struct with a field the struct does not have.
    """
    found: dict[str, list[ResultType]] = {}
    files = sorted(CRATES.glob("*/src/**/*.rs"))
    if not files:
        sys.exit(f"rust_index: no Rust sources under {CRATES}")

    for path in files:
        impls, structs = _parse_file(path)
        for rust_name, (calc_id, fields) in impls.items():
            rust_fields = structs.get(rust_name)
            if rust_fields is None:
                sys.exit(f"rust_index: {path}: `{rust_name}` has no `pub struct` body")
            if len(fields) != len(rust_fields):
                names = [name for name, _ in rust_fields]
                sys.exit(
                    f"rust_index: {path}: `{rust_name}` declares {len(fields)} `FIELDS` "
                    f"{fields} against {len(names)} struct field(s) {names}; they pair "
                    "by index and must be the same list in two spellings"
                )
            found.setdefault(calc_id, []).append(
                ResultType(
                    calc_id=calc_id,
                    rust_name=rust_name,
                    item_path=_item_path(path, rust_name),
                    module_file=path,
                    fields=tuple(
                        (public, rust)
                        for public, (rust, _) in zip(fields, rust_fields, strict=True)
                    ),
                    rust_field_types=tuple(rust_type for _, rust_type in rust_fields),
                )
            )

    for calc_id, types in sorted(found.items()):
        if len(types) > 1:
            names = ", ".join(result.rust_name for result in types)
            sys.exit(f"rust_index: {calc_id} is claimed by more than one result type: {names}")

    return tuple(found[key][0] for key in sorted(found))


_IMPL_BLOCK = re.compile(r"impl\s+([A-Za-z_][A-Za-z0-9_]*)\s*\{")


@lru_cache(maxsize=1)
def _string_accessors() -> dict[str, str]:
    """Each type's own accessor for the string the spec spells an enum as.

    Seven of the eight enums a result can carry expose `as_str`; `PitzerDataset` exposes
    `name`. Both are the type's own method rather than a shared trait, so the name is read
    from the source rather than guessed - a convention written into a generator would be a
    second place to change when a type picks differently.
    """
    out: dict[str, str] = {}
    for path in sorted(CRATES.glob("*/src/**/*.rs")):
        source = _blank_comments(path.read_text(encoding="utf-8"))
        where = str(path.relative_to(ROOT))
        for match in _IMPL_BLOCK.finditer(source):
            body = _block(source, match.end() - 1, where)
            for accessor in ("as_str", "name"):
                if re.search(rf"pub\s+fn\s+{accessor}\s*\(", body):
                    out.setdefault(match.group(1), accessor)
    return out


def string_accessor(type_name: str) -> str:
    """How `type_name` renders its own string, refused where it has no such method."""
    accessor = _string_accessors().get(type_name)
    if accessor is None:
        sys.exit(f"rust_index: {type_name} is a spec enum with no `as_str`/`name` accessor")
    return accessor


#: The crate a spec's `implementations.rust` names, to the module of the binding that carries
#: the `#[pyfunction]` wrapping it. Exhaustive by intent: a seventh crate is a refusal, because a
#: calculation this map could not place is one the extension does not register.
NAMESPACES = {
    "azoth_hydraulics": "hydraulics",
    "azoth_eos": "eos",
    "azoth_thermal": "thermal",
    "azoth_reactions": "reactions",
    "azoth_standards": "standards",
    "azoth_process": "process",
}


@lru_cache(maxsize=1)
def implementations() -> tuple[tuple[str, str, str], ...]:
    """Every registered id as `(calc_id, module, function)`, in id order.

    The function name is the last segment of the spec's own Rust path, which is the name the
    `#[pyfunction]` carries - written once in the spec and once in the module, with this the
    reader that asserts they are the same word rather than assuming a naming rule.
    """
    out: list[tuple[str, str, str]] = []
    for namespace in ("calcs", "models"):
        for path in sorted((ROOT / "specs" / namespace).rglob("*.toml")):
            document = tomllib.loads(path.read_text(encoding="utf-8"))
            calc_id = document["id"]
            crate, _, rest = document["implementations"]["rust"].partition("::")
            module = NAMESPACES.get(crate)
            if module is None:
                sys.exit(
                    f"rust_index: {calc_id} names the crate {crate!r}, which is not one of the "
                    f"six namespaces {sorted(NAMESPACES)}"
                )
            if not rest:
                sys.exit(f"rust_index: {calc_id}'s rust path has no item after the crate")
            out.append((calc_id, module, rest.split("::")[-1]))
    return tuple(sorted(out))


@lru_cache(maxsize=1)
def spec_ids() -> tuple[str, ...]:
    """Every id the specs register, which is what `result_types()` must account for exactly."""
    ids: list[str] = []
    for namespace in ("calcs", "models"):
        for path in sorted((ROOT / "specs" / namespace).glob("**/*.toml")):
            ids.append(tomllib.loads(path.read_text(encoding="utf-8"))["id"])
    return tuple(ids)


def check() -> None:
    """Fail if the Rust index and the specs do not name the same ids."""
    rust = {result.calc_id for result in result_types()}
    specs = set(spec_ids())
    if rust != specs:
        sys.exit(
            "rust_index: the Rust result types and the specs disagree - "
            f"in Rust and not in the specs: {sorted(rust - specs)}; "
            f"in the specs and not in Rust: {sorted(specs - rust)}"
        )


@lru_cache(maxsize=1)
def _table() -> dict[str, object]:
    return tomllib.loads(VOCABULARY.read_text(encoding="utf-8"))


@lru_cache(maxsize=1)
def _units() -> dict[str, dict[str, object]]:
    table = _table()
    return {unit["id"]: unit for unit in table["units"]}  # type: ignore[union-attr]


@lru_cache(maxsize=1)
def _dimensions() -> dict[str, tuple[int, ...]]:
    table = _table()
    return {
        dimension["id"]: tuple(dimension["exponents"])  # type: ignore[index]
        for dimension in table["dimensions"]  # type: ignore[union-attr]
    }


def _entry(unit: str) -> dict[str, object]:
    entry = _units().get(unit)
    if entry is None:
        sys.exit(f"rust_index: unit {unit!r} is not in {VOCABULARY}")
    return entry


def rust_ctor(unit: str) -> str | None:
    """The `azoth_core::units` constructor for a canonical unit, or `None` if it has none.

    A dimensionless unit has no constructor: a bare ratio is an `f64` through the conversion
    table, so there is nothing for an ascription to go through.
    """
    ctor = _entry(unit).get("rust_ctor")
    return None if ctor is None else str(ctor)


@lru_cache(maxsize=1)
def _uom_types() -> dict[tuple[int, ...], str | None]:
    """`gen_vocabulary`'s dimension-to-`uom` map, loaded by path.

    By path rather than by `import`, because a caller (a test) may not have `tools/` on its
    `sys.path` - and a second copy of the map here would be a second answer to the question the
    vocabulary generator already answers once.
    """
    path = Path(__file__).resolve().parent / "gen_vocabulary.py"
    spec = importlib.util.spec_from_file_location("gen_vocabulary_units", path)
    if spec is None or spec.loader is None:
        sys.exit(f"rust_index: cannot load {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.UOM_TYPES  # type: ignore[no-any-return]


def uom_quantity(unit: str) -> str | None:
    """The `uom::si::f64` quantity a unit's dimension carries, or `None` if uom has none."""
    exponents = _dimensions()[str(_entry(unit)["dimension"])]
    return _uom_types().get(exponents)


if __name__ == "__main__":
    check()
    print(f"rust_index: {len(result_types())} result type(s), {len(spec_ids())} spec id(s)")
