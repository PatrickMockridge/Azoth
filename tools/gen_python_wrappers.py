#!/usr/bin/env python3
"""Generate the `#[pyfunction]` wrappers whose signature is the spec's declared input list.

`crates/azoth-python/src/{hydraulics,eos,thermal,reactions,standards}.rs` carried one
hand-written wrapper per calculation, each converting SI magnitudes into the kernel's own
quantities and mapping the result into its transport. The conversion is a function of two things
already in the tree - the kernel's parameter list, and the unit the spec declares for each input -
so the wrapper is emitted rather than written.

**The coverage rule is read from the kernel, not from a list.** An id is covered when the
kernel's parameter names are the spec's declared input names; that is the plan's own "the
signature *is* the declared input list". A mixture-expanded model takes a `Mixture` the spec does
not declare, so its kernel names do not match and this file leaves it alone - which is why the
`expands_to` declaration the plan wants for those is *not* a prerequisite here.

Two things are reused rather than re-derived: a unit's constructor is `rust_index.rust_ctor`, the
vocabulary's own map, and the reading of which `pub fn` a name is comes from the crate root rather
than from a path. **The signature is in the spec's declared order and not the optional-last order**
- `gen_stub`'s `ordered_parameters` belongs to the stub, which describes a Python `def`; `_core`
is not one, and applying that rule here is what put `heat_exchanger`'s arguments in the wrong
places.

    python tools/gen_python_wrappers.py            # write
    python tools/gen_python_wrappers.py --check    # fail if the tree is not what this emits
    python tools/gen_python_wrappers.py --survey   # what it covers, and what it leaves
"""

from __future__ import annotations

import re
import subprocess
import sys
import tomllib
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import rust_index

ROOT = Path(__file__).resolve().parent.parent
CRATES = ROOT / "crates"
OUT = ROOT / "crates" / "azoth-python" / "src" / "wrappers_gen.rs"

#: Kernel parameter types that cross as themselves, with the type the binding should declare.
PASSTHROUGH = {
    "f64": "f64",
    "bool": "bool",
    "usize": "usize",
    "u32": "u32",
    "u64": "u64",
    "i64": "i64",
    "&str": "&str",
    "Option<f64>": "Option<f64>",
    "Option<bool>": "Option<bool>",
    "Option<usize>": "Option<usize>",
    "Option<u32>": "Option<u32>",
    "Option<&str>": "Option<&str>",
    "&[f64]": "Vec<f64>",
    "&[String]": "Vec<String>",
    "&[Vec<f64>]": "Vec<Vec<f64>>",
    "&[[f64; 3]]": "Vec<[f64; 3]>",
    "Option<&[f64]>": "Option<Vec<f64>>",
}

#: Kernel parameters that need a local binding before the call, and what it is.
NEEDS_LOCAL = {"&[&str]": ("Vec<String>", "Vec<&str>", "iter().map(String::as_str).collect()")}


def uom_aliases() -> set[str]:
    """The `uom` quantity aliases `azoth-core` re-exports, read from its own `pub use`.

    Read rather than listed, for the reason every other index here is: a type that is a `uom`
    quantity converts through a constructor, and one that is an enum is parsed from a string.
    Deciding which by convention would be a second place to change when the re-export does.
    """
    source = rust_index._blank_comments(
        (CRATES / "azoth-core" / "src" / "units.rs").read_text(encoding="utf-8")
    )
    match = re.search(r"pub use uom::si::f64::\{([^}]*)\}", source)
    if match is None:
        sys.exit("gen_python_wrappers: azoth-core/src/units.rs no longer re-exports uom aliases")
    return {name.strip() for name in match.group(1).split(",") if name.strip()}


def si_ctors() -> dict[str, str]:
    """Each dimension's SI base constructor, from the vocabulary's own order.

    **The binding takes SI magnitudes**, whatever unit the spec declares - `azoth.hydraulics`
    converts once in Python and hands over the base value. So the constructor is not the declared
    unit's (`millimeters` for an input declared in `mm`, which would scale it a second time) but
    the base unit's for the same dimension. The vocabulary lists the SI base first in every
    dimension, which is the property this reads rather than a name rule kept here.
    """
    table = rust_index._table()
    out: dict[str, str] = {}
    for unit in table["units"]:
        ctor = unit.get("rust_ctor")
        if ctor:
            out.setdefault(str(unit["dimension"]), str(ctor))
    return out


def root_exports(crate: str) -> dict[str, str]:
    """Each name the crate root re-exports, mapped to the module path it comes from.

    A name can be declared twice in one crate - `azoth-process` has `kernels::absorption_column`
    taking a setup record and `models::absorption_column` taking the declared inputs, and the spec
    names the second. Which one the crate root exports is the crate's own answer, so it is read
    rather than guessed at.
    """
    lib = CRATES / crate.replace("_", "-") / "src" / "lib.rs"
    if not lib.exists():
        return {}
    text = rust_index._blank_comments(lib.read_text(encoding="utf-8"))
    out: dict[str, str] = {}
    for group in re.finditer(r"^pub use ([\w:]+)::\{([^}]*)\};", text, re.M | re.S):
        for name in group.group(2).replace("\n", " ").split(","):
            if name.strip():
                out[name.strip()] = group.group(1)
    for single in re.finditer(r"^pub use ([\w:]+)::(\w+);", text, re.M):
        out[single.group(2)] = single.group(1)
    return out


def declared_in(crate: str, path: Path) -> str:
    """The module directory a file sits under, which is the segment a `pub use` names."""
    inside = path.relative_to(CRATES / crate.replace("_", "-") / "src")
    parts = list(inside.with_suffix("").parts)
    if parts and parts[-1] == "mod":
        parts.pop()
    return "::".join(parts[:-1])


def kernel_signatures() -> dict[str, tuple[list[tuple[str, str]], str, str]]:
    """Every `pub fn` the kernels declare, by name, as `(parameters, crate, item path)`.

    The crate is carried so an enum parameter can be named in full where the binding uses it.

    By name rather than by path: a spec's `implementations.rust` names the item, and which module
    it lives in is the crate's business. A name two crates both use is dropped rather than picked,
    so nothing is emitted against a signature that might be the other one.

    **`azoth-python` is not scanned.** It is this binding, and every wrapper it carries has the
    kernel's own name - so including it would make each kernel look like a duplicate of itself.
    """
    found: dict[str, list[tuple[list[tuple[str, str]], str, Path]]] = {}
    for path in sorted(CRATES.glob("*/src/**/*.rs")):
        parts = path.relative_to(CRATES).parts
        if parts[0] == "azoth-python":
            continue
        module = parts[0].replace("-", "_")
        text = rust_index._blank_comments(path.read_text(encoding="utf-8"))
        for name, params in signatures(text):
            found.setdefault(name, []).append((params, module, path))

    out: dict[str, tuple[list[tuple[str, str]], str, Path]] = {}
    exports: dict[str, set[str]] = {}
    for name, candidates in found.items():
        if len(candidates) == 1:
            out[name] = candidates[0]
            continue
        # Declared more than once in one crate: the root's own `pub use` says which is the item the
        # spec's path names. Two that both are, or neither, is a name this reader will not pick.
        for _, module, _ in candidates:
            exports.setdefault(module, root_exports(module))
        matched = [
            candidate
            for candidate in candidates
            if exports[candidate[1]].get(name) == declared_in(candidate[1], candidate[2])
        ]
        if len(matched) == 1:
            out[name] = matched[0]
    return out


def signatures(text: str, visibility: str = r"pub fn") -> list[tuple[str, list[tuple[str, str]]]]:
    """Every top-level visible `fn name(...)` in one file, as `(name, parameters)`.

    `visibility` is the prefix to match - `pub fn` for the kernels, and an optional `pub(crate)`
    for the binding's own helpers, which the mixture builder is.
    """
    out = []
    for match in re.finditer(rf"^{visibility} (\w+)\(", text, re.M):
        depth, i = 0, match.end() - 1
        while i < len(text):
            if text[i] == "(":
                depth += 1
            elif text[i] == ")":
                depth -= 1
                if depth == 0:
                    break
            i += 1
        inner, parts, depth, start = text[match.end() : i], [], 0, 0
        for j, ch in enumerate(inner):
            if ch in "<([":
                depth += 1
            elif ch in ">)]":
                depth -= 1
            elif ch == "," and depth == 0:
                parts.append(inner[start:j])
                start = j + 1
        parts.append(inner[start:])
        params: list[tuple[str, str]] = []
        for part in parts:
            found = re.match(r"(?:mut\s+)?(\w+)\s*:\s*(.+)$", " ".join(part.split()))
            if found:
                params.append((found.group(1), found.group(2)))
        out.append((match.group(1), params))
    return out


def type_path(name: str, source: str, crate: str, file: Path) -> str:
    """Where a kernel names a type it takes, read from that file's own `use`.

    A parameter's type is spelled short in the kernel - `AntoineForm`, `Cubic` - and pyo3 needs
    the path. The kernel's `use` line is where the path is written, so it is read rather than
    guessed: `Cubic` comes from `azoth_eos` and `LiquidViscosityLadder` from a submodule, and a
    rule that prefixed the kernel's own crate would name the wrong one for both.
    """
    for match in re.finditer(r"^use\s+([^;]+);", source, re.M):
        path = " ".join(match.group(1).split())
        if path.endswith("}") and "{" in path:
            base, _, inner = path.partition("::{")
            if re.search(rf"\b{re.escape(name)}\b", inner):
                return f"{_root(base, crate)}::{name}"
            continue
        if path.rsplit("::", 1)[-1] == name:
            return f"{_root(path.rsplit('::', 1)[0], crate)}::{name}"
    if re.search(rf"^pub (?:enum|struct) {re.escape(name)}\b", source, re.M):
        # Defined in the kernel's own file: the path is that file's, not a `use`.
        return rust_index._item_path(file, name)
    return f"{crate}::{name}"


def from_str_error(name: str) -> str:
    """How a `FromStr` failure becomes a `PyErr`: `to_pyerr` for `AzothError`, else `PyValueError`.

    Read from the enum's own `impl FromStr`, because the two are not interchangeable - a
    `Result<_, AzothError>` does not convert to a Python object, and a `&'static str` is not an
    `AzothError`. Which one an enum carries is a fact about that enum.
    """
    last = name.rsplit("::", 1)[-1]
    for path in sorted(CRATES.glob("*/src/**/*.rs")):
        text = rust_index._blank_comments(path.read_text(encoding="utf-8"))
        match = re.search(
            rf"^impl\s+(?:std::str::)?FromStr\s+for\s+{re.escape(last)}\s*\{{", text, re.M
        )
        if match is None:
            continue
        body = rust_index._block(text, match.end() - 1, str(path))
        found = re.search(r"type\s+Err\s*=\s*([\w:]+)", body)
        if found is None:
            continue
        error = found.group(1).rsplit("::", 1)[-1]
        return "to_pyerr" if error == "AzothError" else "value"
    return "value"


def _root(base: str, crate: str) -> str:
    """A `use` path's root, with `crate::` resolved to the package the binding names."""
    if base == "crate":
        return crate
    if base.startswith("crate::"):
        return base.replace("crate", crate, 1)
    if base.startswith(("self::", "super::")):
        return f"{crate}::{base.split('::', 1)[1]}"
    return base


def specs_by_id() -> dict[str, dict]:
    out: dict[str, dict] = {}
    for namespace in ("calcs", "models"):
        for path in sorted((ROOT / "specs" / namespace).rglob("*.toml")):
            document = tomllib.loads(path.read_text(encoding="utf-8"))
            out[document["id"]] = document
    return out


def base_ctor(unit: str | None) -> str | None:
    """The SI base unit's constructor for a declared unit's dimension, or `None`."""
    if not unit:
        return None
    dimension = rust_index._units()[unit]["dimension"] if unit in rust_index._units() else None
    return si_ctors().get(str(dimension))


#: How a struct-typed kernel parameter's own fields cross. A record the spec does not declare is
#: transported by reading its fields, the same reading `rust_index` makes of a result type - so a
#: kernel wanting an `IdealGasModel` is handed `cp_a`-`cp_e` and this assembles the struct.
RECORD_FIELDS = {
    "Vec<f64>": "Vec<f64>",
    "Vec<Vec<f64>>": "Vec<Vec<f64>>",
    "f64": "f64",
    "bool": "bool",
    "usize": "usize",
    "String": "String",
    "Vec<String>": "Vec<String>",
}

#: How a parameter of the binding's own mixture builder crosses: the type pyo3 declares, and the
#: expression the call passes. `build_mixture` *is* the boundary - a kernel wanting a `Mixture`
#: takes eight arguments here, and the builder's signature is where that list is written.
BUILDER_TYPES = {
    "&[f64]": ("Vec<f64>", "{name}"),
    "Vec<f64>": ("Vec<f64>", "{name}"),
    "&PyAssociationSpec": ("PyRef<'_, PyAssociationSpec>", "{name}"),
    "&str": ("&str", "{name}"),
    "Option<&[Vec<f64>]>": ("Option<Vec<Vec<f64>>>", "NAME.as_deref()"),
    "&[f32]": ("Vec<f32>", "{name}"),
}


#: The binding's mixture builders, tried in order. A kernel takes the expansion one of them
#: defines; which one is decided by the emitted arity matching the bridge's call, not guessed.
BUILDERS = ("build_mixture", "build_mixture_with_mass")


def builder_parameters(name: str = "build_mixture") -> list[tuple[str, str, str]]:
    """`build_mixture`'s parameters as `(name, declared type, call expression)`.

    Read from the builder rather than declared here, because the builder is the boundary: a kernel
    that wants a `Mixture` is handed these eight arguments, and the list moves when the builder's
    does. The call expression is `&name` for a slice the builder borrows, `name` for one it takes
    by value, and the `Option` unwrap for the alpha parameters.
    """
    builder = _BUILDER
    source = rust_index._blank_comments(
        (CRATES / "azoth-python" / "src" / "eos.rs").read_text(encoding="utf-8")
    )
    for name, params in signatures(source, visibility=r"(?:pub(?:\s*\([^)]*\))?\s+)?fn"):
        if name != builder:
            continue
        out = []
        for parameter, kind in params:
            if parameter == "py":
                continue
            entry = BUILDER_TYPES.get(kind)
            if entry is None:
                sys.exit(f"gen_python_wrappers: build_mixture's {parameter}: {kind} has no rule")
            declared, call = entry
            call = call.format(name=parameter).replace("NAME", parameter)
            if kind in ("&[f64]", "&PyAssociationSpec"):
                call = f"&{parameter}"
            out.append((parameter, declared, call))
        return out
    sys.exit(f"gen_python_wrappers: eos.rs no longer declares {builder}")


#: The mixture builder's cubic choice: the three parameters the boundary gives Python defaults,
#: and the only part of the bundle that is optional. Named rather than inferred because the
#: builder itself declares no defaults - they are the `Mixture`'s, which is what `boundary_defaults`
#: reads the values from.
CUBIC_CHOICE = ("eos", "alpha", "alpha_params")


def boundary_defaults() -> dict[str, str]:
    """The defaults the cubic choice carries, read from the Python `Mixture` it mirrors.

    `Mixture.cubic` defaults to `PR` and `Mixture.alpha` to `"pr"`, and the boundary takes the
    same two defaults - so the value is read from the dataclass rather than restated here. A
    cubic's boundary spelling is its `name`, which is what the bridge passes.
    """
    mixture = (ROOT / "python" / "src" / "azoth" / "eos" / "mixture.py").read_text(encoding="utf-8")
    cubic_source = (ROOT / "python" / "src" / "azoth" / "eos" / "cubic.py").read_text(
        encoding="utf-8"
    )
    cubic_default = re.search(r"cubic: Cubic = field\(default=(\w+)\)", mixture)
    alpha_default = re.search(r'alpha: str = field\(default="([^"]*)"\)', mixture)
    if cubic_default is None or alpha_default is None:
        sys.exit(
            "gen_python_wrappers: python/src/azoth/eos/mixture.py no longer defaults cubic/alpha"
        )
    named = re.search(
        rf'^{re.escape(cubic_default.group(1))} = Cubic\(\s*name="([^"]*)"', cubic_source, re.M
    )
    if named is None:
        sys.exit(
            f"gen_python_wrappers: {cubic_default.group(1)} in eos/cubic.py carries no name this "
            "reader can resolve"
        )
    return {
        "eos": f'"{named.group(1)}"',
        "alpha": f'"{alpha_default.group(1)}"',
        "alpha_params": "None",
    }


def bridge_call(name: str) -> int | None:
    """How many arguments the bridge passes to `_core.<name>`, or `None` if it does not call it.

    **The bridge is the only caller, so the boundary has to agree with it.** A generated wrapper
    that took a different number of arguments would not fail to compile - it would fail at the
    first call, which is what the arity check here turns into a refusal at generation time.
    """
    source = (ROOT / "python" / "src" / "azoth" / "_rust_bridge.py").read_text(encoding="utf-8")
    # The call's own comments carry commas and would be counted as arguments, so each line's
    # comment is cut first - the reader is counting arguments, not reading prose.
    source = "\n".join(
        line[: line.index("#")]
        if "#" in line and line.count('"', 0, line.index("#")) % 2 == 0
        else line
        for line in source.split("\n")
    )
    match = re.search(rf"_core\.{re.escape(name)}\(", source)
    if match is None:
        return None
    depth, i = 0, match.end() - 1
    while i < len(source):
        if source[i] == "(":
            depth += 1
        elif source[i] == ")":
            depth -= 1
            if depth == 0:
                break
        i += 1
    inner, count, deep, start = source[match.end() : i], 0, 0, 0
    for j, ch in enumerate(inner):
        if ch in "([{":
            deep += 1
        elif ch in ")]}":
            deep -= 1
        elif ch == "," and deep == 0:
            count += 1
            start = j + 1
    if inner[start:].strip():
        count += 1
    return count


def record_fields(type_name: str) -> list[tuple[str, str]] | None:
    """A record's `pub` fields as `(name, type)`, or `None` where no crate declares it."""
    for path in sorted(CRATES.glob("*/src/**/*.rs")):
        text = rust_index._blank_comments(path.read_text(encoding="utf-8"))
        match = re.search(rf"^pub struct {re.escape(type_name)}\s*(?:<[^{{;]*>)?\s*\{{", text, re.M)
        if match is None:
            continue
        return rust_index._field_decls(rust_index._block(text, match.end() - 1, str(path)))
    return None


def bridge_arguments(function: str) -> list[str]:
    """The bridge's call to `_core.<function>`, as its top-level argument expressions."""
    source = (ROOT / "python" / "src" / "azoth" / "_rust_bridge.py").read_text(encoding="utf-8")
    source = "\n".join(
        line[: line.index("#")]
        if "#" in line and line.count('"', 0, line.index("#")) % 2 == 0
        else line
        for line in source.split("\n")
    )
    match = re.search(rf"_core\.{re.escape(function)}\(", source)
    if match is None:
        return []
    depth, i = 0, match.end() - 1
    while i < len(source):
        if source[i] == "(":
            depth += 1
        elif source[i] == ")":
            depth -= 1
            if depth == 0:
                break
        i += 1
    inner, parts, depth, start = source[match.end() : i], [], 0, 0
    for j, ch in enumerate(inner):
        if ch in "([{":
            depth += 1
        elif ch in ")]}":
            depth -= 1
        elif ch == "," and depth == 0:
            parts.append(inner[start:j])
            start = j + 1
    parts.append(inner[start:])
    return [part.strip() for part in parts if part.strip()]


def bridge_order(function: str, names: list[str]) -> bool:
    """Whether the bridge passes `names` in the order they are declared.

    **The bridge is the only caller and it calls positionally**, so the order it passes is the
    contract - and it is not the spec's for every id: `heat_exchanger`'s hand-written wrapper took
    the two component lists first and the fields after, which the spec does not say. Each argument
    is matched to the one declared name it mentions, and a sequence that is not the declared order
    is a signature this file would emit in the wrong places.
    """
    carried = set(names)
    seen: list[str] = []
    for argument in bridge_arguments(function):
        mentioned = [name for name in names if re.search(rf"\b{re.escape(name)}\b", argument)]
        if len(mentioned) == 1:
            seen.append(mentioned[0])
    return [name for name in names if name in carried and name in seen] == seen


def pairing(kernel: list[tuple[str, str]], inputs: dict) -> dict[str, str]:
    """Each kernel parameter to the declared input it carries, or `{}` where that cannot be read.

    **By name, case-insensitively.** A `uom` quantity is a `t` or a `p` by the Rust convention
    where the spec writes the chemistry's `T` and `P`, and the call is positional - so the case is
    the kernel's to choose rather than a second name to declare.

    **Then the single leftover.** Sixteen kernels call the component-name list `names`, `params` or
    `coeffs` where the spec declares `components`, and two call the temperature `temperature`. When
    exactly one name is left unmatched on each side they are the same parameter, and pairing them
    is a reading of the two lists rather than a rename table kept here. Two or more leftovers are a
    refusal: that shape is a re-ordered or a genuinely different signature.
    """
    out: dict[str, str] = {}
    # A parameter the boundary transports *structurally* is not paired with a declared input:
    # the mixture is the spec's `components` read as arrays, and a record is a type the spec does
    # not name at all. Neither is a leftover, and the mixture consumes the components input.
    structural = {
        name
        for name, kind in kernel
        if "Mixture" in kind
        or (
            kind.startswith("&")
            and kind[1:] not in uom_aliases()
            and not kind.startswith(("&[", "&str"))
        )
    }
    declared = dict(inputs)
    if any("Mixture" in kind for _, kind in kernel):
        for spec_name, decl in inputs.items():
            if decl.get("type") == "components":
                declared.pop(spec_name, None)
                break
    casefold = {name.lower(): (name, kind) for name, kind in kernel if name not in structural}
    leftovers = []
    for spec_name in declared:
        found = casefold.pop(spec_name.lower(), None)
        if found is None:
            leftovers.append(spec_name)
        else:
            out[found[0]] = spec_name
    remaining = [name for name, _ in casefold.values()]
    if not remaining and not leftovers:
        return out
    if len(remaining) == 1 and len(leftovers) == 1:
        out[remaining[0]] = leftovers[0]
        return out
    return {}


class Refusal(Exception):
    """This file cannot derive the wrapper, so it stays hand-written."""


def convert(kind: str, unit: str | None, name: str, aliases: set[str]) -> tuple[str, str]:
    """`(declared type, argument)` for one kernel parameter, or a refusal.

    The kernel's own type decides the conversion: a `uom` alias goes through the vocabulary's
    constructor for the unit the spec declares, a bare `f64` crosses as it is, and an enum is
    parsed from the `&str` pyo3 hands over.
    """
    if kind in PASSTHROUGH:
        declared = PASSTHROUGH[kind]
        # A slice crosses as an owned `Vec` and the kernel borrows it; an optional one borrows a
        # slice of the `Vec` the caller gave, which is what `as_deref` is for.
        if kind.startswith("Option<&["):
            return declared, f"{name}.as_deref()"
        return declared, (f"&{name}" if declared.startswith("Vec<") else name)
    if kind in NEEDS_LOCAL:
        return NEEDS_LOCAL[kind][0], f"&{name}_refs"
    if kind.startswith("Option<") and kind.endswith(">"):
        inner = kind[len("Option<") : -1]
        if inner in aliases:
            ctor = base_ctor(unit)
            if ctor is None:
                raise Refusal(f"Option<{inner}> has no constructor for unit {unit!r}")
            return "Option<f64>", f"{name}.map({ctor})"
        raise Refusal(f"Option<{inner}> is neither a uom alias nor a passthrough")
    if kind.startswith("&[") and kind.endswith("]"):
        inner = kind[2:-1]
        if inner not in aliases:
            raise Refusal(f"&[{inner}] is not a uom alias")
        ctor = base_ctor(unit)
        if ctor is None:
            raise Refusal(f"&[{inner}] has no constructor for unit {unit!r}")
        return "Vec<f64>", f"&{name}.iter().map(|v| {ctor}(*v)).collect::<Vec<_>>()"
    if kind in aliases:
        ctor = base_ctor(unit)
        if ctor is None:
            raise Refusal(f"{kind} has no base constructor for unit {unit!r}")
        if rust_index.uom_quantity(str(unit)) != kind:
            # A unit's constructor yields the quantity its *dimension* names, and a kernel that
            # wants a different one - `TemperatureInterval` for a `dT` declared in `K` - needs a
            # conversion this file has no rule for.
            raise Refusal(
                f"unit {unit!r} constructs {rust_index.uom_quantity(str(unit))}, not {kind}"
            )
        return "f64", f"{ctor}({name})"
    if any(ch in kind for ch in "&<"):
        # A reference or a generic this file has no rule for - `&GeNrtlPhaseParameters` is a
        # parameter record the spec does not declare, the same boundary the mixture sits on.
        raise Refusal(f"{kind} is neither a `uom` alias nor an enum this file can parse")
    # A bare capitalised name that is not a `uom` alias is an enum: parsed from the string pyo3
    # hands over.
    return "&str", name


def locals_for(kind: str, name: str) -> list[str]:
    """The `let` statements a parameter needs, before the call."""
    if kind == "&[&str]":
        return [f"    let {name}_refs: Vec<&str> = {name}.iter().map(String::as_str).collect();"]
    return []


def _parse_err(error: str) -> str:
    """The closure a parse failure goes through, from the enum's own `FromStr`."""
    return (
        "|e| to_pyerr(py, e)" if error == "to_pyerr" else "pyo3::exceptions::PyValueError::new_err"
    )


def enum_local(kind: str, name: str, resolved: str, error: str) -> tuple[list[str], str]:
    """An enum parameter, parsed into a local. `(statements, argument)`."""
    if kind.startswith("Option<"):
        return (
            [
                f"    let {name}_parsed: Option<{resolved}> = {name}",
                f"        .map(|value| value.parse().map_err({_parse_err(error)}))",
                "        .transpose()?;",
            ],
            f"{name}_parsed",
        )
    return (
        [
            f"    let {name}_parsed: {resolved} = {name}",
            "        .parse()",
            f"        .map_err({_parse_err(error)})?;",
        ],
        f"{name}_parsed",
    )


def body(
    entry,
    kernel,
    crate: str,
    path: str,
    source: str,
    file: Path,
    aliases: set[str],
    builder: str = "build_mixture",
) -> str:
    global _BUILDER
    _BUILDER = builder
    calc_id, _, function, spec = entry
    inputs = spec["inputs"]
    # **The kernel's own spelling, read case-insensitively.** A `uom` quantity is a `t` or a `p`
    # by the Rust convention where the spec writes the chemistry's `T` and `P`, and the call is
    # positional - so the case is the kernel's to choose and not a second name to declare.
    pairs = pairing(kernel, inputs)
    if not pairs:
        raise Refusal("the kernel's parameter names are not the spec's")

    # What each parameter is, and whether anything about it was refused.
    mixture = next((name for name, kind in kernel if "Mixture" in kind), None)
    #: A struct-typed parameter's own fields, in the order the record declares them.
    record: list[str] = []
    plan: list[tuple[str, str, str, bool]] = []  # (parameter, kind, argument, is_enum)
    for name, kind in kernel:
        if "Mixture" in kind:
            # The builder is bound before the call, so the argument names that binding.
            plan.append(("mixture", kind, "&mixture", False))
            continue
        if kind.startswith("&") and kind[1:] not in aliases and not kind.startswith(("&[", "&str")):
            fields = record_fields(kind[1:])
            if fields is None:
                raise Refusal(f"{kind} is a reference to a type no crate declares")
            for field, field_type in fields:
                declared = RECORD_FIELDS.get(field_type)
                if declared is None:
                    raise Refusal(f"{kind[1:]}'s field {field}: {field_type} has no rule")
                record.append(field)
            plan.append(
                (
                    name,
                    kind,
                    f"&{type_path(kind[1:], source, crate, file)} {{ {', '.join(record)} }}",
                    False,
                )
            )
            continue
        if "Mixture" in kind:
            # The builder is bound before the call, so the argument names that binding.
            plan.append(("mixture", kind, "&mixture", False))
            continue
        is_enum = (
            kind not in PASSTHROUGH
            and kind not in aliases
            and not kind.startswith(("&[", "Option<"))
        )
        if kind.startswith("Option<") and kind.endswith(">"):
            is_enum = kind[len("Option<") : -1] not in aliases and kind not in PASSTHROUGH
        public = pairs.get(name)
        if public is None:
            raise Refusal(f"kernel parameter {name!r} names no declared input")
        try:
            # The *declared* name is what the emitted function binds, so the argument expression
            # reads that; the kernel's own name only decided the kind and the conversion.
            _declared, argument = convert(kind, inputs[public].get("unit"), public, aliases)
        except Refusal as refusal:
            raise Refusal(str(refusal)) from None
        plan.append((public, kind, argument if not is_enum else public, is_enum))

    statements: list[str] = []
    arguments: list[str] = []
    if mixture is not None:
        call = ", ".join(expression for _, _, expression in builder_parameters(builder))
        statements.append(f"    let mixture = {builder}(py, {call})?;")
    for name, kind, argument, is_enum in plan:
        statements.extend(locals_for(kind, name))
        if is_enum:
            inner = kind[len("Option<") : -1] if kind.startswith("Option<") else kind
            # A kernel may write the type inline as `crate::chemical_equilibrium::X` rather than
            # importing it, so a path is resolved by its root and a bare name by the `use` line.
            resolved = (
                _root(inner, crate) if "::" in inner else type_path(inner, source, crate, file)
            )
            parsed, argument = enum_local(kind, name, resolved, from_str_error(resolved))
            statements.extend(parsed)
        arguments.append(argument)

    # **The signature is in the spec's declared order, not the optional-last order.** `_core` is
    # the bridge's private door and the bridge calls it positionally in this order; the
    # optional-last rule belongs to the *stub*, which describes the public function, and `gen_stub`
    # applies it there. `_core`'s signature is not a Python `def`, so it is under no such
    # constraint - and inventing one here is what put `heat_exchanger`'s arguments in the wrong
    # places.
    declared_order = list(inputs.items())
    declarations = []
    defaults = boundary_defaults() if mixture is not None else {}
    if mixture is not None:
        declarations.extend(
            f"{parameter}: {declared}"
            for parameter, declared, _ in builder_parameters(builder)
            if parameter not in CUBIC_CHOICE
        )
    if record:
        for field in record:
            found = next(
                (t for n, nk in kernel for f, t in (record_fields(nk[1:]) or []) if f == field),
                None,
            )
            declarations.append(f"{field}: {RECORD_FIELDS[found]}")
    for name, decl in declared_order:
        kind = next((kind for key, kind in kernel if pairs.get(key) == name), None)
        if kind is not None and decl.get("optional") and not kind.startswith("Option<"):
            # The spec lets the caller omit it and the kernel wants a value, so a default has to
            # be chosen - which is a decision this file will not invent.
            raise Refusal(f"{calc_id}: {name} is optional in the spec and the kernel takes {kind}")
        if kind is None:
            # An input the boundary carries structurally - the mixture's `components`, which the
            # builder's arrays stand for - has no parameter of its own.
            continue
        if "Mixture" in kind:
            # The mixture's own input is carried by the builder's bundle, not by a parameter of
            # the calc's name - the spec declares a component *list* and the boundary takes arrays.
            continue
        declared, _ = convert(kind, decl.get("unit"), name, aliases)
        if (
            decl.get("optional")
            and not declared.startswith("Option<")
            and declared != "&str"
            and not declared.startswith("Vec<")
        ):
            declared = f"Option<{declared}>"
        declarations.append(f"{name}: {declared}")
    if mixture is not None:
        declarations.extend(
            f"{parameter}: {declared}"
            for parameter, declared, _ in builder_parameters(builder)
            if parameter in CUBIC_CHOICE
        )

    rust_name = next(r.rust_name for r in rust_index.result_types() if r.calc_id == calc_id)
    names = [name for name, _ in declared_order]
    if mixture is not None:
        # The mixture's declared input is not a parameter: the bundle carries it.
        paired = set(pairs.values())
        names = [name for name in names if name in paired]
        bundle = [parameter for parameter, _, _ in builder_parameters(builder)]
        names = (
            [n for n in bundle if n not in CUBIC_CHOICE]
            + names
            + [n for n in bundle if n in CUBIC_CHOICE]
        )
    if record:
        # The record's fields sit between the mixture's required half and the calc's own inputs.
        head = 0
        if mixture is not None:
            bundle = [parameter for parameter, _, _ in builder_parameters(builder)]
            head = len([n for n in bundle if n not in CUBIC_CHOICE])
        names = names[:head] + record + names[head:]
    signature_list = ", ".join(f"{n} = {defaults[n]}" if n in defaults else n for n in names)
    lines = [
        f"/// {spec.get('name', function)}.",
        "///",
        "/// All arguments are SI magnitudes. See the module documentation for why.",
        "#[pyfunction]",
        f"#[pyo3(signature = ({signature_list}))]",
        f'#[pyo3(text_signature = "({signature_list.replace(chr(34), chr(92) + chr(34))})")]',
    ]
    if any(name != name.lower() for name in names):
        lines.append("#[allow(non_snake_case)] // symbols from the published equation")
    if len(declarations) > 6:
        lines.append(
            "#[allow(clippy::too_many_arguments)] // The signature is the spec's declared inputs."
        )
    lines.append(f"pub fn {function}(")
    lines.append("    py: Python<'_>,")
    lines.extend(f"    {declaration}," for declaration in declarations)
    lines.append(f") -> PyResult<crate::transport_gen::Py{rust_name}> {{")
    lines.extend(statements)
    lines.append(f"    {path}({', '.join(arguments)})")
    lines.append(f"        .map(|r| crate::transport_gen::Py{rust_name}::from(&r))")
    lines.append("        .map_err(|e| to_pyerr(py, e))")
    lines.append("}")
    return "\n".join(lines)


def covered() -> list[tuple[tuple[str, str, str, dict], str]]:
    """Every id this file can emit, with its rendered wrapper."""
    kernels = kernel_signatures()
    aliases = uom_aliases()
    specs = specs_by_id()
    out = []
    for calc_id, module, function in rust_index.implementations():
        found = kernels.get(function)
        spec = specs[calc_id]
        entry = (calc_id, module, function, spec)
        if found is None:
            continue
        kernel, kernel_module, kernel_file = found
        kernel_path = rust_index._item_path(kernel_file, function)
        if not pairing(kernel, spec["inputs"]):
            continue
        source = rust_index._blank_comments(kernel_file.read_text(encoding="utf-8"))
        called = bridge_call(function)
        # **The builder is searched, not assumed.** A kernel wants *a* mixture expansion, and
        # which one is decided by the emitted arity matching the only caller's call - that is
        # what tells `eos.viscosity` (which needs `molar_mass` at the boundary too) from
        # `eos.pt_flash`. No builder that fits means a refusal, not a guess.
        choices = BUILDERS if any("Mixture" in kind for _, kind in kernel) else ("build_mixture",)
        text = None
        for builder in choices:
            try:
                candidate = body(
                    entry, kernel, kernel_module, kernel_path, source, kernel_file, aliases, builder
                )
            except Refusal:
                continue
            # `py` is injected by pyo3 and is not an argument the bridge passes.
            declared = len(re.findall(r"^    (?:\w+): ", candidate, re.M)) - 1
            if called is None or called == declared:
                text = candidate
                break
        if text is None:
            continue
        if any("Mixture" not in kind for _, kind in kernel):
            # The mixture's own parameters are the bundle's, in the bundle's order; the rest are
            # the spec's, and the bridge has to pass them in the order the spec declares them.
            carried = set(pairing(kernel, spec["inputs"]).values())
            if not bridge_order(function, [name for name in spec["inputs"] if name in carried]):
                continue
        out.append((entry, text))
    return out


def ctors(text: str) -> list[str]:
    """The vocabulary constructors the emitted wrappers name, in call or in `map`.

    A word scan rather than a call scan: a constructor reaches the call as `pascals(P)` but also as
    `duty.map(watts)`, where the name is followed by a bracket rather than a letter.
    """
    return sorted(
        name
        for name in _ALL_CTORS
        # A constructor reaches the call either directly - `pascals(P)` - or through a `map`,
        # where it is the whole argument: `duty.map(watts)`. Matched in those two positions
        # rather than as a bare word, because a parameter can be *named* after a unit.
        if re.search(rf"(?:^|[^\w.]){name}\(", text, re.M)
        or re.search(rf"\.map\(\s*{name}\s*\)", text)
    )


def _all_ctors() -> set[str]:
    table = rust_index._table()
    return {str(u["rust_ctor"]) for u in table["units"] if u.get("rust_ctor")}


_ALL_CTORS = _all_ctors()


def covered_ids() -> set[tuple[str, str]]:
    """The `(module, function)` pairs this file emits, which `gen_python_register` points at."""
    return {(entry[1], entry[2]) for entry, _ in covered()}


def emit() -> str:
    bodies = "\n\n".join(text for _, text in covered())
    used = ctors(bodies)
    imports = "use crate::errors::to_pyerr;\n"
    if "build_mixture(" in bodies or "build_mixture_with_mass(" in bodies:
        # The mixture boundary the expanded wrappers share: the builder and the record it takes.
        builders = sorted({b for b in BUILDERS if f"{b}(" in bodies})
        imports += f"use crate::eos::{{PyAssociationSpec, {', '.join(builders)}}};\n"
    if used:
        imports += f"use azoth_core::units::{{{', '.join(used)}}};\n"
    imports += "use pyo3::prelude::*;\n"
    header = (
        "//! GENERATED FILE - DO NOT EDIT BY HAND.\n"
        "//!\n"
        "//! Generated by `tools/gen_python_wrappers.py` from each spec's declared inputs\n"
        "//! and the kernel's own signature, which agree by name on every id this file covers.\n"
        "//!\n"
        "//! Regenerate with `python tools/gen_python_wrappers.py`; CI runs `--check`\n"
        "//! and fails on any difference.\n"
        "\n"
    )
    return header + imports + "\n" + bodies + "\n"


def format_rust(text: str) -> str:
    """Run the emitted Rust through rustfmt, which the committed file is compared against."""
    try:
        proc = subprocess.run(
            ["rustfmt", "--edition", "2024"],
            input=text,
            capture_output=True,
            text=True,
            check=False,
        )
    except FileNotFoundError:
        sys.exit(
            "gen_python_wrappers: rustfmt is not on PATH, and the generated file is compared "
            "against rustfmt's own output - install the Rust toolchain"
        )
    if proc.returncode != 0:
        sys.exit(f"gen_python_wrappers: rustfmt refused the output:\n{proc.stderr}")
    return proc.stdout


def main() -> None:
    if "--survey" in sys.argv:
        taken = covered()
        total = len(rust_index.implementations())
        print(f"covered: {len(taken)} of {total}")
        names = {cid for (cid, _, _, _), _ in taken}
        for calc_id, _, _ in rust_index.implementations():
            if calc_id not in names:
                print(f"   left hand-written: {calc_id}")
        return
    check = "--check" in sys.argv
    rendered = format_rust(emit())
    existing = OUT.read_text(encoding="utf-8") if OUT.exists() else None
    if existing == rendered:
        print(f"gen_python_wrappers: {OUT.relative_to(ROOT)} is current")
        return
    if check:
        sys.exit(
            f"gen_python_wrappers: {OUT.relative_to(ROOT)} is out of date; run "
            "tools/gen_python_wrappers.py"
        )
    OUT.write_text(rendered, encoding="utf-8")
    print(f"gen_python_wrappers: wrote {OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
