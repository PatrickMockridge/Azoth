#!/usr/bin/env python3
"""Generate the domain each card parameter has to lie in.

`docs/src/calculus/values.md` states the claim this exists for: *a datum a bank holds is inside the
domain the layer declares for it, or the lookup is refused*. Until this generator there was no
domain declared anywhere - `specs/schema/component.schema.json` carried a parameter as
`{value, unit, citation}` with no bound, and both readers checked a value's unit and its dimension
and never its magnitude, so a card could state `Tc = -1e9 K` and the flash would run on it.

**The schema is the one declaration.** Each numeric parameter carries either a bound
(`x-azoth-min`, `x-azoth-max`, the two `-inclusive` flags and `x-azoth-rationale`) or
`x-azoth-unbounded` with a reason, and a parameter carrying neither is **refused** rather than
skipped - so the table stays total as the format grows, and "nothing checks this" is a recorded
decision rather than an omission.

**A bound is stated in the parameter's own declared unit**, beside `x-azoth-unit`, because that is
the unit the schema reader is looking at. Each reader converts the value to the same unit before
comparing, which is why the bound is not emitted in SI: a generator converting it here would have
to write a conversion factor down, and `gen_vocabulary.py` refuses to emit one for the reason
`docs/src/calculus/vocabulary.md` gives - a number written in a third place is a number that can
disagree with both of the ones that know it.

**Why a bound is not a proof.** `Tc > 0` is a statement about the domain the arithmetic requires,
not about the number being right. `Tc = 190.564 K` is empirical and nothing here bears on it. The
two are kept apart everywhere, including on the page this serves.

    python tools/gen_parameter_bounds.py            # write
    python tools/gen_parameter_bounds.py --check    # fail if the tree is not what this emits
    python tools/gen_parameter_bounds.py --survey   # what each parameter's decision is
"""

from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent.parent
SCHEMA = ROOT / "specs" / "schema" / "component.schema.json"
PY_OUT = ROOT / "python" / "src" / "azoth" / "core" / "_bounds_gen.py"
RUST_OUT = ROOT / "crates" / "azoth-core" / "src" / "parameter_bounds_gen.rs"

#: The extension keys a bound is read from, and what each becomes.
BOUND_KEYS = ("min", "max", "min_inclusive", "max_inclusive")


class Refusal(Exception):
    """A parameter whose domain this file cannot state, so it stops rather than defaulting."""


def decisions() -> list[dict[str, Any]]:
    """Every numeric parameter, with its bound or its reason - in schema order."""
    schema = json.loads(SCHEMA.read_text(encoding="utf-8"))
    out: list[dict[str, Any]] = []
    for name, spec in schema["properties"].items():
        # `ion` is a class rather than a number with a unit, and `$defs.parameter` is the
        # shape a *number* has. Anything else is a property this file has no business in.
        if "$ref" not in spec:
            continue
        unit = spec.get("x-azoth-unit")
        if unit is None:
            raise Refusal(
                f"{name}: a parameter with no unit, so a bound could not be stated in one"
            )
        reason = spec.get("x-azoth-unbounded")
        stated = [key for key in BOUND_KEYS if f"x-azoth-{key}" in spec]
        if reason and stated:
            raise Refusal(f"{name}: both a bound and `x-azoth-unbounded`, which contradict")
        if reason:
            out.append({"name": name, "unit": unit, "reason": reason})
            continue
        if not stated:
            raise Refusal(
                f"{name}: neither a bound nor `x-azoth-unbounded`. Every parameter states which, "
                f"because a parameter nothing checks is the thing this generator exists to stop"
            )
        bound = {key: spec.get(f"x-azoth-{key}") for key in BOUND_KEYS}
        if bound["min"] is None and bound["max"] is None:
            raise Refusal(f"{name}: a bound with neither end, which would do nothing")
        if bound["min"] is None and "x-azoth-min-inclusive" in spec:
            raise Refusal(f"{name}: `x-azoth-min-inclusive` with no `x-azoth-min`")
        if bound["max"] is None and "x-azoth-max-inclusive" in spec:
            raise Refusal(f"{name}: `x-azoth-max-inclusive` with no `x-azoth-max`")
        if bound["min"] is not None and bound["max"] is not None and bound["min"] > bound["max"]:
            raise Refusal(f"{name}: its minimum is above its maximum")
        rationale = spec.get("x-azoth-rationale")
        if not rationale:
            raise Refusal(
                f"{name}: a bound with no `x-azoth-rationale`, so a reader cannot tell "
                f"what the domain is a statement about"
            )
        out.append({"name": name, "unit": unit, "bound": bound, "rationale": rationale})
    if not out:
        raise Refusal("no numeric parameter was read; this is looking in the wrong place")
    return out


def emit_python(rows: list[dict[str, Any]]) -> str:
    bounded = [row for row in rows if "bound" in row]
    unbounded = [row for row in rows if "reason" in row]
    lines = [
        '"""The domain each card parameter has to lie in, emitted from the component schema.',
        "",
        "Generated by `tools/gen_parameter_bounds.py` - see its module documentation for why",
        "a bound is stated in the parameter's own declared unit and why it is not a proof.",
        '**Do not edit by hand.**"""',
        "",
        "from typing import NamedTuple",
        "",
        "",
        "class ParameterBound(NamedTuple):",
        '    """One parameter\'s domain, in the unit `COMPONENT_PARAMETERS` declares it in."""',
        "",
        "    min: float | None",
        "    max: float | None",
        "    min_inclusive: bool",
        "    max_inclusive: bool",
        "    #: What the domain is a statement about - carried so both readers refuse in the same",
        "    #: words rather than each writing its own sentence.",
        "    rationale: str",
        "",
        "",
        "#: Each bounded parameter. A parameter absent here is in `UNBOUNDED`, with a reason.",
        "PARAMETER_BOUNDS: dict[str, ParameterBound] = {",
    ]
    for row in bounded:
        bound = row["bound"]
        lines.append(f'    "{row["name"]}": ParameterBound(')
        lines.append(f"        min={bound['min']!r},")
        lines.append(f"        max={bound['max']!r},")
        lines.append(f"        min_inclusive={bool(bound['min_inclusive'])!r},")
        lines.append(f"        max_inclusive={bool(bound['max_inclusive'])!r},")
        lines.append(f"        rationale={row['rationale']!r},")
        lines.append("    ),")
    lines += [
        "}",
        "",
        "#: Every parameter with no stated domain, and why - a decision, not an omission.",
        "UNBOUNDED: dict[str, str] = {",
    ]
    for row in unbounded:
        lines.append(f'    "{row["name"]}": {row["reason"]!r},')
    lines += ["}", ""]
    return "\n".join(lines)


def emit_rust(rows: list[dict[str, Any]]) -> str:
    lines = [
        "//! The domain each card parameter has to lie in, emitted from the component schema.",
        "//!",
        "//! Generated by `tools/gen_parameter_bounds.py` - see its module documentation for why a",
        "//! bound is stated in the parameter's own declared unit and why it is not a proof.",
        "//! **Do not edit by hand.**",
        "",
        "/// A bound on one parameter, in the unit `COMPONENT_PARAMETERS` declares it in.",
        "#[derive(Debug, Clone, Copy)]",
        "pub struct ParameterBound {",
        "    pub min: Option<f64>,",
        "    pub max: Option<f64>,",
        "    pub min_inclusive: bool,",
        "    pub max_inclusive: bool,",
        "    /// What the domain is a statement about, carried so both readers refuse in the same",
        "    /// words rather than each writing its own sentence.",
        "    pub rationale: &'static str,",
        "}",
        "",
        "/// Each bounded parameter, in the unit `COMPONENT_PARAMETERS` declares it in.",
        "pub const PARAMETER_BOUNDS: &[(&str, ParameterBound)] = &[",
    ]
    for row in rows:
        if "bound" not in row:
            continue
        bound = row["bound"]
        lines += [
            f'    ("{row["name"]}", ParameterBound {{',
            f"        min: {_rust_float(bound['min'])},",
            f"        max: {_rust_float(bound['max'])},",
            f"        min_inclusive: {str(bool(bound['min_inclusive'])).lower()},",
            f"        max_inclusive: {str(bool(bound['max_inclusive'])).lower()},",
            f"        rationale: {json.dumps(row['rationale'])},",
            "    }),",
        ]
    lines += [
        "];",
        "",
        "/// Every parameter with no stated domain, and why - a decision, not an omission.",
        "pub const UNBOUNDED: &[(&str, &str)] = &[",
    ]
    for row in rows:
        if "reason" not in row:
            continue
        lines.append(f'    ("{row["name"]}", "{row["reason"]}"),')
    lines += [
        "];",
        "",
        "/// The bound on `name`, or `None` where the domain is not stated.",
        "#[must_use]",
        "pub fn parameter_bound(name: &str) -> Option<&'static ParameterBound> {",
        "    PARAMETER_BOUNDS",
        "        .iter()",
        "        .find(|(parameter, _)| *parameter == name)",
        "        .map(|(_, bound)| bound)",
        "}",
        "",
    ]
    return "\n".join(lines)


def _rust_float(value: object) -> str:
    if value is None:
        return "None"
    # A bound is written as an f64 literal, and `0` is `0.0` in Rust rather than in Python.
    return f"Some({float(value)!r})"


def emit() -> str:
    rows = decisions()
    PY_OUT.write_text(emit_python(rows), encoding="utf-8")
    RUST_OUT.write_text(emit_rust(rows), encoding="utf-8")
    return f"{len(rows)} parameter(s), {sum(1 for r in rows if 'bound' in r)} bounded"


def main() -> None:
    if "--survey" in sys.argv:
        for row in decisions():
            if "bound" in row:
                print(f"{row['name']:28} bounded    {row['unit']}   {row['rationale']}")
            else:
                print(f"{row['name']:28} unbounded  {row['unit']}   {row['reason']}")
        return
    if "--check" in sys.argv:
        before = (PY_OUT.read_text(encoding="utf-8"), RUST_OUT.read_text(encoding="utf-8"))
        summary = emit()
        after = (PY_OUT.read_text(encoding="utf-8"), RUST_OUT.read_text(encoding="utf-8"))
        if before != after:
            raise SystemExit(f"gen_parameter_bounds: the tree is not what this emits ({summary})")
        print(f"gen_parameter_bounds: current ({summary})")
        return
    print(f"gen_parameter_bounds: {emit()}")


if __name__ == "__main__":
    main()
