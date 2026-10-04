#!/usr/bin/env python3
"""The guard-discharge obligations, as the generator and the checker both read them.

`docs/src/calculus/numerics.md` requires every partial function to be total, clamped or
deliberately `NaN`, and the code states which at the call site. Eleven of those statements are
*guarded*: a check above the arithmetic makes the domain condition hold. Nothing read the
sentence, so this is the reader for it.

**The hypotheses are derived and the conclusion is declared.** A guard's hypotheses are its
spec's `severity = "error"` rows on the inputs it names, turned into the predicate that *passes* -
which is `gen_registry.emit_range_check` read backwards, since the Rust emits the predicate that
fails. The *conclusion* is the composite expression the kernel evaluates, and only a person can
say which bounds discharge which expression, so it is declared in `lean/guards.toml` beside the
proof.

One reader rather than one per tool, the reason `vocabulary_table.py` exists: the generator and
the checker must agree about what an obligation is, and two readings of one manifest would drift.

**`equals` is an exclusion.** `RangeCheck::violated` checks it first and returns unconditionally
("an exclusion rather than an interval"), so a row `equals = 1, severity = "error"` refuses
`x = 1` and the hypothesis it yields is `x ≠ 1`. Read the other way it generates false theorems.
"""

from __future__ import annotations

import pathlib
import re
import tomllib
from dataclasses import dataclass

ROOT = pathlib.Path(__file__).resolve().parent.parent
MANIFEST = ROOT / "lean" / "guards.toml"
CRATES = ROOT / "crates"
SPECS = ROOT / "specs"

#: What may own an edge-case class or a guard site. A closed vocabulary, because an owner outside
#: it would be a word nobody can check.
OWNERS = ("lean", "check", "construction", "port", "unformalised")

#: The token a call site carries, and how it is told apart from English. `guarded:` also occurs in
#: prose about a *state* rather than a partial function ("Reproduced rather than guarded: it is
#: what NeqSim returns"), so a marker line is one whose payload is a **single token** - and that
#: token must be a dotted key, or the sweep refuses it rather than passing it by.
MARKER = re.compile(r"guarded:\s*(\S+)\s*$")
KEY = re.compile(r"[a-z0-9_]+(?:\.[a-z0-9_]+)+")

#: Words Lean reserves. A parameter or hypothesis named one of these would shadow rather than
#: declare, so it is refused rather than mangled.
KEYWORDS = frozenset(
    [
        "at",
        "by",
        "def",
        "do",
        "else",
        "end",
        "example",
        "export",
        "for",
        "from",
        "fun",
        "have",
        "if",
        "in",
        "import",
        "include",
        "let",
        "match",
        "mutual",
        "namespace",
        "noncomputable",
        "opaque",
        "open",
        "partial",
        "private",
        "protected",
        "public",
        "rec",
        "return",
        "section",
        "show",
        "structure",
        "syntax",
        "theorem",
        "then",
        "type",
        "unsafe",
        "variable",
        "where",
        "with",
    ]
)

_IDENT = re.compile(r"^[A-Za-z_][A-Za-z0-9_']*$")


class Refusal(Exception):
    """The manifest or the tree says something this reader will not guess at."""


@dataclass(frozen=True)
class Guard:
    """One guarded site: the key its marker carries, and the claim declared about it."""

    key: str
    site: str
    owner: str
    spec: str | None = None
    conditions: tuple[str, ...] = ()
    proof: str | None = None
    source: str | None = None
    reason: str | None = None
    checker: str | None = None

    @property
    def name(self) -> str:
        """The Lean identifier: the key with its separator written as the one Lean allows."""
        return self.key.replace(".", "_")


def manifest() -> tuple[list[Guard], list[dict]]:
    """`(guards, classes)` from `lean/guards.toml`, refused where it is malformed."""
    if not MANIFEST.exists():
        raise Refusal(f"{MANIFEST.relative_to(ROOT)} does not exist")
    document = tomllib.loads(MANIFEST.read_text(encoding="utf-8"))

    guards: list[Guard] = []
    for row in document.get("guard", []):
        guard = Guard(
            key=row["key"],
            site=row["site"],
            owner=row["owner"],
            spec=row.get("spec"),
            conditions=tuple(row["conditions"] if "conditions" in row else [row["condition"]])
            if ("conditions" in row or "condition" in row)
            else (),
            proof=row.get("proof"),
            source=row.get("source"),
            reason=row.get("reason"),
            checker=row.get("checker"),
        )
        if guard.owner not in OWNERS:
            raise Refusal(f"{guard.key}: owner {guard.owner!r} is not one of {list(OWNERS)}")
        if guard.owner == "lean" and (not guard.conditions or not guard.proof):
            raise Refusal(f"{guard.key}: owner `lean` needs a condition and a proof")
        if guard.owner == "unformalised" and not guard.reason:
            raise Refusal(f"{guard.key}: owner `unformalised` needs a reason")
        if guard.owner == "port" and not guard.source:
            raise Refusal(f"{guard.key}: owner `port` needs the source it reproduces")
        if not (ROOT / guard.site).exists():
            raise Refusal(f"{guard.key}: site {guard.site!r} does not exist")
        if guard.checker and not (ROOT / guard.checker).exists():
            raise Refusal(f"{guard.key}: checker {guard.checker!r} does not exist")
        guards.append(guard)

    keys = [guard.key for guard in guards]
    if len(set(keys)) != len(keys):
        raise Refusal(f"a key appears twice in {MANIFEST.name}: {sorted(keys)}")

    classes = list(document.get("class", []))
    for entry in classes:
        if entry.get("owner") not in OWNERS:
            raise Refusal(
                f"class {entry.get('key')!r}: owner {entry.get('owner')!r} is not an owner"
            )
        if entry.get("owner") == "unformalised" and not entry.get("reason"):
            raise Refusal(f"class {entry.get('key')!r}: owner `unformalised` needs a reason")
        for field in ("checker", "gate"):
            named = entry.get(field)
            if named and not (ROOT / named).exists():
                raise Refusal(f"class {entry.get('key')!r}: {field} {named!r} does not exist")
    class_keys_seen = [entry.get("key") for entry in classes]
    if len(set(class_keys_seen)) != len(class_keys_seen):
        raise Refusal(f"a class appears twice in {MANIFEST.name}: {sorted(class_keys_seen)}")
    return guards, classes


def class_keys() -> set[str]:
    """The shapes a `valid_range` row in `specs/` actually takes.

    `<severity>-<equals|interval>-<outside|inside>`, plus `computed-from` for a row naming the
    inputs its quantity is derived from. Derived, so a new shape in a spec is a key with no class
    entry rather than a declaration nobody notices.
    """
    out: set[str] = set()
    for spec in specs().values():
        for row in spec.get("valid_range", []):
            kind = "equals" if "equals" in row else "interval"
            out.add(f"{row.get('severity')}-{kind}-{row.get('when', 'outside')}")
            if "computed_from" in row:
                out.add("computed-from")
    return out


def specs() -> dict[str, dict]:
    """Every calculation and model spec, by id."""
    out: dict[str, dict] = {}
    for namespace in ("calcs", "models"):
        for path in sorted((SPECS / namespace).rglob("*.toml")):
            document = tomllib.loads(path.read_text(encoding="utf-8"))
            out[document["id"]] = document
    return out


def _identifier(name: str, where: str) -> str:
    if not _IDENT.match(name) or name in KEYWORDS:
        raise Refusal(f"{where}: {name!r} is not an identifier this file can declare in Lean")
    return name


def hypotheses(spec: dict, where: str) -> list[tuple[str, str]]:
    """The spec's error rows as `(hypothesis name, Lean type)`, in declaration order.

    One per row, so two rows on one quantity give two hypotheses - `rachford_rice_binary`'s `K1`
    carries both `0 < K1` and `K1 ≠ 1`, and both are needed. A warning row is never one of these:
    assuming it would claim a region the code does not enforce.
    """
    out: list[tuple[str, str]] = []
    seen: dict[tuple[str, str], int] = {}
    for row in spec.get("valid_range", []):
        if row.get("severity") != "error":
            continue
        quantity = _identifier(row["quantity"], where)
        if row.get("when") == "inside":
            raise Refusal(
                f'{where}: {quantity} has an `error` row with `when = "inside"`, whose passing '
                "set is the complement of the band; this reader does not implement that"
            )
        kind: str
        if "equals" in row:
            kind, assertion = "ne", f"{quantity} ≠ {_number(row['equals'])}"
        elif "min" in row:
            inclusive = row.get("min_inclusive", True)
            kind = "ge" if inclusive else "gt"
            op = "≤" if inclusive else "<"
            assertion = f"{_number(row['min'])} {op} {quantity}"
        elif "max" in row:
            inclusive = row.get("max_inclusive", True)
            kind = "le" if inclusive else "lt"
            op = "≤" if inclusive else "<"
            assertion = f"{quantity} {op} {_number(row['max'])}"
        else:
            raise Refusal(f"{where}: {quantity} has an error row with no bound this reader knows")
        index = seen.get((quantity, kind), 0)
        seen[(quantity, kind)] = index + 1
        suffix = f"_{index}" if index else ""
        out.append((f"h_{quantity}_{kind}{suffix}", assertion))
    return out


def _number(value: object) -> str:
    """A declared bound as a Lean literal: an integer where it is one, a decimal otherwise."""
    number = float(value)  # type: ignore[arg-type]
    if number.is_integer():
        return str(int(number))
    return repr(number)


def parameters(spec: dict) -> list[str]:
    """The inputs an error row constrains, in declaration order.

    An input with no error row is not a parameter: it is not one of the checks the site performs,
    so it would be a variable the statement mentions and nothing bounds.
    """
    constrained = {
        row["quantity"] for row in spec.get("valid_range", []) if row.get("severity") == "error"
    }
    return [name for name in spec["inputs"] if name in constrained]


def markers() -> dict[str, tuple[str, int]]:
    """Every `guarded:` marker in the kernels, as `{key: (relative path, line)}`."""
    out: dict[str, tuple[str, int]] = {}
    for path in sorted(CRATES.glob("*/src/**/*.rs")):
        text = path.read_text(encoding="utf-8")
        for number, line in enumerate(text.split("\n"), start=1):
            if "guarded:" not in line:
                continue
            found = MARKER.search(line)
            if found is None:
                continue  # prose: the payload is more than one token
            key = found.group(1)
            if KEY.fullmatch(key) is None:
                raise Refusal(
                    f"{path.relative_to(ROOT)}:{number}: `guarded: {key}` is not a dotted key, "
                    "so nothing can bind it"
                )
            if key in out:
                raise Refusal(
                    f"{key!r} is marked twice: {out[key][0]}:{out[key][1]} and "
                    f"{path.relative_to(ROOT)}:{number}"
                )
            out[key] = (str(path.relative_to(ROOT)), number)
    return out


def theorem(guard: Guard, spec: dict) -> str:
    """The Lean theorem for one `owner = "lean"` guard, as source text."""
    where = guard.key
    rows = hypotheses(spec, where)
    if not rows:
        raise Refusal(
            f"{where}: {guard.spec} declares no error-severity bound, so there is nothing to "
            "discharge this guard - a guard with no check above it is not a guard"
        )
    declared = {
        row["quantity"] for row in spec.get("valid_range", []) if row.get("severity") == "error"
    }
    for condition in guard.conditions:
        for token in re.findall(r"[A-Za-z_][A-Za-z0-9_']*", condition):
            if token in KEYWORDS:
                raise Refusal(f"{where}: condition names the Lean keyword {token!r}")
            if token[0].isupper() or token in declared:
                continue
            if token in {"Real", "sqrt", "rpow"}:
                continue
            raise Refusal(
                f"{where}: condition names {token!r}, which is neither a declared input "
                f"({sorted(declared)}) nor a constant this reader knows"
            )

    names = [_identifier(name, where) for name in parameters(spec)]
    claim = " ∧ ".join(f"({condition})" for condition in guard.conditions)
    if len(guard.conditions) == 1:
        claim = guard.conditions[0]
    binder = " ".join(names)
    hypothesis_args = "\n    ".join(f"({name} : {assertion})" for name, assertion in rows)
    return (
        f"/-- `{guard.spec}`: the checks the site performs discharge this guard. -/\n"
        f"theorem {guard.name} ({binder} : ℝ)\n"
        f"    {hypothesis_args} :\n"
        f"    {claim} := by\n"
        f"  {guard.proof}\n"
    )


if __name__ == "__main__":
    guards, classes = manifest()
    print(f"guards: {len(guards)}  classes: {len(classes)}")
    for guard in guards:
        print(f"  {guard.owner:12s} {guard.key}")
