#!/usr/bin/env python3
"""Check that the numbers and paths the normative pages state are the tree's own.

# Why this exists

`check_links.py` checks that a link resolves, and `prose_lint.py` checks how a
page is written - neither checks whether what a page *says* is true. `prose_lint`
excludes `docs/src` on purpose ("the docs are prose by definition"), so until now
nothing read the docs for facts, and a page quoting a count the code moved went
stale silently. It has: `docs/src/agentic/roadmap.md` said the library had 93
models against 117, `docs/src/calculus/index.md` said 41 `#print axioms` lines in
`Azoth/Axioms.lean` against 44, and `docs/src/calculus/rho.md` denied that the
interoperation layer existed a month after it was built.

The failure is not carelessness. It is that a catch-up pass edits the pages it
remembers and misses the ones it does not, and nothing tells it so.

# How it works, and why it is shaped this way

`docs/claims.toml` declares each claim as a span of a page with the measured
values replaced by `{holes}`. **The declaration holds no measured value.** That is
the property that stops the check being tuned until it passes: a failing run cannot
be silenced by editing `docs/claims.toml`, only by changing the page or the tree.
A hole that carries a value is a parse error.

Each hole names a *probe* - a function that measures the tree. Adding a claim is one
TOML stanza; adding a probe is a code change, deliberately, because a probe is the
only part that has to know how to measure something.

Two sweeps run without any declaration, because enumerating them by hand would
defeat the purpose: every inline-code token that begins with a repo directory must
exist. The prefix rule is what keeps the pages' NeqSim vocabulary - `PhaseGEUniquac`,
`thermo/util/leachman/`, `src/main` - out of it without an allowlist that would rot.

**The specs are swept by the same rule**, and by a second one over the two names a
spec cites most: a spec id and a capture file. A spec's `assumptions`, `description`
and `notes` fields are prose in a data file and name repo paths, neighbor ids and
captures exactly as a page does, so the path sweep reads them too. What it does not
reach is a spec's claim about *what the tree carries* - that is a fact about the
code, and it is held by `tools/check_unported.py`: a refused value is a
`[[unported]]` row, and a sentence making the same claim is an error there, because
a sentence is a claim nothing holds.

# What it does not do

- **Not NeqSim.** The upstream line counts are measured against a checkout CI does
  not have, and the page states the command that measures them. Declared as skips.
- **Generated pages and generated blocks.** Between `<!-- BEGIN GENERATED -->`
  markers is `docs-drift`'s to keep current, and a claim there fails for a reason the
  author cannot fix in place.
- **Rust module paths and third-party crates.** Measured: 16 of the 89 `a::b` tokens
  in the docs are honest non-paths (`middleware::command` is crate-relative,
  `uom::si` is not ours), so a sweep there would cry wolf. Only a declared claim
  resolves one.
- **Whether a proof is the *right* theorem.** `python/tests/test_lean_claims.py`
  owns that, and says so.
- **Prose.** No phrase list. `prose_lint.py` makes the argument for its own.

Usage:
    python tools/check_doc_claims.py            # check; 0 ok, 1 drift, 2 broken declaration
    python tools/check_doc_claims.py --list     # every claim and what it matched
"""

from __future__ import annotations

import argparse
import ast
import importlib
import importlib.util
import re
import sys
import tomllib
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent


def claims_path() -> Path:
    """The declaration, resolved at call time so a test can redirect `ROOT`.

    Outside `docs/src` so mdBook does not build it and `check_links.py` does not walk
    it; `docs-drift` covers its being committed.
    """
    return ROOT / "docs" / "claims.toml"


#: A token that begins a repository-relative path. The sweep uses this to tell an
#: azoth path from a NeqSim class or an upstream directory, which is why the
#: NeqSim vocabulary needs no allowlist.
REPO_DIRS = (
    "agents/",
    "crates/",
    "data/",
    "databank/",
    "docs/",
    "lean/",
    "python/",
    "skills/",
    "specs/",
    "tools/",
    "ui/",
    "validation/",
)

#: A path a page may name *because it does not exist* carries this on its line.
ESCAPE = "doc-claims-ok:"

#: Inline code spans. Fenced blocks are deliberately not read: measured, they carry
#: ten illustrative paths in `skills/` and one ASCII diagram.
INLINE_CODE = re.compile(r"`([^`\n]+)`")

#: A glob or template is a pattern, not a path to resolve.
IS_PATTERN = re.compile(r"[*?<>{}()\[\]]")

#: The marker a calculus page carries once, after its last status.
MARKER = re.compile(r"^\*Enforcement: (construction|check|nothing) — (.+)\*$", re.M)

#: A calculus page's status *block*, and the words it holds in bold. The block and not the
#: line: a page is free to wrap an italic span, and `rho.md` does - reading only the first line
#: left its `**proved**` and `**specified**` invisible to this sweep and to
#: `python/tests/test_claim_kinds.py`, which is a claim compared by nothing.
STATUS = re.compile(r"^\*Status: (.*?)\*\s*$", re.M | re.S)
BOLD = re.compile(r"\*\*([a-z]+)\*\*")


def page_files() -> list[Path]:
    """Every page this check reads: the book, plus the root pages `check_links` guards."""
    import check_links  # sibling tool, imported for its own page list

    pages = sorted(p for p in (ROOT / "docs" / "src").rglob("*.md") if p.name != "SUMMARY.md")
    pages += [ROOT / name for name in check_links.ROOT_PAGES if (ROOT / name).exists()]
    return pages


def spec_files() -> list[Path]:
    """Every spec, which carries inline code the same way a page does.

    A spec's `assumptions`, `description` and `notes` fields are prose in a data file, and
    they name repo paths, spec ids and capture files exactly as a page does. Nothing read
    them, so a spec citing a file that does not exist - or an id the registry never
    declared - was a statement with no reader. The two sweeps below give it one.
    """
    return sorted((ROOT / "specs").rglob("*.toml"))


# --- probes -----------------------------------------------------------------
#
# A *measuring* probe returns the number the page should state. Its holes are
# compared to it. A *validating* probe receives the text a hole captured and
# returns a failure reason, or None.


def positive(count: int, what: str) -> int:
    """Refuse a measurement of nothing.

    A regex or a path that stops matching after a refactor reports 0, and 0 compared
    against a page that also says 0 is a comfortable pass - the shape
    `check_json_keys.py` calls a check that did not run. Every count probe goes through
    here, so a moved tree fails loudly instead of silently agreeing.
    """
    if count <= 0:
        raise ProbeError(f"measured no {what}; the tree is not where this probe looks")
    return count


def specs_calcs() -> int:
    """The calculations. Enumerated the way `gen_registry.py` enumerates them."""
    return positive(len(sorted((ROOT / "specs" / "calcs").rglob("*.toml"))), "calculations")


def specs_models() -> int:
    """The models. Enumerated the way `gen_models.py` enumerates them."""
    return positive(len(sorted((ROOT / "specs" / "models").rglob("*.toml"))), "models")


def specs_ids() -> int:
    """Both registries' ids together, refusing a name that appears in both."""
    calcs = {p.stem for p in (ROOT / "specs" / "calcs").rglob("*.toml")}
    models = {p.stem for p in (ROOT / "specs" / "models").rglob("*.toml")}
    if calcs & models:
        raise ProbeError(f"these ids are in both trees: {sorted(calcs & models)}")
    return positive(len(calcs) + len(models), "ids")


def _palette_ids() -> set[str]:
    text = (ROOT / "crates" / "azoth-process" / "src" / "palette_gen.rs").read_text(
        encoding="utf-8"
    )
    ids = set(re.findall(r'include_str!\("\.\./\.\./\.\./specs/unit_ops/([^"]+)\.toml"\)', text))
    if not ids:
        raise ProbeError("no palette entries parsed from palette_gen.rs")
    # The include path carries the family directory (`column/absorption_column`); the
    # unit's own name is the last segment.
    return {i.rsplit("/", 1)[-1] for i in ids}


def unit_ops_declared() -> int:
    """The palette entries declared under `specs/unit_ops/`."""
    declared = len(sorted((ROOT / "specs" / "unit_ops").rglob("*.toml")))
    embedded = len(_palette_ids())
    if declared != embedded:
        raise ProbeError(f"{declared} specs but {embedded} embedded in palette_gen.rs")
    return declared


def unit_ops_kernels() -> int:
    """The entries that carry a kernel, which is not the same set as the dispatchable ones."""
    modules = (ROOT / "crates" / "azoth-process" / "src" / "kernels" / "mod.rs").read_text(
        encoding="utf-8"
    )
    declared = {m for m in re.findall(r"pub mod (\w+);", modules)}
    if not declared:
        raise ProbeError("no kernel modules parsed from kernels/mod.rs")
    return len(_palette_ids() & declared)


#: The two files an `executor::dispatch` table can be in. `UNRUNNABLE` is hand-written and stays
#: where it is; `DISPATCH` is generated, and read here from whichever file carries it rather than
#: from a path that moves when a table does.
_DISPATCH_SOURCES = (
    ROOT / "crates" / "azoth-process" / "src" / "executor" / "dispatch.rs",
    ROOT / "crates" / "azoth-process" / "src" / "executor" / "dispatch_gen.rs",
)


def _table(name: str) -> list[str]:
    """The ids in one `executor::dispatch` table, parsed from its own source."""
    text = "".join(path.read_text(encoding="utf-8") for path in _DISPATCH_SOURCES)
    # `)];` as well as `];`: a one-entry table is rustfmt'd across lines as
    # `&[(\n    "id",\n    "why",\n)];`.
    match = re.search(rf"pub const {name}[^=]*= &\[(.*?)\n\)?\];", text, re.S)
    if match is None:
        raise ProbeError(f"the {name} table is not where this probe looks for it")
    ids = re.findall(r'"(unit_ops\.\w+)"', match.group(1))
    if not ids:
        raise ProbeError(f"the {name} table parsed to no ids")
    return ids


def unit_ops_dispatch() -> int:
    """The entries the executor runs, cross-checked against the palette it runs them from."""
    dispatched, refused = _table("DISPATCH"), _table("UNRUNNABLE")
    if set(dispatched) & set(refused):
        raise ProbeError(f"both dispatched and refused: {sorted(set(dispatched) & set(refused))}")
    if len(dispatched) + len(refused) != unit_ops_declared():
        raise ProbeError("the two tables together are not the palette")
    return len(dispatched)


def unit_ops_refusals() -> int:
    """The entries refused by name - with a kernel and without one alike."""
    return len(_table("UNRUNNABLE"))


def specs_unported_rows() -> int:
    """The refusal rows the model specs declare, counted where the checker counts them.

    Loaded from `tools/check_unported.py` rather than re-derived here. The claim is
    *about* the declaration, and a second count of the same rows would be a second
    thing to keep in step - which is the defect the declaration exists to remove.
    """
    spec = importlib.util.spec_from_file_location(
        "check_unported", ROOT / "tools" / "check_unported.py"
    )
    if spec is None or spec.loader is None:
        raise ProbeError("tools/check_unported.py cannot be loaded")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    rows, failures = module.declared()
    if failures:
        raise ProbeError(f"the declaration is not well formed: {failures[0]}")
    return sum(len(keys) for keys in rows.values())


def vocabulary_units() -> int:
    text = (ROOT / "specs" / "vocabulary" / "vocabulary.toml").read_text(encoding="utf-8")
    return positive(text.count("[[units]]"), "units in vocabulary.toml")


def _gate_lines(name: str) -> list[str]:
    path = ROOT / "lean" / "Azoth" / name
    if not path.exists():
        raise ProbeError(f"lean/Azoth/{name} does not exist")
    return re.findall(r"^#print axioms (\S+)", path.read_text(encoding="utf-8"), re.M)


def lean_axioms_gated() -> int:
    lines = _gate_lines("Axioms.lean")
    if not lines:
        raise ProbeError("no `#print axioms` lines in Axioms.lean")
    return len(lines)


def lean_axioms_dim() -> int:
    return sum(1 for name in _gate_lines("Axioms.lean") if name.startswith("Azoth.Dim."))


def lean_gate_printed() -> int:
    return len(_gate_lines("Gate.lean"))


def lean_guards_gated() -> int:
    return len(_gate_lines("GuardGate.lean"))


def lean_modules() -> int:
    """Every module under `lean/Azoth/`, the gate files included."""
    modules = sorted((ROOT / "lean" / "Azoth").glob("*.lean"))
    if not modules:
        raise ProbeError("no .lean modules under lean/Azoth")
    return len(modules)


def _wrapper_counts() -> tuple[int, int]:
    """How many registered ids carry a generated pyo3 wrapper, and how many are hand-written.

    Read the way `python/tests/test_gen_python_wrappers.py` reads them, and for the same reason:
    the generated count is the generator's own output and the hand-written one is what the
    namespace modules still carry, so a calculation moved from one to the other moves both.
    """
    import sys as _sys

    _sys.path.insert(0, str(ROOT / "tools"))
    import rust_index

    registered = {function for _, _, function in rust_index.implementations()}
    generated = len(
        re.findall(
            r"^pub fn \w+\(",
            (ROOT / "crates" / "azoth-python" / "src" / "wrappers_gen.rs").read_text(
                encoding="utf-8"
            ),
            re.M,
        )
    )
    hand = 0
    # **A new namespace belongs on this list too.** It is not the generator's, so nothing refuses
    # when it is missing: the hand-written count just comes up short, and a claim whose two numbers
    # are supposed to partition the registry quietly stops partitioning it. `characterization` was
    # absent through three of its ids.
    for module in (
        "hydraulics",
        "eos",
        "thermal",
        "reactions",
        "standards",
        "process",
        "characterization",
    ):
        source = (ROOT / "crates" / "azoth-python" / "src" / f"{module}.rs").read_text(
            encoding="utf-8"
        )
        hand += len(
            [
                name
                for name in re.findall(r"^pub fn (\w+)\(\s*py: Python<'_>", source, re.M)
                if name in registered
            ]
        )
    return generated, hand


def _bridge_counts() -> tuple[int, int]:
    """`(generated, hand-written)` adapters, the same way the wrappers' pair is measured.

    The generated half is the emitted file's `__all__`, which is what the generator wrote; the
    hand-written half is the registered functions still written out in `_rust_bridge.py`. A name
    in both is what `test_gen_python_bridge.py` refuses, so the two counts are disjoint.

    **The registry is read from the Rust, not imported.** `spec-validate` runs this with neither
    `azoth` installed nor `python/src` on the path, so importing the generated tables made the probe
    a `ModuleNotFoundError` there while passing in every checkout - the wrappers' pair above reads
    `rust_index` for the same reason, and one reader is what keeps the two honest.
    """
    import sys as _sys

    _sys.path.insert(0, str(ROOT / "tools"))
    import rust_index

    text = (ROOT / "python" / "src" / "azoth" / "_rust_bridge_gen.py").read_text(encoding="utf-8")
    block = text.partition("__all__ = [")[2].partition("]")[0]
    generated = len(re.findall(r'^    "\w+",', block, re.M))
    if generated == 0:
        raise ProbeError("_rust_bridge_gen.py declares no adapters")

    registered = {function for _, _, function in rust_index.implementations()}
    source = (ROOT / "python" / "src" / "azoth" / "_rust_bridge.py").read_text(encoding="utf-8")
    hand = len([name for name in re.findall(r"^def (\w+)\(", source, re.M) if name in registered])
    return generated, hand


def _result_counts() -> tuple[int, int]:
    """`(generated, hand-written)` result dataclasses, the third of the same three pairs.

    The generated half is the `CALC_ID`s `core/result_gen.py` declares, which is what the
    generator wrote - one class per registered id, in `FIELDS` order. The hand-written half is any
    result class still declared in `core/result.py` itself, which is a re-export shim: the number
    is zero, and the claim states it so that a class drifting back into a module moves a count
    rather than going unnoticed.
    """
    generated = len(
        re.findall(
            r'^    CALC_ID: ClassVar\[str\] = "',
            (ROOT / "python" / "src" / "azoth" / "core" / "result_gen.py").read_text(
                encoding="utf-8"
            ),
            re.M,
        )
    )
    if generated == 0:
        raise ProbeError("core/result_gen.py declares no results")
    hand = len(
        re.findall(
            r"^class \w+\(_HasWarnings\):",
            (ROOT / "python" / "src" / "azoth" / "core" / "result.py").read_text(encoding="utf-8"),
            re.M,
        )
    )
    return generated, hand


def layout_results_generated() -> int:
    return _result_counts()[0]


def layout_results_hand() -> int:
    return _result_counts()[1]


def _generated_lengths(path: Path, names: tuple[str, ...]) -> dict[str, int]:
    """The length of each named module-level dict in a generated file.

    **Read, not imported.** `spec-validate` runs these probes with neither `azoth` installed nor
    `python/src` on the path, so a probe that imports the generated Python is a
    `ModuleNotFoundError` there while passing in every checkout - which is what `_bridge_counts`
    below already records, and what this function is the third instance of.
    """
    tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    out: dict[str, int] = {}
    for node in tree.body:
        if (
            isinstance(node, ast.AnnAssign)
            and isinstance(node.target, ast.Name)
            and node.target.id in names
            and isinstance(node.value, ast.Dict)
        ):
            out[node.target.id] = len(node.value.keys)
    missing = sorted(set(names) - set(out))
    if missing:
        raise ProbeError(f"{path.name} declares no {missing}")
    return out


def _bound_counts() -> tuple[int, int]:
    """`(bounded, unbounded)` card parameters, from the table both readers take.

    `docs/src/calculus/values.md` states how many parameters carry a bound and how many state
    why they have none. The two are the generated table's own length, so a parameter added to
    the schema without a decision moves a count rather than leaving the page behind.
    """
    lengths = _generated_lengths(
        ROOT / "python" / "src" / "azoth" / "core" / "_bounds_gen.py",
        ("PARAMETER_BOUNDS", "UNBOUNDED"),
    )
    bounded, unbounded = lengths["PARAMETER_BOUNDS"], lengths["UNBOUNDED"]
    if not bounded:
        raise ProbeError("no card parameter carries a bound, so the page states nothing")
    return bounded, unbounded


def bounds_card_bounded() -> int:
    return _bound_counts()[0]


def bounds_card_unbounded() -> int:
    return _bound_counts()[1]


def bounds_card_total() -> int:
    return _bound_counts()[0] + _bound_counts()[1]


def layout_bridge_generated() -> int:
    return _bridge_counts()[0]


def layout_bridge_hand() -> int:
    return _bridge_counts()[1]


def layout_wrappers_generated() -> int:
    return _wrapper_counts()[0]


def layout_wrappers_hand() -> int:
    return _wrapper_counts()[1]


def pairs_eos_files() -> int:
    eos = ROOT / "validation" / "eos"
    if not eos.is_dir():
        raise ProbeError("validation/eos does not exist")
    return positive(len(sorted(eos.glob("*.json"))), "paired cases")


def _skills() -> list[dict]:
    with (ROOT / "skills.toml").open("rb") as handle:
        return tomllib.load(handle)["skill"]


def skills_total() -> int:
    return positive(len(_skills()), "skills")


def skills_azoth() -> int:
    return sum(1 for s in _skills() if s.get("calculation_basis") == "azoth")


def skills_screening_p11() -> int:
    return sum(
        1
        for s in _skills()
        if s.get("calculation_basis") == "screening" and s.get("tranche") == "P11"
    )


def skills_placeholders() -> int:
    """Every skill that is not `azoth`-basis - a placeholder, by the catalogue's own words."""
    return positive(len(_skills()) - skills_azoth(), "placeholder skills")


def skills_screening() -> int:
    return sum(1 for s in _skills() if s.get("calculation_basis") == "screening")


def skills_advisory() -> int:
    return sum(1 for s in _skills() if s.get("calculation_basis") == "advisory")


def skills_data_retrieval() -> int:
    return sum(1 for s in _skills() if s.get("calculation_basis") == "data-retrieval")


def skills_hybrid() -> int:
    """Skills that call the library for part of their answer and screen for the rest."""
    return sum(1 for s in _skills() if s.get("calculation_basis") == "hybrid")


def _provenance_table() -> dict[str, dict[str, object]]:
    """The generated provenance table, read as data.

    Parsed rather than imported, for the reason the skills catalog is: this tool runs in
    `spec-validate` before the library is built, so reaching for `azoth` would end that. The
    generated file is a literal, so there is nothing importing it would add.
    """
    path = ROOT / "python" / "src" / "azoth" / "_provenance_gen.py"
    for node in ast.parse(path.read_text(encoding="utf-8")).body:
        if isinstance(node, ast.AnnAssign) and getattr(node.target, "id", "") == "PROVENANCE":
            table = ast.literal_eval(node.value)
            if not isinstance(table, dict):
                raise ProbeError("PROVENANCE is not a mapping")
            return table
    raise ProbeError("the generated table declares no PROVENANCE")


def provenance_ids() -> int:
    return positive(len(_provenance_table()), "provenance entries")


def provenance_partially_verified() -> int:
    """Every id with no independent oracle recorded for it.

    The normal case rather than a defect, and the number a page quoting it has to be held to:
    it moves whenever a validation case is added.
    """
    return sum(
        1
        for entry in _provenance_table().values()
        if entry.get("verification") == "partially_verified"
    )


class ProbeError(Exception):
    """A probe could not measure - the tree is not where the probe looks."""


MEASURES = {
    "specs.calcs": specs_calcs,
    "specs.models": specs_models,
    "specs.ids": specs_ids,
    "specs.unported_rows": specs_unported_rows,
    "unit_ops.declared": unit_ops_declared,
    "unit_ops.kernels": unit_ops_kernels,
    "unit_ops.dispatch": unit_ops_dispatch,
    "unit_ops.refusals": unit_ops_refusals,
    "vocabulary.units": vocabulary_units,
    "lean.axioms_gated": lean_axioms_gated,
    "lean.axioms_dim": lean_axioms_dim,
    "lean.gate_printed": lean_gate_printed,
    "lean.guards_gated": lean_guards_gated,
    "lean.modules": lean_modules,
    "layout.bridge_generated": layout_bridge_generated,
    "layout.bridge_hand": layout_bridge_hand,
    "bounds.card_bounded": bounds_card_bounded,
    "bounds.card_total": bounds_card_total,
    "bounds.card_unbounded": bounds_card_unbounded,
    "layout.results_generated": layout_results_generated,
    "layout.results_hand": layout_results_hand,
    "layout.wrappers_generated": layout_wrappers_generated,
    "layout.wrappers_hand": layout_wrappers_hand,
    "pairs.eos_files": pairs_eos_files,
    "skills.total": skills_total,
    "skills.azoth": skills_azoth,
    "skills.placeholders": skills_placeholders,
    "skills.screening": skills_screening,
    "skills.advisory": skills_advisory,
    "skills.data_retrieval": skills_data_retrieval,
    "skills.hybrid": skills_hybrid,
    "skills.screening_p11": skills_screening_p11,
    "provenance.ids": provenance_ids,
    "provenance.partially_verified": provenance_partially_verified,
}

#: Probes whose measurement is a number. A capture that is not one cannot be compared.
INT_PROBES = frozenset(MEASURES)


def validates_path(captured: str) -> str | None:
    """A `path` hole must resolve to a file or a directory."""
    target = captured.strip().split()[0] if captured.strip() else ""
    if not target:
        return "the hole captured nothing"
    resolved = ROOT / target.rstrip("/")
    if not resolved.exists():
        return f"{target!r} does not exist"
    return None


def validates_lean_gated(captured: str) -> str | None:
    """A `lean.gated` hole must name a theorem one of the two gate files prints axioms for."""
    name = captured.strip().strip("`")
    gated = _gate_lines("Axioms.lean") + _gate_lines("Gate.lean") + _gate_lines("GuardGate.lean")
    if name in gated:
        return None
    leaf = name.rsplit(".", 1)[-1]
    hits = [g for g in gated if g.rsplit(".", 1)[-1] == leaf]
    if len(hits) == 1:
        return None
    if not hits:
        return f"{name!r} is not gated in either gate file"
    return f"{name!r} is ambiguous: {hits}"


VALIDATES = {
    "path": validates_path,
    "lean.gated": validates_lean_gated,
}

#: English tens and units, so a page that writes "Twenty-seven" is measured too.
_WORDS = {
    "one": 1,
    "two": 2,
    "three": 3,
    "four": 4,
    "five": 5,
    "six": 6,
    "seven": 7,
    "eight": 8,
    "nine": 9,
    "ten": 10,
    "eleven": 11,
    "twelve": 12,
    "thirteen": 13,
    "fourteen": 14,
    "fifteen": 15,
    "sixteen": 16,
    "seventeen": 17,
    "eighteen": 18,
    "nineteen": 19,
    "twenty": 20,
    "thirty": 30,
    "forty": 40,
    "fifty": 50,
    "sixty": 60,
    "seventy": 70,
    "eighty": 80,
    "ninety": 90,
}


def as_int(captured: str) -> int:
    """Read a captured hole as a number: digits, or an English number word."""
    text = captured.strip().strip("*`_ ").replace(",", "").replace(" ", "-").lower()
    if text.lstrip("-").isdigit():
        return int(text)
    parts = text.split("-")
    if all(p in _WORDS for p in parts):
        return sum(_WORDS[p] for p in parts)
    raise ValueError(f"{captured!r} is not a number")


def escape_literal(literal: str) -> str:
    """Escape a template's literal text, with whitespace the only liberty.

    A span that wraps in the source may wrap in the page, so a run of whitespace matches
    any run. Nothing else is normalised: a page that retypes an em dash is a stale
    declaration rather than a claim that quietly moved.
    """
    out = []
    for part in re.split(r"(\s+)", literal):
        if not part:
            continue
        out.append(r"\s+" if part.isspace() else re.escape(part))
    return "".join(out)


def template_to_regex(template: str) -> re.Pattern[str]:
    """A claim's text, with each `{hole}` a named capture.

    A hole captures **one token**, not a span: every value these claims state is a number,
    a name or a path, and a hole allowed to contain spaces would let the leftmost match
    start in the prose before the sentence - `29 unit operations` reaching back to swallow
    a list marker and a clause. A claim needing a multi-word capture wants a new probe, not
    a looser hole.
    """
    out, index = [], 0
    for match in re.finditer(r"\{(\w+)\}", template):
        out.append(escape_literal(template[index : match.start()]))
        out.append(rf"(?P<{match.group(1)}>\S+?)")
        index = match.end()
    out.append(escape_literal(template[index:]))
    return re.compile("".join(out))


# --- the claims file --------------------------------------------------------


class Claim:
    def __init__(self, raw: dict, index: int) -> None:
        self.index = index
        self.page = str(raw.get("page", ""))
        self.text = str(raw.get("text", ""))
        self.measure = dict(raw.get("measure", {}))
        self.only = int(raw.get("only", 1))

    def check(self) -> list[str]:
        """Return the failures for this claim - empty when it agrees with the tree."""
        failures: list[str] = []
        where = f"docs/claims.toml claim {self.index + 1} (page = {self.page})"
        page = ROOT / self.page
        if not page.exists():
            return [f"{where}: {self.page} does not exist"]
        if not self.text.strip():
            return [f"{where}: the claim has no `text`"]
        if not self.measure:
            return [f"{where}: the claim names no probe"]

        source = page.read_text(encoding="utf-8")
        pattern = template_to_regex(self.text)
        matches = list(pattern.finditer(source))
        if len(matches) != self.only:
            return [
                f"{where}: the claim's text matches its page {len(matches)} time(s), "
                f"not {self.only} - the page was reworded and the declaration is stale\n"
                f"         declared: {self.text.strip()}\n"
                f"         page:     {self.page}"
            ]
        found = matches[0]

        for hole, probe_name in self.measure.items():
            if hole not in found.groupdict():
                failures.append(f"{where}: the template has no hole {{{hole}}}")
                continue
            captured = found.group(hole)
            try:
                if probe_name in INT_PROBES:
                    stated = as_int(captured)
                    measured = MEASURES[probe_name]()
                    if stated != measured:
                        failures.append(
                            f"{where}: the page states {hole} = {stated} and the tree is "
                            f"{measured}\n"
                            f"         claim:    {self.text.strip()}\n"
                            f"         measured: {probe_name} = {measured}\n"
                            f"         measure:  python tools/check_doc_claims.py --probe "
                            f"{probe_name}"
                        )
                elif probe_name in VALIDATES:
                    reason = VALIDATES[probe_name](captured)
                    if reason is not None:
                        failures.append(
                            f"{where}: {reason}\n"
                            f"         claim:    {self.text.strip()}\n"
                            f"         captured: {captured!r} (probe {probe_name})"
                        )
                else:
                    failures.append(
                        f"{where}: {probe_name!r} is not a known probe. Known: "
                        f"{sorted(MEASURES) + sorted(VALIDATES)}"
                    )
            except ProbeError as error:
                failures.append(f"{where}: the probe {probe_name} could not measure: {error}")
        return failures


def load_claims() -> list[Claim]:
    declaration = claims_path()
    if not declaration.exists():
        raise ProbeError(f"{declaration.relative_to(ROOT)} does not exist")
    with declaration.open("rb") as handle:
        try:
            data = tomllib.load(handle)
        except tomllib.TOMLDecodeError as error:
            # A broken declaration is exit 2, not a traceback: the message has to say
            # which file and which line, because that is what the author fixes.
            raise ProbeError(f"{declaration.relative_to(ROOT)} does not parse: {error}") from None
    raw = data.get("claim", [])
    if not raw:
        raise ProbeError(f"{declaration.relative_to(ROOT)} declares no [[claim]]")
    claims = [Claim(entry, i) for i, entry in enumerate(raw)]

    # The property that stops the check being tuned until it passes.
    for claim in claims:
        for hole, probe in claim.measure.items():
            if probe in INT_PROBES and hole.strip().isdigit():
                raise ProbeError(
                    f"{declaration.relative_to(ROOT)} claim {claim.index + 1}: the hole "
                    f"{{{hole}}} carries a value. The declaration holds no measured value - "
                    f"that is what makes a failure unfixable by editing this file."
                )
    return claims


def load_skips() -> list[dict]:
    declaration = claims_path()
    if not declaration.exists():
        return []
    with declaration.open("rb") as handle:
        return list(tomllib.load(handle).get("skip", []))


# --- the path sweep ---------------------------------------------------------


def sweep_paths(pages: list[Path]) -> tuple[int, list[str], list[str]]:
    """Every inline-code repo path must exist. Returns (checked, failures, escapes).

    The escape is line-scoped but *token-exempting*: a line that names a path because it
    is absent carries the marker, and the paths on that line which do exist are simply
    ordinary paths. So an escape with nothing missing on its line is the finding - not an
    existing path sharing a line with an absent one, which is the usual case and was the
    heuristic that cried wolf on its first run.
    """
    failures: list[str] = []
    escapes: list[str] = []
    checked = 0
    for page in pages:
        for lineno, line in enumerate(page.read_text(encoding="utf-8").splitlines(), start=1):
            exempt = ESCAPE in line
            missing: list[str] = []
            for token in INLINE_CODE.findall(line):
                word = token.strip().split()[0] if token.strip() else ""
                if not word.startswith(REPO_DIRS) or IS_PATTERN.search(word):
                    continue
                checked += 1
                target = word.rstrip(".,;:")
                if not (ROOT / target.rstrip("/")).exists():
                    missing.append(target)
            where = f"{page.relative_to(ROOT)}:{lineno}"
            if exempt:
                reason = line.split(ESCAPE, 1)[1].strip()
                reason = reason.split("-->", 1)[0].strip()
                if not reason:
                    failures.append(f"{where}: the `{ESCAPE}` escape carries no reason")
                elif not missing:
                    failures.append(
                        f"{where}: the `{ESCAPE}` escape exempts nothing - every path on the "
                        f"line exists, so the marker is stale"
                    )
                else:
                    escapes += [f"{where} {target!r}" for target in missing]
            else:
                failures += [
                    f"{where}: inline code names {target!r}, which does not exist"
                    for target in missing
                ]
    return checked, failures, escapes


# --- what a spec names ------------------------------------------------------

#: What an id looks like: two lowercase snake_case segments and nothing else. Measured, all
#: 223 the tree declares are exactly that - never three, never a capital - which is what
#: separates an id from a NeqSim package path such as
#: `process.equipment.reactor.KineticReaction` that a spec cites in the same backticks.
SPEC_ID = re.compile(r"^[a-z][a-z_0-9]*\.[a-z][a-z_0-9]*$")


def spec_ids() -> tuple[set[str], set[str]]:
    """Every id the tree declares, and the namespaces those ids begin with.

    The first segment of an id is a namespace rather than a guess: the set is read out of
    the ids themselves, so a token starting with one is a reference to a spec and a token
    starting with anything else is not. That is what makes the sweep need no allowlist -
    `algorithm.tolerance`, `henry.rs` and `test_model_cases.py` are excluded by *shape*,
    and the namespaces are derived rather than listed.
    """
    declared: set[str] = set()
    for path in spec_files():
        try:
            document = tomllib.loads(path.read_text(encoding="utf-8"))
        except tomllib.TOMLDecodeError:
            continue
        if isinstance(document.get("id"), str):
            declared.add(document["id"])
    return declared, {name.split(".")[0] for name in declared}


def spec_strings(value: object):
    """Every string in a parsed spec, at any depth."""
    if isinstance(value, str):
        yield value
    elif isinstance(value, dict):
        for item in value.values():
            yield from spec_strings(item)
    elif isinstance(value, list):
        for item in value:
            yield from spec_strings(item)


def sweep_spec_refs() -> tuple[int, list[str], list[str]]:
    """Every spec id and capture a spec names in inline code must exist.

    A spec cites its neighbours and its captures by name in prose fields nothing read. The
    two token shapes here are unambiguous by *shape*: a dotted token whose first segment is
    a namespace the tree declares, and a token ending `.tsv`. A capture is named either as a
    bare file or as the path under `validation/neqsim/captures/`, and both are read as the
    file they name.

    Kept apart from `sweep_paths` because the corpus is a different one - a spec is read as
    a document rather than line by line, since its prose is one long TOML string.
    """
    declared, namespaces = spec_ids()
    captures = ROOT / "validation" / "neqsim" / "captures"
    failures: list[str] = []
    checked = 0
    for path in spec_files():
        try:
            document = tomllib.loads(path.read_text(encoding="utf-8"))
        except tomllib.TOMLDecodeError:
            continue
        where = path.relative_to(ROOT)
        for text in spec_strings(document):
            for token in INLINE_CODE.findall(text):
                word = token.strip()
                if word.endswith(".tsv"):
                    checked += 1
                    if not (captures / word.rsplit("/", 1)[-1]).exists():
                        failures.append(
                            f"{where}: names the capture {word!r}, which is not under "
                            f"{captures.relative_to(ROOT)}"
                        )
                    continue
                if not SPEC_ID.match(word) or word.split(".")[0] not in namespaces:
                    continue
                checked += 1
                if word not in declared:
                    failures.append(
                        f"{where}: names the id {word!r}, which the tree does not declare"
                    )
    return checked, failures, []


# --- the port's sources -----------------------------------------------------

#: The three spec trees a port has to account for itself in.
SPEC_TREES = ("unit_ops", "models", "calcs")

#: Where a spec says its arithmetic came from. Empty is a refusal rather than a default: the
#: port rule is that a port names the paper, the NeqSim class and its own notes, and a spec
#: whose `[source].standard` is a bare "NeqSim" has named none of them.
SOURCE = "source.standard"

#: How a spec says its source cannot be named - a class deleted upstream, a source that was
#: never published. The marker travels with the file, so the reason is where the next reader
#: looks, and the count of markers is printed rather than being a silence.
SOURCE_MARKER = "source-ok:"

#: What counts as naming a source: a `.java` file, a backticked symbol or path, a work with a
#: year, a standard's number, or an author list. Five shapes because a port's source is one of
#: five kinds of place - a NeqSim class, a NeqSim method, a paper, a standard, or a textbook the
#: spec cites by its authors.
NAMES_A_CLASS = re.compile(
    r"\.java\b"
    r"|`[A-Za-z_][A-Za-z0-9_.()/]*`"
    r"|\b\d{4}\b"
    r"|\b[A-Z]{2,}[ -]?[A-Z]?\d+"
    # An author list: two or more capitalised surnames joined by `&` or `,` - "Smith, Van Ness &
    # Abbott". Measured, the flash family cites its textbook this way, and those specs name the
    # NeqSim class beside it in `edition`.
    r"|\b[A-Z][a-z]+(?:,? ?& ?|, )[A-Z][a-z]+"
)


def sweep_sources() -> tuple[int, list[str], list[str]]:
    """Every spec says where its arithmetic came from, and says it in a checkable shape.

    Returns `(checked, failures, markers)`.

    **A bare \"NeqSim\" is the failure this exists for.** The port rules are that a port names
    the paper, the NeqSim class and its own notes, and that a divergence from the source is a
    finding rather than a repair - so the class matters: `NeqSim's getAntoineVaporPressure` is
    checkable against the pin, where `NeqSim` alone is a claim about a whole library. Measured,
    127 of the 221 specs name a NeqSim source and 94 a paper or a standard, and every one of
    them names a class, a method or a year.

    A spec that cannot name its source - a class deleted upstream, a correlation that was never
    published - carries a `source-ok:` marker with the reason, which is the same idiom the
    numerics lint uses for a literal a source writes.
    """
    failures: list[str] = []
    markers: list[str] = []
    checked = 0
    found = 0
    for tree in SPEC_TREES:
        for path in sorted((ROOT / "specs" / tree).rglob("*.toml")):
            where = str(path.relative_to(ROOT))
            try:
                spec = tomllib.loads(path.read_text(encoding="utf-8"))
            except tomllib.TOMLDecodeError as error:
                failures.append(f"{where}: is not valid TOML: {error}")
                continue
            source = spec.get("source")
            checked += 1
            if not isinstance(source, dict):
                failures.append(
                    f"{where}: has no `[source]`. A port names the paper or the NeqSim class it "
                    f"came from, and a spec with neither is one nobody can check"
                )
                continue
            standard = str(source.get("standard", "")).strip()
            # **The whole `[source]` block is the citation.** `standard` is the work and
            # `edition` is the locus within it - a line number, a method, the paragraphs of a
            # textbook - and several specs name the NeqSim class there rather than here.
            citation = f"{standard} {source.get('edition', '')}"
            found += 1
            text = path.read_text(encoding="utf-8")
            if SOURCE_MARKER in text:
                markers.append(f"{where} - {standard or 'no standard'}")
                continue
            if not standard:
                failures.append(
                    f"{where}: `{SOURCE}` is empty. Refused rather than defaulted: a spec that "
                    f"names no source cannot be held to one"
                )
            elif not NAMES_A_CLASS.search(citation):
                failures.append(
                    f"{where}: its `[source]` is {standard!r}, and neither it nor `edition` "
                    f"names a class, method, paper, standard or author - so nothing says what "
                    f"this spec was ported from. Name the `.java`, the method, the paper with "
                    f"its year, or mark it `{SOURCE_MARKER} <reason>`"
                )
    if not found:
        raise ProbeError("measured no spec with a `[source]`; the tree is not where this looks")
    return checked, failures, markers


# --- the enforcement sweep ---------------------------------------------------


def swept_pages() -> list[Path]:
    """Every page that states where its claims are enforced.

    Two trees. `docs/src/calculus/` is normative for the *types* - twelve layers, each with
    a Lean module and a status. `docs/src/architecture/` is normative for the *surface* a
    front-end binds, and a claim about the surface is the same kind of claim: `dirty`
    means this, a failed run takes the values with it, the doors disagree about nothing.
    Until this sweep reached it, those rested on tests no page named.
    """
    pages: list[Path] = []
    for tree in ("calculus", "architecture"):
        pages += sorted((ROOT / "docs" / "src" / tree).glob("*.md"))
    return pages


def sweep_enforcement() -> tuple[int, list[str], list[str]]:
    """Every page of the two normative trees states where its claims are enforced.

    Returns (checked, failures, nothing-claims).

    What no page said until now is whether the *implementation* half of a claim is
    enforced by a construction, by a check, or by nothing - and that difference is the
    whole review: a claim that is proved and enforced by nothing is one whose
    implementation half nobody would notice breaking.

    So each page carries one marker, after its last status:

        *Enforcement: construction — `crates/azoth-core/src/unit_vocab_gen.rs`*

    Three kinds, and the two that name something must name a path that resolves:

    * `construction` - the types or the generated tables make the error impossible.
    * `check` - code or a test fails when it happens. Production code that refuses at run
      time and a test are one kind here; the path says which it is.
    * `nothing` - nothing enforces it, **and that is legal only where no status claims a
      proof**. A page with a proved claim and no enforcement is the defect this sweep
      exists for; a page with no status line at all may not claim `nothing` either,
      because nothing then says the absence is deliberate rather than unexamined.
    """
    failures: list[str] = []
    nothing_claims: list[str] = []
    checked = 0
    pages = swept_pages()
    if not pages:
        raise ProbeError("no pages under the swept trees; this is looking in the wrong place")
    for page in pages:
        where = str(page.relative_to(ROOT))
        text = page.read_text(encoding="utf-8")
        markers = MARKER.findall(text)
        if len(markers) != 1:
            failures.append(
                f"{where}: {len(markers)} `*Enforcement:*` marker(s); every page states its "
                f"enforcement exactly once"
            )
            continue
        checked += 1
        kind, reason = markers[0]
        statuses = STATUS.findall(text)
        # The word each block opens with: a status block states its status first and argues
        # after, so the later bold words are the prose's emphasis rather than statuses.
        status_words = {bold[0] for block in statuses if (bold := BOLD.findall(block))}
        if "*Status:" in text and text.index("*Enforcement:") < text.rindex("*Status:"):
            failures.append(
                f"{where}: the marker sits above the last `*Status:*`; it summarises the page's "
                f"claims, so it goes after them"
            )
        # A marker that names a status this page never states is one written about another
        # page's claim - the shape a contradiction takes when two pages describe one claim.
        stray = sorted({word for word in BOLD.findall(reason)} - status_words)
        if status_words and stray:
            failures.append(
                f"{where}: the marker names {stray}, and this page's statuses are "
                f"{sorted(status_words)} - so the enforcement stated here is about a claim "
                f"this page does not make"
            )
        if not statuses and kind == "nothing":
            failures.append(
                f"{where}: claims `nothing` and states no `*Status:*`, so nothing says the "
                f"absence is deliberate rather than unexamined"
            )
        if kind == "nothing":
            proved = [word for line in statuses for word in BOLD.findall(line) if word == "proved"]
            if proved:
                failures.append(
                    f"{where}: claims `nothing` while a status says **proved**; a proved claim "
                    f"whose implementation half enforces nothing is the finding this sweep is for"
                )
            nothing_claims.append(f"{where} - {reason}")
            continue
        named = [word.strip("`") for word in INLINE_CODE.findall(reason)]
        resolved = [word for word in named if (ROOT / word.rstrip("/")).exists()]
        if not named:
            failures.append(f"{where}: `{kind}` names no path, so the claim is unreadable")
        elif not resolved:
            failures.append(
                f"{where}: `{kind}` names {named}, and none of them exists - so the enforcement "
                f"the page claims is not where the page says it is"
            )
    return checked, failures, nothing_claims


def check_skips(skips: list[dict]) -> list[str]:
    """A skip must still quote text its page carries, so it goes stale rather than silent."""
    failures: list[str] = []
    for skip in skips:
        page = ROOT / str(skip.get("page", ""))
        text = str(skip.get("text", ""))
        reason = str(skip.get("reason", "")).strip()
        if not page.exists():
            failures.append(f"docs/claims.toml skip: {skip.get('page')} does not exist")
            continue
        if not reason:
            failures.append(f"docs/claims.toml skip for {skip.get('page')}: no reason")
        if text and text not in page.read_text(encoding="utf-8"):
            failures.append(
                f"docs/claims.toml skip for {skip.get('page')}: the skipped text is no longer on "
                f"the page, so the skip is stale: {text!r}"
            )
    return failures


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--list", action="store_true", help="print every claim and its match")
    parser.add_argument("--probe", metavar="NAME", help="print one probe's measurement and exit")
    args = parser.parse_args(argv)

    if args.probe:
        if args.probe not in MEASURES:
            sys.exit(f"check_doc_claims: no probe {args.probe!r}. Known: {sorted(MEASURES)}")
        try:
            print(MEASURES[args.probe]())
        except ProbeError as error:
            sys.exit(f"check_doc_claims: {args.probe} could not measure: {error}")
        return 0

    pages = page_files()
    if not pages:
        sys.exit("check_doc_claims: no pages found to check")

    try:
        claims = load_claims()
        skips = load_skips()
    except ProbeError as error:
        sys.exit(f"check_doc_claims: {error}")

    failures: list[str] = []
    if args.list:
        for claim in claims:
            pattern = template_to_regex(claim.text)
            hit = pattern.search((ROOT / claim.page).read_text(encoding="utf-8"))
            print(f"  {claim.page}: {'matched' if hit else 'NOT MATCHED'}  {claim.text.strip()!r}")
    failures += [f for claim in claims for f in claim.check()]
    checked, sweep_failures, escapes = sweep_paths(pages + spec_files())
    failures += sweep_failures
    refs, ref_failures, _ = sweep_spec_refs()
    failures += ref_failures
    pages_swept, layer_failures, unenforced = sweep_enforcement()
    failures += layer_failures
    specs, spec_failures, source_markers = sweep_sources()
    failures += spec_failures
    failures += check_skips(skips)

    if failures:
        for failure in failures:
            print(f"  ERROR  {failure}", file=sys.stderr)
        print(
            f"\ncheck_doc_claims: FAILED with {len(failures)} problem(s), "
            f"{checked} path(s) and {refs} spec reference(s) checked",
            file=sys.stderr,
        )
        return 1

    print(
        f"check_doc_claims: OK ({len(claims)} claim(s) across "
        f"{len({c.page for c in claims})} page(s), {checked} path(s) and "
        f"{refs} spec reference(s) checked, "
        f"{pages_swept} swept page(s), "
        f"{specs} spec source(s))"
    )
    if source_markers:
        print(f"  {len(source_markers)} source(s) marked unnameable:")
        for entry in source_markers:
            print(f"    {entry}")
    if unenforced:
        # Printed rather than passed over: a claim nothing enforces is a state to look at,
        # and a number a reader can watch is the difference between a judgement and a habit.
        print(f"  nothing enforces {len(unenforced)}:")
        for entry in unenforced:
            print(f"    {entry}")
    if escapes:
        print(f"  escapes: {len(escapes)}")
        for escape in escapes:
            print(f"    {escape}")
    if skips:
        print(f"  skipped: {len(skips)}")
        for skip in skips:
            print(f"    {skip.get('page')}: {str(skip.get('reason'))[:70]}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
