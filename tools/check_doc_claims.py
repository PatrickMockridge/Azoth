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


def page_files() -> list[Path]:
    """Every page this check reads: the book, plus the root pages `check_links` guards."""
    import check_links  # noqa: PLC0415 - sibling tool, imported for its own list

    pages = sorted(p for p in (ROOT / "docs" / "src").rglob("*.md") if p.name != "SUMMARY.md")
    pages += [ROOT / name for name in check_links.ROOT_PAGES if (ROOT / name).exists()]
    return pages


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


def _table(name: str) -> list[str]:
    """The ids in one `executor::dispatch` table, parsed from its own source."""
    text = (ROOT / "crates" / "azoth-process" / "src" / "executor" / "dispatch.rs").read_text(
        encoding="utf-8"
    )
    match = re.search(rf"pub const {name}[^=]*= &\[(.*?)\n\];", text, re.S)
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


def lean_modules() -> int:
    """Every module under `lean/Azoth/`, the gate files included."""
    modules = sorted((ROOT / "lean" / "Azoth").glob("*.lean"))
    if not modules:
        raise ProbeError("no .lean modules under lean/Azoth")
    return len(modules)


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


class ProbeError(Exception):
    """A probe could not measure - the tree is not where the probe looks."""


MEASURES = {
    "specs.calcs": specs_calcs,
    "specs.models": specs_models,
    "specs.ids": specs_ids,
    "unit_ops.declared": unit_ops_declared,
    "unit_ops.kernels": unit_ops_kernels,
    "unit_ops.dispatch": unit_ops_dispatch,
    "unit_ops.refusals": unit_ops_refusals,
    "vocabulary.units": vocabulary_units,
    "lean.axioms_gated": lean_axioms_gated,
    "lean.axioms_dim": lean_axioms_dim,
    "lean.gate_printed": lean_gate_printed,
    "lean.modules": lean_modules,
    "pairs.eos_files": pairs_eos_files,
    "skills.total": skills_total,
    "skills.azoth": skills_azoth,
    "skills.placeholders": skills_placeholders,
    "skills.screening": skills_screening,
    "skills.advisory": skills_advisory,
    "skills.data_retrieval": skills_data_retrieval,
    "skills.screening_p11": skills_screening_p11,
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
    gated = _gate_lines("Axioms.lean") + _gate_lines("Gate.lean")
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
    "one": 1, "two": 2, "three": 3, "four": 4, "five": 5, "six": 6, "seven": 7, "eight": 8,
    "nine": 9, "ten": 10, "eleven": 11, "twelve": 12, "thirteen": 13, "fourteen": 14,
    "fifteen": 15, "sixteen": 16, "seventeen": 17, "eighteen": 18, "nineteen": 19,
    "twenty": 20, "thirty": 30, "forty": 40, "fifty": 50, "sixty": 60, "seventy": 70,
    "eighty": 80, "ninety": 90,
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


def check_skips(skips: list[dict], pages: list[Path]) -> list[str]:
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
    checked, sweep_failures, escapes = sweep_paths(pages)
    failures += sweep_failures
    failures += check_skips(skips, pages)

    if failures:
        for failure in failures:
            print(f"  ERROR  {failure}", file=sys.stderr)
        print(
            f"\ncheck_doc_claims: FAILED with {len(failures)} problem(s), "
            f"{checked} path(s) checked",
            file=sys.stderr,
        )
        return 1

    print(
        f"check_doc_claims: OK ({len(claims)} claim(s) across "
        f"{len({c.page for c in claims})} page(s), {checked} path(s) checked)"
    )
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
