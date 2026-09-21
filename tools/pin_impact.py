#!/usr/bin/env python3
"""Say which oracle-backed tests a NeqSim pin move can have invalidated.

# Why this exists

azoth pins NeqSim in `databank/manifest.toml` and hard-codes the oracle's numbers into
`crates/*/tests/*.rs`. Nothing makes those literals move when the pin does: the crate's
own tests are the only thing that would notice, and their bars are wide enough that a pin
move of `1e-9` passes them. The last refresh found two whole test files still asserting the
*previous* revision's numbers.

The question that has to be asked at every pin move is **which of this library's tests can
this move have invalidated**, and it has three parts:

1. which NeqSim classes changed between the two revisions;
2. which of those this tree cites at all — a probe, a spec, a test, a kernel;
3. of those, which changed a *value* rather than a brace or a line break.

The third is the one that is not guessable. Upstream's history is mostly fixes, and a
fix that reaches a number azoth asserts is exactly the thing this exists to catch — but
so is a formatting sweep, which reaches nothing and must not be re-measured. So each
changed-and-cited file's numeric literals are compared as multisets across the two
revisions, which separates them without reading either diff.

# What it does not do

It does not run a probe and it does not touch the jar, so it cannot say whether a *value*
that changed reaches the specific number a test asserts. It says which files to look at,
and which of those can be dismissed in one line. Exits 0 either way: a move that changes
values is the expected case, not a failure.

Usage:
    python tools/pin_impact.py --from 805cf0f910819a19fdc45702fc44c4a3675d93d8 \\
                               --to   f0c7436c6923766b1e22957b7075f650600457a7

The checkout is a local prerequisite, like the jar: `~/.claude` is not where NeqSim lives.
Set it with `--checkout` or `$NEQSIM_CHECKOUT`; the default is the copy the refresh used.
"""

from __future__ import annotations

import argparse
import os
import re
import subprocess
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent

#: Where the NeqSim checkout lives. Gitignored, like the jar built from it.
DEFAULT_CHECKOUT = Path("/tmp/neqsim-check/neqsim")

#: The trees that are physics. The `process` tree is unit operations and
#: `fluidmechanics` is hydraulics azoth writes its own; neither can invalidate a test here.
TREES = (
    "src/main/java/neqsim/thermo",
    "src/main/java/neqsim/thermodynamicoperations",
    "src/main/java/neqsim/physicalproperties",
)

#: Where this tree cites a NeqSim class. The probes drive them, the specs and the tests
#: record their numbers, and the kernels transcribe their arithmetic.
CITING = (
    "validation/neqsim/*.java",
    "validation/eos/*.json",
    "specs/**/*.toml",
    "crates/*/src/**/*.rs",
    "crates/*/tests/*.rs",
    "python/src/**/*.py",
)

#: A number in Java source, including the forms a formatting pass leaves alone.
LITERAL = re.compile(r"-?\d+\.\d+(?:[eE][-+]?\d+)?|-?\d+[eE][-+]?\d+")


def git(checkout: Path, *args: str) -> str:
    """Run git in the checkout and return stdout, or raise with its stderr."""
    done = subprocess.run(
        ["git", "-C", str(checkout), *args], capture_output=True, text=True, check=False
    )
    if done.returncode != 0 and "does not exist" not in done.stderr:
        raise SystemExit(f"git {' '.join(args)}: {done.stderr.strip()}")
    return done.stdout


def citation_index() -> dict[str, list[str]]:
    """Every NeqSim class name this tree names, against the documents that name it.

    **A name this tree declares itself is not a citation of NeqSim's.** `Component`,
    `Fluid` and `Water` are azoth's own vocabulary and NeqSim's both, and matching on the
    bare word makes them cite every document in the tree; the declarations are collected
    first so they can be subtracted.
    """
    declared: set[str] = set()
    for pattern in ("crates/*/src/**/*.rs", "python/src/**/*.py"):
        for path in ROOT.glob(pattern):
            declared.update(
                re.findall(
                    r"\b(?:pub\s+)?(?:struct|enum|trait|type)\s+([A-Z][A-Za-z0-9]*)",
                    path.read_text(encoding="utf-8", errors="replace"),
                )
            )

    index: dict[str, list[str]] = {}
    for pattern in CITING:
        for path in sorted(ROOT.glob(pattern)):
            text = path.read_text(encoding="utf-8", errors="replace")
            for name in set(re.findall(r"\b([A-Z][A-Za-z0-9]{3,})\b", text)):
                if name not in declared:
                    index.setdefault(name, []).append(str(path.relative_to(ROOT)))
    return index


def literals_at(
    checkout: Path, rev: str, path: str, cache: dict[tuple[str, str], list[str] | None]
) -> list[str] | None:
    """The file's numeric literals at a revision, or None when the file is not there."""
    key = (rev, path)
    if key not in cache:
        text = git(checkout, "show", f"{rev}:{path}")
        cache[key] = None if not text else sorted(m.group() for m in LITERAL.finditer(text))
    return cache[key]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--from", dest="before", required=True, help="the revision pinned before")
    parser.add_argument("--to", dest="after", required=True, help="the revision pinned now")
    parser.add_argument(
        "--checkout",
        type=Path,
        default=Path(os.environ.get("NEQSIM_CHECKOUT", DEFAULT_CHECKOUT)),
    )
    parser.add_argument(
        "--all",
        action="store_true",
        help="also list the cited files with nothing to re-measure, and the uncited ones",
    )
    args = parser.parse_args()

    if not (args.checkout / ".git").exists():
        raise SystemExit(
            f"no NeqSim checkout at {args.checkout}: pass --checkout or set $NEQSIM_CHECKOUT"
        )

    changed = git(
        args.checkout, "diff", "--name-only", args.before, args.after, "--", *TREES
    ).split("\n")
    changed = [path for path in changed if path.endswith(".java")]
    index = citation_index()
    cache: dict[tuple[str, str], list[str] | None] = {}

    cited: list[tuple[str, str, list[str]]] = []
    uncited: list[tuple[str, str, list[str]]] = []
    for path in changed:
        name = Path(path).stem
        hits = sorted(set(index.get(name, [])))
        (cited if hits else uncited).append((path, name, hits))

    quiet: list[tuple[str, bool]] = []
    moved: list[tuple[str, list[str], list[str], list[str]]] = []
    for path, name, hits in cited:
        before = literals_at(args.checkout, args.before, path, cache) or []
        after = literals_at(args.checkout, args.after, path, cache) or []
        gone = sorted((Counter(before) - Counter(after)).elements())
        added = sorted((Counter(after) - Counter(before)).elements())
        if gone:
            moved.append((name, hits, gone, added))
        else:
            # The same literals, or additions only - new code rather than a number that
            # moved. Either way there is nothing here to re-measure.
            quiet.append((name, bool(added)))

    print(f"# {len(changed)} class file(s) changed between the two revisions")
    print(f"# {len(cited)} of them are cited here, {len(uncited)} are not")
    print(f"# {len(moved)} cited file(s) moved a literal; {len(quiet)} did not\n")

    if moved:
        print("## re-measure these - a literal the tree may assert has moved\n")
        for name, hits, gone, added in moved:
            print(f"  {name}")
            if added:
                print(f"    added : {', '.join(added[:10])}")
            print(f"    gone  : {', '.join(gone[:10])}")
            print(f"    cited by: {', '.join(hits[:6])}")
            print()

    if args.all:
        print("## changed, cited, and nothing to re-measure\n")
        for quiet_name, added_only in quiet:
            why = "new code only" if added_only else "a formatting change"
            print(f"  {quiet_name}  ({why})")
        print()
        print("## changed and cited by nothing here - ignore\n")
        for _path, name, _ in uncited:
            print(f"  {name}")
        print()

    return 0


if __name__ == "__main__":
    sys.exit(main())
