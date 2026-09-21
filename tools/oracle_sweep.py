#!/usr/bin/env python3
"""Re-run every probe and say what the oracle moved, and which test literal is now stale.

# Why this exists

`validation/neqsim/*.java` drives NeqSim and their stdout is committed to
`validation/neqsim/captures/*.tsv`. The crate tests hard-code those numbers. When the pin
moves, a capture that changed and a test that was not renumbered are both invisible: the
tests' bars are wide enough that a move of `1e-9` passes them, which is how two whole test
files stayed on the previous revision through a refresh.

**The committed capture is the previous pin.** That is what makes this decidable with one
jar: a test literal is stale exactly when it equals what a probe printed *before* - the
value still in the committed capture - and not what it prints now. So the run is: rebuild
every capture, take each key that moved, and look for a Rust literal holding the old value.

# What it runs

Every probe with no arguments, plus every invocation written as a recipe in
`validation/neqsim/captures/README.md` - which is where a state that needs arguments is
documented, and which is the only registry of them that is not this file.

# What it cannot see

A state a probe prints only when invoked for it and which no recipe documents. Those are
reported as probes with no committed capture, which is the honest form of the gap: the
number exists, and nothing here can say whether it moved.

Usage:
    python tools/oracle_sweep.py [--jar PATH] [--timeout SECONDS]
"""

from __future__ import annotations

import argparse
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import gen_neqsim_cases as gen
import manifest as manifest_tool

ROOT = Path(__file__).resolve().parent.parent
PROBES = ROOT / "validation" / "neqsim"
README = gen.CAPTURES / "README.md"


def pinned_jar() -> Path:
    """The jar built from the revision `databank/manifest.toml` pins."""
    pinned, _ = manifest_tool.read()
    commit = next((u.commit for u in pinned.upstreams if u.id == "neqsim"), "")
    return PROBES / f"neqsim-{commit[:7]}.jar"


#: A `javac`/`java` recipe in the README, which is how an argument-taking state is written
#: down. The command is taken as far as the redirect into the capture it produces.
RECIPE = re.compile(
    r"java\s+-cp\s+[^\s:]*:neqsim-[0-9a-f]+\.jar\s+(?P<run>[^\s>]+(?:\s+[-\w.]+)*)\s*>\s*"
    r"captures/(?P<capture>[a-z0-9_]+\.tsv)"
)

#: A number in Rust source, at the precision a transcription from a probe carries. Nine
#: significant digits is below every oracle value in this tree and above the constants a
#: test legitimately writes for itself.
LITERAL = re.compile(r"-?\d+\.\d+(?:[eE][-+]?\d+)?")
SIGNIFICANT = 10


def java_files() -> list[Path]:
    return sorted(PROBES.glob("*.java"))


def package_of(source: Path) -> str | None:
    found = re.search(r"^package\s+([\w.]+);", source.read_text(encoding="utf-8"), re.M)
    return found.group(1) if found else None


def build(classes: Path, jar: Path) -> bool:
    """Compile every probe into `classes`, which is also where they are run from."""
    done = subprocess.run(
        [
            "javac",
            "-proc:none",
            "-d",
            str(classes),
            "-cp",
            f"{classes}:{jar}",
            *map(str, java_files()),
        ],
        capture_output=True,
        text=True,
        check=False,
    )
    if done.returncode != 0:
        print(done.stderr.strip()[:2000], file=sys.stderr)
    return done.returncode == 0


def recipes() -> list[tuple[str, str]]:
    """Every documented invocation, as `(arguments, capture name)`."""
    return [
        (m.group("run"), m.group("capture"))
        for m in RECIPE.finditer(README.read_text(encoding="utf-8"))
    ]


def probe_names() -> dict[str, str]:
    """Every probe's runnable name, keyed by the class file's stem."""
    out: dict[str, str] = {}
    for source in java_files():
        pkg = package_of(source)
        out[source.stem] = f"{pkg}.{source.stem}" if pkg else source.stem
    return out


def run(jar: Path, classes: Path, name: str, arguments: str, timeout: float) -> str | None:
    try:
        done = subprocess.run(
            ["java", "-cp", f"{classes}:{jar}", name, *arguments.split()],
            capture_output=True,
            text=True,
            timeout=timeout,
            check=False,
        )
    except subprocess.TimeoutExpired:
        return None
    return done.stdout if done.returncode == 0 else None


def body(text: str) -> str:
    """A probe's output without its own header, which is prose and moves on its own."""
    return "\n".join(
        line for line in text.splitlines() if line.strip() and not line.lstrip().startswith("#")
    )


def as_float(value: str) -> float | None:
    try:
        return float(value)
    except ValueError:
        return None


def relative(mine: float, theirs: float) -> float:
    """The convention `neqsim_layer_diff` uses: absolute when the oracle is zero."""
    if theirs == 0.0:
        return abs(mine)
    return abs(mine / theirs - 1.0)


def moved_keys(before: dict[str, str], after: dict[str, str]) -> list[tuple[str, float, float]]:
    """Every key present in both whose value the probe no longer reproduces."""
    changed = []
    for key, old in before.items():
        was, is_now = as_float(old), as_float(after.get(key, ""))
        if was is not None and is_now is not None and relative(is_now, was) > 1e-12:
            changed.append((key, was, is_now))
    return changed


def is_stale(value: float, key: tuple[str, float, float]) -> bool:
    """Whether a number in a test is the value the probe printed before, and not the one
    it prints now - which is the whole definition, since the capture *is* the before."""
    _, was, is_now = key
    return relative(value, was) <= 1e-12 and relative(value, is_now) > 1e-12


def rust_literals() -> list[tuple[str, int, str, float]]:
    out = []
    for path in sorted(ROOT.glob("crates/*/tests/*.rs")):
        for number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), start=1):
            for token in LITERAL.finditer(line):
                text = token.group()
                if len(re.sub(r"[-.]|^0+|0+$", "", text)) >= SIGNIFICANT:
                    out.append((str(path.relative_to(ROOT)), number, text, float(text)))
    return out


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--jar", type=Path, default=None)
    parser.add_argument("--timeout", type=float, default=120.0, help="per probe, in seconds")
    args = parser.parse_args()

    jar = args.jar or pinned_jar()
    if not jar.exists():
        raise SystemExit(f"no jar at {jar}: build it from the checkout, as .gitignore says")

    workdir = Path(tempfile.mkdtemp(prefix="oracle-sweep-"))
    try:
        if not build(workdir, jar):
            return 2
        names = probe_names()

        # Every probe once with no arguments, then every documented invocation, which is
        # the only way a state that needs arguments is reached.
        runs: list[tuple[str, str, str, str]] = [
            (stem, name, "", "") for stem, name in sorted(names.items())
        ]
        for arguments, produced in recipes():
            head, *rest = arguments.split()
            stem = head.split(".")[-1]
            if stem in names:
                runs.append((stem, names[stem], " ".join(rest), produced))

        # `fresh` is keyed by the probe's stem for a bare run and by the capture's name for
        # a recipe, because a recipe is only ever compared with the capture it produces.
        fresh: dict[str, str] = {}
        for stem, name, arguments, produced in runs:
            text = run(jar, workdir, name, arguments, args.timeout)
            if text is None:
                print(f"# {stem} {arguments}: did not run", file=sys.stderr)
                continue
            fresh[produced or stem] = text

        matched: dict[str, str] = {}
        for stem, text in fresh.items():
            if (gen.CAPTURES / f"{stem}.tsv").exists():
                matched[stem] = f"{stem}.tsv"
                continue
            for candidate in sorted(gen.CAPTURES.glob("*.tsv")):
                if body(text) == body(candidate.read_text(encoding="utf-8")):
                    matched[stem] = candidate.name
                    break

        stale: list[tuple[str, int, str, str, float, float]] = []
        moved: list[tuple[str, list[tuple[str, float, float]]]] = []
        for stem, capture_name in sorted(matched.items()):
            committed = gen.CAPTURES / capture_name
            changed = moved_keys(
                gen.pairs(committed.read_text(encoding="utf-8")), gen.pairs(fresh[stem])
            )
            if not changed:
                continue
            moved.append((stem, changed))
            for path, number, text, value in rust_literals():
                for key, was, is_now in changed:
                    if is_stale(value, (key, was, is_now)):
                        stale.append((path, number, text, f"{capture_name}:{key}", was, is_now))

        print(f"# {len(fresh)} run(s), {len(matched)} compared with a committed capture")
        if moved:
            print(f"\n## {len(moved)} capture(s) moved\n")
            for stem, changed in moved:
                print(f"  {stem} -> {matched[stem]}")
                for key, was, is_now in changed[:8]:
                    print(f"    {key}: {was:.15g} -> {is_now:.15g}")
                if len(changed) > 8:
                    print(f"    ... and {len(changed) - 8} more")
                print()
        else:
            print("\n## no capture moved - every probe reproduces what is committed\n")

        if stale:
            print(f"## {len(stale)} stale literal(s) - a test still holds the old value\n")
            seen: set[tuple[str, int]] = set()
            for path, number, text, where, was, is_now in stale:
                if (path, number) in seen:
                    continue
                seen.add((path, number))
                print(f"  {path}:{number}  {text}  ({where})")
                print(f"    was {was:.15g}, is {is_now:.15g}")
            print()

        loose = sorted({stem for stem, _, _, produced in runs if not produced} - set(matched))
        if loose:
            print(f"## {len(loose)} probe(s) with no committed capture\n")
            print("  " + ", ".join(loose))
            print()

        # The mirror: a capture no run reproduces is one nothing can re-measure, and it is
        # either a state that lost its recipe or one that upstream no longer runs at all.
        compared = set(matched.values())
        orphans = sorted(p.name for p in gen.CAPTURES.glob("*.tsv") if p.name not in compared)
        if orphans:
            print(f"## {len(orphans)} capture(s) no run reproduces\n")
            print("  " + ", ".join(orphans))
            print()
            print()

        return 1 if stale else 0
    finally:
        shutil.rmtree(workdir, ignore_errors=True)


if __name__ == "__main__":
    sys.exit(main())
