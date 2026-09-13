#!/usr/bin/env python3
"""Emit provenance.json: what went into a release, and what it hashes to.

For each calculation the record names the spec it came from, the code that
implements it in both languages, and the tests that exercise it, each with a
SHA-256 of the file's committed bytes. Alongside those it records the git commit
and tag, the lock files, and - for a release - the built artifacts.

# Determinism

The output is byte-for-byte reproducible for a given commit. There is no
timestamp: the git commit and tag already identify the release precisely, and a
wall-clock field would mean two people generating provenance for the same commit
got different files. That matters because the file is signed, and a signature
over something non-reproducible is harder to check than it needs to be.

# Why raw bytes, not canonicalised content

Hashing the file as committed avoids inventing a canonical form for YAML and
JSON. A canonicalisation scheme would be a second definition of the content, and
the two would eventually disagree - which is the failure mode this whole file
exists to rule out.

# What it does not cover

Per-calc hashes cover that calc's own spec, code and tests. They are not a
transitive closure: `azoth-core` and the shared solver affect every result, so
those are listed once under `shared` rather than repeated per calc. Read the two
together.

Usage:
    python tools/provenance.py                          # source provenance
    python tools/provenance.py --tag v0.1.0 --artifact dist/*.whl
    python tools/provenance.py --verify provenance.json # check a record
"""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
from pathlib import Path
from typing import Any

try:
    import yaml
except ImportError:  # pragma: no cover
    sys.exit("provenance requires PyYAML")

ROOT = Path(__file__).resolve().parent.parent
SPEC_DIR = ROOT / "specs" / "calcs"

SCHEMA_VERSION = 1

#: Files that affect every calculation's result. Hashed once here rather than
#: once per calc, because repeating them would suggest they are per-calc.
SHARED = (
    "crates/azoth-core/src/range.rs",
    "crates/azoth-core/src/result.rs",
    "crates/azoth-core/src/units.rs",
    "crates/azoth-core/src/warning.rs",
    "crates/azoth-hydraulics/src/solver.rs",
    "crates/azoth-hydraulics/src/fittings.rs",
    "crates/azoth-hydraulics/src/fluids.rs",
    "crates/azoth-hydraulics/src/provenance.rs",
    "crates/azoth-hydraulics/src/spec_gen.rs",
    "crates/azoth-python/src/hydraulics.rs",
    "crates/azoth-python/src/results.rs",
    "python/src/azoth/core/range.py",
    "python/src/azoth/core/units.py",
    "python/src/azoth/_registry_gen.py",
    "python/src/azoth/_rust_bridge.py",
    # The tools are part of the chain: gen_registry.py writes spec_gen.rs, which
    # every calc reads its range checks from, so a change there changes results.
    # provenance.py hashes itself - the hash of the output then depends on the
    # generator, which is the property you want and not a circularity.
    "tools/gen_registry.py",
    "tools/gen_docs.py",
    "tools/spec_lint.py",
    "tools/provenance.py",
    "tools/check_links.py",
)

#: Data files whose contents change results.
DATA = (
    "data/fittings/crane_k_factors.csv",
    "data/fluids/water.csv",
    "data/fluids/air.csv",
)

LOCK_FILES = ("Cargo.lock", "uv.lock")

#: The licence texts, hashed like anything else. A release's legal terms are part
#: of the release: if a licence file changed between this record and the tree you
#: have, you are looking at different terms, and a hash is the only way a reader
#: finds that out without diffing by eye.
LICENCES = ("LICENSE", "LICENSE-CC-BY-4.0")


def sha256_of(path: Path) -> str:
    """SHA-256 of a file's bytes."""
    return hashlib.sha256(path.read_bytes()).hexdigest()


def describe(path: str) -> dict[str, Any]:
    """A path and its hash, or an explicit record that it is missing.

    Missing files are recorded rather than skipped. A calc with no test file is a
    fact a reader of this document needs, and quietly omitting the entry would
    make the absence invisible.
    """
    full = ROOT / path
    if not full.is_file():
        return {"path": path, "sha256": None, "present": False}
    return {"path": path, "sha256": sha256_of(full), "present": True}


def git(*args: str) -> str | None:
    """Run a git command, returning None when git or the repo is unavailable.

    Provenance generated from an unpacked source directory is still useful - the
    file hashes stand on their own - so this degrades rather than failing.
    """
    try:
        result = subprocess.run(
            ["git", *args],
            cwd=ROOT,
            capture_output=True,
            text=True,
            check=False,
        )
    except FileNotFoundError:
        return None
    if result.returncode != 0:
        return None
    return result.stdout.strip() or None


def calc_entry(spec: dict[str, Any]) -> dict[str, Any]:
    """Provenance for one calculation."""
    name = spec["id"].split(".")[-1]
    namespace = spec["id"].split(".")[0]

    tests = [
        describe(f"python/tests/{namespace}/test_{name}.py"),
        describe(f"crates/{'azoth-' + namespace}/tests/{name}.rs"),
    ]
    return {
        "id": spec["id"],
        "name": spec["name"],
        "equation": spec["equation"],
        "verification": spec["verification"]["status"],
        "spec": describe(f"specs/calcs/{namespace}/{name}.yaml"),
        "code": [
            describe(f"python/src/azoth/{namespace}/reference/{name}.py"),
            describe(f"crates/azoth-{namespace}/src/{name}.rs"),
        ],
        "tests": tests,
    }


def build(artifacts: list[str], tag: str | None) -> dict[str, Any]:
    """Assemble the provenance record."""
    paths = sorted(SPEC_DIR.rglob("*.yaml"))
    calcs = [yaml.safe_load(p.read_text(encoding="utf-8")) for p in paths]
    calcs.sort(key=lambda c: c["id"])

    status = git("status", "--porcelain")
    resolved_tag = tag if tag is not None else git("describe", "--tags", "--exact-match")

    return {
        "schema_version": SCHEMA_VERSION,
        "library": {"name": "azoth", "version": version()},
        "git": {
            "commit": git("rev-parse", "HEAD"),
            "tag": resolved_tag,
            # Recorded because provenance for a dirty tree describes something
            # nobody else can reproduce, and the reader deserves to know.
            "dirty": bool(status),
        },
        "calcs": [calc_entry(c) for c in calcs],
        "shared": [describe(p) for p in SHARED],
        "data": [describe(p) for p in DATA],
        "lock_files": [describe(p) for p in LOCK_FILES],
        "licences": [describe(p) for p in LICENCES],
        "artifacts": [describe(a) for a in artifacts],
        "signature": {
            "sigstore_bundle": None,
            "note": (
                "Filled in by .github/workflows/release.yml after keyless signing "
                "with cosign. Null for provenance generated from a source tree, "
                "which is unsigned by construction. See TRUST.md."
            ),
        },
    }


def version() -> str:
    """The library version, read from pyproject so it cannot drift."""
    text = (ROOT / "pyproject.toml").read_text(encoding="utf-8")
    for line in text.splitlines():
        if line.startswith("version = "):
            return line.split("=", 1)[1].strip().strip('"')
    return "unknown"


def render(record: dict[str, Any]) -> str:
    """Serialise deterministically: sorted keys, fixed indent, trailing newline."""
    return json.dumps(record, indent=2, sort_keys=True) + "\n"


def verify(path: Path) -> int:
    """Recompute the hashes in a record and report any that no longer match.

    This is what makes the instructions in TRUST.md executable rather than
    aspirational: a reader can check a published record against the tree they
    have, rather than being told to trust it.
    """
    recorded = json.loads(path.read_text(encoding="utf-8"))
    problems: list[str] = []
    checked = 0

    def compare(entry: dict[str, Any], where: str) -> None:
        nonlocal checked
        if entry.get("sha256") is None:
            return
        checked += 1
        actual = describe(entry["path"])
        if not actual["present"]:
            problems.append(f"{where}: {entry['path']} is missing")
        elif actual["sha256"] != entry["sha256"]:
            problems.append(
                f"{where}: {entry['path']} hashes to {actual['sha256'][:16]}..., "
                f"recorded {entry['sha256'][:16]}..."
            )

    for entry in recorded.get("calcs", []):
        for key in ("spec",):
            compare(entry[key], f"{entry['id']}.{key}")
        for item in entry["code"]:
            compare(item, f"{entry['id']}.code")
        for item in entry["tests"]:
            compare(item, f"{entry['id']}.tests")

    for section in ("shared", "data", "lock_files", "licences"):
        for item in recorded.get(section, []):
            compare(item, section)

    if problems:
        for problem in problems:
            print(f"  MISMATCH  {problem}", file=sys.stderr)
        print(
            f"\nprovenance: FAILED - {len(problems)} of {checked} file(s) do not match "
            f"{path.name}. This tree is not the one the record describes.",
            file=sys.stderr,
        )
        return 1

    print(f"provenance: OK ({checked} file(s) match {path.name})")
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--out", default="provenance.json", help="where to write")
    parser.add_argument(
        "--artifact",
        action="append",
        default=[],
        help="a built artifact to hash and record; repeatable",
    )
    parser.add_argument("--tag", default=None, help="release tag, if git cannot tell")
    parser.add_argument(
        "--verify",
        metavar="PATH",
        default=None,
        help="check an existing record against this tree instead of writing one",
    )
    args = parser.parse_args()

    if args.verify:
        return verify(Path(args.verify))

    record = build(args.artifact, args.tag)
    out = Path(args.out)
    out.write_text(render(record), encoding="utf-8")

    dirty = " (dirty tree)" if record["git"]["dirty"] else ""
    tag = record["git"]["tag"] or "no tag"
    print(
        f"provenance: wrote {out} - {len(record['calcs'])} calc(s), "
        f"{record['git']['commit'] or 'no commit'} at {tag}{dirty}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
