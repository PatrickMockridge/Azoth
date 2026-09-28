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

Hashing the file as committed avoids inventing a canonical form for TOML and
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
    python tools/provenance.py --tag v0.1.1 --artifact dist/*.whl
    python tools/provenance.py --verify provenance.json # check a record
"""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
import tomllib
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent.parent
SPEC_DIR = ROOT / "specs" / "calcs"
DATA_DIR = ROOT / "data"

SCHEMA_VERSION = 1

#: Files that affect every calculation's result. Hashed once here rather than
#: once per calc, because repeating them would suggest they are per-calc.
SHARED = (
    "crates/azoth-core/src/error.rs",
    "crates/azoth-core/src/range.rs",
    "crates/azoth-core/src/result.rs",
    "crates/azoth-core/src/solver.rs",
    "crates/azoth-core/src/spec.rs",
    "crates/azoth-core/src/units.rs",
    "crates/azoth-core/src/warning.rs",
    # One PyO3 results module for every namespace: the transport types and the
    # introspection functions are shared, and only the per-namespace calc wrappers
    # are not - see NAMESPACE_SUPPORT.
    "crates/azoth-python/src/data.rs",
    "crates/azoth-python/src/results.rs",
    # The Python twins of the Rust modules above, and they are listed for the
    # reason the pair exists: this list used to carry `core/result.rs`,
    # `core/warning.rs` and `core/range.rs` while omitting `core/result.py`,
    # `core/warnings.py` and `core/errors.py`. A change to the Python warning
    # codes - which is what the Keycard work does - moved no recorded hash at all,
    # while the identical change on the Rust side moved one. A mirrored file that
    # is hashed on one side and not the other is a hash that records half a change.
    "python/src/azoth/core/errors.py",
    "python/src/azoth/core/range.py",
    "python/src/azoth/core/result.py",
    "python/src/azoth/core/solver.py",
    "python/src/azoth/core/units.py",
    "python/src/azoth/core/warnings.py",
    "python/src/azoth/_models_gen.py",
    "python/src/azoth/_registry_gen.py",
    "python/src/azoth/_rust_bridge.py",
    # The schemas define the vocabularies the generators emit from - units,
    # warning codes, solver kinds, model kinds. A change there changes what a spec
    # may say, and it would otherwise move no hash because the generated files it
    # feeds are hashed and it is not.
    "specs/schema/calc.schema.json",
    "specs/schema/model.schema.json",
    # The tools are part of the chain: gen_registry.py writes spec_gen.rs, which
    # every calc reads its range checks from, so a change there changes results.
    # provenance.py hashes itself - the hash of the output then depends on the
    # generator, which is the property you want and not a circularity.
    "tools/gen_models.py",
    "tools/gen_registry.py",
    "tools/gen_docs.py",
    "tools/spec_lint.py",
    "tools/provenance.py",
    "tools/check_links.py",
)

#: Files a namespace's calculations share, hashed once each rather than once per
#: calc. Keyed by namespace because they genuinely differ: hydraulics carries a
#: fittings registry and fluid tables, and thermal carries neither. The solver is
#: deliberately *not* here - it affects every namespace, so it sits in SHARED.
#:
#: These used to sit in SHARED, including `crates/azoth-hydraulics/src/spec_gen.rs`
#: - which was already a namespace's file rather than a shared one, and became
#: plainly wrong once each namespace got its own generated tables.
NAMESPACE_SUPPORT = {
    "hydraulics": (
        "crates/azoth-hydraulics/src/fittings.rs",
        "crates/azoth-hydraulics/src/fluids.rs",
        "crates/azoth-hydraulics/src/provenance.rs",
        "crates/azoth-hydraulics/src/spec_gen.rs",
        "crates/azoth-python/src/hydraulics.rs",
    ),
    "thermal": (
        "crates/azoth-thermal/src/spec_gen.rs",
        "crates/azoth-python/src/thermal.rs",
    ),
    "eos": (
        # The mixture layer is the model layer's own arithmetic - components, the
        # van der Waals one-fluid mixing rule, and the mixture fugacity coefficient
        # no registered calc covers. It affects every model in this namespace, so it
        # is recorded once here rather than once per model. `pt_flash.rs` and
        # `pure_saturation.rs` are picked up by the per-model naming convention.
        "crates/azoth-eos/src/mixture.rs",
        # The pressure iteration the two phase-boundary models share. Recorded once
        # here rather than twice under their names, for the same reason `mixture.rs`
        # is: it is one piece of code that both answers depend on.
        "crates/azoth-eos/src/phase_boundary.rs",
        "crates/azoth-eos/src/model_gen.rs",
        # The mixture layer's own arithmetic has no spec to be hashed under - the
        # registry is scalar and there is nowhere in it for a composition vector - so
        # its reduction tests are recorded here rather than under a calc or a model.
        # They are what holds the one piece of this crate that no kernel checks.
        "crates/azoth-eos/tests/mixture.rs",
        "python/tests/eos/test_mixture_layer.py",
        "crates/azoth-eos/src/spec_gen.rs",
        "crates/azoth-python/src/eos.rs",
    ),
}


def _data_files() -> tuple[str, ...]:
    """Every shipped data file, found by walking the data directory.

    Walked rather than listed, because the list that was here was wrong: it named the
    fittings and fluid tables and omitted `data/components/components.csv` and
    `data/components/kij.csv` - the two largest shipped data files and the two every
    `eos` calculation depends on. A provenance record that does not hash them does not
    describe the release, and nothing said so.

    Walking cannot go stale. Whether each file *should* be shipped is a different
    question, and `databank/manifest.toml` is where it is answered: a data file the
    manifest does not declare fails `tools/check_manifest.py`.
    """
    return tuple(sorted(str(path.relative_to(ROOT)) for path in DATA_DIR.rglob("*.csv")))


#: Derived at import; see `_data_files`.
DATA = _data_files()

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


#: The verification vocabulary, most-checked first. The spellings are the ones a
#: `validation/*.json` case already uses, so a status here and a status there read the
#: same way rather than being two dialects for one idea.
VERIFICATION_STATUSES = ("verified", "partially_verified", "source_needed", "unverified")

VALIDATION_DIR = ROOT / "validation"

#: A model's instances, one file per model id. A calc keeps its tests in its own spec;
#: a model's live here, which is the asymmetry `model_cases` exists to absorb.
CASE_DIR = ROOT / "specs" / "cases"


def validation_cases() -> dict[str, list[dict[str, Any]]]:
    """Every validation case, grouped by the id its `calc` key names.

    Walked rather than listed. The cases live under two directories today
    (`validation/eos/` and `validation/crane_tp410/`), and a listing that named one of
    them would silently drop the other's contribution to a status - which is the
    failure mode a hand-kept index has and a walk does not.
    """
    grouped: dict[str, list[dict[str, Any]]] = {}
    for path in sorted(VALIDATION_DIR.rglob("*.json")):
        case = json.loads(path.read_text(encoding="utf-8"))
        calc = case.get("calc")
        if calc:
            grouped.setdefault(calc, []).append(case)
    return grouped


def derive_verification(
    cases: list[dict[str, Any]], tests: list[dict[str, Any]]
) -> tuple[str, int, int, int]:
    """`(status, validation_cases, tests_active, tests_skipped)`.

    **Derived rather than declared, and the spec schema is why.** A calc carries no
    status field - it refuses one - so the honest answer has to be measured from the
    two records that already exist: the tests the spec ships, and the external cases
    in `validation/`. A declared status is a status that can be wrong; this one cannot
    disagree with the evidence it is drawn from.

    `tests` is a calc spec's `[[tests]]` or a model's `[[cases]]`, which is why the
    rule reads `status` leniently: a calc's test declares `active` or `skipped`, and a
    model's case declares neither because it is always run. **A record with no status
    is active**, so the two trees are one rule rather than two that must agree.

    The cascade, first match wins:

    1. **No active test** - nothing exercises it, so nothing can be claimed.
    2. **A source was sought and not found** - a case marked `source_needed`, or a
       test skipped with that reason. This outranks good evidence on purpose:
       `hydraulics.darcy_weisbach` has a `source_needed` case beside two skipped tests
       citing Crane TP-410 Example 3-5 and Perry's Eq. 6-42, and its other evidence
       being sound does not make those attributions resolve. The counts travel beside
       the status so a reader sees which is which.
    3. **Every case verified** - and there is at least one, which is what keeps this
       branch from firing on an empty list.
    4. **Otherwise** - `partially_verified`, which is the normal case rather than a
       defect. Those ids are exercised by both kernels against expectations pinned in
       their own tree, and their `source.standard` is required by the schema, but no
       *independent* oracle is recorded for them. It is not `verified`, because
       nothing outside the spec checks them; and it is not `unverified`, because that
       would be false.
    """
    active = sum(1 for test in tests if test.get("status") != "skipped")
    skipped = sum(1 for test in tests if test.get("status") == "skipped")

    if active == 0:
        return "unverified", len(cases), active, skipped

    sought = any(
        (case.get("source") or {}).get("verification") == "source_needed" for case in cases
    ) or any(
        str(test.get("skip_reason", "")).startswith("source_needed")
        for test in tests
        if test.get("status") == "skipped"
    )
    if sought:
        return "source_needed", len(cases), active, skipped

    verdicts = [(case.get("source") or {}).get("verification") for case in cases]
    if verdicts and all(verdict == "verified" for verdict in verdicts):
        return "verified", len(cases), active, skipped

    return "partially_verified", len(cases), active, skipped


def model_cases(namespace: str, name: str) -> list[dict[str, Any]]:
    """A model's cases, from `specs/cases/`.

    **Not from the model's own spec**, which is where a calc keeps its tests. A model
    is a *procedure* and its spec declares the procedure; the instances it is held to
    are a separate tree, one file per model id, and `tools/gen_models.py` refuses a
    model with none. Reading the wrong tree is silent - it yields an empty list, and
    an empty list reads as "nothing exercises this", which would have been false for
    all 121 models.
    """
    path = CASE_DIR / namespace / f"{name}.toml"
    if not path.is_file():
        return []
    return tomllib.loads(path.read_text(encoding="utf-8")).get("cases", [])


def rust_kernel_path(namespace: str, name: str) -> str:
    """Repo-relative path to a calc's or model's Rust kernel.

    **A naming convention is not enough here, and assuming one was a defect.** The
    unit-operation models live under `crates/azoth-process/src/models/`, so the
    obvious `crates/azoth-<namespace>/src/<name>.rs` does not exist for any of the
    28 `process.*` ids. `describe` records a missing file rather than failing - which
    is right for a test file, whose absence is a fact - so every one of those ids
    recorded its Rust implementation as `present: false` in a record that is signed
    and shipped. The spec is not at fault: it declares `azoth_process::mixer`, and
    `crates/azoth-process/src/lib.rs` re-exports it, so the *module* path is correct
    and only the file's location was guessed.

    Resolution is by looking, in the order the tree is actually laid out: beside the
    crate's other kernels, then under `models/`, then anywhere in the crate. A
    search matching more than one file is left to fall through to the convention, so
    `describe` records it missing - two kernels for one id is a defect a reader needs
    to see rather than a coin to flip.
    """
    crate = ROOT / f"crates/azoth-{namespace}"
    for candidate in (crate / "src" / f"{name}.rs", crate / "src" / "models" / f"{name}.rs"):
        if candidate.is_file():
            return str(candidate.relative_to(ROOT))
    found = [p for p in (crate / "src").rglob(f"{name}.rs") if p.is_file()]
    if len(found) == 1:
        return str(found[0].relative_to(ROOT))
    return f"crates/azoth-{namespace}/src/{name}.rs"


def calc_entry(spec: dict[str, Any], cases: dict[str, list[dict[str, Any]]]) -> dict[str, Any]:
    """Provenance for one calculation.

    `cases` is the whole index rather than this id's slice, because a caller building
    the record walks every spec and would otherwise re-walk `validation/` 192 times.
    """
    name = spec["id"].split(".")[-1]
    namespace = spec["id"].split(".")[0]
    status, case_count, active, skipped = derive_verification(
        cases.get(spec["id"], []), spec.get("tests", [])
    )

    tests = [
        describe(f"python/tests/{namespace}/test_{name}.py"),
        describe(f"crates/{'azoth-' + namespace}/tests/{name}.rs"),
    ]
    return {
        "id": spec["id"],
        "name": spec["name"],
        "equation": spec["equation"],
        "spec": describe(f"specs/calcs/{namespace}/{name}.toml"),
        "code": [
            describe(f"python/src/azoth/{namespace}/reference/{name}.py"),
            describe(rust_kernel_path(namespace, name)),
        ],
        "tests": tests,
        "verification": status,
        "validation_cases": case_count,
        "tests_active": active,
        "tests_skipped": skipped,
    }


def model_entry(spec: dict[str, Any], cases: dict[str, list[dict[str, Any]]]) -> dict[str, Any]:
    """Provenance for one model.

    The same fields a calc entry carries, from a different spec tree. A model's spec
    is what pins its *procedure* down, so hashing it is at least as load-bearing as
    hashing a calc's equation: a changed tolerance or bracket rule changes every
    answer the model returns while its inputs and outputs look identical.
    """
    name = spec["id"].split(".")[-1]
    namespace = spec["id"].split(".")[0]
    status, case_count, active, skipped = derive_verification(
        cases.get(spec["id"], []), model_cases(namespace, name)
    )
    return {
        "id": spec["id"],
        "name": spec["name"],
        "spec": describe(f"specs/models/{namespace}/{name}.toml"),
        "code": [
            describe(f"python/src/azoth/{namespace}/reference/{name}.py"),
            describe(rust_kernel_path(namespace, name)),
        ],
        "tests": [
            describe(f"python/tests/models/test_{name}.py"),
            describe(f"crates/azoth-{namespace}/tests/{name}.rs"),
        ],
        "verification": status,
        "validation_cases": case_count,
        "tests_active": active,
        "tests_skipped": skipped,
    }


def build(artifacts: list[str], tag: str | None) -> dict[str, Any]:
    """Assemble the provenance record."""
    paths = sorted(SPEC_DIR.rglob("*.toml"))
    calcs = [tomllib.loads(p.read_text(encoding="utf-8")) for p in paths]
    calcs.sort(key=lambda c: c["id"])

    model_paths = sorted((ROOT / "specs" / "models").rglob("*.toml"))
    models = [tomllib.loads(p.read_text(encoding="utf-8")) for p in model_paths]
    models.sort(key=lambda m: m["id"])

    status = git("status", "--porcelain")
    resolved_tag = tag if tag is not None else git("describe", "--tags", "--exact-match")
    index = validation_cases()

    return {
        "schema_version": SCHEMA_VERSION,
        "library": {"name": project_field("name"), "version": version()},
        "git": {
            "commit": git("rev-parse", "HEAD"),
            "tag": resolved_tag,
            # Recorded because provenance for a dirty tree describes something
            # nobody else can reproduce, and the reader deserves to know.
            "dirty": bool(status),
        },
        "calcs": [calc_entry(c, index) for c in calcs],
        "models": [model_entry(m, index) for m in models],
        "shared": [describe(p) for p in SHARED],
        "namespaces": [
            {"namespace": ns, "files": [describe(p) for p in files]}
            for ns, files in sorted(NAMESPACE_SUPPORT.items())
        ],
        "data": [describe(p) for p in DATA],
        "lock_files": [describe(p) for p in LOCK_FILES],
        "licences": [describe(p) for p in LICENCES],
        "artifacts": [describe(a) for a in artifacts],
        "signature": {
            "sigstore_bundle": None,
            "note": (
                "Filled in by .github/workflows/release.yml after keyless signing "
                "with cosign. Null for provenance generated from a source tree, "
                "which is unsigned by construction."
            ),
        },
    }


def project_field(key: str) -> str:
    """One `[project]` field, read from pyproject so it cannot drift.

    **The name as well as the version.** The distribution is `azoth-engine` and the import is
    `azoth`, so the two are different strings and a record that wrote one by hand would be a
    second place the distribution name is stated - which is the drift this reader exists to
    prevent, one field over from where it started.
    """
    text = (ROOT / "pyproject.toml").read_text(encoding="utf-8")
    for line in text.splitlines():
        if line.startswith(f"{key} = "):
            return line.split("=", 1)[1].strip().strip('"')
    return "unknown"


def version() -> str:
    """The library version, read from pyproject so it cannot drift."""
    return project_field("version")


def render(record: dict[str, Any]) -> str:
    """Serialise deterministically: sorted keys, fixed indent, trailing newline."""
    return json.dumps(record, indent=2, sort_keys=True) + "\n"


def verify(path: Path) -> int:
    """Recompute the hashes in a record and report any that no longer match.

    This is what makes the verification instructions executable rather than
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

    for namespace in recorded.get("namespaces", []):
        for item in namespace.get("files", []):
            compare(item, f"namespaces.{namespace['namespace']}")

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
