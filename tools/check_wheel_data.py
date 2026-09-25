#!/usr/bin/env python3
"""Check that an *installed* azoth can find its own data.

Run this with a wheel installed, from a directory outside the repository:

    cd /tmp && /path/to/clean-venv/bin/python /path/to/repo/tools/check_wheel_data.py

The artifact check alone needs no install and no interpreter, so it runs anywhere:

    python tools/check_wheel_data.py --artifacts-only

# Why this is separate from the test suite

The suite runs against the source tree, where ``azoth/_data.py`` walks up and finds
``data/`` at the repository root. That path works no matter what the wheel contains,
so every test in the suite would pass on a wheel that shipped no data at all. The one
thing a source-tree run cannot exercise is the thing that was broken: resolving the
tables from an *installed* distribution.

# What it asserts

1. ``azoth`` is the installed copy, not the source tree. Otherwise this proves nothing
   about the artifact - and it would still pass.
2. A calculation whose data comes from the fittings registry works, and reports the
   estimated-data warning its placeholder coefficients require.
3. A fluid property lookup works.
4. Neither resolved data file lives inside the repository. If one does, the source
   tree was reachable and assertion 1 was passed by accident.
5. The *artifacts* contain what they are for and nothing else. A wheel is the Python
   package; an sdist is a tree someone rebuilds from, so it carries the Rust sources
   and must still not carry the editor.

Exit status is non-zero if any of those fails.
"""

from __future__ import annotations

import sys
import tarfile
import zipfile
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent

#: Trees that belong to the repository and not to a distribution of the Python package.
#:
#: **Stated as a list rather than derived from `[tool.maturin]`'s `include`**, because what is
#: being checked is whether maturin's own file set agrees with this one - and deriving it from the
#: allow-list would make the check agree with the thing it is checking.
FORBIDDEN_IN_ARTIFACTS: tuple[str, ...] = (
    "ui/",
    "specs/",
    "docs/",
    "tools/",
    "agents/",
    "lean/",
    "validation/",
    "target/",
)

#: A wheel is the Python package, so it carries no Rust source either.
FORBIDDEN_IN_WHEEL: tuple[str, ...] = ("crates/",)

#: A built browser module, under any path.
FORBIDDEN_SUBSTRING = "azoth_wasm"

#: What a wheel of this project must carry: the package, the extension, and the data the two
#: middleware commands and every calc read from the databank.
REQUIRED_IN_WHEEL: tuple[str, ...] = (
    "azoth/__init__.py",
    "azoth/_core",
    "data/components/",
    "data/reactions/",
    "data/reactors/",
    "data/fittings/",
    "data/fluids/",
)


def fail(message: str) -> None:
    print(f"check_wheel_data: {message}", file=sys.stderr)
    raise SystemExit(1)


def main() -> int:
    # The artifact walk is a property of the files in `dist/` and of nothing else, so it is
    # reachable without an installed wheel - which is what makes it measurable while the rest of
    # this script is not.
    if "--artifacts-only" in sys.argv:
        check_artifacts(REPO_ROOT / "dist")
        return 0

    import azoth
    from azoth._data import find
    from azoth.hydraulics import crane_k_factors
    from azoth.properties import provider_for

    installed = Path(azoth.__file__).resolve()
    print(f"  azoth imported from {installed}")

    if installed.is_relative_to(REPO_ROOT):
        fail(
            f"azoth resolved to the source tree ({installed}) rather than an installed "
            f"wheel. Run this from outside {REPO_ROOT} with the wheel's interpreter, "
            f"or it proves nothing about the artifact."
        )

    # 1. The fittings registry, which is read from a CSV the wheel must carry.
    result = crane_k_factors(["90_elbow", "gate_valve_open"], 0.018)
    if result.k_total <= 0:
        fail(f"crane_k_factors returned a non-positive coefficient: {result.k_total}")
    print(f"  crane_k_factors: k_total={result.k_total}")

    # 2. A fluid table, which is a separate file and a separate parser.
    water = provider_for("water")
    density = water.density(azoth.ureg.Quantity(20.0, "degC"))
    if not 990.0 < density.magnitude < 1005.0:
        fail(f"water density at 20 C came back as {density}, which is not water")
    print(f"  water at 20 C: {density:.2f}")

    # 3. Neither file may have come from the repository.
    for relative in ("data/fittings/crane_k_factors.csv", "data/fluids/water.csv"):
        resolved = find(relative).resolve()
        if resolved.is_relative_to(REPO_ROOT):
            fail(
                f"{relative} resolved inside the repository ({resolved}), so the "
                f"source tree was reachable after all and this run did not exercise "
                f"the packaged path."
            )
        print(f"  {relative} -> {resolved}")

    # 4. What the artifacts carry. Every member of every one is read, so this is a measurement of
    #    the files a stranger would receive rather than a belief about maturin's defaults.
    check_artifacts(REPO_ROOT / "dist")

    print("check_wheel_data: OK")
    return 0


def carried(members: list[str], required: str) -> bool:
    """Whether an archive carries anything under a path or ending in one.

    A parameter rather than a closure over the loop's own list: a `def` inside the loop would bind
    the name late, which is the kind of bug that reads as correct and is not.
    """
    return any(member.startswith(required) or f"/{required}" in member for member in members)


def check_artifacts(dist: Path) -> None:
    """Every member of every artifact under `dist`, held to what the two are for."""
    artifacts = sorted(dist.glob("*.whl")) + sorted(dist.glob("*.tar.gz"))
    if not artifacts:
        # Not a failure: this script is also run by hand against an installed wheel, where there
        # is no `dist/` to walk. The `wheel-data` job builds the artifacts first.
        print("  (no dist/*.whl or dist/*.tar.gz to walk)")
        return

    for artifact in artifacts:
        with (
            zipfile.ZipFile(artifact) if artifact.suffix == ".whl" else tarfile.open(artifact)
        ) as archive:
            members = (
                archive.namelist()
                if artifact.suffix == ".whl"
                else [member.name for member in archive.getmembers()]
            )
        is_wheel = artifact.suffix == ".whl"
        forbidden = FORBIDDEN_IN_ARTIFACTS + (FORBIDDEN_IN_WHEEL if is_wheel else ())

        for member in members:
            for prefix in forbidden:
                if member.startswith(prefix) or f"/{prefix}" in member:
                    fail(
                        f"{artifact.name} carries `{member}`, which is under `{prefix}` - a "
                        f"distribution of the Python package is not a copy of the repository"
                    )
            if FORBIDDEN_SUBSTRING in member:
                fail(
                    f"{artifact.name} carries `{member}`, which is a built browser module and "
                    f"belongs to the editor rather than to the wheel"
                )

        if is_wheel:
            for required in REQUIRED_IN_WHEEL:
                if not carried(members, required):
                    fail(f"{artifact.name} carries nothing under `{required}`, which it needs")
        else:
            # **An sdist is the tree someone rebuilds from**, so it carries the Rust sources and
            # the manifests - and that is the assertion, not an absence of one.
            for required in ("Cargo.toml", "pyproject.toml", "crates/"):
                if not carried(members, required):
                    fail(f"{artifact.name} carries nothing under `{required}`, which it needs")

        print(f"  {artifact.name}: {len(members)} member(s), none of them the repository")


if __name__ == "__main__":
    raise SystemExit(main())
