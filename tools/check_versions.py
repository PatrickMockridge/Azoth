#!/usr/bin/env python3
"""The version numbers this repository states, held to each other.

Two of them are load-bearing and neither sees the other:

* `pyproject.toml`'s `[project] version` is what a wheel is published as, and what
  `tools/provenance.py` stamps into the release record.
* `Cargo.toml`'s `[workspace.package] version` is what every binary reports, and every crate
  takes it by `version.workspace = true`. It reaches a client as `azoth mcp`'s
  `serverInfo.version`, which is `env!("CARGO_PKG_VERSION")`.

So a release that bumps one and not the other ships a wheel whose own MCP server tells a client a
different number from the one on the index - and nothing in the build notices, because both numbers
are correct for the thing they were written in. That is what this gate is for.

**Read, not generated.** A version is bumped by hand in two files because a version *is* a
hand-written fact - the generator half of this repository compiles declarations that already exist,
and there is no file a version could be generated *from*. What can be had instead is a check that
the two hands agreed, which is this.

    python tools/check_versions.py            # check
"""

from __future__ import annotations

import sys
import tomllib
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent


def problems(root: Path) -> list[str]:
    """Every disagreement about a version, as sentences.

    A list rather than an exit, so the check can be tested without a repository around it - the
    part that can be *wrong* is the reading, and the exit code is not.
    """
    found: list[str] = []

    cargo = tomllib.loads((root / "Cargo.toml").read_text(encoding="utf-8"))
    workspace = cargo.get("workspace", {}).get("package", {}).get("version")
    if workspace is None:
        found.append("Cargo.toml states no `[workspace.package] version`")
        return found

    pyproject = tomllib.loads((root / "pyproject.toml").read_text(encoding="utf-8"))
    published = pyproject.get("project", {}).get("version")
    if published is None:
        found.append("pyproject.toml states no `[project] version`")
        return found

    if workspace != published:
        found.append(
            f"the wheel is published as {published} and every binary reports {workspace}; a "
            f"release would ship a server that names a different version from the index"
        )

    # A crate that states its own version is one the pair above cannot speak for: it would ship a
    # number neither file names, and `version.workspace = true` is what makes the pair exhaustive.
    for manifest in sorted((root / "crates").glob("*/Cargo.toml")):
        package = tomllib.loads(manifest.read_text(encoding="utf-8")).get("package", {})
        version = package.get("version")
        if isinstance(version, dict) and "workspace" in version:
            continue
        found.append(
            f"{manifest.relative_to(root)} states its own version ({version!r}) rather than taking "
            f"the workspace's; use `version.workspace = true`"
        )
    return found


def main() -> None:
    found = problems(ROOT)
    if found:
        for problem in found:
            print(f"check_versions: {problem}", file=sys.stderr)
        raise SystemExit(1)
    cargo = tomllib.loads((ROOT / "Cargo.toml").read_text(encoding="utf-8"))
    print(
        f"check_versions: OK ({cargo['workspace']['package']['version']}, and every crate takes it)"
    )


if __name__ == "__main__":
    main()
