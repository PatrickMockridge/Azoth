"""Is the built extension the Rust that is in the tree?

**A stale build is not "Rust is unavailable".** `conftest.py` turns a *missing* extension into
a skip, or into a strict xfail under `AZOTH_REQUIRE_RUST=1` - and both are right for a
developer who has not built it yet. A build that exists and is *older than the sources* is a
different thing: the tests run, they pass, and they are testing code nobody wrote. That is the
trap recorded as `worktree-maturin-after-merge`, where one stale build accounted for
thirty-five of thirty-six sweep failures.

mtime is the right measure here and not a guess: a merge or a checkout rewrites the `.rs`
files *now*, so the direction of the comparison is always the safe one - the sources can only
be newer than the build when something has happened to them since it was made.

The decision is a pure function of two floats so the three states can be pinned without a
build; only `state()` reads the tree.
"""

from __future__ import annotations

from pathlib import Path

#: The three answers, as words rather than booleans, because "absent" and "stale" call for
#: different behaviour and a caller that collapsed them would turn a stale build into a skip.
ABSENT, STALE, CURRENT = "absent", "stale", "current"

#: Where the extension lands. `maturin develop` puts it beside the package's own modules.
EXTENSION = Path("python") / "src" / "azoth" / "_core.abi3.so"

#: What makes a build out of date. `Cargo.toml` as well as the sources: a dependency or a
#: feature flag changes the build as surely as a `fn` body does.
#:
#: **`crates/*/src/`, not all of `crates/`.** An integration test under `crates/*/tests/` is a
#: separate binary that links the library; editing one does not change what the extension is,
#: and counting it made this check fire on the very change that added it.
SOURCES = ("crates/*/src/**/*.rs", "crates/*/Cargo.toml", "Cargo.toml")

#: Generator output under `crates/`, which does **not** make a build out of date in the sense
#: this checks. `gen_registry` rewrites `provenance_gen.rs` on any spec edit, and calling that
#: a stale build would mean a `maturin develop` after every regeneration - noise that would
#: teach a reader to ignore the message. What the check is for is a hand-written kernel that a
#: merge rewrote while an old binary sat in place.
GENERATED = "_gen.rs"


def staleness(newest_source: float | None, extension: float | None) -> str:
    """The pure decision: two mtimes, one word.

    `None` for the extension means it is not there - a different state with a different
    meaning, and the reason this is not a boolean.
    """
    if extension is None:
        return ABSENT
    if newest_source is None:
        # No sources in the tree at all: a checkout of the Python package alone. Nothing to be
        # stale against.
        return CURRENT
    return STALE if newest_source > extension else CURRENT


def newest_source_mtime(root: Path) -> float | None:
    """The newest mtime over every file a build depends on."""
    newest: float | None = None
    for pattern in SOURCES:
        for path in root.glob(pattern):
            if not path.is_file() or path.name.endswith(GENERATED):
                continue
            modified = path.stat().st_mtime
            newest = modified if newest is None else max(newest, modified)
    return newest


def extension_mtime(root: Path) -> float | None:
    """The extension's mtime, or `None` when it has not been built here."""
    path = root / EXTENSION
    return path.stat().st_mtime if path.exists() else None


def state(root: Path) -> str:
    """Which of the three this checkout is in."""
    return staleness(newest_source_mtime(root), extension_mtime(root))
