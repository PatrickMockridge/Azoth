"""Locating the shared data files.

The data lives at the repository root (``data/``) rather than inside the Python
package, because the Rust side embeds the same bytes with ``include_str!``. One
file, two languages.

One file is where the arrangement stops. The two parsers are not compared against
each other: the compiled extension exposes the calcs but not the tables it parses,
so a test on this side has no Rust-parsed value to check. That is a known gap, not
a checked property, and it is worth stating because the whole reason this library
keeps two implementations is to have them disagree when one is wrong.

That arrangement means the Python side has to *find* it, and where it lives
depends on how azoth was obtained:

1. **Installed wheel** - the files are packaged alongside the module. Not wired
   up yet; see the note below.
2. **Source checkout or editable install** - walk up from this file until the
   repository root, identified by the workspace ``Cargo.toml``.

``TODO``: the wheel path (1) needs a maturin ``include`` entry so the CSV ships
inside the distribution. Until then this resolves only from a source tree, which
is what the test suite and the CLI both run from. Recorded here rather than left
implicit because a library that silently cannot find its own data after
installation is a nastier failure than one that says so.
"""

from __future__ import annotations

from functools import cache
from importlib import resources
from pathlib import Path


class DataFileNotFoundError(FileNotFoundError):
    """A shared data file could not be located."""


def _packaged(name: str) -> Path | None:
    """A copy shipped inside the package, if one was installed."""
    try:
        candidate = resources.files("azoth") / "_data" / name
    except (ModuleNotFoundError, TypeError):  # pragma: no cover - packaging edge
        return None
    return Path(str(candidate)) if candidate.is_file() else None


@cache
def find(relative: str) -> Path:
    """Locate a repository data file by its repo-relative path.

    Cached: the files are immutable within a process, and every call would
    otherwise re-walk the tree.
    """
    packaged = _packaged(Path(relative).name)
    if packaged is not None:
        return packaged

    # Walk up looking for the repository root. The workspace Cargo.toml is the
    # marker, so the search stops at the repo boundary instead of wandering up
    # through the filesystem and matching some unrelated directory of the same
    # name.
    for parent in Path(__file__).resolve().parents:
        candidate = parent / relative
        if candidate.is_file():
            return candidate
        if (parent / "Cargo.toml").is_file():
            break

    raise DataFileNotFoundError(
        f"could not locate {relative!r}. This build of azoth resolves data files "
        f"from a source checkout; if it was installed as a wheel, the data was not "
        f"packaged - see azoth/_data.py."
    )
