"""Locating the shared data files.

The data lives at the repository root (``data/``) rather than inside the Python
package, because the Rust side embeds the same bytes with ``include_str!``. One
file, two languages - and the two sides *are* compared:
``python/tests/test_data_agreement.py`` checks every parsed field row by row and
checks that the two sides read byte-identical files. The byte check is the stronger
of the two: parsed values agreeing does not prove the same file was read, and a stale
copy bundled into a wheel would look exactly like an agreement that proves nothing.

That arrangement means the Python side has to *find* it, and where it lives depends
on how azoth was obtained:

1. **Installed wheel** - the files are packaged by a maturin ``include`` entry. Note
   where they land: ``include`` copies a file to its path relative to the project
   root and cannot be redirected, so a wheel carries ``data/fluids/water.csv`` at the
   *distribution* root, not under the ``azoth`` package. That is why the search below
   looks in two places rather than one.
2. **Source checkout or editable install** - walk up from this file until the
   repository root, identified by the workspace ``Cargo.toml``.

Records were once shipped under ``azoth/_data/``, and that path is still searched
first in case a future packaging arrangement puts them there.
"""

from __future__ import annotations

from functools import cache
from importlib import resources
from pathlib import Path


class DataFileNotFoundError(FileNotFoundError):
    """A shared data file could not be located."""


def _packaged(relative: str) -> Path | None:
    """A copy shipped inside the distribution, if one was installed.

    Two locations, because the two ways of shipping data land in different places.
    ``importlib.resources.files("azoth")`` is the package directory and its parent is
    the distribution root - ``site-packages`` for an installed wheel - so a wheel
    built with a maturin ``include`` carries ``data/fluids/water.csv`` there. A file
    placed inside the package would be under ``azoth/_data/`` instead, which is where
    an earlier arrangement expected them.
    """
    try:
        package = resources.files("azoth")
    except (ModuleNotFoundError, TypeError):  # pragma: no cover - packaging edge
        return None
    # `Traversable` has no `parent`, so the distribution root is reached through a
    # real path. For an ordinary install that path is the filesystem, which is the
    # only case where there is a root to reach: an imported-from-archive layout has
    # no directory above the package and falls through to the walk-up below.
    root = Path(str(package)).parent
    for candidate in (Path(str(package)) / "_data" / Path(relative).name, root / relative):
        if candidate.is_file():
            return candidate
    return None


@cache
def find(relative: str) -> Path:
    """Locate a repository data file by its repo-relative path.

    Cached: the files are immutable within a process, and every call would
    otherwise re-walk the tree.
    """
    packaged = _packaged(relative)
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
