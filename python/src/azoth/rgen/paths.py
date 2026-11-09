"""Where the bundle's evidence lives, and what is absent in an installed wheel.

The external validation cases are a repository tree that a distribution does not carry
(the wheel check forbids it). They are optional: a bundle built from a wheel is still a
bundle, and its traceability score records the absence rather than inventing a commit or
a case count.

**The checkout is marked by its source tree, not by a release record.** ``provenance.json``
is generated per release and gitignored, so no clone and no CI checkout carries it; a walk
that asked for it found nothing there and counted four cases instead of five. The markers
are the package source, the ``validation/`` tree and the root ``pyproject.toml``, which a
clone has and a wheel's ``site-packages`` does not. ``runner.build_identity`` still reads
``provenance.json`` when it is present, and records its absence when it is not.

**The module's own path is not always in the checkout.** ``maturin develop`` can leave the
import resolving to a copy under ``site-packages``, while pytest is still started from the
checkout that holds ``validation/``. The walk therefore tries the module first and the
working directory second, so the same test sees the same case set in CI and locally.
"""

from __future__ import annotations

from pathlib import Path

__all__ = ["repo_root"]


def _is_checkout(path: Path) -> bool:
    """Whether `path` is an azoth source checkout.

    All three together: the root `pyproject.toml`, the external `validation/` tree the
    bundle's cases come from, and the package source. A project that merely has a
    `validation/` directory, or a wheel's `site-packages`, cannot satisfy the set.
    """
    return (
        (path / "pyproject.toml").is_file()
        and (path / "validation").is_dir()
        and (path / "python" / "src" / "azoth" / "__init__.py").is_file()
    )


def _walk(start: Path) -> Path | None:
    """The nearest directory at or above `start` that carries the checkout's markers."""
    for parent in (start, *start.parents):
        if _is_checkout(parent):
            return parent
    return None


def repo_root() -> Path | None:
    """The checkout this process is using, or ``None`` in an installed wheel."""
    found = _walk(Path(__file__).resolve().parent)
    if found is not None:
        return found
    return _walk(Path.cwd().resolve())
