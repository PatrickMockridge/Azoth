"""Where the bundle's evidence lives, and what is absent in an installed wheel.

The external validation cases and ``provenance.json`` are repository trees that a
distribution does not carry (the wheel check forbids them). Both are optional: a bundle
built from a wheel is still a bundle, and its traceability score records the absence
rather than inventing a commit or a case count.

**The module's own path is not always in the checkout.** ``maturin develop`` can leave the
import resolving to a copy under ``site-packages``, where no ``provenance.json`` is above
it, while pytest is still started from the checkout that holds ``validation/``. The walk
therefore tries the module first and the working directory second, so the same test sees
the same case set in CI and locally.
"""

from __future__ import annotations

from pathlib import Path

__all__ = ["repo_root"]


def _walk(start: Path) -> Path | None:
    """The nearest directory at or above `start` that carries the checkout's markers."""
    for parent in (start, *start.parents):
        if (parent / "provenance.json").is_file() and (parent / "pyproject.toml").is_file():
            return parent
    return None


def repo_root() -> Path | None:
    """The checkout this process is using, or ``None`` in an installed wheel."""
    found = _walk(Path(__file__).resolve().parent)
    if found is not None:
        return found
    return _walk(Path.cwd().resolve())
