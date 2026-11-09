"""Where the bundle's evidence lives, and what is absent in an installed wheel.

The external validation cases and ``provenance.json`` are repository trees that a
distribution does not carry (the wheel check forbids them). Both are optional: a bundle
built from a wheel is still a bundle, and its traceability score records the absence
rather than inventing a commit or a case count.
"""

from __future__ import annotations

from pathlib import Path

__all__ = ["repo_root"]


def repo_root() -> Path | None:
    """The checkout this module is running from, or ``None`` in a wheel."""
    for parent in Path(__file__).resolve().parents:
        if (parent / "provenance.json").is_file() and (parent / "pyproject.toml").is_file():
            return parent
    return None
