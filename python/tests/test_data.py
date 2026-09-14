"""Locating the shared data files, and what happens when they are not there.

The happy path is covered everywhere: every calculation that reads a fitting table or
a fluid goes through :func:`azoth._data.find`. What is not covered is the failure, and
the failure is the one a user meets - an installed wheel whose data was not packaged
raises from here, and the message is the whole of what that user has to go on.

The search also has a boundary that nothing asserted: it walks up until it finds the
workspace `Cargo.toml`, so it stops at the repository root rather than continuing
through the filesystem and matching an unrelated directory of the same name. That
boundary is tested below by asking for a path that genuinely exists one level *above*
the checkout.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from azoth import _data

REPO_ROOT = Path(__file__).resolve().parents[2]


def test_an_existing_data_file_is_found() -> None:
    found = _data.find("data/fluids/water.csv")

    assert found.is_file(), f"{found} is not a file"
    assert found.read_text(encoding="utf-8"), "the water table is empty"


def test_the_lookup_is_cached() -> None:
    """Two calls are the same path, so the tree is walked once per process."""
    assert _data.find("data/fluids/water.csv") is _data.find("data/fluids/water.csv")


def test_a_missing_file_raises_with_an_actionable_message() -> None:
    with pytest.raises(_data.DataFileNotFoundError) as caught:
        _data.find("data/fluids/this-substance-does-not-exist.csv")

    message = str(caught.value)
    assert "this-substance-does-not-exist.csv" in message, message
    assert "wheel" in message, (
        f"the message must name the likeliest cause - an installed wheel whose data "
        f"was not packaged - and says only: {message}"
    )


def test_the_error_is_a_filenotfounderror() -> None:
    """A caller catching `FileNotFoundError` catches this, which is the point of it."""
    assert issubclass(_data.DataFileNotFoundError, FileNotFoundError)


def test_the_walk_stops_at_the_repository_root() -> None:
    """A path that exists one level above the checkout must not be found.

    `<repo>/../<repo name>/pyproject.toml` is a real file - it is this repository's own
    manifest, seen from its parent. The search must stop at the workspace `Cargo.toml`
    rather than wander up to where that path resolves, or a data lookup would silently
    read from outside the checkout.
    """
    above = f"{REPO_ROOT.name}/pyproject.toml"
    assert (REPO_ROOT.parent / above).is_file(), (
        "this check needs a directory above the checkout that contains a file at the "
        "path being asked for; it found none, so it proves nothing"
    )

    with pytest.raises(_data.DataFileNotFoundError):
        _data.find(above)
