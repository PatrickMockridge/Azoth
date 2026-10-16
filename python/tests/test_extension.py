"""The three states of the built extension, and which one the tree is in.

`conftest.py` skips, or strictly xfails, when the extension is *absent* - both correct for a
developer who has not built it. It stops the run when the extension is **stale**, because a
build older than the sources means every Rust-backed test is checking code nobody wrote. The
difference is the whole reason `_extension` exists, and it is one comparison, so it is pinned
here rather than discovered in a session.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from _extension import ABSENT, CURRENT, STALE, newest_source_mtime, staleness, state

ROOT = Path(__file__).resolve().parents[2]


def test_an_absent_extension_is_absent_rather_than_stale() -> None:
    """The distinction the whole module exists for: no build is a skip, an old build is a stop."""
    assert staleness(newest_source=100.0, extension=None) == ABSENT


def test_a_build_older_than_the_sources_is_stale() -> None:
    assert staleness(newest_source=200.0, extension=100.0) == STALE


def test_a_build_newer_than_the_sources_is_current() -> None:
    assert staleness(newest_source=100.0, extension=200.0) == CURRENT


def test_a_build_exactly_as_old_as_its_newest_source_is_current() -> None:
    """The boundary is `>` rather than `>=`: a build made *from* a source is not stale against
    it, and a `>=` would stop the run on every build whose clock ticked the same second."""
    assert staleness(newest_source=100.0, extension=100.0) == CURRENT


def test_sources_that_are_not_there_are_not_a_reason_to_stop() -> None:
    """A checkout of the Python half alone has nothing to be stale against."""
    assert staleness(newest_source=None, extension=100.0) == CURRENT


@pytest.mark.requires_rust
def test_this_tree_is_not_stale() -> None:
    """**The check the suite makes of itself.** If this fails, the Rust-backed tests above it
    were run against a build nobody wrote - run `maturin develop` and run it again."""
    assert state(ROOT) != STALE, (
        "the extension is older than the newest Rust source; run `maturin develop`"
    )


def test_the_sources_it_reads_are_the_ones_a_build_depends_on() -> None:
    """A guard on the guard: if the globs stop matching, every comparison silently passes."""
    newest = newest_source_mtime(ROOT)
    assert newest is not None, "no Rust sources matched; the SOURCES globs have drifted"
