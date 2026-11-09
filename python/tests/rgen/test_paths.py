"""The checkout marker is found from the module or, failing that, the working directory."""

from __future__ import annotations

from pathlib import Path

from azoth.rgen import paths


def _markers(root: Path) -> None:
    (root / "provenance.json").write_text("{}", encoding="utf-8")
    (root / "pyproject.toml").write_text("", encoding="utf-8")


def test_the_walk_finds_the_nearest_marked_directory(tmp_path: Path) -> None:
    _markers(tmp_path)
    nested = tmp_path / "a" / "b"
    nested.mkdir(parents=True)
    assert paths._walk(nested) == tmp_path


def test_a_tree_without_the_markers_has_no_root(tmp_path: Path) -> None:
    nested = tmp_path / "a"
    nested.mkdir()
    assert paths._walk(nested) is None


def test_this_checkout_is_found_and_carries_the_validation_tree() -> None:
    root = paths.repo_root()
    assert root is not None
    assert (root / "validation" / "eos").is_dir()
