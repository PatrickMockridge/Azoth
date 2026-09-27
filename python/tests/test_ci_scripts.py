"""The checks that only ever run in CI.

`tools/check_links.py` and `tools/provenance.py` are both wired into workflows and
neither has a test, which is the shape of a check that breaks silently: it is not run
by the suite, so the first person to see it fail is whoever pushes next.

What is tested here is the part that can be *wrong* rather than the part that is
obviously right - link resolution, which has three cases that are easy to conflate
(a page that exists, a page that does not, and a link that is not a page at all) -
and the shape of the provenance record, which is consumed by a release rather than by
a reader.
"""

from __future__ import annotations

import importlib
import json
import subprocess
import sys
import tomllib
from pathlib import Path
from types import ModuleType
from typing import Any

REPO_ROOT = Path(__file__).resolve().parents[2]
TOOLS = REPO_ROOT / "tools"


def link_checker() -> ModuleType:
    """`tools/check_links.py`, imported by name.

    The tools are not a package and are outside mypy's `files`, so the import is
    dynamic - the same arrangement `test_prose_lint.py` uses.
    """
    sys.path.insert(0, str(TOOLS))
    try:
        return importlib.import_module("check_links")
    finally:
        sys.path.pop(0)


def check(tmp_path: Path, text: str) -> tuple[int, list[str]]:
    """Run the link checker over one page, rooted at a temporary directory.

    The root is redirected because the error message is relative to it; a page under
    `/tmp` against the real root would raise from `relative_to` rather than report.
    """
    # `Any` at the point the tool is reached, rather than a widened signature:
    # `ModuleType` declares no `ROOT` or `check_file`, and the tools are outside
    # mypy's `files`.
    tool: Any = link_checker()
    tool.ROOT = tmp_path
    page = tmp_path / "page.md"
    page.write_text(text, encoding="utf-8")
    errors: list[str] = []
    checked: int = tool.check_file(page, errors)
    return checked, errors


def test_a_link_to_a_page_that_exists_is_checked_and_passes(tmp_path: Path) -> None:
    (tmp_path / "target.md").write_text("# Target\n", encoding="utf-8")

    checked, errors = check(tmp_path, "see [the target](./target.md)\n")

    assert checked == 1
    assert errors == []


def test_a_link_to_a_page_that_does_not_exist_is_reported(tmp_path: Path) -> None:
    checked, errors = check(tmp_path, "see [nothing](./nowhere.md)\n")

    assert checked == 1, "a broken link is still a link the checker examined"
    assert len(errors) == 1, errors
    assert "nowhere.md" in errors[0], errors[0]
    assert "page.md:1" in errors[0], f"the line number is part of the report: {errors[0]}"


def test_a_link_into_a_directory_resolves_to_its_index(tmp_path: Path) -> None:
    """`./guide/` is a page in a book, and mdBook renders its `index.md`."""
    (tmp_path / "guide").mkdir()
    (tmp_path / "guide" / "index.md").write_text("# Guide\n", encoding="utf-8")

    checked, errors = check(tmp_path, "see [the guide](./guide/)\n")

    assert checked == 1
    assert errors == []


def test_external_links_and_bare_anchors_are_not_resolved(tmp_path: Path) -> None:
    """Two different cases, and the count distinguishes them.

    An external link is skipped outright - fetching it would make the docs build
    depend on someone else's server. A bare anchor is *examined* and then needs no
    file, because it points at a heading on the page it is written on. Neither can
    produce an error, and only the second is counted as a link looked at.
    """
    checked, errors = check(
        tmp_path,
        "see [the docs](https://example.invalid/x) and [below](#a-heading)\n",
    )

    assert checked == 1, "the same-page anchor is examined; the external link is skipped"
    assert errors == []


def test_a_link_with_an_anchor_resolves_on_its_file(tmp_path: Path) -> None:
    (tmp_path / "target.md").write_text("# Target\n", encoding="utf-8")

    checked, errors = check(tmp_path, "see [a heading](./target.md#a-heading)\n")

    assert checked == 1
    assert errors == []


def version_checker() -> ModuleType:
    """`tools/check_versions.py`, imported by name, for the reason `link_checker` is."""
    sys.path.insert(0, str(TOOLS))
    try:
        return importlib.import_module("check_versions")
    finally:
        sys.path.pop(0)


def _version_tree(tmp_path: Path, cargo: str, pyproject: str, crates: dict[str, str]) -> Path:
    """A repository just large enough for the check to read: two manifests and some crates."""
    (tmp_path / "Cargo.toml").write_text(cargo, encoding="utf-8")
    (tmp_path / "pyproject.toml").write_text(pyproject, encoding="utf-8")
    for name, text in crates.items():
        crate = tmp_path / "crates" / name
        crate.mkdir(parents=True)
        (crate / "Cargo.toml").write_text(text, encoding="utf-8")
    return tmp_path


CARGO = '[workspace]\n[workspace.package]\nversion = "0.1.0"\n'
PYPROJECT = '[project]\nname = "azoth-engine"\nversion = "0.1.0"\n'
TAKES_IT = '[package]\nname = "azoth-core"\nversion.workspace = true\n'


def test_two_versions_that_agree_are_no_problems(tmp_path: Path) -> None:
    root = _version_tree(tmp_path, CARGO, PYPROJECT, {"azoth-core": TAKES_IT})
    assert version_checker().problems(root) == []


def test_a_bumped_wheel_and_an_unbumped_binary_are_reported(tmp_path: Path) -> None:
    """**The release this exists to stop**: a wheel published as one version whose own MCP server
    reports another, with nothing in the build able to notice."""
    root = _version_tree(
        tmp_path,
        CARGO,
        '[project]\nname = "azoth-engine"\nversion = "0.2.0"\n',
        {"azoth-core": TAKES_IT},
    )
    (problems,) = version_checker().problems(root)
    assert "0.2.0" in problems and "0.1.0" in problems


def test_a_crate_that_pins_its_own_version_is_reported(tmp_path: Path) -> None:
    """The pair above cannot speak for a crate that states its own number, so it must not."""
    root = _version_tree(
        tmp_path,
        CARGO,
        PYPROJECT,
        {"azoth-core": '[package]\nname = "azoth-core"\nversion = "0.1.0"\n'},
    )
    (problems,) = version_checker().problems(root)
    assert "version.workspace = true" in problems


def test_a_manifest_with_no_version_is_reported_rather_than_assumed(tmp_path: Path) -> None:
    root = _version_tree(tmp_path, "[workspace]\n", PYPROJECT, {})
    (problems,) = version_checker().problems(root)
    assert "no `[workspace.package] version`" in problems


def test_the_version_gate_passes_on_this_tree() -> None:
    """The gate itself, run where it runs in CI."""
    result = subprocess.run(
        [sys.executable, str(TOOLS / "check_versions.py")],
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 0, result.stderr
    assert "OK" in result.stdout


def test_the_link_checker_passes_on_this_tree() -> None:
    """The gate itself, run where it runs in CI."""
    result = subprocess.run(
        [sys.executable, str(TOOLS / "check_links.py")],
        capture_output=True,
        text=True,
        check=False,
    )

    assert result.returncode == 0, f"{result.stdout}{result.stderr}"


def test_provenance_records_a_verifiable_shape(tmp_path: Path) -> None:
    """The record a release is made from, written and read back."""
    out = tmp_path / "provenance.json"

    result = subprocess.run(
        [sys.executable, str(TOOLS / "provenance.py"), "--out", str(out)],
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 0, f"{result.stdout}{result.stderr}"

    record = json.loads(out.read_text(encoding="utf-8"))
    for key in ("schema_version", "git", "calcs", "models", "data", "licences"):
        assert key in record, f"the record has no {key!r}: {sorted(record)}"
    assert record["git"]["commit"], "a record with no commit names no revision"
    assert record["calcs"], "the record lists no calculations"

    # And it must verify against the tree it was written from, or the check the
    # release job runs would fail on the record that job just produced.
    verified = subprocess.run(
        [sys.executable, str(TOOLS / "provenance.py"), "--verify", str(out)],
        capture_output=True,
        text=True,
        check=False,
    )
    assert verified.returncode == 0, f"{verified.stdout}{verified.stderr}"


def provenance_tool() -> ModuleType:
    """`tools/provenance.py`, imported by name, for the same reason as the link checker."""
    sys.path.insert(0, str(TOOLS))
    try:
        return importlib.import_module("provenance")
    finally:
        sys.path.pop(0)


def test_every_registered_id_resolves_to_its_spec_and_kernels() -> None:
    """The record's paths are real, for every id.

    `describe` records a missing file rather than failing, which is right for a test
    file whose absence is a fact - but it also means a *derivation* that is wrong
    reads as an absence rather than as a defect. One was: the Rust kernel was assumed
    to sit beside the crate's other kernels, and the 28 unit-operation models live
    under `crates/azoth-process/src/models/`, so every `process.*` id recorded its
    Rust implementation as `present: false` in a record that is signed and shipped.
    Nothing caught it, because nothing asserted that a path resolves.
    """
    provenance = provenance_tool()
    index = provenance.validation_cases()
    entries = (("specs/calcs", provenance.calc_entry), ("specs/models", provenance.model_entry))
    checked = 0
    for tree, entry in entries:
        for path in sorted((REPO_ROOT / tree).rglob("*.toml")):
            spec = tomllib.loads(path.read_text(encoding="utf-8"))
            record = entry(spec, index)
            checked += 1
            assert record["spec"]["present"], (
                f"{record['id']}: the record names a spec that is not there: "
                f"{record['spec']['path']}"
            )
            missing = [code["path"] for code in record["code"] if not code["present"]]
            assert not missing, f"{record['id']}: the record names no kernel at {missing}"
    # Non-vacuity: a walk that found nothing would pass every assertion above without
    # having checked anything, which is the failure this file exists to catch. Held to
    # the registry rather than to a number typed here, so it cannot go stale.
    from azoth._models_gen import MODELS
    from azoth._registry_gen import CALCS

    assert checked == len(CALCS) + len(MODELS), (
        f"walked {checked} specs but the registry names {len(CALCS) + len(MODELS)}"
    )


def test_the_derived_status_follows_the_evidence() -> None:
    """The four branches of the cascade, each from its own evidence.

    A pure function, so it is worth pinning directly rather than only through the
    tree: the interesting failures are the branches, and a tree walk exercises one
    path per id.

    The leniency of `status` is part of the contract, not an accident of the data. A
    calc's test declares `active` or `skipped`; a model's case declares neither,
    because it is always run. Reading that wrong is silent - it yields "no active
    test" for all 121 models, which is what the first version of this did.
    """
    provenance = provenance_tool()
    derive = provenance.derive_verification
    verified = [{"source": {"verification": "verified"}}]

    assert derive([], [])[0] == "unverified"
    assert derive(verified, [{"id": "a", "status": "skipped"}])[0] == "unverified"

    # A model's case carries no status and is therefore run.
    assert derive([], [{"id": "a"}])[:1] == ("partially_verified",)
    assert derive(verified, [{"id": "a"}])[0] == "verified"

    # `source_needed` outranks good evidence, from either tree. An active test is
    # present in both, so rule 1 does not decide it first.
    assert derive([{"source": {"verification": "source_needed"}}], verified)[0] == "source_needed"
    assert derive(
        verified,
        [
            {"id": "a"},
            {"id": "b", "status": "skipped", "skip_reason": "source_needed: not located"},
        ],
    )[0] == ("source_needed")

    # And an unrelated skip reason is not evidence of a missing source.
    assert (
        derive(
            verified,
            [{"id": "a"}, {"id": "b", "status": "skipped", "skip_reason": "dimension mismatch"}],
        )[0]
        == "verified"
    )

    # The counts travel with the status.
    both = [{"id": "a"}, {"id": "b", "status": "skipped"}]
    assert derive(verified, both) == ("verified", 1, 1, 1)
