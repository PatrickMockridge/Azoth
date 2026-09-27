"""The skills catalog's basis gate: a declaration held to the prose it qualifies.

The catalog says how a skill's numbers are produced and the `SKILL.md` says it in prose, and
only one of the two is read by an agent deciding what to do. These tests pin the rule that
holds them together, and the two ways it is allowed to be silent: a basis that makes no claim
either way, and a skill that names an id *in order to disclaim it*.

The rule is pure over a file, so it is tested directly rather than through the catalog - the
catalog is the data it runs on, and a test that went through it would be asserting today's
skills rather than the rule.
"""

from __future__ import annotations

import importlib
import sys
from pathlib import Path
from types import ModuleType
from typing import Any

REPO_ROOT = Path(__file__).resolve().parents[2]
TOOLS = REPO_ROOT / "tools"

#: Two ids, one a strict prefix of the other, which is a real pair in this registry and the
#: case a substring test gets wrong.
IDS = ["eos.pt_flash", "eos.pt_flash_saft", "hydraulics.darcy_weisbach"]


def validator() -> ModuleType:
    """`tools/validate_skills.py`, imported by name - the tools are not a package."""
    sys.path.insert(0, str(TOOLS))
    try:
        return importlib.import_module("validate_skills")
    finally:
        sys.path.pop(0)


def problems(basis: str, body: str, tmp_path: Path, **extra: Any) -> list[str]:
    """The rule, run over one synthetic skill."""
    page = tmp_path / "SKILL.md"
    page.write_text(body, encoding="utf-8")
    entry = {"name": "azoth-something", "calculation_basis": basis, **extra}
    # Annotated rather than returned directly: the tool is outside mypy's `files`, so reaching
    # it by name resolves to `Any`, and returning that from a typed function is the implicit
    # `Any` mypy's strict mode refuses at exactly the boundary worth keeping typed.
    found: list[str] = validator()._basis_problems(entry, page, "skills.toml skill[0]", IDS)
    return found


def test_a_screening_skill_citing_a_calculation_is_reported(tmp_path: Path) -> None:
    found = problems("screening", "It calls `eos.pt_flash` and returns its answer.\n", tmp_path)

    assert len(found) == 1, found
    assert "screening" in found[0] and "eos.pt_flash" in found[0]


def test_an_azoth_skill_naming_no_calculation_is_reported(tmp_path: Path) -> None:
    found = problems("azoth", "It drives the library.\n", tmp_path)

    assert len(found) == 1, found
    assert "azoth" in found[0]


def test_a_screening_skill_that_disclaims_an_id_is_exempt(tmp_path: Path) -> None:
    """The honest case, and the reason the exemption is declared rather than detected.

    "azoth's `process.compressor` is an isentropic step over an efficiency" names an id in
    order to say the skill is *not* it. No syntactic rule separates that from a citation of a
    calculation the skill does call, so the judgement is declared - with a reason - beside the
    basis it qualifies.
    """
    body = "azoth's `eos.pt_flash` is an equilibrium flash, which is not this skill.\n"

    assert problems("screening", body, tmp_path) != []
    assert problems("screening", body, tmp_path, id_references_ok="disclaims it") == []


def test_a_basis_that_claims_nothing_is_not_checked(tmp_path: Path) -> None:
    """`data-retrieval` returns data and `hybrid` calls the library among other things."""
    body = "It returns `hydraulics.darcy_weisbach`'s answer and some data.\n"

    for basis in ("data-retrieval", "hybrid"):
        assert problems(basis, body, tmp_path) == [], basis


def test_an_id_that_is_a_prefix_of_another_is_not_mis_reported(tmp_path: Path) -> None:
    """`eos.pt_flash` is a strict prefix of `eos.pt_flash_saft`, which is a live pair here.

    A substring test reports the shorter id on a line naming the longer, which would send a
    reader to edit a skill that is already correct - and a gate that cries wolf is a gate that
    gets switched off.
    """
    found = problems("screening", "It calls `eos.pt_flash_saft`.\n", tmp_path)

    assert len(found) == 1, found
    assert "eos.pt_flash_saft" in found[0]
    assert "eos.pt_flash'" not in found[0], "the prefix was reported as well"


def test_the_prefix_of_a_package_path_is_still_matched(tmp_path: Path) -> None:
    """The form a skill writes is `azoth.hydraulics.darcy_weisbach`, id and all.

    A leading guard that rejected a preceding dot would miss every reference a skill actually
    makes, which is the shape this rule was first written with.
    """
    found = problems("screening", "It calls `azoth.hydraulics.darcy_weisbach`.\n", tmp_path)

    assert len(found) == 1, found
    assert "hydraulics.darcy_weisbach" in found[0]


def test_the_catalog_passes_its_own_gate() -> None:
    """The gate, run where it runs in CI."""
    import subprocess

    result = subprocess.run(
        [sys.executable, str(TOOLS / "validate_skills.py")],
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 0, f"{result.stdout}{result.stderr}"
    assert "OK" in result.stdout
