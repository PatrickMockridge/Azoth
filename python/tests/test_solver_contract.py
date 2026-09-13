"""The solver vocabulary is one set, agreed on by three separate places.

A calc spec that is implicit in its unknown names the *scheme* to solve it with, in
its ``solver.kind``. That name is a cross-language contract, in exactly the way the
warning vocabulary and the unit vocabulary are, and it has to hold across:

* ``specs/schema/calc.schema.json``'s ``solver.kind`` enum, which is what a spec is
  validated against;
* :class:`azoth.core.solver.SolverKind`, which the Python reference dispatches on;
* ``crates/azoth-core/src/solver.rs``'s ``SolverKind::ALL``, the Rust side, reachable
  from here through ``azoth._core.solver_kinds()``.

Why this file had to exist before a second kind could be added: the schema's own
description of ``solver.kind`` said so. It reads *"Nothing does that for solver kinds
today, which is why this enum is narrow rather than merely unchecked"* - the enum was
kept to one value not because a second was unwanted but because nothing would have
caught the three lists drifting apart if one were added. That is the same shape as the
``K`` defect the unit contract was written for: a vocabulary with no comparator is not
a contract, it is three opinions.

Whether the schemes *run* identically is a different question with a different test -
``test_cross_impl.py`` compares the numbers and the iteration counts. This file is
about the vocabulary being one set, and about an unknown kind failing loudly rather
than falling back to the only one anybody implemented.
"""

from __future__ import annotations

import importlib
import json
from pathlib import Path
from types import ModuleType
from typing import Any

import pytest

from azoth._registry_gen import CALCS
from azoth.core.errors import InvalidInputError
from azoth.core.solver import SolverKind

REPO_ROOT = Path(__file__).resolve().parents[2]
SCHEMA_PATH = REPO_ROOT / "specs" / "schema" / "calc.schema.json"


def _extension() -> ModuleType:
    """The compiled extension, imported by name.

    The same accessor ``test_units_contract.py`` and ``test_cross_impl.py`` use, and
    for the same reason: the module does not exist until the bindings are built, and
    an attribute mypy cannot resolve is a worse trade than a lookup that fails
    clearly.
    """
    try:
        return importlib.import_module("azoth._core")
    except ImportError as exc:  # pragma: no cover - the marker skips these
        raise AssertionError("azoth._core is not built; run `maturin develop`") from exc


def schema_solver_kinds() -> set[str]:
    """The solver kind strings the spec schema permits."""
    schema: dict[str, Any] = json.loads(SCHEMA_PATH.read_text(encoding="utf-8"))
    return set(schema["properties"]["solver"]["properties"]["kind"]["enum"])


def python_solver_kinds() -> set[str]:
    """The kinds the Python reference can dispatch on."""
    return {kind.value for kind in SolverKind}


def test_the_schema_and_the_python_vocabulary_agree() -> None:
    """The schema's enum and ``SolverKind`` are the same set.

    A name the schema allows but this class lacks would fail at runtime in
    ``SolverKind.parse``, and only for whichever calc happened to declare it. This
    is the cheap version of that check, and it runs without the extension built.
    """
    schema, python_side = schema_solver_kinds(), python_solver_kinds()
    assert schema == python_side, (
        f"solver kind vocabularies differ\n"
        f"  only in the schema: {sorted(schema - python_side)}\n"
        f"  only in SolverKind: {sorted(python_side - schema)}"
    )


@pytest.mark.requires_rust
def test_the_rust_vocabulary_agrees_with_both() -> None:
    """All three lists are one set.

    The Rust leg is the one that was unchecked, and it is the leg that decides what
    a spec's ``solver.kind`` actually *runs*: a kind present in the schema and in
    Python but absent from ``SolverKind::ALL`` would lint clean and then fail in the
    Rust implementation alone. ``solver_kinds()`` exists so this can be asserted
    from Python rather than by reading Rust source.
    """
    rust_side = set(_extension().solver_kinds())
    assert schema_solver_kinds() == rust_side, (
        f"the Rust solver vocabulary differs from the schema\n"
        f"  only in the schema: {sorted(schema_solver_kinds() - rust_side)}\n"
        f"  only in Rust: {sorted(rust_side - schema_solver_kinds())}"
    )
    assert rust_side == python_solver_kinds()


def test_every_solver_kind_a_spec_declares_is_in_the_vocabulary() -> None:
    """No spec may name a solver kind the vocabulary does not carry.

    The schema already enforces this, so this test is a guard on the guard: it
    proves the schema is actually being applied to the specs in the registry, and
    that the generated registry is current.
    """
    assert CALCS, "no calcs in the registry - this test would pass vacuously"

    declared = {calc["solver"]["kind"] for calc in CALCS if calc.get("solver") is not None}
    unknown = declared - schema_solver_kinds()
    assert not unknown, (
        f"spec(s) declare solver kind(s) {sorted(unknown)} which the vocabulary lacks"
    )


def test_at_least_one_spec_declares_a_solver() -> None:
    """The check above must not be vacuous.

    ``every_solver_kind_a_spec_declares_is_in_the_vocabulary`` iterates the specs
    that declare one, so it passes trivially if none do. Colebrook is implicit and
    must declare ``fixed_point``; asserting that here means the vocabulary test
    cannot quietly stop testing anything.
    """
    declaring = [calc["id"] for calc in CALCS if calc.get("solver") is not None]
    assert declaring, (
        "no spec declares a solver, so the vocabulary contract is untested - "
        "hydraulics.friction_factor_colebrook is implicit and should be listed"
    )


def test_every_kind_parses_and_round_trips_through_its_name() -> None:
    """Each kind survives ``parse(kind.value)``, and an unknown name is an error.

    Unknown names must raise rather than default. Falling back to the one
    implemented kind would run a scheme the spec did not ask for, which is a wrong
    number that looks reasonable - the failure this library is organised against.
    """
    for kind in SolverKind:
        assert SolverKind.parse(kind.value) is kind

    with pytest.raises(InvalidInputError, match="bisection"):
        SolverKind.parse("bisection")
