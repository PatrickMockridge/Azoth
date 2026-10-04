"""The Lean development's gate names theorems that exist, and a theorem exists.

`tools/check_lean_axioms.py` proves something about every theorem `Azoth/Axioms.lean`
names: that its proof has no gap and rests on nothing this project did not agree to.
Two things it cannot check, because both are about the *file* rather than about what
Lean reports, are here.

**A name that is not a theorem.** `#print axioms Foo.bar` for an unknown `Foo.bar`
is an error and the gate fails - but only because Lean refuses to elaborate it, and
the message says "unknown constant" rather than which claim went missing, and a
misspelled name fails the same way. This asserts the two sides agree, so the reason
is stated rather than inferred from a Lean error.

**A gate that gates nothing.** `check_lean_axioms.py` refuses an empty
`Axioms.lean` for the same reason, in the language that runs it; asserting it here
as well means a change that emptied the file fails in the suite a contributor runs,
not only in CI.

What neither checks, and what nothing can, is whether a theorem is the one
`docs/src/calculus/` states - a proof of a weakened hypothesis has no gap and is
still not the claim. The pages carry the statement in prose and each names its
theorem, so the two are held together by a reader rather than by a machine, and
saying so here is better than a check that appears to do it.
"""

from __future__ import annotations

import re
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
LEAN_DIR = REPO_ROOT / "lean"

#: The files carrying `#print axioms` lines: hand-written general theorems, one generated
#: line per canonical unit, and one per guard theorem. See `tools/check_lean_axioms.py`.
GATES = (
    LEAN_DIR / "Azoth" / "Axioms.lean",
    LEAN_DIR / "Azoth" / "Gate.lean",
    LEAN_DIR / "Azoth" / "GuardGate.lean",
)
LEAN_SOURCE = LEAN_DIR / "Azoth"

#: `#print axioms Azoth.Dim.ofExponents_nil`
_PRINTED = re.compile(r"^\s*#print axioms\s+(?P<name>[A-Za-z0-9_.]+)\s*$", re.MULTILINE)

#: `theorem ofExponents_nil`, `def exponents`, and the same inside a namespace.
#: A declaration's name, **dot included**: `def Shift.toSi` declares `Shift.toSi`, which is
#: the namespace a structure's own name opens, and a pattern without the dot reads it as
#: `Shift` and then reports the gate as naming a theorem that does not exist.
_DECLARED = re.compile(
    r"^\s*(?:theorem|lemma|def|abbrev)\s+(?P<name>[A-Za-z0-9_'À-ÿ.]+)",
    re.MULTILINE,
)


def gated_names() -> list[str]:
    """Every theorem the gate files ask Lean to report the axioms of."""
    names: list[str] = []
    for gate in GATES:
        names.extend(_PRINTED.findall(gate.read_text(encoding="utf-8")))
    return names


def declared_names() -> set[str]:
    """Every theorem or definition declared under `Azoth/`, fully qualified.

    A stack rather than a variable, because `Dim.lean` nests `Azoth` and `Dim` -
    and a single variable would attribute `ofExponents_nil` to `Azoth` and then
    report the gate as naming a theorem that does not exist.
    """
    names: set[str] = set()
    for path in sorted(LEAN_SOURCE.glob("*.lean")):
        open_namespaces: list[str] = []
        for line in path.read_text(encoding="utf-8").splitlines():
            stripped = line.strip()
            if stripped.startswith("namespace "):
                open_namespaces.append(stripped.removeprefix("namespace ").strip())
                continue
            if stripped.startswith("end ") and open_namespaces:
                open_namespaces.pop()
                continue
            match = _DECLARED.match(line)
            if match:
                leaf = match.group("name")
                names.add(".".join([*open_namespaces, leaf]))
    return names


def test_the_gate_names_something() -> None:
    """An empty gate is not a gate, and must fail here as well as in CI."""
    name = "Azoth.Dim.ofExponents_nil"
    gates = ", ".join(g.relative_to(REPO_ROOT).as_posix() for g in GATES)
    assert gated_names(), (
        f"the gate files ({gates}) name no theorem, so nothing is gated. A Lean file "
        f"that prints nothing passes every axiom check by having nothing to check - "
        f"the vacuous-gate failure the gate exists to catch."
    )
    assert name in gated_names(), (
        f"the gate does not name {name}, which is one of the theorems this "
        f"development proves. Add its `#print axioms` line in the same commit as the "
        f"theorem."
    )


def test_every_gated_name_is_a_declaration_that_exists() -> None:
    """A name in the gate that no file declares is a claim about nothing.

    Lean reports this as `unknown constant`, which fails the build but does not say
    which claim went missing - so the two sides are compared here, where the message
    can.
    """
    declared = declared_names()
    missing = [name for name in gated_names() if name not in declared]
    assert not missing, (
        f"the gate files name {missing}, which no file under "
        f"lean/Azoth/ declares. Either the name is misspelled or the theorem was "
        f"removed; in both cases the gate is naming something that is not there."
    )


# --- the input table ---------------------------------------------------------


def test_the_input_table_has_one_row_per_generated_dimension() -> None:
    """**The theorem is proved over rows, so a generator that emitted none would be vacuous.**

    `Azoth.Inputs.rows` is one entry per (palette entry, input) whose declaration carries a
    unit, and the Rust table beside it is the same set. The two are compared here by count
    because the Lean side cannot read the Rust one: `every_input_dimension_is_its_units`
    closes over whatever rows exist, and a table that collapsed to nothing would close just as
    quietly as one that is right.
    """
    lean = (REPO_ROOT / "lean" / "Azoth" / "Inputs.lean").read_text(encoding="utf-8")
    rust = (REPO_ROOT / "crates" / "azoth-process" / "src" / "model_inputs_gen.rs").read_text(
        encoding="utf-8"
    )

    lean_rows = len(re.findall(r'^    \("unit_ops\.', lean, flags=re.MULTILINE))
    rust_rows = len(re.findall(r"dimension: Some\(", rust))
    assert lean_rows > 0, "the Lean input table is empty"
    assert lean_rows == rust_rows, (
        f"the Lean table has {lean_rows} dimensioned input(s) and the Rust one {rust_rows}; "
        f"both are emitted by `tools/gen_model_inputs.py`, so one of them is stale"
    )
    assert "def rowCount : Nat := rows.length" in lean
    assert "theorem every_input_dimension_is_its_units" in lean
    # And the gate covers it, which is what makes the theorem checked rather than written.
    gate = (REPO_ROOT / "lean" / "Azoth" / "Gate.lean").read_text(encoding="utf-8")
    assert "#print axioms Azoth.Inputs.every_input_dimension_is_its_units" in gate
