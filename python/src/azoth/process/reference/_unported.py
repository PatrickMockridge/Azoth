"""What this port does not carry, as a lookup rather than a sentence.

**A refusal is a fact about the declaration, so it is declared once and read here.** Every
input a spec declares carries an enumeration, and the members of that enumeration reach
the solve; the ones that do not are rows in the spec's ``[[unported]]`` array, which
:mod:`azoth._models_gen` carries whole. A call site names its row by key and gets back the
NeqSim symbol that would close the gap.

**The shape is what makes the two implementations comparable.** A refusal written as a
sentence is a sentence in Python and a different sentence in Rust, and the pair drifts -
measured, ``billet_schultes_1999`` was refused by this half and carried by the Rust half
while the spec said both were ported, and no test could see it. A call site is now one
string literal naming the row, in one language and the other alike, so
``tools/check_unported.py`` compares two sets of keys against the declaration instead of
reading two sets of prose.

A miss is an error as well, and it names the key that is missing: a refusal whose row was
never declared has to be a message somebody reads rather than a silent pass.
"""

from __future__ import annotations

from typing import Any

from azoth import _models_gen
from azoth.core.errors import InvalidInputError


def key_of(row: dict[str, Any]) -> str:
    """The canonical key a call site names a row by.

    ``parameter=value`` for a row that refuses one member of an enumeration, and
    ``parameter@other=member`` for a row that refuses a state - the conditions sorted and
    joined by ``&``, so two spellings of one condition are one key.
    """
    parameter = str(row["parameter"])
    if "value" in row:
        return f"{parameter}={row['value']}"
    conditions = sorted(f"{c['parameter']}={c['value']}" for c in row["when"])
    return f"{parameter}@{'&'.join(conditions)}"


def rows() -> dict[str, dict[str, Any]]:
    """Every declared row, by key."""
    return {key_of(row): row for model in _models_gen.MODELS for row in model.get("unported", ())}


def parameter_of(key: str) -> str:
    """The input a key names - the part before the ``=`` or the ``@``."""
    return key.replace("@", "=").split("=", 1)[0]


def refuse(key: str) -> InvalidInputError:
    """Refuse a value the declaration marks as not carried.

    ``key`` is exactly the row's own key: ``parameter=value``, or
    ``parameter@other=member`` for a row that depends on another input's value. It is a
    literal at every call site on purpose - a key the reader does not have to compute is
    a key the checker can read without following control flow, in this language and in
    the other.
    """
    parameter = parameter_of(key)
    row = rows().get(key)
    if row is None:
        return InvalidInputError(
            parameter,
            f"`{key}` is refused, and no `[[unported]]` row declares it. A refusal is a "
            f"claim about the declaration, so it belongs in the `[[unported]]` array of the "
            f"spec that carries this input, and `tools/check_unported.py` is what holds the "
            f"two together.",
        )
    close = f"`{row['class']}` is the class that would close it"
    measured = f" `{row['capture']}` measures what it does." if row.get("capture") else ""
    return InvalidInputError(parameter, f"`{key}` is not ported: {close}.{measured}")
