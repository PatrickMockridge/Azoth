"""What this namespace does not carry, as a lookup rather than a sentence.

The same shape as :mod:`azoth.process.reference._unported`, and for the same reason: a refusal
written as a sentence is a sentence in Python and a different one in Rust, and the pair drifts
unobserved. A call site names its row by key and ``tools/check_unported.py`` compares the two
key sets against the declaration.

**Local rather than shared**, because the process namespace's copy is reachable only through
:mod:`azoth.process`, whose ``__init__`` imports the compiled extension - and this namespace's
reference implementation must import with or without it.
"""

from __future__ import annotations

from typing import Any

from azoth import _models_gen
from azoth.core.errors import InvalidInputError


def key_of(row: dict[str, Any]) -> str:
    """The canonical key a call site names a row by: ``parameter=value``."""
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

    ``key`` is exactly the row's own key, and it is a literal at every call site on purpose - a
    key the reader does not have to compute is a key the checker can read without following
    control flow, in this language and in the other.
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
