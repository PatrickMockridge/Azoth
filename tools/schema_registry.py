"""Every JSON Schema under `specs/schema/`, and a `referencing` registry over them.

The schemas cross-reference each other by `$id`: `keycard.schema.json` points at
`calc.schema.json` for what a citation and a unit are, `model.schema.json` points
at it for quantities and for what an identifier is, and `calc.schema.json` points
at `unit.schema.json` for the unit enum. A `$ref` across documents resolves only
if both documents are in the registry, so each tool that validates against one
schema needs all of them.

That list was hand-maintained, and adding a schema meant finding and editing every
tool that validates - which is the shape of mistake this repository keeps removing.
So the registry is built by reading the directory: a new schema resolves
everywhere the moment it exists, and no tool has to be told about it.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent.parent
SCHEMA_DIR = ROOT / "specs" / "schema"

#: The documents, keyed by their `$id`. Built by reading the directory rather than
#: by naming files, so a schema added to `specs/schema/` is resolvable everywhere
#: without an edit in each tool.
SCHEMAS: dict[str, dict[str, Any]] = {
    document["$id"]: document
    for document in (
        json.loads(path.read_text(encoding="utf-8")) for path in sorted(SCHEMA_DIR.glob("*.json"))
    )
}


def load(name: str) -> dict[str, Any]:
    """One schema by filename, e.g. `"calc.schema.json"`.

    Raises:
        KeyError: if no schema in `specs/schema/` has that `$id`. A tool asking for
            a schema that is not there is a programming error, not a user
            condition - and a lookup that returned `None` would show up as an
            obscure failure inside a validator instead.
    """
    wanted = f"/specs/schema/{name}"
    for schema_id, document in SCHEMAS.items():
        if schema_id.endswith(wanted):
            return document
    raise KeyError(f"no schema {name!r} in {SCHEMA_DIR}; found {sorted(SCHEMAS)}")


def registry() -> Any:
    """A `referencing` registry carrying every schema, so cross-document `$ref`s resolve."""
    try:
        from referencing import Registry, Resource
    except ImportError:  # pragma: no cover
        import sys

        sys.exit("schema_registry requires referencing: pip install referencing")

    reg = Registry()
    for schema_id, document in SCHEMAS.items():
        reg = reg.with_resource(schema_id, Resource.from_contents(document))
    return reg
