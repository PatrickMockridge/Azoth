"""Write a document as TOML, for the tests that need a synthetic spec.

`tomllib` reads and nothing in this repository writes, so a test that mutates a parsed
spec has no way to put it back. This is the smallest writer that covers what a spec
holds: scalars, arrays of scalars, tables and arrays of tables. Everything a spec's
prose needs is an ordinary quoted string.

TOML forbids adding a key to a table once a sub-table of it has been opened, so a
table's own values are written before any of its children.
"""

from __future__ import annotations

import re
from typing import Any

_BARE = re.compile(r"[A-Za-z0-9_-]+")


def _key(name: str) -> str:
    """A key, quoted where it must be. A dotted key is one key and not a path."""
    return name if _BARE.fullmatch(name) else _string(name)


def _string(value: str) -> str:
    escaped = value.replace("\\", "\\\\").replace('"', '\\"')
    escaped = escaped.replace("\n", "\\n").replace("\t", "\\t").replace("\r", "\\r")
    return f'"{escaped}"'


def _scalar(value: Any) -> str:
    """One TOML value that is not a table or a list of them."""
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, int):
        return str(value)
    if isinstance(value, float):
        if value != value:
            return "nan"
        if value == float("inf"):
            return "inf"
        if value == float("-inf"):
            return "-inf"
        return repr(float(value))
    if isinstance(value, str):
        return _string(value)
    if isinstance(value, list):
        return "[" + ", ".join(_scalar(item) for item in value) + "]"
    raise TypeError(f"a fixture has no {type(value).__name__}")


def _is_table(value: Any) -> bool:
    return isinstance(value, dict)


def _is_table_list(value: Any) -> bool:
    return isinstance(value, list) and bool(value) and all(isinstance(i, dict) for i in value)


def _emit(table: dict[str, Any], at: tuple[str, ...], out: list[str], *, header: bool) -> None:
    """One table: its own values, then each child table under its own header."""
    if header:
        out.append(f"[{'.'.join(_key(part) for part in at)}]")

    for name, value in table.items():
        if not _is_table(value) and not _is_table_list(value):
            out.append(f"{_key(name)} = {_scalar(value)}")

    for name, value in table.items():
        child = (*at, name)
        if _is_table(value):
            _emit(value, child, out, header=True)
        elif _is_table_list(value):
            for record in value:
                out.append(f"[[{'.'.join(_key(part) for part in child)}]]")
                _emit(record, child, out, header=False)


def dump(document: dict[str, Any]) -> str:
    """A document as TOML text, which `tomllib.loads` reads back unchanged."""
    out: list[str] = []
    _emit(document, (), out, header=False)
    return "\n".join(out) + "\n"
