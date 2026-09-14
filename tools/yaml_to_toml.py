#!/usr/bin/env python3
"""Convert a YAML document to TOML, preserving every value exactly.

Written for one migration - this repository's data files moving from YAML to TOML - and
the migration is the only reason it exists. What it is *for* is the check beside it in
`python/tests/test_yaml_to_toml.py`: every YAML file in the repository is converted, the
result is parsed back with `tomllib`, and the two trees are compared. Two parsers, one
document, real data - the arrangement this repository uses wherever a contract has two
implementations.

# What this has to get right, and why a stock writer does not

**The specs are prose.** A `description:` or a bound's `rationale:` is a paragraph, and
YAML's folded block is what keeps it readable. `tomli-w` "avoids multi-line strings by
default", so using it would emit every paragraph as one enormous line - correct, and
unreviewable, which for a file whose purpose is to be reviewed is a real loss.

TOML's multi-line basic string is exact rather than a compromise. A trailing unescaped
backslash trims the newline and the next line's indentation, so a wrapped line rejoins
with the space YAML folded it to; a newline that is *in* the value is written as the
two-character escape followed by a backslash; and the closing delimiter follows the last
character, because TOML trims a newline immediately before it and YAML's clip mode kept
one. So a folded block and its conversion parse to the same string, trailing newline and
all - see `python/tests/test_yaml_to_toml.py`, which is where that claim is checked over
every file in the repository rather than asserted here.

**A table cannot be re-opened.** TOML forbids adding a key to a table once a sub-table of
it has been opened, so the emitter puts every scalar of a table before any of its
children. That means output key order is not input key order - a deliberate consequence,
not an oversight.

**A continuation eats the indentation it lands on.** So a wrapped line can carry the
prose it was wrapped from, and cannot carry a value whose own lines are indented: that
indentation has to be written as escapes to stop being whitespace. It is a short rule and
it is the difference between a value that survives and a value that quietly loses its
layout, which is why both spellings are exercised by name in the test rather than left to
the corpus to notice.

**A record with a newline cannot be an inline table**, because TOML's inline tables are
single-line. So a list of records is emitted as `[[path]]` blocks when any record needs a
multi-line string, and as a compact inline array otherwise - which is what keeps
`databank/manifest.yaml`'s 1,496 one-line column rows from becoming 9,000 lines.

Usage:
    python tools/yaml_to_toml.py specs/calcs/eos/pr_kappa.yaml     # print the TOML
    python tools/yaml_to_toml.py --check specs/calcs/eos/*.yaml   # round-trip, write nothing
"""

from __future__ import annotations

import argparse
import datetime
import json
import re
import sys
from pathlib import Path
from typing import Any

try:
    import yaml
except ImportError:  # pragma: no cover
    sys.exit("yaml_to_toml requires PyYAML: pip install pyyaml")

#: Where prose wraps. The width this repository's own prose is written to, so a
#: converted file reads at the same measure as the one it replaces.
WIDTH = 88

#: The indentation of a continued line, inside a multi-line string.
INDENT = "    "

_BARE_KEY = re.compile(r"[A-Za-z0-9_-]+")

#: What a basic string has to escape: a backslash, a quote, and the control
#: characters TOML has no literal form for. A tab is the one that occurs.
_BASIC_UNSAFE = re.compile(r'[\\"\x00-\x1f]')


def toml_key(name: str) -> str:
    """A TOML key: bare where it can be, quoted where it cannot.

    `hydraulics.orifice_flow` is one key and not a path, which is the case that makes
    this necessary rather than decorative - an unquoted dot would silently become two
    levels of table.
    """
    if _BARE_KEY.fullmatch(name):
        return name
    return json.dumps(name, ensure_ascii=False)


def toml_path(parts: tuple[str, ...]) -> str:
    """A dotted key path, each segment quoted as it needs."""
    return ".".join(toml_key(part) for part in parts)


def _wrap(text: str, width: int) -> list[str]:
    """Break `text` at spaces so no piece exceeds `width`.

    Splitting on a single space and rejoining with one is lossless - a run of two
    spaces becomes an empty piece that rejoins to two - which matters because the
    continuation this feeds joins with a space of its own.
    """
    if len(text) <= width:
        return [text]
    lines: list[str] = []
    current = ""
    for word in text.split(" "):
        if current and len(current) + 1 + len(word) > width:
            lines.append(current)
            current = word
        elif current:
            current = f"{current} {word}"
        else:
            current = word
    lines.append(current)
    return lines


def _indented(value: str) -> bool:
    """Whether any line after the first begins with whitespace.

    The question is whether the layout carries meaning, which is what indentation
    inside a worked substitution or a block of pseudo-code does and what wrapped prose
    does not. The first line is excluded because a line begins after a delimiter rather
    than after a newline, so its own indentation is not layout - and because this
    decides readability only: `_escaped_layout` below makes the wrapped spelling exact
    whatever the lines look like, so this is free to be a judgement.
    """
    return any(line[:1] in (" ", "\t") for line in value.split("\n")[1:])


def _escaped_quotes(value: str) -> str:
    """A value as the body of a multi-line basic string.

    A backslash and a run of three quotes are the two things TOML would read as
    structure, so both are escaped; everything else travels untouched. One or two
    adjacent quotes are legal in a multi-line basic string and stay as themselves,
    which is what keeps ordinary prose's quotation marks readable.
    """
    return value.replace("\\", "\\\\").replace('"""', '\\"\\"\\"')


def _escaped_every_quote(value: str) -> str:
    """A value as the body of a one-line basic string.

    Every quote is escaped here, because the delimiter is one quote wide: a run of two
    is harmless in the multi-line spelling above and is a second delimiter in this one.
    """
    return value.replace("\\", "\\\\").replace('"', '\\"')


def _escaped_layout(line: str) -> str:
    """A line's leading whitespace, escaped so that TOML cannot eat it.

    A backslash at the end of a line "will be trimmed along with all whitespace up to
    the next non-whitespace character", so a wrapped spelling of a value whose second
    line is indented would lose that indentation - silently, and only for the values
    somebody meant to lay out. Writing the indentation as `\\t` and `\\u0020` stops it
    being whitespace, which is the only way to keep it in a basic string.
    """
    stripped = line.lstrip(" \t")
    prefix = line[: len(line) - len(stripped)]
    if not prefix:
        return line
    return "".join("\\t" if character == "\t" else "\\u0020" for character in prefix) + (stripped)


def toml_string(value: str) -> str:
    """A TOML string, in whichever spelling needs the least escaping.

    **Verbatim where the layout means something, wrapped where it does not.** Both
    halves of that rule have a reason:

    * A value whose later lines are indented - a derivation, a block of pseudo-code - is
      a **literal multi-line string**, written exactly as it stands. TOML trims one
      newline after the opening delimiter and nothing else, so this reproduces the value
      byte for byte and needs no escaping at all: backslashes and indentation both
      survive as themselves. Readability is the reason to prefer it, not correctness.
    * A value that is one long line, or several unindented ones, is a **basic multi-line
      string**, wrapped at spaces with a trailing backslash. That continuation trims the
      next line's leading whitespace, so indentation that is part of the value is written
      as escapes - `_escaped_layout` - which makes this spelling exact for any value.
    * A short single-line value is a **basic string** - the ordinary `"..."` - unless it
      holds something a basic string would have to escape. Which it usually does not, and
      when it does it is the one case a literal `'''` string is worth the unusual
      spelling for: a spec's `latex` is LaTeX and its `references` are citations with
      quoted titles, and doubling every backslash in those is how a typo gets in.

    A `'''` in the value, or a value ending in a quote, falls back to the escaped basic
    form, because either would collide with the literal delimiter.
    """
    if "\n" not in value:
        literal = (
            _BASIC_UNSAFE.search(value) is not None
            and len(value) <= WIDTH
            and "'''" not in value
            and not value.endswith("'")
        )
        if literal:
            return f"'''{value}'''"
        return _basic_single(value)

    if "'''" not in value and not value.endswith("'") and _indented(value):
        return "'''\n" + value + "'''"

    escaped = _escaped_quotes(value)
    pieces = [
        f" \\\n{INDENT}".join(_wrap(_escaped_layout(segment), WIDTH))
        for segment in escaped.split("\n")
    ]
    return '"""\\\n' + INDENT + f"\\n\\\n{INDENT}".join(pieces) + '"""'


def _basic_single(value: str) -> str:
    """A one-line basic string, escaped and wrapped if it must be."""
    escaped = _escaped_every_quote(value)
    lines = _wrap(escaped, WIDTH)
    if len(lines) == 1:
        return f'"{lines[0]}"'
    # Wrapped, so this is a multi-line basic string after all, and its first piece sits
    # after a continuation backslash: leading whitespace has to be escaped to survive.
    lines = _wrap(_escaped_layout(escaped), WIDTH)
    return '"""\\\n' + INDENT + f" \\\n{INDENT}".join(lines) + '"""'


def _is_record_list(value: Any) -> bool:
    """Whether a list holds records, which are emitted as tables rather than inline."""
    return isinstance(value, list) and any(isinstance(item, dict) for item in value)


def _is_leaf(value: Any) -> bool:
    """Whether a value is written as `key = value` inside the table it belongs to.

    A record list is a leaf when it fits on one line as an inline array, and **not** a
    leaf when it needs `[[blocks]]`. That difference is not cosmetic: a `[[header]]`
    opens a scope and a bare key does not, so a list emitted inline has to go with the
    other keys of its table and cannot be written after a sub-table has been opened.
    Getting this wrong puts the key inside whichever table was opened last, which parses
    without complaint and loses it from its own.
    """
    if isinstance(value, dict):
        return False
    if _is_record_list(value):
        return not _needs_block(value)
    return True


def _wrapped(value: Any) -> bool:
    """Whether a value holds a string too long, or too multi-line, for one line.

    Recursive because a record's field can be a nested table - a case's `inputs` is one
    - and an inline table is single-line all the way down.
    """
    if isinstance(value, str):
        return len(value) > WIDTH or "\n" in value
    if isinstance(value, dict):
        return any(_wrapped(item) for item in value.values())
    if isinstance(value, list):
        return any(_wrapped(item) for item in value)
    return False


def _needs_block(records: list[dict[str, Any]]) -> bool:
    """Whether a record list must be `[[blocks]]` rather than one line of inline tables.

    An inline table is single-line by definition, so a record holding a string that
    wrapped cannot be one. This is the rule that keeps the manifest compact without
    ever producing output TOML forbids.
    """
    return any(_wrapped(record) for record in records)


def toml_value(value: Any) -> str:
    """A scalar or a list of scalars, as TOML."""
    if isinstance(value, bool):
        # Before `int`: a bool *is* an int in Python, and `true` is not `1`.
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
        # `repr` round-trips a float exactly, which is the whole requirement.
        return repr(value)
    if isinstance(value, str):
        return toml_string(value)
    if isinstance(value, (datetime.datetime, datetime.date)):
        # Bare, so TOML reads it back as the same date or datetime PyYAML produced.
        return value.isoformat()
    if value is None:
        raise ValueError(
            "TOML has no null. A field that is absent must be omitted from the "
            "document rather than written as a null."
        )
    if isinstance(value, list):
        return "[" + ", ".join(toml_value(item) for item in value) + "]"
    raise ValueError(f"cannot write {type(value).__name__} as TOML")


def _emit_scope(
    table: dict[str, Any], at: tuple[str, ...], out: list[str], *, declared: bool
) -> None:
    """One table and everything under it, scalars first so no parent is re-opened.

    `declared` says whether this scope needs a `[header]` of its own. A record inside
    `[[blocks]]` does not - the array-table line above it is the header - but a table
    nested in one does, and it is `[blocks.nested]`, relative to the record it sits in.

    **A child here is always a table or an array of them.** A list of scalars is a leaf,
    and so is a record list that fits on one line, so neither reaches the loop below.
    """
    if declared:
        out.append(f"[{toml_path(at)}]")

    leaves = {name: value for name, value in table.items() if _is_leaf(value)}
    children = {name: value for name, value in table.items() if not _is_leaf(value)}
    for name, value in leaves.items():
        # `_inline_value`, not `toml_value`: a record list that fits on one line is a
        # leaf, and the scalar writer cannot write a table.
        out.append(f"{toml_key(name)} = {_inline_value(value)}")

    for name, value in children.items():
        child = (*at, name)
        if isinstance(value, dict):
            _emit_scope(value, child, out, declared=True)
        else:
            for record in value:
                out.append(f"[[{toml_path(child)}]]")
                # A record's own fields, then its sub-tables - `[tests.inputs]` opens
                # inside the `[[tests]]` above it, which is how TOML nests them.
                _emit_scope(record, child, out, declared=False)


def _inline_value(value: Any) -> str:
    """A field of an inline table: a nested inline table, or a plain value.

    An inline table may hold another one, which a case's `inputs` is - so this recurses
    rather than assuming a record's fields are scalars.
    """
    if isinstance(value, dict):
        return _inline_record(value)
    if _is_record_list(value):
        return _inline(value)
    return toml_value(value)


def _inline_record(record: dict[str, Any]) -> str:
    """One record as an inline table."""
    fields = ", ".join(f"{toml_key(f)} = {_inline_value(v)}" for f, v in record.items())
    return "{" + fields + "}"


def _inline(records: list[dict[str, Any]]) -> str:
    """A list of records as one line of inline tables, which may nest."""
    return "[" + ", ".join(_inline_record(record) for record in records) + "]"


def convert(document: Any) -> str:
    """A YAML document, as TOML.

    Raises:
        ValueError: if the document is not a mapping - a TOML file is a table - or
            holds a null, which TOML has no way to write.
    """
    if not isinstance(document, dict):
        raise ValueError(f"a TOML document is a table, and this is a {type(document).__name__}")
    out: list[str] = []
    _emit_scope(document, (), out, declared=False)
    return "\n".join(out) + "\n"


def convert_text(text: str) -> str:
    """A YAML document as TOML, from its source text."""
    return convert(yaml.safe_load(text))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("paths", nargs="+", type=Path)
    parser.add_argument(
        "--check",
        action="store_true",
        help=(
            "convert and parse back, and fail if a value moved. Writes nothing. "
            "This is the same comparison `python/tests/test_yaml_to_toml.py` makes "
            "over the whole repository; it is here to point at one file."
        ),
    )
    args = parser.parse_args()

    if not args.check:
        for path in args.paths:
            print(convert_text(path.read_text(encoding="utf-8")), end="")
        return 0

    import tomllib

    failed = False
    for path in args.paths:
        original = yaml.safe_load(path.read_text(encoding="utf-8"))
        try:
            converted = convert(original)
        except ValueError as exc:
            print(f"yaml_to_toml: {path}: {exc}", file=sys.stderr)
            failed = True
            continue
        back = tomllib.loads(converted)
        if back != original:
            print(f"yaml_to_toml: {path} does not survive the conversion", file=sys.stderr)
            failed = True
    if failed:
        return 1
    print(f"yaml_to_toml: {len(args.paths)} file(s) convert losslessly")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
