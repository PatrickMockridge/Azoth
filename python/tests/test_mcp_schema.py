"""The MCP pin gate: what it reads, and the ways its declaration can rot.

`tools/check_mcp_schema.py` runs in CI over a tree that is currently clean, which is the shape of a
check that can stop working without anyone noticing. So the readings are tested against a synthetic
tree rather than only against the repository.

**Three things matter and none of them is obvious.**

The pin is read out of `NOTICE` rather than restated in the tool, so a record whose lines are split
across a label and a continuation has to parse. A record the tool failed to parse is worse than a
missing one: the check would compare nothing, find nothing, and pass - which is what it did when the
first version of this parser required one token count for two different line shapes.

The exit codes are a contract, because the scheduled run is the only caller that can see the
difference: *these bytes moved* and *this check could not decide* have to stay apart, or a week of
unreachable upstream reads as a revision that moved.

And the declaration is checked in **both** directions: a file with no record, and a record with no
file. The second is the one a loop over the files never reaches.

The last test runs the real tool over the real tree. A check that is only tested against fakes is a
check whose own subject is unverified.
"""

from __future__ import annotations

import hashlib
import importlib
import sys
from pathlib import Path
from types import ModuleType

REPO_ROOT = Path(__file__).resolve().parents[2]

#: A record laid out the way `NOTICE` lays one out: the first line carries the label, the second is
#: a continuation. Both shapes have to parse, which is what the third test below is about.
NOTICE_HEAD = """  Project:   Model Context Protocol - https://example.test/spec
  Vendored:  python/tests/validation/vendor/mcp/, unmodified
  Source:    {first} https://example.test/{first}
             {second} https://example.test/{second}
  Digest:    {first} sha256 {first_digest}
             {second} sha256 {second_digest}
"""


def gate() -> ModuleType:
    """`tools/check_mcp_schema.py`, imported by name. The tools are not a package."""
    sys.path.insert(0, str(REPO_ROOT / "tools"))
    try:
        return importlib.import_module("check_mcp_schema")
    finally:
        sys.path.pop(0)


def rooted(tmp_path: Path, *files: str, notice: str) -> Path:
    """A synthetic root holding the vendored files named in `files`, each `path\\ntext`."""
    for entry in files:
        name, _, text = entry.partition("\n")
        path = tmp_path / "python/tests/validation/vendor/mcp" / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(text, encoding="utf-8")
    (tmp_path / "NOTICE").write_text(notice, encoding="utf-8")
    return tmp_path


def digest(text: str) -> str:
    """Computed here rather than by the tool, so a wrong digest is not one the tool agrees with."""
    return hashlib.sha256(text.encode()).hexdigest()


TWO = ("schema-2026-07-28.json", "schema-2025-11-25.json")


def notice_for(first: str, second: str, first_digest: str, second_digest: str) -> str:
    return NOTICE_HEAD.format(
        first=first, second=second, first_digest=first_digest, second_digest=second_digest
    )


def test_both_line_shapes_of_a_record_are_read() -> None:
    records = gate().pin(notice_for(TWO[0], TWO[1], digest("a"), digest("b")))
    assert records[TWO[0]] == {"digest": digest("a"), "url": f"https://example.test/{TWO[0]}"}
    assert records[TWO[1]] == {"digest": digest("b"), "url": f"https://example.test/{TWO[1]}"}


def test_a_file_that_matches_its_digest_is_not_drift(tmp_path: Path) -> None:
    tool = gate()
    body = '{"$defs": {}}'
    root = rooted(
        tmp_path,
        f"{TWO[0]}\n{body}",
        f"{TWO[1]}\n{body}",
        notice=notice_for(TWO[0], TWO[1], digest(body), digest(body)),
    )
    records = tool.pin((root / "NOTICE").read_text(encoding="utf-8"))
    assert tool.declaration(root, records) == []
    assert tool.drift(root, records) == []


def test_a_file_that_does_not_match_its_digest_is_drift(tmp_path: Path) -> None:
    tool = gate()
    root = rooted(
        tmp_path,
        f"{TWO[0]}\n" + '{"$defs": {}}',
        f"{TWO[1]}\n" + '{"$defs": {}}',
        # A digest of something else: the file is not what was vendored.
        notice=notice_for(TWO[0], TWO[1], digest("elsewhere"), digest('{"$defs": {}}')),
    )
    records = tool.pin((root / "NOTICE").read_text(encoding="utf-8"))
    found = tool.drift(root, records)
    assert len(found) == 1
    assert TWO[0] in found[0]


def test_a_vendored_file_with_no_record_is_a_broken_declaration(tmp_path: Path) -> None:
    tool = gate()
    # The second file is in the tree and NOTICE says nothing about it, so there is no digest to
    # hold it to: the check has learned nothing, which is not the same as having found drift.
    root = rooted(
        tmp_path,
        f"{TWO[0]}\n" + "body",
        f"{TWO[1]}\n" + "body",
        notice=notice_for(TWO[0], "schema-1970-01-01.json", digest("body"), digest("body")),
    )
    records = tool.pin((root / "NOTICE").read_text(encoding="utf-8"))
    found = tool.declaration(root, records)
    assert any(TWO[1] in problem for problem in found), found


def test_a_record_with_no_file_is_a_broken_declaration(tmp_path: Path) -> None:
    tool = gate()
    # A file deleted without its record: the loop over the files never reaches this one.
    root = rooted(
        tmp_path,
        f"{TWO[0]}\n" + "body",
        notice=notice_for(TWO[0], TWO[1], digest("body"), digest("body")),
    )
    records = tool.pin((root / "NOTICE").read_text(encoding="utf-8"))
    found = tool.declaration(root, records)
    assert any(TWO[1] in problem for problem in found), found


def test_a_record_with_no_url_is_a_broken_declaration(tmp_path: Path) -> None:
    tool = gate()
    # **The state this gate was written for.** A digest with no URL pins the bytes and nothing
    # about whether upstream still publishes them, so a new revision would be invisible.
    root = rooted(
        tmp_path,
        f"{TWO[0]}\n" + "body",
        notice=f"  Digest:    {TWO[0]} sha256 {digest('body')}\n",
    )
    records = tool.pin((root / "NOTICE").read_text(encoding="utf-8"))
    found = tool.declaration(root, records)
    assert any("no URL" in problem for problem in found), found


def test_a_missing_vendor_directory_decides_nothing(tmp_path: Path) -> None:
    tool = gate()
    root = rooted(tmp_path, notice=notice_for(TWO[0], TWO[1], digest("a"), digest("b")))
    records = tool.pin((root / "NOTICE").read_text(encoding="utf-8"))
    found = tool.declaration(root, records)
    assert any("is not there" in problem for problem in found), found


def test_an_empty_vendor_directory_decides_nothing(tmp_path: Path) -> None:
    tool = gate()
    root = rooted(tmp_path, notice=notice_for(TWO[0], TWO[1], digest("a"), digest("b")))
    (root / "python/tests/validation/vendor/mcp").mkdir(parents=True)
    records = tool.pin((root / "NOTICE").read_text(encoding="utf-8"))
    found = tool.declaration(root, records)
    # Both directions fire as well, because neither record has a file to claim it.
    assert any("decided nothing" in problem for problem in found), found


def test_the_real_tree_is_pinned_and_unchanged() -> None:
    tool = gate()
    records = tool.pin((REPO_ROOT / "NOTICE").read_text(encoding="utf-8"))
    assert tool.declaration(REPO_ROOT, records) == []
    assert tool.drift(REPO_ROOT, records) == []
    # Two lanes, and neither is a record this test could pass by not reading.
    assert sorted(records) == sorted(TWO)
