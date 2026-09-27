"""The front-end lint: the three rules, and the exemption that keeps them honest.

`tools/check_ui.py` runs in CI and over a tree that is currently clean — zero findings on the
first run, once the eleven it found were fixed rather than muted. That is the shape of a check
that can stop working without anyone noticing, so each rule is tested against synthetic
TypeScript rather than only against the tree.

**The rules were verified by breaking the thing they guard**, which is the tree's own standard:
each case below is a file that trips a rule, and the tree is the evidence that the rule reads a
real shape rather than an imagined one — the cast rule's first run found eleven sites in
`ui/src`, none of which had been looked at.
"""

from __future__ import annotations

import importlib
import sys
from pathlib import Path
from types import ModuleType
from typing import Any

REPO_ROOT = Path(__file__).resolve().parents[2]


def lint() -> ModuleType:
    """`tools/check_ui.py`, imported by name. The tools are not a package."""
    sys.path.insert(0, str(REPO_ROOT / "tools"))
    try:
        return importlib.import_module("check_ui")
    finally:
        sys.path.pop(0)


def offenders(tmp_path: Path, name: str, text: str) -> list[str]:
    """The lint's verdict on one synthetic file, with the root redirected for the message.

    A synthetic file under `/tmp` against the real root would raise from `relative_to` rather
    than report. `Any` where the tool is reached: a `ModuleType` declares neither `ROOT` nor
    `offenders`, and the tools are outside mypy's `files`.
    """
    path = tmp_path / "ui" / "src" / name
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")
    tool: Any = lint()
    original = tool.ROOT
    tool.ROOT = tmp_path
    try:
        found: list[str] = tool.offenders(path)
        return found
    finally:
        tool.ROOT = original


def test_a_cast_on_the_boundary_is_reported(tmp_path: Path) -> None:
    """The defect the decoder was written for: an assertion where a parse belongs."""
    messages = offenders(tmp_path, "wire/client.ts", "const e = JSON.parse(json) as Envelope;\n")
    assert len(messages) == 1, messages
    assert "asserts a shape" in messages[0]

    # The double cast, which is the same assertion made twice.
    assert len(offenders(tmp_path, "wire/http.ts", "const c = x as unknown as Catalogue;\n")) == 1

    # A decoder's own cast is what the file marker is for.
    marked = (
        "// ui-ok-file: this module is the decoder\nconst r = value as Record<string, unknown>;\n"
    )
    assert offenders(tmp_path, "wire/decode.ts", marked) == []


def test_a_rename_is_not_a_cast(tmp_path: Path) -> None:
    """`export type { EditorCommand as Command }` renames a binding, and neither asserts."""
    text = 'export type { EditorCommand as Command } from "./commands";\n'
    assert offenders(tmp_path, "wire/types.ts", text) == []
    text = 'import { Handle, Position, type NodeProps } from "@xyflow/react";\n'
    assert offenders(tmp_path, "components/UnitOpNode.tsx", text) == []


def test_a_measurement_written_in_the_front_end_is_reported(tmp_path: Path) -> None:
    """`6894.757` is a unit's worth written in the fourth place the vocabulary page allows none."""
    messages = offenders(tmp_path, "state/units.ts", "const psi = 6894.757;\n")
    assert len(messages) == 1, messages
    assert "7 significant digits" in messages[0]

    # **Significant digits, not decimal places.** `273.15` is five of the former and two of the
    # latter, and the first form of this pattern missed it for exactly that reason.
    assert len(offenders(tmp_path, "wire/field.ts", "const offset = 273.15;\n")) == 1
    assert len(offenders(tmp_path, "wire/field.ts", "const inch = 0.0254;\n")) == 1
    assert len(offenders(tmp_path, "state/units.ts", "const seconds = 3600;\n")) == 1


def test_structure_is_not_a_measurement(tmp_path: Path) -> None:
    """`0`, `1`, `4` and `0.5` are a count, an identity, a width — not a unit's worth."""
    for text in (
        "const same = { unit, factor: 1, offset: 0 };\n",
        'const first = sets[0]?.id ?? "";\n',
        "const digits = 4;\n",
        "const half = 0.5;\n",
        "return { unit: wanted, factor, offset: target?.offset ?? 0 };\n",
    ):
        assert offenders(tmp_path, "state/units.ts", text) == [], text


def test_a_part_of_the_document_in_local_state_is_reported(tmp_path: Path) -> None:
    """A graph in state is the second copy the editor's claim says does not exist."""
    messages = offenders(tmp_path, "App.tsx", "const [graph, setGraph] = useState<Graph>(...);\n")
    assert len(messages) == 1, messages
    assert "already carries it" in messages[0]
    rows = "const [rows, setRows] = useState<GraphNode[]>([]);\n"
    assert len(offenders(tmp_path, "App.tsx", rows)) == 1


def test_the_envelope_itself_in_state_is_the_design(tmp_path: Path) -> None:
    """**A rule that flagged this would flag the editor.** The whole document is the state."""
    for text in (
        "const [envelope, setEnvelope] = useState<Envelope | null>(null);\n",
        "const [session, setSession] = useState<Session | null>(null);\n",
        "const [catalogue, setCatalogue] = useState<Catalogue | null>(null);\n",
    ):
        assert offenders(tmp_path, "App.tsx", text) == [], text


def test_a_comment_is_not_a_finding(tmp_path: Path) -> None:
    """The pages and the doc comments name these forms on purpose."""
    for text in (
        "// a cast here would say `as Envelope`, and 6894.757 is a factor\n",
        " * The mirrors, as TypeScript, are held to a fixture.\n",
        "/* useState<Graph> is what this replaced */\n",
    ):
        assert offenders(tmp_path, "App.tsx", text) == [], text


def test_a_marker_exempts_its_whole_block(tmp_path: Path) -> None:
    """The marker covers the block it opens, so a wrapped expression is exempt once."""
    text = (
        "const value = entity\n"
        "  .field  // ui-ok: the projection guarantees the prefix\n"
        "  .other as Role;\n"
    )
    assert offenders(tmp_path, "wire/nodes.ts", text) == []

    # And a marker below a blank line does not reach up to the next block.
    text = "// ui-ok: something else entirely\n\nconst value = entity as Role;\n"
    assert len(offenders(tmp_path, "wire/nodes.ts", text)) == 1


def test_a_tree_with_no_typescript_is_not_a_pass(tmp_path: Path, monkeypatch: Any) -> None:
    """A rule that reads nothing passes, so the empty tree is a failure rather than an OK."""
    tool: Any = lint()
    monkeypatch.setattr(tool, "ROOT", tmp_path)
    assert tool.main() == 1


def test_the_shipped_tree_passes() -> None:
    """The real rules against the real editor — the check's own subject."""
    tool: Any = lint()
    assert tool.main() == 0
