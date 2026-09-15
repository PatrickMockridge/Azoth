"""The skill exporter writes one frontmatter file per catalog entry, per target."""

import subprocess
import sys
import tomllib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
TOOL = ROOT / "tools" / "export_skills.py"
CATALOG = ROOT / "skills.toml"


def _names() -> list[str]:
    document = tomllib.loads(CATALOG.read_text(encoding="utf-8"))
    return [str(entry["name"]) for entry in document.get("skill", [])]


def test_export_writes_one_skill_per_entry_per_target(tmp_path: Path) -> None:
    result = subprocess.run(
        [sys.executable, str(TOOL), "--target", "all", "--root", str(tmp_path)],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stderr
    names = _names()
    dsh = sorted(p.name for p in (tmp_path / ".dsh" / "skills").glob("*.md"))
    assert dsh == sorted(f"{n}.md" for n in names)
    for tool in ("claude", "opencode"):
        written = sorted(
            p.parent.name for p in (tmp_path / f".{tool}" / "skills").glob("*/SKILL.md")
        )
        assert written == sorted(names)


def test_frontmatter_names_and_describes_the_skill(tmp_path: Path) -> None:
    subprocess.run(
        [sys.executable, str(TOOL), "--target", "claude", "--root", str(tmp_path)],
        capture_output=True,
        check=True,
    )
    document = tomllib.loads(CATALOG.read_text(encoding="utf-8"))
    entry = document["skill"][0]
    path = tmp_path / ".claude" / "skills" / entry["name"] / "SKILL.md"
    text = path.read_text(encoding="utf-8")
    assert text.startswith("---\n")
    assert f"name: {entry['name']}\n" in text
    assert "description: " in text
