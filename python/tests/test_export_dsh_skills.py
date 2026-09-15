"""The dsh adapter emits one skill per catalog entry with the right frontmatter."""

import subprocess
import sys
import tomllib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
TOOL = ROOT / "tools" / "export_dsh_skills.py"
CATALOG = ROOT / "skills.toml"


def _names() -> list[str]:
    document = tomllib.loads(CATALOG.read_text(encoding="utf-8"))
    return [str(entry["name"]) for entry in document.get("skill", [])]


def test_export_writes_one_skill_per_entry(tmp_path: Path) -> None:
    result = subprocess.run(
        [sys.executable, str(TOOL), "--out", str(tmp_path)], capture_output=True, text=True
    )
    assert result.returncode == 0, result.stderr
    assert sorted(p.name for p in tmp_path.glob("*.md")) == sorted(f"{n}.md" for n in _names())


def test_frontmatter_names_and_describes_the_skill(tmp_path: Path) -> None:
    subprocess.run(
        [sys.executable, str(TOOL), "--out", str(tmp_path)], capture_output=True, check=True
    )
    document = tomllib.loads(CATALOG.read_text(encoding="utf-8"))
    entry = document["skill"][0]
    text = (tmp_path / f"{entry['name']}.md").read_text(encoding="utf-8")
    assert text.startswith("---\n")
    assert f"name: {entry['name']}\n" in text
    assert "description: " in text
