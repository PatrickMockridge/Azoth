"""The port tooling `SKILL.md` documents must still exist."""

from pathlib import Path

GENERATORS = (
    "gen_registry.py",
    "gen_models.py",
    "gen_vocabulary.py",
    "gen_docs.py",
    "gen_stub.py",
)


def test_the_generators_exist() -> None:
    tools = Path(__file__).resolve().parents[4] / "tools"
    missing = [name for name in GENERATORS if not (tools / name).is_file()]
    assert not missing, f"the port's generators are missing: {missing}"
