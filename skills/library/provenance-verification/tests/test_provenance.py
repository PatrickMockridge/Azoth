"""The provenance `SKILL.md` points at must still exist."""

from pathlib import Path


def test_notice_exists_and_credits_neqsim() -> None:
    notice = Path(__file__).resolve().parents[4] / "NOTICE"
    assert notice.is_file()
    assert "NeqSim" in notice.read_text(encoding="utf-8")
