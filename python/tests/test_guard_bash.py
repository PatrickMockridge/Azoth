"""The Bash guard's decisions, held in the suite as well as in ``--selftest``.

``guard_bash.py --selftest`` is the hook's own gate; this keeps the security-critical
decisions visible to ``test_scope`` and to CI, where the selftest is not run. The
regressions it pins are the ones that made an hour-long command look sanctioned: the
wrapper skipped the wide check, and ``timeout`` hid the program behind it.
"""

from __future__ import annotations

import importlib.util
from pathlib import Path
from types import ModuleType

HOOK = Path(__file__).resolve().parents[2] / ".claude" / "hooks" / "guard_bash.py"


def _guard() -> ModuleType:
    spec = importlib.util.spec_from_file_location("guard_bash", HOOK)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _allowed(command: str) -> bool:
    allow, _ = _guard().decide(command)
    return bool(allow)


def test_a_wrapped_whole_suite_is_still_refused() -> None:
    """The wrapper satisfies the memory cap; it does not make a wide run scoped."""
    assert not _allowed("tools/gated.sh pytest python/tests")
    assert not _allowed("tools/gated.sh cargo test --workspace -j 10")
    assert not _allowed("tools/gated.sh uv run pytest python/tests")
    assert not _allowed("systemd-run --user --scope -q pytest python/tests")


def test_a_timeout_does_not_hide_the_program() -> None:
    """`timeout 300 pytest python/tests` is the same sweep with a word in front."""
    assert not _allowed("timeout 300 pytest python/tests")
    assert not _allowed("tools/gated.sh timeout 300 pytest python/tests")
    assert _allowed(
        "timeout 300 tools/gated.sh python -m pytest "
        "python/tests/models/test_model_cases.py -k flash"
    )


def test_the_whole_tree_needs_the_visible_token() -> None:
    """The token is read before its assignment is stripped, and only where it was written."""
    assert _allowed("AZOTH_FULL_GATE=1 tools/gated.sh cargo test --workspace -j 10")
    assert not _allowed("cargo test --workspace")
    assert not _allowed("echo AZOTH_FULL_GATE=1; cargo test --workspace")


def test_the_scoped_runner_refuses_its_own_whole_tree_modes() -> None:
    """A local `--run` is the edit lane; the whole tree is nightly and needs the token."""
    assert _allowed("python tools/test_scope.py --run")
    assert _allowed("python tools/test_scope.py --lane tranche --run")
    assert _allowed("python tools/test_scope.py --lane nightly --no-gate")
    assert not _allowed("python tools/test_scope.py --lane nightly --run")
    assert not _allowed("python tools/test_scope.py --full --run")
    assert _allowed("AZOTH_FULL_GATE=1 python tools/test_scope.py --full --run")


def test_a_scope_without_a_memory_cap_is_not_the_wrapper() -> None:
    """A bare `systemd-run --scope` runs uncapped, which is what the wrapper prevents."""
    assert not _allowed("systemd-run --user --scope -q java -cp .:neqsim.jar Probe")
    assert _allowed(
        "systemd-run --user --scope -q -p MemoryMax=12G -p MemorySwapMax=0 "
        "java -cp .:neqsim.jar Probe"
    )


def test_a_package_runner_resolves_to_the_program_it_runs() -> None:
    assert not _allowed("npx cargo test")
    assert not _allowed("npx -y pytest python/tests")


def test_the_lane_budget_can_only_be_lowered() -> None:
    """The 2026-10-22 bypass: `--max-minutes 35` against a 15-minute push lane."""
    assert not _allowed("python tools/test_scope.py --lane push --run --max-minutes 35")
    assert not _allowed("python tools/test_scope.py --lane tranche --run --max-minutes 21")
    assert not _allowed("python tools/test_scope.py --lane push --run --max-minutes 0")
    assert _allowed("python tools/test_scope.py --lane push --run --max-minutes 15")
    assert _allowed("python tools/test_scope.py --lane tranche --run --max-minutes 20")


def test_an_hour_long_heavy_command_is_refused() -> None:
    """The per-command ceiling is 30 minutes, whatever a `timeout` in front says."""
    assert not _allowed("timeout 3600 java -cp .:neqsim.jar Probe")
    assert not _allowed(
        "tools/gated.sh timeout 3600 python -m pytest "
        "python/tests/models/test_model_cases.py -k flash"
    )
    assert _allowed("tools/gated.sh timeout 900 java -cp .:neqsim.jar Probe")


def test_a_script_that_runs_a_heavy_program_is_refused(tmp_path: Path) -> None:
    heavy = tmp_path / "probe.sh"
    heavy.write_text("java -cp x Probe\n", encoding="utf-8")
    assert not _allowed(f"bash {heavy}")
    clean = tmp_path / "clean.sh"
    clean.write_text("echo hello\n", encoding="utf-8")
    assert _allowed(f"bash {clean}")


def test_read_only_commands_that_name_tests_are_allowed() -> None:
    """A quoted pattern is not an operator; a mention is not a run."""
    assert _allowed('grep -n -E "cargo test|pytest" .github/workflows/ci.yml')
    assert _allowed("grep -n oracle_sweep tools/")
    assert _allowed("git log --grep pytest")


def test_the_hooks_own_selftest_passes() -> None:
    assert _guard().selftest() == 0
