"""The orchestration layer's DeepSeek Harness boundary.

azoth owns the orchestration *specification* - the rho-calculus and the `agents/`
definitions; DeepSeek Harness is the *executor*. This module is the one place in
`python/src/azoth` that may import `deepseek_harness_sdk`, so the dev-preview runtime
stays behind a thin boundary and nothing else in the library depends on it.

The adapter that makes azoth's skills loadable in dsh is
`tools/export_dsh_skills.py`; the swarm/subagent orchestration API lands with the
first agent, not here.
"""

from __future__ import annotations

from typing import Any


def runtime() -> Any:
    """The DeepSeek Harness runtime, imported on first use.

    Raises:
        ImportError: if `deepseek-harness-sdk` is not installed. The skill export
            (`tools/export_dsh_skills.py`) does not need it; only running an
            orchestration does, which is why it is the optional `agent` extra rather
            than a dependency.
    """
    try:
        import deepseek_harness_sdk as dsh
    except ImportError as exc:
        raise ImportError(
            "azoth orchestration needs `deepseek-harness-sdk`; install `azoth[agent]`"
        ) from exc
    return dsh
