"""What the Rust backend's adapters need that is not an adapter.

`python/src/azoth/_rust_bridge.py` holds one function per registered id and a handful of helpers
they all call. The functions are generated now - `tools/gen_python_bridge.py` emits
`azoth/_rust_bridge_gen.py` - and a generator cannot emit a helper that is itself logic, so the
helpers live here and both modules import them.

**The move is what makes the generated module importable.** A generated module that imported
`_warnings` from `_rust_bridge` would be importing the file that re-exports it, which is the cycle
this file exists to break.
"""

from __future__ import annotations

from collections.abc import Sequence

from azoth import _core
from azoth.core.warnings import Warning, WarningCode

__all__ = ["_warnings"]


def _warnings(raw: Sequence[_core.Warning]) -> tuple[Warning, ...]:
    """Convert transported warnings to the real ones.

    Codes become the enum rather than staying strings, so a caller can compare
    ``warning.code == WarningCode.TRANSITIONAL_FLOW`` regardless of which backend
    produced it. A code the Python side does not know is impossible - the two
    sets are asserted equal by a test - but the ValueError would surface here if
    it ever happened, which is better than a string quietly failing to match.
    """
    return tuple(Warning(WarningCode(w.code), w.message, w.field) for w in raw)
