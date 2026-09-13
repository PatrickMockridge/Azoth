"""Pure-Python reference implementations of the equations-of-state calculations.

Every calculation in azoth has two implementations, and this is the Python one. It
is not a fallback for a missing extension: it is the second opinion the Rust core
is checked against, and the test suite runs both against every case in the specs.
See ``README.md`` for why a single implementation was rejected.

The bodies here are written to mirror the Rust line for line, so a reviewer can
read the two side by side against the published equation. Where they diverge, it is
a bug in one of them - and the cross-implementation tests are what catch it.
"""

from __future__ import annotations

from azoth.eos.reference.pr_kappa import pr_kappa

__all__ = [
    "pr_kappa",
]
