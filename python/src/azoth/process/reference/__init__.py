"""Pure-Python reference implementations of the unit operations.

Every calculation in azoth has two implementations, and this is the Python one. It is
not a fallback for a missing extension: it is the second opinion the Rust core is
checked against, and the test suite runs both against every case in the specs. See
``README.md`` for why a single implementation was rejected.

A unit operation is a *model*, so its reference has the same standing as
``azoth.eos.reference``'s and is held to the same contract - except for one thing worth
knowing before reading one: these take a whole stream rather than a mixture state, so
they carry a molar flow and return two streams where the equation-of-state models carry
one and return a state.

The bodies here are written to mirror the Rust line for line, so a reviewer can read
the two side by side against the NeqSim source they were ported from. Where they
diverge, it is a bug in one of them - and the cross-implementation tests are what catch
it.
"""

from __future__ import annotations

from azoth.process.reference.separator import separator

__all__ = [
    "separator",
]
