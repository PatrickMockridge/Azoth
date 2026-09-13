"""azoth - open, validated, citable chemical engineering calculations.

Every calculation ships with its equation, a source, the range in which it is
validated, its assumptions, a worked example and tests. The docs are generated
from the same machine-readable specs the code is generated from, so they cannot
drift apart.

# The idea worth knowing before anything else

Two things this library does differently from most engineering code, both of
which exist to stop a plausible-looking wrong number from being used:

**Warnings are not errors.** A value outside the range in which a correlation was
validated is still a value, and refusing to return it would be less useful than
returning it with a warning. What the library must never do is return it
silently. Check ``result.warnings``, or ``result.is_clean`` if you are willing to
see every caveat at once.

**A check that could not run is not a check that passed.** When an optional input
is missing, the range check that depends on it reports ``RANGE_CHECK_SKIPPED``
rather than quietly succeeding. "Checked and fine" and "never checked" are
different states and the API keeps them different.

# Start here

>>> import azoth
>>> q = azoth.ureg.Quantity
>>> r = azoth.hydraulics.reynolds_number(
...     q(998.0, "kg/m**3"), q(1.5, "m/s"), q(0.1, "m"), q(1.002e-3, "Pa*s")
... )
>>> round(r.re, 1)
149401.2
>>> str(r.regime)
'turbulent'

# Warnings before you use a result

>>> r = azoth.hydraulics.darcy_weisbach(
...     0.02, q(100.0, "m"), q(0.1, "m"), q(998.0, "kg/m**3"), q(1.5, "m/s")
... )
>>> str(r.warnings[0].code)
'RANGE_CHECK_SKIPPED'

# Not for design work yet

The fitting coefficients in ``data/fittings/crane_k_factors.csv`` are estimated
placeholder values, not engineering data. See that file's header; no warning on a
result announces the placeholder status.
"""

from __future__ import annotations

from typing import Any

from azoth import core
from azoth._dispatch import Backend, available, describe, extension, select, use_backend
from azoth.core.errors import (
    AzothError,
    InvalidInputError,
    OutOfRangeError,
    PropertyUnavailableError,
    SolverNotConvergedError,
    UnitMismatchError,
    UnknownFittingError,
    UnverifiedCalculationError,
)
from azoth.core.result import FlowRegime
from azoth.core.units import Q, ureg
from azoth.core.warnings import Warning, WarningCode

__version__ = "0.1.0"

# Imported last: `azoth.hydraulics` pulls in the dispatch layer, which imports
# this package's submodules. Keeping it at the bottom means everything it needs
# is already bound.
from azoth import batch, eos, hydraulics, properties, thermal

__all__ = [
    "AzothError",
    "Backend",
    "FlowRegime",
    "InvalidInputError",
    "OutOfRangeError",
    "PropertyUnavailableError",
    "Q",
    "SolverNotConvergedError",
    "UnitMismatchError",
    "UnknownFittingError",
    "UnverifiedCalculationError",
    "Warning",
    "WarningCode",
    "__version__",
    "available",
    "backends",
    "batch",
    "core",
    "describe",
    "eos",
    "extension",
    "hydraulics",
    "properties",
    "select",
    "thermal",
    "ureg",
    "use_backend",
]


def backends() -> frozenset[Backend]:
    """Which implementations are available in this installation.

    >>> import azoth
    >>> sorted(azoth.backends())
    ['python']
    """
    return available()


def active_backend() -> Backend:
    """Which implementation calculations will actually use."""
    return select()


def backend_info() -> dict[str, Any]:
    """Diagnostics: what is available, what is selected, and whether the compiled
    extension loaded. Useful in a bug report and in provenance metadata."""
    return describe()
