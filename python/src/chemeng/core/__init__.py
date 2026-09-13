"""Core types shared by every calculation.

Deliberately contains no engineering calculations, mirroring ``chemeng-core`` in
Rust. It holds the vocabulary the calcs share - units, warnings, errors, results,
range checks - so every calculation reports failures, caveats and units the same
way.

The two ideas worth understanding before reading anything else:

1. **Warnings are not errors.** A value outside the range in which a correlation
   was validated is still a value, and refusing to return it would be less useful
   than returning it with a warning. What the library must never do is return it
   *silently*.
2. **A check that could not run is not a check that passed.** When an optional
   input is missing, the range check that depends on it reports
   ``RANGE_CHECK_SKIPPED`` rather than quietly succeeding.
"""

from __future__ import annotations

from chemeng.core.errors import (
    ChemEngError,
    InvalidInputError,
    OutOfRangeError,
    PropertyUnavailableError,
    SolverNotConvergedError,
    UnitMismatchError,
    UnknownFittingError,
    UnverifiedCalculationError,
)
from chemeng.core.range import Band, RangeCheck, Severity, SpecChecks, apply_checks, checks_for
from chemeng.core.result import (
    RESULT_TYPES,
    ColebrookResult,
    DarcyWeisbachResult,
    FlowRegime,
    KComponent,
    KFactorsResult,
    ReynoldsNumberResult,
    SwameeJainResult,
)
from chemeng.core.units import CANONICAL_UNITS, Q, quantity, to_si, unit_for, ureg
from chemeng.core.warnings import Warning, WarningCode

__all__ = [
    "CANONICAL_UNITS",
    "RESULT_TYPES",
    "Band",
    "ChemEngError",
    "ColebrookResult",
    "DarcyWeisbachResult",
    "FlowRegime",
    "InvalidInputError",
    "KComponent",
    "KFactorsResult",
    "OutOfRangeError",
    "PropertyUnavailableError",
    "Q",
    "RangeCheck",
    "ReynoldsNumberResult",
    "Severity",
    "SolverNotConvergedError",
    "SpecChecks",
    "SwameeJainResult",
    "UnitMismatchError",
    "UnknownFittingError",
    "UnverifiedCalculationError",
    "Warning",
    "WarningCode",
    "apply_checks",
    "checks_for",
    "quantity",
    "to_si",
    "unit_for",
    "ureg",
]
