"""Core types shared by every calculation.

Deliberately contains no engineering calculations, mirroring ``azoth-core`` in
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

:mod:`azoth.core.solver` is the one module here that is *machinery* rather than
vocabulary, and it lives here for the same reason: a named solution scheme is
shared the way a unit is, because two implementations of one equation must run the
same scheme to agree numerically. Like :mod:`azoth.core.range`, it is not a
calculation - it has no equation and no spec - and it grows only when the spec
schema's ``solver.kind`` grows.
"""

from __future__ import annotations

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
from azoth.core.range import Band, RangeCheck, Severity, SpecChecks, apply_checks, checks_for
from azoth.core.result import (
    RESULT_TYPES,
    BubblePressureResult,
    ColebrookResult,
    DarcyWeisbachResult,
    DewPressureResult,
    FlowRegime,
    IdealGasCpResult,
    KComponent,
    KFactorsResult,
    MolarEnthalpyEntropyResult,
    Phase,
    PrAlphaAbResult,
    PrDepartureResult,
    PrKappaResult,
    PrMassDensityResult,
    PrMolarVolumeResult,
    PrsvKappaResult,
    PrZFactorResult,
    PtFlashResult,
    PureSaturationResult,
    ReynoldsNumberResult,
    RootStructure,
    StabilityTestResult,
    StabilityVerdict,
    SwameeJainResult,
    Vdw1fMixBinaryResult,
)
from azoth.core.solver import (
    Convergence,
    CubicRootsOutcome,
    SolverKind,
    SolverOutcome,
    cubic_roots,
    fixed_point,
    require_converged,
    require_cubic_converged,
)
from azoth.core.units import CANONICAL_UNITS, Q, from_si, quantity, to_si, unit_for, ureg
from azoth.core.warnings import Warning, WarningCode

__all__ = [
    "CANONICAL_UNITS",
    "RESULT_TYPES",
    "AzothError",
    "Band",
    "BubblePressureResult",
    "ColebrookResult",
    "Convergence",
    "CubicRootsOutcome",
    "DarcyWeisbachResult",
    "DewPressureResult",
    "FlowRegime",
    "IdealGasCpResult",
    "InvalidInputError",
    "KComponent",
    "KFactorsResult",
    "MolarEnthalpyEntropyResult",
    "OutOfRangeError",
    "Phase",
    "PrAlphaAbResult",
    "PrDepartureResult",
    "PrKappaResult",
    "PrMassDensityResult",
    "PrMolarVolumeResult",
    "PrZFactorResult",
    "PropertyUnavailableError",
    "PrsvKappaResult",
    "PtFlashResult",
    "PureSaturationResult",
    "Q",
    "RachfordRiceBinaryResult",
    "RangeCheck",
    "ReynoldsNumberResult",
    "RootStructure",
    "Severity",
    "SolverKind",
    "SolverNotConvergedError",
    "SolverOutcome",
    "SpecChecks",
    "StabilityTestResult",
    "StabilityVerdict",
    "SwameeJainResult",
    "UnitMismatchError",
    "UnknownFittingError",
    "UnverifiedCalculationError",
    "Vdw1fMixBinaryResult",
    "Warning",
    "WarningCode",
    "apply_checks",
    "checks_for",
    "cubic_roots",
    "fixed_point",
    "from_si",
    "quantity",
    "require_converged",
    "require_cubic_converged",
    "to_si",
    "unit_for",
    "ureg",
]
