"""Result shapes.

Mirrors ``azoth_hydraulics::results`` in Rust, field for field. The field names
are a cross-language contract: ``test_registry_contract.py`` asserts that each
Python result dataclass has exactly the fields the corresponding Rust
``CalcResult::FIELDS`` lists. A result whose Python and Rust shapes differ is a
bug no numerical test can catch, because the numbers agree perfectly and only the
attribute name is wrong.

# Why ``eq=False``

The dataclasses are frozen but not comparable with ``==``. Generating
``__eq__`` would compare ``pint`` quantities, whose equality semantics are not
what a numerical test wants anyway: results should be compared within a
tolerance, not bit-wise. :func:`azoth.testing.assert_results_equal` does that
explicitly, which also makes the tolerance visible at every call site.

# Why warnings are tuples

Immutable, and they match ``Vec<Warning>`` semantically. A result whose warnings
can be edited after the fact is a result whose caveats can be removed by whoever
finds them inconvenient.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import StrEnum

from azoth.core.units import Q
from azoth.core.warnings import Warning, WarningCode

#: Upper bound of the laminar band, exclusive. Mirrors
#: ``FlowRegime::LAMINAR_MAX`` in Rust.
#:
#: Module level, not class attributes: any plain assignment inside an ``Enum``
#: body becomes an enum member, and a float member would make this a broken
#: ``StrEnum``.
LAMINAR_MAX: float = 2000.0

#: Lower bound of the turbulent band, exclusive.
TURBULENT_MIN: float = 4000.0


class FlowRegime(StrEnum):
    """Flow regime of a pipe flow, under the Crane/Moody boundaries.

    A note on terminology, because this genuinely trips people up. Crane and
    Moody use "transition zone" to mean something *different* from the 2000-4000
    regime: for them it is the region on the Moody chart between the smooth-flow
    line and complete turbulence, where the friction factor still depends on
    Reynolds number. This enum uses the more common modern sense.

    The boundaries are approximate engineering guidance, not exact physical
    transitions, which is why ``TRANSITIONAL`` is a warning condition everywhere
    it appears rather than a clean category.
    """

    LAMINAR = "laminar"
    TRANSITIONAL = "transitional"
    TURBULENT = "turbulent"

    @classmethod
    def from_reynolds_number(cls, re: float) -> FlowRegime:
        """Classify a Reynolds number.

        Both band edges fall inside ``TRANSITIONAL``, so exactly 2000 or exactly
        4000 is transitional rather than laminar or turbulent. That is a
        convention: the boundaries are not sharp in reality, and the conservative
        reading is the useful one.
        """
        if re < LAMINAR_MAX:
            return cls.LAMINAR
        if re <= TURBULENT_MIN:
            return cls.TRANSITIONAL
        return cls.TURBULENT

    @property
    def is_indeterminate(self) -> bool:
        """Whether the friction factor is indeterminate in this regime."""
        return self is FlowRegime.TRANSITIONAL


class _HasWarnings:
    """Shared warning accessors, mirroring ``CalcResult`` on the Rust side.

    A mixin rather than five copies. ``__slots__`` is empty so that
    ``slots=True`` dataclasses inheriting from it still get real slots - a base
    with a ``__dict__`` would quietly give every result instance one, undoing the
    point of declaring slots.
    """

    __slots__ = ()

    warnings: tuple[Warning, ...]

    @property
    def is_clean(self) -> bool:
        """True when the result carries no warnings at all."""
        return not self.warnings

    def has_warning(self, code: WarningCode) -> bool:
        """True when any warning carries the given code."""
        return any(w.code == code for w in self.warnings)


@dataclass(frozen=True, slots=True, eq=False)
class ReynoldsNumberResult(_HasWarnings):
    """Result of ``hydraulics.reynolds_number``."""

    #: Reynolds number. Dimensionless.
    re: float
    #: Flow regime under the Crane/Moody boundaries.
    regime: FlowRegime
    #: Caveats. Carries ``TRANSITIONAL_FLOW`` when the flow sits in 2000-4000.
    warnings: tuple[Warning, ...]


@dataclass(frozen=True, slots=True, eq=False)
class ColebrookResult(_HasWarnings):
    """Result of ``hydraulics.friction_factor_colebrook``.

    Carries the solver's own report because the answer is only meaningful
    alongside it: ``f`` is the last iterate, and if ``converged`` is false it is
    not a solution to the equation at all.
    """

    #: Darcy friction factor. Dimensionless.
    f: float
    #: Iterations performed.
    iterations: int
    #: Whether the iteration met its tolerance.
    converged: bool
    #: Final change between iterates.
    residual: float
    #: Caveats.
    warnings: tuple[Warning, ...]


@dataclass(frozen=True, slots=True, eq=False)
class SwameeJainResult(_HasWarnings):
    """Result of ``hydraulics.friction_factor_swamee_jain``.

    No solver report: the Swamee-Jain equation is explicit, and reporting an
    iteration count of zero for it would imply a similarity to Colebrook that
    does not exist.
    """

    #: Darcy friction factor. Dimensionless.
    f: float
    #: Caveats.
    warnings: tuple[Warning, ...]


@dataclass(frozen=True, slots=True, eq=False)
class HaalandResult(_HasWarnings):
    """Result of ``hydraulics.friction_factor_haaland``.

    Shaped identically to :class:`SwameeJainResult`, and for the same reason: the
    Haaland equation is explicit, so there is no iteration to report. The two are
    separate classes rather than one shared type because the spec-to-result
    contract is asserted per calc id, and a shared class would make it impossible
    to tell which calc a result came from.
    """

    #: Darcy friction factor. Dimensionless.
    f: float
    #: Caveats.
    warnings: tuple[Warning, ...]


@dataclass(frozen=True, slots=True, eq=False)
class KComponent:
    """One fitting's contribution to the total resistance coefficient."""

    #: Fitting id as it appears in the registry.
    fitting_id: str
    #: Equivalent length ratio (L_eq / D) for this fitting.
    n_ld: float
    #: This fitting's resistance coefficient, ``f_t * n_ld``.
    k: float


@dataclass(frozen=True, slots=True, eq=False)
class KFactorsResult(_HasWarnings):
    """Result of ``hydraulics.crane_k_factors``."""

    #: Total resistance coefficient for all listed fittings. Dimensionless.
    k_total: float
    #: The friction factor the coefficients were based on.
    f_t: float
    #: Per-fitting breakdown, in the order the fittings were supplied.
    components: tuple[KComponent, ...]
    #: Caveats. Carries ``ESTIMATED_DATA`` while the registry holds placeholders.
    warnings: tuple[Warning, ...]


@dataclass(frozen=True, slots=True, eq=False)
class DarcyWeisbachResult(_HasWarnings):
    """Result of ``hydraulics.darcy_weisbach``."""

    #: Pressure drop over the pipe length.
    dp: Q
    #: The friction factor the drop was computed with.
    f: float
    #: Reynolds number, present only when viscosity was supplied.
    re: float | None
    #: Flow regime, present only when viscosity was supplied.
    regime: FlowRegime | None
    #: Caveats. Carries ``RANGE_CHECK_SKIPPED`` when viscosity was omitted,
    #: because the regime then went unchecked.
    warnings: tuple[Warning, ...]


#: Calc id -> the result dataclass it produces. Used by the contract test to
#: check each shape against the Rust side without importing every name by hand.
RESULT_TYPES: dict[str, type[object]] = {
    "hydraulics.reynolds_number": ReynoldsNumberResult,
    "hydraulics.friction_factor_colebrook": ColebrookResult,
    "hydraulics.friction_factor_swamee_jain": SwameeJainResult,
    "hydraulics.friction_factor_haaland": HaalandResult,
    "hydraulics.crane_k_factors": KFactorsResult,
    "hydraulics.darcy_weisbach": DarcyWeisbachResult,
}
