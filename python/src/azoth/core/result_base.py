"""What the result dataclasses are built from, written by hand.

:mod:`azoth.core.result_gen` carries the 194 dataclasses, one per registered id, and it is
generated - but three of the things they need are decisions rather than derivations, and they live
here so that the generated module does not have to import the module that re-exports it.

* the **enums** a result field may be: which values each has, and which of them is a root selector.
* the **``_HasWarnings`` mixin** every result inherits, which is what makes a result self-describing
  - ``provenance``, ``to_dict``, ``to_json``, ``is_clean`` and ``has_warning``.
* the **one record that is not a result**: `KComponent`, a fitting's contribution to a resistance
  coefficient. It is a frozen dataclass with the same decorator and no ``CALC_ID``, so it is not in
  the registry and the generator has nothing to read for it.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import StrEnum
from typing import Any, ClassVar, cast

from azoth._provenance_gen import provenance as provenance_table
from azoth.core import serialise
from azoth.core.provenance import PROVENANCE_KEY, Provenance
from azoth.core.warnings import Warning, WarningCode

LAMINAR_MAX: float = 2000.0

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


class RootStructure(StrEnum):
    """How many admissible real roots a cubic equation of state had.

    Reported so a caller can tell a single-root state from one where ``z_min`` and
    ``z_max`` are genuinely two different roots, without comparing floats.

    **There is no ``TWO_ROOTS``, and that is a theorem rather than an omission.**
    For the Peng-Robinson cubic, ``f(b_reduced)`` is exactly ``-2*b_reduced**2`` -
    the algebra is in ``specs/calcs/eos/pr_z_factor.toml`` - so ``b_reduced`` lies
    either below all three roots or between the middle and the largest one. The
    admissible count is therefore 1 or 3 and never 2. A variant that cannot occur
    would be a value a caller branches on and never sees, which is worse than an
    absent one.
    """

    ONE_ROOT = "one_root"
    THREE_ROOTS = "three_roots"


class Phase(StrEnum):
    """What a converged flash turned out to be.

    A separate type from :class:`RootStructure`, which counts the roots of a *pure*
    component's cubic and says nothing about phases. The two are easy to confuse and
    mean different things: three roots is a mathematical fact about a polynomial,
    ``TWO_PHASE`` is a physical claim about a mixture.

    All four values are reachable, and each is covered by a test - a variant a caller
    branches on and never sees is worse than an absent one, which is the same
    argument that removed ``TWO_ROOTS``.
    """

    #: A genuine split: ``beta`` in ``[0, 1]`` and the two compositions differ.
    TWO_PHASE = "two_phase"

    #: The feed is subcooled liquid. Reached two ways, and ``beta`` distinguishes
    #: them: either the converged Rachford-Rice root is negative, in which case
    #: ``beta`` is present as the negative-flash value; or every K-value is below
    #: one, in which case no root exists at all and ``beta`` is absent.
    ALL_LIQUID = "all_liquid"

    #: The feed is superheated vapour. The same two routes as ``ALL_LIQUID``.
    ALL_VAPOUR = "all_vapour"

    #: The iteration converged to ``x = y = z``.
    #:
    #: The feed is single phase, and **this does not say which one** - the K-values
    #: straddled one throughout, so nothing in the model ever proved which phase the
    #: feed is. That is what a tangent-plane stability analysis decides, and this
    #: model has none. ``beta`` is absent, because at the trivial solution it is
    #: indeterminate rather than out of range.
    TRIVIAL = "trivial"


class _HasWarnings:
    """Shared warning accessors, mirroring ``CalcResult`` on the Rust side.

    A mixin rather than five copies. ``__slots__`` is empty so that
    ``slots=True`` dataclasses inheriting from it still get real slots - a base
    with a ``__dict__`` would quietly give every result instance one, undoing the
    point of declaring slots.
    """

    __slots__ = ()

    warnings: tuple[Warning, ...]

    #: The spec id this result is the answer to, mirroring ``CalcResult::CALC_ID`` in Rust.
    #:
    #: A `ClassVar` and not a field: ``dataclasses.fields`` ignores it, so the field contract
    #: ``test_cross_impl`` holds the two languages to is unaffected and ``slots=True`` does
    #: not give it a slot. Declared with no value here so that a subclass omitting its own is
    #: a type error under strict mypy rather than an ``AttributeError`` at the first call to
    #: :attr:`provenance` - "every result is self-describing" is then a property of the type
    #: system rather than of the day the codemod ran.
    CALC_ID: ClassVar[str]

    @property
    def provenance(self) -> Provenance:
        """What this answer rests on: the generated half, and what this call reported.

        A property rather than a method, so it reads beside ``warnings`` and ``is_clean``.

        Raises:
            KeyError: if ``CALC_ID`` is not a registered id. Unreachable for a result this
                library ships, and a refusal rather than an empty block - a block that
                existed but said nothing would let a number be presented as though it
                carried its provenance when it carried nothing.
        """
        return Provenance.of(provenance_table(type(self).CALC_ID), self.warnings)

    def to_dict(self) -> dict[str, Any]:
        """This result as JSON-shaped data.

        **Inherited rather than written 187 times.** The walk is :mod:`azoth.core.serialise`'s, over
        ``dataclasses.fields``, so a result is serialisable by being a frozen dataclass and nothing
        else - which is what makes "a result as JSON" one writer rather than a codec per model.
        """
        # The walk's return is `Any` - it writes whatever a field holds - and the dataclass branch
        # is the one a result takes, which is a mapping.
        document = cast("dict[str, Any]", serialise.to_dict(self))
        # Through the same walk rather than beside it: the block is a frozen dataclass, so
        # `serialise` writes it too and a result as JSON is still one writer.
        document[PROVENANCE_KEY] = self.provenance.to_dict()
        return document

    def to_json(self) -> str:
        """This result as the JSON document, for a file, a cell or a wire.

        Through :meth:`to_dict`, so the document and the mapping are the same answer. Written
        as ``serialise.to_json(self)`` it would walk the dataclass's fields directly and the
        provenance block - which is not a field - would be in one and not the other.
        """
        return serialise.to_json(self.to_dict())

    @property
    def is_clean(self) -> bool:
        """True when the result carries no warnings at all."""
        return not self.warnings

    def has_warning(self, code: WarningCode) -> bool:
        """True when any warning carries the given code."""
        return any(w.code == code for w in self.warnings)


class HydrateStructure(StrEnum):
    """Which of the two hydrate structures a state's cages are.

    The structure is an *output* of the same comparison the formation temperature is:
    NeqSim evaluates both and keeps the lower water fugacity coefficient, so a caller
    reading an occupancy needs it and a model that reported only the temperature would
    leave the cages unnamed.
    """

    #: Two small ``5^12`` and six large ``5^12 6^2`` cavities per forty-six waters.
    STRUCTURE_I = "structure_i"
    #: Sixteen small ``5^12`` and eight large ``5^12 6^4`` per a hundred and thirty-six.
    STRUCTURE_II = "structure_ii"


@dataclass(frozen=True, slots=True, eq=False)
class KComponent:
    """One fitting's contribution to the total resistance coefficient."""

    #: Fitting id as it appears in the registry.
    fitting_id: str
    #: Equivalent length ratio (L_eq / D) for this fitting.
    n_ld: float
    #: This fitting's resistance coefficient, ``f_t * n_ld``.
    k: float


class StabilityVerdict(StrEnum):
    """Whether a feed is stable as a single phase.

    Two values rather than a boolean because the *asymmetry* between them is the
    point: `UNSTABLE` is a proof - a trial reached a stationary point below the
    tangent plane, so a single phase is not the Gibbs minimum - while `STABLE` is
    the absence of one, from two trials that were placed by a gas-liquid
    correlation. A caller who reads `STABLE` as "no split exists" has read it wrong,
    and a bare `True` invites exactly that.
    """

    STABLE = "stable"
    UNSTABLE = "unstable"


class TpMultiflashSeed(StrEnum):
    """Which of the two things decided the phase count.

    A boolean would say a phase was added; this says whether the answer *is* the two-phase
    flash's or is something the tangent-plane trial found, which is the distinction between
    "this model is `eos.pt_flash`" and "this model found a third phase".
    """

    TWO_PHASE_FLASH = "two_phase_flash"
    STABILITY_SEEDED = "stability_seeded"


class HenryStatus(StrEnum):
    """Whether the guideline's Henry constant was evaluated inside the fitted range.

    A returned number rather than a refusal: the equation is defined over the whole of
    liquid water and the guideline's own extrapolating entry point exists for the rest.
    What an extrapolation means is the consumer's decision - NeqSim's consumers turn it
    into the insoluble limit - and this says which case applies so that the decision is
    visible at the call site.
    """

    #: The temperature is inside the range the gas's row was fitted over.
    WITHIN_FITTED_RANGE = "within_fitted_range"

    #: The equation is defined, but the temperature is outside the row's fitted range.
    GUIDELINE_EXTRAPOLATION = "guideline_extrapolation"
