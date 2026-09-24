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
class PumpPowerResult(_HasWarnings):
    """Result of ``hydraulics.pump_power``."""

    #: Shaft power the pump must be supplied with.
    power: Q
    #: Caveats.
    warnings: tuple[Warning, ...]


@dataclass(frozen=True, slots=True, eq=False)
class OrificeFlowResult(_HasWarnings):
    """Result of ``hydraulics.orifice_flow``."""

    #: Volumetric flow rate through the orifice.
    q: Q
    #: Caveats.
    warnings: tuple[Warning, ...]


@dataclass(frozen=True, slots=True, eq=False)
class ControlValveCvResult(_HasWarnings):
    """Result of ``hydraulics.control_valve_cv``."""

    #: Volumetric flow rate through the valve.
    q: Q
    #: Caveats.
    warnings: tuple[Warning, ...]


@dataclass(frozen=True, slots=True, eq=False)
class ChokedFlowAreaResult(_HasWarnings):
    """Result of ``hydraulics.choked_flow_area``."""

    #: Throat area required for the choked flow.
    a: Q
    #: Caveats.
    warnings: tuple[Warning, ...]


@dataclass(frozen=True, slots=True, eq=False)
class ConductionPlaneWallResult(_HasWarnings):
    """Result of ``thermal.conduction_plane_wall``.

    The first result in the registry to carry a dimensioned output from a
    namespace other than hydraulics, which is what made a second domain worth
    having: nothing about this class is hydraulics-shaped.
    """

    #: Heat flow rate through the wall. Signed, and it follows the sign of the
    #: temperature difference rather than being reported as a magnitude.
    q: Q
    #: Caveats.
    warnings: tuple[Warning, ...]


@dataclass(frozen=True, slots=True, eq=False)
class PrKappaResult(_HasWarnings):
    """Result of ``eos.pr_kappa``.

    Both the input and the output are dimensionless, so this is the first result
    in the registry with no ``pint`` quantity in it at all - which is what an
    equation of state written in reduced variables looks like at a boundary. There
    is no conversion to get wrong because there is no unit to convert.
    """

    #: The Peng-Robinson alpha-function coefficient. Dimensionless, and a property
    #: of the substance alone: no temperature, no pressure.
    kappa: float
    #: Caveats.
    warnings: tuple[Warning, ...]


@dataclass(frozen=True, slots=True, eq=False)
class EquilibriumConstantResult(_HasWarnings):
    """Result of ``reactions.equilibrium_constant``.

    The first result from a namespace that is not about a fluid's own state: its
    inputs are a name and a temperature, and the arithmetic behind them is a fitted
    correlation read from the databank rather than a constitutive equation.
    """

    #: ``ln K`` at the caller's temperature, from the row's four coefficients.
    ln_k: float
    #: ``K``, the exponential of ``ln_k``. Dimensionless, because the correlation's
    #: standard state is what makes the constant so.
    k: float
    #: ``d(ln K)/dT``, van 't Hoff's derivative. The heat is this times ``R T**2``.
    ln_k_derivative: Q
    #: ``d(ln K)/dT * R * T**2``, in J/mol. **The sign is about the reaction as the
    #: table writes it**: a negative value is exothermic in that direction and says
    #: nothing about the reverse.
    reaction_heat: Q
    #: The row's own literature citation, which differs between the sources - so which
    #: fit answered is visible rather than only in the table.
    reference: str
    #: Caveats.
    warnings: tuple[Warning, ...]


@dataclass(frozen=True, slots=True, eq=False)
class ChemicalEquilibriumResult(_HasWarnings):
    """Result of ``reactions.chemical_equilibrium``.

    **An unconverged solve is a result and not an exception.** NeqSim returns its last
    iterate whatever happens and its callers read it, so ``converged`` is what says
    whether to believe the composition rather than the composition being withheld - and
    one of the two captured fluids never converges.
    """

    #: The moles of each species at the answer, or where the solve gave up.
    moles: tuple[Q, ...]
    #: Passes taken.
    iterations: int
    #: ``sum(|dn_i| / n_i)`` over the species that moved on the last pass.
    error: float
    #: Whether the error came in under the tolerance.
    converged: bool
    #: Caveats.
    warnings: tuple[Warning, ...]


@dataclass(frozen=True, slots=True, eq=False)
class HeaterResult(_HasWarnings):
    """Result of ``process.heater``.

    The same record as :class:`PumpResult` with a duty beside it, which is the one number a
    heater adds that its outlet record does not carry. **The duty is the state's and not the
    caller's**: ``Heater.run`` overwrites its own ``energyInput`` with the enthalpy the flash
    reached, in every branch, so it is ``outlet_n * (outlet_h - inlet_h)`` however the outlet
    was specified.
    """

    #: Molar flow out, which is the inlet's.
    outlet_n: Q
    #: Outlet composition, one entry per component.
    outlet_z: tuple[float, ...]
    #: Outlet pressure: the inlet's less ``pressure_drop``.
    outlet_p: Q
    #: Outlet temperature: the stated one, or the one the shifted enthalpy reaches.
    outlet_t: Q
    #: Outlet molar enthalpy.
    outlet_h: Q
    #: The duty moved, as ``outlet_n * (outlet_h - inlet_h)``.
    outlet_duty: Q
    #: Caveats.
    warnings: tuple[Warning, ...]


@dataclass(frozen=True, slots=True, eq=False)
class PumpResult(_HasWarnings):
    """Result of ``process.pump``.

    **A unit operation's result is its outlet's record**, spelled out field by field: the
    flow, composition, pressure, temperature and molar enthalpy the port carries. There is
    no ``Stream`` object, because a model's result has to be a flat set of named fields on
    both sides of the language boundary - the same shape every other model's result has.

    An inlet carries four of the record's five fields and an outlet all five: ``h`` is a
    state function of ``(T, P, z)``, so an inlet's is computed rather than accepted and an
    outlet's is reported because a machine changes it.
    """

    #: Molar flow out, which is the inlet's.
    outlet_n: Q
    #: Outlet composition, one entry per component.
    outlet_z: tuple[float, ...]
    #: Outlet pressure.
    outlet_p: Q
    #: Outlet temperature, solved from the shifted enthalpy.
    outlet_t: Q
    #: Outlet molar enthalpy: the inlet's plus the isentropic head over the efficiency.
    outlet_h: Q
    #: Caveats.
    warnings: tuple[Warning, ...]


@dataclass(frozen=True, slots=True, eq=False)
class ThrottlingValveResult(_HasWarnings):
    """Result of ``process.throttling_valve``.

    The same five fields as :class:`PumpResult`, and for the same reason: a two-port unit
    operation's result is its outlet's record. What differs is the physics - a pump adds
    work and this adds nothing, so the enthalpy crosses unchanged and the temperature is
    what the outlet pressure makes of it.
    """

    #: Outlet molar flow, which is the inlet's.
    outlet_n: Q
    #: Outlet composition, which is the inlet's.
    outlet_z: tuple[float, ...]
    #: Outlet pressure, which is the parameter the valve drops the stream to.
    outlet_p: Q
    #: Outlet temperature, solved from the unchanged enthalpy at the outlet pressure.
    outlet_t: Q
    #: Outlet molar enthalpy, which is the inlet's.
    outlet_h: Q
    #: Caveats.
    warnings: tuple[Warning, ...]


@dataclass(frozen=True, slots=True, eq=False)
class SeparatorResult(_HasWarnings):
    """Result of ``process.separator``.

    **Two outlets of single multiplicity, so ten fields.** The port rule gives a `one`
    port the record's five fields as five result fields, and a separator has two such
    ports, so this is the pair written out under the ports' own names - which is what
    ``unit_ops.separator``'s ``vapour`` and ``liquid`` are called, and not a list of two
    records.
    """

    #: Vapour outlet molar flow.
    vapour_n: Q
    #: Vapour outlet composition.
    vapour_z: tuple[float, ...]
    #: Vapour outlet pressure.
    vapour_p: Q
    #: Vapour outlet temperature.
    vapour_t: Q
    #: Vapour outlet molar enthalpy.
    vapour_h: Q
    #: Liquid outlet molar flow.
    liquid_n: Q
    #: Liquid outlet composition.
    liquid_z: tuple[float, ...]
    #: Liquid outlet pressure.
    liquid_p: Q
    #: Liquid outlet temperature.
    liquid_t: Q
    #: Liquid outlet molar enthalpy.
    liquid_h: Q
    #: Caveats.
    warnings: tuple[Warning, ...]


@dataclass(frozen=True, slots=True, eq=False)
class HeatExchangerResult(_HasWarnings):
    """Result of ``process.heat_exchanger``.

    Two single-multiplicity ports, so ten fields under this unit operation's own port
    names - the shape :class:`SeparatorResult` has, with the two sides carrying different
    fluids rather than two phases of one.
    """

    #: Hot outlet molar flow.
    hot_out_n: Q
    #: Hot outlet composition.
    hot_out_z: tuple[float, ...]
    #: Hot outlet pressure.
    hot_out_p: Q
    #: Hot outlet temperature.
    hot_out_t: Q
    #: Hot outlet molar enthalpy.
    hot_out_h: Q
    #: Cold outlet molar flow.
    cold_out_n: Q
    #: Cold outlet composition.
    cold_out_z: tuple[float, ...]
    #: Cold outlet pressure.
    cold_out_p: Q
    #: Cold outlet temperature.
    cold_out_t: Q
    #: Cold outlet molar enthalpy.
    cold_out_h: Q
    #: Caveats.
    warnings: tuple[Warning, ...]


@dataclass(frozen=True, slots=True, eq=False)
class MixerResult(_HasWarnings):
    """Result of ``process.mixer``.

    **A ``many`` *inlet* port's counterpart to :class:`SplitterResult`**, and the two
    together are the whole of what the process layer's port rule says: the feeds cross as
    vectors and a matrix on the way in, the outlet as the record's five fields on the way
    out. A feed carries four of the five - ``h`` is what the flash computes from the other
    three - because an inlet has no enthalpy to report independently of its state.
    """

    #: Molar flow out: the feeds' sum.
    product_n: Q
    #: Outlet composition, the feeds' weighted by their flows.
    product_z: tuple[float, ...]
    #: Outlet pressure: ``outlet_pressure`` when given, else the lowest feed pressure.
    product_p: Q
    #: Outlet temperature, solved from the joined enthalpy at the outlet pressure.
    product_t: Q
    #: Outlet molar enthalpy, the flow-weighted average of the feeds'.
    product_h: Q
    #: Caveats.
    warnings: tuple[Warning, ...]


@dataclass(frozen=True, slots=True, eq=False)
class SplitterResult(_HasWarnings):
    """Result of ``process.splitter``.

    **A ``many`` port's result, which is the shape the process layer's port rule gives
    it**: one field per record field, each carrying one entry per outlet, and ``z`` a
    matrix with one row per outlet. There is no list of ``Stream`` objects, because a
    result has to be a flat set of named fields on both sides of the language boundary.

    A splitter changes no state, so every entry of every field but ``products_n`` is the
    feed's, written out once per outlet.
    """

    #: Molar flow of each outlet.
    products_n: tuple[Q, ...]
    #: Composition of each outlet, one row per outlet.
    products_z: tuple[tuple[float, ...], ...]
    #: Pressure of each outlet.
    products_p: tuple[Q, ...]
    #: Temperature of each outlet.
    products_t: tuple[Q, ...]
    #: Molar enthalpy of each outlet.
    products_h: tuple[Q, ...]
    #: Caveats.
    warnings: tuple[Warning, ...]


@dataclass(frozen=True, slots=True, eq=False)
class ReactivePhaseEquilibriumResult(_HasWarnings):
    """Result of ``reactions.reactive_phase_equilibrium``.

    **A skip is an answer and not a failure.** ``skipped`` is NeqSim's
    ``getReactivePhaseIndex`` returning ``-1``: the phase is neither aqueous nor liquid
    nor oil, so there is no phase for a water-based equilibrium to be solved in, and the
    composition comes back as it went in. A solve that ran and did not converge is the
    other case, and ``skipped`` is what separates them.
    """

    #: Whether the phase was one the solve runs in.
    skipped: bool
    #: The element matrix: one row per element the components carry, sorted by symbol,
    #: with the electroneutrality row last. Built whether or not the solve ran.
    a_matrix: tuple[tuple[float, ...], ...]
    #: The element amounts the solve conserves, one per row. **The last entry is the
    #: charge correction** and is zero only when every ion in the phase is reactive.
    b: tuple[Q, ...]
    #: Each component's standard-state reference potential, from the independent basis.
    #: Built before the phase is classified, so a skip still carries them.
    chem_ref: tuple[Q, ...]
    #: The composition the phase is left holding.
    moles: tuple[Q, ...]
    #: The solve's passes. Zero when skipped.
    iterations: int
    #: The solve's final error. Zero when skipped.
    error: float
    #: **``solveChemEq``'s return**: the solve converged *and* the three residuals came in
    #: under their tolerances. False when skipped, false for a solve that ran and did not
    #: converge, and **false for one that converged and was not certified**.
    converged: bool
    #: NeqSim's refinement loop: how many ran. One wherever a solve ran, zero on a skip -
    #: the second needs ``useAdaptiveDerivatives`` and a live phase.
    refinements: int
    #: Whether all three residuals came in under their tolerances: ``2e-6`` on the reaction
    #: log residual, ``1e-8`` mol on the net charge, ``1e-8`` mol on the element balance.
    certified: bool
    #: ``max |ln Q - ln K|`` over the reactions the fluid can run, at the answer.
    max_reaction_log_residual: float
    #: The phase's net charge over **every** component it holds. ``nan`` on a skip, which is
    #: what NeqSim returns there.
    net_charge_moles: Q
    #: ``max |A n - b|`` over the element rows, at the answer - the charge row excluded.
    max_element_residual: Q
    #: **Whether the linear program's estimate became the starting composition.** False is
    #: a state and not a failure: NeqSim reads the estimate's absence as keeping the phase's
    #: own composition, which is what ``seed_moles`` reports then.
    seed_applied: bool
    #: The composition the solve started from: the estimate where one exists, the caller's
    #: own ``moles`` where none does.
    seed_moles: tuple[Q, ...]
    #: Caveats.
    warnings: tuple[Warning, ...]


@dataclass(frozen=True, slots=True, eq=True)
class ReactivePhFlashResult(_HasWarnings):
    """Result of ``reactions.reactive_ph_flash``.

    The temperature at which a reactive fluid's thermochemical enthalpy matches a
    specification, and the cost of finding it. **The pass counts are path quantities**: the
    loop is a secant, so two implementations - and NeqSim itself - reach the same temperature
    in a different number of steps.
    """

    #: ``getEquilibriumTemperature``: the temperature the loop stopped at.
    temperature: Q
    #: ``isConverged``. **Also true where the *bracket* closed rather than the residual**,
    #: which the flag does not distinguish.
    converged: bool
    #: ``getOuterIterations``: the temperature steps taken.
    outer_iterations: int
    #: ``getTotalInnerIterations``: every inner flash's passes, summed.
    total_inner_iterations: int
    #: Caveats.
    warnings: tuple[Warning, ...]


@dataclass(frozen=True, slots=True, eq=False)
class ReactiveTpFlashResult(_HasWarnings):
    """Result of ``reactions.reactive_tp_flash``.

    **The composition is the answer and the split is not.** Where two phases converge to
    the *same* composition the Gibbs energy is flat along the direction that trades moles
    between them, so every split satisfies the equilibrium conditions and the one a run
    reports is decided by its path. The rows below are that run's own; the moles summed
    over them are what the element balance and the equilibrium fix.
    """

    #: How many phases the driver stopped on.
    phase_count: int
    #: Each phase's mole numbers, one row per phase and one column per component, in the
    #: driver's own order. **No phase type is promised**: NeqSim's types are its system's
    #: bookkeeping - `gas`, `oil`, `aqueous` - and not a state this model computes.
    phase_moles: tuple[tuple[Q, ...], ...]
    #: Each phase's share of the fluid, and **the fraction the driver weighs it by**: on a
    #: neutral fluid that is the one its bookkeeping left rather than the solve's own.
    phase_fraction: tuple[float, ...]
    #: ``isConverged``. **True on two branches that accept an answer without a converged
    #: solve** - the single-phase one and the ``NR = 0`` fallback - which is the class's
    #: own overwriting of the solve's flag.
    converged: bool
    #: Every solve's passes summed over the outer iterations. **Zero on the ``NR = 0``
    #: fallback**, whose successive substitution the driver never counts.
    total_iterations: int
    #: ``getEquilibriumTotalMoles``, the last solve's total. One mole of feed does not fix
    #: it: a reaction that splits one species into two moves it.
    equilibrium_total_moles: Q
    #: ``computeGibbsEnergy``, the beta-weighted ``sum_i x_i (ln x_i + ln phi_i)``. **Zero
    #: on the single-phase branch**, which returns before the driver computes one, and
    #: twice one phase's worth where the weights are the constructor's ``(1.0, 1.0)``.
    gibbs_energy: float
    #: ``getFinalResidual``: what the solve stopped on. A relaxed multiphase stop is
    #: ``1e-4`` rather than ``1e-9``, and the answer is only as tight as this.
    residual: float
    #: ``getFinalElementResidual`` on its own.
    element_residual: float
    #: Caveats.
    warnings: tuple[Warning, ...]


@dataclass(frozen=True, slots=True, eq=False)
class ReferencePotentialsResult(_HasWarnings):
    """Result of ``reactions.reference_potentials``.

    The first result carrying a *mask* rather than a list of names: which components the
    basis solved for and which reactions survived are positions in the caller's component
    order and the table's row order, and a position is what the basis actually chose over.
    """

    #: The standard-state reference potentials, one per component, in the caller's order.
    #: A potential solved for directly and one propagated are the same quantity;
    #: ``independent`` says which is which.
    potentials: tuple[Q, ...]
    #: A mask over the components: 1.0 where the basis solved directly, 0.0 where the
    #: potential was propagated from the stoichiometry.
    independent: tuple[float, ...]
    #: A mask over the source's **loaded** reactions - the rows whose ``USEREACTION`` is
    #: 1, in the table's physical order. That order is part of the answer, because the
    #: reducer is greedy over it.
    survivors: tuple[float, ...]
    #: The rank the reaction basis reached.
    rank: int
    #: Caveats.
    warnings: tuple[Warning, ...]


@dataclass(frozen=True, slots=True, eq=False)
class KineticRateLawResult(_HasWarnings):
    """Result of ``reactions.kinetic_rate_law``."""

    #: The reaction's rate factor at ``T``, by the selected law - a bare number, as the
    #: class stores and returns it.
    rate_factor: float
    #: Caveats.
    warnings: tuple[Warning, ...]


@dataclass(frozen=True, slots=True, eq=False)
class KineticsResult(_HasWarnings):
    """Result of ``reactions.kinetics``.

    One entry per component of the phase, in the order it was given.
    """

    #: ``reacCoef`` per component: the pseudo-first-order coefficient, the sum of every
    #: reaction's own contribution.
    coefficient: tuple[float, ...]
    #: ``getPhiInfinite`` per component, **zero where no reaction produced one** - the
    #: value the class's field is constructed with, which is what a fresh object reads
    #: back.
    phi_infinite: tuple[float, ...]
    #: A mask: 1.0 where any reaction's scaled ``1/K`` came in under ``1e-3`` while this
    #: component's row was built.
    irreversible: tuple[float, ...]
    #: Caveats.
    warnings: tuple[Warning, ...]


@dataclass(frozen=True, slots=True, eq=False)
class PureSaturationResult(_HasWarnings):
    """Result of ``eos.pure_saturation``.

    The first *model* result: its spec fixes a procedure rather than an equation, and
    it composes the kernels rather than reimplementing them. The shape is an ordinary
    result shape, deliberately - a model's answer is an answer like any other, and
    the difference is in how it was reached.
    """

    #: The saturation pressure.
    p_sat: Q
    #: The common value of ``ln phi_L`` and ``ln phi_V`` at the converged pressure.
    ln_phi: float
    #: Bisection steps taken.
    iterations: int
    #: The dimensionless half-width of the final bracket - the relative uncertainty in
    #: the reduced pressure, not the fugacity residual.
    residual: float
    #: Caveats.
    warnings: tuple[Warning, ...]


@dataclass(frozen=True, slots=True, eq=False)
class FreezingPointResult(_HasWarnings):
    """Result of ``eos.freezing_point``."""

    #: The temperature at which the solid and the fluid meet.
    temperature: Q
    #: Which substance's freezing point this is, of the fluid's candidates.
    component: str
    #: Bracket expansions and bisection steps together, as NeqSim counts them.
    iterations: int
    #: The dimensionless residual at the reported temperature: the Gibbs difference on the
    #: Helmholtz route and the multiphase appearance condition on the tabulated one.
    residual: float
    #: Caveats.
    warnings: tuple[Warning, ...]


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
class HydrateFormationTemperatureResult(_HasWarnings):
    """Result of ``eos.hydrate_formation_temperature``."""

    #: The temperature at which the hydrate's water fugacity meets the fluid's.
    temperature: Q
    #: The structure that comparison found stable.
    structure: HydrateStructure
    #: Flash evaluations taken, the scan's and the bisection's together.
    iterations: int
    #: ``f_w^hydrate / f_w^fluid - 1`` at the reported temperature.
    residual: float
    #: Caveats.
    warnings: tuple[Warning, ...]


@dataclass(frozen=True, slots=True, eq=False)
class TpSolidFlashResult(_HasWarnings):
    """Result of ``eos.tp_solid_flash``."""

    #: The fraction of the feed's moles in the pure solid, and zero where none forms.
    solid_fraction: float
    #: How many phases the feed splits into, the solid counted when it is there.
    phase_count: int
    #: The mole fraction of the feed in each phase, summing to one, the solid last.
    beta: tuple[float, ...]
    #: The composition of each phase, one tuple per phase.
    x: tuple[tuple[float, ...], ...]
    #: The pure solid's fugacity coefficient, the number the phase exists at all is answered
    #: from.
    solid_fugacity_coefficient: float
    #: Fraction-solve steps taken by the last solve.
    iterations: int
    #: The norm of the last fraction correction.
    residual: float
    #: Whether the residual met the tolerance.
    converged: bool
    #: Caveats.
    warnings: tuple[Warning, ...]


@dataclass(frozen=True, slots=True, eq=False)
class TpMultiflashWaxResult(_HasWarnings):
    """Result of ``eos.tp_multiflash_wax``."""

    #: The fraction of the feed's moles in the wax phase, and zero where none forms.
    wax_fraction: float
    #: How many phases the feed splits into: two, or three where the wax survives.
    phase_count: int
    #: The mole fraction of the feed in each phase, summing to one, the wax last.
    beta: tuple[float, ...]
    #: The composition of each phase, one tuple per phase.
    x: tuple[tuple[float, ...], ...]
    #: Fraction-solve steps taken.
    iterations: int
    #: The norm of the last fraction correction.
    residual: float
    #: Whether the residual met the tolerance. **False means the two-phase answer was
    #: returned**, so this is the flag that says which of the two was solved.
    converged: bool
    #: Caveats.
    warnings: tuple[Warning, ...]


@dataclass(frozen=True, slots=True, eq=False)
class WaxSolidFugacityResult(_HasWarnings):
    """Result of ``eos.wax_solid_fugacity``."""

    #: The wax phase's fugacity coefficient for one component. Dimensionless.
    fugacity_coefficient: float
    #: Caveats.
    warnings: tuple[Warning, ...]


@dataclass(frozen=True, slots=True, eq=False)
class SolidFugacityResult(_HasWarnings):
    """Result of ``eos.solid_fugacity``."""

    #: The solid's fugacity coefficient for one component. Dimensionless.
    fugacity_coefficient: float
    #: Caveats.
    warnings: tuple[Warning, ...]


@dataclass(frozen=True, slots=True, eq=False)
class SaltPrecipitationResult(_HasWarnings):
    """Result of ``eos.salt_precipitation``."""

    #: The solid taken, in moles per mole of feed.
    precipitated_moles: float
    #: ``IAP/Ksp`` before anything was taken.
    initial_saturation_ratio: float
    #: ``IAP/Ksp`` at the answer: one where the mineral precipitated, and the initial ratio
    #: where it did not.
    final_saturation_ratio: float
    #: Bisection steps taken.
    iterations: int
    #: The extent reached against the bracket's upper end: one where an ion was exhausted.
    extent_of_maximum: float
    #: Caveats.
    warnings: tuple[Warning, ...]


@dataclass(frozen=True, slots=True, eq=False)
class ScaleSaturationRatioResult(_HasWarnings):
    """Result of ``eos.scale_saturation_ratio``."""

    #: ``IAP/Ksp``, one at saturation, clamped at ``exp(69)`` and floored at zero.
    saturation_ratio: float
    #: The ion activity product, so a divergence can be read as either half.
    ion_activity_product: float
    #: ``Ksp(T, P)`` after the override, the hydrogen-ion correction and the pressure term.
    solubility_product: float
    #: Caveats.
    warnings: tuple[Warning, ...]


@dataclass(frozen=True, slots=True, eq=False)
class TbpFractionPropertiesResult(_HasWarnings):
    """Result of ``eos.tbp_fraction_properties``."""

    #: Critical temperature.
    tc: Q
    #: Critical pressure.
    pc: Q
    #: Normal boiling point.
    boiling_temperature: Q
    #: Pitzer's acentric factor.
    acentric_factor: float
    #: The ``m`` of the cubic's alpha function, which replaces the acentric-factor
    #: correlation for a pseudo-component.
    attraction_exponent: float
    #: Caveats.
    warnings: tuple[Warning, ...]


@dataclass(frozen=True, slots=True, eq=False)
class HydrateFormationPressureResult(_HasWarnings):
    """Result of ``eos.hydrate_formation_pressure``."""

    #: The pressure at which the hydrate's water fugacity meets the fluid's.
    pressure: Q
    #: The structure that comparison found stable.
    structure: HydrateStructure
    #: Flash evaluations taken, the scan's and the bisection's together.
    iterations: int
    #: ``f_w^hydrate / f_w^fluid - 1`` at the reported pressure.
    residual: float
    #: Caveats.
    warnings: tuple[Warning, ...]


@dataclass(frozen=True, slots=True, eq=False)
class HydrateEquilibriumLineResult(_HasWarnings):
    """Result of ``eos.hydrate_equilibrium_line``.

    A hydrate curve: the formation temperature at each of ten equally spaced pressures from a
    minimum to a maximum. The grid is NeqSim's own - ten points, whatever the bounds - and the
    arithmetic at each point is ``eos.hydrate_formation_temperature``'s.
    """

    #: The formation temperatures, in grid order, in kelvin.
    temperature: tuple[float, ...]
    #: The pressures the points were solved at, in pascals, parallel to ``temperature``.
    pressure: tuple[float, ...]
    #: Caveats, from every point of the grid.
    warnings: tuple[Warning, ...]


@dataclass(frozen=True, slots=True, eq=False)
class HydrateInhibitorWtResult(_HasWarnings):
    """Result of ``eos.hydrate_inhibitor_wt``."""

    #: The inhibitor's moles at the answer, the feed's own included.
    inhibitor_moles: float
    #: The **aqueous phase's** inhibitor mass fraction, which is what the secant drives to the
    #: target. Not the feed's own: a gas phase takes some inhibitor with it.
    weight_fraction: float
    #: How many phases the answer's state has.
    phases: int
    #: Secant steps taken, with a floor of three.
    iterations: int
    #: ``-(wtp - wt_target)`` at the answer.
    residual: float
    #: Caveats.
    warnings: tuple[Warning, ...]


@dataclass(frozen=True, slots=True, eq=False)
class HydrateInhibitorConcentrationResult(_HasWarnings):
    """Result of ``eos.hydrate_inhibitor_concentration``."""

    #: The inhibitor's moles at the answer, the feed's own included.
    inhibitor_moles: float
    #: The inhibitor's mass fraction of the inhibitor-and-water pair, which is what a dosing
    #: figure reports - the hydrocarbons are not in its denominator.
    weight_fraction: float
    #: The hydrate temperature the answer's composition gives, which should be the target to
    #: the tolerance.
    hydrate_temperature: Q
    #: Secant steps taken, with a floor of three.
    iterations: int
    #: ``T_hydrate - T_target`` at the answer, in kelvin.
    residual: float
    #: Caveats.
    warnings: tuple[Warning, ...]


@dataclass(frozen=True, slots=True, eq=False)
class HydrateFractionResult(_HasWarnings):
    """Result of ``eos.hydrate_fraction``."""

    #: The fraction of the feed's moles that is hydrate, on the feed's own basis.
    beta: float
    #: The structure the cages at the answer are.
    structure: HydrateStructure
    #: The largest ``|sum_p beta_p x_ip - z_i|`` over the components at the answer.
    #:
    #: The invariant this model exists to keep: zero here by construction, and ``7.3e-13``
    #: on NeqSim's own state at the same feed.
    balance_error: float
    #: Flash evaluations taken: the feed's, the fixed point's and the answer's own.
    iterations: int
    #: ``ln(f_w^hydrate/f_w^fluid)`` at the **feed**, which is what decides between no
    #: hydrate and all of the water. Not read at the answer: the fluid there holds no
    #: water, so the ratio has no finite value at it.
    residual: float
    #: Caveats.
    warnings: tuple[Warning, ...]


@dataclass(frozen=True, slots=True, eq=False)
class PrMolarVolumeResult(_HasWarnings):
    """Result of ``eos.pr_molar_volume``.

    The only result in this namespace carrying a ``pint`` quantity. Everything else
    here is dimensionless, so this is where the boundary rule starts doing work
    again.
    """

    #: Molar volume.
    v: Q
    #: Caveats.
    warnings: tuple[Warning, ...]


@dataclass(frozen=True, slots=True, eq=False)
class PrMassDensityResult(_HasWarnings):
    """Result of ``eos.pr_mass_density``."""

    #: Mass density. Positive whenever the molar mass and volume are, which the
    #: bounds ensure.
    rho: Q
    #: Caveats.
    warnings: tuple[Warning, ...]


@dataclass(frozen=True, slots=True, eq=False)
class PrPenelouxShiftResult(_HasWarnings):
    """Result of ``eos.pr_peneloux_shift``."""

    #: The Peneloux volume-translation parameter.
    c: Q
    #: Caveats.
    warnings: tuple[Warning, ...]


@dataclass(frozen=True, slots=True, eq=False)
class SrkPenelouxShiftResult(_HasWarnings):
    """Result of ``eos.srk_peneloux_shift``."""

    #: The Peneloux volume-translation parameter.
    c: Q
    #: Caveats.
    warnings: tuple[Warning, ...]


@dataclass(frozen=True, slots=True, eq=False)
class HeatOfVaporizationResult(_HasWarnings):
    """Result of ``eos.heat_of_vaporization``."""

    #: The pure-component heat of vaporisation.
    hov: Q
    #: Caveats.
    warnings: tuple[Warning, ...]


@dataclass(frozen=True, slots=True, eq=False)
class LiquidHeatCapacityResult(_HasWarnings):
    """Result of ``eos.liquid_heat_capacity``."""

    #: The pure-component liquid heat capacity.
    cp: Q
    #: Caveats.
    warnings: tuple[Warning, ...]


@dataclass(frozen=True, slots=True, eq=False)
class AntoineVaporPressureResult(_HasWarnings):
    """Result of ``eos.antoine_vapor_pressure``."""

    #: The pure-component vapour pressure.
    p_sat: Q
    #: Caveats.
    warnings: tuple[Warning, ...]


@dataclass(frozen=True, slots=True, eq=False)
class NitricSulfuricAcidVaporPressureResult(_HasWarnings):
    """Result of ``eos.nitric_sulfuric_acid_vapor_pressure``.

    Three pressures from one temperature, because a phase over this system needs all
    three at one state and three calculations would be three chances to evaluate them at
    different ones. Each comes from its own correlation, and the nitric-acid one is
    NeqSim's refit of Pennington's Antoine pair rather than the paper's.
    """

    #: The pure-component saturation vapour pressure of water.
    p_water: Q
    #: The pure-component saturation vapour pressure of nitric acid.
    p_nitric_acid: Q
    #: The pure-component saturation vapour pressure of sulfuric acid.
    p_sulfuric_acid: Q
    #: Caveats.
    warnings: tuple[Warning, ...]


@dataclass(frozen=True, slots=True, eq=False)
class RackettMolarVolumeResult(_HasWarnings):
    """Result of ``eos.rackett_molar_volume``."""

    #: The saturated liquid molar volume.
    v: Q
    #: Caveats.
    warnings: tuple[Warning, ...]


@dataclass(frozen=True, slots=True, eq=False)
class CostaldMolarVolumeResult(_HasWarnings):
    """Result of ``eos.costald_molar_volume``."""

    #: The saturated liquid molar volume.
    v: Q
    #: Caveats.
    warnings: tuple[Warning, ...]


@dataclass(frozen=True, slots=True, eq=False)
class ChungViscosityResult(_HasWarnings):
    """Result of ``eos.chung_viscosity``."""

    #: The gas dynamic viscosity.
    mu: Q
    #: Caveats.
    warnings: tuple[Warning, ...]


@dataclass(frozen=True, slots=True, eq=False)
class ChungConductivityResult(_HasWarnings):
    """Result of ``eos.chung_conductivity``."""

    #: The gas thermal conductivity.
    k: Q
    #: Caveats.
    warnings: tuple[Warning, ...]


@dataclass(frozen=True, slots=True, eq=False)
class MasonSaxenaConductivityResult(_HasWarnings):
    """Result of ``eos.mason_saxena_conductivity``."""

    #: The gas mixture thermal conductivity.
    k: Q
    #: Caveats.
    warnings: tuple[Warning, ...]


@dataclass(frozen=True, slots=True, eq=False)
class TynCalusDiffusivityResult(_HasWarnings):
    """Result of ``eos.tyn_calus_diffusivity``."""

    #: The binary diffusion coefficient.
    d: Q
    #: Caveats.
    warnings: tuple[Warning, ...]


@dataclass(frozen=True, slots=True, eq=False)
class UmrprAlphaResult(_HasWarnings):
    """Result of ``eos.umrpr_alpha``."""

    #: The temperature-dependent alpha function. Dimensionless.
    alpha: float
    #: Caveats.
    warnings: tuple[Warning, ...]


@dataclass(frozen=True, slots=True, eq=False)
class WilkeChangDiffusivityResult(_HasWarnings):
    """Result of ``eos.wilke_chang_diffusivity``."""

    #: The binary diffusion coefficient.
    d: Q
    #: Caveats.
    warnings: tuple[Warning, ...]


@dataclass(frozen=True, slots=True, eq=False)
class HaydukMinhasDiffusivityResult(_HasWarnings):
    """Result of ``eos.hayduk_minhas_diffusivity``."""

    #: The binary diffusion coefficient.
    d: Q
    #: Caveats.
    warnings: tuple[Warning, ...]


@dataclass(frozen=True, slots=True, eq=False)
class SchwartzentruberAlphaResult(_HasWarnings):
    """Result of ``eos.schwartzentruber_alpha``."""

    #: The temperature-dependent alpha function. Dimensionless.
    alpha: float
    #: Caveats.
    warnings: tuple[Warning, ...]


@dataclass(frozen=True, slots=True, eq=False)
class SoreideWhitsonAlphaResult(_HasWarnings):
    """Result of ``eos.soreide_whitson_alpha``."""

    #: The temperature-dependent alpha function. Dimensionless.
    alpha: float
    #: Caveats.
    warnings: tuple[Warning, ...]


@dataclass(frozen=True, slots=True, eq=False)
class SiddiqiLucasDiffusivityResult(_HasWarnings):
    """Result of ``eos.siddiqi_lucas_diffusivity``."""

    #: The binary diffusion coefficient.
    d: Q
    #: Caveats.
    warnings: tuple[Warning, ...]


@dataclass(frozen=True, slots=True, eq=False)
class Co2WaterDiffusivityResult(_HasWarnings):
    """Result of ``eos.co2_water_diffusivity``."""

    #: The binary diffusion coefficient.
    d: Q
    #: Caveats.
    warnings: tuple[Warning, ...]


@dataclass(frozen=True, slots=True, eq=False)
class ParachorSurfaceTensionResult(_HasWarnings):
    """Result of ``eos.parachor_surface_tension``."""

    #: The surface tension, in N/m.
    sigma: Q
    #: Caveats.
    warnings: tuple[Warning, ...]


@dataclass(frozen=True, slots=True, eq=False)
class ViscosityResult(_HasWarnings):
    """Result of ``eos.viscosity``."""

    #: The liquid dynamic viscosity, in Pa*s.
    mu: Q
    #: Caveats.
    warnings: tuple[Warning, ...]


@dataclass(frozen=True, slots=True, eq=False)
class ThermalConductivityResult(_HasWarnings):
    """Result of ``eos.thermal_conductivity``."""

    #: The liquid thermal conductivity, in W/(m*K).
    k: Q
    #: Caveats.
    warnings: tuple[Warning, ...]


@dataclass(frozen=True, slots=True, eq=False)
class Matcop5PrumrAlphaResult(_HasWarnings):
    """Result of ``eos.matcop5_prumr_alpha``."""

    #: The temperature-dependent alpha function. Dimensionless.
    alpha: float
    #: Caveats.
    warnings: tuple[Warning, ...]


@dataclass(frozen=True, slots=True, eq=False)
class MatcopAlphaResult(_HasWarnings):
    """Result of ``eos.matcop_alpha``."""

    #: The temperature-dependent alpha function. Dimensionless.
    alpha: float
    #: Caveats.
    warnings: tuple[Warning, ...]


@dataclass(frozen=True, slots=True, eq=False)
class MatcopPrAlphaResult(_HasWarnings):
    """Result of ``eos.matcop_pr_alpha``."""

    #: The temperature-dependent alpha function. Dimensionless.
    alpha: float
    #: Caveats.
    warnings: tuple[Warning, ...]


@dataclass(frozen=True, slots=True, eq=False)
class MatcopPrumrAlphaResult(_HasWarnings):
    """Result of ``eos.matcop_prumr_alpha``."""

    #: The temperature-dependent alpha function. Dimensionless.
    alpha: float
    #: Caveats.
    warnings: tuple[Warning, ...]


@dataclass(frozen=True, slots=True, eq=False)
class MatcopPrumrNewAlphaResult(_HasWarnings):
    """Result of ``eos.matcop_prumr_new_alpha``."""

    #: The temperature-dependent alpha function. Dimensionless.
    alpha: float
    #: Caveats.
    warnings: tuple[Warning, ...]


@dataclass(frozen=True, slots=True, eq=False)
class WilkeViscosityResult(_HasWarnings):
    """Result of ``eos.wilke_viscosity``."""

    #: The gas mixture dynamic viscosity.
    mu: Q
    #: Caveats.
    warnings: tuple[Warning, ...]


@dataclass(frozen=True, slots=True, eq=False)
class NrtlActivityCoefficientsResult(_HasWarnings):
    """Result of ``eos.nrtl_activity_coefficients``."""

    #: The natural logarithm of each activity coefficient.
    ln_gamma: tuple[float, ...]
    #: The activity coefficient of each component.
    gamma: tuple[float, ...]
    #: Caveats.
    warnings: tuple[Warning, ...]


@dataclass(frozen=True, slots=True, eq=False)
class UnifacActivityCoefficientsResult(_HasWarnings):
    """Result of ``eos.unifac_activity_coefficients``."""

    #: The natural logarithm of each activity coefficient.
    ln_gamma: tuple[float, ...]
    #: The activity coefficient of each component.
    gamma: tuple[float, ...]
    #: Caveats.
    warnings: tuple[Warning, ...]


@dataclass(frozen=True, slots=True, eq=False)
class UnifacUmrpruActivityCoefficientsResult(_HasWarnings):
    """Result of ``eos.unifac_umrpru_activity_coefficients``."""

    #: The natural logarithm of each activity coefficient.
    ln_gamma: tuple[float, ...]
    #: The activity coefficient of each component.
    gamma: tuple[float, ...]
    #: Caveats.
    warnings: tuple[Warning, ...]


@dataclass(frozen=True, slots=True, eq=False)
class UnifacPsrkActivityCoefficientsResult(_HasWarnings):
    """Result of ``eos.unifac_psrk_activity_coefficients``."""

    #: The natural logarithm of each activity coefficient.
    ln_gamma: tuple[float, ...]
    #: The activity coefficient of each component.
    gamma: tuple[float, ...]
    #: Caveats.
    warnings: tuple[Warning, ...]


@dataclass(frozen=True, slots=True, eq=False)
class VanLaarAcidActivityCoefficientsResult(_HasWarnings):
    """Result of ``eos.van_laar_acid_activity_coefficients``."""

    #: The natural logarithm of each activity coefficient.
    ln_gamma: tuple[float, ...]
    #: The activity coefficient of each component.
    gamma: tuple[float, ...]
    #: Caveats.
    warnings: tuple[Warning, ...]


@dataclass(frozen=True, slots=True, eq=False)
class UniquacActivityCoefficientsResult(_HasWarnings):
    """Result of ``eos.uniquac_activity_coefficients``."""

    #: The natural logarithm of each activity coefficient.
    ln_gamma: tuple[float, ...]
    #: The activity coefficient of each component.
    gamma: tuple[float, ...]
    #: Caveats.
    warnings: tuple[Warning, ...]


@dataclass(frozen=True, slots=True, eq=False)
class WilsonActivityCoefficientsResult(_HasWarnings):
    """Result of ``eos.wilson_activity_coefficients``."""

    #: The natural logarithm of each activity coefficient.
    ln_gamma: tuple[float, ...]
    #: The activity coefficient of each component.
    gamma: tuple[float, ...]
    #: Caveats.
    warnings: tuple[Warning, ...]


@dataclass(frozen=True, slots=True, eq=False)
class VuFlashSingleCompResult(_HasWarnings):
    """Result of ``eos.vu_flash_single_comp``."""

    #: The saturation temperature at the pressure asked for - the state's temperature.
    T: Q
    #: The vapour fraction, from the lever rule on the two saturated internal energies.
    beta: float
    #: The molar volume the split implies.
    V: Q
    #: Always ``TWO_PHASE``; the other values are reachable only as a refusal's diagnosis.
    phase: Phase
    #: Caveats.
    warnings: tuple[Warning, ...]


@dataclass(frozen=True, slots=True, eq=False)
class PvfFlashResult(_HasWarnings):
    """Result of ``eos.pvf_flash``."""

    #: The temperature at which the feed's vapour fraction is the one asked for.
    T: Q
    #: The vapour fraction at the answer, as the iteration measured it.
    beta: float
    #: Which phase the feed is in at the answer.
    phase: Phase
    #: Liquid-phase composition at the answer.
    x: tuple[float, ...]
    #: Vapour-phase composition.
    y: tuple[float, ...]
    #: K-values at the answer.
    k: tuple[float, ...]
    #: Liquid root of the cubic at the answer.
    z_liquid: float
    #: Vapour root.
    z_vapour: float
    #: Illinois steps taken.
    iterations: int
    #: ``|beta(T) - beta_spec|`` at the answer.
    residual: float
    #: Caveats.
    warnings: tuple[Warning, ...]


@dataclass(frozen=True, slots=True, eq=False)
class VsFlashResult(_HasWarnings):
    """Result of ``eos.vs_flash``.

    A closed vessel at fixed volume and entropy: both the pressure and the
    temperature are answers, reported alongside the phase split.
    """

    #: The pressure that satisfies the volume and entropy.
    P: Q
    #: The temperature that satisfies the volume and entropy.
    T: Q
    #: The vapour fraction at the answer, or ``None`` for a single-phase feed.
    beta: float | None
    #: Liquid-phase composition at the answer.
    x: tuple[float, ...]
    #: Vapour-phase composition.
    y: tuple[float, ...]
    #: K-values at the answer.
    k: tuple[float, ...]
    #: Which phase the feed is in at the answer.
    phase: Phase
    #: Liquid root of the cubic at the answer.
    z_liquid: float
    #: Vapour root.
    z_vapour: float
    #: Newton steps taken.
    iterations: int
    #: The larger of the relative volume and enthalpy residuals at the answer.
    residual: float
    #: Caveats, deduplicated.
    warnings: tuple[Warning, ...]


@dataclass(frozen=True, slots=True, eq=False)
class TvFractionFlashResult(_HasWarnings):
    """Result of ``eos.tv_fraction_flash``."""

    #: The pressure at which the volume fraction is the one asked for.
    P: Q
    #: The temperature the flash was taken at, echoed.
    T: Q
    #: The vapour fraction at the answer, or ``None`` for a single-phase feed.
    beta: float | None
    #: The gas phase's volume share at the answer.
    volume_fraction: float
    #: Which phase the feed is in at the answer.
    phase: Phase
    #: Liquid-phase composition at the answer.
    x: tuple[float, ...]
    #: Vapour-phase composition.
    y: tuple[float, ...]
    #: K-values at the answer.
    k: tuple[float, ...]
    #: Liquid root of the cubic at the answer.
    z_liquid: float
    #: Vapour root.
    z_vapour: float
    #: Newton steps taken.
    iterations: int
    #: ``|volume_fraction(P) - fraction|`` at the answer.
    residual: float
    #: Caveats.
    warnings: tuple[Warning, ...]


@dataclass(frozen=True, slots=True, eq=False)
class VuFlashResult(_HasWarnings):
    """Result of ``eos.vu_flash``.

    A closed vessel at fixed volume and internal energy: both the pressure and the
    temperature are answers, reported alongside the phase split.
    """

    #: The pressure that satisfies the volume and internal energy.
    P: Q
    #: The temperature that satisfies the volume and internal energy.
    T: Q
    #: The vapour fraction at the answer, or ``None`` for a single-phase feed.
    beta: float | None
    #: Liquid-phase composition at the answer.
    x: tuple[float, ...]
    #: Vapour-phase composition.
    y: tuple[float, ...]
    #: K-values at the answer.
    k: tuple[float, ...]
    #: Which phase the feed is in at the answer.
    phase: Phase
    #: Liquid root of the cubic at the answer.
    z_liquid: float
    #: Vapour root.
    z_vapour: float
    #: Newton steps taken.
    iterations: int
    #: The larger of the relative volume and enthalpy residuals at the answer.
    residual: float
    #: Caveats, deduplicated.
    warnings: tuple[Warning, ...]


@dataclass(frozen=True, slots=True, eq=False)
class PvRefluxFlashResult(_HasWarnings):
    """Result of ``eos.pv_reflux_flash``."""

    #: The temperature at which the phase ratio is the one asked for.
    T: Q
    #: The vapour fraction at the answer, or ``None`` for a single-phase feed.
    beta: float | None
    #: Which phase the feed is in at the answer.
    phase: Phase
    #: Liquid-phase composition at the answer.
    x: tuple[float, ...]
    #: Vapour-phase composition.
    y: tuple[float, ...]
    #: K-values at the answer.
    k: tuple[float, ...]
    #: Liquid root of the cubic at the answer.
    z_liquid: float
    #: Vapour root.
    z_vapour: float
    #: Secant steps taken.
    iterations: int
    #: ``|reflux_ratio(T) - reflux_spec|`` at the answer.
    residual: float
    #: Caveats.
    warnings: tuple[Warning, ...]


@dataclass(frozen=True, slots=True, eq=False)
class VhFlashResult(_HasWarnings):
    """Result of ``eos.vh_flash``.

    A closed vessel at fixed volume and enthalpy: both the pressure and the
    temperature are answers, reported alongside the phase split.
    """

    #: The pressure that satisfies the volume and enthalpy.
    P: Q
    #: The temperature that satisfies the volume and enthalpy.
    T: Q
    #: The vapour fraction at the answer, or ``None`` for a single-phase feed.
    beta: float | None
    #: Liquid-phase composition at the answer.
    x: tuple[float, ...]
    #: Vapour-phase composition.
    y: tuple[float, ...]
    #: K-values at the answer.
    k: tuple[float, ...]
    #: Which phase the feed is in at the answer.
    phase: Phase
    #: Liquid root of the cubic at the answer.
    z_liquid: float
    #: Vapour root.
    z_vapour: float
    #: Newton steps taken.
    iterations: int
    #: The larger of the relative volume and enthalpy residuals at the answer.
    residual: float
    #: Caveats, deduplicated.
    warnings: tuple[Warning, ...]


@dataclass(frozen=True, slots=True, eq=False)
class Vdw1fMixBinaryResult(_HasWarnings):
    """Result of ``eos.vdw1f_mix_binary``.

    The mixture's van der Waals one-fluid parameters. Both dimensionless, like
    everything else in this namespace.
    """

    #: The mixture's attraction parameter.
    a_mix: float
    #: The mixture's repulsion parameter - the mole-fraction-weighted mean of the
    #: pure ones, with no binary correction. vdW1f has no ``l12``.
    b_mix: float
    #: Caveats.
    warnings: tuple[Warning, ...]


@dataclass(frozen=True, slots=True, eq=False)
class RachfordRiceResult(_HasWarnings):
    """Result of ``eos.rachford_rice``."""

    #: The root of the Rachford-Rice equation, returned as the equation gives it:
    #: outside ``[0, 1]`` it is the negative flash rather than a phase split, and the
    #: caller decides what that means. A feed whose K-values do not straddle one has no
    #: root, and comes back as NeqSim's ``1e-12`` clamp at the end it lies towards.
    beta: float
    #: Caveats.
    warnings: tuple[Warning, ...]


@dataclass(frozen=True, slots=True, eq=False)
class RachfordRiceBinaryResult(_HasWarnings):
    """Result of ``eos.rachford_rice_binary``."""

    #: The vapour fraction that solves the Rachford-Rice equation. Outside ``[0, 1]``
    #: the feed is single phase and this is the tangent-plane value rather than a
    #: phase split; the result carries ``OUT_OF_VALID_RANGE`` when so. ``beta < 0``
    #: means subcooled liquid and ``beta > 1`` superheated vapour.
    beta: float
    #: Caveats.
    warnings: tuple[Warning, ...]


@dataclass(frozen=True, slots=True, eq=False)
class PrDepartureResult(_HasWarnings):
    """Result of ``eos.pr_departure``.

    Four dimensionless outputs: the logarithm of the fugacity coefficient, and the
    departure enthalpy, entropy and heat capacity made dimensionless as ``h_dep_rt``,
    ``s_dep_r`` and ``cp_dep_r``. The multiplication by ``R`` and ``T`` happens where
    those live - the model layer - so this namespace stays unit-free end to end.
    """

    #: The logarithm of the fugacity coefficient. Returned as a logarithm rather
    #: than as ``phi`` because the logarithm is what the algebra produces, what the
    #: equilibrium condition equates, and what makes ``ln(K) = ln_phi_l - ln_phi_v``
    #: a subtraction rather than a division.
    ln_phi: float
    #: The departure enthalpy over ``R*T`` - the departure from ideal-gas behaviour
    #: at the same temperature and pressure.
    h_dep_rt: float
    #: The departure entropy over ``R``.
    s_dep_r: float
    #: The departure heat capacity over ``R`` - ``Cp`` relative to the ideal-gas
    #: value at the same state. It is ``d(H_dep/RT)/d(ln T)`` at constant pressure,
    #: which is the derivative the isentropic and isenthalpic flashes need.
    cp_dep_r: float
    #: Caveats.
    warnings: tuple[Warning, ...]


@dataclass(frozen=True, slots=True, eq=False)
class PrsvKappaResult(_HasWarnings):
    """Result of ``eos.prsv_kappa``.

    One dimensionless number, like :class:`PrKappaResult` - but a *different* number
    with a different property: this one varies with temperature, where
    Peng-Robinson's coefficient does not. The two are not interchangeable and a
    caller who substitutes one for the other gets a plausible answer from the wrong
    correlation.
    """

    #: The PRSV alpha-function coefficient. Dimensionless.
    kappa: float
    #: Caveats.
    warnings: tuple[Warning, ...]


@dataclass(frozen=True, slots=True, eq=False)
class PrAlphaAbResult(_HasWarnings):
    """Result of ``eos.pr_alpha_ab``.

    Three dimensionless outputs and no ``pint`` quantity, like
    :class:`PrKappaResult` - this namespace works in reduced variables throughout,
    so there is no unit anywhere in it to convert or to get wrong.

    ``a_reduced`` and ``b_reduced`` are the ``A`` and ``B`` of the cubic. They are
    deliberately *not* named ``A`` and ``B``: output names are snake_case, which
    `spec_lint` enforces, and the conventional symbols appear in the spec's `latex`
    field instead. That division of labour is the schema's own - `equation` is
    machine-evaluable, `latex` is typeset for a reader.
    """

    #: The alpha function, where Peng-Robinson's temperature dependence lives.
    alpha: float
    #: ``A = a*alpha*P/(R**2*T**2)``, the dimensionless attraction parameter.
    a_reduced: float
    #: ``B = b*P/(R*T)``, the dimensionless repulsion parameter.
    b_reduced: float
    #: Caveats.
    warnings: tuple[Warning, ...]


@dataclass(frozen=True, slots=True, eq=False)
class PrZFactorResult(_HasWarnings):
    """Result of ``eos.pr_z_factor``.

    The first result in the registry to carry a solver report *and* an enum output:
    the cubic's roots are the answer, and ``root_structure`` says how many of them
    were admissible without making the caller compare two floats.

    ``z_min`` and ``z_max`` are named for their position in the ordered root set
    rather than as "liquid" and "vapour" on purpose. Above the critical temperature
    there is one root and neither name is true, and a field called ``z_liquid``
    holding a supercritical compressibility factor is a wrong answer that looks
    like a right one.
    """

    #: The smallest admissible root.
    z_min: float
    #: The largest admissible root. Equal to ``z_min`` when only one is admissible.
    z_max: float
    #: How many admissible roots there were.
    root_structure: RootStructure
    #: Newton steps the polish took, summed over the roots. Carried because the
    #: answer alone does not say whether the solver did any work, and because the
    #: cross-language agreement test compares iteration counts as the sharpest
    #: cheap check that both implementations ran the same scheme.
    iterations: int
    #: Whether the polish met its stopping rule.
    converged: bool
    #: The largest ``|x_k - x_{k-1}|`` at the final polish step.
    residual: float
    #: Caveats.
    warnings: tuple[Warning, ...]


@dataclass(frozen=True, slots=True, eq=False)
class SrkKappaResult(_HasWarnings):
    """Result of ``eos.srk_kappa``."""

    #: The Soave-Redlich-Kwong alpha-function coefficient. Dimensionless.
    kappa: float
    #: Caveats.
    warnings: tuple[Warning, ...]


@dataclass(frozen=True, slots=True, eq=False)
class Pr78KappaResult(_HasWarnings):
    """Result of ``eos.pr78_kappa``."""

    #: The 1978 Peng-Robinson alpha-function coefficient. Dimensionless.
    kappa: float
    #: Caveats.
    warnings: tuple[Warning, ...]


@dataclass(frozen=True, slots=True, eq=False)
class TwuKappaResult(_HasWarnings):
    """Result of ``eos.twu_kappa``."""

    #: Twu's alpha-function coefficient. Dimensionless.
    kappa: float
    #: Caveats.
    warnings: tuple[Warning, ...]


@dataclass(frozen=True, slots=True, eq=False)
class TwucoonAlphaResult(_HasWarnings):
    """Result of ``eos.twucoon_alpha``."""

    #: The temperature-dependent alpha function. Dimensionless.
    alpha: float
    #: Caveats.
    warnings: tuple[Warning, ...]


@dataclass(frozen=True, slots=True, eq=False)
class TwucoonParamAlphaResult(_HasWarnings):
    """Result of ``eos.twucoon_param_alpha``."""

    #: The temperature-dependent alpha function. Dimensionless.
    alpha: float
    #: Caveats.
    warnings: tuple[Warning, ...]


@dataclass(frozen=True, slots=True, eq=False)
class TwucoonStatoilAlphaResult(_HasWarnings):
    """Result of ``eos.twucoon_statoil_alpha``."""

    #: The temperature-dependent alpha function. Dimensionless.
    alpha: float
    #: Caveats.
    warnings: tuple[Warning, ...]


@dataclass(frozen=True, slots=True, eq=False)
class PrDaneshAlphaResult(_HasWarnings):
    """Result of ``eos.pr_danesh_alpha``."""

    #: The temperature-dependent alpha function. Dimensionless.
    alpha: float
    #: Caveats.
    warnings: tuple[Warning, ...]


@dataclass(frozen=True, slots=True, eq=False)
class PrDelft1998AlphaResult(_HasWarnings):
    """Result of ``eos.pr_delft1998_alpha``."""

    #: The temperature-dependent alpha function. Dimensionless.
    alpha: float
    #: Caveats.
    warnings: tuple[Warning, ...]


@dataclass(frozen=True, slots=True, eq=False)
class PrGassem2001AlphaResult(_HasWarnings):
    """Result of ``eos.pr_gassem2001_alpha``."""

    #: The temperature-dependent alpha function. Dimensionless.
    alpha: float
    #: Caveats.
    warnings: tuple[Warning, ...]


@dataclass(frozen=True, slots=True, eq=False)
class PrLeeKeslerAlphaResult(_HasWarnings):
    """Result of ``eos.pr_lee_kesler_alpha``."""

    #: The temperature-dependent alpha function. Dimensionless.
    alpha: float
    #: Caveats.
    warnings: tuple[Warning, ...]


@dataclass(frozen=True, slots=True, eq=False)
class MollerupAlphaResult(_HasWarnings):
    """Result of ``eos.mollerup_alpha``."""

    #: The temperature-dependent alpha function. Dimensionless.
    alpha: float
    #: Caveats.
    warnings: tuple[Warning, ...]


@dataclass(frozen=True, slots=True, eq=False)
class SrkAlphaAbResult(_HasWarnings):
    """Result of ``eos.srk_alpha_ab``."""

    #: The Soave alpha function.
    alpha: float
    #: ``A = a*alpha*P/(R**2*T**2)``, the dimensionless attraction parameter.
    a_reduced: float
    #: ``B = b*P/(R*T)``, the dimensionless repulsion parameter.
    b_reduced: float
    #: Caveats.
    warnings: tuple[Warning, ...]


@dataclass(frozen=True, slots=True, eq=False)
class SrkZFactorResult(_HasWarnings):
    """Result of ``eos.srk_z_factor``."""

    #: The smallest admissible root.
    z_min: float
    #: The largest admissible root. Equal to ``z_min`` when only one is admissible.
    z_max: float
    #: How many admissible roots there were.
    root_structure: RootStructure
    #: Newton steps the polish took, summed over the roots.
    iterations: int
    #: Whether the polish met its stopping rule.
    converged: bool
    #: The largest ``|x_k - x_{k-1}|`` at the final polish step.
    residual: float
    #: Caveats.
    warnings: tuple[Warning, ...]


@dataclass(frozen=True, slots=True, eq=False)
class SrkDepartureResult(_HasWarnings):
    """Result of ``eos.srk_departure``."""

    #: The logarithm of the fugacity coefficient.
    ln_phi: float
    #: The departure enthalpy over ``R*T``.
    h_dep_rt: float
    #: The departure entropy over ``R``.
    s_dep_r: float
    #: The departure heat capacity over ``R``.
    cp_dep_r: float
    #: Caveats.
    warnings: tuple[Warning, ...]


@dataclass(frozen=True, slots=True, eq=False)
class RkAlphaAbResult(_HasWarnings):
    """Result of ``eos.rk_alpha_ab``."""

    #: The Redlich-Kwong alpha function, ``1/sqrt(Tr)``.
    alpha: float
    #: ``A = a*alpha*P/(R**2*T**2)``, the dimensionless attraction parameter.
    a_reduced: float
    #: ``B = b*P/(R*T)``, the dimensionless repulsion parameter.
    b_reduced: float
    #: Caveats.
    warnings: tuple[Warning, ...]


@dataclass(frozen=True, slots=True, eq=False)
class RkDepartureResult(_HasWarnings):
    """Result of ``eos.rk_departure``."""

    #: The logarithm of the fugacity coefficient.
    ln_phi: float
    #: The departure enthalpy over ``R*T``.
    h_dep_rt: float
    #: The departure entropy over ``R``.
    s_dep_r: float
    #: The departure heat capacity over ``R``.
    cp_dep_r: float
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
    #: Caveats.
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


@dataclass(frozen=True, slots=True, eq=False)
class PtFlashResult(_HasWarnings):
    """Result of ``eos.pt_flash``.

    The first result in this library whose shape is not a fixed set of scalars: a
    flash's answer is a composition vector, and a mixture has as many of them as it
    has components.

    ``beta`` is ``None`` in the two cases where there is genuinely no vapour
    fraction, and that is deliberate rather than a convenience. At a trivial
    solution every K-value is 1, ``g(beta)`` is identically zero, and ``beta`` is
    *indeterminate* - successive substitution approaches the point through
    geometrically growing ``beta``, and where a bisection stops on an identically
    zero function is a ratio of round-off. Measured at ``-7.7e10`` for one feed and
    ``-2.2e11`` for another, neither reproducible across implementations. A number
    there would be a fabrication indistinguishable from a real vapour fraction; the
    same reasoning applies to a feed with no Rachford-Rice root at all.

    Branch on ``phase`` and on presence, never on whether the number looks
    plausible.
    """

    #: The vapour fraction, or ``None`` when there is no vapour fraction to report.
    beta: float | None
    #: Liquid-phase mole fractions.
    x: tuple[float, ...]
    #: Vapour-phase mole fractions.
    y: tuple[float, ...]
    #: ``K_i = y_i / x_i``, the iterate the loop converges on.
    k: tuple[float, ...]
    #: ``ln phi_i`` in the liquid phase.
    ln_phi_liquid: tuple[float, ...]
    #: ``ln phi_i`` in the vapour phase.
    ln_phi_vapour: tuple[float, ...]
    #: The liquid root of the cubic, the smallest admissible one.
    z_liquid: float
    #: The vapour root, the largest admissible one.
    z_vapour: float
    #: The smallest ``T / Tc_i`` over the components.
    min_t_over_tc: float
    #: What the converged state is.
    phase: Phase
    #: Successive-substitution steps taken.
    iterations: int
    #: ``rms_i |ln K_i - ln K_i_previous|`` at the last step the loop completed, or
    #: ``NaN`` when no step completed.
    residual: float
    #: Caveats.
    warnings: tuple[Warning, ...]


@dataclass(frozen=True, slots=True, eq=False)
class GeNrtlFlashResult(_HasWarnings):
    """Result of ``eos.ge_nrtl_flash``.

    The same shape as :class:`PtFlashResult`, and the same rule about ``beta``: it is
    ``None`` where there is genuinely no vapour fraction - a feed whose K-values are all
    on one side of one, so Rachford-Rice has no root - and a number, possibly outside
    ``[0, 1]``, where the flash converged to one.

    What differs is where the two phases come from. ``ln_phi_liquid`` is an NRTL
    activity-coefficient phase's ``ln(gamma_i P0_i / P)`` and ``ln_phi_vapour`` is the
    cubic's, so the two are not the same kind of number and are reported apart for the
    same reason the phase model reports ``p_sat`` beside ``ln_phi``.
    """

    #: The vapour fraction, or ``None`` when there is no vapour fraction to report.
    beta: float | None
    #: Liquid-phase mole fractions.
    x: tuple[float, ...]
    #: Vapour-phase mole fractions.
    y: tuple[float, ...]
    #: ``K_i = y_i / x_i``, the iterate the loop converges on.
    k: tuple[float, ...]
    #: ``ln phi_i`` in the liquid, from the activity-coefficient phase.
    ln_phi_liquid: tuple[float, ...]
    #: ``ln phi_i`` in the vapour, from the cubic.
    ln_phi_vapour: tuple[float, ...]
    #: The vapour root of the cubic, the largest admissible one. There is no liquid
    #: root: the liquid is not a cubic.
    z_vapour: float
    #: The smallest ``T / Tc_i`` over the components.
    min_t_over_tc: float
    #: What the converged state is.
    phase: Phase
    #: Successive-substitution steps taken.
    iterations: int
    #: ``rms_i |ln K_i - ln K_i_previous|`` at the last step the loop completed, or
    #: ``NaN`` when no step completed.
    residual: float
    #: Caveats.
    warnings: tuple[Warning, ...]


@dataclass(frozen=True, slots=True, eq=False)
class PtPhaseEnvelopeResult(_HasWarnings):
    """Result of ``eos.pt_phase_envelope``.

    The dew and bubble branches of the phase envelope, each a parallel pair of
    temperature and pressure arrays traced upward from the starting pressure to the
    critical point, together with the cricondenbar, cricondentherm and the critical
    point itself.
    """

    #: The dew-point temperatures, in trace order.
    dew_temperature: tuple[float, ...]
    #: The dew-point pressures, parallel to ``dew_temperature``.
    dew_pressure: tuple[float, ...]
    #: The bubble-point temperatures, in trace order.
    bubble_temperature: tuple[float, ...]
    #: The bubble-point pressures, parallel to ``bubble_temperature``.
    bubble_pressure: tuple[float, ...]
    #: The temperature at the cricondenbar.
    cricondenbar_temperature: Q
    #: The cricondenbar pressure.
    cricondenbar_pressure: Q
    #: The cricondentherm temperature.
    cricondentherm_temperature: Q
    #: The pressure at the cricondentherm.
    cricondentherm_pressure: Q
    #: The critical temperature.
    critical_temperature: Q
    #: The critical pressure.
    critical_pressure: Q
    #: The number of continuation points traced.
    iterations: int
    #: The temperature gap between the two branches' endpoints, in kelvin.
    residual: float
    #: Caveats.
    warnings: tuple[Warning, ...]


@dataclass(frozen=True, slots=True, eq=False)
class PhFlashResult(_HasWarnings):
    """Result of ``eos.ph_flash``.

    The state a mixture reaches when a duty is applied at a fixed pressure: the
    temperature that satisfies the energy balance, and the phase split at it. The
    split is reported in as much detail as :class:`PtFlashResult` because it *is* one
    - the flash evaluated at the answer - and a caller who needs the compositions
    should not have to run it again to get them.
    """

    #: The temperature that satisfies the enthalpy. This is the model's answer.
    T: Q
    #: The vapour fraction at that temperature, or ``None`` for a single-phase feed.
    #: ``None`` rather than a number outside ``[0, 1]``: the flash extrapolates a
    #: split that does not exist, and reporting it would invite a caller to use it.
    beta: float | None
    #: Liquid-phase composition at the answer.
    x: tuple[float, ...]
    #: Vapour-phase composition.
    y: tuple[float, ...]
    #: K-values at the answer.
    k: tuple[float, ...]
    #: Which phase the feed is in at the answer.
    phase: Phase
    #: Liquid root of the cubic at the answer.
    z_liquid: float
    #: Vapour root.
    z_vapour: float
    #: Bisection steps taken.
    iterations: int
    #: ``|H(T) - H_target| / max(|H_target|, 1)`` at the answer.
    residual: float
    #: Caveats, deduplicated - the search evaluates the flash thousands of times.
    warnings: tuple[Warning, ...]


@dataclass(frozen=True, slots=True, eq=False)
class PsFlashResult(_HasWarnings):
    """Result of ``eos.ps_flash``.

    The state a stream reaches when it is expanded or compressed **isentropically** at a
    fixed pressure. Same shape as :class:`PhFlashResult`, deliberately: the two models
    differ in which property they invert and agree on everything else, so a caller
    reading one already knows how to read the other.
    """

    #: The temperature that satisfies the entropy. This is the model's answer.
    T: Q
    #: The vapour fraction at that temperature, or ``None`` for a single-phase feed.
    beta: float | None
    #: Liquid-phase composition at the answer.
    x: tuple[float, ...]
    #: Vapour-phase composition.
    y: tuple[float, ...]
    #: K-values at the answer.
    k: tuple[float, ...]
    #: Which phase the feed is in at the answer.
    phase: Phase
    #: Liquid root of the cubic at the answer.
    z_liquid: float
    #: Vapour root.
    z_vapour: float
    #: Bisection steps taken.
    iterations: int
    #: ``|S(T) - S_target| / max(|S_target|, 1)`` at the answer.
    residual: float
    #: Caveats, deduplicated.
    warnings: tuple[Warning, ...]


@dataclass(frozen=True, slots=True, eq=False)
class TvFlashResult(_HasWarnings):
    """Result of ``eos.tv_flash``.

    The pressure that satisfies a volume at a fixed temperature, and the phase split at
    it. Same shape as :class:`PhFlashResult` but with the pressure as the answer.
    """

    #: The pressure that satisfies the volume. This is the model's answer.
    P: Q
    #: The vapour fraction at that pressure, or ``None`` for a single-phase feed.
    beta: float | None
    #: Liquid-phase composition at the answer.
    x: tuple[float, ...]
    #: Vapour-phase composition.
    y: tuple[float, ...]
    #: K-values at the answer.
    k: tuple[float, ...]
    #: Which phase the feed is in at the answer.
    phase: Phase
    #: Liquid root of the cubic at the answer.
    z_liquid: float
    #: Vapour root.
    z_vapour: float
    #: Newton steps taken.
    iterations: int
    #: ``V(P) - V_target`` at the answer.
    residual: float
    #: Caveats, deduplicated.
    warnings: tuple[Warning, ...]


@dataclass(frozen=True, slots=True, eq=False)
class PvFlashResult(_HasWarnings):
    """Result of ``eos.pv_flash``.

    The temperature that satisfies a volume at a fixed pressure, and the phase split at
    it. Same shape as :class:`PhFlashResult`.
    """

    #: The temperature that satisfies the volume. This is the model's answer.
    T: Q
    #: The vapour fraction at that temperature, or ``None`` for a single-phase feed.
    beta: float | None
    #: Liquid-phase composition at the answer.
    x: tuple[float, ...]
    #: Vapour-phase composition.
    y: tuple[float, ...]
    #: K-values at the answer.
    k: tuple[float, ...]
    #: Which phase the feed is in at the answer.
    phase: Phase
    #: Liquid root of the cubic at the answer.
    z_liquid: float
    #: Vapour root.
    z_vapour: float
    #: Newton steps taken.
    iterations: int
    #: ``V(T) - V_target`` at the answer.
    residual: float
    #: Caveats, deduplicated.
    warnings: tuple[Warning, ...]


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


@dataclass(frozen=True, slots=True, eq=False)
class ThFlashResult(_HasWarnings):
    """Result of ``eos.th_flash``.

    Same shape as :class:TvFlashResult, with the pressure as the answer.
    """

    #: The pressure that satisfies the property. This is the model's answer.
    P: Q
    #: The vapour fraction at that pressure, or ``None`` for a single-phase feed.
    beta: float | None
    #: Liquid-phase composition at the answer.
    x: tuple[float, ...]
    #: Vapour-phase composition.
    y: tuple[float, ...]
    #: K-values at the answer.
    k: tuple[float, ...]
    #: Which phase the feed is in at the answer.
    phase: Phase
    #: Liquid root of the cubic at the answer.
    z_liquid: float
    #: Vapour root.
    z_vapour: float
    #: Newton steps taken.
    iterations: int
    #: The relative property residual at the answer.
    residual: float
    #: Caveats, deduplicated.
    warnings: tuple[Warning, ...]


@dataclass(frozen=True, slots=True, eq=False)
class TsFlashResult(_HasWarnings):
    """Result of ``eos.ts_flash``.

    Same shape as :class:TvFlashResult, with the pressure as the answer.
    """

    #: The pressure that satisfies the property. This is the model's answer.
    P: Q
    #: The vapour fraction at that pressure, or ``None`` for a single-phase feed.
    beta: float | None
    #: Liquid-phase composition at the answer.
    x: tuple[float, ...]
    #: Vapour-phase composition.
    y: tuple[float, ...]
    #: K-values at the answer.
    k: tuple[float, ...]
    #: Which phase the feed is in at the answer.
    phase: Phase
    #: Liquid root of the cubic at the answer.
    z_liquid: float
    #: Vapour root.
    z_vapour: float
    #: Newton steps taken.
    iterations: int
    #: The relative property residual at the answer.
    residual: float
    #: Caveats, deduplicated.
    warnings: tuple[Warning, ...]


@dataclass(frozen=True, slots=True, eq=False)
class TuFlashResult(_HasWarnings):
    """Result of ``eos.tu_flash``.

    Same shape as :class:TvFlashResult, with the pressure as the answer.
    """

    #: The pressure that satisfies the property. This is the model's answer.
    P: Q
    #: The vapour fraction at that pressure, or ``None`` for a single-phase feed.
    beta: float | None
    #: Liquid-phase composition at the answer.
    x: tuple[float, ...]
    #: Vapour-phase composition.
    y: tuple[float, ...]
    #: K-values at the answer.
    k: tuple[float, ...]
    #: Which phase the feed is in at the answer.
    phase: Phase
    #: Liquid root of the cubic at the answer.
    z_liquid: float
    #: Vapour root.
    z_vapour: float
    #: Newton steps taken.
    iterations: int
    #: The relative property residual at the answer.
    residual: float
    #: Caveats, deduplicated.
    warnings: tuple[Warning, ...]


@dataclass(frozen=True, slots=True, eq=False)
class PuFlashResult(_HasWarnings):
    """Result of ``eos.pu_flash``.

    Same shape as :class:PvFlashResult, with the temperature as the answer.
    """

    #: The temperature that satisfies the property. This is the model's answer.
    T: Q
    #: The vapour fraction at that temperature, or ``None`` for a single-phase feed.
    beta: float | None
    #: Liquid-phase composition at the answer.
    x: tuple[float, ...]
    #: Vapour-phase composition.
    y: tuple[float, ...]
    #: K-values at the answer.
    k: tuple[float, ...]
    #: Which phase the feed is in at the answer.
    phase: Phase
    #: Liquid root of the cubic at the answer.
    z_liquid: float
    #: Vapour root.
    z_vapour: float
    #: Newton steps taken.
    iterations: int
    #: The relative property residual at the answer.
    residual: float
    #: Caveats, deduplicated.
    warnings: tuple[Warning, ...]


@dataclass(frozen=True, slots=True, eq=False)
class StabilityTestResult(_HasWarnings):
    """Result of ``eos.stability_test``.

    Two trials, always, so `tm` and `w` are fixed-length rather than "one per phase
    found": a trial that converges trivially still has a distance, and it is that
    near-zero number which is the evidence it was trivial.
    """

    #: Whether the feed is stable as a single phase.
    verdict: StabilityVerdict
    #: The tangent-plane distance at each trial's stationary point, vapour-like
    #: trial first. Negative means that trial lies below the tangent plane.
    tm: tuple[float, ...]
    #: The stationary-point composition of each trial, in the same order.
    w: tuple[tuple[float, ...], ...]
    #: Iterations each trial took, in the same order.
    iterations: tuple[int, ...]
    #: The smallest ``T / Tc_i`` over the components.
    min_t_over_tc: float
    #: Caveats.
    warnings: tuple[Warning, ...]


class TpMultiflashSeed(StrEnum):
    """Which of the two things decided the phase count.

    A boolean would say a phase was added; this says whether the answer *is* the two-phase
    flash's or is something the tangent-plane trial found, which is the distinction between
    "this model is `eos.pt_flash`" and "this model found a third phase".
    """

    TWO_PHASE_FLASH = "two_phase_flash"
    STABILITY_SEEDED = "stability_seeded"


@dataclass(frozen=True, slots=True, eq=False)
class TpMultiflashResult(_HasWarnings):
    """Result of ``eos.tp_multiflash``.

    The phases are reported in the order the solve kept them, and **no order is promised**:
    upstream sorts its phases by density and this does not, so a caller matching phase *k*
    across two implementations is matching nothing. `z_factor` is what identifies a phase
    physically - two phases on the same side of the cubic differ in composition and in the
    root they sit on, and the root is the compressibility.
    """

    #: How many phases the feed splits into: 1, 2 or 3.
    phase_count: int
    #: The mole fraction of the feed in each phase, summing to one.
    beta: tuple[float, ...]
    #: The composition of each phase, one tuple per phase, each summing to one.
    x: tuple[tuple[float, ...], ...]
    #: The root of the cubic each phase sits on, as the compressibility ``Z = PV/RT``.
    z_factor: tuple[float, ...]
    #: ``ln phi_i`` in each phase, one tuple per phase.
    ln_phi: tuple[tuple[float, ...], ...]
    #: Whether the tangent-plane trial added a phase that survived the merge.
    seeded: TpMultiflashSeed
    #: The tangent-plane distance at each trial's stationary point.
    tm: tuple[float, ...]
    #: Fraction-solve steps taken.
    iterations: int
    #: The larger of the last step's norm and the gradient norm it was judged beside.
    residual: float
    #: The smallest ``T / Tc_i`` over the components.
    min_t_over_tc: float
    #: Caveats.
    warnings: tuple[Warning, ...]


@dataclass(frozen=True, slots=True, eq=False)
class CapillaryDewPointResult(_HasWarnings):
    """Result of ``eos.capillary_dew_point``.

    ``DewTemperatureResult`` with one field more, which is the whole of the difference: the
    Young-Laplace pressure across the curved interface, reported because it is the one quantity
    this model has that the flat one does not.
    """

    #: The dew-point temperature.
    temperature: Q
    #: The composition of the liquid that first appears.
    incipient: tuple[float, ...]
    #: K-values at the converged temperature, **with the Kelvin shift applied**.
    k: tuple[float, ...]
    #: The liquid root of the cubic at the converged state.
    z_liquid: float
    #: The vapour root.
    z_vapour: float
    #: The Young-Laplace pressure across the interface, ``2 sigma cos(theta) / r``.
    capillary_pressure: Q
    #: The smallest ``T / Tc_i`` over the components at the answer.
    min_t_over_tc: float
    #: Temperature updates taken.
    iterations: int
    #: ``|sum_i y_i / K_i - 1|`` at the last completed step.
    residual: float
    #: Caveats.
    warnings: tuple[Warning, ...]


@dataclass(frozen=True, slots=True, eq=False)
class BubblePressureResult(_HasWarnings):
    """Result of ``eos.bubble_pressure``.

    Shares its shape with :class:`DewPressureResult` rather than being one type with
    a switch, because the two differ in *which* composition is the input - ``x`` here
    and ``y`` there - and a shared field would have to be named after neither.
    """

    #: The bubble-point pressure.
    pressure: Q
    #: The composition of the vapour that first appears.
    incipient: tuple[float, ...]
    #: K-values at the converged pressure.
    k: tuple[float, ...]
    #: The liquid root of the cubic at the converged state.
    z_liquid: float
    #: The vapour root.
    z_vapour: float
    #: The smallest ``T / Tc_i`` over the components.
    min_t_over_tc: float
    #: Pressure updates taken.
    iterations: int
    #: ``|sum_i x_i K_i - 1|`` at the last completed step.
    residual: float
    #: Caveats.
    warnings: tuple[Warning, ...]


@dataclass(frozen=True, slots=True, eq=False)
class DewPressureResult(_HasWarnings):
    """Result of ``eos.dew_pressure``."""

    #: The dew-point pressure.
    pressure: Q
    #: The composition of the liquid that first appears.
    incipient: tuple[float, ...]
    #: K-values at the converged pressure.
    k: tuple[float, ...]
    #: The liquid root of the cubic at the converged state.
    z_liquid: float
    #: The vapour root.
    z_vapour: float
    #: The smallest ``T / Tc_i`` over the components.
    min_t_over_tc: float
    #: Pressure updates taken.
    iterations: int
    #: ``|sum_i y_i / K_i - 1|`` at the last completed step.
    residual: float
    #: Caveats.
    warnings: tuple[Warning, ...]


@dataclass(frozen=True, slots=True, eq=False)
class BubbleTemperatureResult(_HasWarnings):
    """Result of ``eos.bubble_temperature``.

    Shares its shape with :class:`DewTemperatureResult` rather than being one type with
    a switch, because the two differ in *which* composition is the input - ``x`` here
    and ``y`` there - and a shared field would have to be named after neither.
    """

    #: The bubble-point temperature.
    temperature: Q
    #: The composition of the vapour that first appears.
    incipient: tuple[float, ...]
    #: K-values at the converged temperature.
    k: tuple[float, ...]
    #: The liquid root of the cubic at the converged state.
    z_liquid: float
    #: The vapour root.
    z_vapour: float
    #: The smallest ``T / Tc_i`` over the components.
    min_t_over_tc: float
    #: Temperature updates taken.
    iterations: int
    #: ``|sum_i x_i K_i - 1|`` at the last completed step.
    residual: float
    #: Caveats.
    warnings: tuple[Warning, ...]


@dataclass(frozen=True, slots=True, eq=False)
class DewTemperatureResult(_HasWarnings):
    """Result of ``eos.dew_temperature``."""

    #: The dew-point temperature.
    temperature: Q
    #: The composition of the liquid that first appears.
    incipient: tuple[float, ...]
    #: K-values at the converged temperature.
    k: tuple[float, ...]
    #: The liquid root of the cubic at the converged state.
    z_liquid: float
    #: The vapour root.
    z_vapour: float
    #: The smallest ``T / Tc_i`` over the components.
    min_t_over_tc: float
    #: Temperature updates taken.
    iterations: int
    #: ``|sum_i x_i / K_i - 1|`` at the last completed step.
    residual: float
    #: Caveats.
    warnings: tuple[Warning, ...]


@dataclass(frozen=True, slots=True, eq=False)
class HybridEosGeFlashResult(_HasWarnings):
    """Result of ``eos.hybrid_eos_ge_flash``.

    The three roles are **fixed before the fractions are solved** - gas, oil and a GE
    aqueous phase - so the order of every vector here is the model's and not the answer's,
    and no ``role`` field is needed to say which phase is which.
    """

    #: The mole fraction of the feed in each role, in ``[gas, oil, aqueous]`` order.
    beta: tuple[float, ...]
    #: The composition of each role, one row per role, each summing to one.
    x: tuple[tuple[float, ...], ...]
    #: ``ln phi_i`` in each role, from the cubic for the first two and ``eos.pitzer_phase``
    #: for the brine.
    ln_phi: tuple[tuple[float, ...], ...]
    #: Newton steps taken, summed over the fixed-topology passes.
    iterations: int
    #: The larger of the last step's norm and the last gradient's.
    residual: float
    #: The worst material-balance residual. The contract holds it to ``1e-7``.
    max_material_balance_residual: float
    #: The worst cross-role ``ln(x_i phi_i P)`` spread. The contract holds it to ``1e-5``.
    max_log_fugacity_residual: float
    #: The smallest ``T / Tc_i`` over the components.
    min_t_over_tc: float
    #: Caveats.
    warnings: tuple[Warning, ...]


@dataclass(frozen=True, slots=True, eq=False)
class ReactiveHybridEosGeFlashResult(_HasWarnings):
    """Result of ``reactions.reactive_hybrid_eos_ge_flash``.

    A coupled state: the role split and the brine's chemistry solved at once, each as the
    other's input. **A loop that does not certify is an error and not an answer**, so every
    returned result is one whose three conditions held - the pass floor, the brine's
    composition and the fraction solve's residual.
    """

    #: The mole fraction of the feed in each of ``[gas, oil, aqueous]``, summing to one. A
    #: role the solve drove to nothing carries the solver's floor rather than zero.
    beta: tuple[float, ...]
    #: Each role's composition at the coupled state, one row per role and one column per
    #: component. An ion's entry is ``1e-50`` in the two EoS roles.
    x: tuple[tuple[float, ...], ...]
    #: The reaction-adjusted overall inventory the last pass solved at. **Not the feed**
    #: wherever the chemistry moved a species, and a spectator's entry is the feed's.
    coupled_moles: tuple[Q, ...]
    #: The brine's species amounts at the answer, in the reactive set's own order - the
    #: order :func:`azoth.reactions.reference.reactive_hybrid_eos_ge_flash.reactive_components`
    #: reports. This is the state a scale potential is computed from.
    aqueous_moles: tuple[Q, ...]
    #: How many coupled passes the loop made: at least three, and the class's own count
    #: rather than the fraction solve's.
    passes: int
    #: The last pass's composition deviation, the sum of ``|x_old - x_new|`` over the brine.
    #: The class's tolerance on it is ``1e-10``.
    chemical_deviation: float
    #: The last fraction solve's own residual, the larger of its step norm and its gradient
    #: norm, held to ``1e-10``.
    residual: float
    #: The worst ``|z_i - sum_k beta_k x_ik|`` over the coupled inventory. The underlying
    #: model's contract holds it to ``1e-7``.
    max_material_balance_residual: float
    #: The worst cross-role ``ln(x_i phi_i P)`` spread, held to ``1e-5``.
    max_log_fugacity_residual: float
    #: The worst departure of an element balance from what the phases hold, in mol. The
    #: chemistry conserves elements by construction, so this measures the split.
    element_residual: float
    #: The net charge the phases fail to account for, in moles of elementary charge.
    charge_residual: float
    #: Caveats.
    warnings: tuple[Warning, ...]


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


@dataclass(frozen=True, slots=True, eq=False)
class IapwsHenryLawResult(_HasWarnings):
    """Result of ``eos.iapws_henry_law``.

    The standard state is limiting ``f/x`` at pure-water saturation, so ``henry`` is a
    pressure per mole fraction and is large for a sparingly soluble gas.
    ``d_ln_henry_d_t`` comes from the guideline's own logarithmic expression rather than
    from differencing the constant.
    """

    #: ``kH``, the Henry constant.
    henry: Q
    #: ``ln kH``, the logarithm of the value above and not of NeqSim's bar figure.
    ln_henry: float
    #: ``d(ln kH)/dT``.
    d_ln_henry_d_t: float
    #: Whether ``T`` is inside the row's fitted range.
    status: HenryStatus
    #: The row's reported root-mean-square residual in ``ln kH``.
    rms_log_residual: float
    #: Caveats.
    warnings: tuple[Warning, ...]


@dataclass(frozen=True, slots=True, eq=False)
class IdealGasCpResult(_HasWarnings):
    """Result of ``eos.ideal_gas_cp``.

    Reports the polynomial's own dimensionless value as well as the dimensioned heat
    capacity, because the dimensionless form is what a reader checking the arithmetic
    by hand computes first - and because the multiplication by ``R`` is then visible
    as the single step it is, rather than folded into the answer.
    """

    #: The ideal-gas heat capacity. Carries ``OUT_OF_VALID_RANGE`` when it is not
    #: positive, which means the polynomial has been evaluated outside its range.
    cp: Q
    #: Caveats.
    warnings: tuple[Warning, ...]


@dataclass(frozen=True, slots=True, eq=False)
class MolarEnthalpyEntropyResult(_HasWarnings):
    """Result of ``eos.molar_enthalpy_entropy``.

    The ideal-gas and departure parts are reported separately as well as summed,
    because the split is what a caller checking the answer needs: ``h_ideal`` carries
    the datum, ``h_departure`` carries the equation of state, and a single total hides
    which of the two a disagreement came from.
    """

    #: The molar enthalpy, ``h_ideal + h_departure``.
    h: Q
    #: The molar entropy, ``s_ideal + s_departure``.
    s: Q
    #: The ideal-gas part of the enthalpy - the reference values and the integrals.
    h_ideal: Q
    #: The ideal-gas part of the entropy.
    s_ideal: Q
    #: The residual enthalpy, ``R*T*h_dep_rt``.
    h_departure: Q
    #: The residual entropy, ``R*s_dep_r``. Does **not** include the entropy of mixing.
    s_departure: Q
    #: The composition-weighted average of the components' ``psi``.
    psi_bar: float
    #: The molar heat capacity at constant pressure, ``cp_ideal + cp_departure``.
    #: The derivative the isentropic and isenthalpic flashes step on: theirs is
    #: ``dS/dT = cp/T`` in the entropy's case, ``dH/d(1/T) = -T**2*cp`` in the
    #: enthalpy's.
    cp: Q
    #: The ideal-gas part of the heat capacity - the polynomial, at ``T``.
    cp_ideal: Q
    #: The residual heat capacity, ``R*cp_dep_r``.
    cp_departure: Q
    #: Caveats.
    warnings: tuple[Warning, ...]


@dataclass(frozen=True, slots=True, eq=False)
class CriticalPointResult(_HasWarnings):
    """Result of ``eos.critical_point``.

    The four state variables of a mixture critical point. ``Z_c`` is here rather than
    derived by the caller because it is the quantity that distinguishes this model from
    the mechanical conditions: a pure component's is ``(1 - omega_b)/3``, a mixture's
    varies with composition, and the mechanical route cannot produce the second.
    """

    #: The critical temperature.
    tc: Q
    #: The critical pressure.
    pc: Q
    #: The critical molar volume.
    vc: Q
    #: ``Pc Vc/(R Tc)``.
    z_c: float
    #: Outer iterations taken.
    iterations: int
    #: ``max(|smallest eigenvalue|, |cubic form|)`` at the returned state.
    residual: float
    #: Caveats.
    warnings: tuple[Warning, ...]


@dataclass(frozen=True, slots=True, eq=False)
class BwrsPhaseResult(_HasWarnings):
    """Result of ``eos.bwrs_phase``.

    The BWRS (MBWR-32) phase state: the compressibility factor, the fugacity
    coefficients and the residual departures. The departures are residual - real minus
    ideal gas at the same state - so a caller composing them with
    :func:`azoth.eos.molar_enthalpy_entropy` gets absolute values.
    """

    #: The compressibility factor ``Z = P/(rho R T)``.
    z_factor: float
    #: The fugacity coefficients, as logarithms, one per component.
    ln_phi: tuple[float, ...]
    #: The residual enthalpy.
    h_res: Q
    #: The residual entropy.
    s_res: Q
    #: The residual isobaric heat capacity.
    cp_res: Q
    #: Caveats.
    warnings: tuple[Warning, ...]


@dataclass(frozen=True, slots=True, eq=False)
class TpFlashSaftResult(_HasWarnings):
    """Result of ``eos.tp_flash_saft``.

    The SAFT-VR-Mie flash: the vapour fraction, both compositions and both phases' fugacity
    coefficients. **``beta`` is absent rather than zero when the flash reports one phase**,
    because a single-phase answer has no vapour fraction - and which phase it is comes from
    the Gibbs comparison, not from the loop.
    """

    #: The vapour fraction, present only when the flash found a split.
    beta: float | None
    #: Liquid-phase mole fractions, the feed itself when there is one phase.
    x: tuple[float, ...]
    #: Vapour-phase mole fractions, the feed itself when there is one phase.
    y: tuple[float, ...]
    #: ``K_i = phi_liquid_i / phi_vapour_i``, the iterate the loop converges on.
    k: tuple[float, ...]
    #: ``ln phi_i`` in the liquid phase.
    ln_phi_liquid: tuple[float, ...]
    #: ``ln phi_i`` in the vapour phase.
    ln_phi_vapour: tuple[float, ...]
    #: The compressibility factor of the liquid solve.
    z_liquid: float
    #: The compressibility factor of the vapour solve.
    z_vapour: float
    #: What the converged state is.
    phase: Phase
    #: Successive-substitution steps taken.
    iterations: int
    #: The largest relative change in a K-value at the last step.
    residual: float
    #: Caveats.
    warnings: tuple[Warning, ...]


@dataclass(frozen=True, slots=True, eq=False)
class SaftVrMiePhaseResult(_HasWarnings):
    """Result of ``eos.saft_vr_mie_phase``.

    The SAFT-VR-Mie phase state: the compressibility factor, the molar volume the solve
    converged to, the fugacity coefficients and the residual departures.

    **There is no ``cp_res``.** It needs ``d^2(A^R/RT)/dT^2``, and the chain contact value's
    second temperature derivative is not derived in either kernel, so a heat capacity from
    the dispersion alone would be a wrong answer rather than a missing one.
    """

    #: The compressibility factor at the chosen root.
    z_factor: float
    #: The fugacity coefficients, as logarithms, one per component.
    ln_phi: tuple[float, ...]
    #: The molar volume at the chosen root.
    v: Q
    #: The residual enthalpy, real minus ideal gas at the same state.
    h_res: Q
    #: The residual entropy, real minus ideal gas at the same state.
    s_res: Q
    #: Caveats.
    warnings: tuple[Warning, ...]


@dataclass(frozen=True, slots=True, eq=False)
class PcsaftRahmatPhaseResult(_HasWarnings):
    """Result of ``eos.pcsaft_rahmat_phase``.

    The PC-SAFT phase state: the compressibility factor, the molar volume the solve
    converged to, and the fugacity coefficients. **No departure function is reported**,
    because the temperature derivative of the Helmholtz energy is not derived - an
    enthalpy taken from somewhere else would be a different fluid's.
    """

    #: The compressibility factor at the chosen root.
    z_factor: float
    #: The fugacity coefficients, as logarithms, one per component.
    ln_phi: tuple[float, ...]
    #: The molar volume at the chosen root.
    v: Q
    #: The residual enthalpy, real minus ideal gas at the same state.
    h_res: Q
    #: The residual entropy, real minus ideal gas at the same state.
    s_res: Q
    #: Caveats.
    warnings: tuple[Warning, ...]


@dataclass(frozen=True, slots=True, eq=False)
class SrkCpaPhaseResult(_HasWarnings):
    """Result of ``eos.srk_cpa_phase``.

    The Soave-Redlich-Kwong CPA phase state: the compressibility factor, the fugacity
    coefficients and the residual departures, at a root the *association* chose rather
    than the cubic - at 300 K and 100 bar this fluid's associating root is 0.10505 where
    its substituted cubic's is 0.15229.

    **There is no ``cp_res``.** The association's second temperature derivative is not
    derived, so an associating mixture's enthalpy and heat capacity are not consistent;
    reporting a departure heat capacity from the cubic alone would be a wrong answer
    rather than a missing one.
    """

    #: The compressibility factor at the chosen root.
    z_factor: float
    #: The fugacity coefficients, as logarithms, one per component.
    ln_phi: tuple[float, ...]
    #: The residual enthalpy, including the association's ``A/(RT) - T d(A/RT)/dT``.
    h_res: Q
    #: The residual entropy.
    s_res: Q
    #: Caveats.
    warnings: tuple[Warning, ...]


@dataclass(frozen=True, slots=True, eq=False)
class PrCpaPhaseResult(_HasWarnings):
    """Result of ``eos.pr_cpa_phase``.

    The Peng-Robinson CPA phase state: the same model as ``eos.srk_cpa_phase`` under a
    different cubic and a different fitted set. Water's ``kappa_AB`` is 0.0692 for SRK
    against 0.046473789 for PR and its fitted covolume is 1.4515 against 1.456360879, so
    the two are separate fits rather than one converted into the other.

    **There is no NeqSim reading for this model.** Against the pinned jar
    ``SystemPrCPA`` builds CPA components whose sites its phase never sums, so its
    association is computed over nothing - see the spec's assumptions. The divergence is
    recorded with the tranche's others; what is not given up is the two-kernel comparison,
    which this model's case runs.
    """

    #: The compressibility factor at the chosen root.
    z_factor: float
    #: The fugacity coefficients, as logarithms, one per component.
    ln_phi: tuple[float, ...]
    #: The residual enthalpy, including the association's ``A/(RT) - T d(A/RT)/dT``.
    h_res: Q
    #: The residual entropy.
    s_res: Q
    #: Caveats.
    warnings: tuple[Warning, ...]


@dataclass(frozen=True, slots=True, eq=False)
class UmrCpaPhaseResult(_HasWarnings):
    """Result of ``eos.umr_cpa_phase``.

    The only model in this library whose attraction is mixed by a universal rule rather
    than an interaction matrix: ``alpha_mix = sum_i x_i (a_i^T/(b_i R T) + hwfc ln
    gamma_i)`` with ``hwfc = -1/0.53`` over UNIFAC-UMR-PRU's activity coefficients, and
    ``A = n B R T alpha_mix``. No ``kij`` column is read.

    **There is no flash for this model.** ``eos.pt_flash``'s Jacobian is ``d ln phi / d
    n``, which every excess-Gibbs rule refuses because that derivative is the excess
    Gibbs energy's own Hessian.
    """

    #: The compressibility factor at the chosen root.
    z_factor: float
    #: The fugacity coefficients, as logarithms, one per component.
    ln_phi: tuple[float, ...]
    #: The residual enthalpy, including the association's ``-T d(A/RT)/dT``.
    h_res: Q
    #: The residual entropy.
    s_res: Q
    #: Caveats.
    warnings: tuple[Warning, ...]


@dataclass(frozen=True, slots=True, eq=False)
class FurstElectrolyteMod2004PhaseResult(_HasWarnings):
    """Result of ``eos.furst_electrolyte_mod2004_phase``.

    The 2004 revision of :class:`FurstElectrolytePhaseResult`: the same kernels with the
    solvent dielectric constant's temperature and composition derivatives zeroed, and a
    component-independent ``FBornD`` added to every ``ln phi`` - which makes that entry
    depend on the phase's *size*, so this model's cases carry a looser tolerance.
    """

    #: The compressibility factor at the chosen root.
    z_factor: float
    #: The fugacity coefficients, as logarithms, one per component.
    ln_phi: tuple[float, ...]
    #: Caveats.
    warnings: tuple[Warning, ...]


@dataclass(frozen=True, slots=True, eq=False)
class FurstElectrolytePhaseResult(_HasWarnings):
    """Result of ``eos.furst_electrolyte_phase``.

    An SRK cubic with Schwartzentruber's attractive term, three additive Helmholtz terms and
    the Huron-Vidal mixing rule - NeqSim's ``SystemFurstElectrolyteEos``. **The salt is a
    component and not a scalar**, so ``ln_phi`` carries an entry per ion: an ion's is large
    and negative, because its attraction is ``1e-35`` and its covolume is fitted.

    **There is no flash.** ``eos.pt_flash``'s Jacobian is ``d ln phi / d n``, which the
    Huron-Vidal rule refuses because that derivative is the excess Gibbs energy's second
    derivative.
    """

    #: The compressibility factor at the chosen root.
    z_factor: float
    #: The fugacity coefficients, as logarithms, one per component.
    ln_phi: tuple[float, ...]
    #: Caveats.
    warnings: tuple[Warning, ...]


@dataclass(frozen=True, slots=True, eq=False)
class SoreideWhitsonPhaseResult(_HasWarnings):
    """Result of ``eos.soreide_whitson_phase``.

    Peng-Robinson 1978 with a brine in it: water's alpha and the water-gas interaction
    parameter of a water-rich phase both read the salinity, and every other component and
    pair is PR78.

    **There is no flash for this model.** ``eos.pt_flash``'s Jacobian is ``d ln phi / d
    n``, and this rule's interaction matrix moves with the composition being
    differentiated against.
    """

    #: The compressibility factor at the chosen root.
    z_factor: float
    #: The fugacity coefficients, as logarithms, one per component.
    ln_phi: tuple[float, ...]
    #: Caveats.
    warnings: tuple[Warning, ...]


@dataclass(frozen=True, slots=True, eq=False)
class AmmoniaPhaseResult(_HasWarnings):
    """Result of ``eos.ammonia_phase``.

    The ammonia reference phase state: the compressibility factor and the Helmholtz
    property set (internal energy, enthalpy, entropy, heat capacities and Gibbs energy).
    """

    #: The compressibility factor ``Z = P/(rho R T)``.
    z_factor: float
    #: The internal energy.
    u: Q
    #: The enthalpy.
    h: Q
    #: The entropy.
    s: Q
    #: The isochoric heat capacity.
    cv: Q
    #: The isobaric heat capacity.
    cp: Q
    #: The Gibbs energy, ``h - T s``.
    g: Q
    #: Caveats.
    warnings: tuple[Warning, ...]


@dataclass(frozen=True, slots=True, eq=False)
class Co2PhaseResult(_HasWarnings):
    """Result of ``eos.co2_phase``.

    The Span-Wagner CO2 phase state: the compressibility factor and the Helmholtz
    property set (internal energy, enthalpy, entropy, heat capacities and Gibbs energy).
    """

    #: The compressibility factor ``Z = P/(rho R T)``.
    z_factor: float
    #: The internal energy.
    u: Q
    #: The enthalpy.
    h: Q
    #: The entropy.
    s: Q
    #: The isochoric heat capacity.
    cv: Q
    #: The isobaric heat capacity.
    cp: Q
    #: The Gibbs energy, ``h - T s``.
    g: Q
    #: Caveats.
    warnings: tuple[Warning, ...]


@dataclass(frozen=True, slots=True, eq=False)
class HydrogenPhaseResult(_HasWarnings):
    """Result of ``eos.hydrogen_phase``.

    The Leachman hydrogen phase state: the compressibility factor and the Helmholtz
    property set (internal energy, enthalpy, entropy, heat capacities and Gibbs energy).
    """

    #: The compressibility factor ``Z = P/(rho R T)``.
    z_factor: float
    #: The internal energy.
    u: Q
    #: The enthalpy.
    h: Q
    #: The entropy.
    s: Q
    #: The isochoric heat capacity.
    cv: Q
    #: The isobaric heat capacity.
    cp: Q
    #: The Gibbs energy, ``h - T s``.
    g: Q
    #: Caveats.
    warnings: tuple[Warning, ...]


@dataclass(frozen=True, slots=True, eq=False)
class WaterPhaseResult(_HasWarnings):
    """Result of ``eos.water_phase``.

    The IAPWS-IF97 water phase state: the compressibility factor and the Gibbs property
    set (internal energy, enthalpy, entropy, heat capacities and Gibbs energy).
    """

    #: The compressibility factor ``Z = P v/(R T)``.
    z_factor: float
    #: The internal energy.
    u: Q
    #: The enthalpy.
    h: Q
    #: The entropy.
    s: Q
    #: The isochoric heat capacity.
    cv: Q
    #: The isobaric heat capacity.
    cp: Q
    #: The Gibbs energy, ``h - T s``.
    g: Q
    #: Caveats.
    warnings: tuple[Warning, ...]


@dataclass(frozen=True, slots=True, eq=False)
class ArgonSolidPhaseResult(_HasWarnings):
    """Result of ``eos.argon_solid_phase``.

    The solid argon phase state: the compressibility factor and the Helmholtz property
    set (internal energy, enthalpy, entropy, heat capacities and Gibbs energy).
    """

    #: The compressibility factor ``Z = P v/(R T)``.
    z_factor: float
    #: The internal energy.
    u: Q
    #: The enthalpy.
    h: Q
    #: The entropy.
    s: Q
    #: The isochoric heat capacity.
    cv: Q
    #: The isobaric heat capacity.
    cp: Q
    #: The Gibbs energy, ``h - T s``.
    g: Q
    #: Caveats.
    warnings: tuple[Warning, ...]


@dataclass(frozen=True, slots=True, eq=False)
class ParahydrogenSolidPhaseResult(_HasWarnings):
    """Result of ``eos.parahydrogen_solid_phase``.

    The solid para-hydrogen phase state: the compressibility factor and the Helmholtz
    property set (internal energy, enthalpy, entropy, heat capacities and Gibbs energy).
    """

    #: The compressibility factor ``Z = P v/(R T)``.
    z_factor: float
    #: The internal energy.
    u: Q
    #: The enthalpy.
    h: Q
    #: The entropy.
    s: Q
    #: The isochoric heat capacity.
    cv: Q
    #: The isobaric heat capacity.
    cp: Q
    #: The Gibbs energy, ``h - T s``.
    g: Q
    #: Caveats.
    warnings: tuple[Warning, ...]


@dataclass(frozen=True, slots=True, eq=False)
class EosCgPhaseResult(_HasWarnings):
    """Result of ``eos.eos_cg_phase``.

    The EOS-CG phase state: the compressibility factor and the Helmholtz property set
    (internal energy, enthalpy, entropy, heat capacities and Gibbs energy).
    """

    #: The compressibility factor ``Z = P v/(R T)``.
    z_factor: float
    #: The internal energy.
    u: Q
    #: The enthalpy.
    h: Q
    #: The entropy.
    s: Q
    #: The isochoric heat capacity.
    cv: Q
    #: The isobaric heat capacity.
    cp: Q
    #: The Gibbs energy, ``h - T s``.
    g: Q
    #: Caveats.
    warnings: tuple[Warning, ...]


@dataclass(frozen=True, slots=True, eq=False)
class GeNrtlPhaseResult(_HasWarnings):
    """Result of ``eos.ge_nrtl_phase``."""

    #: The NRTL activity coefficient of each component.
    gamma: tuple[float, ...]
    #: The natural logarithm of each activity coefficient.
    ln_gamma: tuple[float, ...]
    #: The natural logarithm of each fugacity coefficient, ``ln(gamma_i P0_i / P)``.
    ln_phi: tuple[float, ...]
    #: The pure-component saturation pressure of each component at ``T``.
    p_sat: tuple[Q, ...]
    #: Caveats.
    warnings: tuple[Warning, ...]


@dataclass(frozen=True, slots=True, eq=False)
class GeUnifacPhaseResult(_HasWarnings):
    """Result of ``eos.ge_unifac_phase``."""

    #: The UNIFAC activity coefficient of each component.
    gamma: tuple[float, ...]
    #: The natural logarithm of each activity coefficient.
    ln_gamma: tuple[float, ...]
    #: The natural logarithm of each fugacity coefficient, ``ln(gamma_i P0_i / P)``.
    ln_phi: tuple[float, ...]
    #: The pure-component saturation pressure of each component at ``T``.
    p_sat: tuple[Q, ...]
    #: Caveats.
    warnings: tuple[Warning, ...]


@dataclass(frozen=True, slots=True, eq=False)
class GeWilsonPhaseResult(_HasWarnings):
    """Result of ``eos.ge_wilson_phase``."""

    #: The Wilson activity coefficient of each component.
    gamma: tuple[float, ...]
    #: The natural logarithm of each activity coefficient.
    ln_gamma: tuple[float, ...]
    #: The natural logarithm of each fugacity coefficient, ``ln(gamma_i P0_i / P)``.
    ln_phi: tuple[float, ...]
    #: The pure-component saturation pressure of each component at ``T``.
    p_sat: tuple[Q, ...]
    #: Caveats.
    warnings: tuple[Warning, ...]


@dataclass(frozen=True, slots=True, eq=False)
class PitzerPhaseResult(_HasWarnings):
    """Result of ``eos.pitzer_phase``."""

    #: The activity coefficient of each component.
    gamma: tuple[float, ...]
    #: The natural logarithm of each activity coefficient.
    ln_gamma: tuple[float, ...]
    #: The natural logarithm of each fugacity coefficient, from the three arms the
    #: spec's assumptions name.
    ln_phi: tuple[float, ...]
    #: The Henry coefficient each component's arm read, or zero where it read none.
    henry: tuple[Q, ...]
    #: The infinite-dilution activity coefficient each arm divided by, or one where it
    #: divided by none.
    gamma_inf: tuple[float, ...]
    #: Each component's molality ``n_i / m_water``, in mol/kg of solvent.
    molality: tuple[float, ...]
    #: ``I = 1/2 sum m_i z_i^2``, in mol/kg.
    ionic_strength: float
    #: The Pitzer osmotic coefficient of the water.
    osmotic_coefficient: float
    #: The water activity ``a_w``.
    water_activity: float
    #: Which parameter dataset answered: ``phreeqc`` or ``legacy``.
    dataset: str
    #: Caveats.
    warnings: tuple[Warning, ...]


@dataclass(frozen=True, slots=True, eq=False)
class KentEisenbergPhaseResult(_HasWarnings):
    """Result of ``eos.kent_eisenberg_phase``."""

    #: The activity coefficient of each component, identically one.
    gamma: tuple[float, ...]
    #: The natural logarithm of each activity coefficient, identically zero.
    ln_gamma: tuple[float, ...]
    #: The natural logarithm of each fugacity coefficient.
    ln_phi: tuple[float, ...]
    #: Caveats.
    warnings: tuple[Warning, ...]


@dataclass(frozen=True, slots=True, eq=False)
class DesmukhMatherPhaseResult(_HasWarnings):
    """Result of ``eos.desmukh_mather_phase``."""

    #: The activity coefficient of each component.
    gamma: tuple[float, ...]
    #: The natural logarithm of each activity coefficient.
    ln_gamma: tuple[float, ...]
    #: Each component's molality ``n_i / m_solvent``, in mol/kg of solvent.
    molality: tuple[float, ...]
    #: ``I = 1/2 sum m_i z_i^2``, in mol/kg.
    ionic_strength: float
    #: The mean molar mass of the ``solvent``-reference components, in kg/mol.
    solvent_molar_mass: float
    #: The natural logarithm of each fugacity coefficient.
    ln_phi: tuple[float, ...]
    #: Caveats.
    warnings: tuple[Warning, ...]


@dataclass(frozen=True, slots=True, eq=False)
class GeUniquacPhaseResult(_HasWarnings):
    """Result of ``eos.ge_uniquac_phase``."""

    #: The UNIQUAC activity coefficient of each component.
    gamma: tuple[float, ...]
    #: The natural logarithm of each activity coefficient.
    ln_gamma: tuple[float, ...]
    #: The natural logarithm of each fugacity coefficient, ``ln(gamma_i P0_i / P)``.
    ln_phi: tuple[float, ...]
    #: The pure-component saturation pressure of each component at ``T``.
    p_sat: tuple[Q, ...]
    #: Caveats.
    warnings: tuple[Warning, ...]


@dataclass(frozen=True, slots=True, eq=False)
class GeVanLaarAcidPhaseResult(_HasWarnings):
    """Result of ``eos.ge_van_laar_acid_phase``."""

    #: The Van Laar acid activity coefficient of each component.
    gamma: tuple[float, ...]
    #: The natural logarithm of each activity coefficient.
    ln_gamma: tuple[float, ...]
    #: The natural logarithm of each fugacity coefficient, ``ln(gamma_i P0_i / P)``.
    ln_phi: tuple[float, ...]
    #: The pure-component saturation pressure of each component at ``T``.
    p_sat: tuple[Q, ...]
    #: Caveats.
    warnings: tuple[Warning, ...]


@dataclass(frozen=True, slots=True, eq=False)
class Gerg2008PhaseResult(_HasWarnings):
    """Result of ``eos.gerg2008_phase``.

    The GERG-2008 phase state: the compressibility factor and the Helmholtz property set
    (internal energy, enthalpy, entropy, heat capacities and Gibbs energy).
    """

    #: The compressibility factor ``Z = P v/(R T)``.
    z_factor: float
    #: The internal energy.
    u: Q
    #: The enthalpy.
    h: Q
    #: The entropy.
    s: Q
    #: The isochoric heat capacity.
    cv: Q
    #: The isobaric heat capacity.
    cp: Q
    #: The Gibbs energy, ``h - T s``.
    g: Q
    #: Caveats.
    warnings: tuple[Warning, ...]


@dataclass(frozen=True, slots=True, eq=False)
class HeliumPhaseResult(_HasWarnings):
    """Result of ``eos.helium_phase``.

    The Vega helium phase state: the compressibility factor and the Helmholtz property
    set (internal energy, enthalpy, entropy, heat capacities and Gibbs energy).
    """

    #: The compressibility factor ``Z = P/(rho R T)``.
    z_factor: float
    #: The internal energy.
    u: Q
    #: The enthalpy.
    h: Q
    #: The entropy.
    s: Q
    #: The isochoric heat capacity.
    cv: Q
    #: The isobaric heat capacity.
    cp: Q
    #: The Gibbs energy, ``h - T s``.
    g: Q
    #: Caveats.
    warnings: tuple[Warning, ...]


@dataclass(frozen=True, slots=True, eq=False)
class EffectiveDiffusionResult(_HasWarnings):
    """Result of ``eos.effective_diffusion``."""

    #: One effective coefficient per component, in ``x``'s order, in m**2/s.
    effective_diffusion: tuple[Q, ...]
    #: Caveats.
    warnings: tuple[Warning, ...]


@dataclass(frozen=True, slots=True, eq=False)
class GeFlashResult(_HasWarnings):
    """Result of ``eos.ge_flash``.

    The same shape as :class:`PtFlashResult`, and the same rule about ``beta``: it is
    ``None`` where there is genuinely no vapour fraction - a feed whose K-values are all
    on one side of one, so Rachford-Rice has no root - and a number, possibly outside
    ``[0, 1]``, where the flash converged to one.

    What differs is where the two phases come from. ``ln_phi_liquid`` is whichever
    activity-coefficient phase ``liquid_model`` named, as ``ln(gamma_i P0_i / P)``, and
    ``ln_phi_vapour`` is the cubic's, so the two are not the same kind of number and are
    reported apart for the same reason the phase model reports ``p_sat`` beside
    ``ln_phi``.
    """

    #: The vapour fraction, or ``None`` when there is no vapour fraction to report.
    beta: float | None
    #: Liquid-phase mole fractions.
    x: tuple[float, ...]
    #: Vapour-phase mole fractions.
    y: tuple[float, ...]
    #: ``K_i = y_i / x_i``, the iterate the loop converges on.
    k: tuple[float, ...]
    #: ``ln phi_i`` in the liquid, from the named activity-coefficient phase.
    ln_phi_liquid: tuple[float, ...]
    #: ``ln phi_i`` in the vapour, from the cubic.
    ln_phi_vapour: tuple[float, ...]
    #: The vapour root of the cubic, the largest admissible one. There is no liquid
    #: root: the liquid is not a cubic.
    z_vapour: float
    #: The smallest ``T / Tc_i`` over the components.
    min_t_over_tc: float
    #: What the converged state is.
    phase: Phase
    #: Successive-substitution steps taken.
    iterations: int
    #: ``rms_i |ln K_i - ln K_i_previous|`` at the last step the loop completed, or
    #: ``NaN`` when no step completed.
    residual: float
    #: Caveats.
    warnings: tuple[Warning, ...]
