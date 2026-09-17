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
    #: The larger of the relative volume and internal-energy residuals at the answer.
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
