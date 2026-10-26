"""``ControlValveSizing``, the ``ValveMechanicalDesign`` default: the IEC 60534 simplified sizing.

The Python twin of ``crates/azoth-process/src/valve.rs``. Every step mirrors the Rust module, and
**the two halves are held to the same capture**: the Rust module's own tests and
``python/tests/test_cross_impl.py`` compare the rows of
``validation/neqsim/captures/process_valve_sizing.tsv``.

**The default sizing standard and not the IEC class.** ``setValveSizingStandard`` selects the
others and this module is the default: ``calcKv``, its inverse flow and its outlet-pressure solve,
on the liquid and gas branches, with the opening characteristic and the fouling fraction.

**The gas formula is the class's, unit mixing included.** ``Q`` at standard conditions (273.15 K,
101.325 kPa) comes from the actual volumetric flow through ``P1``, ``T`` and ``Z``, then
``Kv = Q_std/(N9 P1 Y) sqrt(M T Z / x)`` with ``N9 = 24.6``, ``M`` in g/mol, ``T`` in kelvin and
``P1`` in kPa. The liquid form is ``Kv = Q/sqrt(dP/SG)`` with ``SG = rho/1000``. Both divide by
the opening factor, so the reported ``Kv`` is the rated value at 100% opening.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from azoth.core.units import quantity
from azoth.eos import components as _components
from azoth.eos.reference.molar_enthalpy_entropy import molar_enthalpy_entropy
from azoth.eos.reference.pr_molar_volume import pr_molar_volume
from azoth.eos.reference.pt_flash import pt_flash
from azoth.process.reference import capacity

#: ``Cv = Kv * 1.156``.
KV_TO_CV_FACTOR = 1.156
#: The IEC 60534 constant for standard m3/h, kPa, K and g/mol.
N9 = 24.6
#: The standard temperature the gas flow is quoted at, K.
STANDARD_TEMPERATURE_K = 273.15
#: The standard pressure the gas flow is quoted at, kPa.
STANDARD_PRESSURE_KPA = 101.325
#: ``xT``'s own initialiser.
DEFAULT_XT = 0.137
#: ``Cg = Cv * 1360``.
CG_PER_CV = 1360.0
#: The bisection's ``Kv`` tolerance and cap, the class's own.
BISECTION_TOLERANCE = 1.0e-6
MAX_BISECTION_ITERATIONS = 100
#: ``ThermodynamicConstantsInterface.R``, J/(mol*K).
R_J_PER_MOL_K = 8.314
#: ``ThrottlingValve.percentValveOpening``'s initialiser.
DEFAULT_OPENING_PERCENT = 100.0
#: ``EqualPercentageCharacteristic``'s default rangeability.
DEFAULT_RANGEABILITY = 50.0
#: ``ModifiedParabolicCharacteristic``'s default exponent.
DEFAULT_PARABOLIC_EXPONENT = 2.0

#: The class's four ``ValveCharacteristic`` implementations.
CHARACTERISTICS = ("linear", "equal_percentage", "quick_opening", "modified_parabolic")


def opening_factor(characteristic: str, percent_opening: float) -> float:
    """``getOpeningFactor(percentOpening)``: the fraction of the rated ``Kv``."""
    opening = min(100.0, max(0.0, percent_opening))
    x = opening / 100.0
    if characteristic == "linear":
        return x
    if characteristic == "equal_percentage":
        if opening <= 0.0:
            return 1.0 / DEFAULT_RANGEABILITY
        if opening >= 100.0:
            return 1.0
        return float(DEFAULT_RANGEABILITY ** (x - 1.0))
    if characteristic == "quick_opening":
        return math.sqrt(x)
    return float(x**DEFAULT_PARABOLIC_EXPONENT)


@dataclass(frozen=True, slots=True)
class ValveSizingGeometry:
    """The sizing's configuration, as ``ValveMechanicalDesign`` and ``ThrottlingValve`` hold it."""

    kv: float | None = None
    opening_percent: float = DEFAULT_OPENING_PERCENT
    fouling_fraction: float = 0.0
    characteristic: str = "linear"
    xt: float = DEFAULT_XT
    allow_choked: bool = False


@dataclass(frozen=True, slots=True)
class ValveSizingReport:
    """What ``calcValveSize`` and its inverse solve report."""

    kv: float
    cv: float
    cg: float
    choked: bool
    gas_valve: bool
    calculated_molar_flow: float
    calculated_outlet_pressure: float


@dataclass(frozen=True, slots=True)
class _Reads:
    p1_bara: float
    t_k: float
    flow_m3hr: float
    density_kg_m3: float
    z: float
    gamma: float
    molar_mass_gr_mol: float


def gas_valve(components: list[str], t: float, p: float, z: list[float]) -> bool:
    """``isGasValve()``: a gas phase carrying more than half the volume."""
    mixture, labelled, _shares = capacity._labelled_phases(components, t, p, z)
    flash = pt_flash(mixture, quantity(t, "K"), quantity(p, "Pa"), list(z))
    if flash.phase == "all_vapour":
        return True
    if flash.phase in ("all_liquid", "trivial"):
        return False
    y = flash.vapour_fraction if flash.vapour_fraction is not None else 0.0
    # The first labelled phase is the gas and the second the liquid, as the class's own order.
    gas_density = capacity._mass_density(mixture, t, p, labelled[0][1], labelled[0][2])
    liquid_density = capacity._mass_density(mixture, t, p, labelled[1][1], labelled[1][2])
    gas_mass = capacity._molar_mass(mixture, labelled[0][1])
    liquid_mass = capacity._molar_mass(mixture, labelled[1][1])
    gas_volume = y * gas_mass / gas_density
    liquid_volume = (1.0 - y) * liquid_mass / liquid_density
    return gas_volume / (gas_volume + liquid_volume) > 0.5


def _reads(components: list[str], t: float, p: float, z: list[float], n: float) -> _Reads:
    density = capacity.system_mass_density(components, t, p, z)
    flow_m3s = capacity.system_volumetric_flow(components, t, p, z, n)
    mixture, ideal_gas = _components.mixture_of(components, eos="pr")
    molar_mass = capacity._molar_mass(mixture, list(z))
    flash = pt_flash(mixture, quantity(t, "K"), quantity(p, "Pa"), list(z))
    if flash.phase == "all_vapour":
        root = flash.z_vapour
    elif flash.phase in ("all_liquid", "trivial"):
        root = flash.z_liquid
    else:
        root = flash.z_liquid if not gas_valve(components, t, p, z) else flash.z_vapour
    volume = pr_molar_volume(root, quantity(t, "K"), quantity(p, "Pa")).v.to("m**3/mol").magnitude
    state = molar_enthalpy_entropy(
        mixture, ideal_gas, quantity(t, "K"), quantity(p, "Pa"), list(z), root
    )
    cp = float(state.cp.to_base_units().magnitude)
    return _Reads(
        p1_bara=p / 1.0e5,
        t_k=t,
        flow_m3hr=flow_m3s * 3600.0,
        density_kg_m3=density,
        z=p * volume / (R_J_PER_MOL_K * t),
        gamma=cp / (cp - R_J_PER_MOL_K),
        molar_mass_gr_mol=molar_mass * 1000.0,
    )


def adjusted_kv(rated_kv: float, geometry: ValveSizingGeometry) -> float:
    """``getActualKv`` then ``adjustKv``."""
    return (
        rated_kv
        * opening_factor(geometry.characteristic, geometry.opening_percent)
        * (1.0 - geometry.fouling_fraction)
    )


def calc_kv(
    components: list[str],
    t: float,
    p: float,
    z: list[float],
    n: float,
    outlet_bara: float,
    geometry: ValveSizingGeometry,
) -> float:
    """``calcKv``: the rated ``Kv`` the current flow and drop require."""
    read = _reads(components, t, p, z, n)
    gas = gas_valve(components, t, p, z)
    delta_p = read.p1_bara - outlet_bara
    if delta_p <= 0.0:
        from azoth.core.errors import InvalidInputError

        raise InvalidInputError(
            "outlet_pressure",
            f"sizing needs a positive drop, and {outlet_bara} bara against {read.p1_bara} bara "
            f"is {delta_p} bar",
        )
    factor = opening_factor(geometry.characteristic, geometry.opening_percent)
    if factor <= 0.0:
        from azoth.core.errors import InvalidInputError

        raise InvalidInputError(
            "valve_opening",
            f"an opening of {geometry.opening_percent}% gives an opening factor of {factor}, "
            f"and `calcKv` divides by it",
        )
    if not gas:
        sg = read.density_kg_m3 / 1000.0
        kv = read.flow_m3hr / math.sqrt(delta_p / sg)
    else:
        kv = _gas_kv(read, delta_p, geometry)
    return kv / factor


def _gas_kv(read: _Reads, delta_p: float, geometry: ValveSizingGeometry) -> float:
    p1_kpa = read.p1_bara * 100.0
    x = delta_p / read.p1_bara
    f_gamma = read.gamma / 1.40
    choked = x >= f_gamma * geometry.xt
    x_effective = f_gamma * geometry.xt if choked and geometry.allow_choked else x
    y = max(1.0 - x_effective / (3.0 * f_gamma * geometry.xt), 2.0 / 3.0)
    standard_flow = (
        read.flow_m3hr
        * (p1_kpa / STANDARD_PRESSURE_KPA)
        * (STANDARD_TEMPERATURE_K / read.t_k)
        / read.z
    )
    x_scale = geometry.xt * f_gamma if choked and geometry.allow_choked else x
    return (
        standard_flow
        / (N9 * p1_kpa * y)
        * math.sqrt(read.molar_mass_gr_mol * read.t_k * read.z / x_scale)
    )


def calculated_molar_flow(
    components: list[str],
    t: float,
    p: float,
    z: list[float],
    n: float,
    outlet_bara: float,
    adjusted: float,
    geometry: ValveSizingGeometry,
) -> float:
    """``calculateMolarFlow()``: the mole flow the adjusted ``Kv`` and the drop give."""
    read = _reads(components, t, p, z, n)
    delta_p = read.p1_bara - outlet_bara
    if delta_p <= 0.0:
        return 0.0
    if not gas_valve(components, t, p, z):
        sg = read.density_kg_m3 / 1000.0
        flow_m3s = adjusted * math.sqrt(delta_p / sg) / 3600.0
    else:
        p1_kpa = read.p1_bara * 100.0
        x = delta_p / read.p1_bara
        f_gamma = read.gamma / 1.40
        choked = x >= f_gamma * geometry.xt
        x_effective = f_gamma * geometry.xt if choked and geometry.allow_choked else x
        y = max(1.0 - x_effective / (3.0 * f_gamma * geometry.xt), 2.0 / 3.0)
        x_scale = geometry.xt * f_gamma if choked and geometry.allow_choked else x
        standard_flow = (
            adjusted
            * N9
            * p1_kpa
            * y
            / math.sqrt(read.molar_mass_gr_mol * read.t_k * read.z / x_scale)
        )
        actual = (
            standard_flow
            * (STANDARD_PRESSURE_KPA / p1_kpa)
            * (read.t_k / STANDARD_TEMPERATURE_K)
            * read.z
        )
        flow_m3s = actual / 3600.0
    return flow_m3s * read.density_kg_m3 / (read.molar_mass_gr_mol / 1000.0)


def calculated_outlet_pressure(
    components: list[str],
    t: float,
    p: float,
    z: list[float],
    n: float,
    adjusted: float,
    geometry: ValveSizingGeometry,
) -> float:
    """``findOutletPressureForFixedKv``: the outlet pressure an adjusted ``Kv`` gives, bara."""
    from azoth.core.errors import InvalidInputError

    if adjusted <= 0.0:
        raise InvalidInputError("kv", f"an adjusted Kv of {adjusted} gives no outlet pressure")
    read = _reads(components, t, p, z, n)
    if not gas_valve(components, t, p, z):
        sg = read.density_kg_m3 / 1000.0
        return read.p1_bara - (read.flow_m3hr / adjusted) ** 2 * sg
    p1_kpa = read.p1_bara * 100.0
    standard_flow = (
        read.flow_m3hr
        * (p1_kpa / STANDARD_PRESSURE_KPA)
        * (STANDARD_TEMPERATURE_K / read.t_k)
        / read.z
    )
    f_gamma = read.gamma / 1.40
    low = 0.1
    high = read.p1_bara - 0.001
    mid = (low + high) / 2.0
    for _ in range(MAX_BISECTION_ITERATIONS):
        mid = (low + high) / 2.0
        p2_kpa = mid * 100.0
        x = (p1_kpa - p2_kpa) / p1_kpa
        choked = x >= f_gamma * geometry.xt
        x_effective = f_gamma * geometry.xt if choked and geometry.allow_choked else x
        y = max(1.0 - x_effective / (3.0 * f_gamma * geometry.xt), 2.0 / 3.0)
        x_scale = geometry.xt * f_gamma if choked and geometry.allow_choked else x
        kv = (
            standard_flow
            / (N9 * p1_kpa * y)
            * math.sqrt(read.molar_mass_gr_mol * read.t_k * read.z / x_scale)
        )
        if abs(kv - adjusted) / adjusted < BISECTION_TOLERANCE:
            break
        if kv > adjusted:
            high = mid
        else:
            low = mid
    return mid


def sizing_report(
    components: list[str],
    t: float,
    p: float,
    z: list[float],
    n: float,
    outlet_bara: float,
    geometry: ValveSizingGeometry,
) -> ValveSizingReport:
    """The whole report: the sized (or stated) ``Kv`` and the inverse solves."""
    gas = gas_valve(components, t, p, z)
    read = _reads(components, t, p, z, n)
    if geometry.kv is not None:
        kv = geometry.kv
    elif read.p1_bara - outlet_bara > 0.0:
        kv = calc_kv(components, t, p, z, n, outlet_bara, geometry)
    else:
        # **A non-positive drop has no sizing**, and ``calcKv`` answers ``NaN`` there; the
        # ``acceptNegativeDP`` row's stream is still a state.
        kv = math.nan
    if not math.isfinite(kv):
        return ValveSizingReport(
            kv=kv,
            cv=math.nan,
            cg=math.nan,
            choked=False,
            gas_valve=gas,
            calculated_molar_flow=0.0,
            calculated_outlet_pressure=math.nan,
        )
    adjusted = adjusted_kv(kv, geometry)
    choked = (
        gas and (read.p1_bara - outlet_bara) / read.p1_bara >= (read.gamma / 1.40) * geometry.xt
    )
    return ValveSizingReport(
        kv=kv,
        cv=kv * KV_TO_CV_FACTOR,
        cg=kv * KV_TO_CV_FACTOR * CG_PER_CV,
        choked=choked,
        gas_valve=gas,
        calculated_molar_flow=calculated_molar_flow(
            components, t, p, z, n, outlet_bara, adjusted, geometry
        ),
        calculated_outlet_pressure=calculated_outlet_pressure(
            components, t, p, z, n, adjusted, geometry
        ),
    )
