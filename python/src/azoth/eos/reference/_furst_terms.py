"""The three Helmholtz terms a Furst electrolyte phase adds, and the state they read.

The Python reference, the mirror of ``crates/azoth-eos/src/furst_terms.rs``. The state
builder is here rather than in the model for the reason the Rust module gives: the order is
load-bearing - the packing fraction needs the volume, the dielectric constants need the
packing fraction, ``alphaLR2`` needs the dielectric, and **the shielding parameter needs
``alphaLR2``** - so a ``gamma`` read before the dielectric is a ``gamma`` of the wrong brine.

Two of NeqSim's constants beside the ones ``_furst_dielectric`` carries:
``electronCharge`` is ``1.6021917e-19`` and ``vacumPermittivity`` is ``8.85419e-12``, the
values of their era rather than the current CODATA ones, and ``R`` is ``8.3144621``.
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from typing import NamedTuple

from azoth.core.errors import InvalidInputError, OutOfRangeError
from azoth.eos.reference._furst_dielectric import (
    NEQSIM_AVOGADRO,
    NEQSIM_PI,
    MixingRule,
    packing_fraction,
    packing_fraction_dv,
    packing_fraction_dvdv,
    solvent_dielectric,
)

#: NeqSim's ``electronCharge``.
ELECTRON_CHARGE = 1.6021917e-19

#: NeqSim's ``vacumPermittivity``.
VACUUM_PERMITTIVITY = 8.85419e-12

#: NeqSim's ``R``.
R = 8.3144621


class FurstState:
    """Everything the three terms and their derivatives read at one state.

    A plain class rather than a dataclass, because it is built in two passes: the
    volume-dependent half from ``build_state`` and then the shielding-dependent half from the
    solve, which needs the first half's ``alphaLR2``.
    """

    temperature: float
    molar_volume: float
    moles: float
    packing: float
    ionic_packing: float
    solvent_dielectric: float
    solvent_dielectric_dt: float
    dielectric: float
    dielectric_dt: float
    dielectric_dtdt: float
    dielectric_dv: float
    dielectric_dvdv: float
    dielectric_dtdv: float
    shielding: float
    shielding_dt: float
    xlr: float
    xlr_dt: float
    born_x: float
    w: float
    w_dt: float
    w_dtdt: float
    __slots__ = (
        "born_x",
        "dielectric",
        "dielectric_dt",
        "dielectric_dtdt",
        "dielectric_dtdv",
        "dielectric_dv",
        "dielectric_dvdv",
        "ionic_packing",
        "molar_volume",
        "moles",
        "packing",
        "shielding",
        "shielding_dt",
        "solvent_dielectric",
        "solvent_dielectric_dt",
        "temperature",
        "w",
        "w_dt",
        "w_dtdt",
        "xlr",
        "xlr_dt",
    )


class ComponentState:
    """One component's inputs to the state builder."""

    charge: float
    diameter_m: float
    dielectric: float
    dielectric_dt: float
    dielectric_dtdt: float
    critical_volume: float
    __slots__ = (
        "charge",
        "critical_volume",
        "diameter_m",
        "dielectric",
        "dielectric_dt",
        "dielectric_dtdt",
    )

    def __init__(
        self,
        charge: float,
        diameter_m: float,
        dielectric: float,
        dielectric_dt: float,
        dielectric_dtdt: float,
        critical_volume: float,
    ) -> None:
        self.charge = charge
        self.diameter_m = diameter_m
        self.dielectric = dielectric
        self.dielectric_dt = dielectric_dt
        self.dielectric_dtdt = dielectric_dtdt
        self.critical_volume = critical_volume


def alpha_lr2(dielectric: float, temperature: float) -> float:
    """``alphaLR2 = e^2 N_A/(eps0 eps R T)``, ``eps`` the *phase's* constant."""
    return (
        ELECTRON_CHARGE
        * ELECTRON_CHARGE
        * NEQSIM_AVOGADRO
        / (VACUUM_PERMITTIVITY * dielectric * R * temperature)
    )


def alpha_lr2_dt(state: FurstState) -> float:
    """``dalphaLR2/dT``.

    ``alphaLR2 = scale/(eps T)``, so the derivative is ``-scale/(eps T^2) - scale eps'/(eps^2
    T)``: the second term carries a ``1/T`` that is easy to drop and a central difference of
    ``alphaLR2`` catches immediately.
    """
    scale = ELECTRON_CHARGE * ELECTRON_CHARGE * NEQSIM_AVOGADRO / (VACUUM_PERMITTIVITY * R)
    t, eps = state.temperature, state.dielectric
    return -scale / (eps * t * t) - scale * state.dielectric_dt / (eps * eps * t)


def alpha_lr2_dv(state: FurstState) -> float:
    """``dalphaLR2/dV``."""
    return (
        -ELECTRON_CHARGE
        * ELECTRON_CHARGE
        * NEQSIM_AVOGADRO
        / (VACUUM_PERMITTIVITY * state.dielectric**2 * R * state.temperature)
        * state.dielectric_dv
    )


def born_scale(temperature: float) -> float:
    """``FBorn``'s prefactor ``K(T) = N_A e^2/(4 pi eps0 R T)``."""
    return (
        NEQSIM_AVOGADRO
        * ELECTRON_CHARGE
        * ELECTRON_CHARGE
        / (4.0 * NEQSIM_PI * VACUUM_PERMITTIVITY * R * temperature)
    )


def fsr2(state: FurstState) -> float:
    """The short-range term ``FSR2 = W/(V n (1 - eps))``.

    Raises:
        OutOfRangeError: if ``1 - eps`` is not positive.
    """
    denominator = state.molar_volume * state.moles * (1.0 - state.packing)
    if not math.isfinite(denominator) or denominator <= 0.0:
        raise OutOfRangeError(
            "packing",
            state.packing,
            "the short-range term divides by `V n (1 - eps)` and the excluded volume is "
            "larger than the phase's own, so the term has no value there",
        )
    return state.w / denominator


def t_d_fsr2_dt(state: FurstState) -> float:
    """``T dFSR2/dT``."""
    return (
        state.w_dt * state.temperature / (state.molar_volume * state.moles * (1.0 - state.packing))
    )


def fsr2_dv(state: FurstState) -> float:
    """``dFSR2/dV``, against the total volume in NeqSim's scaling."""
    vn = state.molar_volume * state.moles
    one_minus = 1.0 - state.packing
    fsr2_v = -state.w / (vn * vn * one_minus)
    fsr2_eps = state.w / (vn * one_minus * one_minus)
    eps_dv = -state.packing / vn
    return (fsr2_v + fsr2_eps * eps_dv) * 1.0e-5


def fsr2_dvdv(state: FurstState) -> float:
    """``d^2FSR2/dV^2``."""
    vn = state.molar_volume * state.moles
    one_minus = 1.0 - state.packing
    eps_dv = -state.packing / vn
    eps_dvdv = 2.0 * state.packing / (vn * vn)
    fsr2_vv = 2.0 * state.w / (vn * vn * vn * one_minus)
    fsr2_eps_v = -state.w / (vn * vn * one_minus * one_minus)
    fsr2_eps_eps = 2.0 * state.w / (vn * one_minus * one_minus * one_minus)
    fsr2_eps = state.w / (vn * one_minus * one_minus)
    return (
        fsr2_vv + 2.0 * fsr2_eps_v * eps_dv + fsr2_eps_eps * eps_dv * eps_dv + fsr2_eps * eps_dvdv
    ) * 1.0e-10


def fsr2_dtdv(state: FurstState) -> float:
    """``d^2FSR2/dT dV``."""
    vn = state.molar_volume * state.moles
    one_minus = 1.0 - state.packing
    eps_dv = -state.packing / vn
    fsr2_vw = -1.0 / (vn * vn * one_minus)
    fsr2_eps_w = 1.0 / (vn * one_minus * one_minus)
    return (fsr2_vw * state.w_dt + fsr2_eps_w * eps_dv * state.w_dt) * 1.0e-5


def t2_d2_fsr2_dt2(state: FurstState) -> float:
    """``T^2 d^2FSR2/dT^2``."""
    return (
        state.w_dtdt
        * state.temperature**2
        / (state.molar_volume * state.moles * (1.0 - state.packing))
    )


def flr(state: FurstState) -> float:
    """The MSA long-range term."""
    first = -alpha_lr2(state.dielectric, state.temperature) * state.xlr / (4.0 * NEQSIM_PI)
    second = (
        state.moles * state.molar_volume * state.shielding**3 / (3.0 * NEQSIM_PI * NEQSIM_AVOGADRO)
    )
    return first + second


def t_d_flr_dt(state: FurstState) -> float:
    """``T dFLR/dT``."""
    alpha = alpha_lr2(state.dielectric, state.temperature)
    alpha_dt = alpha_lr2_dt(state)
    term1 = -1.0 / (4.0 * NEQSIM_PI) * (alpha_dt * state.xlr + alpha * state.xlr_dt)
    term2 = (
        state.moles
        * state.molar_volume
        / (3.0 * NEQSIM_PI * NEQSIM_AVOGADRO)
        * 3.0
        * state.shielding**2
        * state.shielding_dt
    )
    return (term1 + term2) * state.temperature


def flr_dv(state: FurstState) -> float:
    """``dFLR/dV``."""
    d_f_d_alpha = -state.xlr / (4.0 * NEQSIM_PI)
    flr_v = state.shielding**3 / (3.0 * NEQSIM_PI * NEQSIM_AVOGADRO)
    return (flr_v + d_f_d_alpha * alpha_lr2_dv(state)) * 1.0e-5


def flr_dvdv(state: FurstState) -> float:
    """``d^2FLR/dV^2``."""
    d_f_d_alpha = -state.xlr / (4.0 * NEQSIM_PI)
    e = ELECTRON_CHARGE
    eps, t = state.dielectric, state.temperature
    alpha_dvdv = (
        2.0
        * e
        * e
        * NEQSIM_AVOGADRO
        / (VACUUM_PERMITTIVITY * eps**3 * R * t)
        * state.dielectric_dv
        * state.dielectric_dv
        - e
        * e
        * NEQSIM_AVOGADRO
        / (VACUUM_PERMITTIVITY * eps * eps * R * t)
        * state.dielectric_dvdv
    )
    return d_f_d_alpha * alpha_dvdv * 1.0e-10


def flr_dtdv(state: FurstState) -> float:
    """``d^2FLR/dT dV``."""
    d_f_d_alpha = -state.xlr / (4.0 * NEQSIM_PI)
    e = ELECTRON_CHARGE
    eps, t = state.dielectric, state.temperature
    alpha_dtdv = (
        e
        * e
        * NEQSIM_AVOGADRO
        / (VACUUM_PERMITTIVITY * eps * eps * R * t * t)
        * state.dielectric_dv
        + 2.0
        * e
        * e
        * NEQSIM_AVOGADRO
        / (VACUUM_PERMITTIVITY * eps**3 * R * t)
        * state.dielectric_dt
        * state.dielectric_dv
        - e
        * e
        * NEQSIM_AVOGADRO
        / (VACUUM_PERMITTIVITY * eps * eps * R * t)
        * state.dielectric_dtdv
    )
    return d_f_d_alpha * alpha_dtdv * 1.0e-5


def t2_d2_flr_dt2(state: FurstState) -> float:
    """``T^2 d^2FLR/dT^2``.

    **NeqSim's own simplification**, not the second derivative: its ``dFLRdTdT`` keeps
    ``-1/(4 pi)(alphaLR'' XLR + 2 alphaLR' XLR')`` and the whole of the ``gamma`` chain is
    missing, which the source says. Reproduced, because completing it would be a different
    model with a divergence nothing would catch.
    """
    e, e0, na = ELECTRON_CHARGE, VACUUM_PERMITTIVITY, NEQSIM_AVOGADRO
    t, eps = state.temperature, state.dielectric
    # NeqSim's `alphaLRdTdT`. Its `2 scale eps'/(eps^2 T^2)` term is written upstream as two
    # identical lines in a row; this is the analytic form, which is the same expression.
    scale = e * e * na / (e0 * R)
    alpha_dtdt = (
        2.0 * scale / (eps * t**3)
        + 2.0 * scale * state.dielectric_dt / (eps * eps * t * t)
        - scale * state.dielectric_dtdt / (eps * eps * t)
        + 2.0 * scale * state.dielectric_dt**2 / (eps**3 * t)
    )
    alpha_dt = alpha_lr2_dt(state)
    term1 = -1.0 / (4.0 * NEQSIM_PI) * alpha_dtdt * state.xlr
    cross = -2.0 / (4.0 * NEQSIM_PI) * alpha_dt * state.xlr_dt
    return (term1 + cross) * t * t


def fborn(state: FurstState) -> float:
    """The Born solvation term."""
    return born_scale(state.temperature) * (1.0 / state.solvent_dielectric - 1.0) * state.born_x


def t_d_fborn_dt(state: FurstState) -> float:
    """``T dFBorn/dT``."""
    k, t = born_scale(state.temperature), state.temperature
    d = 1.0 / state.solvent_dielectric - 1.0
    f_born_t = -(k / t) * d * state.born_x
    f_born_d = -k / state.solvent_dielectric**2 * state.born_x
    return (f_born_t + f_born_d * state.solvent_dielectric_dt) * t


def t2_d2_fborn_dt2(state: FurstState) -> float:
    """``T^2 d^2FBorn/dT^2``.

    **Also simplified upstream**: NeqSim keeps ``FBornTT + FBornTD eps'`` and drops both the
    ``eps''`` and the ``(eps')^2`` terms, so this is not the second derivative of the
    expression above it.
    """
    k, t = born_scale(state.temperature), state.temperature
    d = 1.0 / state.solvent_dielectric - 1.0
    f_born_tt = 2.0 * k / (t * t) * d * state.born_x
    f_born_td = k / t / state.solvent_dielectric**2 * state.born_x
    return (f_born_tt + f_born_td * state.solvent_dielectric_dt) * t * t


class FurstSolution(NamedTuple):
    """The Fürst term's three outputs at one volume.

    A record rather than a mapping, because the three always travel together and a caller
    that mistyped a key would get an ``Any`` and not an error.
    """

    #: The three terms' contribution to ``A^R/(R T)``.
    helmholtz_rt: float
    #: Its derivative with respect to the **SI** volume - the pressure the term carries,
    #: divided by ``-R T``. What the root residual takes.
    helmholtz_rt_dv: float
    #: The three terms' contribution to ``dFdN_i``, one per component.
    ln_phi: list[float]


class CompositionContribution:
    """One component's electrolyte contribution to ``dFdN_i``, by term."""

    __slots__ = ("born", "long_range", "short_range")

    def __init__(self, short_range: float, long_range: float, born: float) -> None:
        self.short_range = short_range
        self.long_range = long_range
        self.born = born

    def total(self) -> float:
        """The three added, which is what the fugacity coefficient takes."""
        return self.short_range + self.long_range + self.born


class ComponentDerivatives:
    """One component's inputs to the composition derivatives."""

    __slots__ = ("charge", "diameter_m", "dielectric", "w_i")

    def __init__(self, charge: float, diameter_m: float, dielectric: float, w_i: float) -> None:
        self.charge = charge
        self.diameter_m = diameter_m
        self.dielectric = dielectric
        #: NeqSim's `calcWi`, which is **`-2 sum_j n_j Wij(i, j, T)`** and not the row sum:
        #: `calcW`, `calcWij` and `calcWi` each return their sum negated and doubled.
        self.w_i = w_i


def ln_phi_contributions(
    state: FurstState,
    components: Sequence[ComponentDerivatives],
    mole_numbers: Sequence[float],
    mod2004: bool = False,
) -> list[CompositionContribution]:
    """The three electrolyte contributions to ``dFdN_i``.

    ``ComponentEos.fugcoef`` is ``exp(dFdN - ln(PV/RT))`` and this library's cubic ``ln phi``
    is already ``dFdN - ln Z`` for its part, so these are the *changes*.

    Raises:
        InvalidInputError: if the two sequences are not one per component.
        OutOfRangeError: if no neutral component carries any moles, which the solvent
            dielectric constant's composition derivative divides by.
    """
    if len(components) != len(mole_numbers):
        raise InvalidInputError(
            "components",
            f"{len(components)} components and {len(mole_numbers)} mole numbers; the "
            f"composition derivative is one per component",
        )
    neutral_moles = sum(n for c, n in zip(components, mole_numbers, strict=True) if c.charge == 0.0)
    if not math.isfinite(neutral_moles) or neutral_moles <= 0.0:
        raise OutOfRangeError(
            "mole_numbers",
            neutral_moles,
            "the solvent dielectric constant's composition derivative divides by the "
            "neutral components' total moles, and a phase with none has no such derivative",
        )

    vn = state.molar_volume * state.moles
    one_minus = 1.0 - state.packing
    fsr2_eps = state.w / (vn * one_minus * one_minus)
    fsr2_w = 1.0 / (vn * one_minus)
    flr_xlr = -alpha_lr2(state.dielectric, state.temperature) / (4.0 * NEQSIM_PI)
    d_f_d_alpha = -state.xlr / (4.0 * NEQSIM_PI)
    k = born_scale(state.temperature)
    f_born_x = k * (1.0 / state.solvent_dielectric - 1.0)
    f_born_d = -k / state.solvent_dielectric**2 * state.born_x
    eps_ionic_half = 1.0 + state.ionic_packing / 2.0

    out: list[CompositionContribution] = []
    for component in components:
        scale = NEQSIM_AVOGADRO * NEQSIM_PI / 6.0 * component.diameter_m**3 / vn
        eps_ionic_i = 0.0 if component.charge == 0.0 else scale
        # **Zero in the 2004 revision**, whose `calcSolventdiElectricdn` returns `0.0` with
        # its body commented out.
        solvent_dn = (
            0.0
            if mod2004 or component.charge != 0.0
            else (component.dielectric - state.solvent_dielectric) / neutral_moles
        )
        x = (1.0 - state.ionic_packing) / eps_ionic_half
        y = state.solvent_dielectric - 1.0
        d_x = eps_ionic_i * -1.5 / eps_ionic_half**2
        dielectric_dn = solvent_dn * x + y * d_x
        alpha_i = (
            -ELECTRON_CHARGE
            * ELECTRON_CHARGE
            * NEQSIM_AVOGADRO
            / (VACUUM_PERMITTIVITY * state.dielectric**2 * R * state.temperature)
            * dielectric_dn
        )
        xlr_i = (
            component.charge**2 * state.shielding / (1.0 + state.shielding * component.diameter_m)
        )
        born_i = component.charge**2 / component.diameter_m if component.diameter_m > 0.0 else 0.0
        # **`FBornD` is *added* in the 2004 revision, not weighted by the solvent's
        # composition derivative.** The variant's `dFBorndN` is `FBornX XBorni + FBornD`
        # where the base's is `+ FBornD solventdiElectricdn`, and since that derivative is
        # zero there the base would give `FBornX XBorni` alone. So the variant adds a
        # component-independent term to every `ln phi`.
        born = f_born_x * born_i + (f_born_d if mod2004 else f_born_d * solvent_dn)
        out.append(
            CompositionContribution(
                short_range=fsr2_eps * scale + fsr2_w * component.w_i,
                long_range=flr_xlr * xlr_i + d_f_d_alpha * alpha_i,
                born=born,
            )
        )
    return out


def shielding_parameter(
    state: FurstState, components: Sequence[ComponentState], mole_numbers: Sequence[float]
) -> float:
    """``calcShieldingParameter``: a damped Newton solve for ``gamma``.

    From ``1e10``, ``gamma -= 0.8 f/f'``, with a 1000-iteration cap and a **three-iteration
    floor**. The floor matters: a phase carrying no ion returns exactly zero, and one
    carrying ions at ``1e-43`` mol returns a residue of the solve rather than a number.
    """
    alpha = alpha_lr2(state.dielectric, state.temperature)
    v_total = state.molar_volume * state.moles
    gamma = 1.0e10
    iterations = 0
    while True:
        iterations += 1
        gamma_old = gamma
        f = 4.0 * gamma * gamma / NEQSIM_AVOGADRO
        df = 8.0 * gamma / NEQSIM_AVOGADRO
        ions = 0
        for component, moles in zip(components, mole_numbers, strict=True):
            if component.charge == 0.0:
                continue
            ions += 1
            sigma = component.diameter_m
            denominator = 1.0 + gamma * sigma
            f -= (
                alpha
                * moles
                / v_total
                * (component.charge / denominator)
                * (component.charge / denominator)
            )
            df += 2.0 * alpha * moles / v_total * component.charge**2 * sigma / denominator**3
        gamma = gamma_old - 0.8 * f / df if ions > 0 else 0.0
        if not ((abs(f) > 1.0e-10 and iterations < 1000) or iterations < 3):
            break
    return gamma


def shielding_parameter_dt(
    state: FurstState, components: Sequence[ComponentState], mole_numbers: Sequence[float]
) -> float:
    """``calcShieldingParameterdT``, by implicit differentiation."""
    if state.shielding < 1.0e-10:
        return 0.0
    alpha = alpha_lr2(state.dielectric, state.temperature)
    alpha_dt = alpha_lr2_dt(state)
    v_total = state.molar_volume * state.moles
    dfdgamma = 8.0 * state.shielding / NEQSIM_AVOGADRO
    total = 0.0
    for component, moles in zip(components, mole_numbers, strict=True):
        if component.charge == 0.0:
            continue
        sigma = component.diameter_m
        denominator = 1.0 + state.shielding * sigma
        dfdgamma += 2.0 * alpha * moles / v_total * component.charge**2 * sigma / denominator**3
        total += moles / v_total * component.charge**2 / denominator**2
    dfdt = -alpha_dt * total
    if abs(dfdgamma) < 1.0e-50:
        return 0.0
    return -dfdt / dfdgamma


def xlr(
    state: FurstState, components: Sequence[ComponentState], mole_numbers: Sequence[float]
) -> float:
    """``calcXLR`` over the ions."""
    return sum(
        moles
        * component.charge**2
        * state.shielding
        / (1.0 + state.shielding * component.diameter_m)
        for component, moles in zip(components, mole_numbers, strict=True)
        if component.charge != 0.0
    )


def xlr_dt(
    state: FurstState, components: Sequence[ComponentState], mole_numbers: Sequence[float]
) -> float:
    """``calcXLRdT``, whose ``d/dT`` carries only the shielding parameter's derivative."""
    total = 0.0
    for component, moles in zip(components, mole_numbers, strict=True):
        if component.charge == 0.0:
            continue
        denominator = 1.0 + state.shielding * component.diameter_m
        total += moles * component.charge**2 * state.shielding_dt / (denominator * denominator)
    return total


def born_x(components: Sequence[ComponentState], mole_numbers: Sequence[float]) -> float:
    """``calcBornX = sum_i n_i z_i^2/sigma_i`` over every component with a diameter."""
    return sum(
        moles * component.charge**2 / component.diameter_m
        for component, moles in zip(components, mole_numbers, strict=True)
        if component.diameter_m > 0.0
    )


def build_state(
    temperature: float,
    molar_volume: float,
    mole_numbers: Sequence[float],
    components: Sequence[ComponentState],
    rule: MixingRule,
    short_range_w: float,
    short_range_w_dt: float,
    short_range_w_dtdt: float,
) -> FurstState:
    """Build the state, in ``volInit``'s own order.

    Raises:
        InvalidInputError: if the sequences are not one per component.
        OutOfRangeError: from the packing fractions, or if the phase has no solvent.
    """
    n = len(components)
    if len(mole_numbers) != n:
        raise InvalidInputError(
            "components",
            f"{n} components and {len(mole_numbers)} mole numbers; the state is one per component",
        )
    is_ion = [c.charge != 0.0 for c in components]
    diameters = [c.diameter_m for c in components]
    total_moles = sum(mole_numbers)

    packing = packing_fraction(
        diameters, list(mole_numbers), is_ion, total_moles, molar_volume, ions_only=False
    )
    ionic_packing = packing_fraction(
        diameters, list(mole_numbers), is_ion, total_moles, molar_volume, ions_only=True
    )

    eps_i = [c.dielectric for c in components]
    eps_dt = [c.dielectric_dt for c in components]
    eps_dtdt = [c.dielectric_dtdt for c in components]
    critical_volumes = [c.critical_volume for c in components]
    solvent = solvent_dielectric(rule, list(mole_numbers), eps_i, is_ion, critical_volumes)
    solvent_dt = solvent_dielectric(rule, list(mole_numbers), eps_dt, is_ion, critical_volumes)
    solvent_dtdt = solvent_dielectric(rule, list(mole_numbers), eps_dtdt, is_ion, critical_volumes)

    x = (1.0 - ionic_packing) / (1.0 + ionic_packing / 2.0)
    y = solvent - 1.0
    ionic_dv = packing_fraction_dv(ionic_packing, total_moles, molar_volume)
    ionic_dvdv = packing_fraction_dvdv(ionic_packing, total_moles, molar_volume)
    d_x_dv = ionic_dv * -1.5 / (1.0 + ionic_packing / 2.0) ** 2
    d_x_dvdv = (
        ionic_dvdv * -1.5 / (1.0 + ionic_packing / 2.0) ** 2
        + ionic_dv * ionic_dv * 1.5 / (1.0 + ionic_packing / 2.0) ** 3
    )

    state = FurstState()
    state.temperature = temperature
    state.molar_volume = molar_volume
    state.moles = total_moles
    state.packing = packing
    state.ionic_packing = ionic_packing
    state.solvent_dielectric = solvent
    state.solvent_dielectric_dt = solvent_dt
    state.dielectric = 1.0 + y * x
    state.dielectric_dt = solvent_dt * x
    state.dielectric_dtdt = solvent_dtdt * x
    state.dielectric_dv = y * d_x_dv
    state.dielectric_dvdv = y * d_x_dvdv
    state.dielectric_dtdv = solvent_dt * d_x_dv
    state.shielding = 0.0
    state.shielding_dt = 0.0
    state.xlr = 0.0
    state.xlr_dt = 0.0
    state.born_x = 0.0
    state.w = short_range_w
    state.w_dt = short_range_w_dt
    state.w_dtdt = short_range_w_dtdt

    state.shielding = shielding_parameter(state, components, mole_numbers)
    state.shielding_dt = shielding_parameter_dt(state, components, mole_numbers)
    state.xlr = xlr(state, components, mole_numbers)
    state.xlr_dt = xlr_dt(state, components, mole_numbers)
    state.born_x = born_x(components, mole_numbers)
    return state
