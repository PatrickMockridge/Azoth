"""The mixture arithmetic the `eos` models share.

A private module, and it exists for one reason: the mixture fugacity coefficient is
the only calculation in this library that no registered spec covers - ``eos.pr_departure``
is the pure-component form and the registry is scalar, so it has no composition
vector to hang a mixture version on. Everything that needs it therefore has to use
*the same* function, or the layer acquires a second implementation of the one thing
it cannot check against a spec.

The models that need it are ``eos.pt_flash``, which holds two compositions and
solves for the split, and ``eos.bubble_pressure`` / ``eos.dew_pressure``, which hold
one and solve for the pressure at which the other appears.

What keeps it honest is not an assertion but a *reduction*: at ``N = 1`` the cross-sum
factor collapses to 1 and :func:`phase_state` returns exactly ``eos.pr_departure``'s
``ln phi``, and at ``N = 2`` :func:`mixture_parameters` reproduces
``eos.vdw1f_mix_binary``. Both are asserted in both languages. The
``eos.vdw1f_mix_binary`` spec's notes record why each holds to a couple of ulps rather
than bit-identically.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import TYPE_CHECKING, Any, NamedTuple

from azoth.core.errors import InvalidInputError, OutOfRangeError
from azoth.core.result import Phase
from azoth.core.warnings import Warning, WarningCode
from azoth.eos.alpha_term import (
    Danesh,
    Delft1998,
    Gassem2001,
    MatCop,
    MatCopFallback,
    Mollerup,
    RkAlpha,
    Schwartzentruber,
    Soave,
    SoreideWhitsonWater,
    TwuCoon,
    matcop_kappa,
    umr_kappa,
)
from azoth.eos.cubic import PR, Cubic
from azoth.eos.mixture import Mixture
from azoth.eos.reference._association import (
    NON_ASSOCIATING,
    Association,
    AssociationComponent,
    R,
    SiteDerivatives,
    SiteScheme,
    family_of,
)
from azoth.eos.reference._furst_dielectric import MixingRule as _FurstRule
from azoth.eos.reference._furst_dielectric import (
    component_dielectric as _furst_component_dielectric,
)
from azoth.eos.reference._furst_dielectric import (
    component_dielectric_dt as _furst_component_dielectric_dt,
)
from azoth.eos.reference._furst_dielectric import (
    component_dielectric_dtdt as _furst_component_dielectric_dtdt,
)
from azoth.eos.reference._furst_terms import (
    ComponentDerivatives as _FurstComponentDerivatives,
)
from azoth.eos.reference._furst_terms import (
    ComponentState as _FurstComponentState,
)
from azoth.eos.reference._furst_terms import (
    FurstSolution as _FurstSolution,
)
from azoth.eos.reference._furst_terms import (
    build_state as _furst_build_state,
)
from azoth.eos.reference._furst_terms import (
    fborn as _furst_fborn,
)
from azoth.eos.reference._furst_terms import (
    flr as _furst_flr,
)
from azoth.eos.reference._furst_terms import (
    flr_dv as _furst_flr_dv,
)
from azoth.eos.reference._furst_terms import (
    fsr2 as _furst_fsr2,
)
from azoth.eos.reference._furst_terms import (
    fsr2_dv as _furst_fsr2_dv,
)
from azoth.eos.reference._furst_terms import (
    ln_phi_contributions as _furst_ln_phi_contributions,
)
from azoth.eos.reference._hv_ge import hv_ln_gamma as _hv_ln_gamma
from azoth.eos.reference.pr78_kappa import pr78_kappa
from azoth.eos.reference.pr_alpha_ab import pr_alpha_ab
from azoth.eos.reference.pr_kappa import pr_kappa
from azoth.eos.reference.pr_z_factor import pr_z_factor
from azoth.eos.reference.rk_alpha_ab import rk_alpha_ab
from azoth.eos.reference.srk_alpha_ab import srk_alpha_ab
from azoth.eos.reference.srk_kappa import srk_kappa
from azoth.eos.reference.srk_z_factor import srk_z_factor
from azoth.eos.reference.twu_kappa import twu_kappa

if TYPE_CHECKING:
    from azoth.eos.components import AssociationParameters
    from azoth.eos.mixture import Component

#: NeqSim's ``R``, which the Fürst electrolyte's electrostatics are built on - and it is
#: ``8.3144621``, not the current CODATA value. The reference module carries its own copy
#: because the ion substitution needs it before any term is evaluated.
R_NEQSIM = 8.3144621

#: Wilson's constant. Some sources print 5.37 and the paper is dated 1968 in some
#: and 1969 in others; the discrepancy is recorded in the model specs' references
#: rather than resolved, because a reader meeting the other value needs to know it
#: is the same correlation and not a correction.
WILSON_CONSTANT = 5.373


def wilson_saturation_pressures(components: tuple[Any, ...], temperature: float) -> list[float]:
    """Wilson's saturation-pressure estimate for each component.

    Extrapolated without complaint for a component above its critical temperature,
    where it is a large number with no physical meaning - which is what a starting
    guess needs and why it is not reported to the caller.
    """
    return [
        component.Pc.to_base_units().magnitude
        * math.exp(
            WILSON_CONSTANT
            * (1.0 + component.omega)
            * (1.0 - component.Tc.to_base_units().magnitude / temperature)
        )
        for component in components
    ]


def wilson_k(mixture: Mixture, temperature: float, pressure: float) -> list[float]:
    """Wilson's correlation for the initial K-values of a flash."""
    return [
        (component.Pc.to_base_units().magnitude / pressure)
        * math.exp(
            WILSON_CONSTANT
            * (1.0 + component.omega)
            * (1.0 - component.Tc.to_base_units().magnitude / temperature)
        )
        for component in mixture.components
    ]


class UnifacTables(NamedTuple):
    """The UNIFAC-UMR-PRU tables the UMR mixing rule reads.

    The three interaction matrices are carried **unevaluated**, as the table states them:
    the rule needs both ``a_mn(T)`` and its temperature derivative, and evaluating the
    first here would lose the ``b`` and ``c`` the second is built from.
    """

    #: Per-component group counts, ``N x G`` row-major.
    groups: tuple[float, ...]
    #: The volume ``R`` of each group, length ``G``.
    group_r: tuple[float, ...]
    #: The surface area ``Q`` of each group, length ``G``.
    group_q: tuple[float, ...]
    #: The constant term of the interaction, ``G x G`` row-major, in kelvin.
    aij: tuple[float, ...]
    #: The linear term, ``G x G`` row-major, in kelvin per kelvin.
    bij: tuple[float, ...]
    #: The quadratic term, ``G x G`` row-major, in kelvin per kelvin squared.
    cij: tuple[float, ...]


class ReducedParameters(NamedTuple):
    """The per-component quantities at one state.

    ``A``, ``B`` and ``psi`` are all functions of the temperature and pressure alone,
    not of the composition, which is why they are computed once per solve rather than
    once per iteration: they are the same numbers for the liquid and the vapour, and
    only the composition re-weights them.

    ``psi`` is the logarithmic derivative of the alpha function. It is carried here
    because the mixture's departure functions are the pure form with ``psi`` replaced
    by a composition-weighted average of these.
    """

    #: ``psi`` is the logarithmic derivative of the alpha function. It is carried here
    #: because the mixture's departure functions are the pure form with ``psi`` replaced
    #: by a composition-weighted average of these.
    a: list[float]
    b: list[float]
    psi: list[float]
    #: ``T * dpsi_i/dT`` for every component: the same derivative already multiplied by
    #: the absolute temperature, which is the form the departure heat capacity needs.
    #: Carried because ``kappa`` and ``Tr`` - the two things it is built from - are
    #: local to :func:`reduced_parameters`.
    psi_t: list[float]
    warnings: list[Warning]
    #: The cubic these were reduced under, so the geometry functions downstream read
    #: the right delta pair without a second dispatch. Defaults to Peng-Robinson, which
    #: is what a caller who built a `ReducedParameters` by hand gets.
    cubic: Cubic = PR
    #: The absolute temperature this was reduced at. Carried because the association
    #: contribution is evaluated at a *volume*, and the volume is `z R T / P` - so a
    #: function handed only the reduced parameters cannot reach it. Rust's
    #: `ReducedParameters` carries both for the same reason.
    #:
    #: Defaulted to zero so a fixture that states none - a test working in reduced
    #: variables has no absolute state to give - need not invent one. It is not a silent
    #: default: a mixture that *associates* reaches `Association.solve` with a non-positive
    #: volume or temperature and is refused there, so an unstated state is a missing answer
    #: rather than a plausible wrong one.
    t_kelvin: float = 0.0
    #: The absolute pressure, in Pa. See `t_kelvin`.
    pressure: float = 0.0
    #: The kernel parameters of the association this mixture runs, or `None` when it runs
    #: none. Built here rather than looked up per call because it is a property of the
    #: *mixture* and of neither the temperature nor the pressure.
    association: Association | None = None
    #: The UNIFAC tables the UMR mixing rule reads, or `None` for every mixture that uses
    #: an interaction matrix. Present is what makes :func:`mixture_parameters` mix by the
    #: universal rule instead.
    umr: UnifacTables | None = None
    #: `Tr_i = T / Tc_i`, one per component. The reduced temperature the alpha terms were
    #: evaluated at, carried because the Soreide-Whitson mixing rule reads it again when it
    #: resolves its phase-dependent interaction matrix - and a caller holding the reduced
    #: parameters should not have to recover `Tc_i` from them to do so.
    reduced_temperatures: tuple[float, ...] = ()
    #: Whether the mixture's interaction matrix is a function of the *phase composition*.
    #: True only for the Soreide-Whitson rule, and carried so that
    #: :func:`phase_derivatives` can refuse it: `A_ij` then has a composition derivative
    #: the classical surface does not have.
    composition_dependent_kij: bool = False
    #: The Huron-Vidal rule's fitted matrices, or `None` for every mixture that does not
    #: name that rule. Present is what makes :func:`mixture_parameters` mix by the
    #: co-volume-weighted NRTL.
    huron_vidal: Any = None
    #: The Fürst electrolyte term, or `None`. Present is what makes :func:`phase_state`
    #: solve a root that is not the cubic's and :func:`phase_state_at` add its `ln phi`.
    furst: Any = None


class PhaseState(NamedTuple):
    """One phase's state at a composition."""

    #: The mixture's attraction parameter at this composition.
    a_mix: float
    #: The mixture's repulsion parameter at this composition.
    b_mix: float
    #: The root of the cubic this phase sits on.
    z: float
    #: ``ln phi_i`` for every component.
    ln_phi: list[float]
    #: The departure enthalpy over ``R*T``, for the mixture.
    h_dep_rt: float
    #: The departure entropy over ``R``, for the mixture.
    s_dep_r: float
    #: The composition-weighted average of the components' ``psi``.
    psi_bar: float
    #: The departure heat capacity over ``R``, for the mixture.
    cp_dep_r: float


class PhaseDerivatives(NamedTuple):
    """One phase's ``ln phi`` together with its first derivatives at a state.

    See :func:`phase_derivatives`, which computes it, for what each family is
    differentiated at and why the composition derivative is a mole-number one.
    """

    #: ``ln phi_i``, the same values :class:`PhaseState` carries at this state.
    ln_phi: list[float]
    #: ``d ln phi_i / d n_j`` at constant ``T``, ``P`` and the other mole numbers.
    d_ln_phi_dn: list[list[float]]
    #: ``d ln phi_i / dT`` at constant ``P`` and composition.
    d_ln_phi_dt: list[float]
    #: ``d ln phi_i / dP`` at constant ``T`` and composition.
    d_ln_phi_dp: list[float]


def _soave_kappa(alpha: str, omega: float) -> tuple[float, tuple[Warning, ...]]:
    """The Soave ``m`` for one component, from the alpha correlation named."""
    if alpha == "srk":
        srk_result = srk_kappa(omega)
        return srk_result.kappa, srk_result.warnings
    if alpha == "pr78":
        pr78_result = pr78_kappa(omega)
        return pr78_result.kappa, pr78_result.warnings
    if alpha == "twu":
        twu_result = twu_kappa(omega)
        return twu_result.kappa, twu_result.warnings
    pr_result = pr_kappa(omega)
    return pr_result.kappa, pr_result.warnings


def _non_soave(
    cubic: Cubic, term: Any, reduced_temperature: float, reduced_pressure: float
) -> tuple[float, float, float, float]:
    """The reduced parameters of a non-Soave alpha term, from its `alpha`, `psi` and
    `psi_t`."""
    alpha = term.alpha(reduced_temperature)
    a_reduced = cubic.omega_a * alpha * reduced_pressure / reduced_temperature**2
    b_reduced = cubic.omega_b * reduced_pressure / reduced_temperature
    return a_reduced, b_reduced, term.psi(reduced_temperature), term.psi_t(reduced_temperature)


def reduced_parameters(mixture: Mixture, temperature: float, pressure: float) -> ReducedParameters:
    """``(A_i, B_i, psi_i, T*dpsi_i/dT, warnings)`` for every component at a state."""
    a: list[float] = []
    b: list[float] = []
    psi: list[float] = []
    psi_t: list[float] = []
    reduced_temperatures: list[float] = []
    warnings: list[Warning] = []
    for index, component in enumerate(mixture.components):
        reduced_temperature = temperature / component.Tc.to_base_units().magnitude
        reduced_temperatures.append(reduced_temperature)
        # **A Fürst ion's attraction and covolume are not the cubic's.** NeqSim's
        # `ComponentModifiedFurstElectrolyteEos` overwrites both in its constructor -
        # `b = (p0 d^3 + p1)` from the fitted parameters and `a = 1e-35` - so a port that
        # left them at the cubic's would give every ion a real attraction and a covolume of
        # the wrong size, and the root would be plausible and wrong.
        if mixture.furst is not None:
            covolume = mixture.furst.ion_covolume(index)
            if covolume is not None:
                r_t = R_NEQSIM * temperature
                a.append(1.0e-35 * pressure / (r_t * r_t))
                b.append(covolume * pressure / r_t)
                # The ion's alpha contributes nothing: `psi` is read only through weights
                # carrying `sqrt(A_i A_j)`, and `A` is `1e-35`.
                psi.append(0.0)
                psi_t.append(0.0)
                continue
        reduced_pressure = pressure / component.Pc.to_base_units().magnitude
        # An associating component carries its own attraction and covolume, and NeqSim's
        # `ComponentSrkCPA` substitutes them for the cubic's: `if (|aCPA| > 1e-6) { a =
        # aCPA; b = bCPA; }`. Not a refinement of the critical-constant values - water's
        # fitted covolume is 1.4515e-5 m3/mol against `0.08664 R Tc/Pc`'s 2.11e-5 - and its
        # alpha coefficient is fitted too, so the whole `a_i(T)` differs.
        #
        # **This comes before the alpha dispatch, not inside it**, because the substitution
        # replaces what that dispatch computes: an associating component reads the fitted
        # Soave `m` whatever correlation the mixture named.
        # **UMR-CPA substitutes a third set, and scales it with the model's own alpha.**
        # `ComponentUMRCPA.setAttractiveTerm` installs the five-parameter Mathias-Copeman
        # term over the substituted `a` and `b`, where `ComponentSrkCPA` calls `setm(mCPA)`
        # and gets a Soave coefficient - so which alpha scales the fitted set is the
        # model's choice, and the model states it by asking for `alpha="matcop_5prumr"`.
        umr_cpa_set = (
            component.association.umr_cpa
            if component.association is not None and mixture.alpha == "matcop_5prumr"
            else None
        )
        fitted = _fitted_set(mixture, component)
        if fitted is not None or umr_cpa_set is not None:
            r_t = R * temperature
            if umr_cpa_set is not None:
                # The model's alpha is `matcop_5prumr` here and the component carries its
                # `UMRCPA_MC1..5` - the same term and fallback the dispatch below builds,
                # reached before it because the substitution replaces what it computes.
                umr_cpa_kappa = pr_kappa(component.omega)
                warnings.extend(umr_cpa_kappa.warnings)
                umr_term = MatCop(
                    kappa=umr_cpa_kappa.kappa,
                    params=component.alpha_params,
                    fallback=MatCopFallback.ALL_UNSET,
                )
                alpha_value = umr_term.alpha(reduced_temperature)
                a.append(umr_cpa_set.attraction * alpha_value * pressure / (r_t * r_t))
                b.append(umr_cpa_set.covolume * pressure / r_t)
                psi.append(umr_term.psi(reduced_temperature))
                psi_t.append(umr_term.psi_t(reduced_temperature))
                continue
            family, record = fitted  # type: ignore[misc]
            term = Soave(kappa=record.alpha_m(family))
            a.append(
                record.attraction(family) * term.alpha(reduced_temperature) * pressure / (r_t * r_t)
            )
            b.append(record.covolume(family) * pressure / r_t)
            psi.append(term.psi(reduced_temperature))
            psi_t.append(term.psi_t(reduced_temperature))
            continue
        # The kappa correlation and the reduced parameters belong to the cubic this
        # mixture is: SRK and PR share the Soave alpha form but differ in both the
        # coefficient and the Omega constants; RK is kappa-free.
        if mixture.cubic.name == "rk":
            rk_ab = rk_alpha_ab(reduced_temperature, reduced_pressure)
            warnings.extend(rk_ab.warnings)
            rk_term = RkAlpha()
            a_reduced = rk_ab.a_reduced
            b_reduced = rk_ab.b_reduced
            psi_value = rk_term.psi(reduced_temperature)
            psi_t_value = rk_term.psi_t(reduced_temperature)
        elif mixture.alpha == "soreide_whitson":
            # **Water by role, and the salinity from the rule.** This is the one alpha
            # that is not a function of the component and the temperature alone, and the
            # role vector is the only per-component identity this layer has - so beside
            # another rule there would be nothing to read either from, and this refuses
            # rather than resolving against a zero. Rust's `Alpha::SoreideWhitson` makes
            # the same refusal.
            parameters = mixture.soreide_whitson
            if parameters is None:
                raise InvalidInputError(
                    "mixing_rule",
                    "the Soreide-Whitson alpha takes its salinity from the "
                    "Soreide-Whitson mixing rule and reads which component is water from "
                    "that rule's roles, so it cannot be resolved beside any other rule",
                )
            if parameters.roles[index] == "water":
                a_reduced, b_reduced, psi_value, psi_t_value = _non_soave(
                    mixture.cubic,
                    SoreideWhitsonWater(
                        salinity=parameters.salinity,
                        critical_temperature=component.Tc.to_base_units().magnitude,
                    ),
                    reduced_temperature,
                    reduced_pressure,
                )
            else:
                # Every other component is Peng-Robinson 1978, which is what
                # `AttractiveTermSoreideWhitson.alpha` delegates to.
                soreide_kappa = pr78_kappa(component.omega)
                warnings.extend(soreide_kappa.warnings)
                a_reduced, b_reduced, psi_value, psi_t_value = _non_soave(
                    mixture.cubic,
                    Soave(kappa=soreide_kappa.kappa),
                    reduced_temperature,
                    reduced_pressure,
                )
        elif mixture.alpha == "twucoon":
            a_reduced, b_reduced, psi_value, psi_t_value = _non_soave(
                mixture.cubic, TwuCoon(omega=component.omega), reduced_temperature, reduced_pressure
            )
        elif mixture.alpha == "gassem2001":
            a_reduced, b_reduced, psi_value, psi_t_value = _non_soave(
                mixture.cubic,
                Gassem2001(omega=component.omega),
                reduced_temperature,
                reduced_pressure,
            )
        elif mixture.alpha == "danesh":
            danesh_kappa = pr78_kappa(component.omega)
            warnings.extend(danesh_kappa.warnings)
            a_reduced, b_reduced, psi_value, psi_t_value = _non_soave(
                mixture.cubic,
                Danesh(kappa=danesh_kappa.kappa),
                reduced_temperature,
                reduced_pressure,
            )
        elif mixture.alpha == "schwartzentruber":
            a_reduced, b_reduced, psi_value, psi_t_value = _non_soave(
                mixture.cubic,
                Schwartzentruber(omega=component.omega, params=component.alpha_params),
                reduced_temperature,
                reduced_pressure,
            )
        elif mixture.alpha == "mollerup":
            a_reduced, b_reduced, psi_value, psi_t_value = _non_soave(
                mixture.cubic,
                Mollerup(params=component.alpha_params),
                reduced_temperature,
                reduced_pressure,
            )
        elif mixture.alpha == "matcop":
            a_reduced, b_reduced, psi_value, psi_t_value = _non_soave(
                mixture.cubic,
                MatCop(
                    kappa=matcop_kappa(component.omega),
                    params=component.alpha_params,
                    fallback=MatCopFallback.NONE,
                ),
                reduced_temperature,
                reduced_pressure,
            )
        elif mixture.alpha == "matcop_pr":
            matcop_kappa_value = pr_kappa(component.omega)
            warnings.extend(matcop_kappa_value.warnings)
            a_reduced, b_reduced, psi_value, psi_t_value = _non_soave(
                mixture.cubic,
                MatCop(
                    kappa=matcop_kappa_value.kappa,
                    params=component.alpha_params,
                    fallback=MatCopFallback.SUPERCRITICAL,
                ),
                reduced_temperature,
                reduced_pressure,
            )
        elif mixture.alpha == "matcop_prumr":
            a_reduced, b_reduced, psi_value, psi_t_value = _non_soave(
                mixture.cubic,
                MatCop(
                    kappa=umr_kappa(component.omega),
                    params=component.alpha_params,
                    fallback=MatCopFallback.UNSET,
                ),
                reduced_temperature,
                reduced_pressure,
            )
        elif mixture.alpha == "matcop_5prumr":
            matcop5_kappa_value = pr_kappa(component.omega)
            warnings.extend(matcop5_kappa_value.warnings)
            a_reduced, b_reduced, psi_value, psi_t_value = _non_soave(
                mixture.cubic,
                MatCop(
                    kappa=matcop5_kappa_value.kappa,
                    params=component.alpha_params,
                    fallback=MatCopFallback.ALL_UNSET,
                ),
                reduced_temperature,
                reduced_pressure,
            )
        elif mixture.alpha == "delft1998":
            delft_kappa = pr78_kappa(component.omega)
            warnings.extend(delft_kappa.warnings)
            a_reduced, b_reduced, psi_value, psi_t_value = _non_soave(
                mixture.cubic,
                Delft1998(kappa=delft_kappa.kappa, is_methane=component.alpha_params == (1.0,)),
                reduced_temperature,
                reduced_pressure,
            )
        else:
            # The kappa correlation belongs to the alpha term; the Omega to the cubic.
            kappa_value, kappa_warnings = _soave_kappa(mixture.alpha, component.omega)
            warnings.extend(kappa_warnings)
            if mixture.cubic.name == "pr":
                pr_ab = pr_alpha_ab(kappa_value, reduced_temperature, reduced_pressure)
                warnings.extend(pr_ab.warnings)
                a_reduced = pr_ab.a_reduced
                b_reduced = pr_ab.b_reduced
            else:
                srk_ab = srk_alpha_ab(kappa_value, reduced_temperature, reduced_pressure)
                warnings.extend(srk_ab.warnings)
                a_reduced = srk_ab.a_reduced
                b_reduced = srk_ab.b_reduced
            term = Soave(kappa=kappa_value)
            psi_value = term.psi(reduced_temperature)
            psi_t_value = term.psi_t(reduced_temperature)
        # The alpha term's own two derivatives, so the form lives in one place rather
        # than being restated per call site.
        a.append(a_reduced)
        b.append(b_reduced)
        psi.append(psi_value)
        psi_t.append(psi_t_value)
    return ReducedParameters(
        a=a,
        b=b,
        psi=psi,
        psi_t=psi_t,
        t_kelvin=temperature,
        pressure=pressure,
        association=_association_of(mixture),
        umr=_unifac_tables(mixture),
        cubic=mixture.cubic,
        warnings=warnings,
        reduced_temperatures=tuple(reduced_temperatures),
        composition_dependent_kij=mixture.mixing_rule == "soreide_whitson",
        huron_vidal=mixture.huron_vidal,
        furst=mixture.furst,
    )


def _unifac_tables(mixture: Mixture) -> UnifacTables | None:
    """The UMR rule's UNIFAC tables, or ``None`` for every other rule.

    A property of the fluid and not of the state, so this is a move rather than a
    computation: the interaction is evaluated at a temperature inside :func:`umr_ader`,
    once per state, where the state's temperature is.
    """
    tables = mixture.umr
    if tables is None:
        return None
    return UnifacTables(
        groups=tables.groups,
        group_r=tables.group_r,
        group_q=tables.group_q,
        aij=tables.aij,
        bij=tables.bij,
        cij=tables.cij,
    )


def _fitted_set(mixture: Mixture, component: Component) -> tuple[str, AssociationParameters] | None:
    """The fitted cubic set one component reads, or ``None`` when it reads none.

    Three conditions, and each is a separate refusal rather than one test: the *mixture*
    must be running association (the same substances are a classical fluid under another
    phase model), the *component* must carry a record, and the record must carry a usable
    set at this family - NeqSim's ``|aCPA| > 1e-6`` guard, without which CO2's all-zero
    ``2A`` row gives a covolume of zero and a ``NaN`` reduced pressure.
    """
    if not mixture.associating:
        return None
    record = component.association
    if record is None:
        return None
    family = family_of(mixture.cubic.name)
    if not record.has_fitted_set(family):
        return None
    return family, record


def _association_of(mixture: Mixture) -> Association | None:
    """The association kernel this mixture runs, or ``None`` when it runs none.

    Built from the components rather than stored on them, because the fitted family is the
    *cubic's* and `Cubic` can change after the mixture is built. An associating mixture
    whose every component lacks a usable set at this family is a classical one, which is
    what NeqSim's `|aCPA| > 1e-6` guard amounts to.
    """
    if not mixture.associating:
        return None
    records = [component.association for component in mixture.components]
    if not any(
        record is not None and record.has_fitted_set(family_of(mixture.cubic.name))
        for record in records
    ):
        return None
    family = family_of(mixture.cubic.name)
    # **The UMR-CPA set replaces the family's for a component that carries one**, the same
    # substitution `reduced_parameters` makes for `a` and `b`, and gated the same way: the
    # model states that it is UMR-CPA by asking for `alpha="matcop_5prumr"`. Water is the
    # check - this set's `kappa_AB` is 0.125 and its `eps` 14177 J/mol against the PR
    # family's, so reading the wrong one moves its `ln phi` in the third decimal while
    # leaving the methane beside it almost untouched.
    umr_cpa = mixture.alpha == "matcop_5prumr"
    components: list[AssociationComponent] = []
    for record in records:
        if record is None:
            components.append(NON_ASSOCIATING)
            continue
        scheme = SiteScheme.from_databank_name(record.scheme)
        if umr_cpa and record.umr_cpa is not None and scheme is not None:
            components.append(
                AssociationComponent(
                    scheme,
                    record.umr_cpa.energy_j_per_mol,
                    record.umr_cpa.volume,
                )
            )
            continue
        if not record.has_fitted_set(family) or scheme is None:
            components.append(NON_ASSOCIATING)
            continue
        components.append(
            AssociationComponent(
                scheme,
                record.energy,
                record.volume_srk if family == "srk" else record.volume_pr,
            )
        )
    return Association(components)


#: NeqSim's ``hwfc`` for the ``UNIFAC_UMRPRU`` GE model, ``EosMixingRuleHandler.init``.
#: A property of the pairing rather than of the cubic: the Huron-Vidal rule beside it takes
#: the cubic's own constant, and nothing about a single component distinguishes the two.
UMR_HWFC = -1.0 / 0.53


def _hv_ader(reduced: ReducedParameters, x: list[float]) -> list[float]:
    """The Huron-Vidal rule's ``ader_i = a_i/b_i - ln(gamma_i)/lambda``.

    The excess-Gibbs rule the Fürst model's system installs. Its activity coefficient is a
    co-volume-weighted NRTL over the *reduced* attraction and repulsion, and ``lambda`` is the
    cubic's Huron-Vidal constant - so this is the same shape as the UMR rule's ``ader`` with a
    different activity model behind it and a different constant in front.
    """
    params = reduced.huron_vidal
    lambda_ = reduced.cubic.hv_constant()
    ln_gamma = _hv_ln_gamma(
        x,
        reduced.t_kelvin,
        reduced.a,
        reduced.b,
        list(params.kij),
        list(params.hv_gij),
        list(params.hv_gij_t),
        list(params.hv_alpha),
        list(params.hv_pairs),
        lambda_,
    )
    return [reduced.a[i] / reduced.b[i] - ln_gamma[i] / lambda_ for i in range(len(x))]


def _ge_alpha_mix(reduced: ReducedParameters, x: list[float]) -> float | None:
    """``sum_i x_i ader_i`` for whichever excess-Gibbs rule the mixture runs.

    ``None`` for a classical rule, which is what tells :func:`mixture_parameters` and
    :func:`phase_state_at` they are on the van der Waals one-fluid path. One dispatcher
    rather than a ternary at each of the three callers, because the choice is the *same*
    choice in all three.
    """
    if reduced.umr is not None:
        return umr_ader(reduced, x)[1]
    if reduced.huron_vidal is not None:
        ader = _hv_ader(reduced, x)
        return sum(xi * a for xi, a in zip(x, ader, strict=True))
    return None


def _ge_ader(reduced: ReducedParameters, x: list[float]) -> list[float] | None:
    """The ``ader_i`` themselves, for the ``ln phi`` form the excess-Gibbs rules write."""
    if reduced.umr is not None:
        return umr_ader(reduced, x)[0]
    if reduced.huron_vidal is not None:
        return _hv_ader(reduced, x)
    return None


def _umr_alpha_mix(reduced: ReducedParameters, x: list[float]) -> float | None:
    """``sum_i x_i ader_i``, or ``None`` for every mixture that mixes classically."""
    if reduced.umr is None:
        return None
    return umr_ader(reduced, x)[1]


def umr_ader(
    reduced: ReducedParameters, x: list[float]
) -> tuple[list[float], float, list[float], float]:
    """The UMR rule's ``(ader, alpha_mix, b_der, T d alpha_mix/dT)`` at a composition.

    ``ader[i] = a_i^T/(b_i R T) + hwfc ln gamma_i``, and ``alpha_mix = sum_i x_i ader_i``.
    The mixture's ``A`` is ``B * alpha_mix``, which is why this returns the rule's own
    quantity rather than a single number. ``b_der`` is the linear co-volume sum, the
    ordinary one.

    **The temperature derivative comes back with it because the departure needs it and the
    two must be one expression.** ``T d alpha_mix/dT = sum_i x_i [qPure_i (psi_i - 1) +
    hwfc T d ln gamma_i/dT]`` - the first term from ``a_i^T/(b_i R T)`` carrying ``1/T``,
    the second from the residual only, since the Flory-Huggins combinatorial is built from
    ``r`` and ``q`` and does not move.
    """
    basis = reduced.umr
    if basis is None:
        raise InvalidInputError("mixing_rule", "the UMR rule needs a UNIFAC basis")
    # Imported here rather than at module scope: `components` reaches this module through
    # `reference._association`, so a top-level import of the activity-coefficient kernel
    # would close a cycle.
    from azoth.eos.components import UnifacUmrpruParameters
    from azoth.eos.reference.unifac_umrpru_activity_coefficients import _ln_gamma

    tables = UnifacUmrpruParameters(
        groups=basis.groups,
        group_r=basis.group_r,
        group_q=basis.group_q,
        aij=basis.aij,
        bij=basis.bij,
        cij=basis.cij,
    )
    ln_gamma, t_d_ln_gamma = _ln_gamma(tables, reduced.t_kelvin, x)
    qpure = [reduced.a[i] / reduced.b[i] for i in range(len(x))]
    ader = [qpure[i] + UMR_HWFC * ln_gamma[i] for i in range(len(x))]
    alpha_mix = sum(x[i] * ader[i] for i in range(len(x)))
    t_d_alpha = sum(
        x[i] * (qpure[i] * (reduced.psi[i] - 1.0) + UMR_HWFC * t_d_ln_gamma[i])
        for i in range(len(x))
    )
    return ader, alpha_mix, list(reduced.b), t_d_alpha


def mixture_parameters(
    a: list[float],
    b: list[float],
    kij: tuple[tuple[float, ...], ...],
    x: list[float],
    alpha_mix: float | None = None,
) -> tuple[float, float]:
    """The mixture parameters for a composition.

    ``alpha_mix`` is the UMR rule's own mixing: ``A = B * alpha_mix`` with ``B`` the linear
    co-volume sum, and ``kij`` unread. Passing it is what distinguishes the two, because a
    universal rule has no interaction matrix to average.
    """
    n = len(x)
    b_mix = sum(x[i] * b[i] for i in range(n))
    if alpha_mix is not None:
        return b_mix * alpha_mix, b_mix
    a_mix = 0.0
    for i in range(n):
        for j in range(n):
            a_mix += x[i] * x[j] * (1.0 - kij[i][j]) * math.sqrt(a[i] * a[j])
    return a_mix, b_mix


def phase_state(
    reduced: ReducedParameters,
    kij: tuple[tuple[float, ...], ...],
    x: list[float],
    *,
    liquid: bool,
) -> PhaseState:
    """One phase's state at a composition.

    ``z`` is selected by *ordering* - the smallest admissible root for the liquid,
    the largest for the vapour - never by an initial guess, which is the rule
    ``eos.pr_z_factor`` fixes. A caller who already has a root should use
    :func:`phase_state_at` instead, which does not re-derive it.
    """
    a_mix, b_mix = mixture_parameters(reduced.a, reduced.b, kij, x, _ge_alpha_mix(reduced, x))
    roots = (
        srk_z_factor(a_mix, b_mix)
        if reduced.cubic.name in ("srk", "rk")
        else pr_z_factor(a_mix, b_mix)
    )
    seed = roots.z_min if liquid else roots.z_max
    if reduced.furst is not None:
        # **The electrolyte terms carry a pressure, so the root is not the cubic's.** The
        # residual is `cubic(Z) - P_e Z/P` with `P_e = -R T d(A_e/(R T))/dV`, the same shape
        # the association's has, and the same solver finds it.
        return phase_state_at(reduced, kij, x, _furst_root(reduced, x, seed, a_mix, b_mix, liquid))
    if reduced.association is None:
        return phase_state_at(reduced, kij, x, seed)
    return phase_state_at(
        reduced,
        kij,
        x,
        _associating_root(reduced, x, seed, a_mix, b_mix, liquid=liquid),
    )


def _furst_state(reduced: ReducedParameters, x: list[float], molar_volume: float) -> Any:
    """The Fürst state at one volume, rebuilt from scratch.

    The term is a function of the volume, so the root solve rebuilds it at every step -
    which is what NeqSim's ``volInit`` does inside its own iteration and the reason this is
    a solve and not a stored state.
    """
    term = reduced.furst
    temperature = reduced.t_kelvin
    components = [
        _FurstComponentState(
            charge=s.charge,
            diameter_m=s.diameter_m,
            dielectric=_furst_component_dielectric(s.dielectric_coefficients, temperature),
            dielectric_dt=_furst_component_dielectric_dt(s.dielectric_coefficients, temperature),
            dielectric_dtdt=_furst_component_dielectric_dtdt(
                s.dielectric_coefficients, temperature
            ),
            critical_volume=s.critical_volume,
        )
        for s in term.species
    ]
    short = term.short_range(temperature, x)
    return _furst_build_state(
        temperature,
        molar_volume,
        list(x),
        components,
        _FurstRule[term.rule],
        short.w,
        short.w_dt,
        short.w_dtdt,
    )


def _furst_solve(reduced: ReducedParameters, x: list[float], molar_volume: float) -> _FurstSolution:
    """The Fürst term's three outputs at one volume.

    ``helmholtz_rt_dv`` is ``1e5`` times the sum of the two volume derivatives, because
    NeqSim's ``dFdV`` is against the **total** volume in its own scaled units - the factor
    the identity ``Z = 1 - V 1e5 dFdV`` fixes.
    """
    term = reduced.furst
    state = _furst_state(reduced, x, molar_volume)
    if term.mod2004:
        # The five zeroings the 2004 revision applies: the solvent dielectric constant's
        # temperature derivative, the shielding parameter's and `XLR`'s. Measured: the base's
        # `getSolventDiElectricConstantdT` is `-0.359218709298880` and the variant's is
        # `-0.00000000000000`.
        state.solvent_dielectric_dt = 0.0
        state.dielectric_dt = 0.0
        state.dielectric_dtdt = 0.0
        state.dielectric_dtdv = 0.0
        state.shielding_dt = 0.0
        state.xlr_dt = 0.0
    derivatives = [
        _FurstComponentDerivatives(
            charge=s.charge,
            diameter_m=s.diameter_m,
            dielectric=_furst_component_dielectric(s.dielectric_coefficients, reduced.t_kelvin),
            # `calcWi`: the row sum negated and doubled, as the handler returns it.
            w_i=-2.0 * sum(x[j] * term.table.wij(i, j, reduced.t_kelvin) for j in range(len(x))),
        )
        for i, s in enumerate(term.species)
    ]
    contributions = _furst_ln_phi_contributions(state, derivatives, list(x), term.mod2004)
    return _FurstSolution(
        helmholtz_rt=_furst_fsr2(state) + _furst_flr(state) + _furst_fborn(state),
        helmholtz_rt_dv=1.0e5 * (_furst_fsr2_dv(state) + _furst_flr_dv(state)),
        ln_phi=[c.total() for c in contributions],
    )


def _furst_root(
    reduced: ReducedParameters,
    x: list[float],
    seed: float,
    a_mix: float,
    b_mix: float,
    liquid: bool,
) -> float:
    """The root of a Fürst mixture's equation of state, by the shared solver.

    The residual is the cubic's plus the electrolyte's pressure, exactly as the
    association's is - so this reaches for the same walk, Newton and scan rather than a
    second copy of them.
    """
    r_t = R * reduced.t_kelvin
    pressure = reduced.pressure
    cubic = reduced.cubic

    def residual(z: float) -> float:
        solution = _furst_solve(reduced, x, z * r_t / pressure)
        p_term = -r_t * solution.helmholtz_rt_dv
        return (
            z
            - z / (z - b_mix)
            + a_mix * z / ((z + cubic.delta1 * b_mix) * (z + cubic.delta2 * b_mix))
            - p_term * z / pressure
        )

    root: float = _root_with_extra_pressure(
        residual, seed, b_mix, liquid=liquid, what="electrolyte"
    )
    return root


def _associating_root(
    reduced: ReducedParameters,
    x: list[float],
    seed: float,
    a_mix: float,
    b_mix: float,
    *,
    liquid: bool,
) -> float:
    """The root of an associating mixture's equation of state.

    **The association carries a pressure, so the root is not the cubic's.**
    ``PhaseSrkCPA.molarVolume`` solves ``BonV - (B/n) dFdV() - P B/(n R T) = 0``, and
    ``PhaseSrkCPA.dFdV()`` is ``super.dFdV() + dFCPAdV()`` - the cubic *plus* the
    association. On water/methanol at 300 K and 100 bar that moves the root from 0.15437 to
    0.10505, and NeqSim's own ``A`` and ``B`` give a cubic root of 0.15229, so it is not a
    small correction.

    The residual, with ``P_a = -R T d(A_assoc/(R T))/dV`` evaluated at ``V = Z R T/P``::

        f(Z) = Z - Z/(Z - B) + A Z/((Z + d1 B)(Z + d2 B)) - P_a Z/P

    The first three terms are the cubic's own residual, so a non-associating component of
    the sum reduces to the cubic exactly.
    """
    cubic = reduced.cubic
    r_t = R * reduced.t_kelvin
    pressure = reduced.pressure
    covolumes = [b * r_t / pressure for b in reduced.b]
    association = reduced.association
    assert association is not None  # the caller checked

    def residual(z: float) -> float:
        v = z * r_t / pressure
        state = association.solve(covolumes, x, v, reduced.t_kelvin)
        p_assoc = -r_t * state.d_helmholtz_dv
        return (
            z
            - z / (z - b_mix)
            + a_mix * z / ((z + cubic.delta1 * b_mix) * (z + cubic.delta2 * b_mix))
        ) - p_assoc * z / pressure

    return _root_with_extra_pressure(residual, seed, b_mix, liquid=liquid, what="associating")


def _root_with_extra_pressure(
    residual: Any,
    seed: float,
    b_mix: float,
    *,
    liquid: bool,
    what: str,
) -> float:
    """The walk, Newton and scan a cubic-plus-pressure root solve shares.

    Both terms that carry a pressure need the same three passes - a geometric walk up from
    the covolume for the liquid side, Newton from the cubic's own root, and a linear scan
    with a bisection for whatever neither finds - and they differ only in the residual.
    Factored rather than copied because the *solver* is the part that must not drift: a fix
    to the walk in one is a fix the other needs.

    ``seed`` is the cubic's own root for the side wanted; ``what`` names the term in the
    failure message, so a refusal says which pressure carried the equation away.
    """

    def bisect(lo: float, hi: float) -> float:
        """The zero of `residual` in `[lo, hi]`, by bisection.

        The count is fixed rather than converged on a residual: near the covolume the
        residual is a difference of two large terms, and a bisection stopping at an
        absolute tolerance there would stop early.
        """
        sign = _sign(residual(lo))
        for _ in range(200):
            mid = 0.5 * (lo + hi)
            if _sign(residual(mid)) == sign:
                lo = mid
            else:
                hi = mid
        return 0.5 * (lo + hi)

    floor = b_mix * 1.000_001

    # `liquid` wants the **lowest** zero above `B`, and the seed cannot find it. At low
    # pressure the substituted cubic's attraction is small enough that it has one real
    # root, so `z_min` and `z_max` are the same number - the vapour root - and Newton from
    # it converges there and stops. **The association's pressure is what creates the
    # liquid root**, so it has to be *found* rather than seeded: the residual is `-inf` at
    # `B+`, and the zero is a fraction of a per cent above it, where a linear scan steps
    # straight over it and a cubic root never points.
    #
    # The walk is geometric because the root's *ratio* to `B` is what is bounded, not its
    # distance.
    if liquid:
        previous = floor
        previous_value = residual(previous)
        steps = 64
        for step in range(1, steps + 1):
            z_try = floor * (seed / floor) ** (step / steps)
            value = residual(z_try)
            if _sign(previous_value) != _sign(value):
                return bisect(previous, z_try)
            previous, previous_value = z_try, value

    # Newton, with the step capped to a tenth of the distance to the covolume so it cannot
    # cross the `Z > B` boundary where the cubic term is singular.
    z = max(seed, floor)
    for _ in range(60):
        value = residual(z)
        if abs(value) < 1.0e-12:
            return z
        # Named `offset` rather than `step`, which the sign-change scans above and below
        # already use for their loop counters.
        offset = 1.0e-7 * max(abs(z), 1.0e-6)
        slope = (residual(z + offset) - value) / offset
        if not math.isfinite(slope) or slope == 0.0:
            break
        delta = -value / slope
        ceiling = 0.1 * (z - b_mix)
        if abs(delta) > ceiling:
            delta = math.copysign(ceiling, delta)
        following = z + delta
        if not math.isfinite(following) or following <= b_mix:
            break
        z = following

    # The fallback: scan for a sign change and bisect. NeqSim reaches
    # `molarVolumeChangePhase` on the same failure, so arriving here is upstream's
    # behaviour rather than a port artifact.
    ceiling = max(seed, 1.0) if liquid else max(seed * 8.0, 1.0)
    steps = 400
    previous = floor
    previous_value = residual(previous)
    for step in range(1, steps + 1):
        z_try = floor + (ceiling - floor) * step / steps
        value = residual(z_try)
        if previous_value == 0.0:
            return previous
        if _sign(previous_value) != _sign(value):
            return bisect(previous, z_try)
        previous, previous_value = z_try, value

    raise OutOfRangeError(
        "z",
        seed,
        f"no volume root for this {what} mixture above B = {b_mix}: the residual has "
        f"no sign change between the covolume and {ceiling}",
    )


def _sign(value: float) -> int:
    """`-1`, `0` or `1`, the comparison the sign-change scan is written in."""
    return (value > 0.0) - (value < 0.0)


def phase_state_at(
    reduced: ReducedParameters,
    kij: tuple[tuple[float, ...], ...],
    x: list[float],
    z: float,
) -> PhaseState:
    """One phase's state, at a root the caller has already chosen.

    For a caller holding a compressibility factor - from ``eos.pr_z_factor``, or from
    a flash that solved for one - who does not want it re-derived. The choice of root
    is a *phase*, and a caller holding a ``Z`` has already made it; re-deriving here
    would silently overrule them.
    """
    a, b = reduced.a, reduced.b
    c = reduced.cubic
    n = len(x)
    alpha_mix = _ge_alpha_mix(reduced, x)
    a_mix, b_mix = mixture_parameters(a, b, kij, x, alpha_mix)
    if not z > b_mix:
        raise OutOfRangeError(
            "z",
            z,
            f"the root must exceed the mixture's B = {b_mix}, because `ln(z - B)` is "
            f"otherwise the logarithm of a negative number. `z <= B` is the "
            f"zero-volume limit, which is not a state",
        )

    cross = [
        sum(x[j] * (1.0 - kij[i][j]) * math.sqrt(a[i] * a[j]) for j in range(n)) for i in range(n)
    ]
    i_term = c.i_term(z, b_mix)
    coefficient = c.coefficient(a_mix, b_mix)
    ln_z_minus_b = math.log(z - b_mix)
    # `z - 1` as the equation of state writes it, which is `z - 1` at the cubic's own root
    # and not at an associating mixture's. See `Cubic.eos_z_minus_one`.
    z_minus_one = c.eos_z_minus_one(z, a_mix, b_mix)

    ln_phi = []
    if alpha_mix is not None:
        # **An excess-Gibbs rule writes `ln phi` in a different form.** The classical
        # `factor * i_term` above is the van der Waals one-fluid rule's; a universal rule
        # carries the activity coefficient per component and a volume-derivative term, and
        # the two expressions do not reduce to one another off the pure component - which
        # is the only state at which this library's tests could tell them apart.
        ader = _ge_ader(reduced, x)
        assert ader is not None  # alpha_mix being set is what put us here
        b_der = list(reduced.b)
        delta1, delta2 = c.delta1, c.delta2
        for i in range(n):
            bd = b_der[i]
            fv = alpha_mix * bd * z / ((z + delta1 * b_mix) * (z + delta2 * b_mix))
            ln_phi.append(-ln_z_minus_b + bd / (z - b_mix) - (ader[i] / c.delta_diff) * i_term - fv)
    else:
        for i in range(n):
            b_ratio = b[i] / b_mix
            # The cross-sum factor, which is 1 for a pure component and makes this
            # identical to `eos.pr_departure` at N = 1.
            factor = 2.0 * cross[i] / a_mix - b_ratio
            ln_phi.append(b_ratio * z_minus_one - ln_z_minus_b - coefficient * factor * i_term)

    # The Wertheim association, when the mixture runs it. Its `ln phi` adds to the cubic's
    # and its Helmholtz energy to the departure enthalpy.
    #
    # **The two must move together.** `s_dep_r` below is the Gibbs identity
    # `h_dep_rt - sum_i x_i ln phi_i`, not a second formula - so adding the association to
    # `ln phi` alone would fold it silently into the entropy. The enthalpy's association
    # term is `A/(RT) - T d(A/RT)/dT`, which is why the kernel carries that derivative.
    #
    # `cp_dep_r` deliberately gets **no** association term, matching the Rust kernel: the
    # association's second temperature derivative is not derived, so an associating
    # mixture's `h` and `cp` are not consistent and a model that reads `cp` must refuse
    # rather than report the cubic's.
    assoc_dep_rt = 0.0
    if reduced.association is not None:
        r_t = R * reduced.t_kelvin
        covolumes = [b_i * r_t / reduced.pressure for b_i in b]
        v = z * r_t / reduced.pressure
        kernel = reduced.association.solve(covolumes, x, v, reduced.t_kelvin)
        for i, addition in enumerate(kernel.ln_phi):
            ln_phi[i] += addition
        # **`-T d(A/RT)/dT`, with no `A/(RT)` beside it.** NeqSim's own departure is
        # `Hres = A_res + T SresTV + P V - n R T` with `SresTV = -(dA_res/dT)_V`, so the
        # association contributes `A - T dA/dT` only if the two appearances of `A` cancel -
        # and they do, leaving `-T d(A/RT)/dT`. Carrying the `A/(RT)` as well overshot
        # NeqSim by `FCPA * R T`, which on methane/water at 70 bar is 7.1 J/mol.
        assoc_dep_rt = -reduced.t_kelvin * (
            reduced.association.temperature_derivative(covolumes, x, v, reduced.t_kelvin, kernel)
        )

    # The Fürst electrolyte contribution. **`ln phi` and not the Helmholtz energy**: the
    # term's `dFdN` is that identity's own `dFdN`, and this library's cubic `ln phi` is
    # already `dFdN - ln Z` for its part, so the contributions add.
    if reduced.furst is not None:
        solution = _furst_solve(reduced, x, z * R * reduced.t_kelvin / reduced.pressure)
        for index, addition in enumerate(solution.ln_phi):
            ln_phi[index] += addition

    # The mixture's departure functions. `psi_bar` is the composition-weighted average
    # of the components' `psi`, and the two lines below are then `pr_departure`'s
    # expressions with `psi_bar` in place of `psi` - which is what makes them reduce to
    # it exactly at one component.
    weight_total = 0.0
    weighted_psi = 0.0
    weighted_psi_t = 0.0
    for i in range(n):
        for j in range(n):
            weight = x[i] * x[j] * (1.0 - kij[i][j]) * math.sqrt(a[i] * a[j])
            psi_pair = 0.5 * (reduced.psi[i] + reduced.psi[j])
            weight_total += weight
            weighted_psi += weight * psi_pair
            # ``T*d(weight*psi_pair)/dT`` for this pair. The weight carries
            # ``sqrt(A_i A_j)``, whose logarithmic derivative is ``psi_pair - 2``, and
            # ``psi_pair`` carries the two components' own derivatives.
            weighted_psi_t += weight * (
                (psi_pair - 2.0) * psi_pair + 0.5 * (reduced.psi_t[i] + reduced.psi_t[j])
            )
    psi_bar = weighted_psi / weight_total
    t_dpsi_bar = weighted_psi_t / weight_total - psi_bar * (psi_bar - 2.0)
    # **A universal rule's departure is not the weighted one.** `h_dep_rt`'s cubic term is
    # `(T dA/dT + A)/(delta_diff B) * I`, and for the classical rule that is
    # `A (psi_bar - 1)/(delta_diff B)` because `A` carries `(R T)^-2` through every weight.
    # For the UMR rule `A = B alpha_mix` with **both** factors moving, so `T dA/dT + A`
    # collapses to `B * T d alpha_mix/dT` and the term is `T d alpha_mix/dT / delta_diff * I`
    # - no `psi_bar` in it at all. A pure fluid hides the difference: there `alpha_mix` is
    # one component's `qPure`, and the two expressions agree exactly.
    # **The UMR rule's departure and the Huron-Vidal rule's are not the same expression**,
    # and only the first is assembled here: `t_d_alpha_mix` is the universal rule's own, and
    # for Huron-Vidal this falls back to the classical form - which is what the Rust kernel
    # does, and which is why neither model reports `h_res`.
    if reduced.umr is not None:
        _, _, _, t_d_alpha = umr_ader(reduced, x)
        excess = t_d_alpha / c.delta_diff
    else:
        excess = coefficient * (psi_bar - 1.0)
    h_dep_rt = (z - 1.0) + excess * i_term + assoc_dep_rt
    s_dep_r = h_dep_rt - sum(xi * lp for xi, lp in zip(x, ln_phi, strict=True))

    # The heat-capacity departure. ``a_mix`` moves with temperature exactly as ``A``
    # does in ``eos.pr_departure`` - its logarithmic derivative is ``psi_bar - 2``,
    # because the weights that form it are the same terms - so these lines are that
    # calc's with ``psi_bar`` in place of ``psi``.
    t_da = a_mix * (psi_bar - 2.0)
    t_db = -b_mix
    t_dc = coefficient * (psi_bar - 1.0)
    d_f_dz = c.df_dz(z, a_mix, b_mix)
    t_dfdt = c.t_dfdt(z, a_mix, b_mix, t_da, t_db)
    t_dz = -t_dfdt / d_f_dz
    n_plus = z + c.delta1 * b_mix
    n_minus = z + c.delta2 * b_mix
    t_di = (t_dz + c.delta1 * t_db) / n_plus - (t_dz + c.delta2 * t_db) / n_minus
    cp_dep_r = (
        h_dep_rt
        + t_dz
        + t_dc * (psi_bar - 1.0) * i_term
        + coefficient * t_dpsi_bar * i_term
        + coefficient * (psi_bar - 1.0) * t_di
    )

    return PhaseState(
        a_mix=a_mix,
        b_mix=b_mix,
        z=z,
        ln_phi=ln_phi,
        h_dep_rt=h_dep_rt,
        s_dep_r=s_dep_r,
        psi_bar=psi_bar,
        cp_dep_r=cp_dep_r,
    )


@dataclass(frozen=True, slots=True)
class _AssociationDerivatives:
    """An associating mixture's association terms at one state.

    Carried together because :func:`phase_derivatives` needs all of them at once and none is
    meaningful alone: the root sensitivities say how ``Z`` responds to the state, ``alpha``
    converts one to a molar volume, and the kernel's derivatives are what the association
    adds at constant volume.
    """

    #: ``R T/P``.
    alpha: float
    #: ``dZ/dn_j``, at constant ``T`` and ``P``.
    d_z_n: tuple[float, ...]
    #: ``T dZ/dT``, at constant ``P`` and composition.
    t_d_z: float
    #: ``P dZ/dP``, at constant ``T`` and composition.
    p_d_z: float
    #: The kernel's derivatives at the converged state.
    kernel: SiteDerivatives
    #: The association's ``ln phi_i``, which adds to the cubic's.
    ln_phi: tuple[float, ...]


def _association_derivatives(
    reduced: ReducedParameters,
    x: list[float],
    a_mix: float,
    b_mix: float,
    abar: list[float],
    t_d_a: float,
    t_d_b: float,
    z: float,
) -> _AssociationDerivatives | None:
    """An associating mixture's root sensitivities to ``n_j``, ``T`` and ``P``.

    **The root is the associating residual's, not the cubic's**, so its sensitivity is taken
    from that residual: the cubic's own ``dZ/dn_j`` is ``-0.8`` where the associating one is
    ``-1.7`` at the state the Rust tests measure, and the difference is the association's
    pressure responding to the composition.

    The residual is ``Q(Z) - (Z/P) P_a(x, V(Z))`` where ``Q`` is the cubic's z-cubic and
    ``P_a = -R T d(A_assoc/(R T))/dV`` - the polynomial form ``associating_root`` solves, so
    that ``dQ/dZ`` is the cubic's own ``df_dz`` and the two sensitivities are consistent.
    ``Q`` vanishes at the same ``Z`` as the monic cubic's residual because it *is* that
    residual times ``q(Z) = (Z + delta1 B)(Z + delta2 B)``.
    """
    association = reduced.association
    if association is None:
        return None
    n = len(x)
    cubic = reduced.cubic
    t = reduced.t_kelvin
    p = reduced.pressure
    alpha = R * t / p
    covolumes = [b_i * alpha for b_i in reduced.b]
    v = z * alpha
    state = association.solve(covolumes, x, v, t)
    kernel = association.derivatives(covolumes, x, v, t, state)

    b = b_mix
    ds, dp = cubic.delta_sum, cubic.delta_prod
    q = z * z + ds * b * z + dp * b * b
    s_factor = (z - b) * q
    s_prime = q + (z - b) * (2.0 * z + ds * b)
    d_s_db = -q + (z - b) * (ds * z + 2.0 * dp * b)

    p_assoc = -R * t * state.d_helmholtz_dv
    a_vv = kernel.d2_helmholtz_dv2
    a_vt = kernel.d2_helmholtz_dv_dt
    fa = cubic.df_da(z, b)
    fb = cubic.df_db(z, a_mix, b)
    fz = cubic.df_dz(z, a_mix, b)

    # `P_a` depends on `V`, and `V = Z R T/P` moves with every one of the three state
    # variables - which is why each sensitivity below carries an `A_vv` term beside the
    # explicit one.
    g_z = fz - (s_prime / p) * p_assoc + s_factor * alpha * alpha * a_vv

    d_a_dn = [2.0 * (abar[j] - a_mix) for j in range(n)]
    d_b_dn = [reduced.b[j] - b_mix for j in range(n)]
    # `dA/dn_j` and `dB/dn_j` are the *normalised* partials - the ones that hold the total
    # mole number at one, and the ones NeqSim reports - so the association's composition
    # derivative has to be taken the same way. The kernel's is raw, so the chain rule from
    # `x = n/(sum n)` subtracts the composition-weighted sum of all of them, exactly as
    # `2 abar_j` becomes `2 (abar_j - A)`.
    weighted_a_vn = sum(x[i] * kernel.d2_helmholtz_dv_dn[i] for i in range(n))
    d_z_n = tuple(
        -(
            (fa * d_a_dn[j] + fb * d_b_dn[j])
            - d_s_db * d_b_dn[j] * p_assoc / p
            + s_factor * alpha * (kernel.d2_helmholtz_dv_dn[j] - weighted_a_vn)
        )
        / g_z
        for j in range(n)
    )

    # `T dB/dT = -B`, so `T dS/dT = -B dS/dB`, and
    # `T dP_a/dT = P_a - R T alpha Z A_vv - R T^2 A_vt` at fixed `Z`.
    t_d_p_assoc = p_assoc - R * t * alpha * z * a_vv - R * t * t * a_vt
    # `T d(monic)/dT` is the chain rule over `A` and `B`, which is what `t_d_a`/`t_d_b`
    # already carry.
    t_d_z = (
        -((fa * t_d_a + fb * t_d_b) + (b * d_s_db / p) * p_assoc - (s_factor / p) * t_d_p_assoc)
        / g_z
    )

    # `P dB/dP = B`, and `P dP_a/dP = Z alpha^2 A_vv` at fixed `Z`.
    p_d_z = (
        -(
            fa * a_mix
            + fb * b
            - d_s_db * b * p_assoc / p
            + s_factor * p_assoc / p
            - s_factor * z * alpha * alpha * a_vv
        )
        / g_z
    )

    return _AssociationDerivatives(
        alpha=alpha,
        d_z_n=d_z_n,
        t_d_z=t_d_z,
        p_d_z=p_d_z,
        kernel=kernel,
        ln_phi=state.ln_phi,
    )


def phase_derivatives(
    reduced: ReducedParameters,
    kij: tuple[tuple[float, ...], ...],
    x: list[float],
    z: float,
    *,
    temperature: float,
    pressure: float,
) -> PhaseDerivatives:
    """One phase's ``ln phi`` together with its first derivatives at a state.

    The three families a second-order flash is written in, all analytic:

    * ``d_ln_phi_dn[i][j]`` - ``d ln phi_i / d n_j`` at constant ``T``, ``P`` and the
      other mole numbers.
    * ``d_ln_phi_dt[i]`` - ``d ln phi_i / dT`` at constant ``P`` and composition.
    * ``d_ln_phi_dp[i]`` - ``d ln phi_i / dP`` at constant ``T`` and composition.

    The composition derivative is taken holding the *other mole numbers* fixed, so the
    total moves with ``j`` and the result is not intensive. The values returned are for
    a total mole number of one, which is what a composition summing to one is, and they
    carry the Gibbs-Duhem sum rule ``sum_i x_i d ln phi_i / d n_j = 0`` at that scale -
    the property a constant-composition derivative would not have.

    ``kij`` is the same matrix :func:`phase_state_at` takes, and it is the *constant*
    one: the rules whose interaction parameter moves with temperature are on the Rust
    side alone, so this takes no ``T * dkij/dT`` and none is silently assumed.
    """
    a, b = reduced.a, reduced.b
    c = reduced.cubic
    n = len(x)
    if reduced.composition_dependent_kij:
        # The Soreide-Whitson rule's matrix is a function of the composition this is
        # differentiated against, so `A_ij` below would carry a term nothing here supplies.
        # The same refusal Rust's `Mixture::phase_derivatives` makes.
        raise InvalidInputError(
            "mixing_rule",
            "the Soreide-Whitson rule resolves its interaction matrix against the phase "
            "composition, so its `A_ij` carries a composition derivative this classical "
            "one does not have. The derivative surface covers the rules whose matrix is "
            "fixed by the state",
        )
    if reduced.umr is not None:
        # The activity-coefficient rules write `ln phi` as `ader`, `alpha_mix` and
        # `b_der`, whose composition derivative is the excess Gibbs energy's second
        # derivative rather than this classical one. Asking for it is an error naming
        # that, not a silent return of the wrong matrix - the same refusal the Rust
        # `Mixture::phase_derivatives` makes.
        raise InvalidInputError(
            "mixing_rule",
            "the UMR rule writes ln phi as `ader`, `alpha_mix` and `b_der`, whose "
            "composition derivative is the excess Gibbs energy's second derivative "
            "rather than this classical one. The derivative surface covers the cubic "
            "family, and this rule is outside it",
        )
    a_mix, b_mix = mixture_parameters(a, b, kij, x)
    if not z > b_mix:
        raise OutOfRangeError(
            "z",
            z,
            f"the root must exceed the mixture's B = {b_mix}, because `ln(z - B)` is "
            f"otherwise the logarithm of a negative number",
        )

    delta1, delta2 = c.delta1, c.delta2
    delta_diff = c.delta_diff
    i_term = c.i_term(z, b_mix)
    coefficient = c.coefficient(a_mix, b_mix)
    ln_z_minus_b = math.log(z - b_mix)
    fz = c.df_dz(z, a_mix, b_mix)
    fa = c.df_da(z, b_mix)
    fb = c.df_db(z, a_mix, b_mix)
    z_minus_one, zmo_dz, zmo_da, zmo_db = c.eos_z_minus_one_partials(z, a_mix, b_mix)

    a_ij = [[(1.0 - kij[i][j]) * math.sqrt(a[i] * a[j]) for j in range(n)] for i in range(n)]
    abar = [sum(x[j] * a_ij[i][j] for j in range(n)) for i in range(n)]
    b_ratio = [b[i] / b_mix for i in range(n)]
    factor = [2.0 * abar[i] / a_mix - b_ratio[i] for i in range(n)]
    ln_phi = [
        b_ratio[i] * z_minus_one - ln_z_minus_b - coefficient * factor[i] * i_term for i in range(n)
    ]

    # --- composition, at constant temperature and pressure -------------------
    #
    # `A_i` and `B_i` are functions of `(T, P)` alone, so only the composition moves. At
    # unit total mole number `dA/dn_j = 2 (abar_j - A)` and `dB/dn_j = B_j - B`.
    d_a_dn = [2.0 * (abar[j] - a_mix) for j in range(n)]
    d_b_dn = [b[j] - b_mix for j in range(n)]

    # `A_ij` moves with temperature through the components' alphas, and its logarithmic
    # derivative is `psi_pair - 2`. `B_ij` scales as `1/T`, so `T dB/dT = -B`. Computed
    # here rather than below because the association's root sensitivities need both.
    t_d_a_ij = [[0.0] * n for _ in range(n)]
    t_d_abar = [0.0] * n
    for i in range(n):
        for j in range(n):
            psi_pair = 0.5 * (reduced.psi[i] + reduced.psi[j])
            t_d_a_ij[i][j] = a_ij[i][j] * (psi_pair - 2.0)
            t_d_abar[i] += x[j] * t_d_a_ij[i][j]
    t_d_a = sum(x[i] * t_d_abar[i] for i in range(n))
    t_d_b = -b_mix

    # --- the root's sensitivities --------------------------------------------
    #
    # **An associating mixture's root is not the cubic's, so its sensitivity is not either.**
    # The association carries a pressure that responds to the composition and to the
    # temperature, and the whole of the difference between this surface and the cubic's is
    # here.
    sensitivities = _association_derivatives(reduced, x, a_mix, b_mix, abar, t_d_a, t_d_b, z)
    if sensitivities is None:
        d_z_dn = [-(fa * d_a_dn[j] + fb * d_b_dn[j]) / fz for j in range(n)]
        t_d_z = -(fa * t_d_a + fb * t_d_b) / fz
        p_d_z = -(fa * a_mix + fb * b_mix) / fz
    else:
        d_z_dn = list(sensitivities.d_z_n)
        t_d_z = sensitivities.t_d_z
        p_d_z = sensitivities.p_d_z
        # The association's `ln phi`, which `PhaseState::ln_phi` carries too: the surface is
        # the same values, so a caller comparing the two must find them equal.
        for i in range(n):
            ln_phi[i] += sensitivities.ln_phi[i]

    d_coeff_dn = [
        d_a_dn[j] / (delta_diff * b_mix) - a_mix * d_b_dn[j] / (delta_diff * b_mix * b_mix)
        for j in range(n)
    ]
    d_i_dn = [
        (d_z_dn[j] + delta1 * d_b_dn[j]) / (z + delta1 * b_mix)
        - (d_z_dn[j] + delta2 * d_b_dn[j]) / (z + delta2 * b_mix)
        for j in range(n)
    ]

    d_ln_phi_dn = [[0.0] * n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            # `d(abar_i)/dn_j = A_ij - abar_i`, and the `B_i/B` term's derivative is
            # `-B_i/B**2 dB/dn_j`, which is `factor_i`'s `-b_ratio_i` differentiated.
            d_factor = (
                2.0 * ((a_ij[i][j] - abar[i]) * a_mix - abar[i] * d_a_dn[j]) / (a_mix * a_mix)
                + b_ratio[i] * d_b_dn[j] / b_mix
            )
            # `B_i/B` times the equation of state's `z - 1`, so the composition moves it
            # through `z`, through `A` and through `B`, and the `B` in the denominator
            # moves it again.
            d_z_minus_one = zmo_dz * d_z_dn[j] + zmo_da * d_a_dn[j] + zmo_db * d_b_dn[j]
            d_ln_phi_dn[i][j] = (
                b_ratio[i] * d_z_minus_one
                - (d_z_dn[j] - d_b_dn[j]) / (z - b_mix)
                # `b_ratio_i` is `B_i/B`, and `B` moves with the composition, so the
                # term above is not the whole of it. The `B_i/B` inside `factor_i`
                # carries the same derivative with the opposite sign, and the two do
                # not cancel.
                - b_ratio[i] * z_minus_one * d_b_dn[j] / b_mix
                - d_coeff_dn[j] * factor[i] * i_term
                - coefficient * d_factor * i_term
                - coefficient * factor[i] * d_i_dn[j]
            )

    # The same normalisation as the root's: the kernel's composition derivative holds the
    # other mole numbers, not the total, so the chain rule from `x = n/(sum n)` subtracts
    # the composition-weighted sum of the row.
    if sensitivities is not None:
        for i in range(n):
            kernel_row = sensitivities.kernel.d_ln_phi_dn[i]
            weighted = sum(x[k] * kernel_row[k] for k in range(n))
            volume = sensitivities.alpha * sensitivities.kernel.d_ln_phi_dv[i]
            for j in range(n):
                d_ln_phi_dn[i][j] += kernel_row[j] - weighted + volume * sensitivities.d_z_n[j]

    # --- temperature, at constant pressure and composition --------------------
    t_d_coeff = t_d_a / (delta_diff * b_mix) - a_mix * t_d_b / (delta_diff * b_mix * b_mix)
    t_d_i = (t_d_z + delta1 * t_d_b) / (z + delta1 * b_mix) - (t_d_z + delta2 * t_d_b) / (
        z + delta2 * b_mix
    )
    d_ln_phi_dt = [0.0] * n
    for i in range(n):
        # `b_ratio_i` is `B_i/B`, two quantities that both scale as `1/T`, so it is
        # temperature-independent and drops out of the sum below.
        t_d_factor = 2.0 * (t_d_abar[i] * a_mix - abar[i] * t_d_a) / (a_mix * a_mix)
        t_d_z_minus_one = zmo_dz * t_d_z + zmo_da * t_d_a + zmo_db * t_d_b
        t_d_ln_phi = (
            b_ratio[i] * t_d_z_minus_one
            - (t_d_z - t_d_b) / (z - b_mix)
            - t_d_coeff * factor[i] * i_term
            - coefficient * t_d_factor * i_term
            - coefficient * factor[i] * t_d_i
        )
        d_ln_phi_dt[i] = t_d_ln_phi / temperature

    # And the association's, whose volume moves with the temperature as
    # `dV/dT = (R/P)(Z + T dZ/dT)`.
    if sensitivities is not None:
        volume = sensitivities.alpha / reduced.t_kelvin * (z + sensitivities.t_d_z)
        for i in range(n):
            d_ln_phi_dt[i] += (
                sensitivities.kernel.d_ln_phi_dt[i] + volume * sensitivities.kernel.d_ln_phi_dv[i]
            )

    # --- pressure, at constant temperature and composition --------------------
    #
    # Every `A_i` and `B_i` is linear in `P`, so `A`, `B` and their row sums scale as
    # `P` while `coefficient` and `factor_i` - ratios of two such quantities - do not
    # move at all.
    p_d_i = (p_d_z + delta1 * b_mix) / (z + delta1 * b_mix) - (p_d_z + delta2 * b_mix) / (
        z + delta2 * b_mix
    )
    d_ln_phi_dp = [0.0] * n
    for i in range(n):
        p_d_z_minus_one = zmo_dz * p_d_z + zmo_da * a_mix + zmo_db * b_mix
        p_d_ln_phi = (
            b_ratio[i] * p_d_z_minus_one
            - (p_d_z - b_mix) / (z - b_mix)
            - coefficient * factor[i] * p_d_i
        )
        d_ln_phi_dp[i] = p_d_ln_phi / pressure

    # The association carries no explicit pressure, so its whole contribution is the
    # volume's: `dV/dP = (R T/P)(dZ/dP - Z/P)`.
    if sensitivities is not None:
        volume = sensitivities.alpha * (sensitivities.p_d_z - z) / reduced.pressure
        for i in range(n):
            d_ln_phi_dp[i] += volume * sensitivities.kernel.d_ln_phi_dv[i]

    return PhaseDerivatives(
        ln_phi=ln_phi,
        d_ln_phi_dn=d_ln_phi_dn,
        d_ln_phi_dt=d_ln_phi_dt,
        d_ln_phi_dp=d_ln_phi_dp,
    )


def helmholtz_energy(
    reduced: ReducedParameters,
    kij: tuple[tuple[float, ...], ...],
    n: list[float],
    compressibility: float,
) -> float:
    """``A^R/(R T)`` - the residual Helmholtz energy, at a set of mole numbers.

    The energy whose first composition derivative is the logarithm of fugacity and
    whose second is :func:`helmholtz_hessian`. It is exposed because those two
    relations are what the tests are built on, and a function that cannot be
    evaluated cannot have its derivatives checked.

        ``d(A^R/RT)/dn_i = ln phi_i + ln Z``

    which is asserted against :func:`phase_state_at`'s ``ln phi`` - a different code
    path, reached by differentiating a departure function rather than an energy, so
    the agreement is evidence rather than a tautology.

    It is also the quantity a stability analysis minimises, and the reason the
    critical point is written in Helmholtz terms at all.

    # Mole numbers, not mole fractions

    ``n`` is a set of mole numbers, not a composition, and the argument is named for
    it because the distinction is load-bearing: ``A^R`` is homogeneous of degree one
    in ``(V, n)`` but **not** in ``n`` at fixed ``V``, so the total is part of the
    state and a function of the fractions alone could not be differentiated. A
    caller holding a composition summing to one is already passing mole numbers.
    """
    count = len(n)
    # `a_i/(R T V)` is `A_i/Z` and `b_i/V` is `B_i/Z`. Those two identities are what
    # make the whole construction dimensionless; the scaled constants below are
    # therefore independent of `n`, which is what lets the sums carry all of the
    # composition dependence.
    a_hat = [value / compressibility for value in reduced.a]
    b_hat = [value / compressibility for value in reduced.b]
    a_ij = [
        [(1.0 - kij[i][j]) * math.sqrt(a_hat[i] * a_hat[j]) for j in range(count)]
        for i in range(count)
    ]
    total = sum(n)
    b_sum = sum(n[i] * b_hat[i] for i in range(count))
    if not b_sum > 0.0:
        raise OutOfRangeError(
            "compressibility",
            compressibility,
            f"the mixture's `B/Z` came out as {b_sum}, and the logarithmic terms of "
            f"the Helmholtz energy are written against it. It is positive for any "
            f"admissible root, so this is a composition or a root that is not a state "
            f"rather than a compressibility that is out of range",
        )
    a_sum = sum(n[i] * n[j] * a_ij[i][j] for i in range(count) for j in range(count))
    g = reduced.cubic.helmholtz_g(b_sum)
    cubic = -total * math.log(1.0 - b_sum) - (a_sum / (reduced.cubic.delta_diff * b_sum)) * g
    # **The association's own `A^R/RT`, at the same volume.** Without it this surface is a
    # cubic's, and the two things that read it - `eos.stability_test`'s choice of the feed's
    # root, and every `RootSide` choice behind it - compare two roots by the wrong energy.
    # An associating mixture's two roots at 356 K and 1 bar differ by `0.157 RT`, and the
    # association is most of that.
    if reduced.association is None:
        return cubic
    r_t = R * reduced.t_kelvin
    # **Mole numbers, and the volume `V = Z R T/P`.** This function's `n` are mole numbers
    # and its volume is fixed by `Z` alone - `b_i/V` is `B_i/Z` for every component, whatever
    # the moles - so a perturbation at fixed `Z` *is* a perturbation at fixed volume. The
    # kernel's `A_assoc/(RT)` is `sum_i n_i sum_A (ln X_A - X_A/2 + 1/2)`, extensive in the
    # moles it is given; fractions there were wrong by the total, which the finite difference
    # caught, and a volume scaled by the total was wrong again.
    covolumes = [b_i * r_t / reduced.pressure for b_i in reduced.b]
    state = reduced.association.solve(
        covolumes, list(n), compressibility * r_t / reduced.pressure, reduced.t_kelvin
    )
    return cubic + state.helmholtz_rt


def helmholtz_hessian(
    reduced: ReducedParameters,
    kij: tuple[tuple[float, ...], ...],
    n: list[float],
    compressibility: float,
) -> list[list[float]]:
    """``d2(A^R/RT)/dn_i dn_j`` at constant temperature and **volume**.

    Not at constant pressure, and the difference is the whole reason this function
    exists rather than being one more derivative of :func:`phase_state_at`. The
    Hessian of the Helmholtz energy is the quantity every criticality condition is
    written in, because its vanishing is what separates a stable phase from a
    metastable one, and it is only the Helmholtz Hessian that has that meaning.
    A constant-pressure composition derivative answers a different question, and
    converting between the two frames needs two further derivative families and a
    partial-molar-volume correction that this avoids entirely.

    **Everything is dimensionless**, like the rest of the ``eos`` core:
    ``a_i/(R T V)`` is ``A_i/Z`` and ``b_i/V`` is ``B_i/Z``, so the reduced
    parameters the caller already holds carry the whole construction and no
    dimensioned quantity appears. ``compressibility`` is the cubic's root for the
    phase, which is what turns those two identities.

    ``n`` is a set of mole numbers, as in :func:`helmholtz_energy`. The criticality
    conditions want the composition, which is a set of mole numbers summing to one,
    so a caller passing one gets the Hessian those conditions are written against.

    # What it is checked against

    A finite difference of :func:`helmholtz_energy`, and the identity that ties the
    energy to :func:`phase_state_at`: ``d(A^R/RT)/dn_i`` must be ``ln phi_i + ln Z``.
    Both are data-free, and the second reaches the same quantity by a different
    route - differentiating a departure function rather than an energy - so the
    agreement is evidence rather than a restatement.

    # The form, and why it is written out rather than factored

    With ``b = sum_i n_i B_i/Z`` (the mixture's ``B`` over the root, times the total),
    ``L = ln(1 - b)`` and ``G = ln((1 + (1+sqrt2)b)/(1 + (1-sqrt2)b))``, every term
    below is one of those three differentiated once or twice. They are written out
    because the alternative - a factored form - hides which derivative each term
    came from, and this is the piece of the critical point most likely to be got
    wrong.
    """
    count = len(n)
    # **Refused rather than answered from the cubic.** The energy carries the association
    # and this does not, so for an associating mixture the two would be inconsistent - the
    # Hessian would not be the energy's second derivative, and `criticality_matrix` would
    # place a critical point on a cubic the mixture is not. The same rule the heat-capacity
    # departure follows: report the association's absence rather than the cubic's value.
    if reduced.association is not None:
        raise InvalidInputError(
            "mixture",
            "this mixture associates, and the association's second composition derivative "
            "is not assembled - so the Hessian here would be the cubic's and not the "
            "energy's",
        )
    a_hat = [value / compressibility for value in reduced.a]
    b_hat = [value / compressibility for value in reduced.b]
    a_ij = [
        [(1.0 - kij[i][j]) * math.sqrt(a_hat[i] * a_hat[j]) for j in range(count)]
        for i in range(count)
    ]
    total = sum(n)
    b_sum = sum(n[i] * b_hat[i] for i in range(count))
    if not b_sum > 0.0:
        raise OutOfRangeError(
            "compressibility",
            compressibility,
            f"the mixture's `B/Z` came out as {b_sum}, and the logarithmic terms of "
            f"the Helmholtz energy are written against it. It is positive for any "
            f"admissible root, so this is a composition or a root that is not a state "
            f"rather than a compressibility that is out of range",
        )

    d = 1.0 - b_sum
    # L(b) = ln(1 - b), G(b) = ln((1 + delta1 b)/(1 + delta2 b)).
    l_prime = -1.0 / d
    l_second = -1.0 / (d * d)
    g = reduced.cubic.helmholtz_g(b_sum)
    g_prime = reduced.cubic.helmholtz_g_prime(b_sum)
    g_second = reduced.cubic.helmholtz_g_second(b_sum)
    half_delta_diff = reduced.cubic.half_delta_diff
    delta_diff = reduced.cubic.delta_diff

    a_bar = [sum(n[j] * a_ij[i][j] for j in range(count)) for i in range(count)]
    a_sum = sum(n[i] * n[j] * a_ij[i][j] for i in range(count) for j in range(count))

    hessian = [[0.0] * count for _ in range(count)]
    for i in range(count):
        for j in range(count):
            pair = a_bar[i] * b_hat[j] + a_bar[j] * b_hat[i]
            product = b_hat[i] * b_hat[j]
            hessian[i][j] = (
                -(b_hat[i] + b_hat[j]) * l_prime
                - total * product * l_second
                - g * a_ij[i][j] / (half_delta_diff * b_sum)
                + g * pair / (half_delta_diff * b_sum * b_sum)
                - g * a_sum * product / (half_delta_diff * b_sum * b_sum * b_sum)
                - g_prime * pair / (half_delta_diff * b_sum)
                + g_prime * a_sum * product / (half_delta_diff * b_sum * b_sum)
                - g_second * a_sum * product / (delta_diff * b_sum)
            )
    return hessian


def criticality_matrix(
    reduced: ReducedParameters,
    kij: tuple[tuple[float, ...], ...],
    n: list[float],
    compressibility: float,
) -> list[list[float]]:
    """Heidemann & Khalil's ``Q``, whose smallest eigenvalue vanishes at a critical point.

    ``Q_ij = sqrt(n_i n_j) * d2(A/RT)/dn_i dn_j`` at constant temperature and volume,
    the **total** Helmholtz energy rather than the residual one. The ideal part of
    that Hessian at constant volume is ``delta_ij / n_i`` - not the
    ``delta_ij/n_i - 1/n`` that appears in the constant-pressure frame, where the
    ``-1/n`` is the entropy of mixing. Using the pressure form here shifts every
    diagonal entry by ``-1``, and at a pure component's critical point it turns a
    quantity that should be zero into exactly ``-1``.

    **The scaling by ``sqrt(n_i n_j)`` is not decoration.** A Maxwell relation makes
    the Hessian symmetric already; the scaling makes that symmetry *structural*
    rather than numerical, so the eigenvalues are real and the eigenvectors are
    available in every case rather than almost every case.

    **``Q`` is not singular at ordinary states**, which is worth stating because the
    opposite is easy to assume. ``A(T, V, n)`` is not homogeneous in ``n`` at fixed
    ``V`` - homogeneity needs the volume to scale with it - so there is no
    Euler-theorem null vector, and the ideal part of the Hessian at constant volume
    (``delta_ij/n_i``) is positive definite on its own. The vanishing is therefore
    informative rather than generic, and the direction it vanishes along is the
    critical composition fluctuation and nothing else; the spec's notes record the
    sweep that measures it.

    That is why the critical point solves for the **smallest-magnitude eigenvalue**
    rather than for ``det(Q)``. The two are not the same equation: the determinant
    is the product of every eigenvalue, so it vanishes when *any* of them does -
    including ones whose vanishing is not criticality - and it is a product, so it
    is badly scaled for a Newton step. The eigenvalue is the quantity whose
    vanishing is the condition.
    """
    hessian = helmholtz_hessian(reduced, kij, n, compressibility)
    count = len(n)
    return [
        [
            math.sqrt(n[i] * n[j]) * (hessian[i][j] + (1.0 / n[i] if i == j else 0.0))
            for j in range(count)
        ]
        for i in range(count)
    ]


def rachford_rice_bounds(k: list[float]) -> tuple[float, float] | None:
    """The interval on which Rachford-Rice has its physical root, or ``None``.

    ``g(beta) = sum_i z_i (K_i - 1) / (1 + beta (K_i - 1))`` has a pole at
    ``1/(1 - K_i)`` for every component, so the root a flash wants - the one that
    keeps ``x_i`` and ``y_i`` non-negative - lies between the largest pole below it
    and the smallest above it. It exists only when the K-values straddle one, so
    ``None`` is a proof that the feed has no two-phase solution rather than a failure
    to bracket: ``sum_i K_i x_i = 1`` alongside ``sum_i x_i = 1`` needs every ``K_i``
    below one and above one at once.
    """
    lo = -math.inf
    hi = math.inf
    for value in k:
        if value > 1.0:
            lo = max(lo, 1.0 / (1.0 - value))
        elif value < 1.0:
            hi = min(hi, 1.0 / (1.0 - value))
        else:
            # `K_i = 1` exactly puts a pole at infinity and makes `g` degenerate.
            return None
    if not (math.isfinite(lo) and math.isfinite(hi)) or lo >= hi:
        return None
    return lo, hi


def rachford_rice(
    z: list[float],
    k: list[float],
    bounds: tuple[float, float],
    tolerance: float,
    max_iterations: int,
) -> float:
    """The vapour fraction that solves Rachford-Rice, by bisection.

    ``bounds`` must be the interval :func:`rachford_rice_bounds` returned, on which
    ``g`` is continuous, strictly decreasing and positive at the lower end.
    """
    lo, hi = bounds

    def g(beta: float) -> float:
        return sum(zi * (ki - 1.0) / (1.0 + beta * (ki - 1.0)) for zi, ki in zip(z, k, strict=True))

    for _ in range(max_iterations):
        if hi - lo <= tolerance:
            break
        mid = 0.5 * (lo + hi)
        if g(mid) > 0.0:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def rms_delta(ln_k_new: list[float], k: list[float]) -> float:
    """The rms change in ``ln K`` across one iteration."""
    n = len(ln_k_new)
    total = sum((new - math.log(old)) ** 2 for new, old in zip(ln_k_new, k, strict=True))
    return math.sqrt(total / n)


def trivial_warning() -> Warning:
    """The caveat every trivial solution carries."""
    return Warning(
        code=WarningCode.TRIVIAL_SOLUTION,
        message=(
            "the iteration converged to x = y = z, so the feed is single phase and "
            "there is no vapour fraction. Which single phase it is, this model does "
            "not say - that needs a stability analysis it does not perform. `beta` is "
            "absent rather than zero, and `phase` is `trivial`."
        ),
        field=None,
    )


def single_phase_warning(phase: Phase) -> Warning:
    """The caveat a feed carries when no Rachford-Rice root exists at all."""
    which = "vapour" if phase is Phase.ALL_VAPOUR else "liquid"
    return Warning(
        code=WarningCode.TRIVIAL_SOLUTION,
        message=(
            "every K-value is on the same side of one, so the Rachford-Rice equation "
            "has no root and the feed has no two-phase solution at this temperature "
            f"and pressure. The feed is single-phase {which}, and `beta` is absent "
            "rather than zero because there is no vapour fraction to report."
        ),
        field=None,
    )


def negative_flash_warning(phase: Phase, beta: float) -> Warning:
    """The caveat a converged-but-out-of-range vapour fraction carries."""
    which = "superheated vapour" if phase is Phase.ALL_VAPOUR else "subcooled liquid"
    return Warning(
        code=WarningCode.OUT_OF_VALID_RANGE,
        message=(
            f"the vapour fraction is {beta}, outside [0, 1], so the feed is "
            f"single-phase {which}. It is reported because the negative flash is a "
            f"real reading - it is the amount of the absent phase that would have to "
            f"be added to bring the feed to saturation - but it is not a phase split."
        ),
        field=None,
    )


def compositions(z: list[float], k: list[float], beta: float) -> tuple[list[float], list[float]]:
    """The two phases' compositions at a vapour fraction."""
    x = [zi / (1.0 + beta * (ki - 1.0)) for zi, ki in zip(z, k, strict=True)]
    y = [ki * xi for ki, xi in zip(k, x, strict=True)]
    return x, y


def is_trivial(k: list[float], tolerance: float) -> bool:
    """Whether every K-value has collapsed to 1."""
    return all(abs(math.log(value)) < tolerance for value in k)


def normalise(values: list[float]) -> list[float]:
    """Rescale a composition to sum to one."""
    total = sum(values)
    return [value / total for value in values] if total > 0.0 else values
