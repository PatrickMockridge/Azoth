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

from azoth.core.errors import OutOfRangeError
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
    warnings: list[Warning] = []
    for component in mixture.components:
        reduced_temperature = temperature / component.Tc.to_base_units().magnitude
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
        fitted = _fitted_set(mixture, component)
        if fitted is not None:
            family, record = fitted
            term = Soave(kappa=record.alpha_m(family))
            r_t = R * temperature
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
        cubic=mixture.cubic,
        warnings=warnings,
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
    components: list[AssociationComponent] = []
    for record in records:
        if record is None or not record.has_fitted_set(family):
            components.append(NON_ASSOCIATING)
            continue
        scheme = SiteScheme.from_databank_name(record.scheme)
        if scheme is None:
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


def mixture_parameters(
    a: list[float], b: list[float], kij: tuple[tuple[float, ...], ...], x: list[float]
) -> tuple[float, float]:
    """The van der Waals one-fluid mixture parameters for a composition."""
    n = len(x)
    b_mix = sum(x[i] * b[i] for i in range(n))
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
    a_mix, b_mix = mixture_parameters(reduced.a, reduced.b, kij, x)
    roots = (
        srk_z_factor(a_mix, b_mix)
        if reduced.cubic.name in ("srk", "rk")
        else pr_z_factor(a_mix, b_mix)
    )
    seed = roots.z_min if liquid else roots.z_max
    if reduced.association is None:
        return phase_state_at(reduced, kij, x, seed)
    return phase_state_at(
        reduced,
        kij,
        x,
        _associating_root(reduced, x, seed, a_mix, b_mix, liquid=liquid),
    )


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
        f"no volume root for this associating mixture above B = {b_mix}: the residual has "
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
    a_mix, b_mix = mixture_parameters(a, b, kij, x)
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
        assoc_dep_rt = kernel.helmholtz_rt - reduced.t_kelvin * (
            reduced.association.temperature_derivative(covolumes, x, v, reduced.t_kelvin, kernel)
        )

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
    h_dep_rt = (z - 1.0) + coefficient * (psi_bar - 1.0) * i_term + assoc_dep_rt
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
    return -total * math.log(1.0 - b_sum) - (a_sum / (reduced.cubic.delta_diff * b_sum)) * g


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
