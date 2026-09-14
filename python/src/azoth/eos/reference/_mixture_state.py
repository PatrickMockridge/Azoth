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
from typing import Any, NamedTuple

from azoth.core.errors import OutOfRangeError
from azoth.core.warnings import Warning
from azoth.eos.mixture import Mixture
from azoth.eos.reference.alpha_term import Soave
from azoth.eos.reference.cubic import PR
from azoth.eos.reference.pr_alpha_ab import pr_alpha_ab
from azoth.eos.reference.pr_kappa import pr_kappa
from azoth.eos.reference.pr_z_factor import pr_z_factor

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


def reduced_parameters(mixture: Mixture, temperature: float, pressure: float) -> ReducedParameters:
    """``(A_i, B_i, psi_i, T*dpsi_i/dT, warnings)`` for every component at a state."""
    a: list[float] = []
    b: list[float] = []
    psi: list[float] = []
    psi_t: list[float] = []
    warnings: list[Warning] = []
    for component in mixture.components:
        kappa = pr_kappa(component.omega)
        warnings.extend(kappa.warnings)
        reduced_temperature = temperature / component.Tc.to_base_units().magnitude
        ab = pr_alpha_ab(
            kappa.kappa,
            reduced_temperature,
            pressure / component.Pc.to_base_units().magnitude,
        )
        warnings.extend(ab.warnings)
        # The alpha term's own two derivatives, so the Soave form lives in one place
        # rather than being restated per call site.
        term = Soave(kappa=kappa.kappa)
        a.append(ab.a_reduced)
        b.append(ab.b_reduced)
        psi.append(term.psi(reduced_temperature))
        psi_t.append(term.psi_t(reduced_temperature))
    return ReducedParameters(a=a, b=b, psi=psi, psi_t=psi_t, warnings=warnings)


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
    roots = pr_z_factor(a_mix, b_mix)
    return phase_state_at(reduced, kij, x, roots.z_min if liquid else roots.z_max)


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
    i_term = PR.i_term(z, b_mix)
    coefficient = PR.coefficient(a_mix, b_mix)
    ln_z_minus_b = math.log(z - b_mix)

    ln_phi = []
    for i in range(n):
        b_ratio = b[i] / b_mix
        # The cross-sum factor, which is 1 for a pure component and makes this
        # identical to `eos.pr_departure` at N = 1.
        factor = 2.0 * cross[i] / a_mix - b_ratio
        ln_phi.append(b_ratio * (z - 1.0) - ln_z_minus_b - coefficient * factor * i_term)

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
    h_dep_rt = (z - 1.0) + coefficient * (psi_bar - 1.0) * i_term
    s_dep_r = h_dep_rt - sum(xi * lp for xi, lp in zip(x, ln_phi, strict=True))

    # The heat-capacity departure. ``a_mix`` moves with temperature exactly as ``A``
    # does in ``eos.pr_departure`` - its logarithmic derivative is ``psi_bar - 2``,
    # because the weights that form it are the same terms - so these lines are that
    # calc's with ``psi_bar`` in place of ``psi``.
    t_da = a_mix * (psi_bar - 2.0)
    t_db = -b_mix
    t_dc = coefficient * (psi_bar - 1.0)
    d_f_dz = PR.df_dz(z, a_mix, b_mix)
    t_dfdt = PR.t_dfdt(z, a_mix, b_mix, t_da, t_db)
    t_dz = -t_dfdt / d_f_dz
    n_plus = z + PR.delta1 * b_mix
    n_minus = z + PR.delta2 * b_mix
    t_di = (t_dz + PR.delta1 * t_db) / n_plus - (t_dz + PR.delta2 * t_db) / n_minus
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
    g = PR.helmholtz_g(b_sum)
    return -total * math.log(1.0 - b_sum) - (a_sum / (PR.delta_diff * b_sum)) * g


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
    g = PR.helmholtz_g(b_sum)
    g_prime = PR.helmholtz_g_prime(b_sum)
    g_second = PR.helmholtz_g_second(b_sum)
    half_delta_diff = PR.half_delta_diff
    delta_diff = PR.delta_diff

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
