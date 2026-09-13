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

# What keeps it honest

Not an assertion but a *reduction*: at ``N = 1`` the cross-sum factor collapses to 1
and :func:`phase_state` returns exactly ``eos.pr_departure``'s ``ln phi``, and at
``N = 2`` :func:`mixture_parameters` reproduces ``eos.vdw1f_mix_binary``. Both are
asserted in both languages, and the second is a cross-layer check no single-language
test can replace.

Both reductions are to within a couple of ulps rather than bit-identical, and the
reason is in ``azoth.eos.mixture``'s Rust counterpart: the registered binary kernel
evaluates its three terms longhand while this sums a double loop, so the two
associate differently. A bit-equality claim here would be a claim about summation
order rather than about the mixing rule.
"""

from __future__ import annotations

import math
from typing import Any, NamedTuple

from azoth.core.errors import OutOfRangeError
from azoth.core.warnings import Warning
from azoth.eos.mixture import Mixture
from azoth.eos.reference.pr_alpha_ab import pr_alpha_ab
from azoth.eos.reference.pr_kappa import pr_kappa
from azoth.eos.reference.pr_z_factor import pr_z_factor

#: Wilson's constant. Some sources print 5.37 and the paper is dated 1968 in some
#: and 1969 in others; the discrepancy is recorded in the model specs' references
#: rather than resolved, because a reader meeting the other value needs to know it
#: is the same correlation and not a correction.
WILSON_CONSTANT = 5.373

_SQRT_2 = math.sqrt(2.0)


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

    a: list[float]
    b: list[float]
    psi: list[float]
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


def reduced_parameters(mixture: Mixture, temperature: float, pressure: float) -> ReducedParameters:
    """``(A_i, B_i, psi_i, warnings)`` for every component at a state."""
    a: list[float] = []
    b: list[float] = []
    psi: list[float] = []
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
        sqrt_tr = math.sqrt(reduced_temperature)
        a.append(ab.a_reduced)
        b.append(ab.b_reduced)
        psi.append(-kappa.kappa * sqrt_tr / (1.0 + kappa.kappa * (1.0 - sqrt_tr)))
    return ReducedParameters(a=a, b=b, psi=psi, warnings=warnings)


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
    i_term = math.log((z + (1.0 + _SQRT_2) * b_mix) / (z + (1.0 - _SQRT_2) * b_mix))
    coefficient = a_mix / (2.0 * _SQRT_2 * b_mix)
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
    for i in range(n):
        for j in range(n):
            weight = x[i] * x[j] * (1.0 - kij[i][j]) * math.sqrt(a[i] * a[j])
            weight_total += weight
            weighted_psi += weight * 0.5 * (reduced.psi[i] + reduced.psi[j])
    psi_bar = weighted_psi / weight_total
    h_dep_rt = (z - 1.0) + coefficient * (psi_bar - 1.0) * i_term
    s_dep_r = h_dep_rt - sum(xi * lp for xi, lp in zip(x, ln_phi, strict=True))

    return PhaseState(
        a_mix=a_mix,
        b_mix=b_mix,
        z=z,
        ln_phi=ln_phi,
        h_dep_rt=h_dep_rt,
        s_dep_r=s_dep_r,
        psi_bar=psi_bar,
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
