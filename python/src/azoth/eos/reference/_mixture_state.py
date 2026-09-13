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
from typing import Any

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


def reduced_parameters(
    mixture: Mixture, temperature: float, pressure: float
) -> tuple[list[float], list[float], list[Warning]]:
    """``(A_i, B_i, warnings)`` for every component at a state.

    These depend on ``T`` and ``P`` alone, not on the composition, which is why they
    are computed once per solve rather than once per iteration: ``A_i`` and ``B_i``
    are the same numbers for the liquid and the vapour, and only the composition
    re-weights them.
    """
    a: list[float] = []
    b: list[float] = []
    warnings: list[Warning] = []
    for component in mixture.components:
        kappa = pr_kappa(component.omega)
        warnings.extend(kappa.warnings)
        ab = pr_alpha_ab(
            kappa.kappa,
            temperature / component.Tc.to_base_units().magnitude,
            pressure / component.Pc.to_base_units().magnitude,
        )
        warnings.extend(ab.warnings)
        a.append(ab.a_reduced)
        b.append(ab.b_reduced)
    return a, b, warnings


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
    a: list[float],
    b: list[float],
    kij: tuple[tuple[float, ...], ...],
    x: list[float],
    *,
    liquid: bool,
) -> tuple[float, list[float]]:
    """``(z, ln_phi)`` for one phase at a composition.

    ``z`` is selected by *ordering* - the smallest admissible root for the liquid,
    the largest for the vapour - never by an initial guess, which is the rule
    ``eos.pr_z_factor`` fixes.
    """
    n = len(x)
    a_mix, b_mix = mixture_parameters(a, b, kij, x)
    roots = pr_z_factor(a_mix, b_mix)
    z = roots.z_min if liquid else roots.z_max

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
    return z, ln_phi


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
