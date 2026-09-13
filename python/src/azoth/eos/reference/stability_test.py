"""``eos.stability_test`` - the tangent-plane stability test.

Whether a feed at a fixed temperature and pressure is stable as a single phase. The
question ``eos.pt_flash`` cannot ask: successive substitution finds *a* stationary
point, and a flash that converges to ``x = y = z`` has shown that its own starting
point was not a split, not that the feed is single phase.

Spec: ``specs/models/eos/stability_test.yaml``

# The criterion

Michelsen's tangent-plane distance. With the feed's chemical potentials as the
reference,

```text
d_i   = ln z_i + ln phi_i(z)          the reference, from the feed itself
w_i   = W_i / sum(W)                  a trial composition
ln W_i = d_i - ln phi_i(w)            iterated to a stationary point
tm    = 1 - sum(W)                    the distance at that point
```

and the feed is unstable if any trial reaches ``tm < 0``. A trial that converges to
the feed itself has ``sum(W) = 1`` and ``tm = 0`` - which is exactly why the
threshold is ``-1e-8`` and not zero, and why no separate trivial-solution test is
needed: a trivial trial cannot fall below it.

# Two trials, and what that costs

Both are seeded from Wilson K-values, one vapour-like (``z K``) and one liquid-like
(``z / K``). Wilson assumes gas-liquid equilibrium, so a feed unstable to a
*liquid-liquid* split can come back ``stable``. That is a real false negative,
stated in the spec's assumptions and not detectable from here. NeqSim carries a
supplementary pure-component trial for it, gated on model-family-name heuristics
that have no place in a spec; it is deliberately not ported.

# What the feed's own root is

``ln phi_i(z)`` needs the feed to sit on *a* root of the cubic, and in the
two-phase region there are three. The feed is placed at the **lower Gibbs energy**
of the admissible roots. The ideal part of ``G`` is identical at both - same
composition, same ``T`` and ``P`` - so comparing ``A^R / RT + Z`` decides it, and
that is what :func:`_feed_state` does.

This matters more than it looks. Taking the liquid root unconditionally would make
a superheated vapour's reference state the wrong one, and every ``tm`` would then be
measured from a state the feed is not in.
"""

from __future__ import annotations

import math

from azoth import _models_gen
from azoth.core.errors import InvalidInputError, SolverNotConvergedError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import StabilityTestResult, StabilityVerdict
from azoth.core.units import Q, input_to_si
from azoth.core.warnings import Warning
from azoth.eos.mixture import Mixture
from azoth.eos.reference._mixture_state import (
    PhaseState,
    ReducedParameters,
    helmholtz_energy,
    mixture_parameters,
    normalise,
    phase_state,
    phase_state_at,
    reduced_parameters,
    wilson_k,
)
from azoth.eos.reference.pr_z_factor import pr_z_factor

CALC_ID = "eos.stability_test"

#: Below this the feed is unstable. Negative rather than zero because a trivial
#: trial - one that converges to the feed - reaches `tm` of order 1e-16, and a
#: strict `tm < 0` would call a single-phase feed unstable on rounding alone.
TM_LIMIT = -1.0e-8

#: A trial whose mole numbers leave this range is diverging rather than converging.
#: `exp` of a large `ln W` overflows to infinity and the next normalisation is a
#: division by it, so the guard is what keeps "diverged" from being reported as a
#: number.
_LOG_W_CEILING = 700.0


def _ln(value: float) -> float:
    """``ln``, with zero giving negative infinity rather than raising.

    A component absent from the feed has ``z_i = 0`` and a reference potential of
    ``-inf``, which is the correct value: its trial mole number is zero and stays
    zero. ``math.log(0.0)`` raises instead, which would make an ordinary feed with
    a zero entry an error.
    """
    return math.log(value) if value > 0.0 else -math.inf


def _feed_state(
    reduced: ReducedParameters, kij: tuple[tuple[float, ...], ...], z: list[float]
) -> PhaseState:
    """The feed's phase state, on whichever admissible root has the lower Gibbs energy.

    ``G / RT = A / RT + Z`` with ``A = A^ideal + A^R``, and ``A^ideal`` carries
    ``-n ln V``. **``V`` is not the same at two roots**, so the ideal part is not
    the same either, and comparing ``A^R / RT + Z`` alone picks the wrong phase.

    Reducing it: with ``sum(n) = 1`` and ``V = Z R T / P``,

    ```text
    A^ideal / RT = sum_i n_i ln(n_i / V) ... = sum_i n_i ln n_i - ln Z - ln(R T / P) ...
    ```

    so between roots at one `(T, P, n)` only ``-ln Z`` differs. The comparison is
    therefore ``A^R / RT - ln Z + Z``.

    **Measured, because the wrong form is plausible and quiet.** Pure methane at
    150 K and 1 bar is superheated vapour - its saturation pressure there is about
    10 bar - and the two roots give:

    ```text
    root             A^R/RT      Z          A^R/RT + Z   A^R/RT - ln Z + Z
    liquid-like     -2.558039   0.003344   -2.554695    3.145800
    vapour-like     -0.015548   0.984493    0.968946    0.984574
    ```

    The wrong form picks the liquid and calls a plain vapour unstable; the right
    one picks the vapour. Nothing else about the model changes, which is what made
    it worth measuring rather than reasoning about.

    A single admissible root is the common case and is taken directly.
    """
    roots = pr_z_factor(*mixture_parameters(reduced.a, reduced.b, kij, z))
    candidates = [roots.z_min] if roots.z_min == roots.z_max else [roots.z_min, roots.z_max]

    best = None
    best_gibbs = math.inf
    for compressibility in candidates:
        residual = helmholtz_energy(reduced, kij, z, compressibility)
        gibbs = residual - math.log(compressibility) + compressibility
        if gibbs < best_gibbs:
            best_gibbs = gibbs
            best = compressibility
    assert best is not None  # `candidates` is never empty
    return phase_state_at(reduced, kij, z, best)


def _trial(
    reduced: ReducedParameters,
    kij: tuple[tuple[float, ...], ...],
    d: list[float],
    seed: list[float],
    *,
    liquid: bool,
    tolerance: float,
    cap: int,
) -> tuple[list[float], float, int, float]:
    """Iterate one trial to a stationary point of the tangent-plane distance.

    Returns ``(w, tm, iterations, residual)``.

    Raises:
        SolverNotConvergedError: if the iteration hits its cap, or if a mole number
            leaves the representable range. **Not** a trial discarded and the feed
            called stable: `tm` bounds stability only at a stationary point, so a
            value read from a partial iteration is a fact about the iteration rather
            than about the mixture, and treating it as evidence of stability is the
            failure this whole library is organised against.
    """
    w = normalise(seed)
    ln_w = [_ln(value) for value in w]
    residual = math.inf

    for step in range(1, cap + 1):
        state = phase_state(reduced, kij, w, liquid=liquid)
        ln_w_new = [di - lp for di, lp in zip(d, state.ln_phi, strict=True)]

        if any(value > _LOG_W_CEILING for value in ln_w_new):
            raise SolverNotConvergedError(step, residual, tolerance)

        # The unnormalised mole numbers, which is what `tm` is written against.
        totals = sum(math.exp(value) for value in ln_w_new)
        residual = math.sqrt(
            sum((new - old) ** 2 for new, old in zip(ln_w_new, ln_w, strict=True)) / len(ln_w_new)
        )
        w = normalise([math.exp(value) for value in ln_w_new])
        ln_w = ln_w_new

        if residual <= tolerance:
            return w, 1.0 - totals, step, residual

    raise SolverNotConvergedError(cap, residual, tolerance)


def stability_test(mixture: Mixture, T: Q, P: Q, z: list[float]) -> StabilityTestResult:
    """Whether a mixture at a temperature and pressure is stable as one phase.

    Args:
        mixture: the components and their interaction parameters.
        T: absolute temperature.
        P: absolute pressure.
        z: overall mole fractions. **Checked rather than renormalised**, as
            everywhere in this namespace.

    Returns:
        The verdict, the tangent-plane distance at each trial's stationary point,
        and the trial compositions.

    Raises:
        InvalidInputError: if ``z`` is the wrong length, has a negative entry, or
            does not sum to one.
        OutOfRangeError: if ``T`` or ``P`` is not positive.
        SolverNotConvergedError: if a trial hits its cap.

    Example:
        >>> import azoth
        >>> from azoth.eos import Component, mixture
        >>> q = azoth.ureg.Quantity
        >>> fluid = mixture(
        ...     [Component(q(190.56, "K"), q(4_599_200.0, "Pa"), 0.01142),
        ...      Component(q(425.12, "K"), q(3_796_000.0, "Pa"), 0.2002)],
        ...     kij={(0, 1): 0.05},
        ... )
        >>> r = stability_test(fluid, q(330.0, "K"), q(25.0, "bar"), [0.6, 0.4])
        >>> r.verdict
        'unstable'
    """
    spec = _models_gen.model(CALC_ID)
    algorithm = spec["algorithm"]

    t_si = input_to_si(spec, "T", T)
    p_si = input_to_si(spec, "P", P)

    n = len(mixture.components)
    if len(z) != n:
        raise InvalidInputError("z", f"a feed for {n} components has {len(z)} entries")
    for i, value in enumerate(z):
        if value < 0.0:
            raise InvalidInputError(
                "z", f"z[{i}] is {value} but a mole fraction cannot be negative"
            )
    if abs(sum(z) - 1.0) > 1.0e-09:
        raise InvalidInputError(
            "z",
            f"the feed's mole fractions sum to {sum(z)}, not to one. Renormalising it "
            f"here would make a composition error invisible in every number "
            f"downstream, so it is refused instead",
        )

    warnings: list[Warning] = []
    checks = checks_for(spec)

    min_t_over_tc = min(
        t_si / component.Tc.to_base_units().magnitude for component in mixture.components
    )
    apply_checks(checks.on_input, {"T": t_si, "P": p_si}.get, warnings)
    apply_checks(
        checks.derived,
        lambda name: min_t_over_tc if name == "min_t_over_tc" else None,
        warnings,
    )

    reduced = reduced_parameters(mixture, t_si, p_si)
    warnings.extend(reduced.warnings)
    kij = mixture.kij

    feed = _feed_state(reduced, kij, list(z))
    d = [_ln(zi) + lp for zi, lp in zip(z, feed.ln_phi, strict=True)]

    k = wilson_k(mixture, t_si, p_si)
    vapour_seed = [zi * ki for zi, ki in zip(z, k, strict=True)]
    liquid_seed = [zi / ki for zi, ki in zip(z, k, strict=True)]

    tolerance = algorithm["tolerance"]
    cap = algorithm["max_iterations"]

    w_rows: list[list[float]] = []
    tm: list[float] = []
    iterations: list[int] = []

    # Ordered vapour-like first, then liquid-like, and the order is the contract:
    # `tm` and `w` are positional and a caller reads `tm[0]` as the vapour-like
    # trial. Swapping them would silently relabel the two.
    for seed, is_liquid in ((vapour_seed, False), (liquid_seed, True)):
        row, distance, steps, _ = _trial(
            reduced,
            kij,
            d,
            seed,
            liquid=is_liquid,
            tolerance=tolerance,
            cap=cap,
        )
        w_rows.append(row)
        tm.append(distance)
        iterations.append(steps)

    unstable = any(distance < TM_LIMIT for distance in tm)

    return StabilityTestResult(
        verdict=StabilityVerdict.UNSTABLE if unstable else StabilityVerdict.STABLE,
        tm=tuple(tm),
        w=tuple(tuple(row) for row in w_rows),
        iterations=tuple(iterations),
        min_t_over_tc=min_t_over_tc,
        warnings=tuple(warnings),
    )
