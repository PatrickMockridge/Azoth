"""``eos.stability_test`` - the tangent-plane stability test.

Whether a feed at a fixed temperature and pressure is stable as a single phase. The
question ``eos.pt_flash`` cannot ask: successive substitution finds *a* stationary
point, and a flash converging to ``x = y = z`` has shown that its own starting point
was not a split, not that the feed is single phase.

Spec: ``specs/models/eos/stability_test.yaml``, which carries the criterion, the two
Wilson trials and what they cost, the root the feed is placed on, and the iteration
cap.

Michelsen's tangent-plane distance, at each of two trial compositions: the feed is
unstable if either trial reaches a stationary point below the tangent plane. An
**unconverged trial raises** rather than being discarded - a tangent-plane distance
bounds stability only at a stationary point, and discarding a partial one would turn
"could not tell" into "stable".
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
    ``-n ln V``. **``V`` is not the same at two roots**, so the ideal part is not the
    same either: reducing with ``sum(n) = 1`` and ``V = Z R T / P`` leaves only
    ``-ln Z`` differing at one ``(T, P, n)``, and the comparison is therefore
    ``A^R / RT - ln Z + Z`` rather than ``A^R / RT + Z``. See the spec's correction 2
    for the measurement behind the wrong form.

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

    ``liquid`` names which root the trial's phase claims: the smallest for the
    liquid-like seed, the largest for the vapour-like one. Selecting by *ordering* is
    the rule :func:`azoth.eos.pr_z_factor` fixes, and both implementations follow it
    rather than re-deriving a root per iteration.

    Returns:
        ``(w, tm, iterations, residual)``.

    Raises:
        SolverNotConvergedError: if the iteration hits its cap, or if a mole number
            leaves the representable range. **Not** a trial discarded and the feed
            called stable: see the module documentation.
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
