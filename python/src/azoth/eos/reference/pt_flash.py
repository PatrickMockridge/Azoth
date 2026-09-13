"""``eos.pt_flash`` - the isothermal two-phase flash.

The two-phase split of a mixture at a fixed temperature and pressure, by successive
substitution from Wilson K-value estimates. A *model* rather than a calculation:
what the spec pins down is the procedure, and this module reads the procedure from
``azoth._models_gen`` rather than choosing it.

Spec: ``specs/models/eos/pt_flash.yaml``

# The algorithm

1. **Initialise** every ``K_i`` from Wilson's correlation,
   ``K_i = (Pc_i / P) * exp(5.373 * (1 + omega_i) * (1 - Tc_i / T))``.
2. **Solve Rachford-Rice** for the vapour fraction, by bisection on
   ``g(beta) = sum_i z_i (K_i - 1) / (1 + beta (K_i - 1))``.
3. **Compose** the two phases, ``x_i = z_i / (1 + beta (K_i - 1))`` and
   ``y_i = K_i x_i``.
4. **Evaluate** ``ln phi_i`` in each phase from the mixture form, at the liquid root
   for ``x`` and the vapour root for ``y``.
5. **Update** ``K_i = exp(ln phi_i^L - ln phi_i^V)`` and stop when the rms change in
   ``ln K`` meets the tolerance.
6. **Re-evaluate once** at the converged ``K``, so the returned ``beta``, ``x``, ``y``
   and ``K`` are one consistent state rather than one state and its predecessor's
   vapour fraction.

# The Rachford-Rice bracket, which is not the textbook one

``g`` has poles at ``1/(1 - K_i)``, is strictly decreasing between consecutive
poles, and jumps from ``-inf`` to ``+inf`` across each. So the leftmost and
rightmost intervals have no root, and every bounded interval between poles has
exactly one. Of those, the one with both phases positive is the interval on which
``1 + beta (K_i - 1) > 0`` for every ``i``:

.. code-block:: text

    lower = max over {i : K_i > 1} of 1/(1 - K_i)
    upper = min over {i : K_i < 1} of 1/(1 - K_i)

and **it exists only when the K-values straddle one**. The textbook statement
``1/(1 - K_max) < beta < 1/(1 - K_min)`` is what those two collapse to in that case,
and is an empty inverted interval otherwise - a bisection handed it walks to
whichever end it happens to reach and returns a number. See the spec's correction 1
for what each of the two wrong versions measured.

# The mixture fugacity coefficient is this layer's own arithmetic

``eos.pr_departure`` covers a *pure* component; the mixture form, which carries the
sum over ``x_j a_ij`` and the ``b_i / b_mix`` term, has no registered calculation
behind it because the registry is scalar and has no composition vector to hang one
on. It is checked by composition instead: at ``N = 1`` the cross-sum factor collapses
to 1 and this becomes exactly ``eos.pr_departure``, and at ``N = 2`` the mixture
parameters and the vapour fraction must reproduce ``eos.vdw1f_mix_binary`` and
``eos.rachford_rice_binary``. Both are asserted in the test files.

# No stability test

Successive substitution finds *a* stationary point. Whether the feed was stable is a
different question, and this model cannot ask it. Where the K-values straddle one
throughout and the feed is single phase, the iteration converges to ``x = y = z`` and
the result says ``TRIVIAL`` - the feed is single phase, and which one is not
something this function can tell you.
"""

from __future__ import annotations

import math

from azoth import _models_gen
from azoth.core.errors import InvalidInputError, SolverNotConvergedError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import Phase, PtFlashResult
from azoth.core.units import Q, input_to_si
from azoth.core.warnings import Warning, WarningCode
from azoth.eos.mixture import Mixture
from azoth.eos.reference.pr_alpha_ab import pr_alpha_ab
from azoth.eos.reference.pr_kappa import pr_kappa
from azoth.eos.reference.pr_z_factor import pr_z_factor

MODEL_ID = "eos.pt_flash"

#: The ``|ln K|`` below which the iteration has found the trivial solution.
TRIVIAL_TOLERANCE = 1.0e-08

#: Wilson's constant. Some sources print 5.37 and the paper is dated 1968 in some
#: and 1969 in others; the discrepancy is recorded in the spec's references rather
#: than resolved, because a reader meeting the other value needs to know it is the
#: same correlation and not a correction.
WILSON_CONSTANT = 5.373

_SQRT_2 = math.sqrt(2.0)


def _wilson_k(mixture: Mixture, temperature: float, pressure: float) -> list[float]:
    """Wilson's correlation for the initial K-values."""
    return [
        (component.Pc.to_base_units().magnitude / pressure)
        * math.exp(
            WILSON_CONSTANT
            * (1.0 + component.omega)
            * (1.0 - component.Tc.to_base_units().magnitude / temperature)
        )
        for component in mixture.components
    ]


def _rachford_rice_bounds(k: list[float]) -> tuple[float, float] | None:
    """The interval on which Rachford-Rice has its physical root, or ``None``.

    ``None`` is not a numerical failure: it is a proof that the feed has no
    two-phase solution at these K-values, since ``sum_i y_i = sum_i K_i x_i = 1``
    with ``sum_i x_i = 1`` is impossible when every K is on the same side of one.
    """
    lower = -math.inf
    upper = math.inf
    for value in k:
        if value > 1.0:
            lower = max(lower, 1.0 / (1.0 - value))
        elif value < 1.0:
            upper = min(upper, 1.0 / (1.0 - value))
        else:
            # ``K_i = 1`` exactly puts a pole at infinity and makes ``g`` degenerate.
            return None
    if lower == -math.inf or upper == math.inf:
        return None
    return (lower, upper)


def _rachford_rice(
    z: list[float], k: list[float], bounds: tuple[float, float], tolerance: float, cap: int
) -> float:
    """The vapour fraction that solves Rachford-Rice, by bisection."""

    def g(beta: float) -> float:
        return sum(zi * (ki - 1.0) / (1.0 + beta * (ki - 1.0)) for zi, ki in zip(z, k, strict=True))

    lower, upper = bounds
    for _ in range(cap):
        mid = 0.5 * (lower + upper)
        if upper - lower <= tolerance:
            break
        if g(mid) > 0.0:
            lower = mid
        else:
            upper = mid
    return 0.5 * (lower + upper)


def _mixture_parameters(
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


def _phase_state(
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
    a_mix, b_mix = _mixture_parameters(a, b, kij, x)
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


def _compositions(z: list[float], k: list[float], beta: float) -> tuple[list[float], list[float]]:
    """The two phases' compositions at a vapour fraction."""
    x = [zi / (1.0 + beta * (ki - 1.0)) for zi, ki in zip(z, k, strict=True)]
    y = [ki * xi for ki, xi in zip(k, x, strict=True)]
    return x, y


def _is_trivial(k: list[float]) -> bool:
    """Whether the K-values have converged onto the feed."""
    return all(abs(math.log(value)) < TRIVIAL_TOLERANCE for value in k)


def _rms_delta(ln_k_new: list[float], k: list[float]) -> float:
    """The rms change in ``ln K`` across one iteration."""
    total = sum((new - math.log(old)) ** 2 for new, old in zip(ln_k_new, k, strict=True))
    return math.sqrt(total / len(ln_k_new))


def pt_flash(mixture: Mixture, T: Q, P: Q, z: list[float]) -> PtFlashResult:
    """The isothermal two-phase flash of a mixture at a temperature and pressure.

    Args:
        mixture: the components and their interaction parameters.
        T: absolute temperature.
        P: absolute pressure.
        z: overall mole fractions. **Checked rather than renormalised** - silently
            rescaling a caller's composition would make their error invisible in a
            way that changes every number downstream.

    Returns:
        The split, or a statement of why there is not one. ``beta`` is ``None``
        whenever there is no vapour fraction to report; read ``phase``.

    Raises:
        InvalidInputError: if ``z`` is the wrong length, has a negative entry, or
            does not sum to one.
        OutOfRangeError: if ``T`` or ``P`` is not positive.
        SolverNotConvergedError: if the iteration hits its cap, or if a K-value
            becomes non-finite - which is the iteration diverging rather than a
            state to diagnose.

    Example:
        >>> import azoth
        >>> from azoth.eos import Component, mixture
        >>> q = azoth.ureg.Quantity
        >>> fluid = mixture(
        ...     [Component(q(190.56, "K"), q(4_599_200.0, "Pa"), 0.01142),
        ...      Component(q(425.12, "K"), q(3_796_000.0, "Pa"), 0.2002)],
        ...     kij={(0, 1): 0.05},
        ... )
        >>> r = azoth.eos.pt_flash(fluid, T=q(330.0, "K"), P=q(2.5e6, "Pa"), z=[0.6, 0.4])
        >>> round(r.beta, 9)
        0.844722027
        >>> r.phase.value
        'two_phase'
    """
    spec = _models_gen.model(MODEL_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    temperatures = [component.Tc.to_base_units().magnitude for component in mixture.components]
    t_si = input_to_si(spec, "T", T)
    p_si = input_to_si(spec, "P", P)

    apply_checks(checks.on_input, {"T": t_si, "P": p_si}.get, warnings)

    n = len(mixture)
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

    min_t_over_tc = min(t_si / tc for tc in temperatures)
    apply_checks(
        checks.derived,
        lambda name: min_t_over_tc if name == "min_t_over_tc" else None,
        warnings,
    )

    # The reduced parameters depend on `T` and `P` alone, so they are computed once
    # rather than once per iteration: `A_i` and `B_i` are the same numbers for both
    # phases, and only the composition re-weights them.
    a: list[float] = []
    b: list[float] = []
    for component in mixture.components:
        kappa = pr_kappa(component.omega)
        warnings.extend(kappa.warnings)
        ab = pr_alpha_ab(
            kappa.kappa,
            t_si / component.Tc.to_base_units().magnitude,
            p_si / component.Pc.to_base_units().magnitude,
        )
        warnings.extend(ab.warnings)
        a.append(ab.a_reduced)
        b.append(ab.b_reduced)

    algorithm = spec["algorithm"]
    inner = algorithm["inner"]
    kij = mixture.kij

    k = _wilson_k(mixture, t_si, p_si)
    iterations = 0
    residual = math.nan
    settled: str | None = None

    for step in range(1, algorithm["max_iterations"] + 1):
        iterations = step

        # Both guards run *before* the Rachford-Rice solve, because both are cases
        # where the bracket is a division by zero: `K_i = 1` puts a pole at
        # infinity, and K-values all on one side of one put a bracket end there.
        if _is_trivial(k):
            settled = "trivial"
            break
        # A K-value that is not finite and positive is the iteration diverging rather
        # than a state to diagnose. Checked because every path below treats the
        # K-values as a composition ratio, and `nan > 1.0` being false would
        # otherwise report a diverged iteration as a single-phase feed.
        if any(not math.isfinite(value) or value <= 0.0 for value in k):
            raise SolverNotConvergedError(iterations, residual, algorithm["tolerance"])
        bounds = _rachford_rice_bounds(k)
        if bounds is None:
            settled = "all_vapour" if all(value > 1.0 for value in k) else "all_liquid"
            break

        beta = _rachford_rice(list(z), k, bounds, inner["tolerance"], inner["max_iterations"])
        x, y = _compositions(list(z), k, beta)
        _, ln_phi_liquid = _phase_state(a, b, kij, x, liquid=True)
        _, ln_phi_vapour = _phase_state(a, b, kij, y, liquid=False)

        ln_k_new = [lp - lv for lp, lv in zip(ln_phi_liquid, ln_phi_vapour, strict=True)]
        residual = _rms_delta(ln_k_new, k)
        k = [math.exp(value) for value in ln_k_new]

        if residual <= algorithm["tolerance"]:
            break

    beta_out: float | None
    if settled == "trivial" or (settled is None and _is_trivial(k)):
        # Stated *exactly* - x = y = z and K = 1 - rather than as the last iterate
        # that approached it, which is a function of where the bisection stopped on
        # an identically-zero function and is not reproducible.
        k = [1.0] * n
        x, y = list(z), list(z)
        phase = Phase.TRIVIAL
        beta_out = None
        warnings.append(_trivial_warning())
    elif settled is not None:
        x, y = list(z), list(z)
        phase = Phase.ALL_VAPOUR if settled == "all_vapour" else Phase.ALL_LIQUID
        beta_out = None
        warnings.append(_single_phase_warning(phase))
    elif residual > algorithm["tolerance"]:
        raise SolverNotConvergedError(iterations, residual, algorithm["tolerance"])
    else:
        bounds = _rachford_rice_bounds(k)
        if bounds is None:  # pragma: no cover - guarded above
            raise SolverNotConvergedError(iterations, residual, algorithm["tolerance"])
        beta_out = _rachford_rice(list(z), k, bounds, inner["tolerance"], inner["max_iterations"])
        x, y = _compositions(list(z), k, beta_out)
        if beta_out < 0.0:
            phase = Phase.ALL_LIQUID
        elif beta_out > 1.0:
            phase = Phase.ALL_VAPOUR
        else:
            phase = Phase.TWO_PHASE
        if phase is not Phase.TWO_PHASE:
            warnings.append(_negative_flash_warning(phase, beta_out))

    z_liquid, ln_phi_liquid = _phase_state(a, b, kij, x, liquid=True)
    z_vapour, ln_phi_vapour = _phase_state(a, b, kij, y, liquid=False)

    return PtFlashResult(
        beta=beta_out,
        x=tuple(x),
        y=tuple(y),
        k=tuple(k),
        ln_phi_liquid=tuple(ln_phi_liquid),
        ln_phi_vapour=tuple(ln_phi_vapour),
        z_liquid=z_liquid,
        z_vapour=z_vapour,
        min_t_over_tc=min_t_over_tc,
        phase=phase,
        iterations=iterations,
        residual=residual,
        warnings=tuple(warnings),
    )


def _trivial_warning() -> Warning:
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


def _single_phase_warning(phase: Phase) -> Warning:
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


def _negative_flash_warning(phase: Phase, beta: float) -> Warning:
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
