"""``eos.pt_flash`` - the isothermal two-phase flash.

Spec: ``specs/models/eos/pt_flash.toml``

The two-phase split of a mixture at a fixed temperature and pressure: Wilson
K-value estimates, then successive substitution, with Rachford-Rice bisected each
iteration. A *model* rather than a calculation - what the spec pins down is the
procedure, and this module reads the procedure from ``azoth._models_gen`` rather than
choosing it.
"""

from __future__ import annotations

import math
from typing import Any

from azoth import _models_gen
from azoth.core.errors import InvalidInputError, SolverNotConvergedError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import Phase, PtFlashResult
from azoth.core.units import Q, input_to_si
from azoth.core.warnings import Warning
from azoth.eos.mixture import Mixture
from azoth.eos.reference._flash_newton import split_at as flash_newton_split
from azoth.eos.reference._flash_newton import step as flash_newton_step
from azoth.eos.reference._mixture_state import (
    compositions,
    is_trivial,
    negative_flash_warning,
    phase_state,
    rachford_rice,
    rachford_rice_bounds,
    reduced_parameters,
    rms_delta,
    single_phase_warning,
    trivial_warning,
    wilson_k,
)

MODEL_ID = "eos.pt_flash"

#: The ``ln K`` below which the iteration has found the trivial solution.
#:
#: ``|ln K_i| < 1e-6`` for every ``i`` means the two phases have converged onto the
#: feed. It is compared against the *K-values* and never against ``beta``, which is
#: indeterminate there - see the spec's correction 2 for the feed a ``beta`` test
#: would have mislabelled.
#:
#: **The threshold must not sit on the scale of the answer, and this one did.** The
#: second-order scheme lands on ``|ln K|`` within ``1e-8`` of one at a single-phase state,
#: so ``1e-8`` made the trivial declaration a coin flip on the last bit: a one-ulp change
#: in the cubic's first term moved that state from ``trivial`` at 128 iterations to no
#: convergence at all. A split has ``|ln K|`` of order one, so three further decades cost
#: nothing. **The Rust kernel's ``TRIVIAL_TOLERANCE`` is this number** - the two are one
#: constant in two languages and have to move together, which nothing enforced until they
#: were found differing by two decades.
TRIVIAL_TOLERANCE = 1.0e-06


def settled_tolerance(algorithm: dict[str, Any], from_newton: bool) -> float:
    """The tolerance the loop settled against, which depends on which scheme finished it.

    The two schemes measure different things - the outer one the change in ``ln K`` and
    the fallback the step it took - so they carry their own stopping rules, and the
    comparison at the end has to be against the one that actually ran.
    """
    fallback = algorithm.get("fallback")
    if from_newton and fallback is not None:
        return float(fallback["algorithm"]["tolerance"])
    return float(algorithm["tolerance"])


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
    reduced = reduced_parameters(mixture, t_si, p_si)
    warnings.extend(reduced.warnings)

    algorithm = spec["algorithm"]
    inner = algorithm["inner"]
    kij = mixture.kij

    k = wilson_k(mixture, t_si, p_si)
    iterations = 0
    residual = math.nan
    # The vapour mole numbers the fallback Newton works in, seeded from the last
    # successive-substitution iterate so the handover starts where the outer scheme
    # had got to rather than from the feed.
    u = [0.0] * n
    newton_steps = 0
    from_newton = False
    # Whether the last outer iterate was a two-phase split. The fallback works in the
    # vapour mole numbers `u = beta y`, which describe two phases only while `beta` is
    # inside `(0, 1)` - a negative flash has no such `u`, and every trial the line
    # search makes at one is infeasible and halved through.
    splits = False
    fallback = algorithm.get("fallback")
    settled: str | None = None

    for step in range(1, algorithm["max_iterations"] + 1):
        iterations = step

        # The handover. NeqSim's `TPflash` runs its own scheme until `activeNewtonLimit`
        # steps have passed and replaces it with the second-order one from then on. Two
        # things end the handover without ending the solve: a step whose linear system
        # is singular, and the second-order scheme's own attempt budget. Both hand the
        # iteration back to the outer scheme - the cap that ends a solve is the outer
        # one, never the fallback's.
        if (
            fallback is not None
            and splits
            and step >= fallback["after"]
            and newton_steps < fallback["algorithm"]["max_iterations"]
        ):
            second_order = fallback["algorithm"]
            newton_steps += 1
            try:
                stepped = flash_newton_step(
                    mixture,
                    reduced,
                    list(z),
                    u,
                    second_order,
                    temperature=t_si,
                    pressure=p_si,
                )
            except SolverNotConvergedError:
                from_newton = False
            else:
                u = stepped.u
                residual = stepped.residual
                from_newton = True
                # The K-values follow the second-order iterate, so that a handover back
                # at `max_iterations` resumes from the state the solve is actually in
                # rather than re-running the outer scheme from where it handed over.
                x_newton, y_newton = flash_newton_split(list(z), u)
                k = [yi / xi for yi, xi in zip(y_newton, x_newton, strict=True)]
                if residual <= second_order["tolerance"]:
                    break
                continue

        # Both guards run *before* the Rachford-Rice solve, because both are cases
        # where the bracket is a division by zero: `K_i = 1` puts a pole at
        # infinity, and K-values all on one side of one put a bracket end there.
        if is_trivial(k, TRIVIAL_TOLERANCE):
            settled = "trivial"
            break
        # A K-value that is not finite and positive is the iteration diverging rather
        # than a state to diagnose. Checked because every path below treats the
        # K-values as a composition ratio, and `nan > 1.0` being false would
        # otherwise report a diverged iteration as a single-phase feed.
        if any(not math.isfinite(value) or value <= 0.0 for value in k):
            raise SolverNotConvergedError(iterations, residual, algorithm["tolerance"])
        bounds = rachford_rice_bounds(k)
        if bounds is None:
            settled = "all_vapour" if all(value > 1.0 for value in k) else "all_liquid"
            break

        beta = rachford_rice(list(z), k, bounds, inner["tolerance"], inner["max_iterations"])
        x, y = compositions(list(z), k, beta)
        liquid_state = phase_state(reduced, kij, x, liquid=True)
        vapour_state = phase_state(reduced, kij, y, liquid=False)

        ln_k_new = [
            lp - lv for lp, lv in zip(liquid_state.ln_phi, vapour_state.ln_phi, strict=True)
        ]
        residual = rms_delta(ln_k_new, k)
        k = [math.exp(value) for value in ln_k_new]
        u = [beta * value for value in y]
        splits = 0.0 < beta < 1.0

        if residual <= algorithm["tolerance"]:
            break

    beta_out: float | None
    if settled == "trivial" or (settled is None and is_trivial(k, TRIVIAL_TOLERANCE)):
        # Stated *exactly* - x = y = z and K = 1 - rather than as the last iterate
        # that approached it, which is a function of where the bisection stopped on
        # an identically-zero function and is not reproducible.
        k = [1.0] * n
        x, y = list(z), list(z)
        phase = Phase.TRIVIAL
        beta_out = None
        warnings.append(trivial_warning())
    elif settled is not None:
        x, y = list(z), list(z)
        phase = Phase.ALL_VAPOUR if settled == "all_vapour" else Phase.ALL_LIQUID
        beta_out = None
        warnings.append(single_phase_warning(phase))
    elif residual > settled_tolerance(algorithm, from_newton):
        raise SolverNotConvergedError(
            iterations, residual, settled_tolerance(algorithm, from_newton)
        )
    elif from_newton:
        # The second-order scheme converged, and it converges on its own measure - the
        # relative step - rather than on the outer scheme's. Its answer is a state in
        # `u`, so the phases and the K-values are read back off it.
        beta_out = sum(u)
        x, y = flash_newton_split(list(z), u)
        k = [yi / xi for yi, xi in zip(y, x, strict=True)]
        if beta_out < 0.0:
            phase = Phase.ALL_LIQUID
        elif beta_out > 1.0:
            phase = Phase.ALL_VAPOUR
        else:
            phase = Phase.TWO_PHASE
        if phase is not Phase.TWO_PHASE:
            warnings.append(negative_flash_warning(phase, beta_out))
    else:
        bounds = rachford_rice_bounds(k)
        if bounds is None:  # pragma: no cover - guarded above
            raise SolverNotConvergedError(iterations, residual, algorithm["tolerance"])
        beta_out = rachford_rice(list(z), k, bounds, inner["tolerance"], inner["max_iterations"])
        x, y = compositions(list(z), k, beta_out)
        if beta_out < 0.0:
            phase = Phase.ALL_LIQUID
        elif beta_out > 1.0:
            phase = Phase.ALL_VAPOUR
        else:
            phase = Phase.TWO_PHASE
        if phase is not Phase.TWO_PHASE:
            warnings.append(negative_flash_warning(phase, beta_out))

    liquid_state = phase_state(reduced, kij, x, liquid=True)
    vapour_state = phase_state(reduced, kij, y, liquid=False)

    return PtFlashResult(
        beta=beta_out,
        x=tuple(x),
        y=tuple(y),
        k=tuple(k),
        ln_phi_liquid=tuple(liquid_state.ln_phi),
        ln_phi_vapour=tuple(vapour_state.ln_phi),
        z_liquid=liquid_state.z,
        z_vapour=vapour_state.z,
        min_t_over_tc=min_t_over_tc,
        phase=phase,
        iterations=iterations,
        residual=residual,
        warnings=tuple(warnings),
    )
