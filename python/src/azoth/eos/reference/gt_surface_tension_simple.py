"""``eos.gt_surface_tension_simple`` - NeqSim's ``GTSurfaceTensionSimple``.

Spec: ``specs/models/eos/gt_surface_tension_simple.toml``. The class walks a linear
density path like LGT, but at each step it solves for the composition gradients that
equalise the chemical potentials - an inner Newton on the ``dmu/dn`` matrix.
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from typing import Any

from azoth.core.errors import InvalidInputError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import GtSurfaceTensionSimpleResult
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning, WarningCode
from azoth.eos.reference._mixture_state import (
    helmholtz_gradient,
    helmholtz_hessian,
    reduced_parameters,
)
from azoth.eos.reference.pt_flash import pt_flash

CALC_ID = "eos.gt_surface_tension_simple"

R = 8.3144621
STEPS = 200
INNER_STEPS = 1200
INNER_TOLERANCE = 1.0e-15


def _solve(a: Sequence[Sequence[float]], b: Sequence[float]) -> list[float]:
    """Gaussian elimination with partial pivoting for the small inner system."""
    count = len(b)
    if count == 0:
        return []
    matrix = [list(row) for row in a]
    rhs = list(b)
    for column in range(count):
        pivot = max(range(column, count), key=lambda row: abs(matrix[row][column]))
        if abs(matrix[pivot][column]) < 1.0e-300:
            raise InvalidInputError(
                "influence_parameters",
                "the GT-simple inner system is singular; the class's LU solver throws there",
            )
        matrix[column], matrix[pivot] = matrix[pivot], matrix[column]
        rhs[column], rhs[pivot] = rhs[pivot], rhs[column]
        for row in range(column + 1, count):
            factor = matrix[row][column] / matrix[column][column]
            for k in range(column, count):
                matrix[row][k] -= factor * matrix[column][k]
            rhs[row] -= factor * rhs[column]
    solution = [0.0] * count
    for row in range(count - 1, -1, -1):
        total = rhs[row]
        for k in range(row + 1, count):
            total -= matrix[row][k] * solution[k]
        solution[row] = total / matrix[row][row]
    return solution


def _state(
    mixture: Any,
    reduced: Any,
    a_abs: Sequence[float],
    b_abs: Sequence[float],
    delta1: float,
    delta2: float,
    density: Sequence[float],
    rt: float,
    reference_pressure: float,
) -> tuple[float, list[float], list[list[float]]]:
    count = len(density)
    if any(not math.isfinite(value) or value <= 0.0 for value in density):
        raise InvalidInputError(
            "z",
            "the GT-simple profile's density vector is not a state: a component density is "
            "non-positive or not finite. The class's local system returns a NaN potential "
            "there, its loop breaks, and it returns the partial integral accumulated so far",
        )
    total = sum(density)
    x = [value / total for value in density]
    volume = 1.0 / total
    b_mix = sum(x[i] * b_abs[i] for i in range(count))
    a_mix = sum(
        x[i] * x[j] * (1.0 - mixture.kij[i][j]) * math.sqrt(a_abs[i] * a_abs[j])
        for i in range(count)
        for j in range(count)
    )
    if not math.isfinite(volume) or volume <= b_mix:
        raise InvalidInputError(
            "z",
            f"the GT-simple profile's density vector is not a state: its molar volume "
            f"{volume} is at or below the mixture covolume {b_mix}. The class evaluates its "
            f"local system there anyway and the potentials become NaN",
        )
    pressure = rt / (volume - b_mix) - a_mix / (
        (volume + delta1 * b_mix) * (volume + delta2 * b_mix)
    )
    compressibility = pressure * volume / rt
    ratio = pressure / reference_pressure
    reduced_step = reduced._replace(
        a=[value * ratio for value in reduced.a],
        b=[value * ratio for value in reduced.b],
    )
    gradient = helmholtz_gradient(reduced_step, mixture.kij, x, compressibility)
    hessian = helmholtz_hessian(reduced_step, mixture.kij, x, compressibility)
    mu = [rt * (math.log(density[i] * rt) + gradient[i]) for i in range(count)]
    dmu = [
        [rt * hessian[i][k] / total + (rt / density[i] if i == k else 0.0) for k in range(count)]
        for i in range(count)
    ]
    return pressure, mu, dmu


def _initial_gradient(
    dmu: Sequence[Sequence[float]],
    sqrt_influence: Sequence[float],
    reference: int,
    delta_rho: Sequence[float],
    density: Sequence[float],
) -> list[float]:
    count = len(dmu)
    gradient = [0.0] * count
    gradient[reference] = delta_rho[reference]
    others = [index for index in range(count) if index != reference]
    if not others:
        return gradient
    b = [
        sqrt_influence[reference] * dmu[i][reference]
        - sqrt_influence[i] * dmu[reference][reference]
        for i in others
    ]
    f = [
        [
            sqrt_influence[i] * dmu[reference][k] - sqrt_influence[reference] * dmu[i][k]
            for k in others
        ]
        for i in others
    ]
    # The class's initial step solves the transposed system; the inner loop does not.
    transposed = [[f[j][i] for j in range(len(f))] for i in range(len(f))]
    solved = _solve(transposed, b)
    for position, i in enumerate(others):
        step = solved[position] * delta_rho[reference]
        if abs(step / density[i]) > 0.1:
            gradient[i] = math.copysign(density[i], solved[position])
        else:
            gradient[i] = step
    return gradient


def gt_surface_tension_simple(
    components: Sequence[str],
    eos: str,
    influence_parameters: Sequence[float],
    reference_component: float,
    pair_interactions: Sequence[Sequence[float]],
    T: Q,
    P: Q,
    z: Sequence[float],
) -> GtSurfaceTensionSimpleResult:
    """The simplified gradient theory's interface tension of a flashed fluid.

    Args:
        components: the fluid's substances, by name.
        eos: the cubic the fluid runs.
        influence_parameters: each component's gradient-theory influence parameter.
        reference_component: the 0-based index of the highest-boiling component.
        pair_interactions: the class's water/MEG ``interact`` matrix; its own branch
            writes zero for every pair.
        T: absolute temperature.
        P: absolute pressure.
        z: the overall mole fractions.

    Returns:
        The interface tension, in ``N/m``.

    Raises:
        InvalidInputError: for a single-phase feed, a singular inner system, or a profile
            state that is not a state.
        OutOfRangeError: if an input is outside its declared range.

    Example:
        >>> import azoth
        >>> q = azoth.ureg.Quantity
        >>> r = gt_surface_tension_simple(
        ...     ["methane", "n-butane"], "srk",
        ...     [9.634812734820117e-21, 2.7571945977326497e-19], 1.0,
        ...     [[0.0, 0.0], [0.0, 0.0]], q(300.0, "K"), q(20.0, "bar"),
        ...     [0.5, 0.5])
        >>> abs(r.sigma.magnitude - 0.00951973185903182) < 1e-6
        True
    """
    from azoth._models_gen import model
    from azoth.eos import components as databank

    spec: dict[str, Any] = model(CALC_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    count = len(z)
    if (
        len(components) != count
        or len(influence_parameters) != count
        or len(pair_interactions) != count
        or any(len(row) != count for row in pair_interactions)
    ):
        raise InvalidInputError(
            "influence_parameters",
            f"{count} component(s) against {len(influence_parameters)} influence "
            f"parameter(s), {len(pair_interactions)} interaction row(s) and "
            f"{len(components)} name(s): every one is one entry per component",
        )
    if (
        not math.isfinite(reference_component)
        or reference_component < 0.0
        or reference_component != int(reference_component)
        or int(reference_component) >= count
    ):
        raise InvalidInputError(
            "reference_component",
            f"{reference_component} is not a whole index in 0..{count}",
        )
    reference = int(reference_component)

    temperature = input_to_si(spec, "T", T)
    pressure = input_to_si(spec, "P", P)
    influences = [float(value) for value in influence_parameters]
    apply_checks(
        checks.on_input,
        {
            "T": temperature,
            "P": pressure,
            "influence_parameters": min(influences) if influences else None,
            "z": min(float(value) for value in z) if z else None,
        }.get,
        warnings,
    )

    mixture, _ = databank.mixture_of(list(components), eos=eos)
    flash = pt_flash(
        mixture, from_si(temperature, "K"), from_si(pressure, "Pa"), [float(v) for v in z]
    )
    if flash.phase != "two_phase":
        raise InvalidInputError(
            "z",
            "the feed is single phase at this state, so there is no interface to have a "
            "tension; a gradient-theory profile runs between two phases",
        )
    if flash.z_vapour <= 0.0 or flash.z_liquid <= 0.0:
        raise InvalidInputError(
            "z", "a phase's compressibility root is not positive, so it has no density"
        )
    for fraction in list(flash.y) + list(flash.x):
        if fraction <= 0.0:
            raise InvalidInputError(
                "z",
                "a component is absent from one of the phases, so its chemical potential "
                "would be `ln(0)`",
            )

    rt = R * temperature
    rho_vapour = [flash.y[i] * pressure / (flash.z_vapour * rt) for i in range(count)]
    rho_liquid = [flash.x[i] * pressure / (flash.z_liquid * rt) for i in range(count)]
    reduced = reduced_parameters(mixture, temperature, pressure)
    a_abs = [reduced.a[i] * rt * rt / pressure for i in range(count)]
    b_abs = [reduced.b[i] * rt / pressure for i in range(count)]
    sqrt_influence = [math.sqrt(value) for value in influences]

    _, mu_equi, _ = _state(
        mixture,
        reduced,
        a_abs,
        b_abs,
        mixture.cubic.delta1,
        mixture.cubic.delta2,
        rho_vapour,
        rt,
        pressure,
    )
    delta_rho = [(rho_liquid[i] - rho_vapour[i]) / STEPS for i in range(count)]
    if delta_rho[reference] == 0.0:
        raise InvalidInputError(
            "reference_component",
            "the reference component's density is the same in both phases, so the profile "
            "has no integration variable",
        )

    density = list(rho_vapour)
    sigma = 0.0
    warned = False
    for step in range(1, STEPS):
        _, _, dmu = _state(
            mixture,
            reduced,
            a_abs,
            b_abs,
            mixture.cubic.delta1,
            mixture.cubic.delta2,
            density,
            rt,
            pressure,
        )
        gradient = _initial_gradient(dmu, sqrt_influence, reference, delta_rho, density)
        converged = False
        for _ in range(INNER_STEPS):
            candidate = [density[i] + gradient[i] for i in range(count)]
            _, mu, dmu_inner = _state(
                mixture,
                reduced,
                a_abs,
                b_abs,
                mixture.cubic.delta1,
                mixture.cubic.delta2,
                candidate,
                rt,
                pressure,
            )
            others = [index for index in range(count) if index != reference]
            b = [
                -sqrt_influence[i] * (mu_equi[reference] - mu[reference])
                + sqrt_influence[reference] * (mu_equi[i] - mu[i])
                for i in others
            ]
            f = [
                [
                    -sqrt_influence[i] * dmu_inner[reference][k]
                    + sqrt_influence[reference] * dmu_inner[i][k]
                    for k in others
                ]
                for i in others
            ]
            x = _solve(f, b)
            error = sum(abs(value) for value in x) / sum(candidate)
            for position, i in enumerate(others):
                gradient[i] += x[position]
            if error <= INNER_TOLERANCE:
                converged = True
                break
        if not converged and not warned:
            warned = True
            warnings.append(
                Warning(
                    code=WarningCode.SOLVER_NOT_CONVERGED,
                    message=(
                        f"the inner mu-equality loop reached {INNER_STEPS} steps at step "
                        f"{step} without meeting its {INNER_TOLERANCE} tolerance; the last "
                        f"iterate is integrated, as the class does"
                    ),
                )
            )
        final = [density[i] + gradient[i] for i in range(count)]
        final_pressure, final_mu, _ = _state(
            mixture,
            reduced,
            a_abs,
            b_abs,
            mixture.cubic.delta1,
            mixture.cubic.delta2,
            final,
            rt,
            pressure,
        )
        kappa = 0.0
        for i in range(count):
            ki = 1.0 if i == reference else gradient[i] / delta_rho[reference]
            for k in range(count):
                kk = 1.0 if k == reference else gradient[k] / delta_rho[reference]
                kappa += (
                    sqrt_influence[i]
                    * sqrt_influence[k]
                    * ki
                    * kk
                    * (1.0 - float(pair_interactions[i][k]))
                )
        mu_times_rho = sum(final[i] * (final_mu[i] - mu_equi[i]) for i in range(count)) - (
            final_pressure - pressure
        )
        if kappa <= 0.0 or mu_times_rho <= 0.0:
            raise InvalidInputError(
                "z",
                f"the gradient-theory integrand is not positive at step {step} "
                f"(kappa = {kappa}, mu_times_rho = {mu_times_rho}); the class takes its "
                f"square root and stops at the NaN, and this port refuses a partial integral",
            )
        sigma += math.sqrt(2.0 * kappa * mu_times_rho) * delta_rho[reference]
        density = final

    return GtSurfaceTensionSimpleResult(sigma=from_si(abs(sigma), "N/m"), warnings=tuple(warnings))
