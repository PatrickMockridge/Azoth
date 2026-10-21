"""``eos.gt_surface_tension`` - NeqSim's ``GTSurfaceTensionFullGT``.

Spec: ``specs/models/eos/gt_surface_tension.toml``. The class solves the full gradient
theory on a finite-difference grid: a damped Newton on the Euler-Lagrange equations, a
block-tridiagonal Jacobian, deferred correction, and an outer remeshing loop that doubles
the grid until the integral stops moving.
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from typing import Any

from azoth.core.errors import InvalidInputError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import GtSurfaceTensionResult
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning
from azoth.eos.reference._mixture_state import (
    helmholtz_gradient,
    helmholtz_hessian,
    reduced_parameters,
)
from azoth.eos.reference.tp_multiflash import tp_multiflash

CALC_ID = "eos.gt_surface_tension"

R = 8.3144621
NLEVEL = 3
NREFINEMENTS = 7
NEWTON_STEPS = 20
MAX_REL_CHANGE = 0.5
NEWTON_TOLERANCE = 1.0e-5
OUTER_TOLERANCE = 0.01
DOMAIN_HALF_WIDTH = 3.0
NM_SCALE = 1.0e18


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
) -> tuple[list[float], list[list[float]]]:
    """The T,V chemical potentials and their physical-density derivative."""
    count = len(density)
    if any(not math.isfinite(value) or value <= 0.0 for value in density):
        raise InvalidInputError(
            "z",
            "the full gradient theory's density vector is not a state: a component density "
            "is non-positive or not finite",
        )
    total = sum(density)
    x = [value / total for value in density]
    volume = 1.0 / total
    b_mix = sum(x[i] * b_abs[i] for i in range(count))
    a_mix = sum(
        x[i] * x[j] * (1.0 - float(mixture.kij[i][j])) * math.sqrt(a_abs[i] * a_abs[j])
        for i in range(count)
        for j in range(count)
    )
    if not math.isfinite(volume) or volume <= b_mix:
        raise InvalidInputError(
            "z",
            "the full gradient theory's density vector is not a state: its molar volume is "
            "at or below the mixture covolume",
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
    return mu, dmu


def _mat_mul(a: Sequence[Sequence[float]], b: Sequence[Sequence[float]]) -> list[list[float]]:
    count = len(a)
    return [
        [sum(a[i][k] * b[k][j] for k in range(count)) for j in range(count)] for i in range(count)
    ]


def _mat_vec(a: Sequence[Sequence[float]], v: Sequence[float]) -> list[float]:
    return [sum(a[i][k] * v[k] for k in range(len(v))) for i in range(len(a))]


def _invert(matrix: Sequence[Sequence[float]]) -> list[list[float]]:
    """The inverse of a small dense matrix, by Gaussian elimination with partial pivoting."""
    count = len(matrix)
    augmented = [
        list(matrix[i]) + [1.0 if i == j else 0.0 for j in range(count)] for i in range(count)
    ]
    for column in range(count):
        pivot = max(range(column, count), key=lambda row: abs(augmented[row][column]))
        if abs(augmented[pivot][column]) < 1.0e-300:
            raise InvalidInputError(
                "z",
                "the full gradient theory's block is singular; the class's band solver "
                "would fail on it",
            )
        augmented[column], augmented[pivot] = augmented[pivot], augmented[column]
        divisor = augmented[column][column]
        augmented[column] = [value / divisor for value in augmented[column]]
        for row in range(count):
            if row == column:
                continue
            factor = augmented[row][column]
            if factor == 0.0:
                continue
            augmented[row] = [
                augmented[row][k] - factor * augmented[column][k] for k in range(2 * count)
            ]
    return [row[count:] for row in augmented]


def _linspace(a: float, b: float, count: int) -> list[float]:
    dx = (b - a) / (count - 1)
    return [a + index * dx for index in range(count)]


def _density_gradient(rhomat: Sequence[Sequence[float]], h: float) -> list[list[float]]:
    count = len(rhomat[0])
    return [
        [(rhomat[i + 1][j] - rhomat[i][j]) / h for j in range(count)]
        for i in range(len(rhomat) - 1)
    ]


def _sigma_calc(
    h: float, drhodz: Sequence[Sequence[float]], cij: Sequence[Sequence[float]]
) -> float:
    count = len(cij)
    sigma = 0.0
    for j in range(count):
        for k in range(count):
            sigma += cij[j][k] * sum(row[j] * row[k] for row in drhodz)
    return sigma * h * 1.0e-9


def _std_integral(
    h: float, drhodz: Sequence[Sequence[float]], cij: Sequence[Sequence[float]]
) -> float:
    """The interface width, with the class's single assignment rather than a sum."""
    count = len(cij)
    total = 0.0
    z_sum = 0.0
    z2_sum = 0.0
    for index, row in enumerate(drhodz):
        z = index * h
        value = 0.0
        for m in range(count):
            for n in range(count):
                value = row[m] * cij[m][n] * row[n]
        total += value * h
        z_sum += value * ((z + h) * (z + h) - z * z) / 2.0
        z2_sum += value * ((z + h) ** 3 - z**3) / 3.0
    if total == 0.0:
        return 0.0
    mean = z_sum / total
    return math.sqrt(max(z2_sum / total - mean * mean, 0.0))


def _follow_alias(rhomat: list[list[float]], alias: tuple[int, int] | None) -> None:
    if alias is not None:
        source, shared = alias
        rhomat[shared] = list(rhomat[source])


def _direct_solve(
    residual: Sequence[Sequence[float]],
    jacobian: Sequence[Sequence[Sequence[float]]],
    cij: Sequence[Sequence[float]],
    h: float,
    rhomat: list[list[float]],
    count: int,
) -> None:
    """The block-tridiagonal Newton step, by the block Thomas algorithm, added in place."""
    ngrid = len(rhomat)
    h2 = h * h
    blocks = ngrid - 2
    rhs = [[0.0] * count for _ in range(blocks)]
    diagonal = [[[0.0] * count for _ in range(count)] for _ in range(blocks)]
    for block, i in enumerate(range(1, ngrid - 1)):
        for j in range(count):
            value = residual[i][j]
            for k in range(count):
                value += (
                    jacobian[i][j][k] * rhomat[i][k]
                    - cij[j][k] * (rhomat[i - 1][k] - 2.0 * rhomat[i][k] + rhomat[i + 1][k]) / h2
                )
            rhs[block][j] = -value
            for k in range(count):
                diagonal[block][j][k] = jacobian[i][j][k] + 2.0 / h2 * cij[j][k]
    off_diagonal = [[-value / h2 for value in row] for row in cij]

    modified_diagonal = [list(block) for block in diagonal]
    modified_rhs = [list(row) for row in rhs]
    for block in range(1, blocks):
        inverse = _invert(modified_diagonal[block - 1])
        product = _mat_mul(off_diagonal, inverse)
        correction = _mat_mul(product, off_diagonal)
        for j in range(count):
            for k in range(count):
                modified_diagonal[block][j][k] -= correction[j][k]
        rhs_correction = _mat_vec(product, modified_rhs[block - 1])
        for j in range(count):
            modified_rhs[block][j] -= rhs_correction[j]
    solution = [[0.0] * count for _ in range(blocks)]
    solution[blocks - 1] = _mat_vec(
        _invert(modified_diagonal[blocks - 1]), modified_rhs[blocks - 1]
    )
    for block in range(blocks - 2, -1, -1):
        rhs_block = list(modified_rhs[block])
        sub = _mat_vec(off_diagonal, solution[block + 1])
        for j in range(count):
            rhs_block[j] -= sub[j]
        solution[block] = _mat_vec(_invert(modified_diagonal[block]), rhs_block)

    for block, i in enumerate(range(1, ngrid - 1)):
        for j in range(count):
            rhomat[i][j] += solution[block][j]
            if not math.isfinite(rhomat[i][j]):
                raise InvalidInputError(
                    "z",
                    "the full gradient theory's Newton step left a non-finite density",
                )


def _newton(
    cij: Sequence[Sequence[float]],
    half_width: float,
    rhomat: list[list[float]],
    mixture: Any,
    reduced: Any,
    a_abs: Sequence[float],
    b_abs: Sequence[float],
    delta1: float,
    delta2: float,
    rt: float,
    reference_pressure: float,
    mu_equilibrium: Sequence[float],
    count: int,
) -> None:
    """One Newton solve for a grid, modifying the density matrix in place."""
    ngrid = len(rhomat)
    h = 2.0 * half_width / (ngrid - 1)
    jacobian = [[[0.0] * count for _ in range(count)] for _ in range(ngrid)]
    dmu = [[0.0] * count for _ in range(ngrid)]
    residual = [[0.0] * count for _ in range(ngrid)]
    for i in range(ngrid):
        mu, jac = _state(
            mixture, reduced, a_abs, b_abs, delta1, delta2, rhomat[i], rt, reference_pressure
        )
        for j in range(count):
            dmu[i][j] = mu[j] - mu_equilibrium[j]
            jacobian[i][j] = list(jac[j])
        for j in range(count):
            residual[i][j] = dmu[i][j] - sum(jacobian[i][j][k] * rhomat[i][k] for k in range(count))
    _add_deferred_correction(residual, dmu, ngrid, count)

    xgrid = _linspace(-half_width, half_width, ngrid)
    # The class's centre-of-gravity shift assigns row *references*, so the row it shifts
    # from and the boundary row become one array that the Newton update moves together.
    alias: tuple[int, int] | None = None
    for _ in range(NEWTON_STEPS):
        previous = [list(row) for row in rhomat]
        _direct_solve(residual, jacobian, cij, h, rhomat, count)
        _follow_alias(rhomat, alias)

        max_change = 0.0
        for i in range(1, ngrid - 1):
            for j in range(count):
                change = abs(rhomat[i][j] - previous[i][j])
                if change > max_change * abs(previous[i][j]):
                    max_change = change / abs(previous[i][j])
        damping = min(MAX_REL_CHANGE / max_change, 1.0) if max_change > 0.0 else 1.0
        for i in range(1, ngrid - 1):
            for j in range(count):
                rhomat[i][j] = damping * rhomat[i][j] + (1.0 - damping) * previous[i][j]
        _follow_alias(rhomat, alias)
        if damping * max_change < NEWTON_TOLERANCE:
            break

        drhodz = _density_gradient(rhomat, h)
        sum_z = 0.0
        total = 0.0
        for i, row in enumerate(drhodz):
            tmp = 0.0
            for j in range(count):
                for k in range(count):
                    tmp += row[j] * cij[j][k] * row[k]
            sum_z += 0.5 * (xgrid[i + 1] + xgrid[i]) * tmp
            total += tmp
        if total != 0.0:
            # `Math.round` is `floor(x + 0.5)`.
            correction = math.floor(sum_z / total / h + 0.5)
            if correction > 0:
                for i in range(1, ngrid - correction):
                    rhomat[i] = list(rhomat[i + correction])
                for i in range(ngrid - correction, ngrid - 1):
                    rhomat[i] = list(rhomat[ngrid - 1])
                alias = (ngrid - 1 - correction, ngrid - 1)
            elif correction < 0:
                shift = -correction
                for i in range(ngrid - 2, shift - 1, -1):
                    rhomat[i] = list(rhomat[i - shift])
                for i in range(1, shift):
                    rhomat[i] = list(rhomat[0])
                alias = (shift, 0)

        for i in range(1, ngrid - 1):
            mu, jac = _state(
                mixture, reduced, a_abs, b_abs, delta1, delta2, rhomat[i], rt, reference_pressure
            )
            for j in range(count):
                dmu[i][j] = mu[j] - mu_equilibrium[j]
                jacobian[i][j] = list(jac[j])
            for j in range(count):
                residual[i][j] = dmu[i][j] - sum(
                    jacobian[i][j][k] * rhomat[i][k] for k in range(count)
                )
        _add_deferred_correction(residual, dmu, ngrid, count)


def _add_deferred_correction(
    residual: list[list[float]], dmu: Sequence[Sequence[float]], ngrid: int, count: int
) -> None:
    for i in range(1, ngrid - 1):
        for j in range(count):
            residual[i][j] += (dmu[i - 1][j] - 2.0 * dmu[i][j] + dmu[i + 1][j]) / 12.0


def gt_surface_tension(
    components: Sequence[str],
    eos: str,
    influence_parameters: Sequence[float],
    T: Q,
    P: Q,
    z: Sequence[float],
) -> GtSurfaceTensionResult:
    """The full gradient theory's interface tension of a flashed fluid.

    Args:
        components: the fluid's substances, by name.
        eos: the cubic the fluid runs.
        influence_parameters: each component's gradient-theory influence parameter.
        T: absolute temperature.
        P: absolute pressure.
        z: the overall mole fractions.

    Returns:
        The interface tension, in ``N/m``.

    Raises:
        InvalidInputError: for a feed that is not two phases, a bulk state out of
            equilibrium, a singular block solve, or a profile that does not stay finite.
        OutOfRangeError: if an input is outside its declared range.

    Example:
        >>> import azoth
        >>> q = azoth.ureg.Quantity
        >>> r = gt_surface_tension(
        ...     ["methane", "n-butane"], "srk",
        ...     [9.634812734820117e-21, 2.7571945977326497e-19],
        ...     q(300.0, "K"), q(20.0, "bar"), [0.5, 0.5])
        >>> abs(r.sigma.magnitude - 0.00950471504464287) < 1e-8
        True
    """
    from azoth._models_gen import model
    from azoth.eos import components as databank

    spec: dict[str, Any] = model(CALC_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    count = len(z)
    if len(components) != count or len(influence_parameters) != count:
        raise InvalidInputError(
            "influence_parameters",
            f"{count} mole fraction(s) against {len(components)} component(s) and "
            f"{len(influence_parameters)} influence parameter(s): every one is one entry "
            f"per component",
        )

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
    flash = tp_multiflash(
        mixture, from_si(temperature, "K"), from_si(pressure, "Pa"), [float(v) for v in z]
    )
    if flash.phase_count != 2:
        if flash.phase_count == 1:
            raise InvalidInputError(
                "z",
                "the feed is single phase at this state, so there is no interface to have a "
                "tension; the full gradient theory solves a profile between two phases",
            )
        raise InvalidInputError(
            "z",
            f"the feed splits into {flash.phase_count} phases at this state, and the class "
            f"profiles its first two; this port refuses a state that is not two phases "
            f"rather than dropping one",
        )
    if flash.z_factor[0] <= 0.0 or flash.z_factor[1] <= 0.0:
        raise InvalidInputError(
            "z", "a phase's compressibility root is not positive, so it has no density"
        )
    for fraction in tuple(flash.x[0]) + tuple(flash.x[1]):
        if fraction <= 0.0:
            raise InvalidInputError(
                "z",
                "a component is absent from one of the phases. The class adds a trace amount "
                "and reflashes; this port refuses the state rather than perturbing the feed",
            )

    rt = R * temperature
    rho_ph1 = [flash.x[0][i] * pressure / (flash.z_factor[0] * rt) for i in range(count)]
    rho_ph2 = [flash.x[1][i] * pressure / (flash.z_factor[1] * rt) for i in range(count)]
    reduced = reduced_parameters(mixture, temperature, pressure)
    a_abs = [reduced.a[i] * rt * rt / pressure for i in range(count)]
    b_abs = [reduced.b[i] * rt / pressure for i in range(count)]
    delta1 = mixture.cubic.delta1
    delta2 = mixture.cubic.delta2

    equilibrium_1, _ = _state(mixture, reduced, a_abs, b_abs, delta1, delta2, rho_ph1, rt, pressure)
    equilibrium_2, _ = _state(mixture, reduced, a_abs, b_abs, delta1, delta2, rho_ph2, rt, pressure)
    for i in range(count):
        difference = abs(equilibrium_1[i] - equilibrium_2[i])
        scale = max(abs(equilibrium_1[i]), abs(equilibrium_2[i]))
        if difference > 1.0e-6 * max(scale, rt):
            raise InvalidInputError(
                "z",
                f"the flash's two phases do not have the same chemical potential for "
                f"component {i}: {equilibrium_1[i] - equilibrium_2[i]:.6e} J/mol apart, and "
                f"the class refuses a bulk state that is not solved",
            )

    cij = [
        [NM_SCALE * math.sqrt(influences[i] * influences[j]) for j in range(count)]
        for i in range(count)
    ]

    ngrid = (1 << NLEVEL) + 1
    half_width = DOMAIN_HALF_WIDTH
    half = ngrid // 2
    rhomat = [list(rho_ph1) if index < half else list(rho_ph2) for index in range(ngrid)]

    sigma = 0.0
    sigma_old = 9.9e99
    for refinement in range(NREFINEMENTS):
        _newton(
            cij,
            half_width,
            rhomat,
            mixture,
            reduced,
            a_abs,
            b_abs,
            delta1,
            delta2,
            rt,
            pressure,
            equilibrium_1,
            count,
        )
        h = 2.0 * half_width / (ngrid - 1)
        drhodz = _density_gradient(rhomat, h)
        sigma = _sigma_calc(h, drhodz, cij)
        std = _std_integral(h, drhodz, cij)
        if refinement > 0 and abs(sigma - sigma_old) < OUTER_TOLERANCE * sigma:
            break
        sigma_old = sigma

        n_half = math.ceil(6.0 * std / h)
        n_ideal = 2 * n_half + 1
        new_half_width = n_half * h
        new_ngrid = 2 * n_ideal - 1
        old_grid = _linspace(-half_width, half_width, ngrid)
        new_grid = _linspace(-new_half_width, new_half_width, new_ngrid)
        next_rows: list[list[float]] = []
        kk = 0
        for x in new_grid:
            if x < old_grid[0]:
                next_rows.append(list(rho_ph1))
            elif x > old_grid[ngrid - 1]:
                next_rows.append(list(rho_ph2))
            else:
                while kk + 1 < ngrid and x > old_grid[kk + 1]:
                    kk += 1
                alpha = (old_grid[kk + 1] - x) / (old_grid[kk + 1] - old_grid[kk])
                next_rows.append(
                    [
                        alpha * rhomat[kk][j] + (1.0 - alpha) * rhomat[kk + 1][j]
                        for j in range(count)
                    ]
                )
        next_rows[0] = list(rho_ph1)
        next_rows[-1] = list(rho_ph2)
        half_width = new_half_width
        ngrid = new_ngrid
        rhomat = next_rows

    return GtSurfaceTensionResult(sigma=from_si(abs(sigma), "N/m"), warnings=tuple(warnings))
