"""``eos.tp_solid_flash`` - how much of a feed has frozen out as one pure solid.

```text
seed     the fluid-only flash, plus a solid of the component the caller named
solve    Q(beta) = sum_k beta_k - sum_i z_i ln E_i      E_i = sum_k beta_k / phi_ik
         with   E_solid = z_solid / phi_solid
drop     every phase the solve pinned under `1.01e-12`
```

Spec: ``specs/models/eos/tp_solid_flash.toml``, which carries the seeding, the pinned
component's rule and the states it is checked at.

NeqSim's ``SolidFlash``, reached through ``Flash.solidPhaseFlash()`` from a ``TPflash`` whose
system has ``setSolidPhaseCheck(name)`` on. A ``TPSolidflash()`` reaches ``SolidFlash1``
instead, which is **one equilibrium and two arithmetics**: this follows ``SolidFlash``,
because that is the class this port's other member is the oracle for and because
``SolidFlash1`` leaves its phases unnormalised.

**The precipitating component's fugacity is pinned by the solid**: ``E_solid`` is
``z_solid/phi_solid`` rather than the sum over phases, which makes ``x_solid^k = phi_solid/phi_ik``
in every fluid phase, and the Hessian loses that component's term because a pinned ``E`` does
not move with ``beta``.

The Rust kernel carries the reasoning; this is the same arithmetic.
"""

from __future__ import annotations

import math
from typing import Any, NamedTuple

from azoth import _models_gen
from azoth.core.errors import InvalidInputError, SolverNotConvergedError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import TpSolidFlashResult
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning
from azoth.eos import components as databank
from azoth.eos.mixture import Mixture
from azoth.eos.reference._mixture_state import (
    ReducedParameters,
    phase_state,
    reduced_parameters,
)
from azoth.eos.reference.pt_flash import pt_flash
from azoth.eos.reference.solid_fugacity import R, solid_fugacity

MODEL_ID = "eos.tp_solid_flash"

#: The floor under a fluid phase fraction, NeqSim's ``phaseFractionMinimumLimit``.
FRACTION_FLOOR = 1.0e-12

#: A fluid phase below this is removed and the solve restarted, upstream's ``1.01e-9``.
REMOVE_LIMIT = 1.01e-9

#: A phase below this is dropped from the answer, upstream's own post-flash sweep.
FRACTION_DROP = 1.01e-12

#: The diagonal regulariser, upstream's ``Qmatrix[i][i] = 1.0e-9``.
REGULARISER = 1.0e-9

#: The restarts a removed fluid phase costs.
MAX_RESTARTS = 8

#: The outer loop's cap and its floor, upstream's ``!(iter > 20)`` and ``iter < 4``.
MAX_OUTER = 20
MINIMUM_OUTER = 4

#: The outer loop's stopping rule on the solid's fraction, upstream's ``1e-3``.
SOLID_TOLERANCE = 1.0e-3

#: The heat-capacity difference NeqSim overrides for water, in J/(mol K).
WATER_DELTA_CP = 37.12

#: The factor NeqSim's ``getDensity()`` carries, which its own fallback divides back out.
DENSITY_SCALE = 1.0e5


class _FluidPhase(NamedTuple):
    """One fluid phase of the solve: an amount, a composition and a root."""

    fraction: float
    composition: list[float]
    liquid: bool


def _polyval(coefficients: tuple[float, ...], temperature: float) -> float:
    """``c1 + c2 T + c3 T^2 + ...``, the form both density correlations are written in."""
    total = 0.0
    for value in reversed(coefficients):
        total = total * temperature + value
    return total


def _molar_volume(
    component: Any, coefficients: tuple[float, ...], temperature: float
) -> float | None:
    """``M / (M 1000 poly(T))``, or None where the table states no correlation."""
    if component.molar_mass is None:
        return None
    molar_mass = float(component.molar_mass.to_base_units().magnitude)
    density = molar_mass * 1000.0 * _polyval(coefficients, temperature)
    if density > 1.0e-20:
        return molar_mass / density
    return None


def _tabulated_solid_fugacity(
    mixture: Mixture, index: int, solid: str, t_si: float, p_si: float, T: Q, P: Q, eos: str
) -> float:
    """A pure solid's fugacity coefficient from the named component's own tables."""
    component = mixture.components[index]
    triple_point = component.triple_point_temperature
    if component.heat_of_fusion <= 0.0 or triple_point <= 0.0:
        raise InvalidInputError(
            "solid",
            f"`{solid}` carries no heat of fusion or no triple-point temperature, and the "
            f"solid route reads both rather than defaulting: a zero heat of fusion is a "
            f"substance that does not melt",
        )
    if component.molar_mass is None:
        raise InvalidInputError(
            "solid",
            f"`{solid}` carries no molar mass, and both density correlations scale by it",
        )

    if solid.strip().lower() == "water":
        delta_cp_sl = WATER_DELTA_CP
    else:
        delta_cp_sl = (
            _polyval(component.cp_liquid, triple_point) - _polyval(component.cp_solid, triple_point)
        ) / 1000.0

    liquid = _molar_volume(component, component.liquid_density_coefs, triple_point)
    if liquid is None:
        # The fallback: the reference phase's own molar volume, read the way NeqSim's
        # `getDensity()` is.
        reference, _ = databank.mixture_of([solid], eos=eos)
        reference_reduced = reduced_parameters(reference, t_si, p_si)
        state = phase_state(reference_reduced, reference.kij, [1.0], liquid=True)
        liquid = state.z * R * t_si / p_si / DENSITY_SCALE
    solid_volume = _molar_volume(component, component.solid_density_coefs, triple_point)
    if solid_volume is None:
        solid_volume = liquid

    return solid_fugacity(
        heat_of_fusion=from_si(component.heat_of_fusion, "J/mol"),
        triple_point_temperature=from_si(triple_point, "K"),
        delta_cp_sl=from_si(delta_cp_sl, "J/(mol*K)"),
        delta_solid_volume=from_si(solid_volume - liquid, "m**3/mol"),
        tc=component.Tc,
        pc=component.Pc,
        omega=component.omega,
        T=T,
        P=P,
        eos=eos,
    ).fugacity_coefficient


def _coefficients(
    reduced: ReducedParameters, kij: Any, phases: list[_FluidPhase]
) -> list[list[float]]:
    """The fugacity coefficients of every component in every fluid phase."""
    return [
        [
            math.exp(value)
            for value in phase_state(reduced, kij, phase.composition, liquid=phase.liquid).ln_phi
        ]
        for phase in phases
    ]


def _denominators(
    z: list[float], phi: list[list[float]], phases: list[_FluidPhase], solid: int, phi_solid: float
) -> list[float]:
    """``E_i`` over the fluid phases, with the precipitating component's value pinned."""
    n = len(z)
    e = [0.0] * n
    for k, phase in enumerate(phases):
        for i in range(n):
            e[i] += phase.fraction / phi[k][i]
    e[solid] = z[solid] / phi_solid
    return e


def _solve(h: list[list[float]], g: list[float]) -> list[float] | None:
    """``H x = g`` by Gaussian elimination with partial pivoting."""
    n = len(g)
    a = [[*h[i][:], g[i]] for i in range(n)]
    for column in range(n):
        pivot = max(range(column, n), key=lambda r: abs(a[r][column]))
        a[column], a[pivot] = a[pivot], a[column]
        if a[column][column] == 0.0:
            return None
        for row in range(column + 1, n):
            factor = a[row][column] / a[column][column]
            for c in range(column, n + 1):
                a[row][c] -= factor * a[column][c]
    x = [0.0] * n
    for row in range(n - 1, -1, -1):
        x[row] = (a[row][n] - sum(a[row][c] * x[c] for c in range(row + 1, n))) / a[row][row]
    return x


def tp_solid_flash(
    components: list[str], solid: str, T: Q, P: Q, z: list[float], eos: str = "srk"
) -> TpSolidFlashResult:
    """The fraction of a feed that has frozen out as one pure solid.

    Args:
        components: the substances, by name.
        solid: the one component allowed to precipitate.
        T: absolute temperature.
        P: absolute pressure.
        z: the overall mole fractions. Checked rather than renormalised.
        eos: the cubic the fluid runs, which is also the one the solid's reference liquid is
            built from.

    Returns:
        The solid's share, the split and the solid's own fugacity coefficient.

    Raises:
        InvalidInputError: if ``z`` is the wrong length, negative or does not sum to one; if
            ``solid`` is not one of ``components``; or if the named component carries no melt
            data.
        OutOfRangeError: if ``T`` or ``P`` is not positive.
        SolverNotConvergedError: if the fraction solve's Hessian is singular.
    """
    spec = _models_gen.model(MODEL_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []
    t_si = input_to_si(spec, "T", T)
    p_si = input_to_si(spec, "P", P)
    apply_checks(checks.on_input, {"T": t_si, "P": p_si}.get, warnings)

    n = len(components)
    if len(z) != n:
        raise InvalidInputError("z", f"a feed for {n} components has {len(z)} entries")
    for i, value in enumerate(z):
        if value < 0.0:
            raise InvalidInputError(
                "z", f"z[{i}] is {value} but a mole fraction cannot be negative"
            )
    total = sum(z)
    if abs(total - 1.0) > 1.0e-09:
        raise InvalidInputError(
            "z",
            f"the mole fractions sum to {total}, not to one. Renormalising them here would "
            f"make a composition error invisible in every number downstream, so it is "
            f"refused instead",
        )

    mixture, _ = databank.mixture_of(components, eos=eos)
    try:
        index = [name.strip().lower() for name in components].index(solid.strip().lower())
    except ValueError:
        raise InvalidInputError("solid", f"`{solid}` is not one of the feed's components") from None

    algorithm = spec["algorithm"]
    tolerance = float(algorithm["tolerance"])
    cap = int(algorithm["max_iterations"])
    reduced = reduced_parameters(mixture, t_si, p_si)
    phi_solid = _tabulated_solid_fugacity(mixture, index, solid, t_si, p_si, T, P, eos)

    flash = pt_flash(mixture, T, P, z)
    if flash.beta is not None:
        phases = [
            _FluidPhase(1.0 - flash.beta, list(flash.x), True),
            _FluidPhase(flash.beta, list(flash.y), False),
        ]
    else:
        liquid = flash.phase.name == "ALL_LIQUID"
        phases = [_FluidPhase(1.0, list(z), liquid)]

    def fluid_only(iterations: int, residual: float, converged: bool) -> TpSolidFlashResult:
        split = flash.beta if flash.beta is not None else 0.0
        return TpSolidFlashResult(
            solid_fraction=0.0,
            phase_count=2,
            beta=(1.0 - split, split),
            x=(tuple(flash.x), tuple(flash.y)),
            solid_fugacity_coefficient=phi_solid,
            iterations=iterations,
            residual=residual,
            converged=converged,
            warnings=tuple(warnings),
        )

    # The candidate screen: how much of the solid component the fluid phases do not account
    # for at the solid's own fugacity.
    start = _coefficients(reduced, mixture.kij, phases)
    candidate = z[index]
    for k, phase in enumerate(phases):
        candidate -= phase.fraction * phi_solid / start[k][index]
    if candidate <= 1.0e-20:
        return fluid_only(0, 0.0, True)

    solid_fraction = candidate
    iterations = 0
    residual = math.nan
    converged = False
    restarts = 0

    def solve_beta() -> tuple[int, float, bool, bool]:
        nonlocal solid_fraction
        steps = 0
        left = math.nan
        for step in range(1, cap + 1):
            steps = step
            phi = _coefficients(reduced, mixture.kij, phases)
            e = _denominators(z, phi, phases, index, phi_solid)

            count = len(phases)
            gradient = [1.0] * count
            for k in range(count):
                for i in range(n):
                    gradient[k] -= z[i] / (e[i] * phi[k][i])
            hessian = [[0.0] * count for _ in range(count)]
            for j in range(count):
                for k in range(count):
                    hessian[j][k] = REGULARISER if j == k else 0.0
                    for i in range(n):
                        if i != index:
                            hessian[j][k] += z[i] / (e[i] * e[i] * phi[j][i] * phi[k][i])

            correction = _solve(hessian, gradient)
            if correction is None:
                raise SolverNotConvergedError(steps, math.nan, tolerance)
            left = math.sqrt(sum(value * value for value in correction))

            damping = (step + 1.0) / (10.0 + step)
            for k, phase in enumerate(phases):
                value = phase.fraction - damping * correction[k]
                if value < FRACTION_FLOOR:
                    value = FRACTION_FLOOR
                elif value > 1.0:
                    value = 1.0 - FRACTION_FLOOR
                phases[k] = _FluidPhase(value, phase.composition, phase.liquid)

            phi = _coefficients(reduced, mixture.kij, phases)
            solid_fraction = z[index] - sum(
                phase.fraction * phi_solid / phi[k][index] for k, phase in enumerate(phases)
            )

            for k, phase in enumerate(phases):
                if abs(phase.fraction) < REMOVE_LIMIT:
                    gone = phases.pop(k)
                    kept = 1.0 - gone.fraction
                    for j, other in enumerate(phases):
                        phases[j] = _FluidPhase(
                            other.fraction / kept, other.composition, other.liquid
                        )
                    return steps, left, False, True

            if step >= 2 and left <= tolerance:
                return steps, left, True, False
        return steps, left, False, False

    def set_compositions() -> None:
        phi = _coefficients(reduced, mixture.kij, phases)
        e = _denominators(z, phi, phases, index, phi_solid)
        for k, phase in enumerate(phases):
            row = [z[i] / (e[i] * phi[k][i]) for i in range(n)]
            row_total = sum(row)
            if row_total > 0.0:
                row = [value / row_total for value in row]
            phases[k] = _FluidPhase(phase.fraction, row, phase.liquid)

    # `run()`, restarted when a fluid phase is removed - the same `continue 'restart` the
    # Rust kernel's labelled loop performs.
    while True:
        restart = False
        for outer in range(1, MAX_OUTER + 1):
            previous = solid_fraction
            iterations, residual, converged, removed = solve_beta()
            if removed:
                restart = restarts < MAX_RESTARTS
                restarts += 1
                break
            set_compositions()
            phi = _coefficients(reduced, mixture.kij, phases)
            solid_fraction = z[index] - sum(
                phase.fraction * phi_solid / phi[k][index] for k, phase in enumerate(phases)
            )
            if outer >= MINIMUM_OUTER and abs(solid_fraction - previous) <= SOLID_TOLERANCE:
                break
        if not restart:
            break

    beta = [phase.fraction for phase in phases] + [solid_fraction]
    x = [list(phase.composition) for phase in phases]
    row = [0.0] * n
    row[index] = 1.0
    x.append(row)
    keep = [k for k in range(len(beta)) if beta[k] >= FRACTION_DROP]
    solid_present = keep and keep[-1] == len(phases)

    return TpSolidFlashResult(
        solid_fraction=solid_fraction if solid_present else 0.0,
        phase_count=len(keep),
        beta=tuple(beta[k] for k in keep),
        x=tuple(tuple(x[k]) for k in keep),
        solid_fugacity_coefficient=phi_solid,
        iterations=iterations,
        residual=residual,
        converged=converged,
        warnings=tuple(warnings),
    )


__all__ = ["tp_solid_flash"]
