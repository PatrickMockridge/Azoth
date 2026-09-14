"""``eos.critical_point`` - the mixture critical point.

Spec: ``specs/models/eos/critical_point.yaml``

# What a critical point is, and why it is not where the spinodal closes

For a pure component the critical point is where the cubic's three roots merge, and
``dP/dV = d2P/dV2 = 0`` finds it exactly. That route does **not** generalise, and the
way it fails is worth stating because it is silent: in reduced variables those two
conditions depend only on ``(a_tilde, v)`` and have a *single universal root*, so they
predict the same ``Z_c = 0.3074`` for every mixture - exact for one component, wrong for
all the rest, and wrong in a way its own name conceals. This model therefore implements
the two conditions that do generalise, from Heidemann & Khalil (1980).

# The two conditions

At a fixed composition, with ``Q`` the scaled Helmholtz Hessian of
:func:`azoth.eos.reference._mixture_state.criticality_matrix`:

1. the **algebraically smallest** eigenvalue of ``Q`` vanishes;
2. the **cubic form** - the third derivative of the Helmholtz energy along the
   eigenvector of that eigenvalue - vanishes.

Both are evaluated at constant temperature and volume, so the state is ``(T, V)`` and
this model works in ``(T, V)`` throughout rather than in ``(T, P)``. That is not a
detail: at a critical point the cubic is degenerate and its root is ill-conditioned,
and solving for the pressure *from the volume* through the explicit equation of state
avoids ever asking the cubic a question it cannot answer there.

# Why the smallest eigenvalue and not ``det(Q)``

The determinant is the product of every eigenvalue, so it vanishes when *any* of them
does - including ones whose vanishing is not criticality - and being a product it is
badly scaled for a Newton step. The eigenvalue is the quantity whose vanishing is the
condition. Note also "algebraically smallest" rather than "smallest in magnitude": the
eigenvalue that crosses zero at a critical point is the one that goes from negative to
positive, and below the critical point it is not the one nearest zero.

# The cubic form, and where the ``-1`` comes from

The cubic form is a central difference of the Hessian along the eigenvector, which needs
no third-derivative tensor:

```text
C(u) = sum_ijk u_i u_j u_k d3(A^R/RT)/dn_i dn_j dn_k  -  sum_i u_i**3 / n_i**2
```

**The second sum is not optional.** ``A^R`` is the *residual* Helmholtz energy, and the
ideal part's third derivative at constant volume is ``delta_ijk / n_i**2``. Omitting it
leaves a constant offset that is exactly ``1`` at unit composition - so a pure
component's cubic form reads ``1.000000000`` at its critical point instead of zero, and
the outer Newton converges on that offset rather than on the critical point.

# The iteration

Nested, as the published method describes: an inner Newton on ``T`` drives the
eigenvalue to zero at fixed ``V``, and an outer Newton on ``V`` drives the cubic form to
zero. The eigenvector is re-derived at each state rather than carried, because the
outer derivative is taken with the eigenvector that belongs to each perturbed state.
"""

from __future__ import annotations

import math
from typing import Any, NamedTuple

from azoth import _models_gen
from azoth.core.errors import InvalidInputError, OutOfRangeError, SolverNotConvergedError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import CriticalPointResult
from azoth.core.units import from_si
from azoth.core.warnings import Warning
from azoth.eos.mixture import Mixture
from azoth.eos.reference._mixture_state import (
    ReducedParameters,
    criticality_matrix,
    helmholtz_hessian,
    reduced_parameters,
)
from azoth.eos.reference.pr_molar_volume import MOLAR_GAS_CONSTANT

MODEL_ID = "eos.critical_point"

#: The reference pressure the reduced parameters are un-reduced at.
#:
#: `eos.pr_alpha_ab` returns `A_i = a_i P/(R^2 T^2)` and `B_i = b_i P/(R T)`, so the
#: dimensioned pair is a division by any pressure. A reference rather than a
#: re-derivation of the alpha function, which would be a second expression of it.
_REFERENCE_PRESSURE = 1.0e5

#: Guard for the Jacobi sweeps: an off-diagonal below this is treated as zero.
_EIGEN_TINY = 1.0e-300

#: The composition step the cubic form is a central difference over.
#:
#: Fixed rather than adaptive, which puts a floor of about `1e-10` on the outer
#: residual: the difference quotient's error falls as the step squared and its
#: round-off rises as the step squared inverted, and `1e-4` is near the sum's
#: minimum for a third derivative in `f64`. See the spec's notes.
_CUBIC_FORM_STEP = 1.0e-4


class Eigen(NamedTuple):
    """A symmetric matrix's eigenvalues, ascending, and their unit eigenvectors."""

    values: list[float]
    vectors: list[list[float]]


def symmetric_eigen(matrix: list[list[float]], *, tolerance: float, max_iterations: int) -> Eigen:
    """Eigen-decomposition of a symmetric matrix, by cyclic Jacobi rotations.

    Jacobi rather than a tridiagonal reduction because the matrices here are small -
    one row per component - and because it is unconditionally stable and needs no
    pivoting: it converges for any symmetric input, and the rotations are orthogonal
    by construction.

    The eigenvalues come back **ascending**, so the criticality condition can name
    ``values[0]`` rather than searching for it, and ``vectors[i]`` is the unit
    eigenvector belonging to ``values[i]``.

    # Why not the characteristic polynomial

    The cubic solver in this namespace finds roots of a degree-three polynomial in
    closed form. Doing the same here would mean a polynomial of degree ``N``, whose
    coefficients are sums of minors - a computation whose conditioning collapses long
    before the matrix does. The rotations touch only the matrix itself.
    """
    n = len(matrix)
    if n == 0:
        return Eigen(values=[], vectors=[])
    a = [list(row) for row in matrix]
    # `vectors[k][i]` is component `k` of eigenvector `i`; the identity to start.
    vectors = [[1.0 if i == j else 0.0 for j in range(n)] for i in range(n)]

    for _ in range(max_iterations):
        off = math.sqrt(sum(a[p][q] ** 2 for p in range(n) for q in range(n) if p != q))
        if off <= tolerance:
            break
        for p in range(n - 1):
            for q in range(p + 1, n):
                if a[p][q] == 0.0 or abs(a[p][q]) < _EIGEN_TINY:
                    continue
                # The rotation that annihilates `a[p][q]`. `t` is the smaller root of
                # `t^2 + 2 theta t - 1 = 0`, which is the choice that keeps the
                # rotation angle under 45 degrees and the iteration stable.
                theta = (a[q][q] - a[p][p]) / (2.0 * a[p][q])
                sign = 1.0 if theta >= 0.0 else -1.0
                t = sign / (abs(theta) + math.sqrt(theta * theta + 1.0))
                c = 1.0 / math.sqrt(t * t + 1.0)
                s = t * c
                for k in range(n):
                    akp, akq = a[k][p], a[k][q]
                    a[k][p] = c * akp - s * akq
                    a[k][q] = s * akp + c * akq
                for k in range(n):
                    apk, aqk = a[p][k], a[q][k]
                    a[p][k] = c * apk - s * aqk
                    a[q][k] = s * apk + c * aqk
                for k in range(n):
                    vkp, vkq = vectors[k][p], vectors[k][q]
                    vectors[k][p] = c * vkp - s * vkq
                    vectors[k][q] = s * vkp + c * vkq

    values = [a[i][i] for i in range(n)]
    order = sorted(range(n), key=lambda i: values[i])
    return Eigen(
        values=[values[i] for i in order],
        vectors=[[vectors[k][i] for k in range(n)] for i in order],
    )


def dimensional_parameters(mixture: Mixture, temperature: float) -> tuple[list[float], list[float]]:
    """``(a_i, b_i)`` in SI, un-reduced from ``eos.pr_alpha_ab``.

    The reduced pair is ``A_i = a_i P/(R^2 T^2)`` and ``B_i = b_i P/(R T)``, so any
    pressure inverts the reduction. Going through the registered calc rather than
    writing the alpha function a second time is deliberate: there is one expression of
    ``alpha``, and this is not it.
    """
    reduced = reduced_parameters(mixture, temperature, _REFERENCE_PRESSURE)
    scale_a = MOLAR_GAS_CONSTANT**2 * temperature**2 / _REFERENCE_PRESSURE
    scale_b = MOLAR_GAS_CONSTANT * temperature / _REFERENCE_PRESSURE
    return [value * scale_a for value in reduced.a], [value * scale_b for value in reduced.b]


def pressure_at_volume(
    mixture: Mixture, temperature: float, volume: float, composition: list[float]
) -> float:
    """The pressure of the explicit equation of state at ``(T, V, z)``.

    ``P = R T/(V - b) - a/(V^2 + 2 b V - b^2)``, evaluated directly. No cubic is solved,
    which is the point: at a critical point the cubic is degenerate, and this is where
    that degeneracy lives.
    """
    a, b = dimensional_parameters(mixture, temperature)
    count = len(composition)
    b_mix = sum(composition[i] * b[i] for i in range(count))
    a_mix = sum(
        composition[i] * composition[j] * (1.0 - mixture.kij[i][j]) * math.sqrt(a[i] * a[j])
        for i in range(count)
        for j in range(count)
    )
    repulsion = MOLAR_GAS_CONSTANT * temperature / (volume - b_mix)
    attraction = a_mix / (volume * volume + 2.0 * b_mix * volume - b_mix * b_mix)
    return repulsion - attraction


def _cubic_form(
    mixture: Mixture,
    reduced: ReducedParameters,
    compressibility: float,
    composition: list[float],
    vector: list[float],
    step: float,
) -> float:
    """The third derivative of the total ``A/RT`` along ``vector``.

    A central difference of the Hessian along the direction, which is the identity

    ```text
    u^T H(n + s u) u - u^T H(n - s u) u  =  2 s * (third derivative along u) + O(s^3)
    ```

    evaluated at unit total moles, minus the ideal part's third derivative
    ``sum_i u_i**3 / n_i**2``. See the module docstring for why that term is there.
    """
    count = len(composition)
    values = []
    for sign in (1.0, -1.0):
        shifted = [composition[i] + sign * step * vector[i] for i in range(count)]
        hessian = helmholtz_hessian(reduced, mixture.kij, shifted, compressibility)
        values.append(
            sum(vector[i] * hessian[i][j] * vector[j] for i in range(count) for j in range(count))
        )
    residual = (values[0] - values[1]) / (2.0 * step)
    ideal = sum(vector[i] ** 3 / composition[i] ** 2 for i in range(count))
    return residual - ideal


class CriticalPointState(NamedTuple):
    """The two quantities the conditions are written on, at one ``(T, V)``."""

    eigenvalue: float
    vector: list[float]
    compressibility: float
    reduced: ReducedParameters


def _state(
    mixture: Mixture, temperature: float, volume: float, composition: list[float]
) -> CriticalPointState:
    """``Q``'s smallest eigenpair at ``(T, V)``, and the state it sits in."""
    pressure = pressure_at_volume(mixture, temperature, volume, composition)
    compressibility = pressure * volume / (MOLAR_GAS_CONSTANT * temperature)
    reduced = reduced_parameters(mixture, temperature, pressure)
    matrix = criticality_matrix(reduced, mixture.kij, composition, compressibility)
    eigen = symmetric_eigen(matrix, tolerance=1.0e-15, max_iterations=100)
    return CriticalPointState(
        eigenvalue=eigen.values[0],
        vector=eigen.vectors[0],
        compressibility=compressibility,
        reduced=reduced,
    )


def solve_critical_point(
    mixture: Mixture,
    composition: list[float],
    algorithm: Any,
) -> tuple[float, float, float, float, int, float]:
    """``(Tc, Pc, Vc, Z_c, iterations, residual)`` for a composition.

    Nested Newton from a Kay-rule temperature and four times the co-volume, which are
    the initial guesses the method is normally stated with: ``Tc0 = sum z_i Tc_i`` and
    ``V0 = 4 b``, and for propane ``4b`` is ``2.25e-4`` against a true ``Vc`` of
    ``2.2251e-4``.

    Tolerances come from the spec's ``algorithm`` block rather than from the caller,
    so the two implementations run the same loop; see
    :func:`azoth.eos.reference.critical_point`.
    """
    temperature_tolerance = algorithm["inner"]["tolerance"]
    inner_iterations = algorithm["inner"]["max_iterations"]
    residual_tolerance = algorithm["tolerance"]
    max_iterations = algorithm["max_iterations"]
    count = len(composition)
    temperature = sum(
        composition[i] * mixture.components[i].Tc.to_base_units().magnitude for i in range(count)
    )
    # Four times the co-volume, which is the starting volume the method is normally
    # stated with. `b` is a per-component constant - `omega_b R Tc/Pc`, with no
    # temperature in it - so evaluating it at the starting temperature is a
    # convenience rather than a dependence.
    _, b = dimensional_parameters(mixture, temperature)
    volume = 4.0 * sum(composition[i] * b[i] for i in range(count))

    iterations = 0
    residual = float("inf")
    while iterations < max_iterations:
        iterations += 1
        # Inner: Newton on T drives the smallest eigenvalue to zero at this volume.
        # It stops on the eigenvalue rather than on the step, so `inner.tolerance` is a
        # bound on the quantity being solved for and not on how far the last step went.
        for _ in range(inner_iterations):
            state = _state(mixture, temperature, volume, composition)
            if abs(state.eigenvalue) < temperature_tolerance:
                break
            probe = temperature * 1.0e-6
            ahead = _state(mixture, temperature + probe, volume, composition)
            derivative = (ahead.eigenvalue - state.eigenvalue) / probe
            if derivative == 0.0 or not math.isfinite(derivative):
                break
            # The step is capped at a tenth of the temperature, so a Newton step taken
            # from a poor starting point cannot leave the region where the state is
            # defined - where the equation of state has no critical point to find.
            step = max(-0.1 * temperature, min(0.1 * temperature, -state.eigenvalue / derivative))
            temperature += step

        state = _state(mixture, temperature, volume, composition)
        form = _cubic_form(
            mixture,
            state.reduced,
            state.compressibility,
            composition,
            state.vector,
            _CUBIC_FORM_STEP,
        )
        residual = max(abs(state.eigenvalue), abs(form))
        if residual < residual_tolerance:
            break

        probe = volume * 1.0e-6
        ahead = _state(mixture, temperature, volume + probe, composition)
        ahead_form = _cubic_form(
            mixture,
            ahead.reduced,
            ahead.compressibility,
            composition,
            ahead.vector,
            _CUBIC_FORM_STEP,
        )
        derivative = (ahead_form - form) / probe
        if derivative == 0.0 or not math.isfinite(derivative):
            break
        volume_step = max(-0.1 * volume, min(0.1 * volume, -form / derivative))
        # Damped by a half: the cubic form is a finite difference of a Hessian, and a
        # full step overshoots on the first passes from the coarse starting volume.
        volume += 0.5 * volume_step

    if residual >= residual_tolerance:
        raise SolverNotConvergedError(iterations, residual, residual_tolerance)

    pressure = pressure_at_volume(mixture, temperature, volume, composition)
    compressibility = pressure * volume / (MOLAR_GAS_CONSTANT * temperature)
    return temperature, pressure, volume, compressibility, iterations, residual


def critical_point(mixture: Mixture, z: list[float]) -> CriticalPointResult:
    """The critical point of a mixture of composition ``z``.

    ``z`` is the composition of the single phase whose critical point is wanted, not a
    feed being split - so unlike the flash's, it is not something the model decides the
    fate of. Whether that composition could exist at the returned state is a different
    question and is not asked.

    Raises:
        InvalidInputError: if ``z`` is the wrong length, has a negative entry, or does
            not sum to one. Renormalising it here would make a caller's error invisible
            in every number downstream, so it is refused instead.
        SolverNotConvergedError: if the outer iteration hits its cap. The two
            conditions are solved to ``algorithm.tolerance``, and a state that does not
            reach it is reported rather than returned - a critical point that is nearly
            one is not one.
        OutOfRangeError: if the returned state is not a state. A diverged iteration can
            converge on a negative temperature or a volume inside the co-volume, and
            those are refused rather than reported.

    See :func:`azoth.eos.reference.critical_point`.
    """
    spec = _models_gen.model(MODEL_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    composition = list(z)
    if len(composition) != len(mixture.components):
        raise InvalidInputError(
            "z",
            f"a mixture of {len(mixture.components)} components needs that many mole "
            f"fractions, but z has {len(composition)}",
        )
    if any(value < 0.0 for value in composition):
        raise InvalidInputError("z", "a mole fraction cannot be negative")
    total = sum(composition)
    if abs(total - 1.0) > 1.0e-09:
        raise InvalidInputError(
            "z",
            f"the composition sums to {total}, not to one. Renormalising it here would "
            f"make a caller's error invisible in every number downstream, so it is "
            f"refused instead",
        )

    temperature, pressure, volume, compressibility, iterations, residual = solve_critical_point(
        mixture, composition, spec["algorithm"]
    )

    # The bounds are on the *result*, which is unusual for this tree: every input is a
    # vector or a matrix, so there is no scalar argument whose range a caller could
    # violate. What can go wrong is the iteration, and these are how a diverged one
    # announces itself.
    computed = {"tc": temperature, "pc": pressure, "z_c": compressibility}
    apply_checks(checks.derived, computed.get, warnings)
    for field, value in (("tc", temperature), ("pc", pressure), ("z_c", compressibility)):
        if not value > 0.0:
            raise OutOfRangeError(
                field,
                value,
                "the critical point this composition was found to have is not a state; "
                "the iteration left the region where the equation of state is defined",
            )

    return CriticalPointResult(
        tc=from_si(temperature, "K"),
        pc=from_si(pressure, "Pa"),
        vc=from_si(volume, "m**3/mol"),
        z_c=compressibility,
        iterations=iterations,
        residual=residual,
        warnings=tuple(warnings),
    )
