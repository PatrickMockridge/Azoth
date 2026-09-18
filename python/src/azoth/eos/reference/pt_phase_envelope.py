"""``eos.pt_phase_envelope`` - the Michelsen natural-parameter phase envelope.

Spec: ``specs/models/eos/pt_phase_envelope.toml``

The envelope is traced as two independent natural-parameter continuations - the
bubble branch at a tiny vapour fraction and the dew branch at one minus it - each
bootstrapped at a low pressure where its K-values are far from one, then continued
upward to the critical point. The ``N+2`` unknowns ``(ln K_i, ln T, ln P)`` satisfy
``N`` isofugacity equations, one material balance, and one specification equation.
"""

from __future__ import annotations

import math

from azoth import _models_gen
from azoth.core.errors import InvalidInputError, SolverNotConvergedError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import PtPhaseEnvelopeResult
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning
from azoth.eos.mixture import Mixture
from azoth.eos.reference._mixture_state import (
    WILSON_CONSTANT,
    phase_state,
    reduced_parameters,
    wilson_saturation_pressures,
)

MODEL_ID = "eos.pt_phase_envelope"

BUBBLE_BETA = 1.0e-6
DEW_BETA = 1.0 - 1.0e-6
D_TMAX = 10.0
D_PMAX = 10.0e5
MAX_PRESSURE = 1000.0e5
NEWTON_TOL = 1.0e-5
NEWTON_MAX = 50

#: The most Newton steps the critical refinement takes, upstream's ten.
#: Upstream's ten, and it is a guard rather than a cap: measured, fifty steps let the Newton
#: wander to ``355.78 K`` with a branch gap of ``21.0 K``.
CRIT_MAX = 10
#: Its convergence tolerance, upstream's. Unreachable on this system in double precision, so
#: the refinement stops at its best iterate instead of converging; kept at upstream's value
#: because lowering it changed nothing - measured, `1e-7` returned the same points.
CRIT_TOL = 1.0e-10
#: The largest single correction it accepts, upstream's ``|dx_j| > 0.5`` abort.
CRIT_STEP = 0.5
#: How far a candidate critical point may sit from where the refinement started, relative.
#: Upstream compares against its polynomial's extrapolation and keeps a candidate only within
#: ``0.10`` in temperature and ``0.20`` in pressure of it; the start is used here because the
#: refinement begins where the trace crossed the K-window and the two are within a few per cent
#: of each other there. Without the guard a diverging Newton's wild iterate wins the
#: best-tracking on ``sum (ln K)^2`` alone, which it otherwise does - measured, a step that
#: jumped 373 K to 437 K had the smallest sum of any iterate.
CRIT_T_BAND = 0.10
CRIT_P_BAND = 0.20
PRESSURE_SPEC_POINTS = 5
HISTORY_LEN = 4


def _solve_linear(a: list[list[float]], b: list[float]) -> list[float] | None:
    """``A x = b`` by Gaussian elimination with partial pivoting, in place."""
    n = len(b)
    for col in range(n):
        pivot = col
        for row in range(col + 1, n):
            if abs(a[row][col]) > abs(a[pivot][col]):
                pivot = row
        if abs(a[pivot][col]) < 1.0e-15:
            return None
        a[col], a[pivot] = a[pivot], a[col]
        b[col], b[pivot] = b[pivot], b[col]
        for row in range(col + 1, n):
            factor = a[row][col] / a[col][col]
            for j in range(col, n):
                a[row][j] -= factor * a[col][j]
            b[row] -= factor * b[col]
    x = [0.0] * n
    for row in range(n - 1, -1, -1):
        s = b[row]
        for j in range(row + 1, n):
            s -= a[row][j] * x[j]
        x[row] = s / a[row][row]
    return x


def _compositions(z: list[float], k: list[float], beta: float) -> tuple[list[float], list[float]]:
    n = len(z)
    x = [0.0] * n
    y = [0.0] * n
    for i in range(n):
        d = 1.0 - beta + beta * k[i]
        x[i] = z[i] / d
        y[i] = k[i] * x[i]
    return x, y


def _residual(
    mixture: Mixture,
    u: list[float],
    beta: float,
    z: list[float],
    last: tuple[int, float] | None,
) -> list[float]:
    """The `N+2` residuals.

    ``last`` is ``(index, value)`` for an ordinary continuation point, pinning ``u[index]``, or
    ``None`` for the criticality condition ``sum_i (ln K_i)^2``, whose Jacobian row the central
    difference below produces without being told.
    """
    n = len(z)
    k = [math.exp(lnk) for lnk in u[:n]]
    t = math.exp(u[n])
    p = math.exp(u[n + 1])
    x, y = _compositions(z, k, beta)
    reduced = reduced_parameters(mixture, t, p)
    liquid = phase_state(reduced, mixture.kij, list(x), liquid=True)
    vapour = phase_state(reduced, mixture.kij, list(y), liquid=False)
    f = [0.0] * (n + 2)
    for i in range(n):
        f[i] = u[i] - liquid.ln_phi[i] + vapour.ln_phi[i]
    f[n] = sum(y) - sum(x)
    if last is None:
        f[n + 1] = sum(value * value for value in u[:n])
    else:
        f[n + 1] = u[last[0]] - last[1]
    return f


def _jacobian(
    mixture: Mixture,
    u: list[float],
    beta: float,
    z: list[float],
    last: tuple[int, float] | None,
) -> list[list[float]]:
    m = len(u)
    f0 = _residual(mixture, u, beta, z, last)
    jac = [[0.0] * m for _ in range(m)]
    for j in range(m):
        delta = max(1.0e-6 * abs(u[j]), 1.0e-8)
        up = u.copy()
        up[j] += delta
        fp = _residual(mixture, up, beta, z, last)
        for i in range(m):
            jac[i][j] = (fp[i] - f0[i]) / delta
    return jac


def _norm2(v: list[float]) -> float:
    return math.sqrt(sum(x * x for x in v))


def _clamp_state(u: list[float], n: int) -> None:
    for i in range(n):
        if not math.isfinite(u[i]):
            u[i] = 0.0
        u[i] = max(-20.0, min(20.0, u[i]))
    if not math.isfinite(u[n]):
        u[n] = math.log(300.0)
    if not math.isfinite(u[n + 1]):
        u[n + 1] = math.log(1.0e5)
    t = max(10.0, min(2000.0, math.exp(u[n])))
    p = max(1.0, min(MAX_PRESSURE * 1.5, math.exp(u[n + 1])))
    u[n] = math.log(t)
    u[n + 1] = math.log(p)


def _newton(
    mixture: Mixture,
    u: list[float],
    beta: float,
    z: list[float],
    last: tuple[int, float] | None,
) -> tuple[list[float], int]:
    u = u.copy()
    for iters in range(1, NEWTON_MAX + 1):
        f = _residual(mixture, u, beta, z, last)
        norm = _norm2(f)
        if norm < NEWTON_TOL:
            return u, iters
        jac = _jacobian(mixture, u, beta, z, last)
        b = f.copy()
        dx = _solve_linear(jac, b)
        if dx is None or any(not math.isfinite(d) for d in dx):
            raise SolverNotConvergedError(iters, norm, NEWTON_TOL)
        step = 1.0
        nxt = u.copy()
        for _ in range(24):
            for i, ui in enumerate(u):
                nxt[i] = ui - step * dx[i]
            _clamp_state(nxt, len(z))
            fnext = _residual(mixture, nxt, beta, z, last)
            if _norm2(fnext) < norm:
                break
            step *= 0.5
        u = nxt
    raise SolverNotConvergedError(
        NEWTON_MAX, _norm2(_residual(mixture, u, beta, z, last)), NEWTON_TOL
    )


def _calc_crit(
    mixture: Mixture, u0: list[float], beta: float, z: list[float]
) -> tuple[float, float] | None:
    """The critical point, refined from a state the trace has brought close to it.

    **This is what replaces a test with a computation.** The trace stops its branches where the
    lightest component's K-value falls below ``1.05`` and the heaviest's rises above ``0.95``,
    which is a place a heuristic happened to fire: it put this library's critical point 6.86 K
    and 4.54 bar from NeqSim's. The refinement Newtons ``sum_i (ln K_i)^2 = 0`` alongside the
    isofugacity rows and the material balance, which is the definition rather than a proxy.

    Returns ``None`` when it does not converge, which the caller reports as the unrefined point:
    the trace's own answer is still a boundary. The **best** iterate is kept as well as the last,
    upstream's rule - the Newton can step past the smallest sum and come back with a larger one.

    The rust kernel carries the reasoning; this is the same arithmetic.
    """
    n = len(z)
    u = list(u0)
    start = (math.exp(u0[n]), math.exp(u0[n + 1]))
    best: tuple[float, float, float] | None = None

    for _ in range(CRIT_MAX):
        f = _residual(mixture, u, beta, z, None)
        sum_ln_k2 = sum(value * value for value in u[:n])
        t, p = math.exp(u[n]), math.exp(u[n + 1])
        inside = (
            abs(t - start[0]) <= CRIT_T_BAND * start[0]
            and abs(p - start[1]) <= CRIT_P_BAND * start[1]
        )
        if inside and math.isfinite(sum_ln_k2) and (best is None or sum_ln_k2 < best[0]):
            best = (sum_ln_k2, t, p)
        if _norm2(f) < CRIT_TOL:
            return (t, p)

        jac = _jacobian(mixture, u, beta, z, None)
        # Levenberg-Marquardt, ``flash_newton.py``'s and at the same magnitude. The system is
        # singular at the critical point, so the Newton is ill-conditioned exactly where it is
        # aimed; without this the iteration stops at a knife-edge iterate and two
        # implementations of the same arithmetic return points 0.035 K apart.
        trace = sum(abs(jac[i][i]) for i in range(n + 2))
        lam = 1.0e-8 * trace / (n + 2)
        for i in range(n + 2):
            jac[i][i] += lam
        b = list(f)
        dx = _solve_linear(jac, b)
        if dx is None or any(not math.isfinite(d) or abs(d) > CRIT_STEP for d in dx):
            break
        if _norm2(dx) < CRIT_TOL:
            return (t, p)

        # The step-halving line search the continuation Newton uses, and it is not optional
        # here: upstream restarts from a polynomial extrapolated to the K = 1 point, which is
        # already inside the critical basin, and this starts from the K-window crossing.
        # Measured without it: 376.696 -> 375.501 -> 403.180, the last outside the band.
        norm = _norm2(f)
        step = 1.0
        accepted = False
        for _ in range(24):
            nxt = list(u)
            for i, d in enumerate(dx):
                nxt[i] = u[i] - step * d
            _clamp_state(nxt, n)
            if _norm2(_residual(mixture, nxt, beta, z, None)) < norm:
                u = nxt
                accepted = True
                break
            step *= 0.5
        if not accepted:
            break

    return None if best is None else (best[1], best[2])


def _sensitivity(
    mixture: Mixture, u: list[float], beta: float, z: list[float], speceq: int
) -> list[float]:
    n = len(z)
    jac = _jacobian(mixture, u, beta, z, (speceq, u[speceq]))
    e = [0.0] * (n + 2)
    e[n + 1] = 1.0
    dxds = _solve_linear(jac, e)
    if dxds is None:
        raise SolverNotConvergedError(0, math.nan, NEWTON_TOL)
    return dxds


def _find_spec(dxds: list[float], n: int) -> int:
    return n if abs(dxds[n]) >= abs(dxds[n + 1]) else n + 1


def _cubic_predict(history: list[list[float]], speceq: int, sny: float) -> list[float]:
    m = len(history[0])
    s = [h[speceq] for h in history]
    vandermonde = [[0.0] * HISTORY_LEN for _ in range(HISTORY_LEN)]
    for i, si in enumerate(s):
        vandermonde[i][0] = 1.0
        vandermonde[i][1] = si
        vandermonde[i][2] = si * si
        vandermonde[i][3] = si * si * si
    result = [0.0] * m
    for j in range(m):
        v = [h[j] for h in history]
        coeffs = _solve_linear([row.copy() for row in vandermonde], v.copy())
        if coeffs is None:
            coeffs = [v[3], 0.0, 0.0, 0.0]
        result[j] = coeffs[0] + sny * (coeffs[1] + sny * (coeffs[2] + sny * coeffs[3]))
    return result


def _light_heavy(mixture: Mixture) -> tuple[int, int]:
    lc = 0
    hc = 0
    for i, c in enumerate(mixture.components):
        if c.Tc.to_base_units().magnitude < mixture.components[lc].Tc.to_base_units().magnitude:
            lc = i
        if c.Tc.to_base_units().magnitude > mixture.components[hc].Tc.to_base_units().magnitude:
            hc = i
    return lc, hc


def _wilson_temperature(mixture: Mixture, pressure: float, z: list[float], beta: float) -> float:
    n = len(z)
    lc, hc = _light_heavy(mixture)
    idx = lc if beta < 0.5 else hc
    c = mixture.components[idx]
    ln_p_ratio = math.log(pressure / c.Pc.to_base_units().magnitude)
    t = (
        c.Tc.to_base_units().magnitude
        * WILSON_CONSTANT
        * (1.0 + c.omega)
        / (WILSON_CONSTANT * (1.0 + c.omega) - ln_p_ratio)
    )
    told = 0.0
    for _ in range(1000):
        psat = wilson_saturation_pressures(mixture.components, t)
        kwil = [ps / pressure for ps in psat]
        s = 0.0
        dsdt = 0.0
        for i in range(n):
            ci = mixture.components[i]
            dlnkdt = WILSON_CONSTANT * (1.0 + ci.omega) * ci.Tc.to_base_units().magnitude / (t * t)
            if beta < 0.5:
                s += z[i] * kwil[i]
                dsdt += z[i] * kwil[i] * dlnkdt
            else:
                s += z[i] / kwil[i]
                dsdt -= z[i] / kwil[i] * dlnkdt
        s -= 1.0
        if abs(s / dsdt) > 0.1 * t:
            t -= 0.001 * s / dsdt
        else:
            t -= s / dsdt
        if abs(t - told) < 1e-5:
            break
        told = t
    return t


def _trace_branch(
    mixture: Mixture, pressure: float, z: list[float], beta: float
) -> tuple[
    list[float], list[float], tuple[float, float], tuple[float, float], tuple[float, float] | None
]:
    spec = _models_gen.model(MODEL_ID)
    max_iterations = spec["algorithm"]["max_iterations"]
    n = len(z)
    lc, hc = _light_heavy(mixture)

    start_t = _wilson_temperature(mixture, pressure, z, beta)
    psat = wilson_saturation_pressures(mixture.components, start_t)
    u = [math.log(ps / pressure) for ps in psat] + [math.log(start_t), math.log(pressure)]

    speceq = n + 1
    specval = u[speceq]
    t_vec: list[float] = []
    p_vec: list[float] = []
    cricondenbar = (start_t, pressure)
    cricondentherm = (start_t, pressure)
    critical: tuple[float, float] | None = None
    history: list[list[float]] = []
    dxds: list[float] = [0.0] * (n + 2)
    ds = 0.1
    last_iters = 2

    u, _ = _newton(mixture, u, beta, z, (speceq, specval))

    for _ in range(max_iterations):
        t = math.exp(u[n])
        pv = math.exp(u[n + 1])
        k = [math.exp(lnk) for lnk in u[:n]]

        if t > cricondentherm[0]:
            cricondentherm = (t, pv)
        if pv > cricondenbar[1]:
            cricondenbar = (t, pv)

        if k[lc] < 1.05 and k[hc] > 0.95:
            critical = _calc_crit(mixture, u, beta, z)
            break

        t_vec.append(t)
        p_vec.append(pv)
        history.append(u.copy())
        if len(history) > HISTORY_LEN:
            history.pop(0)

        if pv > MAX_PRESSURE:
            break

        if len(t_vec) <= PRESSURE_SPEC_POINTS:
            speceq = n + 1
            dxds = _sensitivity(mixture, u, beta, z, speceq)
            ds = 0.1 / dxds[n + 1]
            for i, d in enumerate(dxds):
                u[i] += d * ds
            specval = u[n + 1]
        else:
            speceq = _find_spec(dxds, n)
            sign = 1.0 if dxds[speceq] >= 0.0 else -1.0
            ds = sign * abs(ds)
            dxds = _sensitivity(mixture, u, beta, z, speceq)
            if last_iters > 6:
                ds *= 0.5
            elif last_iters < 3:
                ds *= 1.1
            elif last_iters == 4:
                ds *= 0.9
            elif last_iters > 4:
                ds *= 0.7
            ds = (1.0 if ds >= 0.0 else -1.0) * abs(ds)
            t_cur = math.exp(u[n])
            p_cur = math.exp(u[n + 1])
            if abs(dxds[n]) * abs(ds) > math.log(1.0 + D_TMAX / t_cur):
                ds = (1.0 if ds >= 0.0 else -1.0) * math.log(1.0 + D_TMAX / t_cur) / abs(dxds[n])
            if abs(dxds[n + 1]) * abs(ds) > math.log(1.0 + D_PMAX / p_cur):
                ds = (
                    (1.0 if ds >= 0.0 else -1.0) * math.log(1.0 + D_PMAX / p_cur) / abs(dxds[n + 1])
                )
            if len(history) == HISTORY_LEN:
                sny = ds * dxds[speceq] + history[HISTORY_LEN - 1][speceq]
                u = _cubic_predict(history, speceq, sny)
                specval = sny
            else:
                for i, d in enumerate(dxds):
                    u[i] += d * ds
                specval = u[speceq]

        _clamp_state(u, n)
        predicted = u.copy()
        solved = False
        for _ in range(8):
            try:
                u, iters = _newton(mixture, u, beta, z, (speceq, specval))
                last_iters = iters
                solved = True
                break
            except SolverNotConvergedError:
                ds *= 0.5
                prev = history[-1] if history else predicted
                for i, d in enumerate(dxds):
                    u[i] = prev[i] + d * ds
                specval = u[speceq]
        if not solved:
            break

    return t_vec, p_vec, cricondenbar, cricondentherm, critical


def pt_phase_envelope(mixture: Mixture, P: Q, z: list[float]) -> PtPhaseEnvelopeResult:
    """The PT phase envelope of a mixture of composition ``z``, traced from ``P``.

    Raises:
        InvalidInputError: if the mixture has one component, or if ``z`` is not a
            composition.
        OutOfRangeError: if ``P`` is not positive.
        SolverNotConvergedError: if a branch's continuation hits its cap.
    """
    spec = _models_gen.model(MODEL_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    pressure = input_to_si(spec, "P", P)
    apply_checks(checks.on_input, {"P": pressure}.get, warnings)

    n = len(mixture)
    if n < 2:
        raise InvalidInputError(
            "components",
            "a phase envelope needs two phases with different compositions, and one "
            "component cannot have them. A pure component's envelope is its "
            "vapour-pressure curve, which `eos.pure_saturation` computes.",
        )
    if len(z) != n:
        raise InvalidInputError("z", f"a composition for {n} components has {len(z)} entries")

    bubble = _trace_branch(mixture, pressure, list(z), BUBBLE_BETA)
    dew = _trace_branch(mixture, pressure, list(z), DEW_BETA)

    cricondenbar = bubble[2] if bubble[2][1] >= dew[2][1] else dew[2]
    cricondentherm = bubble[3] if bubble[3][0] >= dew[3][0] else dew[3]

    # The critical point is the refined one when either branch reached it, and either will do:
    # they are the same point and each branch computed it alone, so a disagreement is a finding
    # rather than a choice. Falling back to the branch endpoint keeps a trace that never got
    # near criticality reporting where it stopped rather than a NaN.
    critical = bubble[4] or dew[4] or ((bubble[0][-1], bubble[1][-1]) if bubble[0] else None)
    if critical is None:
        critical = (math.nan, math.nan)
    # The residual is the two branches' disagreement about the critical temperature: zero when
    # both refined to the same point, which is the statement that they meet. A branch that did
    # not refine reports the gap between where the branches stopped - the older, weaker claim.
    if bubble[4] is not None and dew[4] is not None:
        residual = abs(bubble[4][0] - dew[4][0])
    else:
        residual = abs(bubble[0][-1] - dew[0][-1]) if bubble[0] and dew[0] else math.nan

    return PtPhaseEnvelopeResult(
        dew_temperature=tuple(dew[0]),
        dew_pressure=tuple(dew[1]),
        bubble_temperature=tuple(bubble[0]),
        bubble_pressure=tuple(bubble[1]),
        cricondenbar_temperature=from_si(cricondenbar[0], "K"),
        cricondenbar_pressure=from_si(cricondenbar[1], "Pa"),
        cricondentherm_temperature=from_si(cricondentherm[0], "K"),
        cricondentherm_pressure=from_si(cricondentherm[1], "Pa"),
        critical_temperature=from_si(critical[0], "K"),
        critical_pressure=from_si(critical[1], "Pa"),
        iterations=len(bubble[0]) + len(dew[0]),
        residual=residual,
        warnings=tuple(warnings),
    )
