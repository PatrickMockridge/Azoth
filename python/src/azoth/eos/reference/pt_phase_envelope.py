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
    mixture: Mixture, u: list[float], beta: float, z: list[float], speceq: int, specval: float
) -> list[float]:
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
    f[n + 1] = u[speceq] - specval
    return f


def _jacobian(
    mixture: Mixture, u: list[float], beta: float, z: list[float], speceq: int, specval: float
) -> list[list[float]]:
    m = len(u)
    f0 = _residual(mixture, u, beta, z, speceq, specval)
    jac = [[0.0] * m for _ in range(m)]
    for j in range(m):
        delta = max(1.0e-6 * abs(u[j]), 1.0e-8)
        up = u.copy()
        up[j] += delta
        fp = _residual(mixture, up, beta, z, speceq, specval)
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
    mixture: Mixture, u: list[float], beta: float, z: list[float], speceq: int, specval: float
) -> tuple[list[float], int]:
    u = u.copy()
    for iters in range(1, NEWTON_MAX + 1):
        f = _residual(mixture, u, beta, z, speceq, specval)
        norm = _norm2(f)
        if norm < NEWTON_TOL:
            return u, iters
        jac = _jacobian(mixture, u, beta, z, speceq, specval)
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
            fnext = _residual(mixture, nxt, beta, z, speceq, specval)
            if _norm2(fnext) < norm:
                break
            step *= 0.5
        u = nxt
    raise SolverNotConvergedError(
        NEWTON_MAX, _norm2(_residual(mixture, u, beta, z, speceq, specval)), NEWTON_TOL
    )


def _sensitivity(
    mixture: Mixture, u: list[float], beta: float, z: list[float], speceq: int
) -> list[float]:
    n = len(z)
    jac = _jacobian(mixture, u, beta, z, speceq, u[speceq])
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
) -> tuple[list[float], list[float], tuple[float, float], tuple[float, float]]:
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
    history: list[list[float]] = []
    dxds: list[float] = [0.0] * (n + 2)
    ds = 0.1
    last_iters = 2

    u, _ = _newton(mixture, u, beta, z, speceq, specval)

    for _ in range(max_iterations):
        t = math.exp(u[n])
        pv = math.exp(u[n + 1])
        k = [math.exp(lnk) for lnk in u[:n]]

        if t > cricondentherm[0]:
            cricondentherm = (t, pv)
        if pv > cricondenbar[1]:
            cricondenbar = (t, pv)

        if k[lc] < 1.05 and k[hc] > 0.95:
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
                u, iters = _newton(mixture, u, beta, z, speceq, specval)
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

    return t_vec, p_vec, cricondenbar, cricondentherm


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

    critical = (bubble[0][-1], bubble[1][-1]) if bubble[0] else (math.nan, math.nan)
    # The two branches stop at the K-value crossing but not exactly at the same point;
    # the residual reports the temperature gap between their endpoints.
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
