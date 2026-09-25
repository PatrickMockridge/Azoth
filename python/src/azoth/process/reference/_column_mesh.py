"""Naphtali-Sandholm in Python: the twin of ``column/naphtali_sandholm.rs``.

Spec: ``specs/models/process/distillation_column.toml``, ``solver_type``.

The simultaneous MESH correction: one Newton step on the whole column at once - ``N`` tray
blocks of ``C + 2`` variables (the ``C`` liquid component flows, the temperature and the vapour
flow), the ``C`` component balances, the energy balance and ``sum(K x) = 1`` - with a
finite-difference Jacobian, a block-tridiagonal solve, a trust region and an Armijo line search.

**Every constant is `NaphtaliSandholmSolver`'s**, as the Rust twin keeps them, and the paper is
the class's own citation: Naphtali & Sandholm (1971), "Multicomponent Separation Calculations by
Linearization", AIChE Journal 17(1) 148-153.

**The Wilson fallback is the class's own 5.37**, not this library's shared 5.373
(``_mixture_state.WILSON_CONSTANT``, whose own doc records the ambiguity): the class writes its
own constant in `wilsonK`, the port keeps it, and the two differ by `0.06` per cent in the
exponent. It reaches only a tray whose previous K-value is unusable, because every other K comes
from the EOS.

**The flows are mol/hr inside the solve and mol/s outside it.** `flowScale`, the trust region's
floors and the perturbation floor are all absolute, so a solve scaled by one library's flow unit
and stopped by the other's tolerance is a different solve.

**The seed is the class's warm start** (`initializeTrayStateFromColumn`): the substitution
core's converged trays, mapped onto the MESH variables. The class's cold seed is not ported, for
the reason the Rust twin's ``initialize_from_column`` records.
"""

from __future__ import annotations

import math
from typing import TYPE_CHECKING

from azoth.core.errors import InvalidInputError, SolverNotConvergedError
from azoth.core.result import Phase
from azoth.core.units import quantity
from azoth.eos import components as _components
from azoth.eos.reference._mixture_state import phase_state, reduced_parameters
from azoth.eos.reference.molar_enthalpy_entropy import molar_enthalpy_entropy

if TYPE_CHECKING:
    from azoth.process.reference.distillation_column import _States

#: The class's own flow unit, which every absolute constant of this solver is written in.
MOL_PER_HOUR = 3600.0
#: Maximum non-descent line-search steps before the best state is restored.
MAX_NON_DESCENT_LINE_SEARCH_STEPS = 3
#: Forced-root fugacity fixed-point sweeps for one tray evaluation.
THERMO_K_VALUE_ITERATIONS = 2
#: Convergence tolerance for the largest absolute logarithmic K-value update.
THERMO_K_VALUE_TOLERANCE = 1.0e-8
#: Maximum Newton steps.
MAX_ITERATIONS = 80
#: The class's own convergence tolerance, and the one `solveNaphtaliSandholm` sets.
TOLERANCE = 1.0e-8
#: Relative perturbation for the numerical Jacobian, and its absolute floor.
PERTURBATION = 1.0e-4
MIN_PERTURBATION = 1.0e-8
#: The largest per-tray, per-component imbalance an early exit accepts, relative to the feed.
COMPONENT_IMBALANCE_TOLERANCE = 5.0e-3
#: `wilsonK`'s own constant, which is not the library's shared `WILSON_CONSTANT` - see the module.
WILSON_CONSTANT = 5.37


class Mesh:
    """The MESH state, its thermodynamics, and the Newton loop over both."""

    def __init__(self, components: list[str], trays: int, pressures: list[float]) -> None:
        self.components = list(components)
        self.n = trays
        self.c = len(components)
        self.m = self.c + 2
        self.t = [0.0] * trays
        self.p = list(pressures)
        self.v = [0.0] * trays
        self.l = [0.0] * trays
        self.liq = [[0.0] * self.c for _ in range(trays)]
        self.vap = [[0.0] * self.c for _ in range(trays)]
        self.k = [[0.0] * self.c for _ in range(trays)]
        self.hl = [0.0] * trays
        self.hv = [0.0] * trays
        self.feed_liq = [[0.0] * self.c for _ in range(trays)]
        self.feed_vap = [[0.0] * self.c for _ in range(trays)]
        self.feed_hl = 0.0
        self.feed_hv = 0.0
        self.feed_l_total = [0.0] * trays
        self.feed_v_total = [0.0] * trays
        self.internal_vapour_fraction = [1.0] * trays
        self.internal_liquid_fraction = [1.0] * trays
        self.fixed_temperature = [float("nan")] * trays
        #: The residual's scaling, `totalFeed / N` floored at one.
        self.flow_scale = 1.0
        #: The temperature residual's scaling, 100 K.
        self.temp_scale = 100.0
        self.mixture, self.ideal_gas = _components.mixture_of(self.components, eos="pr")

    # ---- the variables ----------------------------------------------------------------

    def variable(self, tray: int, k: int) -> float:
        """One variable, in the order the Jacobian is assembled in."""
        if k < self.c:
            return self.liq[tray][k]
        if k == self.c:
            return self.t[tray]
        return self.v[tray]

    def set_variable(self, tray: int, k: int, value: float) -> None:
        """Set one variable."""
        if k < self.c:
            self.liq[tray][k] = value
        elif k == self.c:
            self.t[tray] = value
        else:
            self.v[tray] = value

    # ---- the equations ----------------------------------------------------------------

    def residual(self) -> list[float]:
        """The residual vector, one block per tray."""
        f: list[float] = []
        for j in range(self.n):
            f.extend(self.residual_for_tray(j))
        return f

    def residual_for_tray(self, j: int) -> list[float]:
        """One tray's `C + 2` residuals: the balances, the energy, and `sum(K x) - 1`."""
        f = [0.0] * self.m
        lj = self.l[j]
        for i in range(self.c):
            balance = self.liq[j][i] + self.vap[j][i]
            if j + 1 < self.n:
                balance -= self.internal_liquid_fraction[j + 1] * self.liq[j + 1][i]
            if j > 0:
                balance -= self.internal_vapour_fraction[j - 1] * self.vap[j - 1][i]
            balance -= self.feed_liq[j][i] + self.feed_vap[j][i]
            f[i] = balance / self.flow_scale

        if self.fixed_temperature[j] == self.fixed_temperature[j]:
            f[self.c] = (self.t[j] - self.fixed_temperature[j]) / self.temp_scale
        else:
            h = lj * self.hl[j] + self.v[j] * self.hv[j]
            if j + 1 < self.n:
                h -= self.internal_liquid_fraction[j + 1] * self.l[j + 1] * self.hl[j + 1]
            if j > 0:
                h -= self.internal_vapour_fraction[j - 1] * self.v[j - 1] * self.hv[j - 1]
            h -= self.feed_l_total[j] * self.feed_hl + self.feed_v_total[j] * self.feed_hv
            tray_flow = max(lj + self.v[j], self.flow_scale)
            latent = max(abs(self.hv[j] - self.hl[j]), 1.0e3)
            f[self.c] = h / (tray_flow * latent)

        sum_kx = 0.0
        for i in range(self.c):
            sum_kx += self.k[j][i] * self.liq[j][i] / max(lj, 1.0e-20)
        f[self.c + 1] = sum_kx - 1.0
        return f

    # ---- the thermodynamics -----------------------------------------------------------

    def ln_phi(self, z: list[float], t: float, p: float, *, liquid: bool) -> list[float]:
        """The forced-phase fugacity coefficients' logarithms at `(T, P, z)`.

        The class's `computeSinglePhaseFugacityCoefficients`: one phase, its root **forced** to
        the side asked for rather than chosen by a stability test - `K_i = phi_L_i(x) /
        phi_V_i(y)` at the tray's own compositions, where a flash would move `x` to the split's
        `x*` and leave a summation residual floor.
        """
        reduced = reduced_parameters(self.mixture, t, p)
        state = phase_state(reduced, self.mixture.kij, list(z), liquid=liquid)
        return list(state.ln_phi)

    def single_phase_enthalpy(self, z: list[float], t: float, p: float, *, liquid: bool) -> float:
        """One forced phase's molar enthalpy at `(T, P, z)`, J/mol."""
        reduced = reduced_parameters(self.mixture, t, p)
        state = phase_state(reduced, self.mixture.kij, list(z), liquid=liquid)
        result = molar_enthalpy_entropy(
            self.mixture,
            self.ideal_gas,
            quantity(t, "K"),
            quantity(p, "Pa"),
            list(z),
            state.z,
        )
        return float(result.h.to("J/mol").magnitude)

    def wilson_k(self, i: int, t: float, p_bar: float) -> float:
        """The class's own `wilsonK`, at its own constant."""
        component = self.mixture.components[i]
        tc = component.Tc.to_base_units().magnitude
        pc = component.Pc.to_base_units().magnitude
        return (pc / 1.0e5 / p_bar) * math.exp(
            WILSON_CONSTANT * (1.0 + component.omega) * (1.0 - tc / t)
        )

    def evaluate_thermo(self) -> None:
        """Evaluate every tray's thermodynamics."""
        for j in range(self.n):
            self.evaluate_thermo_for_tray(j)

    def evaluate_thermo_for_tray(self, j: int) -> None:
        """One tray: its K-values, its enthalpies and its derived vapour flows.

        `V[j]` is **not** recomputed - it is a free variable of the solver - and the vapour
        component flows are derived from it and the K-values at the end.
        """
        c = self.c
        self.l[j] = max(sum(self.liq[j]), 1.0e-20)
        x = [self.liq[j][i] / self.l[j] for i in range(c)]
        p = self.p[j]
        t = self.t[j]
        phi_l = self.ln_phi(x, t, p, liquid=True)

        # The vapour guess from the previous K-values, or Wilson where there is none.
        y = [0.0] * c
        sum_guess = 0.0
        for i in range(c):
            guess = self.k[j][i]
            k_guess = guess if 1.0e-20 < guess < 1.0e15 else self.wilson_k(i, t, p / 1.0e5)
            y[i] = k_guess * x[i]
            sum_guess += y[i]
        y = [value / sum_guess for value in y] if sum_guess > 1.0e-20 else list(x)

        for _ in range(THERMO_K_VALUE_ITERATIONS):
            phi_v = self.ln_phi(y, t, p, liquid=False)
            total = 0.0
            max_update = 0.0
            for i in range(c):
                previous = self.k[j][i]
                if not 1.0e-20 < previous < 1.0e15:
                    previous = self.wilson_k(i, t, p / 1.0e5)
                k_new = min(max(math.exp(phi_l[i] - phi_v[i]), 1.0e-15), 1.0e15)
                max_update = max(max_update, abs(math.log(k_new / previous)))
                self.k[j][i] = k_new
                y[i] = k_new * x[i]
                total += y[i]
            y = [value / total for value in y] if total > 1.0e-20 else list(x)
            if max_update <= THERMO_K_VALUE_TOLERANCE and total > 1.0e-20:
                break

        eta = self.tray_eta(j)
        if eta < 1.0 - 1.0e-10:
            for i in range(c):
                if self.k[j][i] > 0.0:
                    self.k[j][i] = self.k[j][i] ** eta

        sum_kx = sum(self.k[j][i] * x[i] for i in range(c))
        for i in range(c):
            y[i] = self.k[j][i] * x[i] / sum_kx if sum_kx > 1.0e-20 else x[i]

        self.hl[j] = self.single_phase_enthalpy(x, t, p, liquid=True)
        self.hv[j] = self.single_phase_enthalpy(y, t, p, liquid=False)
        for i in range(c):
            self.vap[j][i] = self.k[j][i] * x[i] * self.v[j]

    def tray_eta(self, tray: int) -> float:
        """The Murphree proxy's efficiency for one tray: one everywhere.

        The class forces both ends to one and this port refuses the parameter at its boundary,
        so the proxy is the identity - and the parameter is the class's own `trayEta[j]`, kept
        so the call site reads as the class's does.
        """
        del tray
        return 1.0

    # ---- the closure and the Jacobian's base ------------------------------------------

    def mass_balance_error(self) -> float:
        """`computeMassBalanceError`: the products against the feed, relative."""
        total_feed = sum(
            self.feed_liq[j][i] + self.feed_vap[j][i] for j in range(self.n) for i in range(self.c)
        )
        product = (
            self.internal_vapour_fraction[self.n - 1] * self.v[self.n - 1]
            + self.internal_liquid_fraction[0] * self.l[0]
        )
        for j in range(self.n):
            product += (1.0 - self.internal_vapour_fraction[j]) * self.v[j]
            product += (1.0 - self.internal_liquid_fraction[j]) * self.l[j]
        return abs(product - total_feed) / max(total_feed, 1.0e-20)

    def max_component_imbalance(self) -> float:
        """`computeMaxComponentImbalance`: the worst per-tray, per-component leak."""
        total_feed = sum(
            self.feed_liq[j][i] + self.feed_vap[j][i] for j in range(self.n) for i in range(self.c)
        )
        denom = max(total_feed, 1.0e-20)
        worst = 0.0
        for j in range(self.n):
            for i in range(self.c):
                balance = self.liq[j][i] + self.vap[j][i]
                if j + 1 < self.n:
                    balance -= self.internal_liquid_fraction[j + 1] * self.liq[j + 1][i]
                if j > 0:
                    balance -= self.internal_vapour_fraction[j - 1] * self.vap[j - 1][i]
                balance -= self.feed_liq[j][i] + self.feed_vap[j][i]
                worst = max(worst, abs(balance) / denom)
        return worst

    def mesh_closure_acceptable(self) -> bool:
        """`meshClosureAcceptable`: the component balances close, relative to the feed."""
        return self.max_component_imbalance() <= COMPONENT_IMBALANCE_TOLERANCE

    def save(self) -> tuple[list[list[float]], list[float], list[float]]:
        """Save the primary state."""
        return ([row[:] for row in self.liq], list(self.t), list(self.v))

    def restore(self, saved: tuple[list[list[float]], list[float], list[float]]) -> None:
        """Restore the primary state."""
        self.liq = [row[:] for row in saved[0]]
        self.t = list(saved[1])
        self.v = list(saved[2])

    def save_derived(
        self,
    ) -> tuple[list[list[float]], list[list[float]], list[float], list[float], list[float]]:
        """The derived state a Jacobian build holds still."""
        return (
            [row[:] for row in self.k],
            [row[:] for row in self.vap],
            list(self.l),
            list(self.hl),
            list(self.hv),
        )

    def restore_derived_for_tray(
        self,
        tray: int,
        base: tuple[list[list[float]], list[list[float]], list[float], list[float], list[float]],
    ) -> None:
        """One tray's derived state, put back as it was - **without** re-evaluating it.

        `evaluateThermoForTray` advances its own K-value fixed point, so a restore that
        re-evaluated would leave the next Jacobian column starting from a different base.
        """
        self.k[tray] = list(base[0][tray])
        self.vap[tray] = list(base[1][tray])
        self.l[tray] = base[2][tray]
        self.hl[tray] = base[3][tray]
        self.hv[tray] = base[4][tray]

    def jacobian(self, f0: list[float]) -> list[list[list[list[float]]]]:
        """The numerical Jacobian, as the three block bands.

        **Only the perturbed tray's thermodynamics is re-evaluated**, which is what makes the
        matrix block-tridiagonal: a tray's K-values and enthalpies depend on its own
        temperature and composition and on nothing else, so a perturbed variable can reach the
        equations of its own tray and of its two neighbours and no further.

        Returns `(lower, diagonal, upper)`, one dense `m x m` block per tray.
        """
        lower = [[[0.0] * self.m for _ in range(self.m)] for _ in range(self.n)]
        diagonal = [[[0.0] * self.m for _ in range(self.m)] for _ in range(self.n)]
        upper = [[[0.0] * self.m for _ in range(self.m)] for _ in range(self.n)]
        base = self.save_derived()
        for tray in range(self.n):
            for k in range(self.m):
                original = self.variable(tray, k)
                h = max(abs(original) * PERTURBATION, MIN_PERTURBATION)
                self.set_variable(tray, k, original + h)
                self.evaluate_thermo_for_tray(tray)
                for j in range(max(0, tray - 1), min(self.n - 1, tray + 1) + 1):
                    perturbed = self.residual_for_tray(j)
                    block = diagonal[j] if j == tray else (lower[j] if j > tray else upper[j])
                    for equation in range(self.m):
                        block[equation][k] = (perturbed[equation] - f0[j * self.m + equation]) / h
                self.set_variable(tray, k, original)
                self.restore_derived_for_tray(tray, base)
        return [lower, diagonal, upper]

    # ---- the seed and the outcome ------------------------------------------------------

    def initialize_from_column(self, seed: _States, total_feed: float) -> None:
        """`initializeTrayStateFromColumn`: the MESH variables from the converged trays.

        **This is the class's warm start, and this twin always takes it.**
        """
        previous_total = seed.tray_gas_n[-1] + seed.tray_liquid_n[0]
        if not previous_total > 1.0e-12:
            raise InvalidInputError(
                "column",
                f"the substitution core's products carry {previous_total} mol/s, so there is "
                f"no scale to seed the mesh solve from",
            )
        factor = total_feed / (previous_total * MOL_PER_HOUR)
        for j, temperature in enumerate(seed.tray_temperature):
            if not seed.tray_gas_n[j] > 0.0 or not seed.tray_liquid_n[j] > 0.0:
                raise InvalidInputError(
                    "column",
                    f"tray {j} left the substitution core with {seed.tray_gas_n[j]} mol/s of "
                    f"vapour and {seed.tray_liquid_n[j]} of liquid, and the mesh solve cannot "
                    f"start from that",
                )
            self.t[j] = (
                self.fixed_temperature[j]
                if self.fixed_temperature[j] == self.fixed_temperature[j]
                else temperature
            )
            self.v[j] = seed.tray_gas_n[j] * MOL_PER_HOUR * factor
            self.l[j] = seed.tray_liquid_n[j] * MOL_PER_HOUR * factor
            for i in range(self.c):
                self.liq[j][i] = max(seed.tray_liquid_z[j][i] * self.l[j], 1.0e-20)

    def duties(self) -> tuple[float, float]:
        """The two duties in W, from the ends' energy balance against their traffic."""
        top = self.n - 1
        condenser_in = (
            self.internal_vapour_fraction[top - 1] * self.v[top - 1] * self.hv[top - 1]
            if top > 0
            else 0.0
        )
        reboiler_in = (
            self.internal_liquid_fraction[1] * self.l[1] * self.hl[1] if self.n > 1 else 0.0
        )
        condenser_in += (
            self.feed_l_total[top] * self.feed_hl + self.feed_v_total[top] * self.feed_hv
        )
        reboiler_in += self.feed_l_total[0] * self.feed_hl + self.feed_v_total[0] * self.feed_hv
        condenser_out = self.v[top] * self.hv[top] + self.l[top] * self.hl[top]
        reboiler_out = self.v[0] * self.hv[0] + self.l[0] * self.hl[0]
        return (
            (condenser_out - condenser_in) / MOL_PER_HOUR,
            (reboiler_out - reboiler_in) / MOL_PER_HOUR,
        )


def vector_norm(v: list[float]) -> float:
    """The Euclidean norm of a vector, which is the class's `vectorNorm`."""
    return math.sqrt(sum(value * value for value in v))


def _inverse(matrix: list[list[float]]) -> list[list[float]] | None:
    """The inverse of a small dense matrix by Gauss-Jordan, or `None` where it is singular."""
    m = len(matrix)
    augmented = [[0.0] * (2 * m) for _ in range(m)]
    for i in range(m):
        augmented[i][:m] = list(matrix[i][:m])
        augmented[i][m + i] = 1.0
    for column in range(m):
        pivot_row = column
        largest = abs(augmented[column][column])
        for row in range(column + 1, m):
            if abs(augmented[row][column]) > largest:
                largest = abs(augmented[row][column])
                pivot_row = row
        if largest < 1e-30:
            return None
        augmented[column], augmented[pivot_row] = augmented[pivot_row], augmented[column]
        pivot = augmented[column][column]
        for k in range(2 * m):
            augmented[column][k] /= pivot
        for row in range(m):
            if row == column:
                continue
            factor = augmented[row][column]
            for k in range(2 * m):
                augmented[row][k] -= factor * augmented[column][k]
    return [row[m:] for row in augmented]


def _multiply(a: list[list[float]], b: list[list[float]]) -> list[list[float]]:
    """The product of two square blocks."""
    m = len(a)
    return [[sum(a[i][k] * b[k][j] for k in range(m)) for j in range(m)] for i in range(m)]


def _multiply_vector(a: list[list[float]], v: list[float]) -> list[float]:
    """The product of a square block and a vector."""
    return [sum(a[i][k] * v[k] for k in range(len(v))) for i in range(len(a))]


def solve_block_tridiagonal(
    bands: list[list[list[list[float]]]], rhs: list[list[float]]
) -> list[list[float]] | None:
    """The block Thomas algorithm, which is `solveBlockTridiagonal`.

    The forward sweep *inverts* each reduced diagonal block and multiplies, where an LU of the
    block would be cheaper and better conditioned; the port keeps the inverse because the guard
    against a singular block is the class's own test.
    """
    lower, diagonal, upper = bands
    n = len(diagonal)
    m = len(diagonal[0])
    prime = [[list(row) for row in block] for block in diagonal]
    prime_rhs = [list(row) for row in rhs]
    for j in range(1, n):
        inverse = _inverse(prime[j - 1])
        if inverse is None:
            return None
        factor = _multiply(lower[j], inverse)
        correction = _multiply(factor, upper[j - 1])
        for r in range(m):
            for c in range(m):
                prime[j][r][c] = diagonal[j][r][c] - correction[r][c]
        applied = _multiply_vector(factor, prime_rhs[j - 1])
        for r in range(m):
            prime_rhs[j][r] -= applied[r]
    solution = [[0.0] * m for _ in range(n)]
    inverse = _inverse(prime[n - 1])
    if inverse is None:
        return None
    solution[n - 1] = _multiply_vector(inverse, prime_rhs[n - 1])
    for j in range(n - 2, -1, -1):
        carried = _multiply_vector(upper[j], solution[j + 1])
        adjusted = [prime_rhs[j][r] - carried[r] for r in range(m)]
        inverse = _inverse(prime[j])
        if inverse is None:
            return None
        solution[j] = _multiply_vector(inverse, adjusted)
    return solution


def dense_solve(matrix: list[list[float]], rhs: list[float]) -> list[float] | None:
    """`solveDenseLU`, the class's own fallback for a block route that reports singular."""
    size = len(rhs)
    augmented = [[*matrix[row][:size], rhs[row]] for row in range(size)]
    for column in range(size):
        pivot = column
        largest = abs(augmented[column][column])
        for row in range(column + 1, size):
            if abs(augmented[row][column]) > largest:
                largest = abs(augmented[row][column])
                pivot = row
        if largest < 1e-300:
            return None
        augmented[column], augmented[pivot] = augmented[pivot], augmented[column]
        for row in range(column + 1, size):
            factor = augmented[row][column] / augmented[column][column]
            for k in range(column, size + 1):
                augmented[row][k] -= factor * augmented[column][k]
    dx = [0.0] * size
    for i in range(size - 1, -1, -1):
        total = augmented[i][size]
        for j in range(i + 1, size):
            total -= augmented[i][j] * dx[j]
        dx[i] = total / augmented[i][i]
    return dx


def dense_from_bands(bands: list[list[list[list[float]]]]) -> list[list[float]]:
    """`toDense`: the bands copied into one dense matrix."""
    lower, diagonal, upper = bands
    n = len(diagonal)
    m = len(diagonal[0])
    size = n * m
    dense = [[0.0] * size for _ in range(size)]
    for block in range(n):
        for r in range(m):
            for c in range(m):
                dense[block * m + r][block * m + c] = diagonal[block][r][c]
                if block > 0:
                    dense[block * m + r][(block - 1) * m + c] = lower[block][r][c]
                if block + 1 < n:
                    dense[block * m + r][(block + 1) * m + c] = upper[block][r][c]
    return dense


def apply_trust_region(mesh: Mesh, dx: list[float]) -> float:
    """`applyTrustRegion`: the class's per-variable clamp, scaled whole."""
    max_dt = 10.0
    scale = 1.0
    for j in range(mesh.n):
        base = j * mesh.m
        for i in range(mesh.c):
            cap = 0.5 * mesh.liq[j][i] + 1.0e-3 * mesh.flow_scale
            step = abs(dx[base + i])
            if step > cap and cap > 0.0:
                scale = min(scale, cap / step)
        if mesh.fixed_temperature[j] != mesh.fixed_temperature[j]:
            step = abs(dx[base + mesh.c])
            if step > max_dt:
                scale = min(scale, max_dt / step)
        cap = 0.5 * mesh.v[j] + 0.05 * mesh.flow_scale
        step = abs(dx[base + mesh.c + 1])
        if step > cap and cap > 0.0:
            scale = min(scale, cap / step)
    if scale < 1.0:
        for k in range(len(dx)):
            dx[k] *= scale
    return scale


def apply_update(mesh: Mesh, dx: list[float], alpha: float) -> None:
    """`applyUpdate`: the step, with the physical bounds the class enforces."""
    for j in range(mesh.n):
        base = j * mesh.m
        for i in range(mesh.c):
            mesh.liq[j][i] = max(mesh.liq[j][i] - alpha * dx[base + i], 1.0e-20)
        if mesh.fixed_temperature[j] != mesh.fixed_temperature[j]:
            mesh.t[j] = min(max(mesh.t[j] - alpha * dx[base + mesh.c], 100.0), 1000.0)
        mesh.v[j] = max(mesh.v[j] - alpha * dx[base + mesh.c + 1], 0.0)


def line_search(
    mesh: Mesh, dx: list[float], current_norm: float
) -> tuple[float, list[float], float] | None:
    """`lineSearch`: Armijo backtracking on `||F||`, accepting the **best** trial it saw."""
    armijo = 1.0e-4
    backtrack = 0.5
    max_backtrack = 15
    saved = mesh.save()
    alpha = 1.0
    best: tuple[float, list[float], float] | None = None
    last = float("nan")
    for _ in range(max_backtrack):
        mesh.restore(saved)
        apply_update(mesh, dx, alpha)
        mesh.evaluate_thermo()
        residual = mesh.residual()
        norm = vector_norm(residual)
        last = alpha
        if math.isfinite(norm) and (best is None or norm < best[2]):
            best = (alpha, residual, norm)
        if math.isfinite(norm) and norm < (1.0 - armijo * alpha) * current_norm:
            return (alpha, residual, norm)
        if alpha < 0.01:
            break
        alpha *= backtrack
    if best is None:
        mesh.restore(saved)
        return None
    best_alpha = best[0]
    if best_alpha != last:
        mesh.restore(saved)
        apply_update(mesh, dx, best_alpha)
        mesh.evaluate_thermo()
        residual = mesh.residual()
        norm = vector_norm(residual)
        return (best_alpha, residual, norm)
    return best


def mesh_states(
    components: list[str],
    feed_n: float,
    feed_z: list[float],
    feed_t: float,
    feed_p: float,
    feed_stage: int,
    pressures: list[float],
    fixed_temperature: list[float],
    seed: _States,
) -> _States:
    """The whole solve: the seed, the Newton loop, and the outcome.

    `pressures` is one absolute pressure per tray in Pa and `fixed_temperature` one pin per
    tray in K, `NaN` where the energy equation stands - both the caller's, because the
    pressures are the column's rule and the pins are the ends' specifications.
    """
    trays = len(pressures)
    if not any(value == value for value in fixed_temperature):
        raise InvalidInputError(
            "solver_type",
            "naphtali_sandholm without a temperature pin would run `solveBubblePointMethod`, "
            "whose V-cascade this port does not carry: pin the reboiler's temperature, which "
            "is the class's own `useOverallMBClosure` condition",
        )

    mesh = Mesh(components, trays, pressures)
    mesh.fixed_temperature = list(fixed_temperature)

    # The feed's own flash, split into the vapour and liquid it carries onto its tray.
    mixture, _ideal_gas = _components.mixture_of(mesh.components, eos="pr")
    from azoth.eos.reference.pt_flash import pt_flash

    flash = pt_flash(mixture, quantity(feed_t, "K"), quantity(feed_p, "Pa"), list(feed_z))
    if flash.phase == Phase.TWO_PHASE:
        beta = flash.beta
        if beta is None:
            raise InvalidInputError("feed", "the feed's flash is two-phase with no vapour fraction")
        n_vap, z_vap = feed_n * beta, list(flash.y)
        n_liq, z_liq = feed_n * (1.0 - beta), list(flash.x)
    elif flash.phase == Phase.ALL_VAPOUR:
        n_vap, z_vap, n_liq, z_liq = feed_n, list(feed_z), 0.0, list(feed_z)
    elif flash.phase == Phase.ALL_LIQUID:
        n_vap, z_vap, n_liq, z_liq = 0.0, list(feed_z), feed_n, list(feed_z)
    else:
        raise InvalidInputError(
            "feed",
            "the feed's flash is a trivial solution, so nothing proves which phase it is and "
            "the solver has no split to seed from",
        )

    # The feed enters the tray the caller states, which is the substitution core's own stage 0
    # for a column and the top stage for the feed the class adds there.
    stage = feed_stage
    for i in range(len(components)):
        mesh.feed_liq[stage][i] = n_liq * MOL_PER_HOUR * z_liq[i]
        mesh.feed_vap[stage][i] = n_vap * MOL_PER_HOUR * z_vap[i]
    mesh.feed_l_total[stage] = n_liq * MOL_PER_HOUR
    mesh.feed_v_total[stage] = n_vap * MOL_PER_HOUR
    if n_liq > 0.0:
        mesh.feed_hl = mesh.single_phase_enthalpy(z_liq, feed_t, feed_p, liquid=True)
    if n_vap > 0.0:
        mesh.feed_hv = mesh.single_phase_enthalpy(z_vap, feed_t, feed_p, liquid=False)

    total_feed: float = sum(
        mesh.feed_liq[j][i] + mesh.feed_vap[j][i]
        for j in range(trays)
        for i in range(len(components))
    )
    mesh.flow_scale = max(total_feed / trays, 1.0)
    mesh.initialize_from_column(seed, total_feed)
    mesh.evaluate_thermo()

    residual = mesh.residual()
    norm = vector_norm(residual)
    best = mesh.save()
    best_norm = norm
    failed_steps = 0
    stagnation = 0
    non_descent = 0
    completed = 0
    previous_best = best_norm

    for iteration in range(1, MAX_ITERATIONS + 1):
        completed = iteration
        if norm < TOLERANCE:
            return _outcome(mesh, iteration - 1, norm)
        bands = mesh.jacobian(residual)
        rhs = [residual[j * mesh.m : (j + 1) * mesh.m] for j in range(trays)]
        step = solve_block_tridiagonal(bands, rhs)
        if step is None:
            flat = [value for row in rhs for value in row]
            dense = dense_solve(dense_from_bands(bands), flat)
            step = (
                None
                if dense is None
                else [dense[j * mesh.m : (j + 1) * mesh.m] for j in range(trays)]
            )
        if step is None:
            mesh.restore(best)
            mesh.evaluate_thermo()
            residual = mesh.residual()
            norm = vector_norm(residual)
            failed_steps += 1
            if failed_steps > 5:
                raise _not_converged(iteration, norm)
            continue
        dx = [value for block in step for value in block]
        apply_trust_region(mesh, dx)
        result = line_search(mesh, dx, norm)
        if result is None:
            mesh.restore(best)
            mesh.evaluate_thermo()
            break
        _, trial_residual, trial_norm = result
        if not math.isfinite(trial_norm):
            mesh.restore(best)
            mesh.evaluate_thermo()
            residual = mesh.residual()
            norm = vector_norm(residual)
            failed_steps += 1
            if failed_steps > 5:
                raise _not_converged(iteration, norm)
            continue
        if trial_norm > 1.5 * norm and trial_norm > TOLERANCE:
            mesh.restore(best)
            mesh.evaluate_thermo()
            residual = mesh.residual()
            norm = vector_norm(residual)
            failed_steps += 1
            if failed_steps > 10:
                break
            continue
        if trial_norm >= norm:
            non_descent += 1
            if non_descent >= MAX_NON_DESCENT_LINE_SEARCH_STEPS:
                mesh.restore(best)
                mesh.evaluate_thermo()
                break
        residual = trial_residual
        norm = trial_norm
        failed_steps = 0
        if norm < best_norm:
            best_norm = norm
            best = mesh.save()
        if iteration % 10 == 0:
            if best_norm > 0.95 * previous_best:
                stagnation += 1
                if stagnation >= 2:
                    mesh.restore(best)
                    mesh.evaluate_thermo()
                    if mesh.mass_balance_error() < 0.005 and mesh.mesh_closure_acceptable():
                        return _outcome(mesh, iteration, best_norm)
                    break
            else:
                stagnation = 0
            previous_best = best_norm

    mesh.restore(best)
    mesh.evaluate_thermo()
    norm = vector_norm(mesh.residual())
    if norm < TOLERANCE * 100.0 and mesh.mesh_closure_acceptable():
        return _outcome(mesh, completed, norm)
    raise _not_converged(completed, norm)


def _not_converged(iterations: int, residual: float) -> SolverNotConvergedError:
    """A refusal naming the residual the solve stopped at."""
    return SolverNotConvergedError(iterations, residual, TOLERANCE)


def _outcome(mesh: Mesh, iterations: int, residual: float) -> _States:
    """The solution, as the reference's own record: the trays, the products and the closure."""
    from azoth.process.reference._column_stage import stream_at
    from azoth.process.reference.distillation_column import _States

    liquid_n = [mesh.l[j] / MOL_PER_HOUR for j in range(mesh.n)]
    gas_n = [mesh.v[j] / MOL_PER_HOUR for j in range(mesh.n)]
    liquid_z = [[mesh.liq[j][i] / mesh.l[j] for i in range(mesh.c)] for j in range(mesh.n)]
    gas_z = [[mesh.vap[j][i] / mesh.v[j] for i in range(mesh.c)] for j in range(mesh.n)]

    top = mesh.n - 1
    distillate = stream_at(mesh.components, gas_n[top], gas_z[top], mesh.t[top], mesh.p[top])
    bottoms = stream_at(mesh.components, liquid_n[0], liquid_z[0], mesh.t[0], mesh.p[0])
    condenser_duty, reboiler_duty = mesh.duties()
    feed_enthalpy = (
        sum(
            mesh.feed_l_total[j] * mesh.feed_hl + mesh.feed_v_total[j] * mesh.feed_hv
            for j in range(mesh.n)
        )
        / MOL_PER_HOUR
    )
    products = distillate["n"] * distillate["h"] + bottoms["n"] * bottoms["h"]
    energy_residual = (
        abs(feed_enthalpy + condenser_duty + reboiler_duty - products) / abs(feed_enthalpy)
        if abs(feed_enthalpy) > 0.0
        else 0.0
    )
    mass_residual = 0.0
    for i in range(mesh.c):
        supplied = (
            sum(mesh.feed_liq[j][i] + mesh.feed_vap[j][i] for j in range(mesh.n)) / MOL_PER_HOUR
        )
        if abs(supplied) > 1.0e-12:
            delivered = distillate["n"] * distillate["z"][i] + bottoms["n"] * bottoms["z"][i]
            mass_residual = max(mass_residual, abs(supplied - delivered) / abs(supplied))

    return _States(
        tray_temperature=tuple(mesh.t),
        tray_pressure=tuple(mesh.p),
        tray_gas_n=tuple(gas_n),
        tray_liquid_n=tuple(liquid_n),
        tray_gas_z=tuple(tuple(row) for row in gas_z),
        tray_liquid_z=tuple(tuple(row) for row in liquid_z),
        distillate_n=distillate["n"],
        distillate_z=tuple(distillate["z"]),
        distillate_h=distillate["h"],
        distillate_t=float(distillate["t"]),
        distillate_p=float(distillate["p"]),
        bottoms_n=bottoms["n"],
        bottoms_z=tuple(bottoms["z"]),
        bottoms_h=bottoms["h"],
        bottoms_t=float(bottoms["t"]),
        bottoms_p=float(bottoms["p"]),
        condenser_duty=condenser_duty,
        reboiler_duty=reboiler_duty,
        iterations=iterations,
        # **The temperature residual is the scaled MESH norm here**, which is this solver's own
        # convergence measure; NeqSim's column reports a `NaN` in that field under this strategy.
        temperature_residual=residual,
        mass_residual=mass_residual,
        energy_residual=energy_residual,
        # **The mesh solve has no reactive route, and the model refuses the pair** - see
        # `distillation_column`'s own refusal: this solve takes its fugacities from the
        # MESH equations, so the flag could only be ignored.
        warnings=(),
    )
