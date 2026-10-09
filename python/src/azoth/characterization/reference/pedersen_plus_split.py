"""``characterization.pedersen_plus_split`` - a plus fraction divided into carbon-number cuts.

Spec: ``specs/models/characterization/pedersen_plus_split.toml``. Oracle:
``validation/neqsim/captures/plus_fraction_probe.tsv``.

Mirrors ``crates/azoth-characterization/src/pedersen_plus_split.rs`` line for line.

# Two damped Newton solves, not one

``z_i = exp(a + b*CN)`` and ``rho_i = c + d*ln(CN)``; ``(a, b)`` is fixed against the plus
fraction's mole fraction and mean molar mass, written back, and only then ``(c, d)`` against the
tabulated anchor gravity and the plus fraction's own.

Four things in NeqSim's solver are reproduced rather than repaired: the two loops damp by
``iter/(iter + 50)`` and ``iter/(iter + 5)``; both force three steps whatever their residual, and
stop on the residual measured *before* the final step; every abundance iteration also runs an
inner fixed point on the gravity slope with the intercept held at ``0.80``; and the gravity solve
starts from the constructor's own ``(0.80, 0.0408709562)`` rather than from what the abundance
phase left. Three Jacobian entries are not the derivatives they stand for and are not corrected.

# The linear solve is Jama's

``Jac.solve(fvec)`` is a 2x2 LU decomposition with partial pivoting, and the pivot test matters:
the gravity Jacobian's ``(1,0)`` sits near one against an ``(0,0)`` of exactly one. The closed
form agrees to twelve places and then stops agreeing.

# A range ends at C80's table, not at C200's

``PVTsimMolarMass`` runs C6 to C200 and ``PVTsimDensities`` stops at C80, and the class indexes
both as ``[CN - 6]``. The gravity table is read at exactly one row - the anchor below - so the
short table is enough, and it is kept here at its own length so the asymmetry is visible.
"""

from __future__ import annotations

import csv
import math
from functools import cache

from azoth import _models_gen
from azoth._data import find
from azoth.core.errors import InvalidInputError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import PedersenPlusSplitResult
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning

MODEL_ID = "characterization.pedersen_plus_split"

PLUS_FRACTION_CSV = "data/characterization/plus_fraction.csv"

#: The carbon number the tables begin at, which is the class's own ``[i - 6]`` index base.
TABLE_FIRST_CARBON_NUMBER = 6

#: The residual 2-norm both loops stop below, and the inner fixed point's own threshold.
RESIDUAL_TOLERANCE = 1.0e-6

#: The cap on each Newton loop.
MAX_ITERATIONS = 200

#: A loop is not allowed to stop before this many steps, whatever its residual.
MIN_ITERATIONS = 3

#: The abundance step's damping.
ABUNDANCE_DAMPING = 50.0

#: The gravity step's damping, a tenth of the other's.
GRAVITY_DAMPING = 5.0

#: The cap on the inner fixed point, five times the outer loops'.
DENSITY_SLOPE_MAX_ITERATIONS = 1000

#: ``coefs[0]`` when the abundance solve starts, seeded rather than solved.
Z_SEED_INTERCEPT = 0.1

#: ``coefs[2]``'s constructor value, held fixed by the abundance phase.
DENSITY_SEED_INTERCEPT = 0.80

#: ``coefs[3]``'s constructor value, which the gravity solve starts from regardless.
DENSITY_SEED_SLOPE = 0.0408709562


@cache
def _tables() -> tuple[tuple[float, ...], tuple[float, ...]]:
    """The two PVTsim tables: molar masses from C6, and gravities from C6 to the table's end."""
    molar_mass: list[float] = []
    densities: list[float] = []
    previous: int | None = None
    with find(PLUS_FRACTION_CSV).open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle):
            carbon_number = int(row["carbon_number"])
            expected = (
                TABLE_FIRST_CARBON_NUMBER if previous is None else previous + 1
            )
            if carbon_number != expected:
                raise ValueError(
                    f"the plus-fraction table has C{carbon_number} after C{expected}"
                )
            previous = carbon_number
            molar_mass.append(float(row["molar_mass_g_per_mol"]))
            # Past C80 the density cell is empty, which is how the shorter table is carried.
            density = row["density_g_per_cm3"].strip()
            if density:
                if len(densities) != carbon_number - TABLE_FIRST_CARBON_NUMBER:
                    raise ValueError(
                        f"the plus-fraction table has C{carbon_number}'s density out of "
                        f"sequence with the rows above it"
                    )
                densities.append(float(density))
    return tuple(molar_mass), tuple(densities)


def _dimensionless(spec: dict[str, object], name: str, value: float | Q) -> float:
    """A declared dimensionless input as its SI magnitude, quantity or not.

    A count and a mole fraction carry no unit, so a bare number is one - there is no
    dimension to be wrong about. The dimensioned inputs stay strict: a bare number where a
    density is expected is the mistake this library exists to make impossible.
    """
    if isinstance(value, int | float):
        return float(value)
    return input_to_si(spec, name, value)


def _cut_molar_mass(carbon_number: int) -> float:
    """The cut's tabulated molar mass in g/mol."""
    molar_mass, _ = _tables()
    index = carbon_number - TABLE_FIRST_CARBON_NUMBER
    if index < 0 or index >= len(molar_mass):
        raise InvalidInputError(
            "last_carbon_number",
            f"the molar-mass table runs C{TABLE_FIRST_CARBON_NUMBER} to "
            f"C{TABLE_FIRST_CARBON_NUMBER + len(molar_mass) - 1}",
        )
    return molar_mass[index]


def _anchor(carbon_number: int) -> float:
    """The specific gravity the gravity line is anchored to."""
    _, densities = _tables()
    index = carbon_number - TABLE_FIRST_CARBON_NUMBER
    if index < 0 or index >= len(densities):
        raise InvalidInputError(
            "first_carbon_number",
            f"the density table runs C{TABLE_FIRST_CARBON_NUMBER} to "
            f"C{TABLE_FIRST_CARBON_NUMBER + len(densities) - 1}; the gravity line is anchored "
            f"to the first cut's own row, so a plus fraction outside that span has no anchor",
        )
    return densities[index]


def _gravity_line(c: float, d: float, carbon_number: int) -> float:
    """``c + d*ln(CN)``, the gravity line, in g/cm3."""
    return c + d * math.log(carbon_number)


def _abundance_residuals(a: float, b: float, first: int, last: int, z_plus: float,
                         m_plus: float) -> list[float]:
    """``setfvecAB``: the abundance residuals, on the **kg/mol** scale."""
    z_sum = 0.0
    m_sum = 0.0
    for i in range(first, last):
        z = math.exp(a + b * i)
        z_sum += z
        m_sum += z * (_cut_molar_mass(i) / 1000.0)
    return [z_sum - z_plus, m_sum / z_sum - m_plus]


def _abundance_jacobian(a: float, b: float, first: int, last: int) -> list[list[float]]:
    """``setJacAB``, including the two entries that are not derivatives."""
    n_tot = 0.0
    n_tot2 = 0.0
    for i in range(first, last):
        z = math.exp(a + b * i)
        n_tot += z
        n_tot2 += i * z
    m_tot1 = 0.0
    m_tot2 = 0.0
    z_squares = 0.0
    z_sum = 0.0
    z_times_carbon = 0.0
    for i in range(first, last):
        z = math.exp(a + b * i)
        m = _cut_molar_mass(i) / 1000.0
        m_tot1 += m * z
        m_tot2 += i * m * z
        z_squares += z**2.0
        z_sum += z
        z_times_carbon += i * z
    # `(1,0)` is the class's own `mTot1*zSum - mTot1*zSum`, identically zero; `(1,1)` divides by
    # the sum of squares of `z` where the derivative divides by the square of the sum.
    return [
        [n_tot, n_tot2],
        [0.0, (m_tot2 * z_sum - m_tot1 * z_times_carbon) / z_squares],
    ]


def _density_slope_fixed_point(a: float, b: float, c: float, d: float, first: int, last: int,
                               dens_plus: float) -> float:
    """``setCoefs(double[])``'s inner ``do``-``while``: the gravity slope alone."""
    iterations = 0
    while True:
        iterations += 1
        m_sum = 0.0
        dens_sum = 0.0
        for i in range(first, last):
            z = math.exp(a + b * i)
            m = _cut_molar_mass(i) / 1000.0
            m_sum += z * m
            dens_sum += z * m / _gravity_line(c, d, i)
        dens_sum = m_sum / dens_sum
        d += 1.0 * (dens_plus - dens_sum) / dens_sum * d
        if not (
            abs(dens_plus - dens_sum) > RESIDUAL_TOLERANCE
            and iterations < DENSITY_SLOPE_MAX_ITERATIONS
        ):
            return d


def _gravity_residuals(a: float, b: float, c: float, d: float, first: int, last: int,
                       anchor: float, dens_plus: float) -> list[float]:
    """``setfvecCD``: the gravity residuals, on the **g/mol** scale the class uses here."""
    temp = 0.0
    temp2 = 0.0
    for i in range(first, last):
        weight = math.exp(a + b * i) * _cut_molar_mass(i)
        temp += weight
        temp2 += weight / _gravity_line(c, d, i)
    # The anchor is the carbon number **below** the first cut against the first cut's own row.
    return [_gravity_line(c, d, first - 1) - anchor, temp / temp2 - dens_plus]


def _gravity_jacobian(a: float, b: float, c: float, d: float, first: int, last: int
                      ) -> list[list[float]]:
    """``setJacCD``, with its ``(1,1)`` left at zero and its ``(1,0)`` at the negative."""
    temp = 0.0
    temp2 = 0.0
    temp3 = 0.0
    for i in range(first, last):
        weight = math.exp(a + b * i) * _cut_molar_mass(i)
        rho = _gravity_line(c, d, i)
        temp += weight
        temp2 += weight / rho
        temp3 -= weight / rho**2.0
    derivative = temp / (temp2 * temp2) * temp3
    return [[1.0, math.log(first - 1)], [derivative, 0.0]]


def _solve_2x2(matrix: list[list[float]], rhs: list[float]) -> list[float] | None:
    """``Jama.Matrix.LUDecomposition``'s solve, transcribed for a 2x2 system. ``None`` if singular."""
    lu = [list(row) for row in matrix]
    pivot = [0, 1]
    for j in range(2):
        column = [lu[0][j], lu[1][j]]
        for i in range(2):
            total = 0.0
            for k in range(min(i, j)):
                total += lu[i][k] * column[k]
            column[i] -= total
            lu[i][j] = column[i]
        p = j
        for i in range(j + 1, 2):
            if abs(column[i]) > abs(column[p]):
                p = i
        if p != j:
            lu[p], lu[j] = lu[j], lu[p]
            pivot[p], pivot[j] = pivot[j], pivot[p]
        if lu[j][j] != 0.0:
            for i in range(j + 1, 2):
                lu[i][j] /= lu[j][j]
    if lu[0][0] == 0.0 or lu[1][1] == 0.0:
        return None
    x = [rhs[pivot[0]], rhs[pivot[1]]]
    for k in range(2):
        for i in range(k + 1, 2):
            x[i] -= x[k] * lu[i][k]
    for k in range(1, -1, -1):
        x[k] /= lu[k][k]
        for i in range(k):
            x[i] -= x[k] * lu[i][k]
    return x


def _singular(field: str) -> InvalidInputError:
    """NeqSim's ``LUDecomposition`` throws here; this refuses instead of dividing by zero."""
    return InvalidInputError(
        field,
        "the solve's Jacobian is singular at the plus fraction's own numbers, so this split has "
        "no step to take",
    )


def pedersen_plus_split(
    molar_mass: Q,
    density: Q,
    mole_fraction: float,
    first_carbon_number: int,
    last_carbon_number: int,
) -> PedersenPlusSplitResult:
    """A plus fraction divided into carbon-number cuts, by Pedersen's two Newton solves.

    ``first_carbon_number`` is the cut the range opens at and ``last_carbon_number`` is one past
    the last cut, which is the class's own half-open convention. ``80`` is ``PedersenPlusModel``'s
    default range and ``200`` is ``PedersenHeavyOilPlusModel``'s - the same two solves, over a
    longer table.

    Raises:
        InvalidInputError: if the plus fraction is lighter than the first cut the table gives it,
            if the range is empty, or if a Jacobian comes out singular.
        OutOfRangeError: if an input is outside its declared range - a carbon number below C6,
            above C80 for the first cut, or above C200 for the last.
    """
    spec = _models_gen.model(MODEL_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    values = {
        "molar_mass": input_to_si(spec, "molar_mass", molar_mass),
        "density": input_to_si(spec, "density", density),
        "mole_fraction": _dimensionless(spec, "mole_fraction", mole_fraction),
        "first_carbon_number": _dimensionless(
            spec, "first_carbon_number", first_carbon_number
        ),
        "last_carbon_number": _dimensionless(spec, "last_carbon_number", last_carbon_number),
    }
    apply_checks(checks.on_input, values.get, warnings)

    first = int(values["first_carbon_number"])
    last = int(values["last_carbon_number"])
    z_plus = values["mole_fraction"]
    m_plus = values["molar_mass"]
    # The boundary: the class carries the specific gravity in g/cm3, a thousand smaller.
    dens_plus = values["density"] / 1000.0

    if last <= first:
        raise InvalidInputError(
            "last_carbon_number",
            f"{last} leaves no cut above {first}: the range is half-open, so at least one carbon "
            f"number has to lie between them",
        )

    # `characterizePlusFraction`'s own test: the table's first cut against the plus fraction's
    # molar mass in g/mol. NeqSim logs and returns without a split; this refuses.
    first_cut_gmol = _cut_molar_mass(first)
    if first_cut_gmol > m_plus * 1000.0:
        raise InvalidInputError(
            "molar_mass",
            f"the plus fraction is lighter than the first cut the table gives it, at "
            f"{m_plus * 1000.0} g/mol against C{first}'s {first_cut_gmol} g/mol",
        )
    anchor = _anchor(first)

    # The abundance solve, each step damped, published, and followed by the inner fixed point on
    # the gravity slope that the class runs from inside `setCoefs`.
    a = Z_SEED_INTERCEPT
    b = math.log(z_plus) / first
    c = DENSITY_SEED_INTERCEPT
    d = DENSITY_SEED_SLOPE
    seed_a, seed_b = a, b
    iterations = 0
    while True:
        iterations += 1
        f = _abundance_residuals(seed_a, seed_b, first, last, z_plus, m_plus)
        jacobian = _abundance_jacobian(seed_a, seed_b, first, last)
        dx = _solve_2x2(jacobian, f)
        if dx is None:
            raise _singular("molar_mass")
        scale = iterations / (iterations + ABUNDANCE_DAMPING)
        seed_a -= dx[0] * scale
        seed_b -= dx[1] * scale
        a, b = seed_a, seed_b
        d = _density_slope_fixed_point(a, b, c, d, first, last, dens_plus)
        norm = math.sqrt(f[0] * f[0] + f[1] * f[1])
        if not ((norm > RESIDUAL_TOLERANCE or iterations < MIN_ITERATIONS)
                and iterations < MAX_ITERATIONS):
            break

    # The gravity solve, from the solver's own seed rather than from what the walk above left.
    seed_c, seed_d = DENSITY_SEED_INTERCEPT, DENSITY_SEED_SLOPE
    iterations = 0
    while True:
        iterations += 1
        f = _gravity_residuals(a, b, c, d, first, last, anchor, dens_plus)
        jacobian = _gravity_jacobian(a, b, c, d, first, last)
        dx = _solve_2x2(jacobian, f)
        if dx is None:
            raise _singular("density")
        scale = iterations / (iterations + GRAVITY_DAMPING)
        seed_c -= dx[0] * scale
        seed_d -= dx[1] * scale
        c, d = seed_c, seed_d
        norm = math.sqrt(f[0] * f[0] + f[1] * f[1])
        if not ((norm > RESIDUAL_TOLERANCE or iterations < MIN_ITERATIONS)
                and iterations < MAX_ITERATIONS):
            break

    cut_z = []
    cut_molar_mass = []
    cut_density = []
    for i in range(first, last):
        cut_z.append(math.exp(a + b * i))
        cut_molar_mass.append(from_si(_cut_molar_mass(i) / 1000.0, "kg/mol"))
        cut_density.append(from_si(_gravity_line(c, d, i) * 1000.0, "kg/m**3"))

    apply_checks(checks.derived, lambda _name: None, warnings)

    return PedersenPlusSplitResult(
        cut_z=tuple(cut_z),
        cut_molar_mass=tuple(cut_molar_mass),
        cut_density=tuple(cut_density),
        z_intercept=a,
        z_slope=b,
        density_intercept=c,
        density_slope=d,
        warnings=tuple(warnings),
    )
