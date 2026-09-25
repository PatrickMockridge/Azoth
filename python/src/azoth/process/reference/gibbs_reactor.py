"""``process.gibbs_reactor`` - the equilibrium composition of a fluid at its own state.

Spec: ``specs/models/process/gibbs_reactor.toml``

The Python twin of ``crates/azoth-process/src/reactor/gibbs_solver.rs`` over
``crates/azoth-process/src/kernels/gibbs_reactor.rs``, written to mirror them.

# It is the class's solver, not a better one

``GibbsReactor.run`` never calls ``ChemicalEquilibrium``. It carries a Lagrange-multiplier Newton
iteration over the element balances, against its own species database, whose formation properties
are not the databank's - ``GIBBSENERGYOFFORMATION`` is ``used`` from P10 and its values are not
these. The residual is

    F_i = Gf0_i(T) + RT ln phi_i + RT ln y_i + RT ln(P/1 bara) - sum_j lambda_j a_ij

with one element-balance row per *active* element, and the composition takes ``damping_composition``
of each Newton step while the multipliers take the whole of it.

# Four behaviours that are the class's and are reproduced rather than tidied

**The convergence test is on the undamped step.** ``deltaX = -J^-1 F`` and the norm tested is its,
before ``alpha`` is applied. **``min_iterations`` is a floor that bites** - three of the capture's
six rows stop at exactly 100. **Reaching ``max_iterations`` returns success** while ``converged`` is
false. And **a component the database does not carry is still a variable**: a composition row, no
element entry, an objective reading as zero, and a skipped update.

# The element columns are shifted by one, and the eighth is a charge

``elementNames`` is ``{"O","N","C","H","S","Ar","Z"}`` - seven names for an eight-column table - so
the element the class calls ``Ar`` is the file's ``Na`` column and the one it calls ``Z`` is the
file's ``Ar`` column, which is zero on every row and is therefore never active. The eighth column
is a **charge** (``OH-`` is -1) and no loop in the class reaches it: **the solve does not balance
charge**.

# It is ill-conditioned, which is why the history is the oracle

The ``1e-6`` floor puts an ``RT/n`` entry of about ``1e6`` in the Jacobian against entries near one,
so ``cond(J)`` is ``6.1e6`` on the ammonia row, and the residual is a cancelling sum. One iteration
from the same state agrees between the two implementations to ``1e-13`` on the major species; over
the hundred damped steps the class takes, that is enough to separate two runs.
"""

from __future__ import annotations

import csv
import io
import math
from collections.abc import Sequence
from dataclasses import dataclass
from functools import cache
from typing import Any, Final

from azoth._data import find
from azoth.core.errors import InvalidInputError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import GibbsReactorResult
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning
from azoth.eos import components as _databank
from azoth.eos.reference._mixture_state import (
    phase_derivatives,
    phase_state,
    reduced_parameters,
)

#: ``GibbsReactor.java``'s own ``elementNames``, which has **seven** entries for an eight-column
#: table. Index 5 is named ``Ar`` and reads the file's ``Na`` column; index 6 is named ``Z`` and
#: reads the file's ``Ar`` column.
ELEMENT_NAMES: Final = ("O", "N", "C", "H", "S", "Ar", "Z")

#: ``REFERENCE_TEMPERATURE``, the 298.15 K the formation properties are stated at.
REFERENCE_TEMPERATURE: Final = 298.15

#: ``R_KJ``, the gas constant the class's Gibbs arithmetic uses, kJ/(mol*K). **A fourth gas
#: constant in this workspace** and not interchangeable with ISO 6976's 8.314510, the Sm3
#: conversions' 8.3144621 or ``KineticReaction``'s 8.31446.
R_KJ: Final = 8.314462618e-3

#: ``MIN_MOLES``, the floor the class puts under a database species before the solve.
MIN_MOLES: Final = 1.0e-6

#: ``MIN_JACOBIAN_MOLES``, the same floor where the Jacobian divides by it.
MIN_JACOBIAN_MOLES: Final = 1.0e-6

#: ``ELEMENT_ZERO_THRESHOLD``, the class's test for "this element is not present".
ELEMENT_ZERO_THRESHOLD: Final = 1.0e-6

#: ``performIterationUpdate``'s floor on a composition after the damped step.
MIN_COMPOSITION: Final = 1.0e-15

#: The four heat-capacity coefficients the element subtraction reads, per element, in the order
#: ``O, N, C, H, S`` - ``calculateCorrectedHeatCapacityCoeffs``' own five arrays.
ELEMENT_CP: Final = (
    (12.73, 7.60e-3, -3.58e-6, 6.56e-10),
    (14.4415, -7.85e-4, 4.04e-6, -1.44e-9),
    (8.43, 0.0, 0.0, 0.0),
    (14.544, -9.60e-4, 2.00e-6, -4.35e-10),
    (17.815, 0.001, 0.000, 0.000),
)


def _table(relative: str) -> list[list[str]]:
    """One compiled reactor table, header removed."""
    text = find(relative).read_text(encoding="utf-8")
    records = list(csv.reader(io.StringIO(text)))
    return records[1:]


@dataclass(frozen=True, slots=True)
class Species:
    """One row of the class's database."""

    name: str
    #: The file's eight columns: O, N, C, H, S, Na, Ar, and a charge.
    elements: tuple[float, ...]
    hf298: float
    gf298: float
    sf298: float
    #: The two order-5 polynomials, or ``None`` for a species the coefficient file omits -
    #: the class's ``Double.isNaN`` test for its fallback branch.
    gibbs_polynomial: tuple[float, ...] | None
    enthalpy_polynomial: tuple[float, ...] | None

    def balance_elements(self) -> tuple[float, ...]:
        """The eight counts as the class's *element balance* reads them: seven, under
        ``O, N, C, H, S, Ar, Z``.

        **This drops the charge and shifts the last two names**, which is the class's own
        arithmetic: index 5 is the file's ``Na`` count and index 6 the file's ``Ar`` count.
        """
        return self.elements[:7]

    def corrected_heat_capacity(self, cp: tuple[float, ...]) -> tuple[float, ...]:
        """``calculateCorrectedHeatCapacityCoeffs``: the component's Cp less its elements'."""
        corrected = list(cp)
        for element, contribution in enumerate(ELEMENT_CP):
            for power, coefficient in enumerate(contribution):
                corrected[power] -= self.elements[element] * coefficient
        return tuple(corrected)

    def j_term(self, cp: tuple[float, ...]) -> float:
        """``calculateJ``: **the code, not the javadoc.**

        The comment subtracts every term; the code subtracts the first and adds the other
        three. The two differ for most species, so the comment cannot be followed.
        """
        da, db, dc, dd = self.corrected_heat_capacity(cp)
        tr = REFERENCE_TEMPERATURE
        return (
            self.hf298 - (da * tr - db / 2.0 * tr**2 - dc / 3.0 * tr**3 - dd / 4.0 * tr**4) / 1000.0
        )

    def i_term(self, cp: tuple[float, ...]) -> float:
        """``calculateI``, the reference-temperature term of the fallback branch."""
        da, db, dc, dd = self.corrected_heat_capacity(cp)
        tr = REFERENCE_TEMPERATURE
        return (1.0 / R_KJ) * (
            -(self.j_term(cp) / tr)
            + (da * math.log(tr) + db / 2.0 * tr + dc / 6.0 * tr**2 + dd / 12.0 * tr**3) / 1000.0
        )

    def gibbs_energy(self, temperature: float, cp: tuple[float, ...]) -> float:
        """``calculateGibbsEnergy``, kJ/mol.

        **The polynomial branch is taken if the species has one at all**, so a species whose
        six coefficients are all zero - ``oxygen``, ``hydrogen``, ``nitrogen`` - returns exactly
        ``0.0`` at every temperature rather than taking the fallback.
        """
        if self.gibbs_polynomial is not None:
            return _evaluate(self.gibbs_polynomial, temperature)
        da, db, dc, dd = self.corrected_heat_capacity(cp)
        t = temperature
        ratio = self.gf298 / (R_KJ * REFERENCE_TEMPERATURE)
        rt = (
            ratio
            + self.i_term(cp)
            + (1.0 / R_KJ)
            * (
                self.j_term(cp) / t
                + (-da * math.log(t) - db / 2.0 * t - dc / 6.0 * t**2 - dd / 12.0 * t**3) / 1000.0
            )
        )
        return rt * R_KJ * t

    def enthalpy(self, temperature: float, cp: tuple[float, ...]) -> float:
        """``calculateEnthalpy``, kJ/mol, on the same two branches."""
        if self.enthalpy_polynomial is not None:
            return _evaluate(self.enthalpy_polynomial, temperature)
        da, db, dc, dd = self.corrected_heat_capacity(cp)
        t = temperature
        tr = REFERENCE_TEMPERATURE
        return (
            self.hf298
            + (
                da * (t - tr)
                + db / 2.0 * (t**2 - tr**2)
                + dc / 3.0 * (t**3 - tr**3)
                + dd / 4.0 * (t**4 - tr**4)
            )
            / 1000.0
        )


def _evaluate(coefficients: tuple[float, ...], temperature: float) -> float:
    """``c0*T**5 + ... + c5``, the class's own expansion."""
    return sum(
        coefficient * temperature ** (5 - index) for index, coefficient in enumerate(coefficients)
    )


@cache
def database() -> dict[str, Species]:
    """The class's species database, the two compiled tables joined on the lowercased name.

    **The join lowercases both sides** because the class does: ``componentMap`` is keyed by
    ``molecule.toLowerCase()`` and ``extraCoeffMap`` by ``parts[0].trim().toLowerCase()``. The
    two files do not agree on case - the coefficient file writes ``co2`` where the species file
    writes ``CO2`` - so a case-sensitive join keeps 8 of the 15 rows and says nothing.
    """
    coefficients: dict[str, tuple[tuple[float, ...], tuple[float, ...]]] = {}
    for row in _table("data/reactors/gibbs_reactor_coeffs.csv"):
        values = tuple(float(value) for value in row[1:13])
        coefficients[row[0].lower()] = (values[:6], values[6:])

    species: dict[str, Species] = {}
    for row in _table("data/reactors/gibbs_reactor.csv"):
        name = row[0]
        polynomials = coefficients.get(name.lower())
        species[name.lower()] = Species(
            name=name,
            elements=tuple(float(value) for value in row[1:9]),
            hf298=float(row[9]),
            gf298=float(row[10]),
            sf298=float(row[11]),
            gibbs_polynomial=None if polynomials is None else polynomials[0],
            enthalpy_polynomial=None if polynomials is None else polynomials[1],
        )
    return species


@dataclass(frozen=True, slots=True)
class Solved:
    """What the solve reached."""

    moles: list[float]
    temperature: float
    lambdas: tuple[float, ...]
    converged: bool
    iterations: int
    final_error: float
    element_in: tuple[float, ...]
    element_out: tuple[float, ...]
    gibbs_history: list[float]


def _solve_linear(jacobian: list[list[float]], rhs: list[float]) -> list[float] | None:
    """``J x = rhs`` by Gaussian elimination with partial pivoting, or ``None`` if singular.

    ``solveNewtonSystem`` uses EJML's ``solve`` and falls back to ``pseudoInverse``; the fallback
    is **not carried**, so a singular system is refused rather than answered differently. It is
    reachable: two species against three active elements is over-determined, and the capture's
    ``steam_methane_hot`` row is exactly that.
    """
    n = len(rhs)
    a = [[*jacobian[i], rhs[i]] for i in range(n)]
    for column in range(n):
        pivot = max(range(column, n), key=lambda row: abs(a[row][column]))
        a[column], a[pivot] = a[pivot], a[column]
        if a[column][column] == 0.0:
            return None
        for row in range(column + 1, n):
            factor = a[row][column] / a[column][column]
            for k in range(column, n + 1):
                a[row][k] -= factor * a[column][k]
    x = [0.0] * n
    for row in range(n - 1, -1, -1):
        total = sum(a[row][k] * x[k] for k in range(row + 1, n))
        x[row] = (a[row][n] - total) / a[row][row]
    return x


def solve(
    components: list[str],
    inlet_moles: list[float],
    temperature: float,
    pressure_bara: float,
    *,
    adiabatic: bool,
    damping: float,
    max_iterations: int,
    tolerance: float,
    min_iterations: int,
) -> Solved:
    """Minimise the Gibbs energy of a fluid at a temperature and pressure."""
    table = database()
    mixture, ideal_gas = _databank.mixture_of(components, eos="pr")
    n = len(mixture)
    species = [table.get(name.lower()) for name in components]
    cp = [
        (ideal_gas.cp_a[i], ideal_gas.cp_b[i], ideal_gas.cp_c[i], ideal_gas.cp_d[i])
        for i in range(n)
    ]

    # The element balance on the *feed*, before any floor - both exclusions below use it.
    element_in = _element_balance(species, inlet_moles)

    # `determineFeedExcludedComponents`: a species that needs an element the feed has none of
    # cannot form, so it is frozen at its feed amount and dropped from the matrix.
    excluded = [
        entry is not None
        and any(
            count != 0.0
            and abs(count) > ELEMENT_ZERO_THRESHOLD
            and abs(element_in[j]) <= ELEMENT_ZERO_THRESHOLD
            for j, count in enumerate(entry.balance_elements())
        )
        for entry in species
    ]

    # The floor the class applies twice - `performGibbsMinimization` then
    # `enforceMinimumConcentrations` - and which a feed-excluded species does not get.
    moles = [
        max(inlet_moles[i], MIN_MOLES)
        if species[i] is not None and not excluded[i]
        else inlet_moles[i]
        for i in range(n)
    ]
    # **Everything the database carries and the feed does not exclude**, including a component
    # the database does *not* carry - the class's own handling of a substance it has no species
    # for, and what `testComponentNotInDatabaseMolesUnchanged` asserts.
    variables = [i for i in range(n) if not excluded[i]]
    active = _active_elements(species, element_in)

    lambdas = [0.0] * 7
    temperature_value = temperature
    history: list[float] = []
    converged = False
    iterations = 0
    final_error = math.inf
    enthalpy_old = 0.0

    for iteration in range(1, max_iterations + 1):
        iterations = iteration
        total = sum(moles)
        x = [value / total for value in moles]
        reduced = reduced_parameters(mixture, temperature_value, pressure_bara * 1.0e5)
        state = phase_state(reduced, mixture.kij, x, liquid=False)
        derivatives = phase_derivatives(
            reduced,
            mixture.kij,
            x,
            state.z,
            temperature=temperature_value,
            pressure=pressure_bara * 1.0e5,
        )

        objective = _objective(
            species, cp, moles, total, temperature_value, pressure_bara, lambdas, derivatives
        )
        jacobian = _jacobian(
            species, moles, variables, active, total, temperature_value, derivatives
        )
        vector = [objective[i] for i in variables] + [
            _element_balance(species, moles)[j] - element_in[j] for j in active
        ]
        delta = _solve_linear(jacobian, [-value for value in vector])
        if delta is None:
            raise InvalidInputError(
                "jacobian",
                f"the Newton system is singular at iteration {iteration}, and the class's "
                "fallback for that is `SimpleMatrix.pseudoInverse`, which this port does not "
                "carry - a feed with fewer species than active elements reaches it",
            )

        delta_norm = math.sqrt(sum(value * value for value in delta))
        history.append(_mixture_gibbs(species, moles, temperature_value, cp))

        if adiabatic:
            enthalpy = _mixture_enthalpy(species, moles, temperature_value, cp)
            if iteration == 1:
                enthalpy_old = enthalpy
            else:
                dh = enthalpy - enthalpy_old
                enthalpy_old = enthalpy
                cp_total = total * _molar_cp(
                    mixture, ideal_gas, temperature_value, pressure_bara, x, state.z
                )
                t_out = temperature_value - dh * 1000.0 / cp_total
                if abs(t_out - temperature_value) > 1000.0:
                    raise InvalidInputError(
                        "temperature",
                        "the adiabatic step moves the temperature by more than 1000 K, which "
                        "the class refuses",
                    )
                temperature_value = t_out

        # **The convergence test, on the undamped step, under the class's floor.**
        if (delta_norm < tolerance and iteration >= min_iterations) or iteration == max_iterations:
            converged = delta_norm < tolerance
            final_error = delta_norm
            break

        for slot, index in enumerate(variables):
            moles[index] = max(moles[index] + delta[slot] * damping, MIN_COMPOSITION)
        for slot, element in enumerate(active):
            lambdas[element] += delta[len(variables) + slot]

    element_out = _element_balance(species, moles)
    return Solved(
        moles=moles,
        temperature=temperature_value,
        lambdas=tuple(lambdas),
        converged=converged,
        iterations=iterations,
        final_error=final_error,
        element_in=tuple(element_in),
        element_out=tuple(element_out),
        gibbs_history=history,
    )


def _objective(
    species: list[Species | None],
    cp: Sequence[tuple[float, ...]],
    moles: list[float],
    total: float,
    temperature: float,
    pressure_bara: float,
    lambdas: list[float],
    derivatives: Any,
) -> list[float]:
    """``calculateObjectiveFunctionValues``, with an absent entry reading as zero."""
    rt = R_KJ * temperature
    values = [0.0] * len(moles)
    for i, entry in enumerate(species):
        if entry is None:
            continue
        y = moles[i] / total
        lagrange = sum(lambdas[j] * value for j, value in enumerate(entry.balance_elements()))
        values[i] = (
            entry.gibbs_energy(temperature, cp[i])
            + rt * derivatives.ln_phi[i]
            + rt * math.log(y)
            + rt * math.log(pressure_bara)
            - lagrange
        )
    return values


def _jacobian(
    species: list[Species | None],
    moles: list[float],
    variables: list[int],
    active: list[int],
    total: float,
    temperature: float,
    derivatives: Any,
) -> list[list[float]]:
    """``calculateJacobian``, in the class's own block order."""
    size = len(variables) + len(active)
    matrix = [[0.0] * size for _ in range(size)]
    rt = R_KJ * temperature
    for i, component_i in enumerate(variables):
        ni = max(moles[component_i], MIN_JACOBIAN_MOLES)
        for j, component_j in enumerate(variables):
            dfugdn = derivatives.d_ln_phi_dn[component_i][component_j]
            composition = 1.0 / ni if i == j else 0.0
            matrix[i][j] = rt * (composition - 1.0 / total + dfugdn)
        entry = species[component_i]
        if entry is not None:
            elements = entry.balance_elements()
            for k, element in enumerate(active):
                matrix[i][len(variables) + k] = -elements[element]
    for i, element in enumerate(active):
        for j, component_j in enumerate(variables):
            entry = species[component_j]
            if entry is not None:
                matrix[len(variables) + i][j] = entry.balance_elements()[element]
    return matrix


def _active_elements(species: list[Species | None], element_in: list[float]) -> list[int]:
    """``findActiveElements``: non-zero in the feed *and* carried by some species."""
    return [
        j
        for j in range(7)
        if abs(element_in[j]) > ELEMENT_ZERO_THRESHOLD
        and any(
            entry is not None and abs(entry.balance_elements()[j]) > ELEMENT_ZERO_THRESHOLD
            for entry in species
        )
    ]


def _element_balance(species: list[Species | None], moles: list[float]) -> list[float]:
    """``calculateElementMoleBalance``, over the class's seven names."""
    balance = [0.0] * 7
    for i, entry in enumerate(species):
        if entry is None:
            continue
        for j, value in enumerate(entry.balance_elements()):
            balance[j] += value * moles[i]
    return balance


def _mixture_gibbs(
    species: list[Species | None],
    moles: list[float],
    temperature: float,
    cp: Sequence[tuple[float, ...]],
) -> float:
    """``calculateMixtureGibbsEnergy``, kJ/mol - the total, not the molar one."""
    return sum(
        moles[i] * entry.gibbs_energy(temperature, cp[i])
        for i, entry in enumerate(species)
        if entry is not None
    )


def _mixture_enthalpy(
    species: list[Species | None],
    moles: list[float],
    temperature: float,
    cp: Sequence[tuple[float, ...]],
) -> float:
    """``calculateMixtureEnthalpy``, kJ/mol."""
    return sum(
        moles[i] * entry.enthalpy(temperature, cp[i])
        for i, entry in enumerate(species)
        if entry is not None
    )


def _molar_cp(
    mixture: Any,
    ideal_gas: Any,
    temperature: float,
    pressure_bara: float,
    x: list[float],
    z: float,
) -> float:
    """The fluid's molar heat capacity in J/(mol*K), which the class reads through
    ``getCp("J/K")`` divided by its moles."""
    from azoth.eos.reference.molar_enthalpy_entropy import molar_enthalpy_entropy

    props = molar_enthalpy_entropy(
        mixture,
        ideal_gas,
        from_si(temperature, "K"),
        from_si(pressure_bara * 1.0e5, "Pa"),
        x,
        z,
    )
    return props.cp.to("J/(mol*K)").magnitude


def gibbs_reactor(
    *,
    components: list[str],
    feed_n: Q,
    feed_z: list[float],
    feed_p: Q,
    feed_t: Q,
    energy_mode: str,
    damping_composition: float,
    max_iterations: float,
    convergence_tolerance: float,
    min_iterations: float,
) -> GibbsReactorResult:
    """Bring a feed to its Gibbs equilibrium at its own temperature and pressure.

    Spec: ``specs/models/process/gibbs_reactor.toml``.
    """
    spec = _spec()
    checks = checks_for(spec)
    warnings: list[Warning] = []
    apply_checks(
        checks.on_input,
        {
            "feed_n": input_to_si(spec, "feed_n", feed_n),
            "damping_composition": damping_composition,
            "max_iterations": max_iterations,
            "convergence_tolerance": convergence_tolerance,
        }.get,
        warnings,
    )

    total_flow = input_to_si(spec, "feed_n", feed_n)
    pressure_pa = input_to_si(spec, "feed_p", feed_p)
    temperature_k = input_to_si(spec, "feed_t", feed_t)
    if total_flow <= 0.0:
        raise InvalidInputError("feed_n", "a reactor needs a positive molar flow")

    inlet_moles = [value * total_flow for value in feed_z]
    state = solve(
        list(components),
        inlet_moles,
        temperature_k,
        pressure_pa / 1.0e5,
        adiabatic=energy_mode == "adiabatic",
        damping=damping_composition,
        max_iterations=int(max_iterations),
        tolerance=convergence_tolerance,
        min_iterations=int(min_iterations),
    )

    outlet_total = sum(state.moles)
    if outlet_total <= 0.0:
        raise InvalidInputError("moles", "the solve left no moles, so there is no outlet state")

    from azoth.eos.reference.ph_flash import enthalpy_at

    mixture, ideal_gas = _databank.mixture_of(list(components), eos="pr")
    composition = [value / outlet_total for value in state.moles]
    enthalpy, _ = enthalpy_at(mixture, ideal_gas, temperature_k, pressure_pa, composition)

    apply_checks(
        checks.derived,
        {"iterations": float(state.iterations), "final_error": state.final_error}.get,
        warnings,
    )

    return GibbsReactorResult(
        product_n=from_si(outlet_total, "mol/s"),
        product_z=tuple(composition),
        product_p=from_si(pressure_pa, "Pa"),
        product_t=from_si(state.temperature, "K"),
        product_h=from_si(enthalpy, "J/mol"),
        converged=state.converged,
        iterations=float(state.iterations),
        final_error=state.final_error,
        lagrange_multipliers=tuple(value * 1000.0 for value in state.lambdas),
        element_balance_difference=tuple(
            state.element_out[j] - state.element_in[j] for j in range(7)
        ),
        gibbs_energy_history=tuple(value * 1000.0 for value in state.gibbs_history),
        warnings=tuple(warnings),
    )


def _spec() -> dict[str, Any]:
    """The model spec, which both implementations apply."""
    from azoth._models_gen import model

    return model("process.gibbs_reactor")
