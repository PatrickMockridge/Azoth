"""``reactions.chemical_equilibrium`` - the Smith-Missen reactive solve.

Spec: ``specs/models/reactions/chemical_equilibrium.toml``

The Python twin of ``crates/azoth-reactions/src/chemical_equilibrium.rs``, written to
mirror it line for line.

# The iteration

```text
M[i][k]     = delta_ik / n_i                        (the ideal form)
mu[i]       = mu_ref[i] + ln(n_i) - ln(n_t) + ln(gamma_i)
AMA         = A M^-1 A^T ,   AMU = A M^-1 mu
[AMA  c^T] [lambda]   [AMU + correction]
[c     0 ] [ tau  ] = [ sum(n_i mu_i)   ]
dn          = M^-1 (A^T lambda - mu) + n tau
```

# The three places a naive port goes wrong

**NeqSim keeps two compositions, not one.** ``committed`` is what the phase holds and the
trial is derived from it; a pass whose error got worse is *not* written back, so the next
pass re-derives from the same committed state. Collapsing the two makes the solve wander.

**A trial that goes below zero is not floored.** It is reached on the first pass: flooring
it puts ``1/MIN_MOLES = 1e60`` into ``G_1`` and drives the step to zero. See
:func:`_inner_step`.

**``step`` mixes units.** Its potentials are ``R T`` times the reduced ones while the
``A^T lambda`` it subtracts is reduced, so the factor does not cancel and ``T`` is needed.
"""

from __future__ import annotations

import math

from azoth.core.errors import InvalidInputError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import ChemicalEquilibriumResult
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning
from azoth.reactions.reference import _linalg
from azoth.reactions.reference.equilibrium_constant import GAS_CONSTANT

#: The floor every mole number is raised to before a logarithm, from
#: `ChemicalEquilibrium.MIN_MOLES`.
MIN_MOLES = 1e-60

#: The relative tolerance the conservation-correction check uses.
CONSERVATION_CORRECTION_TOLERANCE = 1e-8

#: Consecutive non-improving iterations before the solve gives up.
STAGNATION_LIMIT = 10


def chemical_equilibrium(
    a_matrix: list[list[float]],
    b: list[float],
    moles: list[float],
    chem_ref: list[float],
    log_activity: list[float],
    T: Q,
    max_iterations: float,
    tolerance: float,
) -> ChemicalEquilibriumResult:
    """The reactive equilibrium composition of a phase.

    Args:
        a_matrix: the element matrix, one row per element and one column per component,
            with the electroneutrality row last.
        b: the element amounts the solve conserves. The charge row's entry is zero.
        moles: the starting composition.
        chem_ref: each component's reduced standard-state potential, ``mu_ref / (R T)``.
        log_activity: each component's ``ln(gamma)``, held fixed for the solve.
        T: absolute temperature, read by the step search.
        max_iterations: the pass cap.
        tolerance: the error the solve is trying to reach.

    Returns:
        The composition at the answer, or at the point the solve gave up, with the pass
        count, the final error and whether it converged. **An unconverged solve is a
        result and not an exception**: one of the two captured fluids never converges.

    Raises:
        InvalidInputError: where the shapes disagree.

    Example:
        >>> import azoth
        >>> q = azoth.ureg.Quantity
        >>> r = chemical_equilibrium(
        ...     [[1.0, 0.0, 0.0, 0.0], [0.0, 2.0, 1.0, 1.0], [0.0, 0.0, -1.0, 1.0]],
        ...     [1.0, 2.0, 0.0],
        ...     [0.5, 1.0, 1e-10, 1e-10],
        ...     [-10.0, -5.0, 3.0, -3.0],
        ...     [0.0, 0.0, 0.0, 0.0],
        ...     q(298.15, "K"),
        ...     100.0,
        ...     1e-8,
        ... )
        >>> r.iterations >= 2
        True
    """
    spec = _spec()
    checks = checks_for(spec)
    warnings: list[Warning] = []
    temperature = input_to_si(spec, "T", T)

    apply_checks(
        checks.on_input,
        {"T": temperature, "tolerance": tolerance, "max_iterations": max_iterations}.get,
        warnings,
    )

    # **Every declared input goes through `input_to_si`, element by element.** The vectors
    # carry units - `b` and `moles` are moles - so the case hands them over as pint
    # quantities, and doing arithmetic on one against a bare float raises rather than
    # converts. The matrix is dimensionless and goes through for the same reason.
    a_matrix = [[_si(spec, "a_matrix", value) for value in row] for row in a_matrix]
    b = [_si(spec, "b", value) for value in b]
    moles = [_si(spec, "moles", value) for value in moles]
    chem_ref = [_si(spec, "chem_ref", value) for value in chem_ref]
    log_activity = [_si(spec, "log_activity", value) for value in log_activity]

    n_elements = len(b)
    species = len(moles)

    if len(a_matrix) != n_elements:
        raise InvalidInputError(
            "a_matrix",
            f"the element matrix has {len(a_matrix)} row(s) against {n_elements} element(s)",
        )
    for row, values in enumerate(a_matrix):
        if len(values) != species:
            raise InvalidInputError(
                "a_matrix",
                f"row {row} has {len(values)} entries against {species} component(s)",
            )
    if len(chem_ref) != species or len(log_activity) != species:
        raise InvalidInputError(
            "chem_ref",
            f"{len(chem_ref)} reference potential(s) and {len(log_activity)} activity "
            f"coefficient(s) against {species} component(s)",
        )

    # Two compositions, because NeqSim has two: `committed` is what the phase holds and
    # `trial` is what this pass computed. A pass whose error did not improve is not
    # written back, so the next pass re-derives from the committed state.
    committed = list(moles)
    trial = list(moles)

    error = 1.0e10
    err_old = 1.0e10
    max_error = tolerance
    iterations = 0
    stagnation = 0
    best_error = math.inf

    while True:
        iterations += 1
        err_old = error
        error = 0.0

        trial = list(committed)
        dn, a_lambda = _chem_solve(a_matrix, b, committed, chem_ref, log_activity)
        step = _step_of(committed, chem_ref, log_activity, dn, a_lambda, temperature)

        for i in range(species):
            if committed[i] < MIN_MOLES:
                continue
            if not math.isfinite(dn[i]):
                error = math.nan
                break
            if abs(dn[i] / committed[i]) > 1e-15:
                error += abs(dn[i] / committed[i])
                trial[i] = dn[i] * step + committed[i]

        if not math.isfinite(error):
            break

        if error < best_error:
            best_error = error
            stagnation = 0
        else:
            stagnation += 1
        if stagnation >= STAGNATION_LIMIT:
            break

        if error <= err_old:
            committed = list(trial)

        continuing = (
            err_old > max_error and abs(error) > max_error and iterations < max_iterations
        ) or iterations < 2
        if not continuing:
            break

        if iterations > 15:
            max_error *= 1.5

    return ChemicalEquilibriumResult(
        moles=tuple(from_si(value, "mol") for value in committed),
        iterations=iterations,
        error=error,
        converged=math.isfinite(error) and error < max_error,
        warnings=tuple(warnings),
    )


def _si(spec: dict[str, object], name: str, value: float | Q) -> float:
    """A declared input as its SI magnitude, quantity or not.

    **The boundary is mixed by design.** A vector whose spec declares a unit arrives as a
    pint quantity and has to be converted; a dimensionless one - ``chem_ref``,
    ``log_activity``, the element matrix - arrives as the bare number it is, and
    ``input_to_si`` refuses that rather than passing it through.
    """
    # A quantity or a bare number. `Q` is a type alias rather than a class, so the check
    # is made the other way round: the numeric branch is the one `isinstance` can name.
    if isinstance(value, int | float):
        return float(value)
    return input_to_si(spec, name, value)


def _chem_solve(
    a_matrix: list[list[float]],
    b: list[float],
    n_mol: list[float],
    chem_ref: list[float],
    log_activity: list[float],
) -> tuple[list[float], list[float]]:
    """One `chemSolve`: the bordered Lagrange system and the Newton direction."""
    n_elements = len(b)
    species = len(n_mol)
    n_t = max(sum(n_mol), MIN_MOLES)

    # `M` is diagonal in the ideal form, so `M^-1 v` is `n_i * v_i` - but it is built and
    # solved rather than shortcut, because the full form is not diagonal and the port
    # would not survive its arrival.
    m = [[0.0] * species for _ in range(species)]
    chem_pot = [0.0] * species
    for i in range(species):
        n_i = max(n_mol[i], MIN_MOLES)
        m[i][i] = 1.0 / n_i
        chem_pot[i] = chem_ref[i] + math.log(n_i) - math.log(n_t) + log_activity[i]

    a_transpose = [[a_matrix[e][i] for e in range(n_elements)] for i in range(species)]
    m_inv_at = _linalg.solve_columns(m, a_transpose)

    ama = [
        [sum(a_matrix[e][i] * m_inv_at[i][f] for i in range(species)) for f in range(n_elements)]
        for e in range(n_elements)
    ]

    mu_column = [[value] for value in chem_pot]
    m_inv_mu = _linalg.solve_columns(m, mu_column)
    amu = [sum(a_matrix[e][i] * m_inv_mu[i][0] for i in range(species)) for e in range(n_elements)]

    nmu = sum(n_mol[i] * chem_pot[i] for i in range(species))

    # The conservation coupling: `b` unless the phase's own element amounts disagree with
    # it, in which case the amounts are what is conserved and `b - A n` is carried.
    coupling = list(b)
    correction = [0.0] * n_elements
    for e in range(n_elements):
        current = sum(a_matrix[e][i] * n_mol[i] for i in range(species))
        if abs(b[e] - current) > CONSERVATION_CORRECTION_TOLERANCE * max(1.0, abs(b[e])):
            coupling = [
                sum(a_matrix[f][i] * n_mol[i] for i in range(species)) for f in range(n_elements)
            ]
            correction = [b[f] - coupling[f] for f in range(n_elements)]
            break

    size = n_elements + 1
    larger = [[0.0] * size for _ in range(size)]
    for e in range(n_elements):
        for f in range(n_elements):
            larger[e][f] = ama[e][f]
        larger[e][n_elements] = coupling[e]
        larger[n_elements][e] = coupling[e]
    rhs = [amu[e] + correction[e] for e in range(n_elements)] + [nmu]

    solved = _linalg.solve_lu(larger, rhs)
    tau = solved[n_elements]

    a_lambda = [0.0] * species
    lambda_rhs = [[0.0] for _ in range(species)]
    for i in range(species):
        total = sum(a_matrix[e][i] * solved[e] for e in range(n_elements))
        a_lambda[i] = total
        lambda_rhs[i][0] = total - chem_pot[i]

    direction = _linalg.solve_columns(m, lambda_rhs)
    dn = [direction[i][0] + n_mol[i] * tau for i in range(species)]
    return dn, a_lambda


def _inner_step(n_mol: list[float], dn: list[float], first: int, species: int) -> float:
    """`innerStep`: the largest step keeping every species non-negative, 3% short of it."""
    shortest = (-n_mol[first] / dn[first]) * 0.97
    for i in range(first, species):
        if n_mol[i] + dn[i] < 0.0:
            candidate = (-n_mol[i] / dn[i]) * 0.97
            if candidate < shortest:
                shortest = candidate
    return 1.0 if shortest > 1.0 else shortest


def _step_of(
    n_mol: list[float],
    chem_ref: list[float],
    log_activity: list[float],
    dn: list[float],
    a_lambda: list[float],
    temperature: float,
) -> float:
    """`step()`: the damped step length, from the two Gibbs measures."""
    species = len(n_mol)
    n_t = max(sum(n_mol), MIN_MOLES)
    r_t = GAS_CONSTANT * temperature

    n_omega = [n_mol[i] + dn[i] for i in range(species)]

    # The negative-moles branch, which is reached on the first pass. See the module
    # docstring for what flooring it costs.
    negative = next((i for i, value in enumerate(n_omega) if value < 0.0), None)
    if negative is not None:
        return _inner_step(n_mol, dn, negative, species)

    def potential(value: float, index: int) -> float:
        return r_t * (
            chem_ref[index] + math.log(max(value, MIN_MOLES)) - math.log(n_t) + log_activity[index]
        )

    g_1 = sum(
        (potential(n_omega[i], i) - a_lambda[i])
        * dn[i]
        * (1.0 / max(n_omega[i], MIN_MOLES) - 1.0 / n_t)
        for i in range(species)
    )

    step = 1.0
    if g_1 > 0.0:
        g_0 = sum(
            (potential(n_mol[i], i) - a_lambda[i])
            * dn[i]
            * (1.0 / max(n_mol[i], MIN_MOLES) - 1.0 / n_t)
            for i in range(species)
        )
        denominator = g_0 - g_1
        if abs(denominator) > 1e-30:
            step = g_0 / denominator

    if not 0.0 <= step <= 1.0 or not math.isfinite(step):
        step = 1.0
    return step


def _spec() -> dict[str, object]:
    """The generated spec, imported here so the module has no import-time cycle."""
    from azoth._models_gen import model

    return model("reactions.chemical_equilibrium")
