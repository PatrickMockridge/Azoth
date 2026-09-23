"""``eos.hybrid_eos_ge_flash`` - NeqSim's fixed-role EoS-gas / EoS-oil / GE-aqueous flash.

Spec: ``specs/models/eos/hybrid_eos_ge_flash.toml``, which records the fixed topology, the
two rules that make an ion aqueous and the acceptance contract the answer is held to.

This is the pure-Python reference: a second, independent expression of the same physics as
the Rust kernel, held to it by ``test_cross_impl.py``. The two share no solver - each
carries its own dense elimination, which is what makes their agreement a statement about
the arithmetic rather than about a copied loop.

The three roles are **fixed before the fractions are solved**, so this is a fraction Newton
over a topology the caller already knows and not a stability analysis: no phase is added or
removed, and a fluid whose phases are not already known cannot be answered by it.
"""

from __future__ import annotations

import math
from typing import Any

from azoth import _models_gen
from azoth.core.errors import (
    InvalidInputError,
    SolverNotConvergedError,
)
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import HybridEosGeFlashResult
from azoth.core.units import Q, input_to_si
from azoth.core.warnings import Warning
from azoth.eos import components as _components
from azoth.eos.reference._mixture_state import phase_state, reduced_parameters
from azoth.eos.reference.pitzer_phase import pitzer_phase

MODEL_ID = "eos.hybrid_eos_ge_flash"

#: A phase fraction is kept strictly inside ``(0, 1)``, upstream's ``phaseFractionMinimumLimit``.
FRACTION_FLOOR = 1.0e-12

#: ``E_i``'s floor, upstream's.
MINIMUM_E = 1.0e-100

#: The ion concentration a non-aqueous role is given instead of a computed one.
ION_CONFINEMENT = 1.0e-50

#: The cap on the beta step where the inventory carries an ion.
IONIC_BETA_STEP_SCALE = 0.1

#: The ratio a phase fraction may move by in one projection.
AQUEOUS_STEP_LIMIT = 2.0

#: The aqueous fraction's margin over the ionic inventory.
ION_CAPACITY_MARGIN = 100.0 * FRACTION_FLOOR

#: The diagonal regulariser, ``TPHybridEosGeFlash.calcQ``'s own.
REGULARISER = 1.0e-10

#: The gradient norm the iteration must also reach.
GRADIENT_TOLERANCE = 1.0e-10

#: A fraction at or below this is a phase that is not there, in the acceptance contract.
MATERIAL_PHASE_FRACTION = 1.0e-10

#: The fixed-topology solve's own cap of passes.
MAXIMUM_HYBRID_ITERATIONS = 50

#: The fewest outer passes, upstream's ``iteration >= 4``.
MINIMUM_OUTER_ITERATIONS = 4

#: The roles, in the solver's own order, which is what every vector here is indexed by.
ROLE_ORDER = ("gas", "oil", "aqueous")

#: The polar solvents `SystemEosGE.isAqueousComponent` lists beside water.
AQUEOUS_SOLVENTS = frozenset({"water", "meg", "teg", "deg", "methanol", "ethanol"})

#: The component class the two EoS roles' seeding reads a hydrocarbon from.
HYDROCARBON = "hc"


class _Phase:
    """One role's fraction and composition."""

    __slots__ = ("composition", "fraction")

    def __init__(self, fraction: float, composition: list[float]) -> None:
        self.fraction = fraction
        self.composition = composition


def _solve(hessian: list[list[float]], gradient: list[float]) -> list[float] | None:
    """``H x = g`` by Gaussian elimination with partial pivoting, or ``None`` if singular.

    The pivot takes the **last** largest magnitude, which is what Rust's ``max_by`` does:
    two kernels that break a tie differently take different elimination paths and a case
    that pins a step count would disagree for a reason that is not the physics.
    """
    n = len(gradient)
    a = [[*hessian[i][:], gradient[i]] for i in range(n)]
    for column in range(n):
        pivot = column
        for row in range(column + 1, n):
            if abs(a[row][column]) >= abs(a[pivot][column]):
                pivot = row
        a[column], a[pivot] = a[pivot], a[column]
        if a[column][column] == 0.0:
            return None
        pivot_value = a[column][column]
        for row in range(column + 1, n):
            factor = a[row][column] / pivot_value
            for entry in range(column, n + 1):
                a[row][entry] -= factor * a[column][entry]
    x = [0.0] * n
    for row in range(n - 1, -1, -1):
        total = sum(a[row][k] * x[k] for k in range(row + 1, n))
        x[row] = (a[row][n] - total) / a[row][row]
    return x


def _resolved(names: list[str], card: Any = None) -> list[Any]:
    """The databank entries, in the feed's order."""
    return [_components.entry(name, card=card) for name in names]


def _is_ion(entry: Any) -> bool:
    """Whether the component may not leave the aqueous role.

    The databank's own classification and not a second input: ``COMPTYPE == ion`` is the
    compiled table's marker, and NeqSim tests the charge and the flag beside it.
    """
    return bool(entry.component_type == "ion" or entry.ionic_charge != 0.0)


def _seed(entries: list[Any], z: list[float], ion: list[bool]) -> list[_Phase]:
    """`SystemEosGE.prepareHybridEosGeFlash`'s per-class seed.

    A component's *class* and not its state: an ion, water and the polar solvents start in
    the brine, a heavy hydrocarbon in the oil and a light one in the gas, each at the feed's
    own share of its group floored at ``1e-5``.
    """
    n = len(z)
    gas = [0.0] * n
    oil = [0.0] * n
    aqueous = [0.0] * n
    gas_feed = oil_feed = aqueous_feed = 0.0

    for i in range(n):
        entry = entries[i]
        fraction = max(z[i], ION_CONFINEMENT)
        heavy = (entry.molar_mass is not None) and (entry.molar_mass.to("kg/mol").magnitude > 0.045)
        if ion[i]:
            gas[i] = ION_CONFINEMENT
            oil[i] = ION_CONFINEMENT
            aqueous[i] = fraction
            aqueous_feed += max(z[i], 0.0)
        elif entry.name.strip().lower() in AQUEOUS_SOLVENTS:
            gas[i] = min(fraction * 1.0e-3, 1.0e-8)
            oil[i] = min(fraction * 1.0e-4, 1.0e-10)
            aqueous[i] = fraction
            aqueous_feed += max(z[i], 0.0)
        elif entry.component_type == HYDROCARBON:
            if heavy:
                gas[i] = max(fraction * 1.0e-2, 1.0e-16)
                oil[i] = fraction
                oil_feed += max(z[i], 0.0)
            else:
                gas[i] = fraction
                oil[i] = max(fraction * 5.0e-2, 1.0e-16)
                gas_feed += max(z[i], 0.0)
            aqueous[i] = max(fraction * 1.0e-8, 1.0e-30)
        else:
            gas[i] = fraction
            oil[i] = max(fraction * 1.0e-1, 1.0e-16)
            aqueous[i] = max(fraction * 1.0e-2, 1.0e-20)
            gas_feed += max(z[i], 0.0)

    for values in (gas, oil, aqueous):
        total = sum(values)
        if total > 0.0:
            for index in range(n):
                values[index] /= total

    floor = 1.0e-5
    gas_feed = max(gas_feed, floor)
    oil_feed = max(oil_feed, floor)
    aqueous_feed = max(aqueous_feed, floor)
    total = gas_feed + oil_feed + aqueous_feed
    return [
        _Phase(gas_feed / total, gas),
        _Phase(oil_feed / total, oil),
        _Phase(aqueous_feed / total, aqueous),
    ]


def _enforce_aqueous_fraction_bounds(
    phases: list[_Phase], ion: list[bool], z: list[float], aqueous: int, previous: float
) -> None:
    """`enforceAqueousFractionBounds`: the aqueous fraction projected onto its band.

    The band's floor is the ionic inventory plus a margin - a brine smaller than the ions it
    must hold is not a state - and both ends are held within a factor of
    ``AQUEOUS_STEP_LIMIT`` of the fraction the correction started from.

    Raises:
        SolverNotConvergedError: if the other roles cannot give up enough fraction to reach
            the floor. NeqSim throws the same state.
    """
    inventory = sum(max(z[i], 0.0) for i in range(len(z)) if ion[i])
    if not (inventory > 0.0) or not math.isfinite(inventory):
        return

    maximum = 1.0 - 2.0 * FRACTION_FLOOR
    stepped_minimum = (
        previous / AQUEOUS_STEP_LIMIT if previous > 0.0 else inventory + ION_CAPACITY_MARGIN
    )
    minimum = min(maximum, max(inventory + ION_CAPACITY_MARGIN, stepped_minimum))
    stepped_maximum = min(previous * AQUEOUS_STEP_LIMIT, maximum)
    current = phases[aqueous].fraction
    if minimum <= current <= stepped_maximum:
        return

    others = [k for k in range(3) if k != aqueous]
    if current > stepped_maximum:
        excess = current - stepped_maximum
        share_total = sum(phases[k].fraction for k in others)
        for k in others:
            share = phases[k].fraction / share_total if share_total > 0.0 else 0.5
            phases[k].fraction += excess * share
        phases[aqueous].fraction = stepped_maximum
        return

    required = minimum - current
    adjustable = sum(max(phases[k].fraction - FRACTION_FLOOR, 0.0) for k in others)
    if adjustable + FRACTION_FLOOR < required:
        raise SolverNotConvergedError(0, inventory, minimum)
    remaining = required
    last = None
    for k in others:
        available = max(phases[k].fraction - FRACTION_FLOOR, 0.0)
        if available <= 0.0:
            continue
        last = k
        transfer = min(required * available / adjustable, remaining)
        phases[k].fraction -= transfer
        remaining -= transfer
    if remaining > 0.0 and last is not None:
        phases[last].fraction -= remaining
    non_aqueous = sum(phases[k].fraction for k in others)
    phases[aqueous].fraction = 1.0 - non_aqueous


def _set_compositions(
    phases: list[_Phase],
    z: list[float],
    e: list[float],
    inverted: list[list[float]],
    ion: list[bool],
    aqueous: int,
) -> None:
    """`TPHybridEosGeFlash.setXY`: the compositions the material balance gives.

    Three rules and not one. An **ion** is ``z_i / beta`` in the brine and
    ``ION_CONFINEMENT`` in the other two roles. A **neutral** is ``z_i / (E_i phi_ik)``,
    floored. And the brine's neutrals are then scaled so that they and the ions together sum
    to one, which is what makes the ions' share exact rather than diluted.

    Raises:
        SolverNotConvergedError: if the brine cannot hold the ionic inventory at this
            fraction. NeqSim throws the same state.
    """
    n = len(z)
    for k, phase in enumerate(phases):
        fraction = max(phase.fraction, FRACTION_FLOOR)
        ion_sum = 0.0
        neutral_sum = 0.0
        row = [0.0] * n
        for i in range(n):
            if ion[i]:
                value = z[i] / fraction if k == aqueous else ION_CONFINEMENT
            else:
                value = z[i] / e[i] * inverted[k][i]
            if not math.isfinite(value) or value <= 0.0:
                value = ION_CONFINEMENT
            row[i] = value
            if ion[i]:
                ion_sum += value
            else:
                neutral_sum += value
        if k == aqueous:
            if not (ion_sum < 1.0) or not (neutral_sum > 0.0):
                raise SolverNotConvergedError(0, ion_sum, 1.0)
            neutral_total = 1.0 - ion_sum
            for i in range(n):
                if not ion[i]:
                    row[i] = neutral_total * row[i] / neutral_sum
        else:
            row_total = sum(row)
            if row_total > 0.0:
                row = [value / row_total for value in row]
        phase.composition = row


def _acceptance(
    beta: list[float],
    composition: list[list[float]],
    ln_phi: list[list[float]],
    z: list[float],
    ion: list[bool],
    pressure: float,
) -> tuple[float, float]:
    """The acceptance contract's two residuals.

    The balance is ``z_i - sum_k beta_k x_ik``; the equilibrium is the spread of
    ``ln(x_i phi_i P)`` over the roles that are there, taken over the components the contract
    admits - an empty component or an ion takes no part.
    """
    n = len(z)
    balance = 0.0
    for i in range(n):
        split = sum(beta[k] * composition[k][i] for k in range(len(beta)))
        balance = max(balance, abs(z[i] - split))

    fugacity = 0.0
    for i in range(n):
        if z[i] <= 1.0e-30 or ion[i]:
            continue
        first: float | None = None
        for k in range(len(beta)):
            if beta[k] <= MATERIAL_PHASE_FRACTION:
                continue
            value = composition[k][i] * math.exp(ln_phi[k][i]) * pressure
            if not (value > 0.0) or not math.isfinite(value):
                continue
            logged = math.log(value)
            if first is None:
                first = logged
            else:
                fugacity = max(fugacity, abs(first - logged))
    if not math.isfinite(balance) or not math.isfinite(fugacity):
        raise SolverNotConvergedError(0, max(balance, fugacity), 0.0)
    return balance, fugacity


def hybrid_eos_ge_flash(
    components: list[str],
    cubic: str,
    T: Q,
    P: Q,
    moles: list[float],
) -> HybridEosGeFlashResult:
    """The isothermal flash of a fixed gas-oil-brine topology.

    Args:
        components: the substances the feed is made of, by name.
        cubic: the equation of state **both** EoS roles are.
        T: absolute temperature.
        P: absolute pressure.
        moles: the feed's mole numbers, one per component. Moles and not mole fractions:
            the material balance and the ionic inventory are both mole sums.

    Returns:
        The fraction and composition of each of ``[gas, oil, aqueous]``, the coefficients
        they were solved against, and the two residuals the acceptance contract is stated
        in.

    Raises:
        InvalidInputError: if a name is unknown, or a component has no heat-capacity
            coefficients.
        OutOfRangeError: if ``T`` or ``P`` is not positive.
        SolverNotConvergedError: if a Newton correction cannot be solved at an iterate, or
            the brine cannot hold the ionic inventory.

    Example:
        >>> import azoth
        >>> q = azoth.ureg.Quantity
        >>> r = hybrid_eos_ge_flash(
        ...     ["methane", "n-heptane", "water", "Na+", "Cl-"],
        ...     "srk",
        ...     q(313.15, "K"),
        ...     q(50.0e5, "Pa"),
        ...     [5.0, 2.0, 55.5, 1.0, 1.0],
        ... )
        >>> [round(beta, 6) for beta in r.beta]
        [0.070088, 0.038593, 0.89132]
    """
    spec = _models_gen.model(MODEL_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    t_si = input_to_si(spec, "T", T)
    p_si = input_to_si(spec, "P", P)
    apply_checks(checks.on_input, {"T": t_si, "P": p_si}.get, warnings)

    # **The boundary is mixed by design.** `moles` carries a unit, so the case hands it over
    # as pint quantities and doing arithmetic on one against a bare float raises rather than
    # converts. `T` and `P` are scalars and go through `input_to_si` below.
    moles = [_si(spec, "moles", value) for value in moles]

    n = len(components)
    if len(moles) != n:
        raise InvalidInputError("moles", f"a feed for {n} components has {len(moles)} entries")
    if any(value < 0.0 for value in moles):
        raise InvalidInputError("moles", "a mole number cannot be negative")
    total = sum(moles)
    if not (total > 0.0) or not math.isfinite(total):
        raise InvalidInputError("moles", f"the feed sums to {total}, which is not a mole count")
    z = [value / total for value in moles]

    entries = _resolved(components)
    ion = [_is_ion(entry) for entry in entries]

    mixture, _ = _components.mixture_of_with_ions(components, eos=cubic)
    reduced = reduced_parameters(mixture, t_si, p_si)
    warnings.extend(reduced.warnings)
    kij = mixture.kij

    algorithm = spec["algorithm"]
    tolerance = float(algorithm["tolerance"])
    cap = int(algorithm["max_iterations"])

    phases = _seed(entries, z, ion)
    aqueous = 2
    iterations, residual, gradient_norm = solve_fixed_topology_from(
        reduced,
        kij,
        phases,
        moles,
        ion,
        components,
        t_si,
        p_si,
        tolerance,
        cap,
    )

    beta = [phase.fraction for phase in phases]
    composition = [list(phase.composition) for phase in phases]
    ln_phi = [
        [math.log(value) for value in row]
        for row in _coefficients(reduced, kij, phases, components, t_si, p_si, aqueous)
    ]

    balance, fugacity = _acceptance(beta, composition, ln_phi, z, ion, p_si)
    min_t_over_tc = min(reduced.reduced_temperatures)

    return HybridEosGeFlashResult(
        beta=tuple(beta),
        x=tuple(tuple(row) for row in composition),
        ln_phi=tuple(tuple(row) for row in ln_phi),
        iterations=iterations,
        residual=max(residual, gradient_norm),
        max_material_balance_residual=balance,
        max_log_fugacity_residual=fugacity,
        min_t_over_tc=min_t_over_tc,
        warnings=tuple(warnings),
    )


def solve_fixed_topology_from(
    reduced: Any,
    kij: Any,
    phases: list[_Phase],
    feed: list[float],
    ion: list[bool],
    components: list[str],
    t_si: float,
    p_si: float,
    tolerance: float,
    cap: int,
) -> tuple[int, float, float]:
    """The fixed-topology solve **from a split the caller already has**.

    Split out of :func:`hybrid_eos_ge_flash` because the reactive coupling re-enters it: that
    loop re-equilibrates the brine's chemistry, writes the adjusted species inventory back into
    the roles and solves the fractions again, and each pass starts from the fractions the last
    one converged on rather than from the seed. NeqSim's ``run`` is arranged the same way - its
    ``system`` carries the split between passes - and re-seeding each pass would throw away the
    only thing the passes have in common.

    ``feed`` is the **reaction-adjusted overall inventory and not the feed constant**: the
    chemistry changes species amounts, and what the material balance is stated over moves with
    them. ``phases`` is in and out, in ``[gas, oil, aqueous]`` order.

    Returns the steps taken, the last step's norm and the last gradient's, which is what the
    caller's own convergence test reads.
    """
    aqueous = 2
    n = len(feed)
    z = [value / sum(feed) for value in feed]
    iterations = 0
    residual = math.nan
    gradient_norm = math.inf

    for outer in range(MAXIMUM_HYBRID_ITERATIONS):
        for step in range(1, cap + 1):
            iterations += 1
            phi = _coefficients(reduced, kij, phases, components, t_si, p_si, aqueous)
            # An ion is excluded from the EoS roles **exactly**: its inverse coefficient is
            # zero there, so it never enters a gradient, and a coefficient the model could
            # not answer with is excluded the same way.
            inverted = [
                [
                    0.0 if (ion[i] and k != aqueous) or not (phi[k][i] > 0.0) else 1.0 / phi[k][i]
                    for i in range(n)
                ]
                for k in range(3)
            ]

            e = [0.0] * n
            for k, phase in enumerate(phases):
                for i in range(n):
                    e[i] += phase.fraction * inverted[k][i]
            for i in range(n):
                if not math.isfinite(e[i]) or e[i] < MINIMUM_E:
                    e[i] = MINIMUM_E

            gradient = [1.0] * 3
            for k in range(3):
                for i in range(n):
                    gradient[k] -= z[i] * inverted[k][i] / e[i]
            hessian = [[0.0] * 3 for _ in range(3)]
            for j in range(3):
                for k in range(3):
                    for i in range(n):
                        hessian[j][k] += z[i] * inverted[j][i] * inverted[k][i] / (e[i] * e[i])
                    if j == k:
                        hessian[j][k] += REGULARISER

            gradient_norm = math.sqrt(sum(value * value for value in gradient))
            correction = _solve(hessian, gradient)
            if correction is None:
                raise SolverNotConvergedError(iterations, math.nan, tolerance)
            residual = math.sqrt(sum(value * value for value in correction))

            # The step is damped by `n/(n+3)`, capped at a tenth while the inventory carries
            # an ion: an unconstrained proposal can push the brine below its own ions.
            scale = step / (step + 3.0)
            if any(ion):
                scale = min(scale, IONIC_BETA_STEP_SCALE)
            # `calcQ` captures the settled fraction *before* the proposal is written.
            previous_aqueous = phases[aqueous].fraction
            total_fraction = 0.0
            for k, phase in enumerate(phases):
                candidate = phase.fraction - scale * correction[k]
                phase.fraction = min(max(candidate, FRACTION_FLOOR), 1.0 - FRACTION_FLOOR)
                total_fraction += phase.fraction
            for phase in phases:
                phase.fraction /= total_fraction

            _enforce_aqueous_fraction_bounds(phases, ion, z, aqueous, previous_aqueous)
            _set_compositions(phases, z, e, inverted, ion, aqueous)

            if step >= 2 and residual <= tolerance and gradient_norm <= GRADIENT_TOLERANCE:
                break
        if (
            outer >= MINIMUM_OUTER_ITERATIONS
            and residual <= tolerance
            and gradient_norm <= GRADIENT_TOLERANCE
        ):
            break
    return iterations, residual, gradient_norm


def _coefficients(
    reduced: Any,
    kij: Any,
    phases: list[_Phase],
    components: list[str],
    t_si: float,
    p_si: float,
    aqueous: int,
) -> list[list[float]]:
    """The fugacity coefficients of every component in every role.

    The two EoS roles are the cubic's, at each role's own composition and on its own root;
    the brine's are ``eos.pitzer_phase``'s, which is the model whose *activity* surface the
    phase is and whose fugacity surface this flash is the first caller of.
    """
    out: list[list[float]] = []
    for k, phase in enumerate(phases):
        if k == aqueous:
            result = pitzer_phase(components, _kelvin(t_si), _pascal(p_si), list(phase.composition))
            out.append([math.exp(value) for value in result.ln_phi])
        else:
            state = phase_state(reduced, kij, list(phase.composition), liquid=k != 0)
            out.append([math.exp(value) for value in state.ln_phi])
    if any(not (value > 0.0) or not math.isfinite(value) for row in out for value in row):
        # A coefficient a role cannot answer with is not a state this solve can carry, and
        # the caller's reading of it is an equilibrium quantity. The Rust kernel excludes it
        # the same way; both then leave the material balance to say so.
        pass
    return out


def _si(spec: Any, name: str, value: float | Q) -> float:
    """A declared input as its SI magnitude, quantity or not."""
    if isinstance(value, int | float):
        return float(value)
    return input_to_si(spec, name, value)


def _kelvin(value: float) -> Q:
    """A bare kelvin quantity, for the phase model's own argument."""
    import azoth

    return azoth.ureg.Quantity(value, "K")


def _pascal(value: float) -> Q:
    """A bare pascal quantity, for the phase model's own argument."""
    import azoth

    return azoth.ureg.Quantity(value, "Pa")
