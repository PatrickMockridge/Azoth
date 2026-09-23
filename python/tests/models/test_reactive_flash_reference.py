"""The Python twin of the reactive flash stack, against the same capture as the Rust one.

The oracle is ``validation/neqsim/ReactiveFlashProbe.java`` and its capture
``validation/neqsim/captures/reactive_flash_probe.tsv``. These tests exist because the
stack's pieces are private machinery rather than a registered model: the accelerator, the
solve and the driver are all reached through ``reactions.reactive_tp_flash`` from outside,
and pinning them here is what says the two kernels agree before the model lands on top.

**The composition is the pin and the split is not.** Where two phases converge to the same
composition the Gibbs energy is flat along the direction that trades moles between them, so
the split a run reports is decided by its path - the driver's own tests record that NeqSim
and this port land on different points of that line.
"""

from __future__ import annotations

from typing import Any

from azoth.reactions.reference import _diis

#: The capture's ``diis-accelerator`` block: seven scripted entries, a history of four, and
#: the combination at each step. Only the six steps past the first combine.
CAPTURED = {
    1: [1.592_150_170_648_464_9, 3.184_300_341_296_929_7, 4.776_450_511_945_395],
    2: [0.499_999_999_999_399_6, 0.999_999_999_998_799_2, 1.499_999_999_998_198_8],
    3: [0.500_000_000_000_008, 1.000_000_000_000_016, 1.500_000_000_000_021_3],
    4: [0.499_999_999_999_989_8, 0.999_999_999_999_979_6, 1.499_999_999_999_977],
    5: [0.500_000_000_000_852_7, 1.000_000_000_001_705_3, 1.500_000_000_002_572_2],
    6: [0.500_000_000_001_357_1, 1.000_000_000_002_714_3, 1.500_000_000_004_064_3],
}


def _entry(step: int) -> tuple[list[float], list[float]]:
    """The sequence the probe feeds: three components wide, seven entries long."""
    iterate = [0.5 * (step + 1) * (i + 1) for i in range(3)]
    residual = [(i + 1) / (step + 1) + 0.1 * step for i in range(3)]
    return iterate, residual


def test_the_pulay_combination_is_the_classes() -> None:
    """The combination on a buffer that **wraps**, which is what makes ``bufferIndex`` real.

    Two pairs is the class's minimum and the history holds four, so the last three steps
    extrapolate over a rolling window rather than a growing one. Every value is the
    capture's, to ``1e-9``.
    """
    diis = _diis.DiisAccelerator(3, 4)
    assert diis.count() == 0
    assert not diis.can_extrapolate()

    for step in range(7):
        iterate, residual = _entry(step)
        diis.add_entry(iterate, residual)
        assert diis.count() == min(step + 1, 4)
        assert diis.can_extrapolate() is (step > 0)
        extrapolated = diis.extrapolate()
        if step == 0:
            assert extrapolated is None
            continue
        assert extrapolated is not None
        for component, (got, want) in enumerate(zip(extrapolated, CAPTURED[step], strict=True)):
            assert abs(got - want) / abs(want) < 1.0e-9, (step, component, got, want)


def test_identical_residuals_are_refused() -> None:
    """Two equal residuals make every overlap entry the same and the elimination leaves a
    zero pivot: the class says it *can* extrapolate and then returns ``null``."""
    singular = _diis.DiisAccelerator(3, 4)
    singular.add_entry([1.0, 1.0, 1.0], [1.0, 2.0, 3.0])
    singular.add_entry([2.0, 2.0, 2.0], [1.0, 2.0, 3.0])
    assert singular.can_extrapolate()
    assert singular.extrapolate() is None


def test_reset_discards_the_history() -> None:
    """``reset`` empties the buffer, which is what the capture's three lines after the loop
    show."""
    diis = _diis.DiisAccelerator(3, 4)
    for step in (0, 1):
        iterate, residual = _entry(step)
        diis.add_entry(iterate, residual)
    assert diis.extrapolate() is not None
    diis.reset()
    assert diis.count() == 0
    assert not diis.can_extrapolate()
    assert diis.extrapolate() is None


# --- the driver, against the same capture -----------------------------------------------

NAMES = ["CO", "water", "CO2", "hydrogen"]
FEED = [0.25, 0.25, 0.25, 0.25]


def _element_matrix(names: list[str]) -> tuple[list[list[float]], list[str]]:
    """The element matrix in the order the Rust side collects it: first-encountered.

    `FormulaMatrix.build` walks the components and appends an element the first time it
    sees it, so the WGS fluid's rows are `C O H` and not a sorted `C H O`.
    """
    from azoth.reactions.reference import _tables

    element_names: list[str] = []
    per_component: list[list[tuple[str, float]]] = []
    for name in names:
        found = _tables.element_composition(name)
        assert found is not None, name
        per_component.append(list(found))
        for element, _count in found:
            if element not in element_names:
                element_names.append(element)

    matrix = [[0.0] * len(names) for _ in element_names]
    for column, entries in enumerate(per_component):
        for element, count in entries:
            matrix[element_names.index(element)][column] = count
    return matrix, element_names


def _formation(names: list[str]) -> list[Any]:
    """The thermo data `standard_potentials` reads, from the databank and the component row."""
    from azoth.eos.components import entry
    from azoth.reactions.reference import _tables
    from azoth.reactions.reference._rand_solver import ThermoData

    out = []
    for name in names:
        formation = _tables.formation_properties(name)
        assert formation is not None, name
        cp = entry(name).cp
        assert cp is not None, name
        out.append(
            ThermoData(
                enthalpy_of_formation=formation.enthalpy_of_formation,
                absolute_entropy=formation.absolute_entropy,
                gibbs_energy_of_formation=formation.gibbs_energy_of_formation,
                cp=tuple(cp),  # type: ignore[arg-type]
            )
        )
    return out


def _driver(temperature: float, phases: list[Any], max_phases: int = 2) -> Any:
    """One `run` at the state, with SRK and the vapour root for every phase."""
    from azoth.eos.components import from_names
    from azoth.eos.reference._mixture_state import phase_state, reduced_parameters
    from azoth.reactions.reference import _reactive_flash as driver
    from azoth.reactions.reference._rand_solver import solve_single_phase, standard_potentials
    from azoth.reactions.reference._reactive_stability import CriticalConstants

    mixture = from_names(NAMES, eos="srk")
    reduced = reduced_parameters(mixture, temperature, 1.0e5)
    matrix, _names = _element_matrix(NAMES)
    b = [sum(row[i] * FEED[i] for i in range(len(FEED))) for row in matrix]
    g0 = standard_potentials(_formation(NAMES), temperature, 1.0)

    constants = [
        CriticalConstants(
            tc=component.Tc.to("K").magnitude,
            pc=component.Pc.to("bar").magnitude,
            omega=component.omega,
        )
        for component in mixture.components
    ]

    def ln_phi_at(x: list[float]) -> list[float]:
        return list(phase_state(reduced, mixture.kij, x, liquid=False).ln_phi)

    def phase_ln_phi(_index: int, x: list[float]) -> list[float]:
        return ln_phi_at(x)

    def ce(x: list[float]) -> list[float]:
        solved = solve_single_phase(matrix, g0, b, x, ln_phi_at)
        if not solved.converged:
            return list(x)
        total = sum(solved.moles)
        return [moles / total for moles in solved.moles]

    return driver.run(
        driver.DriverState(
            feed_moles=FEED,
            a_matrix=matrix,
            g0=g0,
            b=b,
            total_moles=1.0,
            constants=constants,
            charges=[0.0] * len(NAMES),
            temperature=temperature,
            pressure=1.0,
            max_phases=max_phases,
            phases=phases,
        ),
        phase_ln_phi,
        ln_phi_at,
        ce,
    )


def test_the_driver_reproduces_the_captured_state() -> None:
    """The class's own fluid at 600 K, two phases at the feed: the outer loop runs, the
    composition is the capture's, and the Gibbs measure is **twice** one phase's worth."""
    from azoth.reactions.reference._rand_solver import PhaseFeed

    outcome = _driver(600.0, [PhaseFeed(fractions=list(FEED), beta=1.0)] * 2)
    assert outcome.converged
    assert len(outcome.phases) == 2
    assert outcome.solution is not None
    overall = [sum(phase[i] for phase in outcome.solution.phase_moles) for i in range(len(NAMES))]
    captured = [
        0.079_140_502_548_023_4,
        0.079_140_502_560_693_7,
        0.420_859_061_906_977_14,
        0.420_859_061_880_687_67,
    ]
    for index, (got, want) in enumerate(zip(overall, captured, strict=True)):
        assert abs(got - want) / abs(want) < 1.0e-5, (index, got, want)

    captured_gibbs = -2.259_535_542_715_054_7
    assert abs(outcome.gibbs_energy - captured_gibbs) / abs(captured_gibbs) < 1.0e-5


def test_the_single_phase_branch_reports_no_gibbs_energy() -> None:
    """Forced to one phase the VLE initialisation declines (`V = 1`), the analysis finds the
    fluid stable, and the branch returns **before** the driver computes an energy."""
    from azoth.reactions.reference._rand_solver import PhaseFeed

    outcome = _driver(600.0, [PhaseFeed(fractions=list(FEED), beta=1.0)])
    assert outcome.converged
    assert outcome.gibbs_energy == 0.0


def test_the_vle_initialisation_sets_the_captured_betas() -> None:
    """The 300 K state forced to one phase: the VLE initialisation adds a phase at the
    captured `(0.2106, 0.7894)`, and the driver reports one phase's worth of energy."""
    from azoth.reactions.reference._rand_solver import PhaseFeed

    outcome = _driver(300.0, [PhaseFeed(fractions=list(FEED), beta=1.0)])
    assert outcome.converged
    assert len(outcome.phases) == 2
    assert abs(outcome.phases[0].beta - 0.210_596_317_973_390_96) < 1.0e-12
    assert abs(outcome.phases[1].beta - 0.789_403_682_026_609) < 1.0e-12
    captured_gibbs = -0.716_211_225_797_898_8
    assert abs(outcome.gibbs_energy - captured_gibbs) / abs(captured_gibbs) < 1.0e-4


# --- the PH flash's outer loop ----------------------------------------------------------


def _ideal_gas() -> Any:
    """The mixture's heat-capacity polynomial, from the same component rows `_formation` reads."""
    from azoth.eos import IdealGasModel
    from azoth.eos.components import entry

    coefficients = []
    for name in NAMES:
        cp = entry(name).cp
        assert cp is not None, name
        coefficients.append(cp)
    return IdealGasModel(
        cp_a=tuple(c[0] for c in coefficients),
        cp_b=tuple(c[1] for c in coefficients),
        cp_c=tuple(c[2] for c in coefficients),
        cp_d=tuple(c[3] for c in coefficients),
        cp_e=tuple(c[4] for c in coefficients),
    )


def _thermochemical(temperature: float, outcome: Any) -> float:
    """The state's thermochemical enthalpy: the phases' sensible enthalpies plus the formation
    inventory, which is what the class's own residual compares."""
    import azoth
    from azoth.eos import molar_enthalpy_entropy
    from azoth.eos.components import from_names
    from azoth.eos.reference._mixture_state import phase_state, reduced_parameters
    from azoth.reactions.reference import _tables

    mixture = from_names(NAMES, eos="srk")
    reduced = reduced_parameters(mixture, temperature, 1.0e5)
    ideal = _ideal_gas()
    total = 0.0
    inventory = 0.0
    assert outcome.solution is not None
    for row in outcome.solution.phase_moles:
        held = sum(row)
        if held <= 0.0:
            continue
        composition = [moles / held for moles in row]
        z = phase_state(reduced, mixture.kij, composition, liquid=False).z
        state = molar_enthalpy_entropy(
            mixture,
            ideal,
            azoth.ureg.Quantity(temperature, "K"),
            azoth.ureg.Quantity(1.0, "bar"),
            composition,
            z,
        )
        total += held * float(state.h.to("J/mol").magnitude)
        for index, moles in enumerate(row):
            formation = _tables.formation_properties(NAMES[index])
            assert formation is not None
            inventory += moles * formation.enthalpy_of_formation
    return total + inventory


def _heat_capacity(temperature: float, outcome: Any) -> float:
    """The system's heat capacity, `sum_j n_j cp_j` over the phases."""
    import azoth
    from azoth.eos import molar_enthalpy_entropy
    from azoth.eos.components import from_names
    from azoth.eos.reference._mixture_state import phase_state, reduced_parameters

    mixture = from_names(NAMES, eos="srk")
    reduced = reduced_parameters(mixture, temperature, 1.0e5)
    ideal = _ideal_gas()
    total = 0.0
    assert outcome.solution is not None
    for row in outcome.solution.phase_moles:
        held = sum(row)
        if held <= 0.0:
            continue
        composition = [moles / held for moles in row]
        z = phase_state(reduced, mixture.kij, composition, liquid=False).z
        state = molar_enthalpy_entropy(
            mixture,
            ideal,
            azoth.ureg.Quantity(temperature, "K"),
            azoth.ureg.Quantity(1.0, "bar"),
            composition,
            z,
        )
        total += held * float(state.cp.to("J/K/mol").magnitude)
    return total


def test_the_ph_loop_finds_its_temperature_back() -> None:
    """The class's own round trip, in Python: the WGS equilibrium at 600 K sets the
    specification, the search starts from 500 K and has to come back.

    **The recovered temperature is the capture's to `1e-5` relative** (`600.0000398 K`), and
    the outer pass count is *reported* rather than pinned: the loop is a secant, so its step
    count is a path quantity.
    """
    from azoth.reactions.reference._reactive_ph_flash import PhState, reactive_ph_flash

    def round_trip(flash_temperature: float, perturbed: float) -> tuple[float, int]:
        from azoth.reactions.reference._rand_solver import PhaseFeed

        # The specification: the reactive equilibrium at the flash temperature.
        reference = _driver(flash_temperature, [PhaseFeed(fractions=list(FEED), beta=1.0)] * 2)
        specified = _thermochemical(flash_temperature, reference)

        def inner(temperature: float) -> PhState:
            outcome = _driver(temperature, [PhaseFeed(fractions=list(FEED), beta=1.0)] * 2)
            return PhState(
                iterations=outcome.total_iterations,
                thermochemical_enthalpy=_thermochemical(temperature, outcome),
                cp=_heat_capacity(temperature, outcome),
            )

        result = reactive_ph_flash(perturbed, specified, inner)
        assert result.converged
        return result.temperature, result.outer_iterations

    recovered, _outer = round_trip(600.0, 500.0)
    captured = 600.000_039_758_328_5
    assert abs(recovered - captured) / captured < 1.0e-5, (recovered, captured)


def test_a_trace_ion_fluid_is_answered_by_the_split_it_has() -> None:
    """**The trace-ion short circuit**: a multiphase fluid whose ions are all traces is
    answered by the split it already has, and no solve runs at all.

    The oracle is ``TraceIonProbe`` and ``captures/trace_ion_probe.tsv``: two ``SystemSrkEos``
    fluids differing only in how much salt they carry, where ``1e-12`` mol takes the circuit
    (``total_iterations = 0``, the split untouched) and a mole does not (69 iterations).

    **The sum is over the ions' overall mole *fractions*, not over their charges**, which is why
    the trace fluid here carries a trace *amount* of the charged component.
    """
    from azoth.eos.components import from_names
    from azoth.eos.reference._mixture_state import phase_state, reduced_parameters
    from azoth.reactions.reference import _reactive_flash as driver
    from azoth.reactions.reference._rand_solver import (
        PhaseFeed,
        solve_single_phase,
        standard_potentials,
    )
    from azoth.reactions.reference._reactive_stability import CriticalConstants

    mixture = from_names(NAMES, eos="srk")
    reduced = reduced_parameters(mixture, 600.0, 1.0e5)
    matrix, _names = _element_matrix(NAMES)
    g0 = standard_potentials(_formation(NAMES), 600.0, 1.0)
    constants = [
        CriticalConstants(
            tc=component.Tc.to("K").magnitude,
            pc=component.Pc.to("bar").magnitude,
            omega=component.omega,
        )
        for component in mixture.components
    ]

    def ln_phi_at(x: list[float]) -> list[float]:
        return list(phase_state(reduced, mixture.kij, x, liquid=False).ln_phi)

    def phase_ln_phi(_index: int, x: list[float]) -> list[float]:
        return ln_phi_at(x)

    def one(feed: list[float], charges: list[float]) -> Any:
        b = [sum(row[i] * feed[i] for i in range(len(feed))) for row in matrix]

        def ce(x: list[float]) -> list[float]:
            solved = solve_single_phase(matrix, g0, b, x, ln_phi_at)
            if not solved.converged:
                return list(x)
            total = sum(solved.moles)
            return [moles / total for moles in solved.moles]

        phases = [
            PhaseFeed(fractions=list(feed), beta=1.0),
            PhaseFeed(fractions=list(feed), beta=1.0),
        ]
        return driver.run(
            driver.DriverState(
                feed_moles=feed,
                a_matrix=matrix,
                g0=g0,
                b=b,
                total_moles=1.0,
                constants=constants,
                charges=charges,
                temperature=600.0,
                pressure=1.0,
                max_phases=2,
                phases=phases,
            ),
            phase_ln_phi,
            ln_phi_at,
            ce,
        )

    trace_feed = [0.25, 0.25, 0.25, 1.0e-12]
    trace_charges = [0.0, 0.0, 0.0, 1.0]
    ion_fraction = sum(
        abs(fraction)
        for fraction, charge in zip(trace_feed, trace_charges, strict=True)
        if charge != 0.0
    )
    assert ion_fraction < driver.TRACE_ION_Z_TOLERANCE
    # The fractions are the feed's own and not normalised by the driver, so `0.25` here is the
    # component's share of the *feed* - which for the trace fluid sums to 1.000000000001.
    traced = one([0.25, 0.25, 0.25, 1.0e-12], trace_charges)
    assert traced.converged
    assert traced.total_iterations == 0, "no solve runs"
    assert traced.solution is None, "and none is reported"
    assert traced.phases[0].fractions == trace_feed, "the split is untouched"
    assert traced.phases[1].fractions == trace_feed

    # A quarter mole of the same charged component is not a trace, and the driver runs.
    molar = one(FEED, [0.0, 0.0, 0.0, 1.0])
    assert molar.total_iterations > 0, "a molar ion takes the solve"
