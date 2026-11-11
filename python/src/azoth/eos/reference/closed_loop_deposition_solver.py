"""``eos.closed_loop_deposition_solver`` - NeqSim's ``ClosedLoopDepositionSolver``.

The pipe is re-run at a shrinking effective diameter and the accumulator evaluated at each
run; the worst segment's thickness sets the next diameter until the diameter stops moving,
the iteration cap is reached, or the line blocks.

Spec: ``specs/models/eos/closed_loop_deposition_solver.toml``.
"""

from __future__ import annotations

from typing import Any

from azoth.core.range import apply_checks, checks_for
from azoth.core.result import ClosedLoopDepositionSolverResult
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning
from azoth.eos.reference.pipe_beggs_and_brills import pipe_beggs_and_brills
from azoth.eos.reference.scale_deposition_accumulator import scale_deposition_accumulator

CALC_ID = "eos.closed_loop_deposition_solver"

MINIMUM_TOLERANCE = 1.0e-7


def closed_loop_deposition_solver(
    components: list[str],
    inlet_n: Q,
    inlet_z: list[float],
    inlet_p: Q,
    inlet_t: Q,
    length: Q,
    diameter: Q,
    elevation: Q,
    number_of_increments: int,
    roughness: Q,
    calcium: Q,
    bicarbonate: Q,
    sulphate: Q,
    barium: Q,
    sodium: Q,
    total_dissolved_solids: Q,
    ph: float,
    service_time: Q,
    scale_density: Q,
    base_rate: Q,
    inhibitor_efficiency: float,
    tolerance: Q,
    max_iterations: int,
    eos: str | None = None,
) -> ClosedLoopDepositionSolverResult:
    """The closed-loop deposited-diameter solve over a Beggs-and-Brills line.

    Args:
        components: the fluid's substances, by name.
        inlet_n: the inlet's molar flow.
        inlet_z: the inlet composition.
        inlet_p: the inlet pressure.
        inlet_t: the inlet temperature, which the isothermal march holds.
        length: the pipe's total length.
        diameter: the pipe's original internal diameter.
        elevation: the total elevation change, positive uphill.
        number_of_increments: the number of pipe segments.
        roughness: the pipe wall's absolute roughness.
        calcium: the calcium ion's concentration.
        bicarbonate: the bicarbonate ion's concentration.
        sulphate: the sulphate ion's concentration.
        barium: the barium ion's concentration.
        sodium: the sodium ion's concentration.
        total_dissolved_solids: TDS.
        ph: the aqueous phase's pH.
        service_time: the exposure time the deposit masses are scaled by.
        scale_density: the assumed deposit density.
        base_rate: the deposition rate at SI = 1.
        inhibitor_efficiency: the uniform fractional efficiency, 0..1.
        tolerance: the convergence tolerance on the effective diameter.
        max_iterations: the iteration cap.
        eos: the cubic the pipe's flash runs; absent takes ``srk``.

    Returns:
        The iteration count, convergence, the final diameter and the three histories.

    Raises:
        OutOfRangeError: if an input is outside the spec's hard range.
        InvalidInputError: if the pipe refuses a segment's negative outlet pressure.

    Example:
        >>> import azoth
        >>> q = azoth.ureg.Quantity
        >>> r = closed_loop_deposition_solver(
        ...     ["methane", "CO2", "water"], q(80.0, "mol/s"), [0.85, 0.05, 0.10],
        ...     q(1.0e7, "Pa"), q(333.15, "K"), q(1000.0, "m"), q(0.1, "m"),
        ...     q(0.0, "m"), 10, q(1.0e-5, "m"),
        ...     q(1.0, "kg/m**3"), q(0.5, "kg/m**3"), q(0.1, "kg/m**3"),
        ...     q(0.01, "kg/m**3"), q(20.0, "kg/m**3"), q(50.0, "kg/m**3"),
        ...     6.5, q(31_557_600.0, "s"), q(2700.0, "kg/m**3"),
        ...     q(1.5844043908606944e-8, "kg/(m**2*s)"), 0.0,
        ...     q(1.0e-4, "m"), 10, "srk")
        >>> r.iterations_taken
        1
    """
    from azoth._models_gen import model as _model_for

    spec: dict[str, Any] = _model_for(CALC_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    values: dict[str, float] = {
        name: input_to_si(spec, name, value)
        for name, value in {
            "inlet_n": inlet_n,
            "inlet_p": inlet_p,
            "inlet_t": inlet_t,
            "length": length,
            "diameter": diameter,
            "elevation": elevation,
            "roughness": roughness,
            "calcium": calcium,
            "bicarbonate": bicarbonate,
            "sulphate": sulphate,
            "barium": barium,
            "sodium": sodium,
            "total_dissolved_solids": total_dissolved_solids,
            "service_time": service_time,
            "scale_density": scale_density,
            "base_rate": base_rate,
            "tolerance": tolerance,
        }.items()
    }
    values["ph"] = ph
    values["inhibitor_efficiency"] = inhibitor_efficiency
    apply_checks(checks.on_input, values.get, warnings)

    tolerance_m = max(values["tolerance"], MINIMUM_TOLERANCE)
    cap = max(int(max_iterations), 1)
    original = values["diameter"]

    current = original
    previous = original
    converged = False
    iterations_taken = 0
    final_diameter = 0.0
    diameter_history: list[float] = []
    thickness_history: list[float] = []
    velocity_history: list[float] = []

    for index in range(cap):
        iterations_taken = index + 1
        run = pipe_beggs_and_brills(
            components,
            from_si(values["inlet_n"], "mol/s"),
            list(inlet_z),
            inlet_p,
            inlet_t,
            length,
            from_si(current, "m"),
            elevation,
            number_of_increments,
            roughness,
            eos,
        )
        deposition = scale_deposition_accumulator(
            calcium,
            bicarbonate,
            sulphate,
            barium,
            sodium,
            total_dissolved_solids,
            ph,
            from_si(current, "m"),
            length,
            list(run.temperature_profile),
            list(run.pressure_profile),
            service_time,
            scale_density,
            base_rate,
            inhibitor_efficiency,
        )

        thickness = deposition.max_thickness.to("m").magnitude
        diameter_history.append(current)
        thickness_history.append(thickness)
        velocities = run.mixture_superficial_velocity_profile
        velocity_history.append(
            0.0
            if not velocities
            else sum(v.to("m/s").magnitude for v in velocities) / len(velocities)
        )

        new_diameter = original - 2.0 * thickness
        if new_diameter <= 0.0:
            final_diameter = 0.0
            converged = True
            break
        if abs(new_diameter - previous) < tolerance_m:
            final_diameter = new_diameter
            converged = True
            break
        previous = current
        current = new_diameter
    if not converged:
        final_diameter = current

    return ClosedLoopDepositionSolverResult(
        iterations_taken=iterations_taken,
        converged=converged,
        final_effective_diameter=from_si(final_diameter, "m"),
        diameter_history=tuple(from_si(value, "m") for value in diameter_history),
        max_thickness_history=tuple(from_si(value, "m") for value in thickness_history),
        velocity_history=tuple(from_si(value, "m/s") for value in velocity_history),
        warnings=tuple(warnings),
    )
