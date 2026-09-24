"""``process.absorption_column`` - the tray absorber, as a registered id.

Spec: ``specs/models/process/absorption_column.toml``.

The Python twin of ``crates/azoth-process/src/models/absorption_column.rs``, whose arithmetic is
the column's own: `AbsorptionColumn extends DistillationColumn` and overrides **no `run`**, so
what this module adds is a shape - no condenser, no reboiler, the gas at stage 0 and the solvent
at the top stage - and the refusals.

**The class's own tests solve an isothermal column, and that is a one-sweep one.** They pin every
stage with `SimpleTray.setOutletTemperature`, which makes the base's own gate - the mean
tray-temperature change - exactly zero, so the solve stops after its first sweep. The capture
holds those rows as evidence and the unpinned columns as the oracles.
"""

from __future__ import annotations

from azoth.core.errors import InvalidInputError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import AbsorptionColumnResult
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning
from azoth.process.reference.distillation_column import (
    UNPORTED_SOLVERS,
    _feed,
    _si,
    _states,
    unported_solver_class,
)

#: The class's own `DEFAULT_MAX_ALLOWABLE_GAS_LOAD_FACTOR`.
DEFAULT_MAX_ALLOWABLE_GAS_LOAD_FACTOR = 0.15


def absorption_column(
    gas_components: list[str],
    gas_n: Q,
    gas_z: list[float],
    gas_p: Q,
    gas_t: Q,
    solvent_components: list[str],
    solvent_n: Q,
    solvent_z: list[float],
    solvent_p: Q,
    solvent_t: Q,
    number_of_stages: int,
    top_pressure: Q,
    bottom_pressure: Q,
    tray_temperatures: list[float] | None = None,
    temperature_tolerance: float = 1.0e-2,
    max_iterations: int = 80,
    murphree_efficiency: float | None = None,
    component_murphree_efficiency: list[float] | None = None,
    max_allowable_gas_load_factor: float | None = None,
    solver_type: str | None = None,
) -> AbsorptionColumnResult:
    """Solve a tray absorber.

    Args:
        gas_components: the gas's substances, by name.
        gas_n: the gas feed's molar flow.
        gas_z: the gas feed's composition.
        gas_p: the gas feed's pressure.
        gas_t: the gas feed's temperature.
        solvent_components: the solvent's substances, in the gas's own order.
        solvent_n: the solvent feed's molar flow.
        solvent_z: the solvent feed's composition.
        solvent_p: the solvent feed's pressure.
        solvent_t: the solvent feed's temperature.
        number_of_stages: the trays, numbered 0 at the bottom where the gas enters.
        top_pressure: the pressure at the top stage.
        bottom_pressure: the pressure at the bottom stage.
        tray_temperatures: one outlet-temperature pin per tray, `NaN` where a tray has none.
            **A pinned column stops after its first sweep**, because the gate is the
            tray-temperature change it made zero.
        temperature_tolerance: the gate on the mean tray-temperature change.
        max_iterations: the iteration cap.
        murphree_efficiency: **not ported**; `SimpleTray.setMurphreeEfficiency` would close it.
        component_murphree_efficiency: **not ported**; the per-component override
            `applyMurphreeCorrection` applies.
        max_allowable_gas_load_factor: a design limit the solve does not read.
        solver_type: the base's strategy, by `process.distillation_column`'s own names.

    Returns:
        The tray profile, the treated gas and the loaded solvent.

    Raises:
        InvalidInputError: for a declared parameter whose arithmetic is not ported, or an inlet
            pair that does not carry the same substances in the same order.
        SolverNotConvergedError: when the solve misses its gate.

    See :func:`azoth.process.reference.distillation_column._states`.
    """
    _refuse_unported(
        murphree_efficiency, component_murphree_efficiency, max_allowable_gas_load_factor
    )
    if solver_type is not None and solver_type not in ("direct_substitution", "naphtali_sandholm"):
        raise InvalidInputError(
            "solver_type",
            f"solver_type = {solver_type} is not ported: `ColumnSolverFactory."
            f"{unported_solver_class(solver_type)}` is the class that would close it. The "
            f"strategies are path variants of one another rather than different physics - "
            f"`process.distillation_column`'s `solver_type` carries the measurement.",
        )
    if list(gas_components) != list(solvent_components):
        raise InvalidInputError(
            "solvent_components",
            "the two inlets of one column carry the same substances in the same order: the "
            "class's own feed handling indexes both by the same component list, and a solvent "
            "that names different substances in a different order is a different fluid, not an "
            "inlet",
        )

    spec = _spec()
    checks = checks_for(spec)
    warnings: list[Warning] = []

    gas_pressure = input_to_si(spec, "gas_p", gas_p)
    gas_temperature = input_to_si(spec, "gas_t", gas_t)
    solvent_pressure = input_to_si(spec, "solvent_p", solvent_p)
    solvent_temperature = input_to_si(spec, "solvent_t", solvent_t)
    top = input_to_si(spec, "top_pressure", top_pressure)
    bottom = input_to_si(spec, "bottom_pressure", bottom_pressure)
    stages = int(_si(spec, "number_of_stages", number_of_stages))
    tolerance = _si(spec, "temperature_tolerance", temperature_tolerance)
    iterations_cap = int(_si(spec, "max_iterations", max_iterations))
    gas_flow = input_to_si(spec, "gas_n", gas_n)
    solvent_flow = input_to_si(spec, "solvent_n", solvent_n)

    apply_checks(
        checks.on_input,
        {
            "number_of_stages": float(stages),
            "top_pressure": top,
            "bottom_pressure": bottom,
            "temperature_tolerance": tolerance,
            "gas_t": gas_temperature,
            "solvent_t": solvent_temperature,
        }.get,
        warnings,
    )

    # **The gas is the column's own `feed` and the solvent its `top_feed`**, which is the
    # position `addGasInStream` and `addSolventInStream` give them: stage 0 and the top stage.
    states = _states(
        list(gas_components),
        gas_flow,
        list(gas_z),
        gas_temperature,
        gas_pressure,
        stages,
        0,
        False,
        False,
        top,
        bottom,
        None,
        None,
        tolerance,
        iterations_cap,
        None,
        None,
        solver_type,
        top_feed=_feed(
            list(solvent_components),
            solvent_flow,
            list(solvent_z),
            solvent_temperature,
            solvent_pressure,
        ),
        tray_temperatures=None if tray_temperatures is None else tuple(tray_temperatures),
    )

    return AbsorptionColumnResult(
        tray_temperature=[from_si(value, "K") for value in states.tray_temperature],
        tray_pressure=[from_si(value, "Pa") for value in states.tray_pressure],
        tray_gas_n=[from_si(value, "mol/s") for value in states.tray_gas_n],
        tray_liquid_n=[from_si(value, "mol/s") for value in states.tray_liquid_n],
        gas_out_n=from_si(states.distillate_n, "mol/s"),
        gas_out_z=states.distillate_z,
        gas_out_p=from_si(states.distillate_p, "Pa"),
        gas_out_t=from_si(states.distillate_t, "K"),
        gas_out_h=from_si(states.distillate_h, "J/mol"),
        liquid_out_n=from_si(states.bottoms_n, "mol/s"),
        liquid_out_z=states.bottoms_z,
        liquid_out_p=from_si(states.bottoms_p, "Pa"),
        liquid_out_t=from_si(states.bottoms_t, "K"),
        liquid_out_h=from_si(states.bottoms_h, "J/mol"),
        iterations=states.iterations,
        temperature_residual=states.temperature_residual,
        mass_residual=states.mass_residual,
        energy_residual=states.energy_residual,
        warnings=tuple(warnings),
    )


def _spec() -> dict[str, object]:
    """This model's generated spec."""
    import azoth._models_gen as models

    return models.model("process.absorption_column")


def _refuse_unported(
    murphree_efficiency: float | None,
    component_murphree_efficiency: list[float] | None,
    max_allowable_gas_load_factor: float | None,
) -> None:
    """Refuse every parameter the palette declares and this stage does not implement."""
    if murphree_efficiency is not None:
        raise InvalidInputError(
            "murphree_efficiency",
            f"a Murphree efficiency of {murphree_efficiency} is not ported: "
            f"`SimpleTray.setMurphreeEfficiency` and the per-tray correction "
            f"`applyMurphreeCorrection` applies are the classes that would close it. Omitted "
            f"means the ideal stage, which is the class's own default of one",
        )
    if component_murphree_efficiency is not None:
        raise InvalidInputError(
            "component_murphree_efficiency",
            f"{len(component_murphree_efficiency)} component Murphree efficiencies are not "
            f"ported: `AbsorptionColumn.setComponentMurphreeEfficiency(int, String, double)` "
            f"and the `applyMurphreeCorrection` override it feeds are the classes that would "
            f"close them",
        )
    # The design limit: read by `isGasLoadFactorWithinDesignLimit`,
    # `getGasLoadFactorUtilization` and `getMinimumDiameterForGasLoadLimit`, and by nothing on
    # the run path - so it is accepted and the separation is indifferent to it.
    _ = max_allowable_gas_load_factor


__all__ = ["DEFAULT_MAX_ALLOWABLE_GAS_LOAD_FACTOR", "UNPORTED_SOLVERS", "absorption_column"]
