"""``process.gas_scrubber`` - the separator's kernel under the other entry, and the one
number the class adds to it.

Spec: ``specs/models/process/gas_scrubber.toml``

The Python twin of ``crates/azoth-process/src/models/gas_scrubber.rs`` over
``crates/azoth-process/src/kernels/gas_scrubber.rs``, written to mirror them: the steady state
is the separator's flash, and ``getCapacityUtilization`` is the whole of the difference.

# The vessel holds the temperature, and that is not a throttling

``Separator.run`` sets the pressure to ``P - pressureDrop`` and runs a ``TPflash``, so the
outlets come out at the **feed's** temperature - not at a temperature a caller states, and
not at the one an isenthalpic drop would give. Measured on the captured fluid, a 2 bar drop
takes -2705.35 J/mol in to a weighted -2364.95 out. That is the class's own model of the
vessel, and reproducing it means the balance is open whenever a pressure drop is asked for.
A ``heat_input`` is the one thing that moves the flash, and it moves it by ``Q / n`` J/mol -
``Heater.run``'s one-second convention, where the duty is W and the flow is mol/s.

# Entrainment moves material and not energy

``gas_in_liquid`` carries that fraction of the vapour's moles into the liquid outlet,
component by component, which is NeqSim's ``addPhaseFractionToPhase`` on its
``feed``/``mole`` basis. Both outlets are then at the same temperature, so the transfer
does not close the balance either - and it acts only where both phases exist, because a
from-phase or a to-phase that is not there has nothing to move.

# The capacity metric reads the *entrained* system

``getCapacityUtilization`` is ``Q / (K sqrt((rho_l - rho_g)/rho_g) A)`` with ``A = pi d**2/4``.
``Separator.run`` applies the entrainment to the system and only then does the getter read it, so
``Q`` is the vapour that is *left* and ``rho_l`` is the oil that has absorbed the carried moles -
measured on the capture's pair, ``0.0033609430039118087`` against the dry ``0.003341057093792472``.
Two rules of the class's own are worth naming because the columns' capacity families do neither:
``rho_l`` is the ``oil`` phase's density, ``aqueous``'s where there is no oil and ``1000.0`` where
there is neither (no ``10.0`` floor), and a stream with no gas phase answers ``0.0``. NeqSim
answers ``NaN`` when the liquid is no denser than the gas; this refuses that state rather than
carrying it.
"""

from __future__ import annotations

import math

from azoth.core.errors import InvalidInputError, OutOfRangeError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import GasScrubberResult, Phase
from azoth.core.units import Q, from_si, input_to_si, quantity
from azoth.core.warnings import Warning
from azoth.eos import components as _components
from azoth.eos.reference.hydrate_inhibitor_wt import GAS
from azoth.eos.reference.ph_flash import enthalpy_at
from azoth.eos.reference.ph_flash import ph_flash as ph_flash_solve
from azoth.eos.reference.pt_flash import pt_flash
from azoth.process.reference._phase_outlet import phase_enthalpy
from azoth.process.reference.capacity import _labelled_phases, _mass_density, _molar_mass
from azoth.process.reference.separator import MechanicalDesign

#: ``Separator.getMaxAllowableGasVelocity``'s own substitution where the system has no liquid
#: phase at all - a bare ``1000.0`` beside the ``50.0`` it puts on an absent gas phase, which a
#: stream with no gas never reaches because ``getCapacityUtilization`` answers zero first.
DEFAULT_LIQUID_DENSITY = 1000.0


def gas_scrubber(
    components: list[str],
    feed_n: Q,
    feed_z: list[float],
    feed_p: Q,
    feed_t: Q,
    pressure_drop: Q,
    gas_in_liquid: float,
    internal_diameter: Q | None = None,
    design_gas_load_factor: Q | None = None,
    heat_input: Q | None = None,
    mechanical_design: bool | None = None,
    orientation: str | None = None,
    volume_safety_factor: float | None = None,
    design_liquid_level_fraction: float | None = None,
    max_operation_pressure: Q | None = None,
    allowable_stress: Q | None = None,
    joint_efficiency: float | None = None,
    corrosion_allowance: Q | None = None,
    safety_factor: float | None = None,
    mesh_pad_area: float | None = None,
    mesh_pad_thickness: Q | None = None,
    mesh_pad_elevation: Q | None = None,
) -> GasScrubberResult:
    """Flash a stream into vapour and liquid outlets.

    Args:
        components: the fluid's substances, by name.
        feed_n: the inlet's molar flow.
        feed_z: the inlet composition.
        feed_p: the inlet pressure.
        feed_t: the inlet temperature, which the vessel holds.
        pressure_drop: the pressure drop across the vessel.
        gas_in_liquid: the fraction of the vapour's moles carried into the liquid outlet.
        internal_diameter: the vessel's internal diameter, or ``None`` for no capacity metric.
        design_gas_load_factor: the Souders-Brown ``K`` the vessel is designed to, or ``None``.
        heat_input: heat added before the flash, or ``None`` for a flash at the feed's
            temperature.

    Returns:
        Both outlets' records.

    Raises:
        InvalidInputError: where the pressure drop takes the outlet below zero, the
            entrainment fraction is outside ``[0, 1]``, only one of the two mechanical
            inputs is stated, or the state makes the capacity metric undefined.
        OutOfRangeError: for a non-positive feed flow, or a pressure drop or fraction the
            spec bounds.

    Example:
        >>> import azoth
        >>> q = azoth.ureg.Quantity
        >>> r = separator(
        ...     ["methane", "n-butane"],
        ...     q(1.0, "mol/s"),
        ...     [0.7, 0.3],
        ...     q(20.0, "bar"),
        ...     q(300.0, "K"),
        ...     q(0.0, "Pa"),
        ... )
        >>> round(r.vapour_n.to("mol/s").magnitude, 6)
        0.818221
    """
    spec = _spec()
    checks = checks_for(spec)
    warnings: list[Warning] = []

    n = input_to_si(spec, "feed_n", feed_n)
    t = input_to_si(spec, "feed_t", feed_t)
    p = input_to_si(spec, "feed_p", feed_p)
    drop = input_to_si(spec, "pressure_drop", pressure_drop)

    diameter = (
        None
        if internal_diameter is None
        else input_to_si(spec, "internal_diameter", internal_diameter)
    )
    load_factor = (
        None
        if design_gas_load_factor is None
        else input_to_si(spec, "design_gas_load_factor", design_gas_load_factor)
    )
    design = bool(mechanical_design) if mechanical_design is not None else False
    volume_factor = None if volume_safety_factor is None else float(volume_safety_factor)
    level_fraction = (
        None if design_liquid_level_fraction is None else float(design_liquid_level_fraction)
    )
    max_pressure = (
        None
        if max_operation_pressure is None
        else input_to_si(spec, "max_operation_pressure", max_operation_pressure)
    )
    stress = (
        None
        if allowable_stress is None
        else input_to_si(spec, "allowable_stress", allowable_stress)
    )
    joint = None if joint_efficiency is None else float(joint_efficiency)
    corrosion = (
        None
        if corrosion_allowance is None
        else input_to_si(spec, "corrosion_allowance", corrosion_allowance) * 1000.0
    )
    apply_checks(
        checks.on_input,
        {
            "pressure_drop": drop,
            "gas_in_liquid": gas_in_liquid,
            "internal_diameter": diameter,
            "design_gas_load_factor": load_factor,
            "volume_safety_factor": volume_factor,
            "design_liquid_level_fraction": level_fraction,
            "max_operation_pressure": max_pressure,
            "allowable_stress": stress,
            "joint_efficiency": joint,
            "corrosion_allowance": corrosion,
        }.get,
        warnings,
    )

    if not 0.0 <= gas_in_liquid <= 1.0:
        raise InvalidInputError(
            "gas_in_liquid",
            f"an entrainment fraction is a fraction, and {gas_in_liquid} is not in [0, 1]",
        )
    # **A stated diameter needs its design basis; a design basis does not need a diameter.**
    if diameter is not None and load_factor is None:
        raise InvalidInputError(
            "design_gas_load_factor",
            "`internal_diameter` was stated and `design_gas_load_factor` was not, and the "
            "Souders-Brown capacity metric needs both",
        )
    if not design and load_factor is not None and diameter is None:
        raise InvalidInputError(
            "internal_diameter",
            "`design_gas_load_factor` was stated without `internal_diameter` and without "
            "`mechanical_design`, and the design basis has no vessel to be used on",
        )
    if design and load_factor is None:
        raise InvalidInputError(
            "design_gas_load_factor",
            "`mechanical_design` was stated without `design_gas_load_factor`, and the design's "
            "Souders-Brown basis is that number",
        )
    if not design:
        for stated, value in (
            ("orientation", orientation),
            ("volume_safety_factor", volume_safety_factor),
            ("design_liquid_level_fraction", design_liquid_level_fraction),
            ("max_operation_pressure", max_operation_pressure),
            ("allowable_stress", allowable_stress),
            ("joint_efficiency", joint_efficiency),
            ("corrosion_allowance", corrosion_allowance),
            ("mesh_pad_area", mesh_pad_area),
            ("mesh_pad_thickness", mesh_pad_thickness),
            ("mesh_pad_elevation", mesh_pad_elevation),
        ):
            if value is not None:
                raise InvalidInputError(
                    stated, f"`{stated}` was stated without `mechanical_design`"
                )
    kind = "vertical" if orientation is None else orientation
    if kind not in ("horizontal", "vertical"):
        raise InvalidInputError(
            "orientation",
            f"a vessel is `horizontal` or `vertical`, and {kind!r} is neither",
        )
    if p <= 0.0:
        raise OutOfRangeError("P", p, "an absolute pressure; zero and below are not states")
    p_out = p - drop
    if p_out <= 0.0:
        raise InvalidInputError(
            "pressure_drop",
            f"a pressure drop of {drop} Pa takes a {p} Pa inlet to {p_out} Pa, which is "
            f"not a pressure",
        )
    if not math.isfinite(n) or not n > 0.0:
        raise OutOfRangeError("feed_n", n, "a unit operation needs a positive molar flow")

    mixture, ideal_gas = _components.mixture_of(components, eos="pr")

    if heat_input is None:
        # The vessel holds the feed's temperature, so the flash is a TP one and the
        # outlets share the temperature they came in at.
        flash = pt_flash(mixture, feed_t, from_si(p_out, "Pa"), feed_z)
        temperature = feed_t
        x, y = list(flash.x), list(flash.y)
        beta = _vapour_fraction(flash.phase, flash.vapour_fraction)
    else:
        # W over (mol/s) is J/mol, for the reason `Heater.run` gives.
        h_in, _ = enthalpy_at(mixture, ideal_gas, t, p, feed_z)
        duty = input_to_si(spec, "heat_input", heat_input)
        moved = ph_flash_solve(
            mixture,
            ideal_gas,
            from_si(p_out, "Pa"),
            from_si(h_in + duty / n, "J/mol"),
            feed_z,
        )
        temperature = moved.T
        x, y = list(moved.x), list(moved.y)
        beta = _vapour_fraction(moved.phase, moved.vapour_fraction)

    n_vapour = n * beta
    n_liquid = n * (1.0 - beta)
    # Only where both phases exist: a from-phase or a to-phase that is not there has
    # nothing to move, which is NeqSim's own guard.
    carried = n_vapour * gas_in_liquid if n_vapour > 0.0 and n_liquid > 0.0 else 0.0

    liquid_n = n_liquid + carried
    liquid_z = (
        [(n_liquid * x[i] + carried * y[i]) / liquid_n for i in range(len(x))]
        if carried > 0.0
        else x
    )

    # **The outlet is a phase, except where the class re-runs it as a stream.** `GasScrubber`
    # inherits `Separator.run`, so the same branch decides the liquid's enthalpy: a phase root
    # while nothing was carried into it, a re-flash of the carried composition once something
    # was. See `_phase_outlet`, which is where the difference is measured.
    t_out = temperature.to("K").magnitude
    vapour_h = phase_enthalpy(mixture, ideal_gas, t_out, p_out, y, liquid=False)
    if gas_in_liquid != 0.0:
        liquid_h, _ = enthalpy_at(mixture, ideal_gas, t_out, p_out, list(liquid_z))
    else:
        liquid_h = phase_enthalpy(mixture, ideal_gas, t_out, p_out, list(liquid_z), liquid=True)

    # The class reads the system `run` left, so the metric is at the *outlet* state - the feed's
    # temperature unless a `heat_input` moved the flash.
    capacity = (
        None
        if diameter is None or load_factor is None
        else _capacity_utilization(
            components,
            t_out,
            p_out,
            list(feed_z),
            n_vapour - carried,
            list(liquid_z),
            diameter,
            load_factor,
        )
    )
    required_diameter = None
    if diameter is not None and load_factor is not None and safety_factor is not None:
        from azoth.process.reference.separator import _sizing_getters

        sizing = _sizing_getters(
            components,
            t_out,
            p_out,
            list(feed_z),
            n_vapour - carried,
            list(liquid_z),
            diameter,
            load_factor,
            kind,
            0.8 if level_fraction is None else level_fraction,
            diameter * 5.0,
            0.5,
            float(safety_factor),
        )
        if sizing is not None and sizing.required_diameter is not None:
            required_diameter = from_si(sizing.required_diameter, "m")
    gas_fraction = (
        1.0 if kind != "horizontal" else 1.0 - (0.8 if level_fraction is None else level_fraction)
    )
    vessel = (
        None
        if not design or load_factor is None
        else _scrubber_design(
            components,
            t_out,
            p_out,
            list(y),
            n_vapour - carried,
            list(liquid_z),
            liquid_n,
            kind == "horizontal",
            load_factor,
            gas_fraction,
            1.1 if volume_factor is None else volume_factor,
            100.0e5 if max_pressure is None else max_pressure,
            (425.0e6 if stress is None else stress) / 1.0e6,
            1.0 if joint is None else joint,
            0.0 if corrosion is None else corrosion,
        )
    )

    return GasScrubberResult(
        vapour_n=from_si(n_vapour - carried, "mol/s"),
        vapour_z=tuple(y),
        vapour_p=from_si(p_out, "Pa"),
        vapour_t=temperature,
        vapour_h=from_si(vapour_h, "J/mol"),
        liquid_n=from_si(liquid_n, "mol/s"),
        liquid_z=tuple(liquid_z),
        liquid_p=from_si(p_out, "Pa"),
        liquid_t=temperature,
        liquid_h=from_si(liquid_h, "J/mol"),
        capacity_utilization=capacity,
        required_diameter=required_diameter,
        has_mesh_pad=(
            None
            if not design
            else (
                mesh_pad_area is not None
                and input_to_si(spec, "mesh_pad_area", mesh_pad_area) > 0.0
            )
        ),
        mesh_pad_area=(
            None
            if not design
            else from_si(
                0.0 if mesh_pad_area is None else input_to_si(spec, "mesh_pad_area", mesh_pad_area),
                "m**2",
            )
        ),
        mesh_pad_thickness=(
            None
            if not design
            else from_si(
                0.1
                if mesh_pad_thickness is None
                else input_to_si(spec, "mesh_pad_thickness", mesh_pad_thickness),
                "mm",
            )
        ),
        mesh_pad_elevation=(
            None
            if not design
            else from_si(
                0.0
                if mesh_pad_elevation is None
                else input_to_si(spec, "mesh_pad_elevation", mesh_pad_elevation),
                "m",
            )
        ),
        has_inlet_cyclones=None if not design else False,
        has_demisting_cyclones=None if not design else False,
        design_inner_diameter=(None if vessel is None else from_si(vessel.inner_diameter, "m")),
        design_tantan_length=(None if vessel is None else from_si(vessel.tantan_length, "m")),
        design_outer_diameter=(None if vessel is None else from_si(vessel.outer_diameter, "m")),
        design_wall_thickness=(None if vessel is None else from_si(vessel.wall_thickness, "m")),
        max_design_gas_volume_flow=(
            None
            if vessel is None
            else from_si(vessel.max_design_gas_volume_flow_m3_hr / 3600.0, "m**3/h")
        ),
        vessel_shell_weight=(
            None if vessel is None else from_si(vessel.vessel_shell_weight_kg, "kg")
        ),
        piping_weight=None if vessel is None else from_si(vessel.piping_weight_kg, "kg"),
        structural_weight=(None if vessel is None else from_si(vessel.structural_weight_kg, "kg")),
        electrical_weight=(None if vessel is None else from_si(vessel.electrical_weight_kg, "kg")),
        internals_weight=(None if vessel is None else from_si(vessel.internals_weight_kg, "kg")),
        total_weight=None if vessel is None else from_si(vessel.total_weight_kg, "kg"),
        module_width=None if vessel is None else from_si(vessel.module_width, "m"),
        module_length=None if vessel is None else from_si(vessel.module_length, "m"),
        module_height=None if vessel is None else from_si(vessel.module_height, "m"),
        vessel_volume=(None if vessel is None else from_si(vessel.vessel_volume_m3, "m**3")),
        hhll=None if vessel is None else from_si(vessel.hhll, "m"),
        hll=None if vessel is None else from_si(vessel.hll, "m"),
        nll=None if vessel is None else from_si(vessel.nll, "m"),
        lll=None if vessel is None else from_si(vessel.lll, "m"),
        llll=None if vessel is None else from_si(vessel.llll, "m"),
        weir_height=None if vessel is None else from_si(vessel.weir_height, "m"),
        nil=None if vessel is None else from_si(vessel.nil, "m"),
        effective_length_liquid=(
            None if vessel is None else from_si(vessel.effective_length_liquid, "m")
        ),
        effective_length_gas=(
            None if vessel is None else from_si(vessel.effective_length_gas, "m")
        ),
        inlet_nozzle_id=(None if vessel is None else from_si(vessel.inlet_nozzle_id, "m")),
        gas_outlet_nozzle_id=(
            None if vessel is None else from_si(vessel.gas_outlet_nozzle_id, "m")
        ),
        oil_outlet_nozzle_id=(
            None if vessel is None else from_si(vessel.oil_outlet_nozzle_id, "m")
        ),
        surge_volume=(None if vessel is None else from_si(vessel.surge_volume_m3, "m**3")),
        surge_time=None if vessel is None else from_si(vessel.surge_time_s, "s"),
        liquid_retention_time=(
            None if vessel is None else from_si(vessel.liquid_retention_time_s, "s")
        ),
        warnings=tuple(warnings),
    )


def _scrubber_design(
    components: list[str],
    t: float,
    p: float,
    gas_z: list[float],
    gas_n: float,
    liquid_z: list[float],
    liquid_n: float,
    horizontal: bool,
    gas_load_factor: float,
    fg: float,
    volume_safety_factor: float,
    max_operation_pressure_pa: float,
    allowable_stress_mpa: float,
    joint_efficiency: float,
    corrosion_allowance_mm: float,
) -> MechanicalDesign:
    """The scrubber's mechanical design: the separator's arithmetic as a `Scrubber` vessel."""
    from azoth.process.reference.separator import _mechanical_design

    return _mechanical_design(
        components,
        t,
        p,
        gas_z,
        gas_n,
        liquid_z,
        liquid_n,
        horizontal,
        gas_load_factor,
        fg,
        volume_safety_factor,
        60.0,
        max_operation_pressure_pa,
        allowable_stress_mpa,
        joint_efficiency,
        corrosion_allowance_mm,
        vessel="scrubber",
    )


def _capacity_utilization(
    components: list[str],
    t: float,
    p: float,
    z: list[float],
    gas_n: float,
    liquid_z: list[float],
    diameter_m: float,
    design_gas_load_factor: float,
) -> float:
    """``Separator.getCapacityUtilization``, on the system the run left behind.

    ``z`` is the feed the flash started from, because the phase *types* the class asks for are
    the ones that flash assigned - ``hasPhaseType`` reads a label that ``addPhaseFractionToPhase``
    does not revisit - while ``liquid_z`` is what the oil phase *holds* now, because
    ``initPhysicalProperties`` does recompute its density. The two are the same on a run with no
    entrainment, which is what the capture's dry row fixes.
    """
    mixture, labelled, _shares = _labelled_phases(components, t, p, z, from_feed=True)
    gas = next((entry for entry in labelled if entry[0] == GAS), None)
    if gas is None:
        # The class's own first answer: a stream with no gas phase uses none of its capacity.
        return 0.0
    _gas_label, gas_z, gas_root = gas
    gas_density = _mass_density(mixture, t, p, gas_z, gas_root)
    gas_molar_mass = _molar_mass(mixture, gas_z)

    if any(entry[0] != GAS for entry in labelled):
        # The oil's density is the *entrained* composition's, at the root that composition's own
        # flash puts it on - `Separator.run` adds the carried moles before the getter reads.
        entrained = pt_flash(mixture, quantity(t, "K"), quantity(p, "Pa"), list(liquid_z))
        liquid_density = _mass_density(mixture, t, p, list(liquid_z), entrained.z_liquid)
    else:
        liquid_density = DEFAULT_LIQUID_DENSITY

    if math.isnan(liquid_density) or math.isnan(gas_density) or liquid_density <= gas_density:
        raise InvalidInputError(
            "capacity_utilization",
            f"the liquid phase is {liquid_density} kg/m3 against the gas's {gas_density}, so "
            f"`rho_l - rho_g` is not positive and the Souders-Brown velocity is the square root "
            f"of a negative number. NeqSim's `getCapacityUtilization` answers NaN there; this "
            f"refuses the state rather than answering it",
        )

    area = math.pi * diameter_m * diameter_m / 4.0
    velocity = design_gas_load_factor * math.sqrt((liquid_density - gas_density) / gas_density)
    volumetric_flow = gas_n * gas_molar_mass / gas_density
    return volumetric_flow / (velocity * area)


def _vapour_fraction(phase: Phase, beta: float | None) -> float:
    """The vapour fraction a flash's phase and ``beta`` imply.

    ``beta`` is absent for a single-phase answer - the flash reports ``None`` rather than a
    number outside ``[0, 1]``, because there is no split to report - so the phase is what
    says which side of the interval the outlet is on.
    """
    if phase in (Phase.ALL_VAPOUR, Phase.TRIVIAL):
        return 1.0
    return 0.0 if beta is None else float(beta)


def _spec() -> dict[str, object]:
    """The generated spec, imported here so the module has no import-time cycle."""
    from azoth._models_gen import model

    return model("process.gas_scrubber")
