"""``process.filter`` - a pressure drop, flashed at the feed's temperature.

Spec: ``specs/models/process/filter.toml``

The Python twin of ``crates/azoth-process/src/models/filter.rs`` over
``crates/azoth-process/src/kernels/filter.rs``, written to mirror them.

# It holds the temperature, which is what separates it from a throttling valve

``Filter.run`` reduces the pressure and runs a ``TPflash``, so the outlet is at the feed's
temperature and its enthalpy moves with the pressure. A valve is a ``PHflash`` and holds the
enthalpy instead. Measured on the fluid these cases share, one bar of drop moves the enthalpy
``24.8`` J/mol - the same pressure dependence the no-specification branch of
``process.heater`` prices at ``49.6`` over two.

# A drop past the inlet pressure is clamped rather than refused

``run`` applies ``min(max(0, dP), max(0, P_in - 1e-6 bar))`` and logs a warning, so the
outlet lands a millionth of a bar above vacuum. ``applied_drop`` carries what was applied and
the model warns when the two differ, because a result that carried the clamped value silently
would be indistinguishable from one that was asked for it.

# The clean drop and the capture

``FilterPressureDropModel``'s ``FIXED`` returns the configured drop, ``FLOW_SCALED`` scales it
by ``(Q/Q_ref)^n``, ``TABULATED`` interpolates a measured curve, and ``ERGUN`` evaluates the
packed-bed correlation. The capture is ``updateParticleCapturePerformance``: a beta-ratio curve
or a scalar efficiency, the inlet concentration the caller states, and two outputs that are
zero when it does not.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, replace
from itertools import pairwise

from azoth.core.errors import InvalidInputError, OutOfRangeError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import FilterResult
from azoth.core.units import Q, from_si, input_records_to_si, input_to_si
from azoth.core.warnings import Warning, WarningCode
from azoth.eos import components as _components
from azoth.eos.reference.ph_flash import enthalpy_at
from azoth.process.reference import capacity

#: ``Filter``'s media defaults.
DEFAULT_MEDIA_AREA_M2 = 1.0
DEFAULT_MEDIA_BED_DEPTH_M = 1.0
DEFAULT_MEDIA_PARTICLE_DIAMETER_M = 0.001
DEFAULT_MEDIA_VOID_FRACTION = 0.40
#: ``Filter.particleSizeMicrometre``'s initialiser, in metres.
DEFAULT_PARTICLE_SIZE_M = 10.0e-6
#: ``Filter.deltaP``'s clamp ceiling: a millionth of a bar, in Pa.
CLAMP_MARGIN_PA = 0.1
#: ``Filter``'s default terminal differential pressure, 5 bar, in Pa.
DEFAULT_TERMINAL_DROP_PA = 500_000.0
#: ``setBreakthroughStartFraction``'s ceiling; the divisor ``1 - start`` stays positive.
MAX_BREAKTHROUGH_START = 0.999999


def filter(
    components: list[str],
    inlet_n: Q,
    inlet_z: list[float],
    inlet_p: Q,
    inlet_t: Q,
    pressure_drop: Q,
    pressure_drop_model: str | None = None,
    reference_flow: Q | None = None,
    flow_exponent: float | None = None,
    pressure_drop_curve: list[dict[str, Q]] | None = None,
    media_area: Q | None = None,
    media_bed_depth: Q | None = None,
    media_particle_diameter: Q | None = None,
    media_void_fraction: float | None = None,
    viscosity: Q | None = None,
    performance_curve: list[dict[str, Q]] | None = None,
    particle_size: Q | None = None,
    nominal_removal_efficiency: float | None = None,
    inlet_particle_concentration: float | None = None,
    solids_loading: Q | None = None,
    loading_capacity: Q | None = None,
    pressure_drop_increase_at_capacity: Q | None = None,
    breakthrough_start_fraction: float | None = None,
    bypass_cracking_pressure: Q | None = None,
    terminal_pressure_drop: Q | None = None,
    time_step: Q | None = None,
    holdup_volume: Q | None = None,
    solids_loading_rate: Q | None = None,
    backwash_removal_rate: Q | None = None,
    regeneration_removal_rate: Q | None = None,
    backwash_active: bool | None = None,
    regeneration_active: bool | None = None,
) -> FilterResult:
    """Drop a stream's pressure by the selected clean-drop model at a constant temperature.

    Args:
        components: the fluid's substances, by name.
        inlet_n: the inlet's molar flow.
        inlet_z: the inlet composition.
        inlet_p: the inlet pressure.
        inlet_t: the inlet temperature, which is also the outlet's.
        pressure_drop: the configured clean pressure loss.
        pressure_drop_model: ``fixed``, ``flow_scaled``, ``tabulated`` or ``ergun``.
        reference_flow: the ``flow_scaled`` reference flow.
        flow_exponent: the ``flow_scaled`` exponent.
        pressure_drop_curve: the ``tabulated`` flow/drop rows.
        media_area: the Ergun bed's area.
        media_bed_depth: the Ergun bed's depth.
        media_particle_diameter: the Ergun bed's particle diameter.
        media_void_fraction: the Ergun bed's void fraction.
        viscosity: the fluid's dynamic viscosity, required by ``ergun``.
        performance_curve: the beta-ratio rows.
        particle_size: the size the beta curve is evaluated at.
        nominal_removal_efficiency: the scalar efficiency without a beta curve.
        inlet_particle_concentration: the inlet concentration as a mass fraction; stating it
            enables the capture.

    Returns:
        The outlet's record, the clean and applied drops, and the capture's answers.

    Raises:
        InvalidInputError: where the shapes disagree or a model's own settings are not one.
        OutOfRangeError: for a temperature the spec bounds.

    Example:
        >>> import azoth
        >>> q = azoth.ureg.Quantity
        >>> r = filter(
        ...     ["methane", "n-butane"],
        ...     1.0,
        ...     q([0.9, 0.1]),
        ...     q(30.0, "bar"),
        ...     q(320.0, "K"),
        ...     q(1.0, "bar"),
        ... )
        >>> round(r.outlet_t.to("K").magnitude, 1)
        320.0
    """
    spec = _spec()
    checks = checks_for(spec)
    warnings: list[Warning] = []

    n = input_to_si(spec, "inlet_n", inlet_n)
    t = input_to_si(spec, "inlet_t", inlet_t)
    p = input_to_si(spec, "inlet_p", inlet_p)
    drop = input_to_si(spec, "pressure_drop", pressure_drop)
    reference_si = (
        None if reference_flow is None else input_to_si(spec, "reference_flow", reference_flow)
    )
    viscosity_si = None if viscosity is None else input_to_si(spec, "viscosity", viscosity)
    solids_si = (
        None if solids_loading is None else input_to_si(spec, "solids_loading", solids_loading)
    )
    capacity_si = (
        None
        if loading_capacity is None
        else input_to_si(spec, "loading_capacity", loading_capacity)
    )
    increase_si = (
        None
        if pressure_drop_increase_at_capacity is None
        else input_to_si(
            spec, "pressure_drop_increase_at_capacity", pressure_drop_increase_at_capacity
        )
    )
    bypass_si = (
        None
        if bypass_cracking_pressure is None
        else input_to_si(spec, "bypass_cracking_pressure", bypass_cracking_pressure)
    )
    terminal_si = (
        None
        if terminal_pressure_drop is None
        else input_to_si(spec, "terminal_pressure_drop", terminal_pressure_drop)
    )
    step_si = None if time_step is None else input_to_si(spec, "time_step", time_step)
    holdup_si = None if holdup_volume is None else input_to_si(spec, "holdup_volume", holdup_volume)
    rate_si = (
        None
        if solids_loading_rate is None
        else input_to_si(spec, "solids_loading_rate", solids_loading_rate)
    )
    backwash_si = (
        None
        if backwash_removal_rate is None
        else input_to_si(spec, "backwash_removal_rate", backwash_removal_rate)
    )
    regen_si = (
        None
        if regeneration_removal_rate is None
        else input_to_si(spec, "regeneration_removal_rate", regeneration_removal_rate)
    )

    apply_checks(
        checks.on_input,
        {
            "inlet_t": t,
            "pressure_drop": drop,
            "flow_exponent": flow_exponent,
            "media_void_fraction": media_void_fraction,
            "nominal_removal_efficiency": nominal_removal_efficiency,
            "solids_loading": solids_si,
            "pressure_drop_increase_at_capacity": increase_si,
            "breakthrough_start_fraction": breakthrough_start_fraction,
            "bypass_cracking_pressure": bypass_si,
            "terminal_pressure_drop": terminal_si,
            "time_step": step_si,
            "holdup_volume": holdup_si,
            "solids_loading_rate": rate_si,
            "backwash_removal_rate": backwash_si,
            "regeneration_removal_rate": regen_si,
        }.get,
        warnings,
    )

    model = "fixed" if pressure_drop_model is None else pressure_drop_model
    if model not in ("fixed", "flow_scaled", "tabulated", "ergun"):
        raise InvalidInputError(
            "pressure_drop_model",
            f"{model} is not one of fixed, flow_scaled, tabulated or ergun",
        )
    exponent = 1.0 if flow_exponent is None else float(flow_exponent)
    if not exponent > 0.0:
        raise InvalidInputError(
            "flow_exponent",
            f"the scaling exponent is {exponent}, and a non-positive one is not a curve",
        )
    if reference_si is not None and not reference_si > 0.0:
        raise InvalidInputError(
            "reference_flow",
            f"the reference flow is {reference_si} m3/s, and the scaling divides by it",
        )
    drop_curve = _drop_curve(spec, pressure_drop_curve)
    beta_curve = _beta_curve(spec, performance_curve)

    media = (
        DEFAULT_MEDIA_AREA_M2
        if media_area is None
        else input_to_si(spec, "media_area", media_area),
        DEFAULT_MEDIA_BED_DEPTH_M
        if media_bed_depth is None
        else input_to_si(spec, "media_bed_depth", media_bed_depth),
        DEFAULT_MEDIA_PARTICLE_DIAMETER_M
        if media_particle_diameter is None
        else input_to_si(spec, "media_particle_diameter", media_particle_diameter),
        DEFAULT_MEDIA_VOID_FRACTION if media_void_fraction is None else float(media_void_fraction),
    )
    if model == "ergun":
        if viscosity_si is None:
            raise InvalidInputError(
                "viscosity",
                "the Ergun model reads the fluid's viscosity, and this filter was not given one",
            )
        area, depth, diameter, void = media
        if area <= 0.0 or depth <= 0.0 or diameter <= 0.0:
            raise InvalidInputError(
                "media_area",
                "an Ergun bed needs a positive area, depth and particle diameter",
            )
        if not 0.0 < void < 1.0:
            raise InvalidInputError(
                "media_void_fraction",
                f"the void fraction is {void}, and an Ergun bed's is in (0, 1)",
            )

    efficiency = 0.0 if nominal_removal_efficiency is None else float(nominal_removal_efficiency)
    if not 0.0 <= efficiency <= 1.0:
        raise InvalidInputError(
            "nominal_removal_efficiency",
            f"the removal efficiency is {efficiency}, and a fraction is in [0, 1]",
        )
    if inlet_particle_concentration is not None and (
        not float(inlet_particle_concentration) >= 0.0
    ):
        raise InvalidInputError(
            "inlet_particle_concentration",
            f"the inlet concentration is {inlet_particle_concentration}, not a non-negative "
            f"fraction",
        )
    for name, value in (
        ("solids_loading", solids_si),
        ("loading_capacity", capacity_si),
        ("pressure_drop_increase_at_capacity", increase_si),
        ("breakthrough_start_fraction", breakthrough_start_fraction),
        ("time_step", step_si),
        ("holdup_volume", holdup_si),
        ("solids_loading_rate", rate_si),
        ("backwash_removal_rate", backwash_si),
        ("regeneration_removal_rate", regen_si),
    ):
        if value is not None and not math.isfinite(float(value)):
            raise InvalidInputError(name, f"{value} is not a finite number")
    if solids_loading_rate is not None and inlet_particle_concentration is not None:
        raise InvalidInputError(
            "solids_loading_rate",
            "the generic captured rate and the concentration model each disable the other "
            "in `Filter`, so the last setter decides; the port refuses both rather than "
            "choosing",
        )
    if bypass_si is not None and not (math.isfinite(bypass_si) and bypass_si > 0.0):
        raise InvalidInputError(
            "bypass_cracking_pressure",
            f"the cracking pressure is {bypass_si} Pa, and `setBypassCrackingDeltaP` needs a "
            f"finite positive one",
        )
    if terminal_si is not None and not (math.isfinite(terminal_si) and terminal_si > 0.0):
        raise InvalidInputError(
            "terminal_pressure_drop",
            f"the terminal drop is {terminal_si} Pa, and `setTerminalDeltaP` needs a finite "
            f"positive one",
        )

    size_m = (
        DEFAULT_PARTICLE_SIZE_M
        if particle_size is None
        else input_to_si(spec, "particle_size", particle_size)
    )

    loading = LoadingRequest(
        solids_loading_kg=solids_si,
        loading_capacity_kg=capacity_si,
        pressure_drop_increase_at_capacity_pa=increase_si,
        breakthrough_start_fraction=breakthrough_start_fraction,
        bypass_cracking_pressure_pa=bypass_si,
        terminal_pressure_drop_pa=terminal_si,
        time_step_s=step_si,
        holdup_volume_m3=holdup_si,
        solids_loading_rate_kg_h=None if rate_si is None else rate_si * 3600.0,
        backwash_removal_rate_kg_h=None if backwash_si is None else backwash_si * 3600.0,
        regeneration_removal_rate_kg_h=None if regen_si is None else regen_si * 3600.0,
        backwash_active=bool(backwash_active),
        regeneration_active=bool(regeneration_active),
    )

    states = _route(
        components,
        n,
        inlet_z,
        p,
        t,
        drop,
        model,
        reference_si,
        exponent,
        drop_curve,
        media,
        None if viscosity_si is None else float(viscosity_si),
        beta_curve,
        size_m,
        efficiency,
        inlet_particle_concentration,
        loading,
    )

    if states.applied_drop != states.restricted_drop:
        warnings.append(
            Warning(
                code=WarningCode.OUT_OF_VALID_RANGE,
                field="pressure_drop",
                message=(
                    f"a restricted drop of {states.restricted_drop} Pa was computed and "
                    f"{states.applied_drop} Pa was applied: the outlet is held a millionth "
                    f"of a bar above vacuum, which is the clamp `Filter.run` applies"
                ),
            )
        )

    return FilterResult(
        outlet_n=from_si(n, "mol/s"),
        outlet_z=tuple(inlet_z),
        outlet_p=from_si(states.p_out, "Pa"),
        outlet_t=from_si(t, "K"),
        outlet_h=from_si(states.h_out, "J/mol"),
        clean_drop=from_si(states.clean_drop, "Pa"),
        applied_drop=from_si(states.applied_drop, "Pa"),
        nominal_removal_efficiency=states.efficiency,
        outlet_particle_concentration=states.outlet_concentration,
        captured_rate=from_si(states.captured_rate_kg_s, "kg/h"),
        unrestricted_drop=from_si(states.unrestricted_drop, "Pa"),
        bypass_fraction=states.bypass_fraction,
        loading_fraction=states.loading_fraction,
        breakthrough_fraction=states.breakthrough_fraction,
        current_removal_efficiency=states.current_efficiency,
        solids_loading=from_si(states.solids_loading_kg, "kg"),
        differential_pressure_utilization=states.differential_pressure_utilization,
        replacement_required=states.replacement_required,
        holdup_residence_time=(
            None if states.holdup_residence_s is None else from_si(states.holdup_residence_s, "s")
        ),
        warnings=tuple(warnings),
    )


@dataclass(frozen=True, slots=True)
class LoadingRequest:
    """``Filter``'s carried loading state, as the model id's optional inputs state it."""

    solids_loading_kg: float | None = None
    loading_capacity_kg: float | None = None
    pressure_drop_increase_at_capacity_pa: float | None = None
    breakthrough_start_fraction: float | None = None
    bypass_cracking_pressure_pa: float | None = None
    terminal_pressure_drop_pa: float | None = None
    time_step_s: float | None = None
    holdup_volume_m3: float | None = None
    solids_loading_rate_kg_h: float | None = None
    backwash_removal_rate_kg_h: float | None = None
    regeneration_removal_rate_kg_h: float | None = None
    backwash_active: bool = False
    regeneration_active: bool = False


@dataclass(frozen=True, slots=True)
class LoadingState:
    """``updateHydraulicState``'s arithmetic and the capture's current efficiency."""

    solids_loading_kg: float
    loading_fraction: float
    breakthrough_fraction: float
    efficiency: float
    unrestricted_drop: float
    restricted_drop: float
    bypass_fraction: float
    holdup_residence_s: float | None = 0.0


@dataclass(frozen=True, slots=True)
class FilterStates:
    """The filter's intermediates, in the order it computes them, in SI."""

    #: The clean drop the model computed, the pre-clamp restricted one and the applied one.
    clean_drop: float
    restricted_drop: float
    applied_drop: float
    #: The outlet pressure.
    p_out: float
    #: The outlet molar enthalpy, which is not the inlet's: the drop is isothermal.
    h_out: float
    #: The nominal removal efficiency, the current one, and the capture's two answers.
    efficiency: float
    current_efficiency: float
    outlet_concentration: float
    captured_rate_kg_s: float
    #: The carried loading state and the replacement check.
    unrestricted_drop: float
    bypass_fraction: float
    loading_fraction: float
    breakthrough_fraction: float
    solids_loading_kg: float
    differential_pressure_utilization: float
    replacement_required: bool
    holdup_residence_s: float | None


def _route(
    components: list[str],
    inlet_n: float,
    inlet_z: list[float],
    inlet_p: float,
    inlet_t: float,
    pressure_drop: float,
    model: str,
    reference_flow: float | None,
    exponent: float,
    drop_curve: list[tuple[float, float]],
    media: tuple[float, float, float, float],
    viscosity: float | None,
    beta_curve: list[tuple[float, float]],
    particle_size_m: float,
    nominal_efficiency: float,
    inlet_concentration: float | None,
    loading: LoadingRequest,
) -> FilterStates:
    """The filter's intermediates, in SI, in the order it computes them."""
    if not inlet_p > 0.0:
        # The flash the outlet passes through names the absolute pressure `P`, and the Rust
        # entry reaches that refusal rather than a model-level one; the two are one error.
        raise OutOfRangeError("P", inlet_p, "an absolute pressure; zero and below are not states")

    flow_m3h = (
        capacity.system_volumetric_flow(components, inlet_t, inlet_p, inlet_z, inlet_n) * 3600.0
    )
    clean_pa = _clean_drop(
        components,
        inlet_t,
        inlet_p,
        inlet_z,
        pressure_drop,
        model,
        reference_flow,
        exponent,
        drop_curve,
        media,
        viscosity,
        flow_m3h,
    )
    nominal = _nominal_efficiency(beta_curve, particle_size_m, nominal_efficiency)
    pre = _loading_state(clean_pa, nominal, loading)

    mixture, ideal_gas = _components.mixture_of(components, eos="pr")
    mass_flow = inlet_n * capacity._molar_mass(mixture, list(inlet_z))

    if (loading.time_step_s or 0.0) > 0.0:
        captured_kg_s = (
            mass_flow * inlet_concentration * pre.efficiency
            if inlet_concentration is not None
            else max(0.0, loading.solids_loading_rate_kg_h or 0.0) / 3600.0
        )
        state = _advance(pre, clean_pa, nominal, loading, flow_m3h, captured_kg_s)
    else:
        state = pre

    ceiling = max(0.0, inlet_p - CLAMP_MARGIN_PA)
    applied = min(max(0.0, state.restricted_drop), ceiling)
    p_out = inlet_p - applied
    h_out, _ = enthalpy_at(mixture, ideal_gas, inlet_t, p_out, inlet_z)

    if inlet_concentration is None:
        outlet_concentration = 0.0
        captured_rate = 0.0
    else:
        outlet_concentration = inlet_concentration * (1.0 - state.efficiency)
        captured_rate = mass_flow * inlet_concentration * state.efficiency

    terminal = (
        DEFAULT_TERMINAL_DROP_PA
        if loading.terminal_pressure_drop_pa is None
        else loading.terminal_pressure_drop_pa
    )
    utilization = state.unrestricted_drop / terminal
    replacement_required = utilization >= 1.0 or state.loading_fraction >= 1.0

    return FilterStates(
        clean_drop=clean_pa,
        restricted_drop=state.restricted_drop,
        applied_drop=applied,
        p_out=p_out,
        h_out=h_out,
        efficiency=nominal,
        current_efficiency=state.efficiency,
        outlet_concentration=outlet_concentration,
        captured_rate_kg_s=captured_rate,
        unrestricted_drop=state.unrestricted_drop,
        bypass_fraction=state.bypass_fraction,
        loading_fraction=state.loading_fraction,
        breakthrough_fraction=state.breakthrough_fraction,
        solids_loading_kg=state.solids_loading_kg,
        differential_pressure_utilization=utilization,
        replacement_required=replacement_required,
        holdup_residence_s=state.holdup_residence_s,
    )


def _loading_state(clean_pa: float, nominal: float, request: LoadingRequest) -> LoadingState:
    """``updateHydraulicState`` and the capture's current efficiency, before the clamp."""
    return _loading_state_at(clean_pa, nominal, request.solids_loading_kg, request)


def _loading_state_at(
    clean_pa: float, nominal: float, loading_kg: float | None, request: LoadingRequest
) -> LoadingState:
    solids_loading_kg = max(0.0, 0.0 if loading_kg is None else loading_kg)
    capacity_kg = request.loading_capacity_kg
    loading_fraction = (
        solids_loading_kg / capacity_kg
        if capacity_kg is not None and math.isfinite(capacity_kg) and capacity_kg > 0.0
        else 0.0
    )
    start = (
        0.8 if request.breakthrough_start_fraction is None else request.breakthrough_start_fraction
    )
    start = min(max(0.0, start), MAX_BREAKTHROUGH_START)
    breakthrough_fraction = (
        0.0 if loading_fraction <= start else min(1.0, (loading_fraction - start) / (1.0 - start))
    )
    increase = max(0.0, request.pressure_drop_increase_at_capacity_pa or 0.0)
    unrestricted_drop = clean_pa + increase * max(0.0, loading_fraction)
    cracking = request.bypass_cracking_pressure_pa
    if cracking is not None and unrestricted_drop > cracking:
        bypass_fraction = min(max(1.0 - math.sqrt(cracking / unrestricted_drop), 0.0), 1.0)
        restricted_drop = cracking
    else:
        bypass_fraction = 0.0
        restricted_drop = unrestricted_drop
    efficiency = min(
        max(nominal * (1.0 - breakthrough_fraction) * (1.0 - bypass_fraction), 0.0), 1.0
    )
    return LoadingState(
        solids_loading_kg=solids_loading_kg,
        loading_fraction=loading_fraction,
        breakthrough_fraction=breakthrough_fraction,
        efficiency=efficiency,
        unrestricted_drop=unrestricted_drop,
        restricted_drop=restricted_drop,
        bypass_fraction=bypass_fraction,
    )


def _advance(
    pre: LoadingState,
    clean_pa: float,
    nominal: float,
    request: LoadingRequest,
    flow_m3h: float,
    captured_kg_s: float,
) -> LoadingState:
    """``updateDynamicLoading`` and ``updateHoldupResidenceTime``, one step."""
    step = max(0.0, request.time_step_s or 0.0)
    removed_kg_s = 0.0
    if request.backwash_active:
        removed_kg_s += max(0.0, request.backwash_removal_rate_kg_h or 0.0) / 3600.0
    if request.regeneration_active:
        removed_kg_s += max(0.0, request.regeneration_removal_rate_kg_h or 0.0) / 3600.0
    loading_kg = max(0.0, pre.solids_loading_kg + (captured_kg_s - removed_kg_s) * step)
    state = _loading_state_at(clean_pa, nominal, loading_kg, request)
    holdup = max(0.0, request.holdup_volume_m3 or 0.0)
    residence: float | None
    if holdup <= 0.0:
        residence = 0.0
    else:
        flow_m3s = flow_m3h / 3600.0
        residence = None if not math.isfinite(flow_m3s) or flow_m3s <= 0.0 else holdup / flow_m3s
    return replace(state, holdup_residence_s=residence)


def _clean_drop(
    components: list[str],
    t: float,
    p: float,
    z: list[float],
    configured: float,
    model: str,
    reference_flow: float | None,
    exponent: float,
    drop_curve: list[tuple[float, float]],
    media: tuple[float, float, float, float],
    viscosity: float | None,
    flow_m3h: float,
) -> float:
    """``calculateCleanPressureDrop``, in Pa."""
    if model == "fixed":
        return configured
    if model == "flow_scaled":
        reference = flow_m3h if reference_flow is None else reference_flow * 3600.0
        if reference <= 0.0:
            return 0.0
        return float(configured * (flow_m3h / reference) ** exponent)
    if model == "tabulated":
        if not drop_curve:
            return configured
        return _curve_drop(drop_curve, flow_m3h)

    if flow_m3h <= 0.0:
        return 0.0
    density = capacity.system_mass_density(components, t, p, z)
    if viscosity is None:
        return configured
    if not density > 0.0 or not viscosity > 0.0:
        return configured
    area, depth, diameter, void = media
    u = flow_m3h / 3600.0 / area
    solid = 1.0 - void
    void_cube = void**3
    viscous = 150.0 * viscosity * solid * solid / (void_cube * diameter * diameter) * u
    inertial = 1.75 * density * solid / (void_cube * diameter) * u * u
    return max(0.0, depth * (viscous + inertial))


def _curve_drop(curve: list[tuple[float, float]], flow_m3h: float) -> float:
    """``FilterPressureDropCurve.getPressureDrop``, Pa."""
    if flow_m3h <= 0.0:
        return 0.0
    first_flow, first_drop = curve[0]
    if flow_m3h <= first_flow:
        return max(0.0, first_drop * flow_m3h / first_flow)
    for (lower_flow, lower_drop), (upper_flow, upper_drop) in pairwise(curve):
        if flow_m3h <= upper_flow:
            fraction = (flow_m3h - lower_flow) / (upper_flow - lower_flow)
            return max(0.0, lower_drop + fraction * (upper_drop - lower_drop))
    if len(curve) == 1:
        return max(0.0, first_drop * flow_m3h / first_flow)
    lower_flow, lower_drop = curve[-2]
    upper_flow, upper_drop = curve[-1]
    fraction = (flow_m3h - lower_flow) / (upper_flow - lower_flow)
    return max(0.0, lower_drop + fraction * (upper_drop - lower_drop))


def _nominal_efficiency(
    curve: list[tuple[float, float]], particle_size_m: float, scalar: float
) -> float:
    """``Filter.getNominalRemovalEfficiency``: the curve's value, else the scalar."""
    if not curve:
        return min(1.0, max(0.0, scalar))
    beta = _beta_ratio(curve, particle_size_m)
    return min(1.0, max(0.0, 1.0 - 1.0 / beta))


def _beta_ratio(curve: list[tuple[float, float]], particle_size_m: float) -> float:
    """``FilterPerformanceCurve.getBetaRatio``: linear in size, logarithmic in the ratio."""
    first_size, first_beta = curve[0]
    if particle_size_m <= first_size:
        return first_beta
    last_size, last_beta = curve[-1]
    if particle_size_m >= last_size:
        return last_beta
    for (lower_size, lower_beta), (upper_size, upper_beta) in pairwise(curve):
        if particle_size_m <= upper_size:
            fraction = (particle_size_m - lower_size) / (upper_size - lower_size)
            lower_log = math.log(lower_beta)
            upper_log = math.log(upper_beta)
            return math.exp(lower_log + fraction * (upper_log - lower_log))
    return last_beta


def _drop_curve(
    spec: dict[str, object], rows: list[dict[str, Q]] | None
) -> list[tuple[float, float]]:
    """The `pressure_drop_curve` rows as `(flow m3/h, drop Pa)`, validated."""
    if rows is None:
        return []
    converted: list[dict[str, float]] = input_records_to_si(spec, "pressure_drop_curve", rows)
    curve: list[tuple[float, float]] = []
    for index, row in enumerate(converted):
        flow_m3h = float(row["flow"]) * 3600.0
        drop_pa = float(row["drop"])
        if not flow_m3h > 0.0:
            raise InvalidInputError(
                "pressure_drop_curve",
                f"point {index}'s flow is {flow_m3h} m3/hr, and the curve is over positive flows",
            )
        if curve and flow_m3h <= curve[-1][0]:
            raise InvalidInputError(
                "pressure_drop_curve", "the curve's flows are not strictly increasing"
            )
        if drop_pa < 0.0:
            raise InvalidInputError(
                "pressure_drop_curve",
                f"point {index}'s drop is {drop_pa} Pa, and a drop is non-negative",
            )
        if curve and drop_pa < curve[-1][1]:
            raise InvalidInputError(
                "pressure_drop_curve", "the curve's drops are not non-decreasing"
            )
        curve.append((flow_m3h, drop_pa))
    return curve


def _beta_curve(
    spec: dict[str, object], rows: list[dict[str, Q]] | None
) -> list[tuple[float, float]]:
    """The `performance_curve` rows as `(particle size m, beta)`, validated."""
    if rows is None:
        return []
    converted: list[dict[str, float]] = input_records_to_si(spec, "performance_curve", rows)
    curve: list[tuple[float, float]] = []
    for index, row in enumerate(converted):
        size = float(row["particle_size"])
        beta = float(row["beta_ratio"])
        if not size > 0.0:
            raise InvalidInputError(
                "performance_curve",
                f"point {index}'s particle size is {size} m, and a size is positive",
            )
        if curve and size <= curve[-1][0]:
            raise InvalidInputError(
                "performance_curve", "the curve's particle sizes are not strictly increasing"
            )
        if beta < 1.0:
            raise InvalidInputError(
                "performance_curve",
                f"point {index}'s beta ratio is {beta}, and a ratio is at least one",
            )
        curve.append((size, beta))
    return curve


def _spec() -> dict[str, object]:
    """The generated spec, imported here so the module has no import-time cycle."""
    from azoth._models_gen import model

    return model("process.filter")
