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
from dataclasses import dataclass
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

    apply_checks(
        checks.on_input,
        {
            "inlet_t": t,
            "pressure_drop": drop,
            "flow_exponent": flow_exponent,
            "media_void_fraction": media_void_fraction,
            "nominal_removal_efficiency": nominal_removal_efficiency,
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

    size_m = (
        DEFAULT_PARTICLE_SIZE_M
        if particle_size is None
        else input_to_si(spec, "particle_size", particle_size)
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
    )

    if states.applied_drop != states.clean_drop:
        warnings.append(
            Warning(
                code=WarningCode.OUT_OF_VALID_RANGE,
                field="pressure_drop",
                message=(
                    f"a clean drop of {states.clean_drop} Pa was computed and "
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
        warnings=tuple(warnings),
    )


@dataclass(frozen=True, slots=True)
class FilterStates:
    """The filter's intermediates, in the order it computes them, in SI."""

    #: The clean drop the model computed and the one applied after the clamp.
    clean_drop: float
    applied_drop: float
    #: The outlet pressure.
    p_out: float
    #: The outlet molar enthalpy, which is not the inlet's: the drop is isothermal.
    h_out: float
    #: The nominal removal efficiency, and the capture's two answers.
    efficiency: float
    outlet_concentration: float
    captured_rate_kg_s: float


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
    ceiling = max(0.0, inlet_p - CLAMP_MARGIN_PA)
    applied = min(max(0.0, clean_pa), ceiling)
    p_out = inlet_p - applied

    mixture, ideal_gas = _components.mixture_of(components, eos="pr")
    h_out, _ = enthalpy_at(mixture, ideal_gas, inlet_t, p_out, inlet_z)

    efficiency = _nominal_efficiency(beta_curve, particle_size_m, nominal_efficiency)
    if inlet_concentration is None:
        outlet_concentration = 0.0
        captured_rate = 0.0
    else:
        mass_flow = inlet_n * capacity._molar_mass(mixture, list(inlet_z))
        outlet_concentration = inlet_concentration * (1.0 - efficiency)
        captured_rate = mass_flow * inlet_concentration * efficiency

    return FilterStates(
        clean_drop=clean_pa,
        applied_drop=applied,
        p_out=p_out,
        h_out=h_out,
        efficiency=efficiency,
        outlet_concentration=outlet_concentration,
        captured_rate_kg_s=captured_rate,
    )


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
