"""``process.separator_performance`` - NeqSim's droplet-based separation calculator.

Spec: ``specs/models/process/separator_performance.toml``. The Python twin of
``crates/azoth-process/src/kernels/separator_performance.rs`` over
``crates/azoth-process/src/models/separator_performance.rs``.

The class's standard path - the Rosin-Rammler/log-normal droplet distribution, the grade
efficiency curves, the Schiller-Naumann settling, and the six fractions and five
efficiencies they produce - **and its enhanced chain**: the flow-regime map, the inlet
device, the vessel geometry and the API 12J record.

**The class's quirks are reproduced**: the three-phase gas entrainment writes the whole
``1 - overall`` to both liquid fractions, the gas carry-under writes one fraction twice, and
the liquid-liquid stage uses a fixed 300 s. The efficiency integrals are the midpoint sums
over 50 log-spaced classes without renormalising, and ``erf``/``inverse_erf`` are the
class's own approximations rather than the ``math`` module's exact ones.
"""

from __future__ import annotations

import math
import sys
from collections.abc import Callable, Sequence
from dataclasses import dataclass
from itertools import pairwise

from azoth.core.errors import InvalidInputError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import SeparatorPerformanceResult
from azoth.core.units import Q, from_si, input_records_to_si, input_to_si
from azoth.core.warnings import Warning

MODEL_ID = "process.separator_performance"

_G = 9.81
_MAX_ITERATIONS = 50
_TOLERANCE = 1.0e-6
_NUMBER_OF_CLASSES = 50
_LIQUID_LIQUID_RESIDENCE_TIME = 300.0
_MIN_LIQUID_VELOCITY = 0.0001
_MIN_PHASE_FRACTION = 0.01

_WIRE_MESH = "wire_mesh"
_VANE_PACK = "vane_pack"
_AXIAL_CYCLONE = "axial_cyclone"
_PLATE_PACK = "plate_pack"
_CUSTOM = "custom"
_GRAVITY = "gravity"

#: The class's per-type defaults: `(d_50, sharpness, eta_max)`.
_DEFAULTS = {
    _WIRE_MESH: (5.0e-6, 2.5, 0.998),
    _VANE_PACK: (12.0e-6, 2.0, 0.995),
    _AXIAL_CYCLONE: (3.0e-6, 4.0, 0.998),
}
_PLATE_PACK_SHARPNESS = 2.0


@dataclass(frozen=True, slots=True)
class _Distribution:
    """``DropletSizeDistribution``'s closed form and its discrete classes."""

    kind: str
    characteristic_diameter: float
    spread: float

    def cumulative(self, diameter: float) -> float:
        """``cumulativeFraction``: the volume fraction smaller than ``diameter``."""
        if diameter <= 0.0:
            return 0.0
        if self.kind == "rosin_rammler":
            return 1.0 - math.exp(-((diameter / self.characteristic_diameter) ** self.spread))
        z = (math.log(diameter) - math.log(self.characteristic_diameter)) / self.spread
        return 0.5 * (1.0 + _erf(z / math.sqrt(2.0)))

    def inverse(self, fraction: float) -> float:
        """``inverseCDF``."""
        if fraction <= 0.0:
            return 0.0
        if fraction >= 1.0:
            return self.characteristic_diameter * 10.0
        if self.kind == "rosin_rammler":
            return float(
                self.characteristic_diameter * (-math.log(1.0 - fraction)) ** (1.0 / self.spread)
            )
        return math.exp(
            math.log(self.characteristic_diameter)
            + self.spread * math.sqrt(2.0) * _inverse_erf(2.0 * fraction - 1.0)
        )

    def classes(self) -> list[tuple[float, float, float]]:
        """``getDiscreteClasses``: ``(lower, midpoint, volume fraction)`` over 50 bins."""
        lower_bound = self.inverse(0.001)
        if lower_bound <= 0.0:
            lower_bound = self.characteristic_diameter * 0.01
        upper_bound = self.inverse(0.999)
        log_min = math.log(lower_bound)
        log_max = math.log(upper_bound)
        step = (log_max - log_min) / _NUMBER_OF_CLASSES
        out: list[tuple[float, float, float]] = []
        for index in range(_NUMBER_OF_CLASSES):
            lower = math.exp(log_min + index * step)
            upper = math.exp(log_min + (index + 1) * step)
            middle = math.exp(log_min + (index + 0.5) * step)
            out.append((lower, middle, self.cumulative(upper) - self.cumulative(lower)))
        return out


@dataclass(frozen=True, slots=True)
class _LookupPoint:
    """One row of a custom curve's table, in SI."""

    diameter: float
    efficiency: float


@dataclass(frozen=True, slots=True)
class _CurveSettings:
    """One declared internals curve, before its defaults are resolved."""

    internals: str
    cut_diameter: float | None
    sharpness: float | None
    max_efficiency: float | None
    points: tuple[_LookupPoint, ...] | None


@dataclass(frozen=True, slots=True)
class _Curve:
    """``GradeEfficiencyCurve``, resolved."""

    internals: str
    cut_diameter: float
    sharpness: float
    max_efficiency: float
    min_efficiency: float
    points: tuple[_LookupPoint, ...]

    @classmethod
    def gravity(cls, cut_diameter: float) -> _Curve:
        return cls(_GRAVITY, cut_diameter, 2.0, 1.0, 0.0, ())

    @classmethod
    def resolve(cls, role: str, settings: _CurveSettings) -> _Curve:
        """Resolve one declared curve, with the class's per-type defaults."""
        if settings.internals == _GRAVITY:
            raise InvalidInputError(
                role, "gravity is the settling section, not a configured internals curve"
            )
        if settings.internals == _CUSTOM:
            points = settings.points
            if points is None:
                raise InvalidInputError(
                    "points", f"{role} = custom needs a lookup table of at least two points"
                )
            if len(points) < 2:
                raise InvalidInputError(
                    "points", f"{role} = custom needs at least two points, not {len(points)}"
                )
            clamped = tuple(
                _LookupPoint(point.diameter, min(1.0, max(0.0, point.efficiency)))
                for point in points
            )
            return cls(_CUSTOM, 0.0, 2.0, 1.0, 0.0, clamped)
        if settings.internals == _PLATE_PACK:
            if settings.cut_diameter is None:
                raise InvalidInputError("cut_diameter", f"{role} = plate_pack needs a cut diameter")
            if settings.max_efficiency is None:
                raise InvalidInputError(
                    "max_efficiency", f"{role} = plate_pack needs a maximum efficiency"
                )
            return cls(
                _PLATE_PACK,
                settings.cut_diameter,
                _PLATE_PACK_SHARPNESS if settings.sharpness is None else settings.sharpness,
                settings.max_efficiency,
                0.0,
                (),
            )
        if settings.points is not None:
            raise InvalidInputError(
                "points", f"{role} = {settings.internals} does not read a lookup table"
            )
        default_cut, default_sharpness, default_max = _DEFAULTS[settings.internals]
        return cls(
            settings.internals,
            default_cut if settings.cut_diameter is None else settings.cut_diameter,
            default_sharpness if settings.sharpness is None else settings.sharpness,
            default_max if settings.max_efficiency is None else settings.max_efficiency,
            0.0,
            (),
        )

    def efficiency(self, diameter: float) -> float:
        """``getEfficiency``: the fractional capture at one diameter."""
        if diameter <= 0.0:
            return self.min_efficiency
        if self.internals == _CUSTOM:
            return _interpolate_custom(self.points, diameter)
        ratio = diameter / self.cut_diameter
        if self.internals == _GRAVITY:
            eta = ratio * ratio
        else:
            eta = self.max_efficiency * (1.0 - math.exp(-0.693 * ratio**self.sharpness))
        return min(self.max_efficiency, max(self.min_efficiency, eta))

    def overall(self, distribution: _Distribution) -> float:
        """``calcOverallEfficiency``: the class's midpoint sum, not renormalised."""
        return _sum_classes(distribution.classes(), lambda middle: self.efficiency(middle))


def _sum_classes(
    classes: Sequence[tuple[float, float, float]], capture: Callable[[float], float]
) -> float:
    total = sum(capture(middle) * fraction for _lower, middle, fraction in classes)
    return min(1.0, max(0.0, total))


def _combined(distribution: _Distribution, first: _Curve, second: _Curve) -> float:
    """``calcCombinedEfficiency``: the two stages in series, per size class."""

    def capture(middle: float) -> float:
        return 1.0 - (1.0 - first.efficiency(middle)) * (1.0 - second.efficiency(middle))

    return _sum_classes(distribution.classes(), capture)


def _interpolate_custom(points: Sequence[_LookupPoint], diameter: float) -> float:
    """``interpolateCustom``."""
    if len(points) < 2:
        return 0.0
    if diameter <= points[0].diameter:
        return points[0].efficiency
    if diameter >= points[-1].diameter:
        return points[-1].efficiency
    for first, second in pairwise(points):
        if diameter >= first.diameter and diameter < second.diameter:
            fraction = (diameter - first.diameter) / (second.diameter - first.diameter)
            return first.efficiency + fraction * (second.efficiency - first.efficiency)
    return points[-1].efficiency


def _erf(x: float) -> float:
    """Abramowitz and Stegun 7.1.26, the class's own approximation."""
    negative = x < 0.0
    absolute = abs(x)
    t = 1.0 / (1.0 + 0.3275911 * absolute)
    poly = t * (
        0.254829592 + t * (-0.284496736 + t * (1.421413741 + t * (-1.453152027 + t * 1.061405429)))
    )
    result = 1.0 - poly * math.exp(-absolute * absolute)
    return -result if negative else result


def _inverse_erf(x: float) -> float:
    """Winitzki (2008), the class's own approximation."""
    if x <= -1.0:
        return -6.0
    if x >= 1.0:
        return 6.0
    a = 0.147
    log_term = math.log(1.0 - x * x)
    pi_inverse = 2.0 / (math.pi * a)
    half_log = log_term / 2.0
    root_term = math.sqrt((pi_inverse + half_log) ** 2 - log_term / a)
    result = math.sqrt(root_term - pi_inverse - half_log)
    return result if x >= 0.0 else -result


def _drag_coefficient(reynolds: float) -> float:
    """``calcDragCoefficient``: Schiller-Naumann below 1000, Newton above."""
    if reynolds <= 0.0:
        return 1.0e10
    if reynolds < 1000.0:
        return float((24.0 / reynolds) * (1.0 + 0.15 * reynolds**0.687))
    return 0.44


def _reynolds(diameter: float, velocity: float, density: float, viscosity: float) -> float:
    if viscosity <= 0.0:
        return 0.0
    return density * abs(velocity) * diameter / viscosity


def _terminal_velocity(
    diameter: float, continuous: float, dispersed: float, viscosity: float
) -> float:
    """``calcTerminalVelocity``: positive when the dispersed phase is heavier."""
    if diameter <= 0.0 or continuous <= 0.0 or viscosity <= 0.0:
        return 0.0
    delta = abs(dispersed - continuous)
    if delta < 1.0e-6:
        return 0.0
    velocity = diameter * diameter * delta * _G / (18.0 * viscosity)
    for _ in range(_MAX_ITERATIONS):
        coefficient = _drag_coefficient(_reynolds(diameter, velocity, continuous, viscosity))
        updated = math.sqrt(4.0 * _G * diameter * delta / (3.0 * coefficient * continuous))
        if abs(updated - velocity) / (velocity + 1.0e-30) < _TOLERANCE:
            velocity = updated
            break
        velocity = 0.5 * (velocity + updated)
    return velocity if dispersed > continuous else -velocity


def _critical_diameter(
    available_height: float,
    residence_time: float,
    continuous: float,
    dispersed: float,
    viscosity: float,
) -> float:
    """``calcCriticalDiameter``: the smallest droplet that settles the height in the time."""
    if residence_time <= 0.0 or available_height <= 0.0:
        return sys.float_info.max
    required = available_height / residence_time
    delta = abs(dispersed - continuous)
    if delta < 1.0e-6:
        return sys.float_info.max
    diameter = math.sqrt(18.0 * viscosity * required / (delta * _G))
    for _ in range(_MAX_ITERATIONS):
        velocity = abs(_terminal_velocity(diameter, continuous, dispersed, viscosity))
        if velocity <= 0.0:
            break
        updated = diameter * math.sqrt(required / velocity)
        if abs(updated - diameter) / (diameter + 1.0e-30) < _TOLERANCE:
            diameter = updated
            break
        diameter = 0.5 * (diameter + updated)
    return diameter


def _estimate_liquid_residence_time(
    liquid_level_fraction: float, length: float, gas_velocity: float
) -> float:
    """``estimateLiquidResidenceTime``."""
    gas_fraction = max(_MIN_PHASE_FRACTION, 1.0 - liquid_level_fraction)
    liquid_fraction = max(_MIN_PHASE_FRACTION, liquid_level_fraction)
    liquid_velocity = gas_velocity * gas_fraction / liquid_fraction
    return length / max(liquid_velocity, _MIN_LIQUID_VELOCITY)


def _gas_section(
    orientation: str,
    diameter: float,
    length: float,
    liquid_level_fraction: float,
    gas_velocity: float,
) -> tuple[float, float]:
    gas_fraction = 1.0 - liquid_level_fraction
    height = length * gas_fraction if orientation == "vertical" else diameter * gas_fraction
    residence = length / gas_velocity if gas_fraction > 0.0 and gas_velocity > 0.0 else 0.0
    return height, residence


def _liquid_height(
    orientation: str, diameter: float, length: float, liquid_level_fraction: float
) -> float:
    if orientation == "vertical":
        return length * liquid_level_fraction
    return diameter * liquid_level_fraction


def _clamp01(value: float) -> float:
    return min(1.0, max(0.0, value))


def _distribution_kind(field: str, value: str) -> str:
    if value in ("rosin_rammler", "log_normal"):
        return value
    raise InvalidInputError(field, f"'{value}' is not one of the two declared distributions")


def _optional_distribution(
    role: str,
    distribution: str | None,
    characteristic_diameter: float | None,
    spread: float | None,
) -> _Distribution | None:
    if distribution is None and characteristic_diameter is None and spread is None:
        return None
    if distribution is None or characteristic_diameter is None or spread is None:
        raise InvalidInputError(
            f"{role}_distribution",
            "a distribution needs its type, characteristic diameter and spread together",
        )
    return _Distribution(
        _distribution_kind(f"{role}_distribution", distribution),
        characteristic_diameter,
        spread,
    )


def _lookup_points(rows: Sequence[dict[str, float]] | None) -> tuple[_LookupPoint, ...] | None:
    if rows is None:
        return None
    return tuple(_LookupPoint(float(row["diameter"]), float(row["efficiency"])) for row in rows)


def _segment_area(diameter: float, height: float) -> float:
    """`SeparatorGeometryCalculator.calcSegmentArea`."""
    if height <= 0.0 or diameter <= 0.0:
        return 0.0
    if height >= diameter:
        return math.pi * diameter * diameter / 4.0
    radius = diameter / 2.0
    theta = 2.0 * math.acos((radius - height) / radius)
    return radius * radius / 2.0 * (theta - math.sin(theta))


@dataclass(frozen=True, slots=True)
class _Geometry:
    gas_area: float
    liquid_area: float
    gas_volume: float
    liquid_volume: float
    effective_gas_settling_height: float
    effective_liquid_settling_height: float
    gas_residence_time: float
    liquid_residence_time: float
    gravity_section_length: float
    oil_pad_thickness: float
    water_layer_height: float


def _geometry_two_phase(
    diameter: float,
    length: float,
    horizontal: bool,
    level: float,
    gas_flow: float,
    liquid_flow: float,
) -> _Geometry:
    total_area = math.pi * diameter * diameter / 4.0
    if horizontal:
        liquid_height = level * diameter
        liquid_area = _segment_area(diameter, liquid_height)
        gas_area = total_area - liquid_area
        gravity_length = length * (0.90 - 0.05)
        gas_volume = gas_area * length
        liquid_volume = liquid_area * length
        gas_residence = (
            gravity_length / (gas_flow / gas_area) if gas_area > 0.0 and gas_flow > 0.0 else 0.0
        )
        liquid_residence = (
            liquid_volume / liquid_flow if liquid_volume > 0.0 and liquid_flow > 0.0 else 0.0
        )
        return _Geometry(
            gas_area,
            liquid_area,
            gas_volume,
            liquid_volume,
            diameter - liquid_height,
            liquid_height,
            gas_residence,
            liquid_residence,
            gravity_length,
            0.0,
            0.0,
        )
    liquid_height = level * length
    gas_height = length - liquid_height
    gravity_length = max(0.0, length * 0.90 - liquid_height)
    gas_volume = total_area * gas_height
    liquid_volume = total_area * liquid_height
    gas_residence = (
        gas_height / (gas_flow / total_area) if total_area > 0.0 and gas_flow > 0.0 else 0.0
    )
    liquid_residence = (
        liquid_volume / liquid_flow if liquid_volume > 0.0 and liquid_flow > 0.0 else 0.0
    )
    return _Geometry(
        total_area,
        total_area,
        gas_volume,
        liquid_volume,
        diameter,
        liquid_height,
        gas_residence,
        liquid_residence,
        gravity_length,
        0.0,
        0.0,
    )


def _geometry_three_phase(
    diameter: float,
    length: float,
    horizontal: bool,
    level: float,
    gas_flow: float,
    oil_flow: float,
    water_flow: float,
) -> _Geometry:
    if not horizontal:
        return _geometry_two_phase(diameter, length, False, level, gas_flow, oil_flow + water_flow)
    total_area = math.pi * diameter * diameter / 4.0
    liquid_height = level * diameter
    total_liquid = oil_flow + water_flow
    water_frac = water_flow / total_liquid if total_liquid > 0.0 else 0.5
    water_layer = liquid_height * water_frac
    oil_pad = liquid_height * (1.0 - water_frac)
    liquid_area = _segment_area(diameter, liquid_height)
    gas_area = total_area - liquid_area
    gravity_length = length * (0.60 - 0.05)
    oil_area = liquid_area - _segment_area(diameter, water_layer)
    gas_volume = gas_area * length
    liquid_volume = liquid_area * length
    gas_residence = (
        gravity_length / (gas_flow / gas_area) if gas_area > 0.0 and gas_flow > 0.0 else 0.0
    )
    oil_volume = oil_area * gravity_length
    liquid_residence = oil_volume / oil_flow if oil_flow > 0.0 and oil_volume > 0.0 else 0.0
    return _Geometry(
        gas_area,
        liquid_area,
        gas_volume,
        liquid_volume,
        diameter - liquid_height,
        oil_pad,
        gas_residence,
        liquid_residence,
        gravity_length,
        oil_pad,
        water_layer,
    )


def _flow_regime(
    horizontal: bool,
    vsg: float,
    vsl: float,
    rho_g: float,
    rho_l: float,
    mu_l: float,
    sigma: float,
    pipe_diameter: float,
) -> tuple[str, float, float]:
    """`MultiphaseFlowRegime.predict`, and its generated Rosin-Rammler parameters."""
    g = 9.81
    if horizontal:
        if vsl < 0.01:
            if vsg < 3.0:
                regime = "STRATIFIED_SMOOTH"
            elif vsg < 8.0:
                regime = "STRATIFIED_WAVY"
            else:
                regime = "ANNULAR"
        elif vsl > 3.0:
            regime = "ANNULAR_MIST" if vsg > 10.0 else "DISPERSED_BUBBLE"
        elif vsg < 0.3:
            regime = "PLUG"
        elif vsg < 3.0:
            regime = "SLUG"
        elif vsg < 15.0:
            regime = "STRATIFIED_WAVY" if vsl < 0.1 else "SLUG"
        else:
            regime = "ANNULAR" if vsl < 0.5 else "ANNULAR_MIST"
    else:
        delta = rho_l - rho_g
        bubble_slug = 0.25 + 0.35 * math.sqrt(g * pipe_diameter * delta / rho_l)
        slug_churn = math.sqrt(g * pipe_diameter)
        kutateladze = 3.1 * (sigma * g * delta / (rho_g * rho_g)) ** 0.25
        re_mix = rho_l * (vsg + vsl) * pipe_diameter / mu_l
        if vsl > 3.0 and re_mix > 8000.0:
            regime = "DISPERSED_BUBBLE"
        elif vsg < bubble_slug:
            regime = "BUBBLE"
        elif vsg < slug_churn:
            regime = "SLUG"
        elif vsg < kutateladze:
            regime = "CHURN"
        elif vsl < 0.1:
            regime = "ANNULAR"
        else:
            regime = "ANNULAR_MIST"
    d = pipe_diameter
    if sigma <= 0.0 or rho_g <= 0.0:
        return regime, 100.0e-6, 2.6
    weber = max(rho_g * vsg * vsg * d / sigma, 1.0)
    re_liquid = max(rho_l * vsl * d / mu_l, 1.0)
    if regime in ("STRATIFIED_SMOOTH", "STRATIFIED_WAVY"):
        d_max = min(max(weber**-0.6 * d, 50.0e-6), 5.0e-3)
        return regime, d_max / 2.5, 2.0
    if regime in ("ANNULAR", "ANNULAR_MIST"):
        re_liquid = max(rho_l * vsl * d / mu_l, 100.0)
        d32 = min(max(0.069 * weber**-0.6 * re_liquid**0.1 * d, 5.0e-6), 1.0e-3)
        return regime, d32 * 1.5, 2.8
    if regime == "SLUG":
        d_max = min(max(0.725 * weber**-0.6 * d, 30.0e-6), 3.0e-3)
        return regime, d_max / 3.0, 1.8
    if regime == "CHURN":
        d_max = min(max(0.5 * weber**-0.6 * d, 20.0e-6), 2.0e-3)
        return regime, d_max / 3.0, 2.2
    if regime == "PLUG":
        return regime, min(d * 0.05, 5.0e-3), 2.0
    v_mix = vsg + vsl
    weber_liquid = max(rho_l * v_mix * v_mix * d / sigma, 1.0)
    d_max = min(max(0.725 * weber_liquid**-0.6 * d, 100.0e-6), 10.0e-3)
    return regime, d_max / 2.0, 2.5


_INLET_PARAMETERS = {
    "none": (0.0, 1.0, 0.0, None),
    "deflector_plate": (0.60, 0.85, 0.5, 1500.0),
    "half_pipe": (0.70, 0.75, 1.0, 3500.0),
    "inlet_vane": (0.80, 0.55, 1.5, 6000.0),
    "inlet_cyclone": (0.95, 0.35, 3.0, 15000.0),
    "schoepentoeter": (0.85, 0.50, 2.5, 9000.0),
    "impingement_plate": (0.65, 0.80, 0.8, 2000.0),
}


def _inlet_device(
    device: str,
    incoming: tuple[float, float, bool],
    rho_g: float,
    rho_l: float,
    gas_flow: float,
    liquid_flow: float,
    nozzle_diameter: float,
) -> tuple[float, float, float, float, tuple[float, float, bool]]:
    """`InletDeviceModel.calculate`: bulk, dp, velocity, momentum, downstream DSD."""
    base, multiplier, coefficient, max_momentum = _INLET_PARAMETERS[device]
    area = math.pi * nozzle_diameter * nozzle_diameter / 4.0
    total = gas_flow + liquid_flow
    velocity = total / area if area > 0.0 else 0.0
    gas_fraction = gas_flow / total if total > 0.0 else 1.0
    mix_density = gas_fraction * rho_g + (1.0 - gas_fraction) * rho_l
    momentum = mix_density * velocity * velocity
    if max_momentum is None:
        bulk = 0.0
    elif momentum > max_momentum:
        degradation = min(1.0, (momentum - max_momentum) / max_momentum)
        bulk = base * (1.0 - 0.5 * degradation)
    elif momentum < max_momentum * 0.6 * 0.3:
        bulk = base * 0.5
    else:
        bulk = base
    pressure_drop = coefficient * 0.5 * mix_density * velocity * velocity
    if device == "none":
        downstream = incoming
    else:
        new_diameter = max(incoming[0] * multiplier, 1.0e-6)
        spread = incoming[1]
        if device in ("inlet_cyclone", "inlet_vane"):
            spread = min(spread * 1.2, 4.0)
        if device in ("half_pipe", "deflector_plate") and momentum > 3000.0:
            spread = max(spread * 0.8, 1.5)
        downstream = (new_diameter, spread, incoming[2])
    return bulk, pressure_drop, velocity, momentum, downstream


def _turbulence_corrected_cut(
    cut: float,
    velocity: float,
    height: float,
    k_factor: float,
    design_k: float,
    rho_g: float,
    rho_l: float,
    mu_g: float,
) -> float:
    if cut <= 0.0 or velocity <= 0.0 or height <= 0.0:
        return cut
    turbulence = 0.05
    if design_k > 0.0:
        load = min(1.0, k_factor / design_k)
        turbulence = 0.05 + 0.15 * load * load
    u_rms = turbulence * velocity
    eddy = height / 5.0
    diffusivity = u_rms * eddy
    terminal = abs(_terminal_velocity(cut, rho_g, rho_l, mu_g))
    if terminal <= 0.0:
        return cut
    argument = min(1.0 + diffusivity / (terminal * height), 81.0)
    return float(cut * argument**0.25)


def _enhanced_result(
    warnings: list[Warning],
    gas_rho: float,
    oil_rho: float,
    water_rho: float,
    gas_mu: float,
    oil_mu: float,
    water_mu: float,
    velocity: float,
    diameter: float,
    length: float,
    orientation: str,
    level: float,
    gas_liquid: _Distribution | None,
    gas_bubble: _Distribution | None,
    oil_in_water: _Distribution | None,
    water_in_oil: _Distribution | None,
    mist: _CurveSettings | None,
    coalescer: _CurveSettings | None,
    include_gravity_section: bool | None,
    liquid_in_gas_factor: float | None,
    gas_carry_under_factor: float | None,
    liquid_liquid_factor: float | None,
    inlet_pipe: float | None,
    sigma: float | None,
    oil_fraction: float | None,
    turbulence: bool | None,
    inlet_device: str | None,
    design_k: float | None,
    ll_residence: float | None,
) -> SeparatorPerformanceResult:
    """`calculateEnhanced`: the regime, device, geometry and API chain, in order."""
    mist_curve = None if mist is None else _Curve.resolve("mist_eliminator", mist)
    coalescer_curve = (
        None if coalescer is None else _Curve.resolve("oil_water_coalescer", coalescer)
    )
    horizontal = orientation == "horizontal"
    has_oil = oil_rho > 0.0
    has_water = water_rho > 0.0
    three_phase = has_oil and has_water
    liquid_density = oil_rho if has_oil else water_rho
    liquid_viscosity = oil_mu if has_oil else water_mu
    inlet_pipe = 0.2 if inlet_pipe is None else inlet_pipe
    sigma = 0.025 if sigma is None else sigma
    oil_fraction = 0.5 if oil_fraction is None else oil_fraction
    turbulence = True if turbulence is None else turbulence
    device = "half_pipe" if inlet_device is None else inlet_device
    if device not in _INLET_PARAMETERS:
        raise InvalidInputError(
            "inlet_device", f"'{device}' is not one of the seven declared devices"
        )

    vsl = velocity * level / (1.0 - level + 1.0e-10)
    if vsl < 0.001:
        vsl = 0.01
    regime, generated_d, generated_spread = _flow_regime(
        horizontal, velocity, vsl, gas_rho, liquid_density, liquid_viscosity, sigma, inlet_pipe
    )
    if gas_liquid is None:
        incoming = (generated_d, generated_spread, False)
    else:
        incoming = (
            gas_liquid.characteristic_diameter,
            gas_liquid.spread,
            gas_liquid.kind == "log_normal",
        )
    pipe_area = math.pi * inlet_pipe * inlet_pipe / 4.0
    bulk, pressure_drop, nozzle_velocity, momentum, downstream = _inlet_device(
        device,
        incoming,
        gas_rho,
        liquid_density,
        velocity * pipe_area,
        vsl * pipe_area,
        inlet_pipe,
    )
    gravity_dsd = _Distribution(
        kind="log_normal" if downstream[2] else "rosin_rammler",
        characteristic_diameter=downstream[0],
        spread=downstream[1],
    )
    vessel_area = math.pi * diameter * diameter / 4.0
    vessel_gas = velocity * vessel_area * (1.0 - level)
    vessel_liquid = vsl * vessel_area * level
    if three_phase:
        fraction = min(0.99, max(0.01, oil_fraction))
        geometry = _geometry_three_phase(
            diameter,
            length,
            horizontal,
            level,
            vessel_gas,
            vessel_liquid * fraction,
            vessel_liquid * (1.0 - fraction),
        )
    else:
        geometry = _geometry_two_phase(
            diameter, length, horizontal, level, vessel_gas, vessel_liquid
        )
    delta = liquid_density - gas_rho
    k_factor = velocity * math.sqrt(gas_rho / delta) if delta > 0.0 and gas_rho > 0.0 else 0.0
    if design_k is None:
        if mist is not None and mist.internals != _WIRE_MESH:
            raise InvalidInputError(
                "design_k_factor",
                "the internals database's design K is refused for this mist eliminator; "
                "state `design_k_factor`",
            )
        design_k_value = 0.107
    else:
        design_k_value = design_k
    k_utilization = k_factor / design_k_value if design_k_value > 0.0 else 0.0
    flooded = k_utilization > 1.0

    cut = _critical_diameter(
        geometry.effective_gas_settling_height,
        geometry.gas_residence_time,
        gas_rho,
        liquid_density,
        gas_mu,
    )
    effective_cut = cut
    if turbulence:
        effective_cut = _turbulence_corrected_cut(
            cut,
            velocity,
            geometry.effective_gas_settling_height,
            k_factor,
            design_k_value,
            gas_rho,
            liquid_density,
            gas_mu,
        )
    gravity_curve = _Curve.gravity(effective_cut)
    include_gravity = True if include_gravity_section is None else include_gravity_section
    gravity_efficiency = gravity_curve.overall(gravity_dsd) if include_gravity else 0.0
    mist_efficiency = 0.0
    overall = 0.0
    if mist_curve is not None:
        if flooded:
            penalty = max(0.0, 1.0 - (k_utilization - 1.0))
            mist_efficiency = mist_curve.overall(gravity_dsd) * penalty * 0.5
            overall = gravity_efficiency + (1.0 - gravity_efficiency) * mist_efficiency
        else:
            mist_efficiency = mist_curve.overall(gravity_dsd)
            overall = (
                _combined(gravity_dsd, gravity_curve, mist_curve)
                if include_gravity
                else mist_efficiency
            )
    elif include_gravity:
        overall = gravity_efficiency

    api_gas = False
    api_liquid = False
    api_cut_um = 0.0
    api_k_utilization = 0.0
    api_max_k = 0.0
    api_min_residence = 0.0
    api_gas_liquid_comment = None
    api_liquid_comment = None
    if cut > 0.0:
        liquid_residence = geometry.liquid_residence_time
        if liquid_residence <= 0.0:
            liquid_residence = _estimate_liquid_residence_time(level, length, velocity)
        api_max_k = 0.25 if mist_curve is not None else (0.120 if horizontal else 0.107)
        api_k_utilization = k_factor / api_max_k
        api_cut_um = cut * 1.0e6
        api_gas = cut <= 100.0e-6 and k_factor <= api_max_k
        api_min_residence = 300.0 if three_phase else 180.0
        api_liquid = liquid_residence >= api_min_residence
        gas_prefix = "COMPLIANT" if api_gas else "NON-COMPLIANT"
        liquid_prefix = "COMPLIANT" if api_liquid else "NON-COMPLIANT"
        api_gas_liquid_comment = (
            f"{gas_prefix}: cut={api_cut_um:.0f} um (max 100 um), "
            f"K={k_factor:.3f} m/s (max {api_max_k:.3f} m/s)"
        )
        api_liquid_comment = (
            f"{liquid_prefix}: liquid residence time = {liquid_residence:.0f} s "
            f"(min {api_min_residence:.0f} s)"
        )

    combined_efficiency = bulk + (1.0 - bulk) * overall
    overall = min(1.0, combined_efficiency)
    entrainment = 1.0 - overall
    oil_in_gas = 0.0
    water_in_gas = 0.0
    if three_phase:
        oil_in_gas = entrainment
        water_in_gas = entrainment
    elif has_oil:
        oil_in_gas = entrainment
    else:
        water_in_gas = entrainment

    gas_in_oil = 0.0
    gas_in_water = 0.0
    if gas_bubble is not None:
        height = geometry.effective_liquid_settling_height
        residence = (
            120.0 if geometry.liquid_residence_time <= 0.0 else geometry.liquid_residence_time
        )
        bubble_cut = _critical_diameter(
            height, residence, liquid_density, gas_rho, liquid_viscosity
        )
        efficiency = _Curve.gravity(bubble_cut).overall(gas_bubble)
        gas_in_oil = 1.0 - efficiency
        gas_in_water = gas_in_oil

    oil_in_water_fraction = 0.0
    water_in_oil_fraction = 0.0
    ll_gravity_efficiency = 0.0
    if three_phase:
        residence = geometry.liquid_residence_time if ll_residence is None else ll_residence
        if residence <= 0.0:
            residence = 300.0
        oil_pad = (
            diameter * level / 2.0
            if geometry.oil_pad_thickness <= 0.0
            else geometry.oil_pad_thickness
        )
        if water_in_oil is not None:
            water_cut = _critical_diameter(oil_pad, residence, oil_rho, water_rho, oil_mu)
            efficiency = _Curve.gravity(water_cut).overall(water_in_oil)
            water_in_oil_fraction = (
                (1.0 - efficiency) * (1.0 - coalescer_curve.overall(water_in_oil))
                if coalescer_curve is not None
                else 1.0 - efficiency
            )
        if oil_in_water is not None:
            water_layer = (
                oil_pad if geometry.water_layer_height <= 0.0 else geometry.water_layer_height
            )
            oil_cut = _critical_diameter(water_layer, residence, water_rho, oil_rho, water_mu)
            efficiency = _Curve.gravity(oil_cut).overall(oil_in_water)
            ll_gravity_efficiency = efficiency
            oil_in_water_fraction = (
                (1.0 - efficiency) * (1.0 - coalescer_curve.overall(oil_in_water))
                if coalescer_curve is not None
                else 1.0 - efficiency
            )

    liquid_in_gas_factor = 1.0 if liquid_in_gas_factor is None else max(0.0, liquid_in_gas_factor)
    gas_carry_under_factor = (
        1.0 if gas_carry_under_factor is None else max(0.0, gas_carry_under_factor)
    )
    liquid_liquid_factor = 1.0 if liquid_liquid_factor is None else max(0.0, liquid_liquid_factor)
    return SeparatorPerformanceResult(
        oil_in_gas=_clamp01(oil_in_gas * liquid_in_gas_factor),
        water_in_gas=_clamp01(water_in_gas * liquid_in_gas_factor),
        gas_in_oil=_clamp01(gas_in_oil * gas_carry_under_factor),
        gas_in_water=_clamp01(gas_in_water * gas_carry_under_factor),
        oil_in_water=_clamp01(oil_in_water_fraction * liquid_liquid_factor),
        water_in_oil=_clamp01(water_in_oil_fraction * liquid_liquid_factor),
        gravity_section_efficiency=gravity_efficiency,
        mist_eliminator_efficiency=mist_efficiency,
        overall_gas_liquid_efficiency=overall,
        liquid_liquid_gravity_efficiency=ll_gravity_efficiency,
        gravity_cut_diameter=from_si(cut, "m"),
        inlet_flow_regime=regime,
        inlet_device_bulk_efficiency=bulk,
        inlet_device_pressure_drop=from_si(pressure_drop, "Pa"),
        inlet_nozzle_velocity=from_si(nozzle_velocity, "m/s"),
        inlet_momentum_flux=from_si(momentum, "Pa"),
        k_factor=from_si(k_factor, "m/s"),
        k_factor_utilization=k_utilization,
        mist_eliminator_flooded=flooded,
        gas_area=from_si(geometry.gas_area, "m**2"),
        liquid_area=from_si(geometry.liquid_area, "m**2"),
        effective_gas_settling_height=from_si(geometry.effective_gas_settling_height, "m"),
        effective_liquid_settling_height=from_si(geometry.effective_liquid_settling_height, "m"),
        gas_residence_time=from_si(geometry.gas_residence_time, "s"),
        liquid_residence_time=from_si(geometry.liquid_residence_time, "s"),
        gravity_section_length=from_si(geometry.gravity_section_length, "m"),
        oil_pad_thickness=from_si(geometry.oil_pad_thickness, "m"),
        water_layer_height=from_si(geometry.water_layer_height, "m"),
        api_gas_liquid_compliant=api_gas,
        api_liquid_compliant=api_liquid,
        api_fully_compliant=api_gas and api_liquid,
        api_gas_liquid_comment=api_gas_liquid_comment,
        api_liquid_comment=api_liquid_comment,
        api_gravity_cut_diameter=from_si(api_cut_um / 1.0e6, "m"),
        api_k_factor_utilization=api_k_utilization,
        api_max_k_factor=from_si(api_max_k, "m/s"),
        api_min_liquid_residence_time=from_si(api_min_residence, "s"),
        post_inlet_device_characteristic_diameter=from_si(downstream[0], "m"),
        post_inlet_device_spread=downstream[1],
        warnings=tuple(warnings),
    )


def separator_performance(
    gas_density: Q,
    oil_density: Q,
    water_density: Q,
    gas_viscosity: Q,
    oil_viscosity: Q,
    water_viscosity: Q,
    gas_velocity: Q,
    vessel_diameter: Q,
    vessel_length: Q,
    orientation: str,
    liquid_level_fraction: float,
    gas_liquid_distribution: str | None = None,
    gas_liquid_characteristic_diameter: Q | None = None,
    gas_liquid_spread: float | None = None,
    gas_bubble_distribution: str | None = None,
    gas_bubble_characteristic_diameter: Q | None = None,
    gas_bubble_spread: float | None = None,
    oil_in_water_distribution: str | None = None,
    oil_in_water_characteristic_diameter: Q | None = None,
    oil_in_water_spread: float | None = None,
    water_in_oil_distribution: str | None = None,
    water_in_oil_characteristic_diameter: Q | None = None,
    water_in_oil_spread: float | None = None,
    mist_eliminator: str | None = None,
    mist_eliminator_cut_diameter: Q | None = None,
    mist_eliminator_sharpness: float | None = None,
    mist_eliminator_max_efficiency: float | None = None,
    mist_eliminator_points: list[dict[str, Q]] | None = None,
    oil_water_coalescer: str | None = None,
    oil_water_coalescer_cut_diameter: Q | None = None,
    oil_water_coalescer_max_efficiency: float | None = None,
    oil_water_coalescer_points: list[dict[str, Q]] | None = None,
    include_gravity_section: bool | None = None,
    liquid_in_gas_calibration_factor: float | None = None,
    gas_carry_under_calibration_factor: float | None = None,
    liquid_liquid_calibration_factor: float | None = None,
    enhanced: bool | None = None,
    inlet_pipe_diameter: Q | None = None,
    surface_tension: Q | None = None,
    oil_volume_fraction: float | None = None,
    apply_turbulence_correction: bool | None = None,
    inlet_device: str | None = None,
    design_k_factor: Q | None = None,
    liquid_liquid_residence_time: Q | None = None,
) -> SeparatorPerformanceResult:
    """The class's standard separation calculation.

    Args:
        gas_density: the gas phase's mass density; zero is the class's absent phase.
        oil_density: the oil phase's mass density; zero is absent.
        water_density: the water phase's mass density; zero is absent.
        gas_viscosity: the gas phase's dynamic viscosity.
        oil_viscosity: the oil phase's dynamic viscosity.
        water_viscosity: the water phase's dynamic viscosity.
        gas_velocity: the superficial gas velocity in the gas section.
        vessel_diameter: the vessel's internal diameter.
        vessel_length: the vessel's tangent-to-tangent length.
        orientation: ``horizontal`` or ``vertical``.
        liquid_level_fraction: the liquid's share of the cross-section or height.
        gas_liquid_distribution: the entrained liquid's distribution; absent leaves the
            gas-liquid branch unrun.
        gas_liquid_characteristic_diameter: its ``d_63.2`` or ``d_50``.
        gas_liquid_spread: its spread.
        gas_bubble_distribution: the entrained gas's distribution.
        gas_bubble_characteristic_diameter: its characteristic diameter.
        gas_bubble_spread: its spread.
        oil_in_water_distribution: the dispersed oil's distribution in water.
        oil_in_water_characteristic_diameter: its characteristic diameter.
        oil_in_water_spread: its spread.
        water_in_oil_distribution: the dispersed water's distribution in oil.
        water_in_oil_characteristic_diameter: its characteristic diameter.
        water_in_oil_spread: its spread.
        mist_eliminator: wire mesh, vane pack, axial cyclone, plate pack or custom.
        mist_eliminator_cut_diameter: ``d_50``; absent takes the type's default.
        mist_eliminator_sharpness: the sigmoid's ``n``; absent takes the type's default.
        mist_eliminator_max_efficiency: ``eta_max``; absent takes the type's default.
        mist_eliminator_points: a custom curve's table.
        oil_water_coalescer: plate pack or custom.
        oil_water_coalescer_cut_diameter: the plate pack's ``d_50``.
        oil_water_coalescer_max_efficiency: the plate pack's maximum efficiency.
        oil_water_coalescer_points: a custom coalescer's table.
        include_gravity_section: whether the gravity stage is included; default true.
        liquid_in_gas_calibration_factor: the multiplier on the two gas entrainments.
        gas_carry_under_calibration_factor: the multiplier on the two carry-unders.
        liquid_liquid_calibration_factor: the multiplier on the two liquid-liquid fractions.

    Returns:
        The six fractions, the four efficiencies and the gravity cut diameter.

    Raises:
        InvalidInputError: for a value outside its declared range, an undeclared name, a
            distribution whose fields are not stated together, a curve missing what its
            type requires, or a lookup table on a type that does not read one.

    See :func:`azoth.process.separator_performance`.
    """
    from azoth._models_gen import model

    spec = model(MODEL_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    gas_rho = input_to_si(spec, "gas_density", gas_density)
    oil_rho = input_to_si(spec, "oil_density", oil_density)
    water_rho = input_to_si(spec, "water_density", water_density)
    gas_mu = input_to_si(spec, "gas_viscosity", gas_viscosity)
    oil_mu = input_to_si(spec, "oil_viscosity", oil_viscosity)
    water_mu = input_to_si(spec, "water_viscosity", water_viscosity)
    velocity = input_to_si(spec, "gas_velocity", gas_velocity)
    diameter = input_to_si(spec, "vessel_diameter", vessel_diameter)
    length = input_to_si(spec, "vessel_length", vessel_length)

    def optional(name: str, value: Q | None) -> float | None:
        return None if value is None else input_to_si(spec, name, value)

    gas_liquid_diameter = optional(
        "gas_liquid_characteristic_diameter", gas_liquid_characteristic_diameter
    )
    gas_bubble_diameter = optional(
        "gas_bubble_characteristic_diameter", gas_bubble_characteristic_diameter
    )
    oil_in_water_diameter = optional(
        "oil_in_water_characteristic_diameter", oil_in_water_characteristic_diameter
    )
    water_in_oil_diameter = optional(
        "water_in_oil_characteristic_diameter", water_in_oil_characteristic_diameter
    )
    mist_cut = optional("mist_eliminator_cut_diameter", mist_eliminator_cut_diameter)
    coalescer_cut = optional("oil_water_coalescer_cut_diameter", oil_water_coalescer_cut_diameter)
    inlet_pipe = optional("inlet_pipe_diameter", inlet_pipe_diameter)
    sigma = optional("surface_tension", surface_tension)
    design_k = optional("design_k_factor", design_k_factor)
    ll_residence = optional("liquid_liquid_residence_time", liquid_liquid_residence_time)
    oil_fraction = None if oil_volume_fraction is None else float(oil_volume_fraction)

    apply_checks(
        checks.on_input,
        {
            "gas_density": gas_rho,
            "oil_density": oil_rho,
            "water_density": water_rho,
            "gas_viscosity": gas_mu,
            "oil_viscosity": oil_mu,
            "water_viscosity": water_mu,
            "gas_velocity": velocity,
            "vessel_diameter": diameter,
            "vessel_length": length,
            "liquid_level_fraction": liquid_level_fraction,
            "gas_liquid_characteristic_diameter": gas_liquid_diameter,
            "gas_liquid_spread": gas_liquid_spread,
            "gas_bubble_characteristic_diameter": gas_bubble_diameter,
            "gas_bubble_spread": gas_bubble_spread,
            "oil_in_water_characteristic_diameter": oil_in_water_diameter,
            "oil_in_water_spread": oil_in_water_spread,
            "water_in_oil_characteristic_diameter": water_in_oil_diameter,
            "water_in_oil_spread": water_in_oil_spread,
            "mist_eliminator_cut_diameter": mist_cut,
            "mist_eliminator_sharpness": mist_eliminator_sharpness,
            "mist_eliminator_max_efficiency": mist_eliminator_max_efficiency,
            "oil_water_coalescer_cut_diameter": coalescer_cut,
            "oil_water_coalescer_max_efficiency": oil_water_coalescer_max_efficiency,
            "inlet_pipe_diameter": inlet_pipe,
            "surface_tension": sigma,
            "oil_volume_fraction": oil_fraction,
            "design_k_factor": design_k,
            "liquid_liquid_residence_time": ll_residence,
        }.get,
        warnings,
    )

    if orientation not in ("horizontal", "vertical"):
        raise InvalidInputError(
            "orientation", f"'{orientation}' is not one of the two declared orientations"
        )

    gas_liquid = _optional_distribution(
        "gas_liquid", gas_liquid_distribution, gas_liquid_diameter, gas_liquid_spread
    )
    gas_bubble = _optional_distribution(
        "gas_bubble", gas_bubble_distribution, gas_bubble_diameter, gas_bubble_spread
    )
    oil_in_water = _optional_distribution(
        "oil_in_water", oil_in_water_distribution, oil_in_water_diameter, oil_in_water_spread
    )
    water_in_oil = _optional_distribution(
        "water_in_oil", water_in_oil_distribution, water_in_oil_diameter, water_in_oil_spread
    )

    mist = _curve_settings(
        "mist_eliminator",
        {"none", _WIRE_MESH, _VANE_PACK, _AXIAL_CYCLONE, _PLATE_PACK, _CUSTOM},
        mist_eliminator,
        mist_cut,
        mist_eliminator_sharpness,
        mist_eliminator_max_efficiency,
        None
        if mist_eliminator_points is None
        else _lookup_points(
            input_records_to_si(spec, "mist_eliminator_points", mist_eliminator_points)
        ),
    )
    coalescer = _curve_settings(
        "oil_water_coalescer",
        {"none", _PLATE_PACK, _CUSTOM},
        oil_water_coalescer,
        coalescer_cut,
        None,
        oil_water_coalescer_max_efficiency,
        None
        if oil_water_coalescer_points is None
        else _lookup_points(
            input_records_to_si(spec, "oil_water_coalescer_points", oil_water_coalescer_points)
        ),
    )
    mist_curve = None if mist is None else _Curve.resolve("mist_eliminator", mist)
    coalescer_curve = (
        None if coalescer is None else _Curve.resolve("oil_water_coalescer", coalescer)
    )

    if enhanced:
        return _enhanced_result(
            warnings,
            gas_rho,
            oil_rho,
            water_rho,
            gas_mu,
            oil_mu,
            water_mu,
            velocity,
            diameter,
            length,
            orientation,
            liquid_level_fraction,
            gas_liquid,
            gas_bubble,
            oil_in_water,
            water_in_oil,
            mist,
            coalescer,
            include_gravity_section,
            liquid_in_gas_calibration_factor,
            gas_carry_under_calibration_factor,
            liquid_liquid_calibration_factor,
            inlet_pipe,
            sigma,
            oil_fraction,
            apply_turbulence_correction,
            inlet_device,
            design_k,
            ll_residence,
        )

    has_liquid = oil_rho > 0.0 or water_rho > 0.0
    liquid_density = oil_rho if oil_rho > 0.0 else water_rho
    liquid_viscosity = oil_mu if oil_rho > 0.0 else water_mu
    oil_in_gas = 0.0
    water_in_gas = 0.0
    gas_in_oil = 0.0
    gas_in_water = 0.0
    oil_in_water_fraction = 0.0
    water_in_oil_fraction = 0.0
    gravity_section_efficiency = 0.0
    mist_efficiency = 0.0
    overall = 0.0
    liquid_liquid_gravity_efficiency = 0.0
    gravity_cut_diameter = 0.0

    if gas_liquid is not None and has_liquid:
        height, residence = _gas_section(
            orientation, diameter, length, liquid_level_fraction, velocity
        )
        gravity_cut_diameter = _critical_diameter(
            height, residence, gas_rho, liquid_density, gas_mu
        )
        gravity = _Curve.gravity(gravity_cut_diameter)
        include_gravity = True if include_gravity_section is None else include_gravity_section
        if include_gravity:
            gravity_section_efficiency = gravity.overall(gas_liquid)
        if mist_curve is not None:
            mist_efficiency = mist_curve.overall(gas_liquid)
        if include_gravity and mist_curve is not None:
            overall = _combined(gas_liquid, gravity, mist_curve)
        elif include_gravity:
            overall = gravity_section_efficiency
        elif mist_curve is not None:
            overall = mist_efficiency
        entrainment = 1.0 - overall
        if oil_rho > 0.0 and water_rho > 0.0:
            oil_in_gas = entrainment
            water_in_gas = entrainment
        elif oil_rho > 0.0:
            oil_in_gas = entrainment
        else:
            water_in_gas = entrainment

    if gas_bubble is not None and has_liquid:
        bubble_height = _liquid_height(orientation, diameter, length, liquid_level_fraction)
        bubble_residence = (
            0.0
            if bubble_height <= 0.0
            else _estimate_liquid_residence_time(liquid_level_fraction, length, velocity)
        )
        cut = _critical_diameter(
            bubble_height, bubble_residence, liquid_density, gas_rho, liquid_viscosity
        )
        gas_in_oil = 1.0 - _Curve.gravity(cut).overall(gas_bubble)
        gas_in_water = gas_in_oil

    if oil_rho > 0.0 and water_rho > 0.0:
        settling_height = _liquid_height(orientation, diameter, length, liquid_level_fraction) / 2.0
        if water_in_oil is not None:
            cut = _critical_diameter(
                settling_height,
                _LIQUID_LIQUID_RESIDENCE_TIME,
                oil_rho,
                water_rho,
                oil_mu,
            )
            efficiency = _Curve.gravity(cut).overall(water_in_oil)
            water_in_oil_fraction = 1.0 - efficiency
            if coalescer_curve is not None:
                water_in_oil_fraction *= 1.0 - coalescer_curve.overall(water_in_oil)
        if oil_in_water is not None:
            cut = _critical_diameter(
                settling_height,
                _LIQUID_LIQUID_RESIDENCE_TIME,
                water_rho,
                oil_rho,
                water_mu,
            )
            efficiency = _Curve.gravity(cut).overall(oil_in_water)
            liquid_liquid_gravity_efficiency = efficiency
            oil_in_water_fraction = 1.0 - efficiency
            if coalescer_curve is not None:
                oil_in_water_fraction *= 1.0 - coalescer_curve.overall(oil_in_water)

    # The class's setters are `Math.max(0.0, factor)`.
    liquid_in_gas_factor = (
        1.0
        if liquid_in_gas_calibration_factor is None
        else max(0.0, liquid_in_gas_calibration_factor)
    )
    gas_carry_under_factor = (
        1.0
        if gas_carry_under_calibration_factor is None
        else max(0.0, gas_carry_under_calibration_factor)
    )
    liquid_liquid_factor = (
        1.0
        if liquid_liquid_calibration_factor is None
        else max(0.0, liquid_liquid_calibration_factor)
    )
    oil_in_gas = _clamp01(oil_in_gas * liquid_in_gas_factor)
    water_in_gas = _clamp01(water_in_gas * liquid_in_gas_factor)
    gas_in_oil = _clamp01(gas_in_oil * gas_carry_under_factor)
    gas_in_water = _clamp01(gas_in_water * gas_carry_under_factor)
    oil_in_water_fraction = _clamp01(oil_in_water_fraction * liquid_liquid_factor)
    water_in_oil_fraction = _clamp01(water_in_oil_fraction * liquid_liquid_factor)

    return SeparatorPerformanceResult(
        oil_in_gas=oil_in_gas,
        water_in_gas=water_in_gas,
        gas_in_oil=gas_in_oil,
        gas_in_water=gas_in_water,
        oil_in_water=oil_in_water_fraction,
        water_in_oil=water_in_oil_fraction,
        gravity_section_efficiency=gravity_section_efficiency,
        mist_eliminator_efficiency=mist_efficiency,
        overall_gas_liquid_efficiency=overall,
        liquid_liquid_gravity_efficiency=liquid_liquid_gravity_efficiency,
        gravity_cut_diameter=from_si(gravity_cut_diameter, "m"),
        inlet_flow_regime=None,
        inlet_device_bulk_efficiency=None,
        inlet_device_pressure_drop=None,
        inlet_nozzle_velocity=None,
        inlet_momentum_flux=None,
        k_factor=None,
        k_factor_utilization=None,
        mist_eliminator_flooded=None,
        gas_area=None,
        liquid_area=None,
        effective_gas_settling_height=None,
        effective_liquid_settling_height=None,
        gas_residence_time=None,
        liquid_residence_time=None,
        gravity_section_length=None,
        oil_pad_thickness=None,
        water_layer_height=None,
        api_gas_liquid_compliant=None,
        api_liquid_compliant=None,
        api_fully_compliant=None,
        api_gas_liquid_comment=None,
        api_liquid_comment=None,
        api_gravity_cut_diameter=None,
        api_k_factor_utilization=None,
        api_max_k_factor=None,
        api_min_liquid_residence_time=None,
        post_inlet_device_characteristic_diameter=None,
        post_inlet_device_spread=None,
        warnings=tuple(warnings),
    )


def _curve_settings(
    role: str,
    allowed: set[str],
    internals: str | None,
    cut_diameter: float | None,
    sharpness: float | None,
    max_efficiency: float | None,
    points: tuple[_LookupPoint, ...] | None,
) -> _CurveSettings | None:
    if internals is None or internals == "none":
        if (
            cut_diameter is not None
            or sharpness is not None
            or max_efficiency is not None
            or points is not None
        ):
            raise InvalidInputError(role, "curve parameters are stated without an internals type")
        return None
    if internals not in allowed:
        raise InvalidInputError(
            role, f"'{internals}' is not one of the declared internals for {role}"
        )
    return _CurveSettings(internals, cut_diameter, sharpness, max_efficiency, points)


__all__ = ["separator_performance"]
