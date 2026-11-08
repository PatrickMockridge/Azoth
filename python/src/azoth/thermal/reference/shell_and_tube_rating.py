"""``thermal.shell_and_tube_rating`` - the shell-and-tube thermal-hydraulic rating.

Spec: ``specs/calcs/thermal/shell_and_tube_rating.toml``

The Python twin of ``crates/azoth-thermal/src/shell_and_tube_rating.rs``, and both mirror
``ThermalDesignCalculator``. **The port implements ``calculateStrict()``, not ``calculate()``**:
a geometry the class would leave at zero is refused here, because a zero coefficient is not
an answer.

# What the class carries that the name does not say

``Js`` is computed with the central baffle spacing three times, so it is 1.0; the class never
varies the inlet or outlet spacing. ``tubeHeating`` is set by the equipment path and read by no
arithmetic here. And the Darcy friction factor is the class's own helper - ``64/Re`` below 2300
and the smooth Swamee-Jain form above - rather than the library's Colebrook.
"""

from __future__ import annotations

import math

from azoth._registry_gen import spec as _spec_for
from azoth.core.errors import InvalidInputError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import ShellAndTubeRatingResult
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning

CALC_ID = "thermal.shell_and_tube_rating"

KERN = "kern"
BELL_DELAWARE = "bell_delaware"


def darcy_friction(re: float) -> float:
    """The class's ``calcDarcyFriction``: ``64/Re`` below 2300, Swamee-Jain above."""
    if re <= 0.0:
        return 0.0
    if re < 2300.0:
        return 64.0 / re
    log_term = math.log10(5.74 / math.pow(re, 0.9))
    return 0.25 / (log_term * log_term)


def gnielinski_nusselt(re: float, pr: float, friction_factor: float) -> float:
    """``HeatTransferCoefficientCalculator.calculateGnielinskiNusselt``."""
    if re < 3000.0 or pr < 0.5:
        return laminar_nusselt()
    f = friction_factor
    numerator = (f / 8.0) * (re - 1000.0) * pr
    denominator = 1.0 + 12.7 * math.sqrt(f / 8.0) * (math.pow(pr, 2.0 / 3.0) - 1.0)
    return max(3.66, numerator / denominator)


def laminar_nusselt() -> float:
    """``calculateLaminarNusselt(true)``, constant wall temperature."""
    return 3.66


def _whole(value: float, field: str) -> int:
    if not math.isfinite(value) or value <= 0.0 or value != math.floor(value):
        raise InvalidInputError(field, f"{value} is not a positive whole count")
    return int(value)


def _optional_whole(value: float | None, field: str, default: int) -> int:
    if value is None:
        return default
    if not math.isfinite(value) or value < 0.0 or value != math.floor(value):
        raise InvalidInputError(field, f"{value} is not a non-negative whole count")
    return int(value)


def _required(raw: dict[str, float | None], name: str) -> float:
    value = raw[name]
    if value is None:
        raise InvalidInputError(name, "a required input is absent")
    return value


def shell_and_tube_rating(
    tube_outside_diameter: Q,
    tube_inside_diameter: Q,
    tube_length: Q,
    tube_count: float,
    tube_passes: float,
    tube_pitch: Q,
    triangular_pitch: bool,
    shell_inside_diameter: Q,
    baffle_spacing: Q,
    baffle_count: float,
    baffle_cut: float,
    tube_wall_conductivity: Q,
    tube_density: Q,
    tube_viscosity: Q,
    tube_cp: Q,
    tube_conductivity: Q,
    tube_mass_flow: Q,
    shell_density: Q,
    shell_viscosity: Q,
    shell_cp: Q,
    shell_conductivity: Q,
    shell_mass_flow: Q,
    fouling_tube: float,
    fouling_shell: float,
    shell_side_method: str,
    tube_fouling_layer_thickness: Q | None = None,
    tube_fouling_layer_conductivity: Q | None = None,
    shell_viscosity_wall: Q | None = None,
    has_sealing: bool | None = None,
    sealing_pairs: float | None = None,
    tube_to_baffle_clearance: Q | None = None,
    shell_to_baffle_clearance: Q | None = None,
    bypass_area: Q | None = None,
) -> ShellAndTubeRatingResult:
    """The thermal-hydraulic rating of a fixed shell-and-tube geometry.

    Raises:
        InvalidInputError: for every state the class's ``validateStrictInputs`` refuses, and
            for a ``shell_side_method`` that is not ``kern`` or ``bell_delaware``.
        OutOfRangeError: for an input outside the spec's declared bounds.
    """
    spec = _spec_for(CALC_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    raw = {
        "tube_outside_diameter": input_to_si(spec, "tube_outside_diameter", tube_outside_diameter),
        "tube_inside_diameter": input_to_si(spec, "tube_inside_diameter", tube_inside_diameter),
        "tube_length": input_to_si(spec, "tube_length", tube_length),
        "tube_count": tube_count,
        "tube_passes": tube_passes,
        "tube_pitch": input_to_si(spec, "tube_pitch", tube_pitch),
        "shell_inside_diameter": input_to_si(spec, "shell_inside_diameter", shell_inside_diameter),
        "baffle_spacing": input_to_si(spec, "baffle_spacing", baffle_spacing),
        "baffle_count": baffle_count,
        "baffle_cut": baffle_cut,
        "tube_wall_conductivity": input_to_si(
            spec, "tube_wall_conductivity", tube_wall_conductivity
        ),
        "tube_density": input_to_si(spec, "tube_density", tube_density),
        "tube_viscosity": input_to_si(spec, "tube_viscosity", tube_viscosity),
        "tube_cp": input_to_si(spec, "tube_cp", tube_cp),
        "tube_conductivity": input_to_si(spec, "tube_conductivity", tube_conductivity),
        "tube_mass_flow": input_to_si(spec, "tube_mass_flow", tube_mass_flow),
        "shell_density": input_to_si(spec, "shell_density", shell_density),
        "shell_viscosity": input_to_si(spec, "shell_viscosity", shell_viscosity),
        "shell_cp": input_to_si(spec, "shell_cp", shell_cp),
        "shell_conductivity": input_to_si(spec, "shell_conductivity", shell_conductivity),
        "shell_mass_flow": input_to_si(spec, "shell_mass_flow", shell_mass_flow),
        "fouling_tube": fouling_tube,
        "fouling_shell": fouling_shell,
        "tube_fouling_layer_thickness": None
        if tube_fouling_layer_thickness is None
        else input_to_si(spec, "tube_fouling_layer_thickness", tube_fouling_layer_thickness),
        "tube_fouling_layer_conductivity": None
        if tube_fouling_layer_conductivity is None
        else input_to_si(spec, "tube_fouling_layer_conductivity", tube_fouling_layer_conductivity),
        "shell_viscosity_wall": None
        if shell_viscosity_wall is None
        else input_to_si(spec, "shell_viscosity_wall", shell_viscosity_wall),
        "sealing_pairs": sealing_pairs,
        "tube_to_baffle_clearance": None
        if tube_to_baffle_clearance is None
        else input_to_si(spec, "tube_to_baffle_clearance", tube_to_baffle_clearance),
        "shell_to_baffle_clearance": None
        if shell_to_baffle_clearance is None
        else input_to_si(spec, "shell_to_baffle_clearance", shell_to_baffle_clearance),
        "bypass_area": None
        if bypass_area is None
        else input_to_si(spec, "bypass_area", bypass_area),
    }
    apply_checks(checks.on_input, raw.get, warnings)

    _validate_strict(raw, shell_side_method, has_sealing, triangular_pitch)

    tube_count_i = _whole(_required(raw, "tube_count"), "tube_count")
    tube_passes_i = _whole(_required(raw, "tube_passes"), "tube_passes")
    baffle_count_i = _whole(_required(raw, "baffle_count"), "baffle_count")
    sealing_pairs_i = _optional_whole(raw["sealing_pairs"], "sealing_pairs", 0)
    sealed = bool(has_sealing)
    thickness = raw["tube_fouling_layer_thickness"] or 0.0
    deposit_k = raw["tube_fouling_layer_conductivity"]
    deposit_k = 1.0 if deposit_k is None else deposit_k
    tube_to_baffle = raw["tube_to_baffle_clearance"]
    tube_to_baffle = 0.0004 if tube_to_baffle is None else tube_to_baffle
    shell_to_baffle = raw["shell_to_baffle_clearance"]
    shell_to_baffle = 0.003 if shell_to_baffle is None else shell_to_baffle
    bypass = raw["bypass_area"] or 0.0

    tube_od = _required(raw, "tube_outside_diameter")
    tube_id = _required(raw, "tube_inside_diameter")
    tube_length_si = _required(raw, "tube_length")
    tube_pitch_si = _required(raw, "tube_pitch")
    shell_id = _required(raw, "shell_inside_diameter")
    baffle_spacing_si = _required(raw, "baffle_spacing")
    baffle_cut = _required(raw, "baffle_cut")
    tube_density_si = _required(raw, "tube_density")
    tube_viscosity_si = _required(raw, "tube_viscosity")
    tube_cp_si = _required(raw, "tube_cp")
    tube_k = _required(raw, "tube_conductivity")
    tube_mass_flow_si = _required(raw, "tube_mass_flow")
    shell_density_si = _required(raw, "shell_density")
    shell_viscosity_si = _required(raw, "shell_viscosity")
    shell_cp_si = _required(raw, "shell_cp")
    shell_k = _required(raw, "shell_conductivity")
    shell_mass_flow_si = _required(raw, "shell_mass_flow")

    effective_id = tube_id - 2.0 * thickness
    outside_area = math.pi * tube_od * tube_length_si * tube_count_i

    area_per_tube = math.pi / 4.0 * effective_id * effective_id
    tubes_per_pass = tube_count_i // tube_passes_i
    tube_flow_area = tubes_per_pass * area_per_tube
    tube_velocity = tube_mass_flow_si / (tube_density_si * tube_flow_area)
    tube_re = tube_density_si * tube_velocity * effective_id / tube_viscosity_si
    tube_pr = tube_cp_si * tube_viscosity_si / tube_k
    if tube_re > 3000.0:
        tube_nu = gnielinski_nusselt(tube_re, tube_pr, darcy_friction(tube_re))
    elif tube_re > 2300.0:
        laminar = laminar_nusselt()
        turbulent = gnielinski_nusselt(3000.0, tube_pr, darcy_friction(3000.0))
        fraction = (tube_re - 2300.0) / 700.0
        tube_nu = laminar + fraction * (turbulent - laminar)
    else:
        tube_nu = laminar_nusselt()
    tube_htc = tube_nu * tube_k / effective_id

    stated_wall = raw["shell_viscosity_wall"]
    wall_viscosity = (
        stated_wall if stated_wall is not None and stated_wall > 0.0 else shell_viscosity_si
    )
    shell_area = _crossflow_area(shell_id, baffle_spacing_si, tube_od, tube_pitch_si)
    shell_mass_flux = shell_mass_flow_si / shell_area
    shell_velocity = shell_mass_flux / shell_density_si

    if shell_side_method == BELL_DELAWARE:
        shell_re = shell_mass_flux * tube_od / shell_viscosity_si
        ideal = _ideal_crossflow_htc(
            shell_mass_flux,
            tube_od,
            shell_viscosity_si,
            shell_cp_si,
            shell_k,
            wall_viscosity,
            triangular_pitch,
        )
        jc = _jc(baffle_cut)
        jl = _jl(
            tube_to_baffle, shell_to_baffle, shell_area, tube_count_i, tube_od, baffle_spacing_si
        )
        if bypass <= 0.0:
            bypass = _bypass_area(shell_id, shell_id * 0.9, baffle_spacing_si)
        rows = _tube_rows_crossflow(shell_id, baffle_cut, tube_pitch_si, triangular_pitch)
        jb = _jb(bypass, shell_area, sealed, sealing_pairs_i, rows)
        js = _js(baffle_spacing_si, baffle_spacing_si, baffle_spacing_si, baffle_count_i)
        jr = _jr(shell_re, rows)
        shell_htc = ideal * jc * jl * jb * js * jr
    else:
        equivalent_diameter = _shell_equivalent_diameter(tube_od, tube_pitch_si, triangular_pitch)
        shell_re = shell_mass_flux * equivalent_diameter / shell_viscosity_si
        shell_htc = _kern_shell_side_htc(
            shell_mass_flux,
            equivalent_diameter,
            shell_viscosity_si,
            shell_cp_si,
            shell_k,
            wall_viscosity,
        )
        jc = jl = jb = js = jr = 1.0

    do_by_di = tube_od / tube_id
    total_resistance = (
        1.0 / shell_htc
        + _required(raw, "fouling_shell")
        + tube_od * math.log(do_by_di) / (2.0 * _required(raw, "tube_wall_conductivity"))
        + do_by_di * _required(raw, "fouling_tube")
        + tube_od / effective_id / tube_htc
    )
    deposit_resistance = 0.0
    if thickness > 0.0:
        deposit_resistance = tube_od * math.log(tube_id / effective_id) / (2.0 * deposit_k)
        total_resistance += deposit_resistance
    overall_u = 1.0 / total_resistance

    friction = darcy_friction(tube_re)
    tube_dp_friction = (
        friction
        * tube_length_si
        / effective_id
        * tube_density_si
        * tube_velocity
        * tube_velocity
        / 2.0
    )
    tube_dp_return = 2.5 * tube_density_si * tube_velocity * tube_velocity / 2.0
    tube_dp = tube_passes_i * (tube_dp_friction + tube_dp_return)

    if shell_side_method == BELL_DELAWARE:
        rows = _tube_rows_crossflow(shell_id, baffle_cut, tube_pitch_si, triangular_pitch)
        ideal_dp = _ideal_crossflow_pressure_drop(
            rows, shell_mass_flux, shell_density_si, shell_re, triangular_pitch
        )
        rl = _rl(
            tube_to_baffle, shell_to_baffle, shell_area, tube_count_i, tube_od, baffle_spacing_si
        )
        if bypass <= 0.0:
            bypass = _bypass_area(shell_id, shell_id * 0.9, baffle_spacing_si)
        rb = _rb(bypass, shell_area, sealed, sealing_pairs_i, rows)
        window_area = math.pi / 4.0 * shell_id * shell_id * baffle_cut * 0.5
        rows_window = max(1, rows // 3)
        window_dp = _window_pressure_drop(
            window_area, shell_mass_flow_si, shell_density_si, rows_window
        )
        shell_dp = (
            ideal_dp * (baffle_count_i - 1) * rl * rb
            + baffle_count_i * window_dp * rl
            + 2.0 * ideal_dp * (1.0 + rows / 5.0) * rb
        )
    else:
        equivalent_diameter = _shell_equivalent_diameter(tube_od, tube_pitch_si, triangular_pitch)
        shell_dp = _kern_shell_side_pressure_drop(
            shell_mass_flux,
            equivalent_diameter,
            shell_id,
            baffle_count_i,
            shell_density_si,
            shell_viscosity_si,
            wall_viscosity,
        )

    return ShellAndTubeRatingResult(
        tube_side_heat_transfer_coefficient=from_si(tube_htc, "W/(m**2*K)"),
        shell_side_heat_transfer_coefficient=from_si(shell_htc, "W/(m**2*K)"),
        overall_heat_transfer_coefficient=from_si(overall_u, "W/(m**2*K)"),
        tube_side_pressure_drop=from_si(tube_dp, "Pa"),
        shell_side_pressure_drop=from_si(shell_dp, "Pa"),
        tube_side_velocity=from_si(tube_velocity, "m/s"),
        shell_side_velocity=from_si(shell_velocity, "m/s"),
        tube_side_reynolds=tube_re,
        shell_side_reynolds=shell_re,
        effective_tube_inside_diameter=from_si(effective_id, "m"),
        outside_heat_transfer_area=from_si(outside_area, "m**2"),
        tube_fouling_layer_resistance_outside=deposit_resistance,
        jc=jc,
        jl=jl,
        jb=jb,
        js=js,
        jr=jr,
        warnings=tuple(warnings),
    )


def _validate_strict(
    raw: dict[str, float | None],
    shell_side_method: str,
    has_sealing: bool | None,
    triangular_pitch: bool,
) -> None:
    if shell_side_method not in (KERN, BELL_DELAWARE):
        raise InvalidInputError(
            "shell_side_method",
            f"{shell_side_method} is not one of {KERN!r} or {BELL_DELAWARE!r}",
        )
    positives = [
        _required(raw, "tube_outside_diameter"),
        _required(raw, "tube_inside_diameter"),
        _required(raw, "tube_length"),
        _required(raw, "tube_pitch"),
        _required(raw, "shell_inside_diameter"),
        _required(raw, "baffle_spacing"),
        _required(raw, "tube_wall_conductivity"),
        _required(raw, "tube_density"),
        _required(raw, "tube_viscosity"),
        _required(raw, "tube_cp"),
        _required(raw, "tube_conductivity"),
        _required(raw, "tube_mass_flow"),
        _required(raw, "shell_density"),
        _required(raw, "shell_viscosity"),
        _required(raw, "shell_cp"),
        _required(raw, "shell_conductivity"),
        _required(raw, "shell_mass_flow"),
    ]
    if any(value is None or not math.isfinite(value) or value <= 0.0 for value in positives):
        raise InvalidInputError(
            "geometry",
            "strict thermal rating requires positive finite geometry and fluid inputs",
        )
    if (
        _required(raw, "tube_outside_diameter") <= _required(raw, "tube_inside_diameter")
        or _required(raw, "tube_pitch") <= _required(raw, "tube_outside_diameter")
        or _required(raw, "shell_inside_diameter") <= _required(raw, "tube_outside_diameter")
        or _required(raw, "baffle_spacing") > _required(raw, "tube_length")
    ):
        raise InvalidInputError(
            "geometry",
            "strict thermal rating requires OD > ID, pitch > OD, shell ID > OD, and baffle "
            "spacing no greater than tube length",
        )
    tube_count = _whole(_required(raw, "tube_count"), "tube_count")
    tube_passes = _whole(_required(raw, "tube_passes"), "tube_passes")
    baffle_count = _whole(_required(raw, "baffle_count"), "baffle_count")
    if tube_count % tube_passes != 0:
        raise InvalidInputError(
            "tube_count",
            "strict thermal rating requires the tube count to divide evenly by the passes",
        )
    if baffle_count < 1:
        raise InvalidInputError(
            "baffle_count", "strict thermal rating requires at least one baffle"
        )
    baffle_cut = _required(raw, "baffle_cut")
    if not math.isfinite(baffle_cut) or baffle_cut <= 0.0 or baffle_cut >= 0.5:
        raise InvalidInputError(
            "baffle_cut",
            "strict thermal rating requires a baffle cut between zero and one half",
        )
    if (
        not math.isfinite(_required(raw, "fouling_tube"))
        or _required(raw, "fouling_tube") < 0.0
        or not math.isfinite(_required(raw, "fouling_shell"))
        or _required(raw, "fouling_shell") < 0.0
    ):
        raise InvalidInputError(
            "fouling", "strict thermal rating requires finite nonnegative fouling resistances"
        )
    stated_wall = raw["shell_viscosity_wall"]
    if stated_wall is not None:
        value = stated_wall
        if not math.isfinite(value) or value == 0.0:
            raise InvalidInputError(
                "shell_viscosity_wall",
                "strict thermal rating requires a finite nonzero wall viscosity; negative "
                "selects the bulk viscosity",
            )
    thickness = raw["tube_fouling_layer_thickness"] or 0.0
    conductivity = raw["tube_fouling_layer_conductivity"]
    conductivity = 1.0 if conductivity is None else conductivity
    if not math.isfinite(thickness) or thickness < 0.0:
        raise InvalidInputError(
            "tube_fouling_layer_thickness",
            "the deposit thickness must be finite and nonnegative",
        )
    if not math.isfinite(conductivity) or conductivity <= 0.0:
        raise InvalidInputError(
            "tube_fouling_layer_conductivity",
            "the deposit conductivity must be finite and positive",
        )
    effective_id = _required(raw, "tube_inside_diameter") - 2.0 * thickness
    if effective_id <= 0.0:
        raise InvalidInputError(
            "tube_fouling_layer_thickness",
            "the deposit must leave a positive tube flow diameter",
        )
    tube_flow_area = math.pi / 4.0 * effective_id * effective_id * (tube_count / tube_passes)
    shell_area = _crossflow_area(
        _required(raw, "shell_inside_diameter"),
        _required(raw, "baffle_spacing"),
        _required(raw, "tube_outside_diameter"),
        _required(raw, "tube_pitch"),
    )
    equivalent_diameter = _shell_equivalent_diameter(
        _required(raw, "tube_outside_diameter"), _required(raw, "tube_pitch"), triangular_pitch
    )
    outside_area = (
        math.pi
        * _required(raw, "tube_outside_diameter")
        * _required(raw, "tube_length")
        * tube_count
    )
    if any(
        not math.isfinite(value) or value <= 0.0
        for value in (tube_flow_area, shell_area, equivalent_diameter, outside_area)
    ):
        raise InvalidInputError(
            "geometry",
            "strict thermal rating requires finite positive flow and heat-transfer geometry",
        )
    if shell_side_method == BELL_DELAWARE:
        tube_to_baffle = raw["tube_to_baffle_clearance"]
        tube_to_baffle = 0.0004 if tube_to_baffle is None else tube_to_baffle
        shell_to_baffle = raw["shell_to_baffle_clearance"]
        shell_to_baffle = 0.003 if shell_to_baffle is None else shell_to_baffle
        bypass = raw["bypass_area"] or 0.0
        sealing_pairs = _optional_whole(raw["sealing_pairs"], "sealing_pairs", 0)
        sealed = bool(has_sealing)
        if (
            not math.isfinite(tube_to_baffle)
            or tube_to_baffle < 0.0
            or tube_to_baffle
            >= _required(raw, "tube_pitch") - _required(raw, "tube_outside_diameter")
            or not math.isfinite(shell_to_baffle)
            or shell_to_baffle < 0.0
            or shell_to_baffle >= _required(raw, "shell_inside_diameter")
            or not math.isfinite(bypass)
            or bypass < 0.0
            or bypass > _required(raw, "shell_inside_diameter") * _required(raw, "baffle_spacing")
            or (sealed and sealing_pairs == 0)
        ):
            raise InvalidInputError(
                "bell_delaware",
                "strict Bell-Delaware rating requires valid clearances, bypass area and "
                "sealing-strip counts",
            )


def _ideal_crossflow_htc(
    mass_flux: float,
    tube_od: float,
    mu: float,
    cp: float,
    k: float,
    mu_wall: float,
    triangular: bool,
) -> float:
    if mu <= 0.0 or k <= 0.0 or cp <= 0.0 or tube_od <= 0.0:
        return 0.0
    re = mass_flux * tube_od / mu
    pr = cp * mu / k
    pr_wall = cp * mu_wall / k
    viscosity_correction = math.pow(pr / pr_wall, 0.25) if pr_wall > 0.0 else 1.0
    if triangular:
        if re < 100.0:
            nu = 0.9 * math.pow(re, 0.4) * math.pow(pr, 0.36) * viscosity_correction
        elif re < 1000.0:
            nu = 0.52 * math.pow(re, 0.5) * math.pow(pr, 0.36) * viscosity_correction
        elif re < 200_000.0:
            nu = 0.27 * math.pow(re, 0.63) * math.pow(pr, 0.36) * viscosity_correction
        else:
            nu = 0.033 * math.pow(re, 0.8) * math.pow(pr, 0.36) * viscosity_correction
    elif re < 100.0:
        nu = 0.9 * math.pow(re, 0.4) * math.pow(pr, 0.36) * viscosity_correction
    elif re < 1000.0:
        nu = 0.52 * math.pow(re, 0.5) * math.pow(pr, 0.36) * viscosity_correction
    elif re < 200_000.0:
        nu = 0.27 * math.pow(re, 0.63) * math.pow(pr, 0.36) * viscosity_correction
    else:
        nu = 0.021 * math.pow(re, 0.84) * math.pow(pr, 0.36) * viscosity_correction
    return nu * k / tube_od


def _kern_shell_side_htc(
    mass_flux: float,
    equivalent_diameter: float,
    mu: float,
    cp: float,
    k: float,
    mu_wall: float,
) -> float:
    if mu <= 0.0 or k <= 0.0 or cp <= 0.0 or equivalent_diameter <= 0.0:
        return 0.0
    re = mass_flux * equivalent_diameter / mu
    pr = cp * mu / k
    viscosity_correction = math.pow(mu / mu_wall, 0.14) if mu_wall > 0.0 else 1.0
    nu = 0.36 * math.pow(re, 0.55) * math.pow(pr, 1.0 / 3.0) * viscosity_correction
    return nu * k / equivalent_diameter


def _kern_shell_side_pressure_drop(
    mass_flux: float,
    equivalent_diameter: float,
    shell_diameter: float,
    baffle_count: int,
    rho: float,
    mu: float,
    mu_wall: float,
) -> float:
    if mu <= 0.0 or rho <= 0.0 or equivalent_diameter <= 0.0:
        return 0.0
    re = mass_flux * equivalent_diameter / mu
    viscosity_correction = math.pow(mu / mu_wall, 0.14) if mu_wall > 0.0 else 1.0
    if re < 10.0:
        friction = 64.0 / max(re, 0.1)
    elif re < 1000.0:
        friction = 5.6 * math.pow(re, -0.344)
    else:
        friction = 0.228 * math.pow(re, -0.1574)
    velocity = mass_flux / rho
    return (
        friction
        * shell_diameter
        * (baffle_count + 1.0)
        * rho
        * velocity
        * velocity
        / (2.0 * equivalent_diameter * viscosity_correction)
    )


def _shell_equivalent_diameter(tube_od: float, tube_pitch: float, triangular: bool) -> float:
    if tube_od <= 0.0 or tube_pitch <= 0.0:
        return 0.0
    if triangular:
        flow_area = (
            math.sqrt(3.0) / 4.0 * tube_pitch * tube_pitch - math.pi / 8.0 * tube_od * tube_od
        )
        wetted_perimeter = math.pi / 2.0 * tube_od
    else:
        flow_area = tube_pitch * tube_pitch - math.pi / 4.0 * tube_od * tube_od
        wetted_perimeter = math.pi * tube_od
    return 4.0 * flow_area / wetted_perimeter


def _crossflow_area(
    shell_id: float, baffle_spacing: float, tube_od: float, tube_pitch: float
) -> float:
    if shell_id <= 0.0 or baffle_spacing <= 0.0 or tube_pitch <= 0.0:
        return 0.0
    return baffle_spacing * shell_id * (1.0 - tube_od / tube_pitch)


def _crossflow_fraction(baffle_cut: float) -> float:
    if baffle_cut <= 0.15:
        return 0.9
    if baffle_cut >= 0.45:
        return 0.3
    return 0.9 - (baffle_cut - 0.15) / (0.45 - 0.15) * 0.6


def _jc(baffle_cut: float) -> float:
    return 0.55 + 0.72 * _crossflow_fraction(baffle_cut)


def _jl(
    tube_to_baffle: float,
    shell_to_baffle: float,
    crossflow_area: float,
    tube_count: int,
    tube_od: float,
    baffle_spacing: float,
) -> float:
    if crossflow_area <= 0.0:
        return 0.7
    tube_baffle_area = tube_count * math.pi * tube_od * tube_to_baffle
    shell_baffle_area = math.pi * baffle_spacing * shell_to_baffle
    total_ratio = (tube_baffle_area + shell_baffle_area) / crossflow_area
    shell_fraction = shell_baffle_area / max(tube_baffle_area + shell_baffle_area, 1e-10)
    value = 0.44 * (1.0 - shell_fraction) + (1.0 - 0.44 * (1.0 - shell_fraction)) * math.exp(
        -2.2 * total_ratio
    )
    return min(max(value, 0.2), 1.0)


def _jb(
    bypass_area: float,
    crossflow_area: float,
    has_sealing: bool,
    sealing_pairs: int,
    tube_rows_crossflow: int,
) -> float:
    if crossflow_area <= 0.0:
        return 0.8
    fbp = bypass_area / crossflow_area
    if has_sealing and tube_rows_crossflow > 0 and sealing_pairs / tube_rows_crossflow >= 0.5:
        return 1.0
    return min(max(math.exp(-1.35 * fbp), 0.5), 1.0)


def _js(
    central_spacing: float, inlet_spacing: float, outlet_spacing: float, baffle_count: int
) -> float:
    if central_spacing <= 0.0 or baffle_count < 1:
        return 1.0
    nb = float(baffle_count)
    inlet_ratio = central_spacing / inlet_spacing
    outlet_ratio = central_spacing / outlet_spacing
    numerator = nb - 1.0 + math.pow(inlet_ratio, 0.6) + math.pow(outlet_ratio, 0.6)
    denominator = nb - 1.0 + inlet_ratio + outlet_ratio
    if denominator <= 0.0:
        return 1.0
    return min(max(numerator / denominator, 0.5), 1.0)


def _jr(re: float, tube_rows_crossflow: int) -> float:
    if re >= 100.0 or tube_rows_crossflow <= 10:
        return 1.0
    laminar = math.pow(10.0 / tube_rows_crossflow, 0.18)
    if re < 20.0:
        return min(max(laminar, 0.4), 1.0)
    fraction = (re - 20.0) / 80.0
    value = laminar + fraction * (1.0 - laminar)
    return min(max(value, 0.4), 1.0)


def _rl(
    tube_to_baffle: float,
    shell_to_baffle: float,
    crossflow_area: float,
    tube_count: int,
    tube_od: float,
    baffle_spacing: float,
) -> float:
    if crossflow_area <= 0.0:
        return 0.5
    tube_baffle_area = tube_count * math.pi * tube_od * tube_to_baffle
    shell_baffle_area = math.pi * baffle_spacing * shell_to_baffle
    total_ratio = (tube_baffle_area + shell_baffle_area) / crossflow_area
    shell_fraction = shell_baffle_area / max(tube_baffle_area + shell_baffle_area, 1e-10)
    value = math.exp(-1.33 * (1.0 + shell_fraction) * math.pow(total_ratio, 0.8))
    return min(max(value, 0.1), 1.0)


def _rb(
    bypass_area: float,
    crossflow_area: float,
    has_sealing: bool,
    sealing_pairs: int,
    tube_rows_crossflow: int,
) -> float:
    if crossflow_area <= 0.0:
        return 0.7
    fbp = bypass_area / crossflow_area
    if has_sealing and tube_rows_crossflow > 0 and sealing_pairs / tube_rows_crossflow >= 0.5:
        return 1.0
    return min(max(math.exp(-3.7 * fbp), 0.3), 1.0)


def _ideal_crossflow_pressure_drop(
    tube_rows_crossflow: int, mass_flux: float, rho: float, re: float, triangular: bool
) -> float:
    if rho <= 0.0 or mass_flux <= 0.0:
        return 0.0
    if triangular:
        friction = 48.0 / max(re, 0.1) if re < 1000.0 else 0.35 * math.pow(re, -0.15)
    else:
        friction = 35.0 / max(re, 0.1) if re < 1000.0 else 0.18 * math.pow(re, -0.15)
    return friction * tube_rows_crossflow * mass_flux * mass_flux / (2.0 * rho)


def _window_pressure_drop(
    window_area: float, mass_flow: float, rho: float, tube_rows_window: int
) -> float:
    if window_area <= 0.0 or rho <= 0.0:
        return 0.0
    velocity = mass_flow / (rho * window_area)
    return (2.0 + 0.6 * tube_rows_window) * rho * velocity * velocity / 2.0


def _tube_rows_crossflow(
    shell_id: float, baffle_cut: float, tube_pitch: float, triangular: bool
) -> int:
    if shell_id <= 0.0 or tube_pitch <= 0.0:
        return 1
    distance_between_tips = shell_id * (1.0 - 2.0 * baffle_cut)
    row_pitch = tube_pitch * math.sqrt(3.0) / 2.0 if triangular else tube_pitch
    return max(1, int(distance_between_tips / row_pitch))


def _bypass_area(shell_id: float, bundle_diameter: float, baffle_spacing: float) -> float:
    if shell_id <= 0.0 or bundle_diameter <= 0.0 or baffle_spacing <= 0.0:
        return 0.0
    return (shell_id - bundle_diameter) * baffle_spacing
