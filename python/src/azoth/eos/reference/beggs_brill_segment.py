"""``eos.beggs_brill_segment`` - NeqSim's ``PipeBeggsAndBrills`` segment correlation.

One pipe segment's Beggs-and-Brill flow regime, liquid holdup with its inclination
correction, no-slip density and viscosity, Reynolds number, friction factor and pressure
drop - without the flash or the march, which are the model layer's.

Spec: ``specs/calcs/eos/beggs_brill_segment.toml``.
"""

from __future__ import annotations

import math

from azoth.core.range import apply_checks, checks_for
from azoth.core.result import BeggsBrillSegmentResult
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning

CALC_ID = "eos.beggs_brill_segment"

#: The class's field-unit constants.
FEET_PER_METER = 3.2808399
FT3_PER_M3 = FEET_PER_METER**3
KG_M3_PER_LB_FT3 = 16.018_463_37
PROFILE_LB_FT3_TO_KG_M3 = 16.018_46
BAR_PER_IMPERIAL_DROP = 1.487_27e-5
CP_PER_PA_S = 1000.0
FROUDE_G = 32.174
HYDROSTATIC_G = 32.2
POUNDS_PER_SQUARE_INCH = 6894.757_293_168
MINIMUM_SURFACE_TENSION_DYNE_CM = 1.0


def _liquid_velocity_number(
    liquid_density_lb_ft3: float,
    liquid_velocity_ft_s: float,
    temperature_k: float,
    pressure_pa: float,
) -> float:
    specific_gravity = liquid_density_lb_ft3 / (1000.0 * 0.062_427_960_6)
    api = 141.5 / specific_gravity - 131.5
    sigma68 = 39.0 - 0.2571 * api
    sigma100 = 37.5 - 0.2571 * api
    temperature_f = (temperature_k - 273.15) * 9.0 / 5.0 + 32.0
    sigma: float
    if temperature_f > 100.0:
        sigma = sigma100
    elif temperature_f < 68.0:
        sigma = sigma68
    else:
        sigma = sigma68 + (temperature_f - 68.0) * (sigma100 - sigma68) / (100.0 - 68.0)
    pressure_psi = pressure_pa / POUNDS_PER_SQUARE_INCH
    sigma *= 1.0 - 0.024 * pressure_psi**0.45
    sigma = (
        MINIMUM_SURFACE_TENSION_DYNE_CM
        if math.isnan(sigma)
        else max(sigma, MINIMUM_SURFACE_TENSION_DYNE_CM)
    )
    return float(1.938 * liquid_velocity_ft_s * (liquid_density_lb_ft3 / sigma) ** 0.25)


def _inclination_corrected_holdup(
    regime: str,
    lam: float,
    froude: float,
    nvl: float,
    sin18: float,
    elevation_ft: float,
) -> float:
    horizontal: float
    if regime == "segregated":
        horizontal = 0.98 * lam**0.4846 / froude**0.0868
    elif regime == "intermittent":
        horizontal = 0.845 * lam**0.5351 / froude**0.0173
    else:
        horizontal = 1.065 * lam**0.5824 / froude**0.0609

    log_arg: float
    if elevation_ft > 0.0:
        if regime == "segregated":
            log_arg = 0.011 * nvl**3.539 / (lam**3.768 * froude**1.614)
        elif regime == "intermittent":
            log_arg = 2.96 * lam**0.305 * froude**0.0978 / nvl**0.4473
        else:
            log_arg = math.nan
    else:
        log_arg = 4.70 * nvl**0.1244 / (lam**0.3692 * froude**0.5056)

    betta: float = (1.0 - lam) * math.log(log_arg) if log_arg > 0.0 else 0.0
    betta = max(0.0, betta)
    return horizontal * (1.0 + betta * (sin18 - sin18**3 / 3.0))


def beggs_brill_segment(
    gas_flow: Q,
    liquid_flow: Q,
    gas_density: Q,
    liquid_density: Q,
    gas_viscosity: Q,
    liquid_viscosity: Q,
    temperature: Q,
    pressure: Q,
    diameter: Q,
    length: Q,
    elevation: Q,
    roughness: Q,
) -> BeggsBrillSegmentResult:
    """One pipe segment's Beggs-and-Brill state.

    Args:
        gas_flow: the gas phase's volumetric flow.
        liquid_flow: the liquid's volumetric flow, a single phase.
        gas_density: the gas phase's density.
        liquid_density: the liquid's density.
        gas_viscosity: the gas phase's dynamic viscosity.
        liquid_viscosity: the liquid's dynamic viscosity.
        temperature: the segment's absolute temperature.
        pressure: the segment's absolute pressure.
        diameter: the pipe's internal diameter.
        length: the segment's length.
        elevation: the segment's elevation change, positive uphill.
        roughness: the pipe wall's absolute roughness.

    Returns:
        The regime, holdup, velocities, no-slip properties, friction and pressure drop.

    Raises:
        OutOfRangeError: if an input is outside the spec's hard range.
        InvalidInputError: if neither phase flows.

    Example:
        >>> import azoth
        >>> q = azoth.ureg.Quantity
        >>> r = beggs_brill_segment(
        ...     q(0.01796694168385081, "m**3/s"), q(1.42993961618865e-4, "m**3/s"),
        ...     q(69.58154569795788, "kg/m**3"), q(970.1200910245294, "kg/m**3"),
        ...     q(1.5089488997328844e-5, "Pa*s"), q(4.6535425555414334e-4, "Pa*s"),
        ...     q(333.15, "K"), q(9999960.703327656, "Pa"), q(0.1, "m"), q(1.0, "m"),
        ...     q(0.0, "m"), q(1.0e-5, "m"))
        >>> r.regime
        'segregated'
    """
    from azoth._registry_gen import spec as _spec_for

    spec = _spec_for(CALC_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    values: dict[str, float] = {
        name: input_to_si(spec, name, value)
        for name, value in {
            "gas_flow": gas_flow,
            "liquid_flow": liquid_flow,
            "gas_density": gas_density,
            "liquid_density": liquid_density,
            "gas_viscosity": gas_viscosity,
            "liquid_viscosity": liquid_viscosity,
            "temperature": temperature,
            "pressure": pressure,
            "diameter": diameter,
            "length": length,
            "roughness": roughness,
        }.items()
    }
    values["elevation"] = input_to_si(spec, "elevation", elevation)
    apply_checks(checks.on_input, values.get, warnings)

    if values["gas_flow"] <= 0.0 and values["liquid_flow"] <= 0.0:
        from azoth.core.errors import InvalidInputError

        raise InvalidInputError(
            "gas_flow",
            "at least one phase must flow; the superficial velocities are undefined otherwise",
        )

    d_ft = values["diameter"] * FEET_PER_METER
    l_ft = values["length"] * FEET_PER_METER
    e_ft = values["elevation"] * FEET_PER_METER
    area = math.pi / 4.0 * d_ft * d_ft
    v_sl = values["liquid_flow"] * FT3_PER_M3 / area
    v_sg = values["gas_flow"] * FT3_PER_M3 / area
    v_mix = v_sl + v_sg

    rho_g = values["gas_density"] / KG_M3_PER_LB_FT3
    rho_l = values["liquid_density"] / KG_M3_PER_LB_FT3
    mu_g = values["gas_viscosity"] * CP_PER_PA_S
    mu_l = values["liquid_viscosity"] * CP_PER_PA_S

    if values["liquid_flow"] <= 0.0:
        regime, lam, el = "single_phase", 0.0, 0.0
        rho_mix, rho_ns, mu_ns, s_multiplier = rho_g, rho_g, mu_g, 0.0
    elif values["gas_flow"] <= 0.0:
        regime, lam, el = "single_phase", 1.0, 1.0
        rho_mix, rho_ns, mu_ns, s_multiplier = rho_l, rho_l, mu_l, 0.0
    else:
        lam = v_sl / v_mix
        froude = v_mix * v_mix / (FROUDE_G * d_ft)
        l1 = 316.0 * lam**0.302
        l2 = 0.000_925_2 * lam**-2.468_4
        l3 = 0.1 * lam**-1.451_6
        l4 = 0.5 * lam**-6.738
        if (lam < 0.01 and froude < l1) or (lam >= 0.01 and froude < l2):
            regime = "segregated"
        elif lam >= 0.01 and l2 <= froude <= l3:
            regime = "transition"
        elif (lam < 0.4 and l3 < froude <= l1) or (lam >= 0.4 and l3 < froude <= l4):
            regime = "intermittent"
        elif (lam < 0.4 and froude >= l1) or (lam >= 0.4 and froude > l4):
            regime = "distributed"
        else:
            from azoth.core.errors import InvalidInputError

            raise InvalidInputError(
                "gas_flow", "the Beggs-and-Brill regime map has no regime for this state"
            )
        a = (l3 - froude) / (l3 - l2)

        nvl = _liquid_velocity_number(rho_l, v_sl, values["temperature"], values["pressure"])
        theta = math.asin(values["elevation"] / values["length"])
        sin18 = math.sin(1.8 * theta)

        if regime == "transition":
            holdup = a * _inclination_corrected_holdup(
                "segregated", lam, froude, nvl, sin18, e_ft
            ) + (1.0 - a) * _inclination_corrected_holdup(
                "intermittent", lam, froude, nvl, sin18, e_ft
            )
        else:
            holdup = _inclination_corrected_holdup(regime, lam, froude, nvl, sin18, e_ft)
        el = max(lam, min(1.0, holdup))
        rho_mix = rho_l * el + rho_g * (1.0 - el)
        rho_ns = rho_l * lam + rho_g * (1.0 - lam)
        mu_ns = mu_l * lam + mu_g * (1.0 - lam)

        y = lam / (el * el)
        if 1.0 < y < 1.2:
            s_multiplier = math.log(2.2 * y - 1.2)
        else:
            ln_y = math.log(y)
            s_multiplier = ln_y / (-0.0523 + 3.182 * ln_y - 0.8725 * ln_y**2 + 0.01853 * ln_y**4)

    reynolds = rho_ns * v_mix * d_ft * (16.0 / (3.28 * 3.28)) / (0.001 * mu_ns)
    relative_roughness = values["roughness"] / values["diameter"]
    if abs(reynolds) < 1.0e-10:
        friction_factor = 0.0
    elif abs(reynolds) < 2300.0:
        friction_factor = 64.0 / reynolds
    elif abs(reynolds) < 4000.0:
        laminar = 64.0 / 2300.0
        turbulent = (
            1.0 / (-1.8 * math.log10((relative_roughness / 3.7) ** 1.11 + 6.9 / 4000.0))
        ) ** 2
        friction_factor = laminar + (turbulent - laminar) * (reynolds - 2300.0) / 1700.0
    else:
        friction_factor = (
            1.0 / (-1.8 * math.log10((relative_roughness / 3.7) ** 1.11 + 6.9 / reynolds))
        ) ** 2
    two_phase = friction_factor * math.exp(s_multiplier)

    hydrostatic_raw = rho_mix * HYDROSTATIC_G * e_ft
    friction_raw = two_phase * v_mix * v_mix * rho_ns * l_ft / (2.0 * d_ft)
    to_pa = BAR_PER_IMPERIAL_DROP * 1.0e5

    return BeggsBrillSegmentResult(
        regime=regime,
        liquid_holdup=el,
        superficial_gas_velocity=from_si(v_sg / FEET_PER_METER, "m/s"),
        superficial_liquid_velocity=from_si(v_sl / FEET_PER_METER, "m/s"),
        superficial_mixture_velocity=from_si(v_mix / FEET_PER_METER, "m/s"),
        no_slip_density=from_si(rho_ns * PROFILE_LB_FT3_TO_KG_M3, "kg/m**3"),
        no_slip_viscosity=from_si(mu_ns / CP_PER_PA_S, "Pa*s"),
        reynolds=reynolds,
        friction_factor=friction_factor,
        two_phase_friction_factor=two_phase,
        hydrostatic_drop=from_si(hydrostatic_raw * to_pa, "Pa"),
        friction_drop=from_si(friction_raw * to_pa, "Pa"),
        total_drop=from_si((hydrostatic_raw + friction_raw) * to_pa, "Pa"),
        warnings=tuple(warnings),
    )
