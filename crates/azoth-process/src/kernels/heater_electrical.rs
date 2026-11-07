//! `HeatExchangerElectricalDesign` - the duty as shaft power, the IEC motor tables, the cable and
//! the switchgear.
//!
//! **The electric heater's "motor" is the element bank.** `getProcessShaftPowerKW` returns the
//! process duty in kW, the base `ElectricalDesign` sizes a standard IEC motor on it, and the
//! motor's tabulated efficiency and power factor are what the electrical input, current and
//! losses come from. A shell-and-tube heater has no motor: its shaft power is zero and only the
//! auxiliary load is left, and it carries no cable. The air-cooler branch sizes fan power at
//! `max(2 x fans, 1% duty)`.
//!
//! `ElectricalCable.sizeCable` picks the smallest tray ampacity that carries the full-load
//! current through the 40 C ambient and tray-grouping derating, then upsizes for a 5 % drop.
//! `Switchgear.sizeSwitchgear` picks the starter from the rated power and the voltage's class, a
//! 1.25x breaker, a 1.6x fuse and the class's short-circuit rating and weight.
//!
//! Not here: the VFD, transformer and hazardous-area components (each its own class), the motor's
//! frame, weight and cost, the NEMA frame table, and the cable's and switchgear's costs.

/// `HeatExchangerInstrumentDesign.HeatExchangerType`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum HeatExchangerType {
    /// An electric resistance heater; the class detects this for a `Heater`.
    ElectricHeater,
    /// An air cooler's fan instrumentation.
    AirCooler,
    /// A shell-and-tube exchanger's utility-side instrumentation.
    ShellAndTube,
}

impl HeatExchangerType {
    /// The class's string spelling.
    #[must_use]
    pub fn parse(spelling: &str) -> Option<Self> {
        Some(match spelling {
            "electric_heater" => Self::ElectricHeater,
            "air_cooler" => Self::AirCooler,
            "shell_and_tube" => Self::ShellAndTube,
            _ => return None,
        })
    }
}

/// `ElectricalMotor.IEC_STANDARD_POWERS_KW`.
pub const IEC_POWERS_KW: [f64; 56] = [
    0.37, 0.55, 0.75, 1.1, 1.5, 2.2, 3.0, 4.0, 5.5, 7.5, 11.0, 15.0, 18.5, 22.0, 30.0, 37.0, 45.0,
    55.0, 75.0, 90.0, 110.0, 132.0, 160.0, 200.0, 250.0, 315.0, 355.0, 400.0, 450.0, 500.0, 560.0,
    630.0, 710.0, 800.0, 900.0, 1000.0, 1120.0, 1250.0, 1400.0, 1600.0, 1800.0, 2000.0, 2240.0,
    2500.0, 2800.0, 3150.0, 3550.0, 4000.0, 4500.0, 5000.0, 5600.0, 6300.0, 7100.0, 8000.0, 9000.0,
    10000.0,
];

/// `ElectricalMotor.efficiencyClass`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EfficiencyClass {
    /// IE1.
    Ie1,
    /// IE2.
    Ie2,
    /// IE3, the class's default.
    Ie3,
    /// IE4.
    Ie4,
}

impl EfficiencyClass {
    /// The class's string spelling.
    #[must_use]
    pub fn parse(spelling: &str) -> Option<Self> {
        Some(match spelling {
            "ie1" => Self::Ie1,
            "ie2" => Self::Ie2,
            "ie3" => Self::Ie3,
            "ie4" => Self::Ie4,
            _ => return None,
        })
    }

    fn offset(self) -> f64 {
        match self {
            Self::Ie1 => -3.0,
            Self::Ie2 => -1.5,
            Self::Ie3 => 0.0,
            Self::Ie4 => 1.0,
        }
    }
}

/// The electrical design's settings, with the class's defaults resolved by the model.
#[derive(Debug, Clone, Copy)]
pub struct ElectricalRequest {
    /// The motor sizing margin over the shaft power.
    pub motor_sizing_margin: f64,
    /// The motor's pole count.
    pub motor_poles: u32,
    /// The efficiency class.
    pub efficiency_class: EfficiencyClass,
    /// The supply frequency, Hz.
    pub frequency_hz: f64,
    /// The supply phases.
    pub phases: u32,
    /// The air cooler's fan count.
    pub number_of_fans: u32,
    /// The auxiliary instrumentation load, W.
    pub instrumentation_power_w: f64,
    /// The shell-and-tube cooling-water pump load, W.
    pub cooling_water_pump_power_w: f64,
    /// The power cable's route length, m.
    pub cable_length_m: f64,
}

/// What the electrical design computed, in SI.
#[derive(Debug, Clone, Copy)]
pub struct ElectricalDesign {
    /// The shaft power the type implies, W.
    pub shaft_power_w: f64,
    /// The electrical input, W.
    pub electrical_input_w: f64,
    /// The apparent power, VA.
    pub apparent_power_va: f64,
    /// The reactive power, VAR.
    pub reactive_power_var: f64,
    /// The power factor.
    pub power_factor: f64,
    /// The selected voltage, V.
    pub rated_voltage_v: f64,
    /// The stated supply frequency, Hz.
    pub frequency_hz: f64,
    /// The full-load current, A.
    pub full_load_current_a: f64,
    /// The starting current, A.
    pub starting_current_a: f64,
    /// The electrical losses, W.
    pub total_electrical_losses_w: f64,
    /// The selected motor's rated power, W.
    pub motor_rated_power_w: f64,
    /// The motor's efficiency at full load, a fraction.
    pub motor_efficiency: f64,
    /// The motor's power factor at full load.
    pub motor_power_factor: f64,
    /// The motor's rated speed, rpm.
    pub motor_rated_speed_rpm: f64,
    /// The motor's locked-rotor current multiplier.
    pub motor_locked_rotor_multiplier: f64,
    /// The auxiliary load, W.
    pub total_auxiliary_w: f64,
    /// The electrical input plus the auxiliary load, W.
    pub total_connected_w: f64,
    /// The power cable's route length, m.
    pub cable_length_m: f64,
    /// The selected cable cross-section, mm2.
    pub cable_cross_section_mm2: f64,
    /// The derated cable ampacity, A.
    pub cable_ampacity_a: f64,
    /// The cable's voltage drop, percent.
    pub cable_voltage_drop_percent: f64,
    /// The cable's 1-second short-circuit withstand, kA.
    pub cable_short_circuit_withstand_ka: f64,
    /// The switchgear's starter type.
    pub switchgear_starter: &'static str,
    /// The switchgear's type.
    pub switchgear_type: &'static str,
    /// The circuit-breaker rating, A.
    pub switchgear_circuit_breaker_a: f64,
    /// The fuse rating, A.
    pub switchgear_fuse_a: f64,
    /// The switchgear's rated current, A.
    pub switchgear_rated_current_a: f64,
    /// The switchgear's short-circuit rating, kA.
    pub switchgear_short_circuit_ka: f64,
    /// The switchgear's estimated weight, kg.
    pub switchgear_weight_kg: f64,
}

fn select_power(required_kw: f64) -> f64 {
    // The comparison carries a 1e-9 relative slack: the duty the flash reports at a standard
    // size lands a few ulps either side of it, and the class's own run on the captured state
    // reports just under 100 kW where this run reports just over. Without the slack the table
    // jumps a frame on a rounding difference that is not a sizing difference.
    for standard in IEC_POWERS_KW {
        if standard >= required_kw * (1.0 - 1.0e-9) {
            return standard;
        }
    }
    (required_kw / 100.0).ceil() * 100.0
}

fn slip(power_kw: f64) -> f64 {
    if power_kw <= 7.5 {
        0.05
    } else if power_kw <= 75.0 {
        0.03
    } else if power_kw <= 375.0 {
        0.02
    } else {
        0.015
    }
}

fn efficiency(power_kw: f64, poles: u32, class: EfficiencyClass) -> f64 {
    let base = if power_kw <= 1.1 {
        82.5
    } else if power_kw <= 4.0 {
        87.7
    } else if power_kw <= 11.0 {
        91.0
    } else if power_kw <= 37.0 {
        93.3
    } else if power_kw <= 110.0 {
        95.0
    } else if power_kw <= 375.0 {
        96.0
    } else {
        96.5
    };
    let pole_offset = match poles {
        2 => -0.5,
        6 => -1.0,
        8 => -1.5,
        _ => 0.0,
    };
    (base + class.offset() + pole_offset).min(99.0)
}

fn power_factor(power_kw: f64) -> f64 {
    if power_kw <= 4.0 {
        0.78
    } else if power_kw <= 15.0 {
        0.83
    } else if power_kw <= 55.0 {
        0.86
    } else if power_kw <= 200.0 {
        0.88
    } else if power_kw <= 1000.0 {
        0.90
    } else {
        0.91
    }
}

fn locked_rotor(power_kw: f64) -> f64 {
    if power_kw <= 7.5 {
        6.0
    } else if power_kw <= 75.0 {
        6.5
    } else if power_kw <= 500.0 {
        7.0
    } else {
        6.5
    }
}

fn current(input_kw: f64, voltage_v: f64, factor: f64, phases: u32) -> f64 {
    if voltage_v <= 0.0 || factor <= 0.0 {
        return 0.0;
    }
    if phases == 3 {
        input_kw * 1000.0 / (3.0_f64.sqrt() * voltage_v * factor)
    } else {
        input_kw * 1000.0 / (voltage_v * factor)
    }
}

/// `ElectricalCable.STANDARD_CROSS_SECTIONS`, mm2.
pub const CABLE_CROSS_SECTIONS_MM2: [f64; 19] = [
    1.5, 2.5, 4.0, 6.0, 10.0, 16.0, 25.0, 35.0, 50.0, 70.0, 95.0, 120.0, 150.0, 185.0, 240.0,
    300.0, 400.0, 500.0, 630.0,
];
/// `BASE_AMPACITY_XLPE_TRAY`, A.
pub const CABLE_BASE_AMPACITY_A: [f64; 19] = [
    18.0, 25.0, 33.0, 43.0, 60.0, 80.0, 106.0, 131.0, 158.0, 200.0, 241.0, 278.0, 315.0, 360.0,
    420.0, 480.0, 560.0, 640.0, 740.0,
];
/// `Switchgear.STANDARD_RATINGS_A`.
pub const SWITCHGEAR_RATINGS_A: [f64; 14] = [
    100.0, 160.0, 250.0, 400.0, 630.0, 800.0, 1000.0, 1250.0, 1600.0, 2000.0, 2500.0, 3150.0,
    4000.0, 5000.0,
];

fn next_cable_size(min_mm2: f64) -> f64 {
    for standard in CABLE_CROSS_SECTIONS_MM2 {
        if standard >= min_mm2 {
            return standard;
        }
    }
    *CABLE_CROSS_SECTIONS_MM2.last().unwrap_or(&630.0)
}

fn voltage_drop_percent(current_a: f64, voltage_v: f64, length_m: f64, section: f64) -> f64 {
    if section <= 0.0 || voltage_v <= 0.0 || length_m <= 0.0 {
        return 0.0;
    }
    let resistance = 0.0175 / section;
    let reactance = 0.000_08;
    let cos_phi = 0.85_f64;
    let sin_phi = (1.0 - cos_phi * cos_phi).sqrt();
    let drop = 3.0_f64.sqrt() * current_a * length_m * (resistance * cos_phi + reactance * sin_phi);
    drop / voltage_v * 100.0
}

/// `ElectricalCable.sizeCable` with the heater path's Tray/40 C/copper defaults.
fn size_cable(current_a: f64, voltage_v: f64, length_m: f64) -> (f64, f64, f64, f64) {
    let ambient = ((90.0_f64 - 40.0) / (90.0 - 30.0)).sqrt().clamp(0.5, 1.0);
    let grouping = 0.85;
    let derating = ambient * grouping;
    let required = if derating > 0.0 {
        current_a / derating
    } else {
        current_a
    };
    let mut section = CABLE_CROSS_SECTIONS_MM2[0];
    let mut base_ampacity = CABLE_BASE_AMPACITY_A[0];
    for (index, standard) in CABLE_CROSS_SECTIONS_MM2.iter().enumerate() {
        if CABLE_BASE_AMPACITY_A[index] >= required {
            section = *standard;
            base_ampacity = CABLE_BASE_AMPACITY_A[index];
            break;
        }
        if index == CABLE_CROSS_SECTIONS_MM2.len() - 1 {
            section = *standard;
            base_ampacity = CABLE_BASE_AMPACITY_A[index];
        }
    }
    let ampacity = base_ampacity * derating;
    let mut drop = voltage_drop_percent(current_a, voltage_v, length_m, section);
    while drop > 5.0 && section < CABLE_CROSS_SECTIONS_MM2[CABLE_CROSS_SECTIONS_MM2.len() - 1] {
        section = next_cable_size(section + 1.0);
        drop = voltage_drop_percent(current_a, voltage_v, length_m, section);
    }
    let withstand = 143.0 * section / 1000.0;
    (section, ampacity, drop, withstand)
}

/// `Switchgear.sizeSwitchgear` with `useVFD = false`.
fn size_switchgear(
    current_a: f64,
    rated_power_kw: f64,
    voltage_v: f64,
) -> (&'static str, &'static str, f64, f64, f64, f64, f64) {
    let starter = if rated_power_kw > 200.0 {
        "Soft Starter"
    } else if rated_power_kw > 11.0 {
        "Star-Delta"
    } else {
        "DOL"
    };
    let switchgear_type = if voltage_v > 1000.0 {
        "MV Switchgear"
    } else {
        "MCC"
    };
    let breaker = {
        let required = current_a * 1.25;
        SWITCHGEAR_RATINGS_A
            .iter()
            .copied()
            .find(|rating| *rating >= required)
            .unwrap_or((required / 100.0).ceil() * 100.0)
    };
    let fuse = {
        let target = current_a * 1.6;
        const FUSES: [f64; 20] = [
            16.0, 20.0, 25.0, 32.0, 40.0, 50.0, 63.0, 80.0, 100.0, 125.0, 160.0, 200.0, 250.0,
            315.0, 400.0, 500.0, 630.0, 800.0, 1000.0, 1250.0,
        ];
        FUSES
            .iter()
            .copied()
            .find(|rating| *rating >= target)
            .unwrap_or(FUSES[FUSES.len() - 1])
    };
    let short_circuit = if voltage_v > 6000.0 {
        40.0
    } else if voltage_v > 1000.0 {
        31.5
    } else {
        25.0
    };
    let weight = if voltage_v > 1000.0 {
        breaker * 0.5 + 200.0
    } else {
        breaker * 0.15 + 50.0
    };
    (
        starter,
        switchgear_type,
        breaker,
        fuse,
        breaker,
        short_circuit,
        weight,
    )
}

/// `HeatExchangerElectricalDesign.calcDesign`.
#[must_use]
pub fn design(
    kind: HeatExchangerType,
    process_duty_w: f64,
    request: &ElectricalRequest,
) -> ElectricalDesign {
    let duty_kw = process_duty_w.abs() / 1000.0;
    let shaft_kw = match kind {
        HeatExchangerType::ElectricHeater => duty_kw,
        HeatExchangerType::AirCooler => {
            (2.0 * f64::from(request.number_of_fans)).max(duty_kw * 0.01)
        }
        HeatExchangerType::ShellAndTube => 0.0,
    };
    let rated_voltage = match kind {
        HeatExchangerType::ElectricHeater => {
            if shaft_kw > 200.0 {
                6600.0
            } else if shaft_kw > 75.0 {
                690.0
            } else {
                400.0
            }
        }
        HeatExchangerType::AirCooler => {
            if shaft_kw > 200.0 {
                690.0
            } else {
                400.0
            }
        }
        HeatExchangerType::ShellAndTube => 400.0,
    };

    let (rated_power, motor_efficiency, motor_power_factor, speed, locked) = if shaft_kw > 0.0 {
        let required = shaft_kw * request.motor_sizing_margin;
        let rated = select_power(required);
        let poles = request.motor_poles.max(1);
        // The motor is sized before the design's frequency is copied to it, so the
        // synchronous speed always uses the motor's own 50 Hz default.
        let synchronous = 120.0 * 50.0 / f64::from(poles);
        (
            rated,
            efficiency(rated, poles, request.efficiency_class),
            power_factor(rated),
            synchronous * (1.0 - slip(rated)),
            locked_rotor(rated),
        )
    } else {
        (0.0, 95.0, 0.86, 0.0, 7.0)
    };
    let instrumentation = request.instrumentation_power_w;
    let auxiliary = instrumentation
        + if kind == HeatExchangerType::ShellAndTube {
            request.cooling_water_pump_power_w
        } else {
            0.0
        };

    if kind == HeatExchangerType::ShellAndTube {
        let input_w = auxiliary;
        let factor = 0.85;
        let apparent = if factor > 0.0 { input_w / factor } else { 0.0 };
        let input_kw = input_w / 1000.0;
        let full_load = current(input_kw, rated_voltage, factor, request.phases);
        return ElectricalDesign {
            shaft_power_w: 0.0,
            electrical_input_w: input_w,
            apparent_power_va: apparent,
            reactive_power_var: apparent * (1.0 - factor * factor).sqrt(),
            power_factor: factor,
            rated_voltage_v: rated_voltage,
            frequency_hz: request.frequency_hz,
            full_load_current_a: full_load,
            starting_current_a: full_load * locked,
            total_electrical_losses_w: input_w,
            motor_rated_power_w: 0.0,
            motor_efficiency: motor_efficiency / 100.0,
            motor_power_factor,
            motor_rated_speed_rpm: speed,
            motor_locked_rotor_multiplier: locked,
            total_auxiliary_w: auxiliary,
            total_connected_w: input_w + auxiliary,
            cable_length_m: request.cable_length_m,
            cable_cross_section_mm2: 0.0,
            cable_ampacity_a: 0.0,
            cable_voltage_drop_percent: 0.0,
            cable_short_circuit_withstand_ka: 0.0,
            switchgear_starter: "DOL",
            switchgear_type: "MCC",
            switchgear_circuit_breaker_a: 0.0,
            switchgear_fuse_a: 0.0,
            switchgear_rated_current_a: 0.0,
            switchgear_short_circuit_ka: 25.0,
            switchgear_weight_kg: 0.0,
        };
    }

    let efficiency_fraction = if motor_efficiency > 0.0 {
        motor_efficiency / 100.0
    } else {
        0.0
    };
    let input_kw = if efficiency_fraction > 0.0 {
        shaft_kw / efficiency_fraction
    } else {
        shaft_kw
    };
    let apparent = if motor_power_factor > 0.0 {
        input_kw / motor_power_factor
    } else {
        0.0
    };
    let full_load = current(input_kw, rated_voltage, motor_power_factor, request.phases);
    let (cable_section, cable_ampacity, cable_drop, cable_withstand) =
        size_cable(full_load, rated_voltage, request.cable_length_m);
    let (starter, switchgear_type, breaker, fuse, rated_current, short_circuit, weight) =
        size_switchgear(full_load, rated_power, rated_voltage);
    ElectricalDesign {
        shaft_power_w: shaft_kw * 1000.0,
        electrical_input_w: input_kw * 1000.0,
        apparent_power_va: apparent * 1000.0,
        reactive_power_var: apparent
            * 1000.0
            * (1.0 - motor_power_factor * motor_power_factor).sqrt(),
        power_factor: motor_power_factor,
        rated_voltage_v: rated_voltage,
        frequency_hz: request.frequency_hz,
        full_load_current_a: full_load,
        starting_current_a: full_load * locked,
        total_electrical_losses_w: (input_kw - shaft_kw) * 1000.0,
        motor_rated_power_w: rated_power * 1000.0,
        motor_efficiency: efficiency_fraction,
        motor_power_factor,
        motor_rated_speed_rpm: speed,
        motor_locked_rotor_multiplier: locked,
        total_auxiliary_w: auxiliary,
        total_connected_w: (input_kw * 1000.0) + auxiliary,
        cable_length_m: request.cable_length_m,
        cable_cross_section_mm2: cable_section,
        cable_ampacity_a: cable_ampacity,
        cable_voltage_drop_percent: cable_drop,
        cable_short_circuit_withstand_ka: cable_withstand,
        switchgear_starter: starter,
        switchgear_type,
        switchgear_circuit_breaker_a: breaker,
        switchgear_fuse_a: fuse,
        switchgear_rated_current_a: rated_current,
        switchgear_short_circuit_ka: short_circuit,
        switchgear_weight_kg: weight,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn request() -> ElectricalRequest {
        ElectricalRequest {
            motor_sizing_margin: 1.10,
            motor_poles: 4,
            efficiency_class: EfficiencyClass::Ie3,
            frequency_hz: 50.0,
            phases: 3,
            number_of_fans: 2,
            instrumentation_power_w: 1500.0,
            cooling_water_pump_power_w: 0.0,
            cable_length_m: 50.0,
        }
    }

    fn close(actual: f64, expected: f64, tolerance: f64, what: &str) {
        let error = (actual - expected).abs() / expected.abs().max(f64::MIN_POSITIVE);
        assert!(
            error <= tolerance,
            "{what}: {actual} against {expected}, relative error {error:.3e} > {tolerance:.3e}"
        );
    }

    #[test]
    fn the_five_kilowatt_heater_is_the_capture() {
        let design = design(
            HeatExchangerType::ElectricHeater,
            4999.9999986086605,
            &request(),
        );
        close(design.electrical_input_w, 5494.50549297655, 1.0e-9, "input");
        close(
            design.apparent_power_va,
            6619.886136116326,
            1.0e-9,
            "apparent",
        );
        close(
            design.reactive_power_var,
            3692.330137297538,
            1.0e-9,
            "reactive",
        );
        close(
            design.full_load_current_a,
            9.55498260672858,
            1.0e-12,
            "current",
        );
        close(
            design.starting_current_a,
            57.329895640371475,
            1.0e-12,
            "starting",
        );
        close(design.motor_rated_power_w, 5500.0, 1.0e-12, "motor");
        close(design.motor_efficiency, 0.91, 1.0e-12, "efficiency");
        close(design.motor_rated_speed_rpm, 1425.0, 1.0e-12, "speed");
        close(
            design.total_connected_w,
            6994.50549297655,
            1.0e-9,
            "connected",
        );
    }

    #[test]
    fn the_hundred_kilowatt_heater_takes_the_690_volt_tap() {
        let design = design(
            HeatExchangerType::ElectricHeater,
            99_999.99997673612,
            &request(),
        );
        close(design.rated_voltage_v, 690.0, 1.0e-12, "voltage");
        close(design.motor_rated_power_w, 110_000.0, 1.0e-12, "motor");
        close(design.motor_efficiency, 0.95, 1.0e-12, "efficiency");
        close(design.motor_power_factor, 0.88, 1.0e-12, "factor");
        close(design.motor_rated_speed_rpm, 1470.0, 1.0e-12, "speed");
        close(design.motor_locked_rotor_multiplier, 7.0, 1.0e-12, "locked");
        close(
            design.full_load_current_a,
            100.08845937440393,
            1.0e-9,
            "current",
        );
        close(
            design.starting_current_a,
            700.6192156208275,
            1.0e-9,
            "starting",
        );
        close(
            design.total_connected_w,
            106_763.15787024854,
            1.0e-9,
            "connected",
        );
    }

    #[test]
    fn the_five_hundred_kilowatt_heater_takes_the_6600_volt_tap() {
        let design = design(
            HeatExchangerType::ElectricHeater,
            499_999.9998608659,
            &request(),
        );
        close(design.rated_voltage_v, 6600.0, 1.0e-12, "voltage");
        close(design.motor_rated_power_w, 560_000.0, 1.0e-12, "motor");
        close(design.motor_efficiency, 0.965, 1.0e-12, "efficiency");
        close(design.motor_rated_speed_rpm, 1477.5, 1.0e-12, "speed");
        close(design.motor_locked_rotor_multiplier, 6.5, 1.0e-12, "locked");
        close(
            design.total_electrical_losses_w,
            18_134.715020860483,
            1.0e-9,
            "losses",
        );
    }

    #[test]
    fn the_shell_and_tube_design_is_the_auxiliary_load() {
        let mut request = request();
        request.cooling_water_pump_power_w = 3000.0;
        let design = design(HeatExchangerType::ShellAndTube, 100_000.0, &request);
        close(design.shaft_power_w, 0.0, 1.0e-12, "shaft");
        close(design.electrical_input_w, 4500.0, 1.0e-12, "input");
        close(design.power_factor, 0.85, 1.0e-12, "factor");
        close(
            design.apparent_power_va,
            5294.117647058823,
            1.0e-12,
            "apparent",
        );
        close(
            design.full_load_current_a,
            7.6414006216274,
            1.0e-12,
            "current",
        );
        close(design.motor_rated_power_w, 0.0, 1.0e-12, "motor");
        close(design.total_auxiliary_w, 4500.0, 1.0e-12, "auxiliary");
        close(design.total_connected_w, 9000.0, 1.0e-12, "connected");
    }

    #[test]
    fn the_air_cooler_sizes_the_fans() {
        let design = design(HeatExchangerType::AirCooler, 100_000.0, &request());
        close(design.shaft_power_w, 4000.0, 1.0e-12, "shaft");
        close(design.motor_rated_power_w, 5500.0, 1.0e-12, "motor");
        close(
            design.electrical_input_w,
            4395.604395604395,
            1.0e-9,
            "input",
        );
        close(
            design.starting_current_a,
            45.863916525059636,
            1.0e-9,
            "starting",
        );
    }

    #[test]
    fn the_cable_and_switchgear_are_the_capture() {
        let five = design(
            HeatExchangerType::ElectricHeater,
            4_999.9999986086605,
            &request(),
        );
        close(five.cable_cross_section_mm2, 1.5, 1.0e-12, "5 kW section");
        close(
            five.cable_ampacity_a,
            13.966925216381735,
            1.0e-12,
            "5 kW ampacity",
        );
        close(
            five.cable_voltage_drop_percent,
            2.0601932342502276,
            1.0e-9,
            "5 kW drop",
        );
        close(
            five.cable_short_circuit_withstand_ka,
            0.2145,
            1.0e-12,
            "5 kW withstand",
        );
        assert_eq!(five.switchgear_starter, "DOL");
        assert_eq!(five.switchgear_type, "MCC");
        close(five.switchgear_circuit_breaker_a, 100.0, 1.0e-12, "breaker");
        close(five.switchgear_fuse_a, 16.0, 1.0e-12, "fuse");
        close(five.switchgear_rated_current_a, 100.0, 1.0e-12, "rated");
        close(five.switchgear_short_circuit_ka, 25.0, 1.0e-12, "sc");
        close(five.switchgear_weight_kg, 65.0, 1.0e-12, "weight");

        let hundred = design(
            HeatExchangerType::ElectricHeater,
            99_999.99997673612,
            &request(),
        );
        close(
            hundred.cable_cross_section_mm2,
            35.0,
            1.0e-12,
            "100 kW section",
        );
        close(
            hundred.cable_ampacity_a,
            101.64817796366708,
            1.0e-12,
            "100 kW ampacity",
        );
        close(
            hundred.cable_voltage_drop_percent,
            0.5868336821889366,
            1.0e-9,
            "100 kW drop",
        );
        close(
            hundred.cable_short_circuit_withstand_ka,
            5.005,
            1.0e-12,
            "100 kW withstand",
        );
        assert_eq!(hundred.switchgear_starter, "Star-Delta");
        close(
            hundred.switchgear_circuit_breaker_a,
            160.0,
            1.0e-12,
            "breaker",
        );
        close(hundred.switchgear_fuse_a, 200.0, 1.0e-12, "fuse");
        close(hundred.switchgear_weight_kg, 74.0, 1.0e-12, "weight");

        let five_hundred = design(
            HeatExchangerType::ElectricHeater,
            499_999.9998608659,
            &request(),
        );
        close(
            five_hundred.cable_cross_section_mm2,
            16.0,
            1.0e-12,
            "500 kW section",
        );
        close(
            five_hundred.cable_ampacity_a,
            62.07522318391882,
            1.0e-12,
            "500 kW ampacity",
        );
        close(
            five_hundred.cable_voltage_drop_percent,
            0.06422034990752382,
            1.0e-9,
            "500 kW drop",
        );
        close(
            five_hundred.cable_short_circuit_withstand_ka,
            2.288,
            1.0e-12,
            "500 kW withstand",
        );
        assert_eq!(five_hundred.switchgear_starter, "Soft Starter");
        assert_eq!(five_hundred.switchgear_type, "MV Switchgear");
        close(
            five_hundred.switchgear_short_circuit_ka,
            40.0,
            1.0e-12,
            "500 kW sc",
        );
        close(
            five_hundred.switchgear_weight_kg,
            250.0,
            1.0e-12,
            "500 kW weight",
        );
    }

    #[test]
    fn the_shell_branch_carries_no_cable() {
        let mut request = request();
        request.cooling_water_pump_power_w = 3000.0;
        let shell = design(HeatExchangerType::ShellAndTube, 100_000.0, &request);
        close(shell.cable_length_m, 50.0, 1.0e-12, "length");
        close(shell.cable_cross_section_mm2, 0.0, 1.0e-12, "section");
        close(shell.cable_ampacity_a, 0.0, 1.0e-12, "ampacity");
        close(shell.cable_voltage_drop_percent, 0.0, 1.0e-12, "drop");
        close(
            shell.cable_short_circuit_withstand_ka,
            0.0,
            1.0e-12,
            "withstand",
        );
        assert_eq!(shell.switchgear_starter, "DOL");
        close(shell.switchgear_circuit_breaker_a, 0.0, 1.0e-12, "breaker");
        close(shell.switchgear_short_circuit_ka, 25.0, 1.0e-12, "sc");
        close(shell.switchgear_weight_kg, 0.0, 1.0e-12, "weight");
    }

    #[test]
    fn the_cable_upsizes_for_voltage_drop() {
        let mut request = request();
        request.cable_length_m = 500.0;
        let half_kilometre = design(
            HeatExchangerType::ElectricHeater,
            99_999.99997673612,
            &request,
        );
        close(
            half_kilometre.cable_cross_section_mm2,
            50.0,
            1.0e-12,
            "500 m section",
        );
        close(
            half_kilometre.cable_ampacity_a,
            101.64817796366708,
            1.0e-12,
            "500 m ampacity",
        );
        close(
            half_kilometre.cable_voltage_drop_percent,
            4.2666569049915735,
            1.0e-9,
            "500 m drop",
        );
        close(
            half_kilometre.cable_short_circuit_withstand_ka,
            7.15,
            1.0e-12,
            "500 m withstand",
        );

        request.cable_length_m = 2000.0;
        let two_kilometres = design(
            HeatExchangerType::ElectricHeater,
            99_999.99997673612,
            &request,
        );
        close(
            two_kilometres.cable_cross_section_mm2,
            300.0,
            1.0e-12,
            "2000 m section",
        );
        close(
            two_kilometres.cable_ampacity_a,
            101.64817796366708,
            1.0e-12,
            "2000 m ampacity",
        );
        close(
            two_kilometres.cable_voltage_drop_percent,
            4.609117155205686,
            1.0e-9,
            "2000 m drop",
        );
        close(
            two_kilometres.cable_short_circuit_withstand_ka,
            42.9,
            1.0e-12,
            "2000 m withstand",
        );
    }

    #[test]
    fn the_air_cooler_cable_is_its_own_current() {
        let air = design(HeatExchangerType::AirCooler, 100_000.0, &request());
        close(air.cable_cross_section_mm2, 1.5, 1.0e-12, "section");
        close(
            air.cable_ampacity_a,
            13.966925216381735,
            1.0e-12,
            "ampacity",
        );
        close(
            air.cable_voltage_drop_percent,
            1.6481545878588109,
            1.0e-9,
            "drop",
        );
    }

    #[test]
    fn the_efficiency_class_and_pole_count_move_the_table() {
        let mut ie4 = request();
        ie4.efficiency_class = EfficiencyClass::Ie4;
        close(
            design(HeatExchangerType::ElectricHeater, 99_999.99997673612, &ie4).motor_efficiency,
            0.96,
            1.0e-12,
            "IE4 efficiency",
        );
        let mut two_pole = request();
        two_pole.motor_poles = 2;
        let design = design(
            HeatExchangerType::ElectricHeater,
            99_999.99997673612,
            &two_pole,
        );
        close(design.motor_efficiency, 0.945, 1.0e-12, "2-pole efficiency");
        close(
            design.motor_rated_speed_rpm,
            2940.0,
            1.0e-12,
            "2-pole speed",
        );
    }
}
