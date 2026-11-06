//! `HeatExchangerElectricalDesign` - the duty as shaft power and the IEC motor tables.
//!
//! **The electric heater's "motor" is the element bank.** `getProcessShaftPowerKW` returns the
//! process duty in kW, the base `ElectricalDesign` sizes a standard IEC motor on it, and the
//! motor's tabulated efficiency and power factor are what the electrical input, current and
//! losses come from. A shell-and-tube heater has no motor: its shaft power is zero and only the
//! auxiliary load is left. The air-cooler branch sizes fan power at `max(2 x fans, 1% duty)`.
//!
//! Not here: the cable, switchgear, VFD, transformer and hazardous-area components (each its own
//! class), the motor's frame, weight and cost, and the NEMA frame table - all named as outside.


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
    0.37, 0.55, 0.75, 1.1, 1.5, 2.2, 3.0, 4.0, 5.5, 7.5, 11.0, 15.0, 18.5, 22.0, 30.0, 37.0,
    45.0, 55.0, 75.0, 90.0, 110.0, 132.0, 160.0, 200.0, 250.0, 315.0, 355.0, 400.0, 450.0,
    500.0, 560.0, 630.0, 710.0, 800.0, 900.0, 1000.0, 1120.0, 1250.0, 1400.0, 1600.0, 1800.0,
    2000.0, 2240.0, 2500.0, 2800.0, 3150.0, 3550.0, 4000.0, 4500.0, 5000.0, 5600.0, 6300.0,
    7100.0, 8000.0, 9000.0, 10000.0,
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
    ElectricalDesign {
        shaft_power_w: shaft_kw * 1000.0,
        electrical_input_w: input_kw * 1000.0,
        apparent_power_va: apparent * 1000.0,
        reactive_power_var: apparent * 1000.0 * (1.0 - motor_power_factor * motor_power_factor).sqrt(),
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
        let design = design(HeatExchangerType::ElectricHeater, 4999.9999986086605, &request());
        close(design.electrical_input_w, 5494.50549297655, 1.0e-9, "input");
        close(design.apparent_power_va, 6619.886136116326, 1.0e-9, "apparent");
        close(design.reactive_power_var, 3692.330137297538, 1.0e-9, "reactive");
        close(design.full_load_current_a, 9.55498260672858, 1.0e-12, "current");
        close(
            design.starting_current_a,
            57.329895640371475,
            1.0e-12,
            "starting",
        );
        close(design.motor_rated_power_w, 5500.0, 1.0e-12, "motor");
        close(design.motor_efficiency, 0.91, 1.0e-12, "efficiency");
        close(design.motor_rated_speed_rpm, 1425.0, 1.0e-12, "speed");
        close(design.total_connected_w, 6994.50549297655, 1.0e-9, "connected");
    }

    #[test]
    fn the_hundred_kilowatt_heater_takes_the_690_volt_tap() {
        let design = design(HeatExchangerType::ElectricHeater, 99_999.99997673612, &request());
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
        close(design.total_connected_w, 106_763.15787024854, 1.0e-9, "connected");
    }

    #[test]
    fn the_five_hundred_kilowatt_heater_takes_the_6600_volt_tap() {
        let design = design(HeatExchangerType::ElectricHeater, 499_999.9998608659, &request());
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
        close(design.apparent_power_va, 5294.117647058823, 1.0e-12, "apparent");
        close(design.full_load_current_a, 7.6414006216274, 1.0e-12, "current");
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
        let design = design(HeatExchangerType::ElectricHeater, 99_999.99997673612, &two_pole);
        close(design.motor_efficiency, 0.945, 1.0e-12, "2-pole efficiency");
        close(design.motor_rated_speed_rpm, 2940.0, 1.0e-12, "2-pole speed");
    }
}
