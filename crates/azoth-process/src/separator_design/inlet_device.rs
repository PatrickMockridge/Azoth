//! `InletDeviceModel`: the inlet device's bulk removal, pressure drop and downstream DSD.
//!
//! The seven device types carry their own efficiency, DSD multiplier and pressure-drop
//! coefficient. The bulk efficiency degrades above a device-specific momentum limit (Bothamley
//! 2013) and is halved below a third of the optimum. The downstream distribution's
//! characteristic diameter is multiplied by the device's factor; a cyclone or vane narrows
//! the spread, and a half-pipe or deflector above `3000 Pa` broadens it through
//! re-atomisation.
//!
//! The source is `process/equipment/separator/entrainment/InletDeviceModel.java`.

/// The class's seven inlet devices, with `(bulk efficiency, DSD multiplier, dp coefficient)`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum InletDevice {
    None,
    DeflectorPlate,
    HalfPipe,
    InletVane,
    InletCyclone,
    Schoepentoeter,
    ImpingementPlate,
}

impl InletDevice {
    /// The class's enum spelling.
    #[must_use]
    pub fn name(self) -> &'static str {
        match self {
            Self::None => "NONE",
            Self::DeflectorPlate => "DEFLECTOR_PLATE",
            Self::HalfPipe => "HALF_PIPE",
            Self::InletVane => "INLET_VANE",
            Self::InletCyclone => "INLET_CYCLONE",
            Self::Schoepentoeter => "SCHOEPENTOETER",
            Self::ImpingementPlate => "IMPINGEMENT_PLATE",
        }
    }

    fn bulk_efficiency(self) -> f64 {
        match self {
            Self::None => 0.0,
            Self::DeflectorPlate => 0.60,
            Self::HalfPipe => 0.70,
            Self::InletVane => 0.80,
            Self::InletCyclone => 0.95,
            Self::Schoepentoeter => 0.85,
            Self::ImpingementPlate => 0.65,
        }
    }

    fn dsd_multiplier(self) -> f64 {
        match self {
            Self::None => 1.0,
            Self::DeflectorPlate => 0.85,
            Self::HalfPipe => 0.75,
            Self::InletVane => 0.55,
            Self::InletCyclone => 0.35,
            Self::Schoepentoeter => 0.50,
            Self::ImpingementPlate => 0.80,
        }
    }

    fn pressure_drop_coefficient(self) -> f64 {
        match self {
            Self::None => 0.0,
            Self::DeflectorPlate => 0.5,
            Self::HalfPipe => 1.0,
            Self::InletVane => 1.5,
            Self::InletCyclone => 3.0,
            Self::Schoepentoeter => 2.5,
            Self::ImpingementPlate => 0.8,
        }
    }

    fn max_momentum(self) -> Option<f64> {
        match self {
            Self::None => None,
            Self::DeflectorPlate => Some(1500.0),
            Self::HalfPipe => Some(3500.0),
            Self::InletVane => Some(6000.0),
            Self::InletCyclone => Some(15000.0),
            Self::Schoepentoeter => Some(9000.0),
            Self::ImpingementPlate => Some(2000.0),
        }
    }
}

/// The device's result, in SI.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct InletResult {
    pub bulk_efficiency: f64,
    pub pressure_drop: f64,
    pub nozzle_velocity: f64,
    pub momentum_flux: f64,
    /// The downstream distribution's `(d63.2, spread, is_log_normal)`.
    pub downstream: (f64, f64, bool),
}

/// `calculate`: the device on the incoming distribution.
#[allow(clippy::too_many_arguments)] // one per declared setting, as the class
#[must_use]
pub fn calculate(
    device: InletDevice,
    incoming: (f64, f64, bool),
    gas_density: f64,
    liquid_density: f64,
    gas_volume_flow: f64,
    liquid_volume_flow: f64,
    nozzle_diameter: f64,
    override_bulk_efficiency: f64,
    override_dsd_multiplier: f64,
) -> InletResult {
    let nozzle_area = std::f64::consts::PI * nozzle_diameter * nozzle_diameter / 4.0;
    let total_flow = gas_volume_flow + liquid_volume_flow;
    let nozzle_velocity = if nozzle_area > 0.0 {
        total_flow / nozzle_area
    } else {
        0.0
    };
    let gas_vol_fraction = if total_flow > 0.0 {
        gas_volume_flow / total_flow
    } else {
        1.0
    };
    let mix_density = gas_vol_fraction * gas_density + (1.0 - gas_vol_fraction) * liquid_density;
    let momentum_flux = mix_density * nozzle_velocity * nozzle_velocity;

    let bulk_efficiency = if override_bulk_efficiency > 0.0 {
        override_bulk_efficiency
    } else {
        match device.max_momentum() {
            None => 0.0,
            Some(maximum) => {
                let base = device.bulk_efficiency();
                if momentum_flux > maximum {
                    let degradation = ((momentum_flux - maximum) / maximum).min(1.0);
                    base * (1.0 - 0.5 * degradation)
                } else if momentum_flux < maximum * 0.6 * 0.3 {
                    base * 0.5
                } else {
                    base
                }
            }
        }
    };

    let pressure_drop =
        device.pressure_drop_coefficient() * 0.5 * mix_density * nozzle_velocity * nozzle_velocity;

    let downstream = if device == InletDevice::None {
        incoming
    } else {
        let multiplier = if override_dsd_multiplier > 0.0 {
            override_dsd_multiplier
        } else {
            device.dsd_multiplier()
        };
        let diameter = (incoming.0 * multiplier).max(1.0e-6);
        let mut spread = incoming.1;
        if matches!(device, InletDevice::InletCyclone | InletDevice::InletVane) {
            spread = (spread * 1.2).min(4.0);
        }
        if matches!(device, InletDevice::HalfPipe | InletDevice::DeflectorPlate)
            && momentum_flux > 3000.0
        {
            spread = (spread * 0.8).max(1.5);
        }
        (diameter, spread, incoming.2)
    };

    InletResult {
        bulk_efficiency,
        pressure_drop,
        nozzle_velocity,
        momentum_flux,
        downstream,
    }
}

#[cfg(test)]
mod tests {
    use super::{InletDevice, calculate};

    /// The enhanced capture's default device is `HALF_PIPE`; its bulk efficiency is the
    /// device's own `0.70` at the state's momentum.
    #[test]
    fn the_half_pipe_bulk_efficiency_matches_the_capture() {
        let result = calculate(
            InletDevice::HalfPipe,
            (100.0e-6, 2.6, false),
            50.0,
            750.0,
            // The inlet-pipe flows: `gasVelocity * pipeArea` each, not the vessel's.
            0.062_831_853_071_795_86,
            0.062_831_853_071_795_86,
            0.2,
            0.0,
            0.0,
        );
        // The capture's two-phase row: velocity 4.0 m/s, momentum 6400 Pa, efficiency
        // 0.4100000001840002 (the class's own floating-point path).
        assert!((result.nozzle_velocity - 4.0).abs() < 1.0e-12);
        assert!((result.momentum_flux - 6400.0).abs() < 1.0e-9);
        assert!((result.bulk_efficiency - 0.41).abs() < 1.0e-12);
        assert!((result.downstream.0 - 75.0e-6).abs() < 1.0e-15);
    }
}
