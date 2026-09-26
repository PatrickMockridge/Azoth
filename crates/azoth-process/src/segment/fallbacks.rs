//! The class's own `DEFAULT_*` constants, and the record of which of them a state took.
//!
//! `RateBasedPackedColumn` substitutes a constant where a physical-property model answers
//! nothing. There are nine of them, and they fall into two kinds that behave differently.

/// The gas diffusivity the reference scaling falls back to, in m²/s.
pub const DEFAULT_GAS_DIFFUSIVITY: f64 = 1.5e-5;

/// The liquid diffusivity, likewise.
pub const DEFAULT_LIQUID_DIFFUSIVITY: f64 = 1.5e-9;

/// The surface tension, in N/m.
pub const DEFAULT_SURFACE_TENSION: f64 = 0.025;

/// The gas thermal conductivity, in W/(m·K).
pub const DEFAULT_GAS_THERMAL_CONDUCTIVITY: f64 = 0.030;

/// The liquid thermal conductivity, in W/(m·K).
pub const DEFAULT_LIQUID_THERMAL_CONDUCTIVITY: f64 = 0.60;

/// The gas heat capacity, in J/(kg·K).
pub const DEFAULT_GAS_HEAT_CAPACITY: f64 = 2200.0;

/// The liquid heat capacity, in J/(kg·K).
pub const DEFAULT_LIQUID_HEAT_CAPACITY: f64 = 4200.0;

/// The floor the snapshot lifts a gas diffusivity to, in m²/s.
pub const MIN_GAS_DIFFUSIVITY: f64 = 1.0e-7;

/// The floor it lifts a liquid diffusivity to, in m²/s.
pub const MIN_LIQUID_DIFFUSIVITY: f64 = 1.0e-12;

/// One substituted property.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Property {
    Diffusivity,
    SurfaceTension,
    ThermalConductivity,
    HeatCapacity,
}

impl Property {
    /// The name the warning carries.
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::Diffusivity => "diffusivity",
            Self::SurfaceTension => "surface tension",
            Self::ThermalConductivity => "thermal conductivity",
            Self::HeatCapacity => "heat capacity",
        }
    }
}

/// Which of the class's substitutions one segment's snapshot took.
///
/// **The diffusivity pair is taken on every state of the class's own tests, and that is
/// measured rather than assumed.** `averageDiffusivity` averages
/// `getEffectiveDiffusionCoefficient(i)` over the components whose value is positive, and
/// **NeqSim's flash path never populates that vector**: on the CO2/water absorber's own state
/// every entry of both phases is `0.0`, so the sum is empty and the constant stands in. What
/// the same state *does* carry is the pair matrix - `getDiffusionCoefficient(0, 1)` answers
/// `9.457085848059925e-7` for the gas, which is the value `eos.phase_transport` is separately
/// held to - so the film model's ratios are real even though their reference is not.
///
/// `Diffusivity.calcEffectiveDiffusionCoefficients` is the class that would close it, and the
/// capture prints the vector it leaves at zero.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct Fallbacks {
    pub diffusivity: bool,
    pub surface_tension: bool,
    pub thermal_conductivity: bool,
    pub heat_capacity: bool,
}

impl Fallbacks {
    /// The properties this segment substituted, in the order the class declares them.
    pub fn taken(&self) -> Vec<Property> {
        let mut taken = Vec::new();
        if self.diffusivity {
            taken.push(Property::Diffusivity);
        }
        if self.surface_tension {
            taken.push(Property::SurfaceTension);
        }
        if self.thermal_conductivity {
            taken.push(Property::ThermalConductivity);
        }
        if self.heat_capacity {
            taken.push(Property::HeatCapacity);
        }
        taken
    }

    /// Merge another segment's record, so a column reports each property once.
    pub fn merge(&mut self, other: Self) {
        self.diffusivity |= other.diffusivity;
        self.surface_tension |= other.surface_tension;
        self.thermal_conductivity |= other.thermal_conductivity;
        self.heat_capacity |= other.heat_capacity;
    }
}
