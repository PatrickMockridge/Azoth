//! DEXPI reference units <-> SI, as one explicit table and a refusal.
//!
//! A DEXPI physical quantity names its unit by an object reference such as
//! `Core/PhysicalQuantities.MassFlowRateUnit.KilogramPerHour`, and the Process profile's
//! writers use kg/h, bar absolute and °C. Azoth's vocabulary is SI, so an exchange converts.
//!
//! **Nothing here is inferred from a name.** A reference not in the table is refused, gauge
//! pressure is refused rather than given an atmosphere, and Fahrenheit and Rankine are refused
//! by name - a unit this port has not reviewed is a value it cannot carry.

/// What a reference unit measures.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum QuantityKind {
    MassFlow,
    PressureAbsolute,
    Temperature,
    Length,
}

/// One reference unit and its two conversions.
struct UnitRef {
    reference: &'static str,
    kind: QuantityKind,
    to_si: fn(f64) -> f64,
    from_si: fn(f64) -> f64,
}

fn identity(value: f64) -> f64 {
    value
}

fn hours_to_seconds(value: f64) -> f64 {
    value / 3600.0
}

fn seconds_to_hours(value: f64) -> f64 {
    value * 3600.0
}

fn bar_to_pascal(value: f64) -> f64 {
    value * 1.0e5
}

fn pascal_to_bar(value: f64) -> f64 {
    value / 1.0e5
}

fn celsius_to_kelvin(value: f64) -> f64 {
    value + 273.15
}

fn kelvin_to_celsius(value: f64) -> f64 {
    value - 273.15
}

fn millimetres_to_metres(value: f64) -> f64 {
    value / 1000.0
}

fn metres_to_millimetres(value: f64) -> f64 {
    value * 1000.0
}

/// The units this port can convert.
static UNITS: &[UnitRef] = &[
    UnitRef {
        reference: "Core/PhysicalQuantities.MassFlowRateUnit.KilogramPerHour",
        kind: QuantityKind::MassFlow,
        to_si: hours_to_seconds,
        from_si: seconds_to_hours,
    },
    UnitRef {
        reference: "Core/PhysicalQuantities.MassFlowRateUnit.KilogramPerSecond",
        kind: QuantityKind::MassFlow,
        to_si: identity,
        from_si: identity,
    },
    UnitRef {
        reference: "Core/PhysicalQuantities.PressureAbsoluteUnit.Bar",
        kind: QuantityKind::PressureAbsolute,
        to_si: bar_to_pascal,
        from_si: pascal_to_bar,
    },
    UnitRef {
        reference: "Core/PhysicalQuantities.PressureAbsoluteUnit.Pascal",
        kind: QuantityKind::PressureAbsolute,
        to_si: identity,
        from_si: identity,
    },
    UnitRef {
        reference: "Core/PhysicalQuantities.TemperatureUnit.DegreeCelsius",
        kind: QuantityKind::Temperature,
        to_si: celsius_to_kelvin,
        from_si: kelvin_to_celsius,
    },
    UnitRef {
        reference: "Core/PhysicalQuantities.TemperatureUnit.Kelvin",
        kind: QuantityKind::Temperature,
        to_si: identity,
        from_si: identity,
    },
    UnitRef {
        reference: "Core/PhysicalQuantities.LengthUnit.Millimetre",
        kind: QuantityKind::Length,
        to_si: millimetres_to_metres,
        from_si: metres_to_millimetres,
    },
];

/// Why a value could not be converted.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum UnitError {
    /// The reference is not one this port has reviewed.
    Unsupported { reference: String },
    /// The reference is a gauge pressure, which needs an atmosphere this exchange does not carry.
    Gauge { reference: String },
    /// The value is not finite.
    NonFinite { reference: String, value: String },
}

impl UnitError {
    /// A one-line message naming what was refused.
    #[must_use]
    pub fn message(&self) -> String {
        match self {
            Self::Unsupported { reference } => format!(
                "`{reference}` is not a unit this port converts; the table in azoth-dexpi's \
                 units.rs names the reviewed references"
            ),
            Self::Gauge { reference } => format!(
                "`{reference}` is a gauge pressure and this exchange carries no atmosphere to \
                 add; write the absolute pressure instead"
            ),
            Self::NonFinite { reference, value } => {
                format!("`{value}` under `{reference}` is not a finite number")
            }
        }
    }
}

/// The kind a reference measures, if it is a reviewed unit.
#[must_use]
pub fn kind(reference: &str) -> Option<QuantityKind> {
    UNITS
        .iter()
        .find(|unit| unit.reference == reference)
        .map(|unit| unit.kind)
}

/// Convert a DEXPI magnitude to the SI base of its kind.
///
/// # Errors
/// [`UnitError`] for an unreviewed reference, a gauge pressure, or a non-finite value.
pub fn to_si(reference: &str, value: f64) -> Result<f64, UnitError> {
    if reference.contains("PressureGaugeUnit") {
        return Err(UnitError::Gauge {
            reference: reference.to_string(),
        });
    }
    if !value.is_finite() {
        return Err(UnitError::NonFinite {
            reference: reference.to_string(),
            value: format!("{value}"),
        });
    }
    let unit = UNITS
        .iter()
        .find(|unit| unit.reference == reference)
        .ok_or_else(|| UnitError::Unsupported {
            reference: reference.to_string(),
        })?;
    Ok((unit.to_si)(value))
}

/// Convert an SI magnitude to a DEXPI reference unit.
///
/// # Errors
/// [`UnitError`] for an unreviewed reference or a non-finite value.
pub fn from_si(reference: &str, value: f64) -> Result<f64, UnitError> {
    if !value.is_finite() {
        return Err(UnitError::NonFinite {
            reference: reference.to_string(),
            value: format!("{value}"),
        });
    }
    let unit = UNITS
        .iter()
        .find(|unit| unit.reference == reference)
        .ok_or_else(|| UnitError::Unsupported {
            reference: reference.to_string(),
        })?;
    Ok((unit.from_si)(value))
}
