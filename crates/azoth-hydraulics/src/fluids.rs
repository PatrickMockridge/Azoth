//! Fluid property tables, read from `data/fluids/`.
//!
//! Mirrors `azoth.properties` on the Python side, and reads the same files -
//! one source of truth, two languages, with a cross-language test comparing the
//! parsed values. The Rust side embeds them with `include_str!`; the Python side
//! locates them at runtime.
//!
//! # Provenance
//!
//! The values are real published figures read from widely used engineering
//! tables that trace back to IAPWS-IF97 / NIST and the CRC Handbook. They are
//! marked [`VerifyStatus::Unverified`] because they have not been checked
//! against a primary formulation. That is a different claim from the Crane
//! coefficients, which are `EstimatedDummy` placeholders - see
//! [`crate::provenance`] for why the distinction is kept.
//!
//! # No extrapolation
//!
//! Outside the tabulated range the lookup raises rather than extending the
//! trend. Water's viscosity changes by a factor of six across 0-100 C, so a
//! straight-line extension past either end is a confident wrong number rather
//! than a small error.

use std::sync::OnceLock;

use azoth_core::units::{DynamicViscosity, MassDensity, kilograms_per_cubic_meter, pascal_seconds};
use azoth_core::{AzothError, Result};

use crate::provenance::VerifyStatus;

const WATER_CSV: &str = include_str!("../../../data/fluids/water.csv");
const AIR_CSV: &str = include_str!("../../../data/fluids/air.csv");

/// One tabulated row.
#[derive(Debug, Clone, PartialEq)]
pub struct FluidPoint {
    /// Temperature in degrees Celsius.
    pub temperature_c: f64,
    /// Density in kg/m^3.
    pub density_kg_m3: f64,
    /// Dynamic viscosity in Pa*s.
    pub dynamic_viscosity_pa_s: f64,
    /// Where the value came from.
    pub citation: String,
    /// How far it can be trusted.
    pub status: VerifyStatus,
}

/// A named fluid's property table.
#[derive(Debug, Clone, PartialEq)]
pub struct FluidTable {
    /// Identifier, e.g. `water`.
    pub name: String,
    /// Rows, sorted by temperature.
    pub points: Vec<FluidPoint>,
}

impl FluidTable {
    /// The tabulated temperature range, in degrees Celsius.
    #[must_use]
    pub fn temperature_range_c(&self) -> (f64, f64) {
        match (self.points.first(), self.points.last()) {
            (Some(low), Some(high)) => (low.temperature_c, high.temperature_c),
            _ => (f64::NAN, f64::NAN),
        }
    }

    /// Linear interpolation between tabulated points, exact at them.
    ///
    /// # Errors
    /// Returns [`AzothError::OutOfRange`] outside the tabulated range.
    fn interpolate(&self, temperature_c: f64, attribute: &str) -> Result<f64> {
        let (low, high) = self.temperature_range_c();
        if temperature_c.is_nan() || temperature_c < low || temperature_c > high {
            return Err(AzothError::out_of_range(
                "temperature",
                temperature_c,
                format!(
                    "outside the tabulated range {low}-{high} C for {}; this provider \
                     does not extrapolate",
                    self.name
                ),
            ));
        }

        // Exact hits first, so a tabulated point returns the tabulated value
        // rather than an interpolation that happens to land on it.
        for point in &self.points {
            if temperature_c == point.temperature_c {
                return Ok(field_of(point, attribute));
            }
        }

        for pair in self.points.windows(2) {
            let (lower, upper) = (&pair[0], &pair[1]);
            if lower.temperature_c < temperature_c && temperature_c < upper.temperature_c {
                let span = upper.temperature_c - lower.temperature_c;
                let fraction = (temperature_c - lower.temperature_c) / span;
                let lo = field_of(lower, attribute);
                let hi = field_of(upper, attribute);
                return Ok(lo + fraction * (hi - lo));
            }
        }

        // Unreachable: the range check above guarantees a bracketing pair.
        Err(AzothError::invalid_input(
            "temperature",
            format!("no bracketing points for {temperature_c} C"),
        ))
    }

    /// Density at a temperature in degrees Celsius.
    ///
    /// # Errors
    /// Returns an error outside the tabulated range.
    pub fn density_at_celsius(&self, temperature_c: f64) -> Result<MassDensity> {
        Ok(kilograms_per_cubic_meter(
            self.interpolate(temperature_c, "density")?,
        ))
    }

    /// Dynamic viscosity at a temperature in degrees Celsius.
    ///
    /// # Errors
    /// Returns an error outside the tabulated range.
    pub fn dynamic_viscosity_at_celsius(&self, temperature_c: f64) -> Result<DynamicViscosity> {
        Ok(pascal_seconds(
            self.interpolate(temperature_c, "dynamic_viscosity")?,
        ))
    }

    /// True when any row is a placeholder rather than a published value.
    #[must_use]
    pub fn has_placeholder_rows(&self) -> bool {
        self.points.iter().any(|p| p.status.is_placeholder())
    }
}

fn field_of(point: &FluidPoint, attribute: &str) -> f64 {
    match attribute {
        "density" => point.density_kg_m3,
        _ => point.dynamic_viscosity_pa_s,
    }
}

fn parse(name: &str, raw: &str) -> Result<FluidTable> {
    let body: String = raw
        .lines()
        .filter(|line| !line.trim_start().starts_with('#') && !line.trim().is_empty())
        .collect::<Vec<_>>()
        .join("\n");

    let mut reader = csv::Reader::from_reader(body.as_bytes());
    let mut points = Vec::new();
    for record in reader.records() {
        let record = record.map_err(|e| {
            AzothError::invalid_input("fluids", format!("malformed CSV row in {name}: {e}"))
        })?;
        let get = |index: usize, column: &str| -> Result<&str> {
            record.get(index).ok_or_else(|| {
                AzothError::invalid_input(
                    "fluids",
                    format!("row is missing column `{column}` in {name}"),
                )
            })
        };
        let number = |index: usize, column: &str| -> Result<f64> {
            get(index, column)?.trim().parse().map_err(|e| {
                AzothError::invalid_input(column, format!("not a number in {name}: {e}"))
            })
        };
        points.push(FluidPoint {
            temperature_c: number(0, "temperature_c")?,
            density_kg_m3: number(1, "density_kg_m3")?,
            dynamic_viscosity_pa_s: number(2, "dynamic_viscosity_pa_s")?,
            citation: get(3, "citation")?.to_string(),
            status: VerifyStatus::parse(get(4, "verify_status")?)?,
        });
    }

    points.sort_by(|a, b| a.temperature_c.total_cmp(&b.temperature_c));
    Ok(FluidTable {
        name: name.to_string(),
        points,
    })
}

/// Water at 1 atm, tabulated 0-100 C.
///
/// # Errors
/// Returns an error if the embedded table is malformed, which is a build-time
/// invariant that a test asserts.
pub fn water() -> Result<&'static FluidTable> {
    static TABLE: OnceLock<Result<FluidTable>> = OnceLock::new();
    cache(&TABLE, || parse("water", WATER_CSV))
}

/// Dry air at 1 atm, tabulated 0-100 C.
///
/// # Errors
/// Returns an error if the embedded table is malformed.
pub fn air() -> Result<&'static FluidTable> {
    static TABLE: OnceLock<Result<FluidTable>> = OnceLock::new();
    cache(&TABLE, || parse("air", AIR_CSV))
}

fn cache(
    slot: &'static OnceLock<Result<FluidTable>>,
    build: impl FnOnce() -> Result<FluidTable>,
) -> Result<&'static FluidTable> {
    match slot.get_or_init(build) {
        Ok(table) => Ok(table),
        Err(e) => Err(e.clone()),
    }
}

/// Look up a built-in fluid by name.
///
/// # Errors
/// Returns [`AzothError::InvalidInput`] for anything but `water` or `air`. An
/// error rather than a default: silently substituting water for an unknown fluid
/// would produce a plausible pressure drop for the wrong substance.
pub fn provider_for(name: &str) -> Result<&'static FluidTable> {
    match name.trim().to_lowercase().as_str() {
        "water" => water(),
        "air" => air(),
        other => Err(AzothError::invalid_input(
            "fluid",
            format!("`{other}` is not a built-in fluid; available: air, water"),
        )),
    }
}

/// Names of the built-in fluids.
#[must_use]
pub const fn available_fluids() -> &'static [&'static str] {
    &["air", "water"]
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn both_tables_parse() {
        assert!(!water().unwrap().points.is_empty());
        assert!(!air().unwrap().points.is_empty());
    }

    #[test]
    fn tables_are_sorted_and_span_the_expected_range() {
        for table in [water().unwrap(), air().unwrap()] {
            let temperatures: Vec<f64> = table.points.iter().map(|p| p.temperature_c).collect();
            assert!(
                temperatures.windows(2).all(|w| w[0] < w[1]),
                "{} is not sorted",
                table.name
            );
            assert_eq!(table.temperature_range_c(), (0.0, 100.0), "{}", table.name);
        }
    }

    #[test]
    fn tabulated_points_are_returned_exactly() {
        let table = water().unwrap();
        let point = table
            .points
            .iter()
            .find(|p| p.temperature_c == 20.0)
            .expect("water has a 20 C row");
        assert_eq!(
            table.density_at_celsius(20.0).unwrap().value,
            point.density_kg_m3
        );
        assert_eq!(
            table.dynamic_viscosity_at_celsius(20.0).unwrap().value,
            point.dynamic_viscosity_pa_s
        );
    }

    #[test]
    fn interpolation_matches_a_hand_computed_midpoint() {
        let table = water().unwrap();
        let low = table
            .points
            .iter()
            .find(|p| p.temperature_c == 20.0)
            .unwrap();
        let high = table
            .points
            .iter()
            .find(|p| p.temperature_c == 40.0)
            .unwrap();
        let expected = (low.density_kg_m3 + high.density_kg_m3) / 2.0;
        let actual = table.density_at_celsius(30.0).unwrap().value;
        assert!((actual - expected).abs() < 1e-12, "{actual} vs {expected}");
    }

    #[test]
    fn water_properties_fall_with_temperature() {
        let table = water().unwrap();
        let densities: Vec<f64> = [20.0, 40.0, 60.0, 80.0]
            .iter()
            .map(|t| table.density_at_celsius(*t).unwrap().value)
            .collect();
        assert!(
            densities.windows(2).all(|w| w[0] > w[1]),
            "water density is not decreasing: {densities:?}"
        );

        let viscosities: Vec<f64> = [20.0, 40.0, 60.0, 80.0]
            .iter()
            .map(|t| table.dynamic_viscosity_at_celsius(*t).unwrap().value)
            .collect();
        assert!(viscosities.windows(2).all(|w| w[0] > w[1]));
        assert!(
            viscosities[0] / viscosities[3] > 2.0,
            "viscosity should vary substantially over this range"
        );
    }

    #[test]
    fn outside_the_table_is_an_error_rather_than_an_extrapolation() {
        let table = water().unwrap();
        for temperature in [-10.0, 150.0, f64::NAN] {
            let err = table.density_at_celsius(temperature).unwrap_err();
            assert!(matches!(err, AzothError::OutOfRange { .. }), "{err}");
        }
    }

    #[test]
    fn unknown_fluid_is_an_error_not_a_default() {
        assert!(provider_for("unobtainium").is_err());
        assert_eq!(provider_for("  WATER ").unwrap().name, "water");
        assert_eq!(available_fluids(), &["air", "water"]);
    }

    #[test]
    fn fluid_rows_are_not_placeholders() {
        // Different claim from the fittings registry: these are real published
        // values, merely unchecked against the primary formulation. If someone
        // replaces them with invented numbers this fails, which is the point.
        for table in [water().unwrap(), air().unwrap()] {
            assert!(
                !table.has_placeholder_rows(),
                "{} has estimated_dummy rows; fluid tables should carry real \
                 published values, not placeholders",
                table.name
            );
        }
        // And they must say they are unverified rather than silently claiming more.
        for table in [water().unwrap(), air().unwrap()] {
            for point in &table.points {
                assert_eq!(point.status, VerifyStatus::Unverified);
                assert!(!point.citation.is_empty());
            }
        }
    }
}
