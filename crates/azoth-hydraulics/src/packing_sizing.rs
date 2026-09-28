//! `hydraulics.packing_sizing` - the internal diameter a packed bed needs to run at a chosen
//! fraction of flood.
//!
//! ```text
//! u_design = u_flood * f
//! q_v      = m_v / max(rho_v, 0.01)
//! A        = q_v / u_design
//! d_raw    = sqrt(4 A / pi)
//! d        = the first standard vessel size >= d_raw
//! ```
//!
//! Spec: `specs/calcs/hydraulics/packing_sizing.toml`
//!
//! # The trial diameter is inert, and that is what makes this a closed form
//!
//! `PackingHydraulicsCalculator.sizeColumnDiameter` sets `columnDiameter = 1.0` and then calls
//! `calculateFloodingVelocity()` - which reads **no diameter and no area**, only the packing's
//! factor, the flows, the densities and the liquid's viscosity. So the trial is set and never
//! used, and sizing is not the iteration the name `trial` suggests: one fit evaluation, one
//! area, one rounding.
//!
//! # The rounding is up, and its floor is reached
//!
//! `roundToStandardDiameter` returns the first size in a fixed table of thirty-one that is
//! `>= d_raw`, and only above `8.0` m falls back to half-metre steps. Nothing smaller than
//! `0.3` m is ever returned, so a requirement of `0.1195` m and one of `0.2973` m both answer
//! `0.3` m - which is why the raw diameter is reported beside the rounded one rather than
//! instead of it.

use azoth_core::units::{DynamicViscosity, MassDensity, MassRate, meters};
use azoth_core::{Result, apply_checks};

use crate::packing::packing_or_default;
use crate::packing_hydraulics::eckert_flooding_velocity;
use crate::results::PackingSizingResult;
use crate::spec_gen;

/// The standard vessel sizes `roundToStandardDiameter` chooses from, in metres.
///
/// Thirty-one of them, in hundredths to `1.2`, then in fifths and tenths to `4.0`, then in
/// halves to `8.0`. Above the table the class falls back to `ceil(2 d)/2`, which is the same
/// half-metre step continuing - so the table is a list of what is *bought* rather than a
/// discretisation of what is possible.
const STANDARD_DIAMETERS_M: [f64; 31] = [
    0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.4, 1.5, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8,
    3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.5, 5.0, 5.5, 6.0, 7.0, 8.0,
];

/// The diameter `sizeColumnDiameter` starts from, and the one it answers with when the design
/// velocity is not positive.
const TRIAL_DIAMETER_M: f64 = 1.0;

/// The floor the vapour density is held at in the volumetric flow, in kg/m**3.
const MIN_VAPOR_DENSITY: f64 = 0.01;

/// The state the sizing reads, which is a *subset* of what the full hydraulics reads.
///
/// Deliberately not the whole [`crate::packing_hydraulics::PackingState`]: this calculation
/// touches no surface tension, no diffusivity and no packed height, and taking them would say
/// it does.
#[derive(Debug, Clone, Copy)]
pub struct PackingSizingState {
    /// The fraction of the flooding velocity the column is designed to run at.
    pub design_flood_fraction: f64,
    /// The vapour's mass flow.
    pub vapor_mass_flow: MassRate,
    /// The liquid's mass flow, which enters only through the flow parameter.
    pub liquid_mass_flow: MassRate,
    /// The vapour's mass density.
    pub vapor_density: MassDensity,
    /// The liquid's mass density.
    pub liquid_density: MassDensity,
    /// The liquid's dynamic viscosity.
    pub liquid_viscosity: DynamicViscosity,
    /// The packing's relative hydraulic capacity factor. One is the class's default.
    pub hydraulic_capacity_factor: f64,
}

/// The internal diameter a packed bed needs at a chosen fraction of flood.
///
/// # Errors
/// * [`azoth_core::AzothError::OutOfRange`] if the design flood fraction, either density or the
///   liquid viscosity is not positive.
///
/// # Example
/// ```
/// use azoth_core::units::{
///     kilograms_per_cubic_meter, kilograms_per_second, pascal_seconds,
/// };
/// use azoth_hydraulics::packing_sizing::{PackingSizingState, packing_sizing};
///
/// let r = packing_sizing("Pall-Ring-50", PackingSizingState {
///     design_flood_fraction: 0.7,
///     vapor_mass_flow: kilograms_per_second(0.35),
///     liquid_mass_flow: kilograms_per_second(3.5),
///     vapor_density: kilograms_per_cubic_meter(45.0),
///     liquid_density: kilograms_per_cubic_meter(990.0),
///     liquid_viscosity: pascal_seconds(6.5e-4),
///     hydraulic_capacity_factor: 1.0,
/// })?;
/// assert!((r.required_diameter.value - 0.371_717_652_150_071_8).abs() < 1e-12);
/// assert!((r.column_diameter.value - 0.4).abs() < 1e-12);
/// # Ok::<(), azoth_core::AzothError>(())
/// ```
pub fn packing_sizing(name: &str, state: PackingSizingState) -> Result<PackingSizingResult> {
    let spec = &spec_gen::PACKING_SIZING_SPEC;
    let mut warnings = Vec::new();

    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "design_flood_fraction" => Some(state.design_flood_fraction),
            "vapor_mass_flow" => Some(state.vapor_mass_flow.value),
            "liquid_mass_flow" => Some(state.liquid_mass_flow.value),
            "vapor_density" => Some(state.vapor_density.value),
            "liquid_density" => Some(state.liquid_density.value),
            "liquid_viscosity" => Some(state.liquid_viscosity.value),
            "hydraulic_capacity_factor" => Some(state.hydraulic_capacity_factor),
            _ => None,
        },
        &mut warnings,
    )?;

    let packing = packing_or_default(name);

    let flooding_velocity = eckert_flooding_velocity(
        packing.packing_factor,
        state.hydraulic_capacity_factor,
        state.vapor_mass_flow.value,
        state.liquid_mass_flow.value,
        state.vapor_density.value,
        state.liquid_density.value,
        state.liquid_viscosity.value,
    );

    let design_velocity = flooding_velocity * state.design_flood_fraction;
    // **The class answers its trial diameter here rather than dividing by zero.** `if (uDesign
    // <= 0) return 1.0;` is a guard, not a physical range, so the four intermediates are the
    // trial's own state and the answer is the trial rather than an error.
    let (vapor_volumetric_flow, required_area, required_diameter, column_diameter) =
        if design_velocity <= 0.0 {
            (0.0, 0.0, TRIAL_DIAMETER_M, TRIAL_DIAMETER_M)
        } else {
            let volumetric =
                state.vapor_mass_flow.value / state.vapor_density.value.max(MIN_VAPOR_DENSITY);
            let area = volumetric / design_velocity;
            let raw = (4.0 * area / std::f64::consts::PI).sqrt();
            (volumetric, area, raw, round_to_standard_diameter(raw))
        };

    Ok(PackingSizingResult {
        packing_name: packing.name.clone(),
        packing_factor: packing.packing_factor,
        flooding_velocity,
        design_velocity,
        vapor_volumetric_flow,
        required_area: uom::si::f64::Area::new::<uom::si::area::square_meter>(required_area),
        required_diameter: meters(required_diameter),
        column_diameter: meters(column_diameter),
        warnings,
    })
}

/// `roundToStandardDiameter`: the first standard size at or above `diameter`.
fn round_to_standard_diameter(diameter: f64) -> f64 {
    for size in STANDARD_DIAMETERS_M {
        if size >= diameter {
            return size;
        }
    }
    (diameter * 2.0).ceil() / 2.0
}

#[cfg(test)]
mod tests {
    use super::*;
    use azoth_core::units::{kilograms_per_cubic_meter, kilograms_per_second, pascal_seconds};

    fn absorber() -> PackingSizingState {
        PackingSizingState {
            design_flood_fraction: 0.7,
            vapor_mass_flow: kilograms_per_second(0.35),
            liquid_mass_flow: kilograms_per_second(3.5),
            vapor_density: kilograms_per_cubic_meter(45.0),
            liquid_density: kilograms_per_cubic_meter(990.0),
            liquid_viscosity: pascal_seconds(6.5e-4),
            hydraulic_capacity_factor: 1.0,
        }
    }

    /// The table's floor is a bound on the answer, not a formality: two raw diameters nearly
    /// three times apart both buy the same vessel.
    #[test]
    fn the_floor_collapses_two_different_requirements() {
        assert!((round_to_standard_diameter(0.119_496_391_078_763_78) - 0.3).abs() < 1e-12);
        assert!((round_to_standard_diameter(0.297_311_141_804_488_14) - 0.3).abs() < 1e-12);
    }

    /// Above the table the step continues at half a metre, which is the `ceil(2d)/2` fallback.
    #[test]
    fn above_the_table_the_step_is_half_a_metre() {
        assert!((round_to_standard_diameter(8.1) - 8.5).abs() < 1e-12);
        assert!((round_to_standard_diameter(8.0) - 8.0).abs() < 1e-12);
    }

    /// **Two degenerate states, two different answers, and the difference is the fit's floor.**
    ///
    /// `FLV` is clamped to `0.005` rather than to zero, so a *zero vapour flow* still gets a
    /// positive flooding velocity: the area is zero, the raw diameter is zero, and the answer is
    /// the table's `0.3` m floor. Only a zero **capacity factor** makes the flood itself zero,
    /// which is the state `if (uDesign <= 0) return 1.0;` exists for - and there the answer is
    /// the trial diameter, not the floor.
    #[test]
    fn a_zero_flood_answers_the_trial_and_a_zero_flow_answers_the_floor() {
        let mut state = absorber();
        state.design_flood_fraction = -1.0;
        assert!(
            packing_sizing("Pall-Ring-50", state).is_err(),
            "the range check is an error"
        );

        let no_vapour = packing_sizing(
            "Pall-Ring-50",
            PackingSizingState {
                vapor_mass_flow: kilograms_per_second(0.0),
                ..absorber()
            },
        )
        .expect("a zero vapour flow is a warning, not an error");
        assert!(
            no_vapour.flooding_velocity > 0.0,
            "the fit's FLV floor keeps the flood alive"
        );
        assert!((no_vapour.required_diameter.value - 0.0).abs() < 1e-12);
        assert!((no_vapour.column_diameter.value - 0.3).abs() < 1e-12);

        let no_flood = packing_sizing(
            "Pall-Ring-50",
            PackingSizingState {
                hydraulic_capacity_factor: 0.0,
                ..absorber()
            },
        )
        .expect("a zero capacity factor is a state, not an error");
        assert!((no_flood.design_velocity - 0.0).abs() < 1e-12);
        assert!((no_flood.column_diameter.value - TRIAL_DIAMETER_M).abs() < 1e-12);
    }
}
