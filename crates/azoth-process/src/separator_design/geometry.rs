//! `SeparatorGeometryCalculator`: the areas, settling heights and residence times the
//! enhanced entrainment chain reads.
//!
//! The calculator is a closed form: a horizontal vessel's liquid area is the circular segment
//! below the normal level, its gas area the rest of the cross-section; a vertical vessel uses
//! the full section for both and its heights are `level * length` and the remainder. The
//! three-phase horizontal branch splits the liquid into an oil pad and a water layer by the
//! **flow ratio**, not by the stated level, and its gravity section ends at the weir.
//!
//! The source is `process/equipment/separator/entrainment/SeparatorGeometryCalculator.java`,
//! through `calculate`, `calculateThreePhase` and `calcSegmentArea`.

/// The vessel the geometry is read on.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GeometryInput {
    /// Internal diameter, m.
    pub diameter: f64,
    /// Tangent-to-tangent length, m.
    pub length: f64,
    /// Whether the vessel is horizontal.
    pub horizontal: bool,
    /// The normal liquid level as a fraction of the diameter (horizontal) or length
    /// (vertical).
    pub level_fraction: f64,
    /// The inlet device's position as a fraction of the length. Default `0.05`.
    pub inlet_position: f64,
    /// The mist eliminator's position as a fraction of the length. Default `0.90`.
    pub mist_position: f64,
    /// The weir's position as a fraction of the length. Default `0.60`.
    pub weir_position: f64,
}

impl Default for GeometryInput {
    fn default() -> Self {
        Self {
            diameter: 2.0,
            length: 8.0,
            horizontal: true,
            level_fraction: 0.5,
            inlet_position: 0.05,
            mist_position: 0.90,
            weir_position: 0.60,
        }
    }
}

/// The calculator's results, in SI.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Geometry {
    pub gas_area: f64,
    pub liquid_area: f64,
    pub total_area: f64,
    pub gas_volume: f64,
    pub liquid_volume: f64,
    pub total_volume: f64,
    pub effective_gas_settling_height: f64,
    pub effective_liquid_settling_height: f64,
    pub gas_residence_time: f64,
    pub liquid_residence_time: f64,
    pub gravity_section_length: f64,
    pub oil_pad_thickness: f64,
    pub water_layer_height: f64,
}

/// `calculate`: the two-phase geometry.
#[must_use]
pub fn two_phase(input: &GeometryInput, gas_volume_flow: f64, liquid_volume_flow: f64) -> Geometry {
    let total_area = std::f64::consts::PI * input.diameter * input.diameter / 4.0;
    if input.horizontal {
        let liquid_height = input.level_fraction * input.diameter;
        let liquid_area = segment_area(input.diameter, liquid_height);
        let gas_area = total_area - liquid_area;
        let gas_settling = input.diameter - liquid_height;
        let liquid_settling = liquid_height;
        let gravity_length = input.length * (input.mist_position - input.inlet_position);
        let gas_volume = gas_area * input.length;
        let liquid_volume = liquid_area * input.length;
        let gas_residence = if gas_area > 0.0 && gas_volume_flow > 0.0 {
            gravity_length / (gas_volume_flow / gas_area)
        } else {
            0.0
        };
        let liquid_residence = if liquid_volume > 0.0 && liquid_volume_flow > 0.0 {
            liquid_volume / liquid_volume_flow
        } else {
            0.0
        };
        Geometry {
            gas_area,
            liquid_area,
            total_area,
            gas_volume,
            liquid_volume,
            total_volume: total_area * input.length,
            effective_gas_settling_height: gas_settling,
            effective_liquid_settling_height: liquid_settling,
            gas_residence_time: gas_residence,
            liquid_residence_time: liquid_residence,
            gravity_section_length: gravity_length,
            oil_pad_thickness: 0.0,
            water_layer_height: 0.0,
        }
    } else {
        let liquid_height = input.level_fraction * input.length;
        let gas_height = input.length - liquid_height;
        let gas_settling = input.diameter;
        let gravity_length = (input.length * input.mist_position - liquid_height).max(0.0);
        let gas_volume = total_area * gas_height;
        let liquid_volume = total_area * liquid_height;
        let gas_residence = if total_area > 0.0 && gas_volume_flow > 0.0 {
            gas_height / (gas_volume_flow / total_area)
        } else {
            0.0
        };
        let liquid_residence = if liquid_volume > 0.0 && liquid_volume_flow > 0.0 {
            liquid_volume / liquid_volume_flow
        } else {
            0.0
        };
        Geometry {
            gas_area: total_area,
            liquid_area: total_area,
            total_area,
            gas_volume,
            liquid_volume,
            total_volume: total_area * input.length,
            effective_gas_settling_height: gas_settling,
            effective_liquid_settling_height: liquid_height,
            gas_residence_time: gas_residence,
            liquid_residence_time: liquid_residence,
            gravity_section_length: gravity_length,
            oil_pad_thickness: 0.0,
            water_layer_height: 0.0,
        }
    }
}

/// `calculateThreePhase`: the horizontal three-phase geometry.
///
/// **The split is the flow ratio, not the stated oil level.** `oil_level_fraction` is read
/// nowhere in this branch; the water layer and oil pad are `level * water_frac` and
/// `level * oil_frac`, and the gravity section ends at the weir. A vertical vessel delegates
/// to the two-phase branch with the two liquid flows summed.
#[must_use]
pub fn three_phase(
    input: &GeometryInput,
    gas_volume_flow: f64,
    oil_volume_flow: f64,
    water_volume_flow: f64,
    _oil_level_fraction: f64,
) -> Geometry {
    if !input.horizontal {
        return two_phase(input, gas_volume_flow, oil_volume_flow + water_volume_flow);
    }
    let total_area = std::f64::consts::PI * input.diameter * input.diameter / 4.0;
    let liquid_height = input.level_fraction * input.diameter;
    let total_liquid = oil_volume_flow + water_volume_flow;
    let water_frac = if total_liquid > 0.0 {
        water_volume_flow / total_liquid
    } else {
        0.5
    };
    let water_layer = liquid_height * water_frac;
    let oil_pad = liquid_height * (1.0 - water_frac);
    let liquid_area = segment_area(input.diameter, liquid_height);
    let gas_area = total_area - liquid_area;
    let gravity_length = input.length * (input.weir_position - input.inlet_position);
    let oil_area = liquid_area - segment_area(input.diameter, water_layer);
    let gas_volume = gas_area * input.length;
    let liquid_volume = liquid_area * input.length;
    let gas_residence = if gas_area > 0.0 && gas_volume_flow > 0.0 {
        gravity_length / (gas_volume_flow / gas_area)
    } else {
        0.0
    };
    let oil_volume = oil_area * gravity_length;
    let liquid_residence = if oil_volume_flow > 0.0 && oil_volume > 0.0 {
        oil_volume / oil_volume_flow
    } else {
        0.0
    };
    Geometry {
        gas_area,
        liquid_area,
        total_area,
        gas_volume,
        liquid_volume,
        total_volume: total_area * input.length,
        effective_gas_settling_height: input.diameter - liquid_height,
        effective_liquid_settling_height: oil_pad,
        gas_residence_time: gas_residence,
        liquid_residence_time: liquid_residence,
        gravity_section_length: gravity_length,
        oil_pad_thickness: oil_pad,
        water_layer_height: water_layer,
    }
}

/// `calcSegmentArea`: the circular segment below `height`.
#[must_use]
pub fn segment_area(diameter: f64, height: f64) -> f64 {
    if height <= 0.0 || diameter <= 0.0 {
        return 0.0;
    }
    if height >= diameter {
        return std::f64::consts::PI * diameter * diameter / 4.0;
    }
    let radius = diameter / 2.0;
    let theta = 2.0 * ((radius - height) / radius).acos();
    radius * radius / 2.0 * (theta - theta.sin())
}

/// `calcKFactor`: the Souders-Brown K at a velocity and two densities.
#[must_use]
pub fn k_factor(gas_velocity: f64, gas_density: f64, liquid_density: f64) -> f64 {
    let delta = liquid_density - gas_density;
    if delta <= 0.0 || gas_density <= 0.0 {
        return 0.0;
    }
    gas_velocity * (gas_density / delta).sqrt()
}

#[cfg(test)]
mod tests {
    use super::{GeometryInput, k_factor, segment_area, two_phase};

    /// The enhanced capture's `two_phase_horizontal` geometry: a 2.4 m x 7 m vessel at a
    /// half-full level.
    #[test]
    fn the_captures_horizontal_geometry_is_reproduced() {
        let input = GeometryInput {
            diameter: 2.4,
            length: 7.0,
            horizontal: true,
            level_fraction: 0.5,
            ..GeometryInput::default()
        };
        let gas_area = std::f64::consts::PI * 2.4 * 2.4 / 4.0 * 0.5;
        let geometry = two_phase(&input, 4.523_893_421_169_302, 4.523_893_421_169_302);
        let close = |actual: f64, expected: f64, what: &str| {
            let error = (actual - expected).abs() / expected.abs().max(f64::MIN_POSITIVE);
            assert!(error <= 1.0e-12, "{what}: {actual} against {expected}");
        };
        close(geometry.gas_area, gas_area, "gas area");
        close(geometry.liquid_area, gas_area, "liquid area");
        close(geometry.gravity_section_length, 5.95, "gravity length");
        close(geometry.effective_gas_settling_height, 1.2, "gas height");
        close(
            geometry.effective_liquid_settling_height,
            1.2,
            "liquid height",
        );
        close(geometry.gas_residence_time, 2.975, "gas residence");
        close(geometry.liquid_residence_time, 3.5, "liquid residence");
        close(segment_area(2.4, 1.2), gas_area, "the half-full segment");
        close(
            k_factor(2.0, 50.0, 750.0),
            // The capture's `0.5345224838248488`.
            0.534_522_483_824_848_8,
            "the K factor",
        );
    }
}
