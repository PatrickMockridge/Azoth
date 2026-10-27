//! `ExpanderChartKhader` - a turbo-expander map over the velocity ratio and the IGV opening.
//!
//! **The abscissa is `U/C`, not flow.** Each measured IGV opening is one curve of isentropic
//! efficiency and stage head against the blade-speed/spouting-velocity ratio; the chart
//! interpolates **bilinearly** - linearly along each curve, then linearly between the two
//! bracketing openings - and the U/C direction is **linearly extrapolated** beyond the
//! measured ends where the IGV direction is **clamped** to the first or last curve. Both
//! rules are the class's and differ on purpose.
//!
//! **The head is stored as a coefficient over the reference sound speed.** `setCurves`
//! divides every stage head (J/kg) by `cs_ref²`, so the stored curve is dimensionless; a read
//! multiplies it by the *local* sound speed squared, which is how the same map scales to a
//! different gas. Azoth has no real-fluid sound-speed model, so both speeds are stated
//! inputs - the reference's phase sound speed normalises the map, and the process's is the
//! caller's at evaluation.
//!
//! **The map is built from unordered points.** The class sorts each U/C curve by its ratio
//! (the head and efficiency travel with it) and leaves the IGV order to the caller's javadoc
//! promise that it is strictly increasing; this refuses an unordered one rather than reading
//! past it.

use azoth_core::{AzothError, Result};

/// The class's degree-2 polynomial fit's parameter count is not used here; the chart is
/// piecewise linear by construction.
const LEGACY_PARABOLA: f64 = -3.56;
/// `getOptimumVelocityRatio`'s default when no map is defined.
const DEFAULT_OPTIMUM_UC: f64 = 0.7;

/// One measured point: an IGV opening, a velocity ratio, and the efficiency and stage head.
///
/// `igv` and `uc` are fractions; `efficiency` is a fraction; `head` is a specific energy in
/// J/kg, the class's kJ/kg multiplied by `1000`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CurvePoint {
    /// The inlet-guide-vane opening, 0..1.
    pub igv: f64,
    /// The blade-speed/spouting-velocity ratio.
    pub uc: f64,
    /// The measured isentropic efficiency, as a fraction.
    pub efficiency: f64,
    /// The measured stage head, J/kg.
    pub head: f64,
}

/// One IGV curve, sorted by U/C, as the chart stores it.
#[derive(Debug, Clone, PartialEq)]
struct Curve {
    uc: Vec<f64>,
    efficiency: Vec<f64>,
    head_coefficient: Vec<f64>,
}

/// The bilinear map and the reference sound speed.
#[derive(Debug, Clone, PartialEq)]
pub struct ExpanderChart {
    reference_sound_speed: f64,
    igv_positions: Vec<f64>,
    curves: Vec<Curve>,
}

impl ExpanderChart {
    /// Build the map from a flattened point list and the reference sound speed.
    ///
    /// **The rows are grouped by `igv` in order**, each curve sorted by U/C; the stage heads
    /// are divided by `reference_sound_speed²` once, here.
    ///
    /// # Errors
    /// [`AzothError::InvalidInput`] for no points, a non-finite value, a non-positive
    /// reference sound speed, or IGV openings that do not strictly increase in order.
    pub fn new(points: &[CurvePoint], reference_sound_speed: f64) -> Result<Self> {
        if points.is_empty() {
            return Err(AzothError::invalid_input(
                "chart",
                "an expander map needs at least one measured point",
            ));
        }
        if !reference_sound_speed.is_finite() || reference_sound_speed <= 0.0 {
            return Err(AzothError::invalid_input(
                "reference_sound_speed",
                format!(
                    "the map's head coefficient divides by the reference sound speed squared, \
                     and {reference_sound_speed} m/s is not positive"
                ),
            ));
        }
        for (index, point) in points.iter().enumerate() {
            for (name, value) in [
                ("igv", point.igv),
                ("uc", point.uc),
                ("efficiency", point.efficiency),
                ("head", point.head),
            ] {
                if !value.is_finite() {
                    return Err(AzothError::invalid_input(
                        "chart",
                        format!("point {index}'s {name} is {value}, not a finite number"),
                    ));
                }
            }
        }

        let mut igv_positions: Vec<f64> = Vec::new();
        let mut groups: Vec<Vec<CurvePoint>> = Vec::new();
        for point in points {
            match igv_positions.last() {
                Some(last) if *last == point.igv => groups
                    .last_mut()
                    .expect("a group per position")
                    .push(*point),
                _ => {
                    if let Some(last) = igv_positions.last() {
                        if point.igv <= *last {
                            return Err(AzothError::invalid_input(
                                "chart",
                                format!(
                                    "the IGV openings are not strictly increasing: \
                                     {} follows {last}",
                                    point.igv
                                ),
                            ));
                        }
                    }
                    igv_positions.push(point.igv);
                    groups.push(vec![*point]);
                }
            }
        }

        let mut curves = Vec::with_capacity(groups.len());
        for group in groups {
            let mut rows = group;
            rows.sort_by(|a, b| a.uc.partial_cmp(&b.uc).expect("finite U/C"));
            curves.push(Curve {
                uc: rows.iter().map(|row| row.uc).collect(),
                efficiency: rows.iter().map(|row| row.efficiency).collect(),
                head_coefficient: rows
                    .iter()
                    .map(|row| row.head / (reference_sound_speed * reference_sound_speed))
                    .collect(),
            });
        }

        Ok(Self {
            reference_sound_speed,
            igv_positions,
            curves,
        })
    }

    /// The reference sound speed the map was normalised with, m/s.
    pub fn reference_sound_speed(&self) -> f64 {
        self.reference_sound_speed
    }

    /// `getEfficiency(uc, igv)`, as a fraction.
    pub fn efficiency(&self, uc: f64, igv: f64) -> f64 {
        self.interpolate(uc, igv, |curve| &curve.efficiency)
    }

    /// `getStageHeadDrop(uc, igv, localCs)`, J/kg.
    ///
    /// The class returns kJ/kg and takes a `SystemInterface`; this takes the local sound
    /// speed the caller read and answers SI.
    pub fn stage_head_drop(&self, uc: f64, igv: f64, local_sound_speed: f64) -> f64 {
        let coefficient = self.interpolate(uc, igv, |curve| &curve.head_coefficient);
        coefficient * local_sound_speed * local_sound_speed
    }

    /// `getOptimumVelocityRatio(igv)`: the U/C at peak efficiency on the nearest curve.
    pub fn optimum_velocity_ratio(&self, igv: f64) -> f64 {
        if self.curves.is_empty() {
            return DEFAULT_OPTIMUM_UC;
        }
        let mut nearest = 0;
        let mut difference = (self.igv_positions[0] - igv).abs();
        for (index, position) in self.igv_positions.iter().enumerate().skip(1) {
            let candidate = (position - igv).abs();
            if candidate < difference {
                difference = candidate;
                nearest = index;
            }
        }
        let curve = &self.curves[nearest];
        let mut best = 0;
        for (index, efficiency) in curve.efficiency.iter().enumerate().skip(1) {
            if *efficiency > curve.efficiency[best] {
                best = index;
            }
        }
        curve.uc[best]
    }

    fn interpolate(&self, uc: f64, igv: f64, grid: impl Fn(&Curve) -> &Vec<f64>) -> f64 {
        if self.curves.is_empty() {
            return LEGACY_PARABOLA * (uc - 1.0) * (uc - 1.0) + 1.0;
        }
        if self.curves.len() == 1 || igv <= self.igv_positions[0] {
            return interpolate_curve(uc, &self.curves[0].uc, grid(&self.curves[0]));
        }
        let last = self.curves.len() - 1;
        if igv >= self.igv_positions[last] {
            return interpolate_curve(uc, &self.curves[last].uc, grid(&self.curves[last]));
        }
        let mut lower = 0;
        for index in 0..last {
            if igv >= self.igv_positions[index] && igv <= self.igv_positions[index + 1] {
                lower = index;
                break;
            }
        }
        let v0 = interpolate_curve(uc, &self.curves[lower].uc, grid(&self.curves[lower]));
        let v1 = interpolate_curve(
            uc,
            &self.curves[lower + 1].uc,
            grid(&self.curves[lower + 1]),
        );
        let t = (igv - self.igv_positions[lower])
            / (self.igv_positions[lower + 1] - self.igv_positions[lower]);
        v0 + t * (v1 - v0)
    }
}

/// `interpolateCurve`: piecewise linear, **extrapolating** beyond the measured U/C ends.
fn interpolate_curve(uc: f64, x: &[f64], y: &[f64]) -> f64 {
    let n = x.len();
    if n == 1 {
        return y[0];
    }
    if uc <= x[0] {
        let slope = (y[1] - y[0]) / (x[1] - x[0]);
        return y[0] + slope * (uc - x[0]);
    }
    if uc >= x[n - 1] {
        let slope = (y[n - 1] - y[n - 2]) / (x[n - 1] - x[n - 2]);
        return y[n - 1] + slope * (uc - x[n - 1]);
    }
    for index in 0..n - 1 {
        if uc >= x[index] && uc <= x[index + 1] {
            let t = (uc - x[index]) / (x[index + 1] - x[index]);
            return y[index] + t * (y[index + 1] - y[index]);
        }
    }
    y[n - 1]
}

#[cfg(test)]
mod tests {
    use super::*;

    /// `TurboExpanderEnhancementsTest.buildChart`'s table, head in kJ/kg.
    fn table() -> Vec<CurvePoint> {
        let igv = [0.6, 1.0];
        let uc = [0.5, 0.6, 0.7, 0.8, 0.9];
        let eta = [
            [0.72, 0.80, 0.84, 0.82, 0.76],
            [0.78, 0.85, 0.88, 0.86, 0.80],
        ];
        let head = [
            [30.0, 32.0, 33.0, 32.5, 31.0],
            [34.0, 36.0, 37.0, 36.5, 35.0],
        ];
        let mut points = Vec::new();
        for i in 0..2 {
            for j in 0..5 {
                points.push(CurvePoint {
                    igv: igv[i],
                    uc: uc[j],
                    efficiency: eta[i][j],
                    head: head[i][j] * 1000.0,
                });
            }
        }
        points
    }

    #[test]
    fn the_capture_rows_reproduce_the_class() {
        // `validation/neqsim/captures/process_expander_chart.tsv`, reference fluid absent:
        // the class's `csRef` falls back to `1.0`, so the head round-trips in J/kg.
        let chart = ExpanderChart::new(&table(), 1.0).expect("a valid map");
        for (uc, igv, expected) in [
            (0.7, 1.0, 0.88),
            (0.7, 0.6, 0.84),
            (0.7, 0.8, 0.86),
            (0.4, 1.0, 0.7100000000000001),
            (0.95, 1.0, 0.7700000000000001),
        ] {
            assert!(
                (chart.efficiency(uc, igv) - expected).abs() < 1.0e-12,
                "efficiency at {uc}/{igv}: {} against {expected}",
                chart.efficiency(uc, igv)
            );
        }
        assert!((chart.optimum_velocity_ratio(1.0) - 0.7).abs() < 1.0e-12);
        assert!((chart.optimum_velocity_ratio(0.6) - 0.7).abs() < 1.0e-12);
        for (uc, igv, expected) in [
            (0.7, 1.0, 37_000.0),
            (0.7, 0.6, 33_000.0),
            (0.7, 0.8, 35_000.0),
            (0.4, 1.0, 32_000.0),
            (0.95, 1.0, 34_250.0),
        ] {
            assert!(
                (chart.stage_head_drop(uc, igv, 1.0) - expected).abs() < 1.0e-6,
                "head at {uc}/{igv}: {} against {expected}",
                chart.stage_head_drop(uc, igv, 1.0)
            );
        }
    }

    #[test]
    fn the_head_scales_with_the_local_sound_speed() {
        // The capture's reference phase sound speed and two system-level ones.
        let chart = ExpanderChart::new(&table(), 415.419842626026).expect("a valid map");
        let reference = chart.stage_head_drop(0.7, 1.0, 830.839685252052);
        let carbon_dioxide = chart.stage_head_drop(0.7, 1.0, 459.52905504799423);
        assert!((reference - 148_000.00000000003).abs() < 1.0e-6);
        assert!((carbon_dioxide - 45_274.45202145469).abs() < 1.0e-6);
    }

    #[test]
    fn an_unordered_igv_list_is_refused() {
        let points = [
            CurvePoint {
                igv: 1.0,
                uc: 0.5,
                efficiency: 0.8,
                head: 10_000.0,
            },
            CurvePoint {
                igv: 0.6,
                uc: 0.5,
                efficiency: 0.78,
                head: 9_000.0,
            },
        ];
        let error = ExpanderChart::new(&points, 400.0).expect_err("the IGV list is backwards");
        assert_eq!(error.field(), Some("chart"));
    }
}
