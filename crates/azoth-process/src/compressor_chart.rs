//! `CompressorChart`, the default fan-law map.
//!
//! **One quadratic in `Q/N` per table, exactly as the pump's chart is.** `setCurves` divides
//! every measured flow by its speed and every head by the speed squared, fits one degree-2
//! polynomial through all the reduced points, and reads it back as
//! `head = f(Q/N) N²`; the polytropic efficiency is fitted against `Q/N` with no speed
//! reduction. **The compressor's chart clamps nothing**: `getPolytropicHead` and
//! `getPolytropicEfficiency` extrapolate the quadratic wherever they are asked, and the
//! compressor clamps the efficiency downstream in `run`.
//!
//! **The efficiency's own flow column is the head's here.** `CompressorChart.setCurves` takes a
//! separate `flowPolytropicEfficiency` matrix; every oracle row states the same flows, so the
//! model carries one `flow` per point and this module fits the efficiency against it.
//!
//! **The surge and stonewall lookups are per curve, and no curve is built from the table.**
//! `setCurves` builds no boundary curve; `getSurgeFlowAtSpeed` and its siblings take the
//! measured curve closest to the requested speed and read its minimum or maximum flow. The
//! spline boundary curves (`SafeSplineSurgeCurve`/`SafeSplineStoneWallCurve`) belong to the
//! anti-surge machinery and are not here.
//!
//! **The dead fan-law correction is not repeated.** `setCurves` fits a third polynomial over
//! `speed/speed[0]` that no getter reads.

use azoth_core::{AzothError, Result};
use azoth_reactions::linalg::solve_least_squares;

/// The class's degree-2 polynomial fit has three coefficients.
const COEFFICIENTS: usize = 3;
/// `getSpeedValue`'s own Newton constants.
const MAX_NEWTON: usize = 50;
const SPEED_TOLERANCE: f64 = 1.0e-6;
const DAMPING: f64 = 0.7;
const MAX_STEP_FRACTION: f64 = 0.3;
/// The fallback's bisection cap and interval stop.
const MAX_BISECTION: usize = 50;
const MIN_SPEED_INTERVAL: f64 = 1.0;
/// `getFlow`'s Newton cap and tolerance.
const MAX_FLOW_ITERATIONS: usize = 100;
const FLOW_TOLERANCE: f64 = 1.0e-6;

/// One measured point of the map: a flow and the head and polytropic efficiency there.
///
/// `head` is a specific energy in J/kg and `efficiency` a fraction; the class's table is
/// normally quoted in kJ/kg and percent, and both conversions are linear, so the fit is the
/// same polynomial. `speed` is rpm and `flow` m3/hr, the class's own coordinates.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CurvePoint {
    /// The compressor speed, rpm.
    pub speed: f64,
    /// The measured flow, m3/hr.
    pub flow: f64,
    /// The measured polytropic head, J/kg.
    pub head: f64,
    /// The measured polytropic efficiency, as a fraction.
    pub efficiency: f64,
}

/// One speed line, kept for the surge and stonewall lookups.
#[derive(Debug, Clone, PartialEq)]
struct SpeedCurve {
    speed: f64,
    flow: Vec<f64>,
    head: Vec<f64>,
}

/// The reduced fits and the speed lines, as `CompressorChart` holds them.
#[derive(Debug, Clone, PartialEq)]
pub struct CompressorChart {
    head: [f64; COEFFICIENTS],
    efficiency: [f64; COEFFICIENTS],
    reference_speed: f64,
    min_speed: f64,
    max_speed: f64,
    curves: Vec<SpeedCurve>,
}

impl CompressorChart {
    /// Build the reduced fits and the speed lines from a flattened map.
    ///
    /// **Three points are the fit's own requirement**, as the pump's chart states: the class
    /// delegates to `PolynomialCurveFitter.create(2)`, which refuses fewer.
    ///
    /// # Errors
    /// [`AzothError::InvalidInput`] for fewer than three points, a non-finite value or a
    /// non-positive speed.
    pub fn new(points: &[CurvePoint]) -> Result<Self> {
        if points.len() < COEFFICIENTS {
            return Err(AzothError::invalid_input(
                "curve",
                format!(
                    "a degree-2 reduced curve needs at least {COEFFICIENTS} measured points, and \
                     this chart has {}",
                    points.len()
                ),
            ));
        }
        for (index, point) in points.iter().enumerate() {
            for (name, value) in [
                ("speed", point.speed),
                ("flow", point.flow),
                ("head", point.head),
                ("efficiency", point.efficiency),
            ] {
                if !value.is_finite() {
                    return Err(AzothError::invalid_input(
                        "curve",
                        format!("point {index}'s {name} is {value}, not a finite number"),
                    ));
                }
            }
            if point.speed <= 0.0 {
                return Err(AzothError::invalid_input(
                    "curve",
                    format!(
                        "point {index}'s speed is {} rpm, and the reduced variables divide by it",
                        point.speed
                    ),
                ));
            }
        }

        let reduced: Vec<f64> = points
            .iter()
            .map(|point| point.flow / point.speed)
            .collect();
        let head = fit(
            &reduced,
            &points
                .iter()
                .map(|point| point.head / (point.speed * point.speed))
                .collect::<Vec<f64>>(),
        )?;
        let efficiency = fit(
            &reduced,
            &points
                .iter()
                .map(|point| point.efficiency)
                .collect::<Vec<f64>>(),
        )?;

        let mut curves: Vec<SpeedCurve> = Vec::new();
        for point in points {
            match curves.last_mut() {
                Some(curve) if curve.speed == point.speed => {
                    curve.flow.push(point.flow);
                    curve.head.push(point.head);
                }
                _ => curves.push(SpeedCurve {
                    speed: point.speed,
                    flow: vec![point.flow],
                    head: vec![point.head],
                }),
            }
        }
        let min_speed = points
            .iter()
            .map(|point| point.speed)
            .fold(f64::INFINITY, f64::min);
        let max_speed = points
            .iter()
            .map(|point| point.speed)
            .fold(f64::NEG_INFINITY, f64::max);

        Ok(Self {
            head,
            efficiency,
            reference_speed: 0.5 * (max_speed + min_speed),
            min_speed,
            max_speed,
            curves,
        })
    }

    /// The midpoint of the measured speed range, rpm.
    pub fn reference_speed(&self) -> f64 {
        self.reference_speed
    }

    /// `getPolytropicHead(flow, speed)`, J/kg. **No clamp**: the class extrapolates.
    pub fn polytropic_head(&self, flow: f64, speed: f64) -> f64 {
        evaluate(&self.head, flow / speed) * speed * speed
    }

    /// `getPolytropicEfficiency(flow, speed)`, as a fraction. **No clamp**.
    pub fn polytropic_efficiency(&self, flow: f64, speed: f64) -> f64 {
        evaluate(&self.efficiency, flow / speed)
    }

    /// `getSpeedValue(flow, head)`: the speed that reaches a head at a flow, rpm.
    ///
    /// The class's fan-law guess, damped Newton (`0.7`, a `30 %` step cap), and its
    /// bisection fallback with the bounds doubled or halved when the target is outside. The
    /// two constants are the class's own and are kept rather than re-tuned.
    pub fn speed_for_head(&self, flow: f64, head: f64) -> f64 {
        let reference_head = self.polytropic_head(flow, self.reference_speed);
        let initial = if reference_head > 0.0 && head > 0.0 {
            self.reference_speed * (head / reference_head).sqrt()
        } else {
            self.reference_speed
        };

        let mut lower = self.min_speed;
        if lower <= 0.0 {
            lower = 100.0;
        }
        let upper = self.max_speed;
        let mut speed = initial.clamp(lower, upper);

        let mut old_speed = speed * 1.01;
        let mut old_error = self.polytropic_head(flow, old_speed) - head;
        for _ in 0..MAX_NEWTON {
            let new_head = self.polytropic_head(flow, speed);
            let error = new_head - head;
            if error.abs() < SPEED_TOLERANCE {
                return speed;
            }
            let mut derivative = (error - old_error) / (speed - old_speed);
            if derivative.abs() < 1.0e-10 {
                derivative = 2.0 * new_head / speed;
                if derivative.abs() < 1.0e-10 {
                    break;
                }
            }
            let update = (DAMPING * error / derivative)
                .clamp(-MAX_STEP_FRACTION * speed, MAX_STEP_FRACTION * speed);
            old_speed = speed;
            old_error = error;
            speed = (speed - update).clamp(lower, upper);
            if (speed - old_speed).abs() < 1.0e-10 {
                break;
            }
        }

        let mut low = lower;
        let mut high = upper;
        let head_at_low = self.polytropic_head(flow, low);
        let head_at_high = self.polytropic_head(flow, high);
        if head < head_at_low && head < head_at_high {
            low *= 0.5;
            high = if head_at_low < head_at_high {
                low * 2.0
            } else {
                high
            };
        } else if head > head_at_low && head > head_at_high {
            high *= 2.0;
        }

        for _ in 0..MAX_BISECTION {
            let mid = 0.5 * (low + high);
            let mid_head = self.polytropic_head(flow, mid);
            if (mid_head - head).abs() < SPEED_TOLERANCE {
                return mid;
            }
            if mid_head < head {
                low = mid;
            } else {
                high = mid;
            }
            if high - low < MIN_SPEED_INTERVAL {
                return mid;
            }
        }
        0.5 * (low + high)
    }

    /// `getSpeed(flow, head)`: the rounded speed, rpm.
    pub fn speed(&self, flow: f64, head: f64) -> i64 {
        self.speed_for_head(flow, head).round() as i64
    }

    /// `getFlow(head, speed, guessFlow)`: the flow at which `head/eta` reaches the target.
    ///
    /// **The residual divides by the efficiency**, which is the class's own counterintuitive
    /// form: `newhead = head(Q)/eta(Q)` against the target. It is the transient path's solver
    /// and this reproduces it rather than the steady branch's head alone.
    pub fn flow_for_head(&self, head: f64, speed: f64, guess: f64) -> f64 {
        let mut flow = guess;
        let mut old_flow = flow * 1.1;
        let mut old_error = self.polytropic_head(old_flow, speed) - head;
        let mut error = 1.0;
        let mut iteration = 1;
        loop {
            iteration += 1;
            let new_head =
                self.polytropic_head(flow, speed) / self.polytropic_efficiency(flow, speed);
            error = new_head - head;
            let derivative = (error - old_error) / (flow - old_flow);
            if derivative.abs() < 1.0e-10 {
                break;
            }
            old_flow = flow;
            old_error = error;
            flow -= error / derivative;
            if flow < 0.0 {
                flow = guess * 0.1;
            }
            if !(error.abs() > FLOW_TOLERANCE && iteration < MAX_FLOW_ITERATIONS) {
                break;
            }
        }
        flow.max(0.0)
    }

    /// `getSurgeFlowAtSpeed`: the minimum measured flow of the closest speed line.
    pub fn surge_flow_at_speed(&self, speed: f64) -> f64 {
        match self.closest(speed) {
            Some(curve) => curve.flow.iter().copied().fold(f64::INFINITY, f64::min),
            None => f64::NAN,
        }
    }

    /// `getStonewallFlowAtSpeed`: the maximum measured flow of the closest speed line.
    pub fn stonewall_flow_at_speed(&self, speed: f64) -> f64 {
        match self.closest(speed) {
            Some(curve) => curve.flow.iter().copied().fold(f64::NEG_INFINITY, f64::max),
            None => f64::NAN,
        }
    }

    /// `getSurgeHeadAtSpeed`: the head at the minimum-flow point of the closest line.
    pub fn surge_head_at_speed(&self, speed: f64) -> f64 {
        match self.closest(speed) {
            Some(curve) => {
                let index = minimum_index(&curve.flow);
                curve.head[index]
            }
            None => f64::NAN,
        }
    }

    /// `getStoneWallHeadAtSpeed`: the head at the maximum-flow point of the closest line.
    pub fn stonewall_head_at_speed(&self, speed: f64) -> f64 {
        match self.closest(speed) {
            Some(curve) => {
                let index = maximum_index(&curve.flow);
                curve.head[index]
            }
            None => f64::NAN,
        }
    }

    /// `getFlowRangeStatus`: `IN_RANGE`, the two extrapolation names, or `NO_CHART`.
    pub fn flow_range_status(&self, flow: f64, speed: f64) -> &'static str {
        if self.curves.is_empty() {
            return "NO_CHART";
        }
        let low = self.surge_flow_at_speed(speed);
        let high = self.stonewall_flow_at_speed(speed);
        if low.is_nan() || high.is_nan() {
            return "NO_CHART";
        }
        let tolerance = 1.0e-6 * high.max(1.0);
        if flow < low - tolerance {
            "EXTRAPOLATED_LOW_FLOW"
        } else if flow > high + tolerance {
            "EXTRAPOLATED_HIGH_FLOW"
        } else {
            "IN_RANGE"
        }
    }

    fn closest(&self, speed: f64) -> Option<&SpeedCurve> {
        let mut best: Option<&SpeedCurve> = None;
        let mut best_difference = f64::INFINITY;
        for curve in &self.curves {
            let difference = (curve.speed - speed).abs();
            if difference < best_difference {
                best_difference = difference;
                best = Some(curve);
            }
        }
        best
    }
}

fn minimum_index(values: &[f64]) -> usize {
    let mut index = 0;
    let mut value = values[0];
    for (candidate, &entry) in values.iter().enumerate().skip(1) {
        if entry < value {
            value = entry;
            index = candidate;
        }
    }
    index
}

fn maximum_index(values: &[f64]) -> usize {
    let mut index = 0;
    let mut value = values[0];
    for (candidate, &entry) in values.iter().enumerate().skip(1) {
        if entry > value {
            value = entry;
            index = candidate;
        }
    }
    index
}

fn fit(reduced_flow: &[f64], values: &[f64]) -> Result<[f64; COEFFICIENTS]> {
    let design: Vec<Vec<f64>> = reduced_flow
        .iter()
        .map(|reduced| vec![1.0, *reduced, reduced * reduced])
        .collect();
    let coefficients = solve_least_squares(&design, values).map_err(|error| {
        AzothError::invalid_input("curve", format!("the reduced fit did not solve: {error}"))
    })?;
    Ok([coefficients[0], coefficients[1], coefficients[2]])
}

fn evaluate(coefficients: &[f64; COEFFICIENTS], x: f64) -> f64 {
    (coefficients[2] * x + coefficients[1]) * x + coefficients[0]
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The eight-speed `CompressorChartTest` table, head in kJ/kg and efficiency in percent.
    fn table() -> Vec<CurvePoint> {
        let speeds = [
            12913.0, 12298.0, 11683.0, 11098.0, 10453.0, 9224.0, 8609.0, 8200.0,
        ];
        let flows = [
            [
                2789.1285, 3174.0375, 3689.2288, 4179.4503, 4570.2768, 4954.7728, 5246.0329,
                5661.0331,
            ],
            [
                2571.1753, 2943.7254, 3440.2675, 3837.4448, 4253.0898, 4668.6643, 4997.1926,
                5387.4952,
            ],
            [
                2415.3793, 2763.0706, 3141.7095, 3594.7436, 4047.6467, 4494.1889, 4853.7353,
                5138.7858,
            ],
            [
                2247.2043, 2799.7342, 3178.3428, 3656.1551, 4102.778, 4394.1591, 4648.3224,
                4840.4998,
            ],
            [
                2072.8397, 2463.9483, 2836.4078, 3202.5266, 3599.6333, 3978.0203, 4257.0022,
                4517.345,
            ],
            [
                1835.9552, 2208.455, 2618.1322, 2940.8034, 3244.7852, 3530.1279, 3753.3738,
                3895.9746,
            ],
            [
                1711.3386, 1965.8848, 2356.9431, 2685.9247, 3008.5154, 3337.2855, 3591.5092, 0.0,
            ],
            [
                1636.5807, 2002.8708, 2338.0319, 2642.1245, 2896.4894, 3113.6264, 3274.8764,
                3411.2977,
            ],
        ];
        let heads = [
            [
                80.0375, 78.8934, 76.2142, 71.8678, 67.0062, 60.6061, 53.0499, 39.728,
            ],
            [
                72.2122, 71.8369, 68.9009, 65.8341, 60.7167, 54.702, 47.2749, 35.7471,
            ],
            [
                65.1576, 64.5253, 62.6118, 59.1619, 54.0455, 47.0059, 39.195, 31.6387,
            ],
            [
                58.6154, 56.9627, 54.6647, 50.4462, 44.4322, 38.4144, 32.9084, 28.8109,
            ],
            [
                52.3295, 51.0573, 49.5283, 46.3326, 42.3685, 37.2502, 31.4884, 25.598,
            ],
            [
                40.6578, 39.6416, 37.6008, 34.6603, 30.9503, 27.1116, 23.2713, 20.4546,
            ],
            [
                35.2705, 34.6359, 32.7228, 31.0645, 27.0985, 22.7482, 18.0113, 0.0,
            ],
            [
                32.192, 31.1756, 29.1329, 26.833, 23.8909, 21.3324, 18.7726, 16.3403,
            ],
        ];
        let efficiencies = [
            [
                77.2452238409573,
                79.4154186459363,
                80.737960012489,
                80.5229826589649,
                79.2210931638144,
                75.4719133864634,
                69.6034181197298,
                58.7322388482707,
            ],
            [
                77.0107837113504,
                79.3069974136389,
                80.8941189021135,
                80.7190194665918,
                79.5313242980328,
                75.5912622896367,
                69.6846136362097,
                60.0043057990909,
            ],
            [
                77.0043065299874,
                79.1690958847856,
                80.8038169975675,
                80.6543975614197,
                78.8532389102705,
                73.6664774270613,
                66.2735600426727,
                57.671664571658,
            ],
            [
                77.0716623789093,
                80.4629750233093,
                81.1390811169072,
                79.6374242667478,
                75.380928428817,
                69.5332969549779,
                63.7997587622339,
                58.8120614497758,
            ],
            [
                76.9705872525642,
                79.8335492585324,
                80.9468133671171,
                80.5806471927835,
                78.0462158225426,
                73.0403707523258,
                66.5572286338589,
                59.8624822515064,
            ],
            [
                77.5063036680357,
                80.2056198362559,
                81.0339108025933,
                79.6085962687939,
                76.3814534404405,
                70.8027503005902,
                64.6437367160571,
                60.5299349982342,
            ],
            [
                77.8175271586685,
                80.065165942218,
                81.0631362122632,
                79.8955051771299,
                76.1983240929369,
                69.289982774309,
                60.8567149372229,
                0.0,
            ],
            [
                78.0924334304045,
                80.9353551568667,
                80.7904437766234,
                78.8639325223295,
                75.2170936751143,
                70.3105081673411,
                65.5507568533569,
                61.0391468300337,
            ],
        ];
        let mut points = Vec::new();
        for (row, speed) in speeds.iter().enumerate() {
            let count = if row == 6 { 7 } else { 8 };
            for column in 0..count {
                points.push(CurvePoint {
                    speed: *speed,
                    flow: flows[row][column],
                    head: heads[row][column] * 1000.0,
                    efficiency: efficiencies[row][column] / 100.0,
                });
            }
        }
        points
    }

    #[test]
    fn the_capture_rows_reproduce_the_class() {
        // `validation/neqsim/captures/process_compressor_curve.tsv`.
        let chart = CompressorChart::new(&table()).expect("a valid chart");
        for (flow, speed, head_kj, efficiency_percent) in [
            (3974.8455, 12913.0, 71.33705849976621, 80.7609351748815),
            (3599.6333, 10453.0, 42.32952036551759, 77.58262465303945),
            (1636.5807, 8200.0, 31.99842270699048, 76.35542784660711),
            (6000.0, 12913.0, 28.035736231331526, 50.487041963431494),
            (100.0, 12913.0, 49.67512659722273, 17.784484294150236),
        ] {
            let head = head_kj * 1000.0;
            let efficiency = efficiency_percent / 100.0;
            assert!(
                (chart.polytropic_head(flow, speed) - head).abs() < 1.0e-9 * head.abs(),
                "head at {flow} m3/hr and {speed} rpm: {} against {head}",
                chart.polytropic_head(flow, speed)
            );
            assert!(
                (chart.polytropic_efficiency(flow, speed) - efficiency).abs() < 1.0e-12,
                "efficiency at {flow} m3/hr and {speed} rpm: {} against {efficiency}",
                chart.polytropic_efficiency(flow, speed)
            );
        }

        let speed = chart.speed_for_head(3974.8455, 70_000.0);
        // The class stops its Newton at `1e-6` on the head, so its reported speed can differ
        // in the last few hundredths of an rpm; the fit itself agrees to `1e-12`.
        assert!(
            (speed - 12822.382970080176).abs() < 1.0e-2,
            "speed: {speed}"
        );
        assert_eq!(chart.speed(3974.8455, 70_000.0), 12822);

        let flow = chart.flow_for_head(70_000.0, 12913.0, 4000.0);
        assert!((flow - 5399.049830646614).abs() < 1.0e-2, "flow: {flow}");
    }

    #[test]
    fn the_surge_and_stonewall_rows_are_the_closest_curve_s() {
        let chart = CompressorChart::new(&table()).expect("a valid chart");
        assert!((chart.surge_flow_at_speed(11683.0) - 2415.3793).abs() < 1.0e-9);
        assert!((chart.stonewall_flow_at_speed(11683.0) - 5138.7858).abs() < 1.0e-9);
        assert!((chart.surge_head_at_speed(11683.0) - 65_157.6).abs() < 1.0e-6);
        assert!((chart.stonewall_head_at_speed(11683.0) - 31_638.7).abs() < 1.0e-6);

        let mid = 0.5 * (2415.3793 + 5138.7858);
        assert_eq!(chart.flow_range_status(mid, 11683.0), "IN_RANGE");
        assert_eq!(
            chart.flow_range_status(2415.3793 - 500.0, 11683.0),
            "EXTRAPOLATED_LOW_FLOW"
        );
        assert_eq!(
            chart.flow_range_status(5138.7858 + 500.0, 11683.0),
            "EXTRAPOLATED_HIGH_FLOW"
        );
    }

    #[test]
    fn a_two_point_chart_is_refused() {
        let points = [
            CurvePoint {
                speed: 1000.0,
                flow: 10.0,
                head: 100_000.0,
                efficiency: 0.7,
            },
            CurvePoint {
                speed: 1000.0,
                flow: 20.0,
                head: 90_000.0,
                efficiency: 0.75,
            },
        ];
        let error = CompressorChart::new(&points).expect_err("two points cannot fit");
        assert_eq!(error.field(), Some("curve"));
    }
}
