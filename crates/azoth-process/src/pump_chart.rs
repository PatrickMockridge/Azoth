//! `PumpChart` - a manufacturer's head and efficiency curve, reduced by the affinity laws.
//!
//! **One quadratic in `Q/N`.** `setCurves` divides every measured flow by its speed and every
//! head by the speed squared, then fits one degree-2 polynomial through all the reduced points;
//! `getHead` evaluates it at the clamped reduced flow and multiplies by the speed squared. The
//! efficiency is fitted and evaluated the same way but **not** divided by the speed: the class
//! stores the measured percentage as its own reduced variable. That is the affinity law the
//! class's tests measure, and it is why a chart at two speeds has one curve, not two.
//!
//! **The reduced flow is clamped to the measured range, and the NPSH flow is not.** `getHead`
//! and `getEfficiency` substitute `minReducedFlow` or `maxReducedFlow` outside the range, so a
//! query cannot leave the fitted interval; `getNPSHRequired` evaluates wherever it is asked and
//! clamps only a negative reduced NPSH to `0.01`.
//!
//! **The head is in metres and the speed in rpm.** `headUnit` has a second value, `kJ/kg`, for
//! the same specific energy divided by `g`; the model's record carries the metre head, and the
//! conversion is a unit conversion rather than a second calculation.
//!
//! **The corrections are the class's own.** `getCorrectedHead` multiplies by
//! `referenceDensity/actualDensity`; `getFullyCorrectedHead` first converts the flow through the
//! Hydraulic Institute's flow factor, evaluates the water curve there and multiplies by the head
//! factor, then applies the density correction. The factors are fitted by
//! [`viscosity_factors`], the class's `calculateViscosityCorrection`.

use azoth_core::{AzothError, Result};
use azoth_reactions::linalg::solve_least_squares;

/// `ThermodynamicConstantsInterface.gravity`, m/s^2. Also `uom`'s `standard_gravity`.
pub const GRAVITY: f64 = 9.80665;
/// The class's `m3/hr -> gpm` factor in the Hydraulic Institute correlation.
pub const GPM_PER_CUBIC_METRE_PER_HOUR: f64 = 4.40287;
/// The class's `m -> ft` factor in the same correlation.
pub const FEET_PER_METRE: f64 = 3.28084;
/// `getEfficiency`'s lower clamp, as a fraction (the class clamps at 0.1 percent).
pub const MIN_EFFICIENCY: f64 = 0.001;
/// `getEfficiency`'s upper clamp, as a fraction (the class clamps at 99 percent).
pub const MAX_EFFICIENCY: f64 = 0.99;
/// `getNPSHRequired`'s floor for a negative reduced NPSH, m per rpm squared.
pub const MIN_REDUCED_NPSH: f64 = 0.01;
/// The class's degree-2 polynomial fit has three coefficients.
const COEFFICIENTS: usize = 3;
/// `getBestEfficiencyFlowRate`'s own scan resolution.
const BEP_INTERVALS: usize = 50;

/// One measured point of a chart: a flow and the head and efficiency there, at one speed.
///
/// The class's `setCurves` takes parallel speed, flow, head and efficiency arrays; the model's
/// record list is the same table flattened, and this is one row of it. `efficiency` is a
/// fraction, not the class's percentage, and `npsh` is the same point's required NPSH where the
/// chart states one.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CurvePoint {
    /// The pump speed, rpm.
    pub speed: f64,
    /// The measured flow, m3/hr.
    pub flow: f64,
    /// The measured head, m.
    pub head: f64,
    /// The measured efficiency, as a fraction.
    pub efficiency: f64,
    /// The required NPSH at this point, m, where the chart states one.
    pub npsh: Option<f64>,
}

/// The Hydraulic Institute's three correction factors for a viscous fluid.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ViscosityFactors {
    /// `Cq`: the flow factor.
    pub flow: f64,
    /// `Ch`: the head factor.
    pub head: f64,
    /// `Ceta`: the efficiency factor.
    pub efficiency: f64,
}

/// The reduced quadratic and the corrections, as `PumpChart` holds them.
#[derive(Debug, Clone, PartialEq)]
pub struct PumpChart {
    head: [f64; COEFFICIENTS],
    efficiency: [f64; COEFFICIENTS],
    npsh: Option<[f64; COEFFICIENTS]>,
    min_reduced_flow: f64,
    max_reduced_flow: f64,
    reference_speed: f64,
    reference_density: Option<f64>,
    use_viscosity_correction: bool,
}

impl PumpChart {
    /// Build the reduced fits from a flattened chart table.
    ///
    /// **Three points are the fit's own requirement.** The class delegates to
    /// `PolynomialCurveFitter.create(2)`, which refuses fewer than three observations, and a
    /// table whose points all share one reduced flow does not define a quadratic; both are
    /// [`AzothError::InvalidInput`] here rather than a fitted line the data does not support.
    ///
    /// **An NPSH curve is all or nothing.** The class's `setNPSHCurve` takes one value per flow
    /// point, so a table that states NPSH at some of its points has no curve; this refuses it
    /// rather than fitting the subset, which would silently drop rows.
    pub fn new(
        points: &[CurvePoint],
        reference_density: Option<f64>,
        use_viscosity_correction: bool,
    ) -> Result<Self> {
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
            if let Some(npsh) = point.npsh {
                if !npsh.is_finite() {
                    return Err(AzothError::invalid_input(
                        "curve",
                        format!("point {index}'s npsh is {npsh}, not a finite number"),
                    ));
                }
            }
        }
        if points.iter().any(|point| point.npsh.is_some())
            && points.iter().any(|point| point.npsh.is_none())
        {
            return Err(AzothError::invalid_input(
                "curve",
                "the chart states a required NPSH at some points and not at others, and the \
                 class's NPSH curve is one value per flow point",
            ));
        }

        let reduced_flow: Vec<f64> = points
            .iter()
            .map(|point| point.flow / point.speed)
            .collect();
        let min_reduced_flow = reduced_flow.iter().copied().fold(f64::INFINITY, f64::min);
        let max_reduced_flow = reduced_flow
            .iter()
            .copied()
            .fold(f64::NEG_INFINITY, f64::max);
        if max_reduced_flow == min_reduced_flow {
            return Err(AzothError::invalid_input(
                "curve",
                format!(
                    "every point has the same reduced flow {min_reduced_flow}, so the reduced \
                     curve has no interval to fit"
                ),
            ));
        }

        let head = fit(
            &reduced_flow,
            &points
                .iter()
                .map(|point| point.head / (point.speed * point.speed))
                .collect::<Vec<f64>>(),
        )?;
        let efficiency = fit(
            &reduced_flow,
            &points
                .iter()
                .map(|point| point.efficiency)
                .collect::<Vec<f64>>(),
        )?;
        let npsh = if points.iter().all(|point| point.npsh.is_some()) {
            Some(fit(
                &reduced_flow,
                &points
                    .iter()
                    .map(|point| point.npsh.unwrap_or(0.0) / (point.speed * point.speed))
                    .collect::<Vec<f64>>(),
            )?)
        } else {
            None
        };

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
            npsh,
            min_reduced_flow,
            max_reduced_flow,
            reference_speed: 0.5 * (max_speed + min_speed),
            reference_density: reference_density.filter(|density| *density > 0.0),
            use_viscosity_correction,
        })
    }

    /// The fitted curve's own midpoint speed, rpm.
    pub fn reference_speed(&self) -> f64 {
        self.reference_speed
    }

    /// The measured reduced-flow interval, m3/hr per rpm.
    pub fn reduced_flow_range(&self) -> (f64, f64) {
        (self.min_reduced_flow, self.max_reduced_flow)
    }

    /// Whether the chart carries an NPSH curve.
    pub fn has_npsh_curve(&self) -> bool {
        self.npsh.is_some()
    }

    /// Whether a reference density makes the density correction active.
    pub fn has_density_correction(&self) -> bool {
        self.reference_density.is_some()
    }

    /// `isUseViscosityCorrection`.
    pub fn uses_viscosity_correction(&self) -> bool {
        self.use_viscosity_correction
    }

    /// `getHead(flow, speed)`, m.
    pub fn head(&self, flow: f64, speed: f64) -> f64 {
        let reduced = self.clamped_reduced_flow(flow, speed);
        evaluate(&self.head, reduced) * speed * speed
    }

    /// `getEfficiency(flow, speed)`, as a fraction.
    pub fn efficiency(&self, flow: f64, speed: f64) -> f64 {
        let reduced = self.clamped_reduced_flow(flow, speed);
        evaluate(&self.efficiency, reduced).clamp(MIN_EFFICIENCY, MAX_EFFICIENCY)
    }

    /// `getNPSHRequired(flow, speed)`, m; zero where the chart has no NPSH curve.
    pub fn npsh_required(&self, flow: f64, speed: f64) -> f64 {
        let Some(coefficients) = self.npsh else {
            return 0.0;
        };
        let reduced = evaluate(&coefficients, flow / speed);
        // **Only a negative reduced NPSH is floored**, which is the class's own branch: it
        // substitutes `0.01` for a negative value and leaves a small positive one alone.
        let reduced = if reduced < 0.0 {
            MIN_REDUCED_NPSH
        } else {
            reduced
        };
        reduced * speed * speed
    }

    /// `getBestEfficiencyFlowRate(speed)`, m3/hr.
    pub fn best_efficiency_flow_rate(&self, speed: f64) -> f64 {
        if speed <= 0.0 {
            return 0.0;
        }
        let mut best_flow = 0.0;
        let mut best_efficiency = 0.0;
        for index in 0..=BEP_INTERVALS {
            let reduced = self.min_reduced_flow
                + (self.max_reduced_flow - self.min_reduced_flow) * index as f64
                    / BEP_INTERVALS as f64;
            let efficiency = evaluate(&self.efficiency, reduced);
            if efficiency > best_efficiency {
                best_efficiency = efficiency;
                best_flow = reduced * speed;
            }
        }
        best_flow
    }

    /// `getCorrectedHead(flow, speed, actualDensity)`, m.
    pub fn corrected_head(&self, flow: f64, speed: f64, actual_density: f64) -> f64 {
        let head = self.head(flow, speed);
        match self.reference_density {
            Some(reference) if actual_density > 0.0 => head * reference / actual_density,
            _ => head,
        }
    }

    /// `getFullyCorrectedHead(flow, speed, actualDensity, actualViscosity)`, m.
    pub fn fully_corrected_head(
        &self,
        flow: f64,
        speed: f64,
        actual_density: f64,
        actual_viscosity: f64,
    ) -> f64 {
        if !(self.use_viscosity_correction && actual_viscosity > 1.0) {
            return self.corrected_head(flow, speed, actual_density);
        }
        let factors = self.factors_at(actual_viscosity, speed);
        let water_head = self.head(flow / factors.flow, speed) * factors.head;
        match self.reference_density {
            Some(reference) if actual_density > 0.0 => water_head * reference / actual_density,
            _ => water_head,
        }
    }

    /// `getCorrectedEfficiency(flow, speed, actualViscosity)`, as a fraction.
    pub fn corrected_efficiency(&self, flow: f64, speed: f64, actual_viscosity: f64) -> f64 {
        let base = self.efficiency(flow, speed);
        if self.use_viscosity_correction && actual_viscosity > 1.0 {
            return base * self.factors_at(actual_viscosity, speed).efficiency;
        }
        base
    }

    /// The class's `calculateViscosityCorrection` at this chart's best-efficiency point.
    fn factors_at(&self, viscosity: f64, speed: f64) -> ViscosityFactors {
        let flow_bep = self.best_efficiency_flow_rate(self.reference_speed);
        let head_bep = self.head(flow_bep, self.reference_speed);
        viscosity_factors(viscosity, flow_bep, head_bep, speed)
    }

    fn clamped_reduced_flow(&self, flow: f64, speed: f64) -> f64 {
        (flow / speed).clamp(self.min_reduced_flow, self.max_reduced_flow)
    }
}

/// `calculateViscosityCorrection(viscosity, flowBEP, headBEP, speed)`: the HI factors.
///
/// The correlation is the class's, in its own units: `Q` in gpm, `H` in feet, `N` in rpm. The
/// piecewise answer is the class's too - unity below `B = 1`, the empirical factors below
/// `B = 40`, and the floor triple above it.
pub fn viscosity_factors(
    viscosity: f64,
    flow_bep: f64,
    head_bep: f64,
    speed: f64,
) -> ViscosityFactors {
    let unity = ViscosityFactors {
        flow: 1.0,
        head: 1.0,
        efficiency: 1.0,
    };
    if viscosity <= 1.0 {
        return unity;
    }
    let flow_gpm = flow_bep * GPM_PER_CUBIC_METRE_PER_HOUR;
    let head_ft = head_bep * FEET_PER_METRE;
    if flow_gpm.is_nan()
        || flow_gpm <= 0.0
        || head_ft.is_nan()
        || head_ft <= 0.0
        || speed.is_nan()
        || speed <= 0.0
    {
        return unity;
    }
    let b =
        26.6 * viscosity.sqrt() * head_ft.powf(0.0625) / (flow_gpm.powf(0.375) * speed.powf(0.25));
    if b <= 1.0 {
        return unity;
    }
    if b <= 40.0 {
        return ViscosityFactors {
            flow: (1.0 - 0.01 * b.powf(0.9)).max(0.6),
            head: (1.0 - 0.008 * b).max(0.7),
            efficiency: (1.0 - 0.015 * b.powf(0.85)).max(0.4),
        };
    }
    ViscosityFactors {
        flow: 0.6,
        head: 0.7,
        efficiency: 0.4,
    }
}

/// One degree-2 least-squares fit over the reduced points.
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

/// Horner's rule over the fitted coefficients, ascending.
fn evaluate(coefficients: &[f64; COEFFICIENTS], x: f64) -> f64 {
    (coefficients[2] * x + coefficients[1]) * x + coefficients[0]
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The eight-speed table `PumpChartTest` states, at the reduced points.
    fn eight_speed_table() -> Vec<CurvePoint> {
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
                    head: heads[row][column],
                    efficiency: efficiencies[row][column] / 100.0,
                    npsh: None,
                });
            }
        }
        points
    }

    fn simple_chart() -> PumpChart {
        let flows = [10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0];
        let heads = [120.0, 118.0, 115.0, 110.0, 103.0, 94.0, 83.0, 70.0];
        let efficiencies = [0.60, 0.70, 0.78, 0.82, 0.81, 0.76, 0.68, 0.55];
        let points: Vec<CurvePoint> = (0..8)
            .map(|index| CurvePoint {
                speed: 1000.0,
                flow: flows[index],
                head: heads[index],
                efficiency: efficiencies[index],
                npsh: None,
            })
            .collect();
        PumpChart::new(&points, Some(998.0), false).expect("a valid chart")
    }

    #[test]
    fn the_capture_rows_reproduce_the_class() {
        // `validation/neqsim/captures/process_pump_curve.tsv`.
        let chart = PumpChart::new(&eight_speed_table(), None, false).expect("a valid chart");

        for (flow, speed, head, efficiency) in [
            (3974.8455, 12913.0, 71.33705849976621, 0.807609351748815),
            (3599.6333, 10453.0, 42.32952036551759, 0.7758262465303945),
            (1636.5807, 8200.0, 31.99842270699048, 0.763554278466071),
            (6000.0, 12913.0, 37.37780855561962, 0.5816069550481612),
            (100.0, 12913.0, 79.33940792723672, 0.7617942861177876),
        ] {
            assert!(
                (chart.head(flow, speed) - head).abs() < 1.0e-9 * head.abs(),
                "head at {flow} m3/hr and {speed} rpm: {} against {head}",
                chart.head(flow, speed)
            );
            assert!(
                (chart.efficiency(flow, speed) - efficiency).abs() < 1.0e-12,
                "efficiency at {flow} m3/hr and {speed} rpm: {} against {efficiency}",
                chart.efficiency(flow, speed)
            );
        }

        assert_eq!(chart.reference_speed(), 10556.5);
        assert!((chart.best_efficiency_flow_rate(10556.5) - 2909.341180877987).abs() < 1.0e-6);
        assert!((chart.best_efficiency_flow_rate(1500.0) - 413.39570608790604).abs() < 1.0e-6);
    }

    #[test]
    fn the_density_correction_is_the_reference_over_the_actual() {
        let chart = simple_chart();
        let base = chart.head(40.0, 1000.0);
        assert!((base - 109.95833333333331).abs() < 1.0e-9);
        assert!(chart.has_density_correction());
        assert!((chart.corrected_head(40.0, 1000.0, 998.0) - base).abs() < 1.0e-9);
        assert!((chart.corrected_head(40.0, 1000.0, 800.0) - 137.1730208333333).abs() < 1.0e-9);
        assert!((chart.corrected_head(40.0, 1000.0, 1200.0) - 91.44868055555554).abs() < 1.0e-9);
    }

    #[test]
    fn the_npsh_curve_reproduces_the_capture() {
        let points: Vec<CurvePoint> = {
            let speeds = [1000.0, 1500.0];
            let flows = [
                [10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0],
                [15.0, 30.0, 45.0, 60.0, 75.0, 90.0, 105.0, 120.0],
            ];
            let heads = [
                [120.0, 118.0, 115.0, 110.0, 103.0, 94.0, 83.0, 70.0],
                [270.0, 265.5, 258.8, 247.5, 231.8, 211.5, 186.8, 157.5],
            ];
            let efficiencies = [
                [0.60, 0.70, 0.78, 0.82, 0.81, 0.76, 0.68, 0.55],
                [0.62, 0.71, 0.79, 0.83, 0.82, 0.77, 0.69, 0.56],
            ];
            let npsh = [
                [2.0, 2.2, 2.5, 3.0, 3.8, 4.8, 6.2, 8.0],
                [4.5, 4.95, 5.625, 6.75, 8.55, 10.8, 13.95, 18.0],
            ];
            let mut points = Vec::new();
            for row in 0..2 {
                for column in 0..8 {
                    points.push(CurvePoint {
                        speed: speeds[row],
                        flow: flows[row][column],
                        head: heads[row][column],
                        efficiency: efficiencies[row][column],
                        npsh: Some(npsh[row][column]),
                    });
                }
            }
            points
        };
        let chart = PumpChart::new(&points, None, false).expect("a valid chart");
        assert!(chart.has_npsh_curve());
        for (flow, speed, expected) in [
            (40.0, 1000.0, 2.97797619047619),
            (60.0, 1500.0, 6.700446428571427),
            (35.0, 1000.0, 2.6635416666666663),
            (150.0, 1000.0, 26.837499999999988),
            (1.0, 1000.0, 2.301279761904759),
        ] {
            assert!(
                (chart.npsh_required(flow, speed) - expected).abs() < 1.0e-9 * expected.abs(),
                "npsh at {flow} m3/hr and {speed} rpm: {} against {expected}",
                chart.npsh_required(flow, speed)
            );
        }
    }

    #[test]
    fn the_viscosity_factors_reproduce_the_capture() {
        for (viscosity, flow, head, efficiency) in [
            (
                50.0,
                0.9376922588128928,
                0.9389178325111737,
                0.9155709225691215,
            ),
            (
                75.0,
                0.9252204110055453,
                0.9251899286345778,
                0.899693095354078,
            ),
            (
                100.0,
                0.9148851172309777,
                0.9136167703181568,
                0.8866478422061969,
            ),
            (
                500.0,
                0.8243929257526138,
                0.8068412263154211,
                0.7753574379039675,
            ),
        ] {
            let factors = viscosity_factors(viscosity, 30.0, 88.0, 1000.0);
            assert!((factors.flow - flow).abs() < 1.0e-12);
            assert!((factors.head - head).abs() < 1.0e-12);
            assert!((factors.efficiency - efficiency).abs() < 1.0e-12);
        }
        let water = viscosity_factors(1.0, 30.0, 88.0, 1000.0);
        assert_eq!(water.flow, 1.0);
        assert_eq!(water.head, 1.0);
        assert_eq!(water.efficiency, 1.0);
    }

    #[test]
    fn the_fully_corrected_head_reproduces_the_capture() {
        let mut chart = simple_chart();
        chart.use_viscosity_correction = true;
        // The chart's own BEP at its reference speed, as `getFullyCorrectedHead` recomputes it.
        assert!(
            (chart.fully_corrected_head(30.0, 1000.0, 850.0, 100.0) - 123.70507324402753).abs()
                < 1.0e-9
        );
        assert!(
            (chart.corrected_efficiency(30.0, 1000.0, 100.0) - 0.7019571848170408).abs() < 1.0e-12
        );
        assert!(
            (chart.fully_corrected_head(30.0, 1000.0, 850.0, 1.0) - 135.17029411764702).abs()
                < 1.0e-9
        );
    }

    #[test]
    fn a_chart_of_two_points_is_refused() {
        let points = [
            CurvePoint {
                speed: 1000.0,
                flow: 10.0,
                head: 100.0,
                efficiency: 0.7,
                npsh: None,
            },
            CurvePoint {
                speed: 1000.0,
                flow: 20.0,
                head: 90.0,
                efficiency: 0.75,
                npsh: None,
            },
        ];
        let error = PumpChart::new(&points, None, false).expect_err("two points cannot fit");
        assert_eq!(error.field(), Some("curve"));
    }

    #[test]
    fn a_partial_npsh_curve_is_refused() {
        let points = simple_chart_points();
        let mut with_partial = points.clone();
        with_partial[0].npsh = Some(2.0);
        let error = PumpChart::new(&with_partial, None, false).expect_err("a partial NPSH curve");
        assert_eq!(error.field(), Some("curve"));
    }

    fn simple_chart_points() -> Vec<CurvePoint> {
        let flows = [10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0];
        let heads = [120.0, 118.0, 115.0, 110.0, 103.0, 94.0, 83.0, 70.0];
        let efficiencies = [0.60, 0.70, 0.78, 0.82, 0.81, 0.76, 0.68, 0.55];
        (0..8)
            .map(|index| CurvePoint {
                speed: 1000.0,
                flow: flows[index],
                head: heads[index],
                efficiency: efficiencies[index],
                npsh: None,
            })
            .collect()
    }
}
