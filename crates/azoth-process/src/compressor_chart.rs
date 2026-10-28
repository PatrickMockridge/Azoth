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
//! **The boundary curves are carried, and `setCurves` still builds none.** `getSurgeFlowAtSpeed`
//! and its siblings take the measured curve closest to the requested speed and read its
//! minimum or maximum flow; the spline boundaries (`SafeSplineSurgeCurve` /
//! `SafeSplineStoneWallCurve`) are set explicitly or built by `generateSurgeCurve` /
//! `generateStoneWallCurve` from those same extreme points, which is what the class does - its
//! `setCurves` does not call the generators either.
//!
//! **The dead fan-law correction is not repeated.** `setCurves` fits a third polynomial over
//! `speed/speed[0]` that no getter reads.

use azoth_core::{AzothError, Result};
use azoth_reactions::linalg::solve_least_squares;

use crate::surge_curve::{Spline, StoneWallCurve, SurgeCurve};

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
/// `head` is in the chart's head unit - the model's record carries **metres**, the class's
/// default `headUnit`, so a curve quoted in kJ/kg is multiplied by `1000/9.81` (the class's
/// own `polytropicHeadMeter`). `efficiency` is a fraction; the class's percent column is
/// divided by `100`. Both conversions are linear, so the fit is the same polynomial either
/// way. `speed` is rpm and `flow` m3/hr, the class's own coordinates.
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

/// One speed line, kept for the surge and stonewall lookups and the alternative maps.
#[derive(Debug, Clone, PartialEq)]
struct SpeedCurve {
    speed: f64,
    flow: Vec<f64>,
    head: Vec<f64>,
    efficiency: Vec<f64>,
}

/// Which map the chart interpolates.
///
/// `FanLaw` is `CompressorChart`, the reduced quadratics. `MapLookup` and
/// `MapLookupExtrapolate` are `CompressorChartAlternativeMapLookup` and its extrapolating
/// subclass: the measured speed lines read directly, with a natural cubic spline along each
/// line.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum ChartKind {
    /// `CompressorChart`: one reduced quadratic in `Q/N` per table.
    #[default]
    FanLaw,
    /// `CompressorChartAlternativeMapLookup`: the mean of the nearest lines' splines.
    MapLookup,
    /// `CompressorChartAlternativeMapLookupExtrapolate`: the same, with linear extrapolation
    /// past the flow knots and a linear speed blend.
    MapLookupExtrapolate,
    /// `CompressorChartKhader2015`: the extrapolating lookup in the map's similarity
    /// coordinates - machine Mach number, corrected flow and head over the sound speed.
    Khader2015,
}

/// `CompressorChartKhader2015`'s map: the speed lines in similarity coordinates.
#[derive(Debug, Clone, PartialEq)]
struct KhaderMap {
    /// The reference fluid's phase sound speed, which fits the map, m/s.
    reference_sound_speed: f64,
    /// The process fluid's system sound speed, which reads it, m/s.
    process_sound_speed: f64,
    /// The impeller's outer diameter, m.
    impeller_diameter: f64,
    /// One line per speed, as machine Mach numbers.
    lines: Vec<KhaderLine>,
}

/// One speed line in the similarity coordinates.
#[derive(Debug, Clone, PartialEq)]
struct KhaderLine {
    /// `N/60 D/cs`, the machine Mach number.
    mach: f64,
    /// `Q/3600/cs/D^2`, the corrected flow **both getters interpolate against**. The class
    /// stores a second column without the `3600` and never reads it on this path.
    head_flow: Vec<f64>,
    /// `h/cs^2`.
    head: Vec<f64>,
    /// The measured efficiency, as a fraction.
    efficiency: Vec<f64>,
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
    khader: Option<KhaderMap>,
    kind: ChartKind,
    surge: Option<SurgeCurve>,
    stonewall: Option<StoneWallCurve>,
}

impl CompressorChart {
    /// Build the reduced fits and the speed lines from a flattened map, the fan-law chart.
    ///
    /// **Three points are the fit's own requirement**, as the pump's chart states: the class
    /// delegates to `PolynomialCurveFitter.create(2)`, which refuses fewer.
    ///
    /// # Errors
    /// [`AzothError::InvalidInput`] for fewer than three points, a non-finite value or a
    /// non-positive speed.
    pub fn new(points: &[CurvePoint]) -> Result<Self> {
        Self::with_kind(points, ChartKind::FanLaw)
    }

    /// Build the same speed lines under one of the alternative map lookups.
    ///
    /// The reduced fits are still built - they cost a few flops and keep one constructor - and
    /// the two lookup kinds interpolate the lines instead of reading them.
    ///
    /// # Errors
    /// The same refusals as [`Self::new`].
    pub fn with_kind(points: &[CurvePoint], kind: ChartKind) -> Result<Self> {
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
                    curve.efficiency.push(point.efficiency);
                }
                _ => curves.push(SpeedCurve {
                    speed: point.speed,
                    flow: vec![point.flow],
                    head: vec![point.head],
                    efficiency: vec![point.efficiency],
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
            khader: None,
            kind,
            surge: None,
            stonewall: None,
        })
    }

    /// Build the Khader 2015 chart: the extrapolating lookup in similarity coordinates.
    ///
    /// **Two sound speeds, because the class uses two.** The map is fitted in the *reference*
    /// fluid's phase sound speed - the default methane/ethane/propane mix at the chart
    /// conditions - and read in the *process* fluid's system sound speed. The head's corrected
    /// flow divides by `3600` and the efficiency's does not; both asymmetries are the class's.
    ///
    /// # Errors
    /// The same refusals as [`Self::new`], plus a non-positive sound speed or diameter.
    pub fn with_khader(
        points: &[CurvePoint],
        reference_sound_speed: f64,
        process_sound_speed: f64,
        impeller_diameter: f64,
    ) -> Result<Self> {
        for (name, value) in [
            ("reference_sound_speed", reference_sound_speed),
            ("process_sound_speed", process_sound_speed),
            ("impeller_diameter", impeller_diameter),
        ] {
            if !value.is_finite() || value <= 0.0 {
                return Err(AzothError::invalid_input(
                    name,
                    format!("the Khader map divides by {name}, and {value} is not positive"),
                ));
            }
        }
        let mut chart = Self::with_kind(points, ChartKind::Khader2015)?;
        let area = impeller_diameter * impeller_diameter;
        let mut lines: Vec<KhaderLine> = Vec::new();
        for curve in &chart.curves {
            let mach = curve.speed / 60.0 * impeller_diameter / reference_sound_speed;
            let head_flow = curve
                .flow
                .iter()
                .map(|flow| flow / 3600.0 / reference_sound_speed / area)
                .collect();
            let head = curve
                .head
                .iter()
                .map(|head| head / (reference_sound_speed * reference_sound_speed))
                .collect();
            lines.push(KhaderLine {
                mach,
                head_flow,
                head,
                efficiency: curve.efficiency.clone(),
            });
        }
        chart.khader = Some(KhaderMap {
            reference_sound_speed,
            process_sound_speed,
            impeller_diameter,
            lines,
        });
        Ok(chart)
    }

    /// The interpolation scheme this chart uses.
    #[must_use]
    pub fn kind(&self) -> ChartKind {
        self.kind
    }

    /// The Khader map, where one is built.
    #[must_use]
    pub fn khader(&self) -> bool {
        self.khader.is_some()
    }

    /// The midpoint of the measured speed range, rpm.
    pub fn reference_speed(&self) -> f64 {
        self.reference_speed
    }

    /// `getPolytropicHead(flow, speed)`, J/kg. **No clamp**: the fan law extrapolates its
    /// quadratic, and the lookups refuse a flow outside their spline unless they extrapolate.
    ///
    /// # Errors
    /// [`AzothError::InvalidInput`] where a non-extrapolating lookup's flow lies outside a
    /// speed line's knots - the class's `SplineInterpolator` throws there and it does not catch
    /// it.
    pub fn polytropic_head(&self, flow: f64, speed: f64) -> Result<f64> {
        match self.kind {
            ChartKind::FanLaw => Ok(self.fan_law_head(flow, speed)),
            ChartKind::Khader2015 => {
                let map = self.khader_map()?;
                let (corrected_flow, mach) = map.coordinates(flow, speed);
                Ok(
                    lookup_lines(&map.map_lines(), corrected_flow, mach, true, true, false)?
                        * map.process_sound_speed
                        * map.process_sound_speed,
                )
            }
            _ => self.lookup(flow, speed, true),
        }
    }

    /// The fan-law head, which cannot fail: the reduced quadratic read back at `Q/N`.
    fn fan_law_head(&self, flow: f64, speed: f64) -> f64 {
        evaluate(&self.head, flow / speed) * speed * speed
    }

    /// `getPolytropicEfficiency(flow, speed)`, as a fraction. **No clamp**.
    ///
    /// # Errors
    /// The same refusal as [`Self::polytropic_head`].
    pub fn polytropic_efficiency(&self, flow: f64, speed: f64) -> Result<f64> {
        if self.kind == ChartKind::FanLaw {
            return Ok(evaluate(&self.efficiency, flow / speed));
        }
        if self.kind == ChartKind::Khader2015 {
            let map = self.khader_map()?;
            let (corrected_flow, mach) = map.coordinates(flow, speed);
            return lookup_lines(&map.map_lines(), corrected_flow, mach, false, true, false);
        }
        self.lookup(flow, speed, false)
    }

    /// The alternative lookups: the nearest speed line(s), splined along `flow`.
    fn lookup(&self, flow: f64, speed: f64, head: bool) -> Result<f64> {
        let lines: Vec<MapLine<'_>> = self
            .curves
            .iter()
            .map(|curve| MapLine {
                speed: curve.speed,
                head_flow: &curve.flow,
                efficiency_flow: &curve.flow,
                head: &curve.head,
                efficiency: &curve.efficiency,
            })
            .collect();
        lookup_lines(
            &lines,
            flow,
            speed,
            head,
            self.kind == ChartKind::MapLookupExtrapolate,
            matches!(self.kind, ChartKind::MapLookup),
        )
    }

    fn khader_map(&self) -> Result<&KhaderMap> {
        self.khader.as_ref().ok_or_else(|| {
            AzothError::invalid_input("chart_type", "the Khader kind has no map built")
        })
    }

    /// `getSpeedValue(flow, head)`: the speed that reaches a head at a flow, rpm.
    ///
    /// **The fan-law inverse.** The alternative classes override `getSpeedValue` with their
    /// own Newton; that is the speed-solve tranche's, and this reads the fan-law quadratic.
    ///
    /// The class's fan-law guess, damped Newton (`0.7`, a `30 %` step cap), and its
    /// bisection fallback with the bounds doubled or halved when the target is outside. The
    /// two constants are the class's own and are kept rather than re-tuned.
    pub fn speed_for_head(&self, flow: f64, head: f64) -> f64 {
        let reference_head = self.fan_law_head(flow, self.reference_speed);
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
        let mut old_error = self.fan_law_head(flow, old_speed) - head;
        for _ in 0..MAX_NEWTON {
            let new_head = self.fan_law_head(flow, speed);
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
            // numerics-ok: `CompressorChart.getSpeedValue` caps its damped Newton step
            // at 30 % of the speed, the class's own bound.
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
        let head_at_low = self.fan_law_head(flow, low);
        let head_at_high = self.fan_law_head(flow, high);
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
            let mid_head = self.fan_law_head(flow, mid);
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
        let mut old_error = self.fan_law_head(old_flow, speed) - head;
        let mut error;
        let mut iteration = 1;
        loop {
            iteration += 1;
            let new_head =
                self.fan_law_head(flow, speed) / evaluate(&self.efficiency, flow / speed);
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

    /// `setSurgeCurve`: the explicit surge boundary, in the chart's head unit.
    pub fn set_surge_curve(&mut self, curve: SurgeCurve) {
        self.surge = Some(curve);
    }

    /// `setStoneWallCurve`: the explicit stonewall boundary, in the chart's head unit.
    pub fn set_stonewall_curve(&mut self, curve: StoneWallCurve) {
        self.stonewall = Some(curve);
    }

    /// The surge curve, where one is set.
    #[must_use]
    pub fn surge_curve(&self) -> Option<&SurgeCurve> {
        self.surge.as_ref()
    }

    /// The stonewall curve, where one is set.
    #[must_use]
    pub fn stonewall_curve(&self) -> Option<&StoneWallCurve> {
        self.stonewall.as_ref()
    }

    /// `getSurgeCurve().getSurgeFlow(head)`: the curve's flow at a head, or `None` without one.
    #[must_use]
    pub fn surge_flow(&self, head: f64) -> Option<f64> {
        self.surge.as_ref().map(|curve| curve.flow_at_head(head))
    }

    /// `getStoneWallCurve().getStoneWallFlow(head)`, or `None` without a curve.
    #[must_use]
    pub fn stonewall_flow(&self, head: f64) -> Option<f64> {
        self.stonewall
            .as_ref()
            .map(|curve| curve.flow_at_head(head))
    }

    /// `getSurgeCurve().getSurgeHead(flow)`, or `None` without a curve.
    #[must_use]
    pub fn surge_head(&self, flow: f64) -> Option<f64> {
        self.surge.as_ref().map(|curve| curve.head_at_flow(flow))
    }

    /// `getStoneWallCurve().getStoneWallHead(flow)`, or `None` without a curve.
    #[must_use]
    pub fn stonewall_head(&self, flow: f64) -> Option<f64> {
        self.stonewall
            .as_ref()
            .map(|curve| curve.head_at_flow(flow))
    }

    /// `isSurge(head, flow)`: false where no surge curve is set, the class's inactive answer.
    #[must_use]
    pub fn is_surge(&self, head: f64, flow: f64) -> bool {
        self.surge
            .as_ref()
            .is_some_and(|curve| curve.is_limit(head, flow))
    }

    /// `isStoneWall(head, flow)`: false where no stonewall curve is set.
    #[must_use]
    pub fn is_stonewall(&self, head: f64, flow: f64) -> bool {
        self.stonewall
            .as_ref()
            .is_some_and(|curve| curve.is_limit(head, flow))
    }

    /// `generateSurgeCurve`: one point per speed row - its lowest flow - deduped by flow.
    ///
    /// # Errors
    /// [`AzothError::InvalidInput`] if the chart carries no speed rows, or if the points do not
    /// make a curve (two distinct rows, which the class's spline refuses).
    pub fn generate_surge_curve(&mut self) -> Result<()> {
        let (flow, head) = self.extreme_points(true)?;
        self.surge = Some(SurgeCurve::new(&flow, &head)?);
        Ok(())
    }

    /// `generateStoneWallCurve`: one point per speed row - its highest flow - deduped by flow.
    ///
    /// # Errors
    /// The same refusals as [`Self::generate_surge_curve`].
    pub fn generate_stonewall_curve(&mut self) -> Result<()> {
        let (flow, head) = self.extreme_points(false)?;
        self.stonewall = Some(StoneWallCurve::new(&flow, &head)?);
        Ok(())
    }

    /// The lowest-flow (surge) or highest-flow (stonewall) point of every speed row, in row
    /// order and deduped by flow - the class's `TreeMap` with `containsKey`'s first wins.
    fn extreme_points(&self, lowest: bool) -> Result<(Vec<f64>, Vec<f64>)> {
        if self.curves.is_empty() {
            return Err(AzothError::invalid_input(
                "curve",
                "the chart has no speed rows, so it has no boundary to generate",
            ));
        }
        let mut seen: Vec<f64> = Vec::new();
        let mut flow = Vec::new();
        let mut head = Vec::new();
        for curve in &self.curves {
            let index = if lowest {
                minimum_index(&curve.flow)
            } else {
                maximum_index(&curve.flow)
            };
            let value = curve.flow[index];
            if seen.contains(&value) {
                continue;
            }
            seen.push(value);
            flow.push(value);
            head.push(curve.head[index]);
        }
        Ok((flow, head))
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

/// One speed line, borrowed for a lookup.
struct MapLine<'a> {
    speed: f64,
    head_flow: &'a [f64],
    efficiency_flow: &'a [f64],
    head: &'a [f64],
    efficiency: &'a [f64],
}

impl KhaderMap {
    /// The corrected flow and machine Mach number a read is taken at.
    fn coordinates(&self, flow: f64, speed: f64) -> (f64, f64) {
        let area = self.impeller_diameter * self.impeller_diameter;
        let corrected = flow / 3600.0 / self.process_sound_speed / area;
        let mach = speed / 60.0 * self.impeller_diameter / self.process_sound_speed;
        (corrected, mach)
    }

    /// The map's lines, with the head's flow axis or the efficiency's.
    fn map_lines(&self) -> Vec<MapLine<'_>> {
        self.lines
            .iter()
            .map(|line| {
                // **Both getters read `curve.flow`**, the head's corrected flow; the class's
                // `flowPolytropicEfficiency` column is stored and never read on this path.
                MapLine {
                    speed: line.mach,
                    head_flow: &line.head_flow,
                    efficiency_flow: &line.head_flow,
                    head: &line.head,
                    efficiency: &line.efficiency,
                }
            })
            .collect()
    }
}

/// The mean-over-lines rule (`MapLookup`) or the extrapolating blend, over any set of lines.
fn lookup_lines(
    lines: &[MapLine<'_>],
    flow: f64,
    speed: f64,
    head: bool,
    extrapolate: bool,
    mean: bool,
) -> Result<f64> {
    if lines.is_empty() {
        return Err(AzothError::invalid_input(
            "curve",
            "the chart has no speed lines, so a map lookup has nothing to read",
        ));
    }
    let mut order: Vec<usize> = (0..lines.len()).collect();
    order.sort_by(|a, b| lines[*a].speed.total_cmp(&lines[*b].speed));
    let mut index = 0;
    while index < order.len() && lines[order[index]].speed < speed {
        index += 1;
    }
    let closest: Vec<usize> = if index < order.len() && lines[order[index]].speed == speed {
        vec![order[index]]
    } else if index == 0 {
        vec![order[0]]
    } else if index == order.len() {
        vec![order[order.len() - 1]]
    } else {
        vec![order[index - 1], order[index]]
    };

    let mut values: Vec<(f64, f64)> = Vec::with_capacity(closest.len());
    for position in closest {
        let line = &lines[position];
        let (knots, ordinates) = if head {
            (line.head_flow, line.head)
        } else {
            (line.efficiency_flow, line.efficiency)
        };
        values.push((
            line.speed,
            line_value(knots, ordinates, flow, extrapolate, line.speed)?,
        ));
    }
    if mean {
        return Ok(values.iter().map(|(_, value)| value).sum::<f64>() / values.len() as f64);
    }
    if values.len() == 1 {
        let (reference, value) = values[0];
        // The extrapolating class scales a single line's head by the speed ratio (not its
        // square), and leaves the efficiency alone.
        return Ok(if head {
            value * speed / reference
        } else {
            value
        });
    }
    let (first, first_value) = values[0];
    let (second, second_value) = values[1];
    let slope = (second_value - first_value) / (second - first);
    let blended = first_value + slope * (speed - first);
    // **The pinned class scales the blended head by `speed / max(s1, s2)`**, and leaves the
    // efficiency unscaled; the later working-tree source drops it and the jar is the oracle.
    Ok(if head {
        blended * speed / first.max(second)
    } else {
        blended
    })
}

/// One speed line's head or efficiency at a flow: the natural spline inside the knots, and
/// either a refusal or a linear extrapolation outside them.
fn line_value(
    knots: &[f64],
    ordinates: &[f64],
    flow: f64,
    extrapolate: bool,
    speed: f64,
) -> Result<f64> {
    let n = knots.len();
    if n < 2 || knots.windows(2).any(|window| window[1] <= window[0]) {
        return Err(AzothError::invalid_input(
            "curve",
            format!("the speed line at {speed} rpm is not strictly increasing in flow"),
        ));
    }
    if flow < knots[0] || flow > knots[n - 1] {
        if !extrapolate {
            // The class's `PolynomialSplineFunction` throws outside its knots and the lookup
            // does not catch it; the refusal is that throw.
            return Err(AzothError::invalid_input(
                "curve",
                format!(
                    "the map lookup's spline at {speed} rpm covers {} to {} m3/hr, and {flow} \
                     is outside it; only the extrapolating chart reads past the knots",
                    knots[0],
                    knots[n - 1]
                ),
            ));
        }
        let (low, high) = if flow < knots[0] {
            (0, 1)
        } else {
            (n - 2, n - 1)
        };
        let slope = (ordinates[high] - ordinates[low]) / (knots[high] - knots[low]);
        return Ok(ordinates[low] + slope * (flow - knots[low]));
    }
    Ok(Spline::new(knots.to_vec(), ordinates.to_vec()).at(flow))
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
                    head: heads[row][column] * KJ_PER_KG_TO_M,
                    efficiency: efficiencies[row][column] / 100.0,
                });
            }
        }
        points
    }

    /// The class's `polytropicHeadMeter`: kJ/kg to metres by `1000/9.81`.
    const KJ_PER_KG_TO_M: f64 = 1000.0 / 9.81;

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
            let head = head_kj * KJ_PER_KG_TO_M;
            let efficiency = efficiency_percent / 100.0;
            assert!(
                (chart
                    .polytropic_head(flow, speed)
                    .expect("fan law never refuses")
                    - head)
                    .abs()
                    < 1.0e-9 * head.abs(),
                "head at {flow} m3/hr and {speed} rpm: {} against {head}",
                chart
                    .polytropic_head(flow, speed)
                    .expect("fan law never refuses")
            );
            assert!(
                (chart
                    .polytropic_efficiency(flow, speed)
                    .expect("fan law never refuses")
                    - efficiency)
                    .abs()
                    < 1.0e-12,
                "efficiency at {flow} m3/hr and {speed} rpm: {} against {efficiency}",
                chart
                    .polytropic_efficiency(flow, speed)
                    .expect("fan law never refuses")
            );
        }

        let speed = chart.speed_for_head(3974.8455, 70.0 * KJ_PER_KG_TO_M);
        // The class stops its Newton at `1e-6` on the head, so its reported speed can differ
        // in the last few hundredths of an rpm; the fit itself agrees to `1e-12`.
        assert!(
            (speed - 12822.382970080176).abs() < 1.0e-2,
            "speed: {speed}"
        );
        assert_eq!(chart.speed(3974.8455, 70.0 * KJ_PER_KG_TO_M), 12822);

        let flow = chart.flow_for_head(70.0 * KJ_PER_KG_TO_M, 12913.0, 4000.0);
        assert!((flow - 5399.049830646614).abs() < 1.0e-2, "flow: {flow}");
    }

    #[test]
    fn the_generated_boundaries_reproduce_the_class() {
        // `process_anti_surge.tsv`'s `generated_boundaries` block: the class's own
        // `generateSurgeCurve`/`generateStoneWallCurve` on this table, read at three heads
        // (stated there in kJ/kg, as the class' head unit is set to) and at speed 11683.
        let mut chart = CompressorChart::new(&table()).expect("a valid chart");
        chart.generate_surge_curve().expect("a surge curve");
        chart.generate_stonewall_curve().expect("a stonewall curve");
        for (head_kj, expected) in [
            (31.0, 1607.6343245574142),
            (50.0, 2019.1008057284437),
            (80.0, 2788.0154491600683),
        ] {
            let got = chart
                .surge_flow(head_kj * KJ_PER_KG_TO_M)
                .expect("the curve is set");
            assert!(
                (got - expected).abs() < 1.0e-8 * expected.abs(),
                "surge flow at {head_kj} kJ/kg: {got} against {expected}"
            );
        }
        for (head_kj, expected) in [
            (31.0, 5078.449674683021),
            (50.0, 6366.848696674117),
            (80.0, 8428.226023409274),
        ] {
            let got = chart
                .stonewall_flow(head_kj * KJ_PER_KG_TO_M)
                .expect("the curve is set");
            assert!(
                (got - expected).abs() < 1.0e-8 * expected.abs(),
                "stonewall flow at {head_kj} kJ/kg: {got} against {expected}"
            );
        }
        // The per-speed lookups still agree with the generated curve's own extreme points.
        assert!((chart.surge_flow_at_speed(11683.0) - 2415.3793).abs() < 1.0e-9);
        assert!((chart.stonewall_flow_at_speed(11683.0) - 5138.7858).abs() < 1.0e-9);
        // And the limit tests read the same curve.
        let head = 50.0 * KJ_PER_KG_TO_M;
        assert!(chart.is_surge(head, 1000.0));
        assert!(!chart.is_surge(head, 3000.0));
        assert!(chart.is_stonewall(head, 9000.0));
        assert!(!chart.is_stonewall(head, 3000.0));
    }

    /// The three-speed table `CompressorChartSpeedInterpolationTest` reads.
    fn interpolation_table() -> Vec<CurvePoint> {
        let mut points = Vec::new();
        for (speed, heads) in [
            (6000.0, [65.0, 60.0, 50.0]),
            (7000.0, [85.0, 80.0, 70.0]),
            (8000.0, [105.0, 100.0, 90.0]),
        ] {
            for (index, flow) in [1000.0, 2000.0, 3000.0].into_iter().enumerate() {
                points.push(CurvePoint {
                    speed,
                    flow,
                    head: heads[index],
                    efficiency: [0.75, 0.80, 0.75][index],
                });
            }
        }
        points
    }

    #[test]
    fn the_extrapolating_lookup_reproduces_the_capture() {
        // `process_compressor_map_lookup.tsv`'s `speed_interpolation_extrapolate` block, on
        // the three-speed table. **The pinned class blends in speed and then scales the head by
        // `speed / max(s1, s2)`** - which is why 6500 rpm reads 65, not the later source's 70.
        let chart =
            CompressorChart::with_kind(&interpolation_table(), ChartKind::MapLookupExtrapolate)
                .expect("the table fits");
        for (speed, expected) in [
            (6000.0, 60.0),
            (6001.0, 51.454_288_571_428_58),
            (6500.0, 65.0),
            (7000.0, 80.0),
            (7500.0, 84.375),
            (8000.0, 100.0),
            (0.0, 0.0),
            (5000.0, 50.0),
            (9000.0, 112.5),
        ] {
            let got = chart
                .polytropic_head(2000.0, speed)
                .expect("inside the knots, or the end line");
            assert!(
                (got - expected).abs() < 1.0e-9,
                "head at {speed} rpm: {got} against {expected}"
            );
        }
        assert!(
            (chart
                .polytropic_efficiency(2000.0, 7000.0)
                .expect("a number")
                - 0.80)
                .abs()
                < 1.0e-12
        );
        // A single line scales by the speed ratio too: 60 at 6000 rpm gives 70 at 7000.
        let single = CompressorChart::with_kind(
            &interpolation_table()
                .into_iter()
                .filter(|point| point.speed == 6000.0)
                .collect::<Vec<CurvePoint>>(),
            ChartKind::MapLookupExtrapolate,
        )
        .expect("one line fits");
        assert!((single.polytropic_head(2000.0, 7000.0).expect("a number") - 70.0).abs() < 1.0e-12);
    }

    #[test]
    fn the_base_lookup_is_the_mean_of_the_bracketing_lines() {
        let table = table();
        let chart = CompressorChart::with_kind(&table, ChartKind::MapLookup).expect("the table");
        let flow = 3000.0;
        let lower = chart
            .polytropic_head(flow, 11683.0)
            .expect("the line's own spline");
        let upper = chart
            .polytropic_head(flow, 12298.0)
            .expect("the line's own spline");
        let bracket = chart
            .polytropic_head(flow, 12000.0)
            .expect("between two lines");
        assert!((bracket - 0.5 * (lower + upper)).abs() < 1.0e-9);
        // Outside the speed range the end line is read alone, with no scaling.
        let below = chart.polytropic_head(flow, 7000.0).expect("the end line");
        let lowest = chart.polytropic_head(flow, 8200.0).expect("the end line");
        assert!((below - lowest).abs() < 1.0e-9);
        let above = chart.polytropic_head(flow, 14000.0).expect("the end line");
        let highest = chart.polytropic_head(flow, 12913.0).expect("the end line");
        assert!((above - highest).abs() < 1.0e-9);
    }

    #[test]
    fn a_base_lookup_refuses_a_flow_outside_the_knots() {
        let chart = CompressorChart::with_kind(&table(), ChartKind::MapLookup).expect("the table");
        assert!(chart.polytropic_head(100.0, 11683.0).is_err());
        let extrapolating = CompressorChart::with_kind(&table(), ChartKind::MapLookupExtrapolate)
            .expect("the table");
        assert!(extrapolating.polytropic_head(100.0, 11683.0).is_ok());
    }

    #[test]
    fn the_surge_and_stonewall_rows_are_the_closest_curve_s() {
        let chart = CompressorChart::new(&table()).expect("a valid chart");
        assert!((chart.surge_flow_at_speed(11683.0) - 2415.3793).abs() < 1.0e-9);
        assert!((chart.stonewall_flow_at_speed(11683.0) - 5138.7858).abs() < 1.0e-9);
        assert!((chart.surge_head_at_speed(11683.0) - 65.1576 * KJ_PER_KG_TO_M).abs() < 1.0e-6);
        assert!((chart.stonewall_head_at_speed(11683.0) - 31.6387 * KJ_PER_KG_TO_M).abs() < 1.0e-6);

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
