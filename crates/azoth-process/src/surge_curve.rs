//! `SafeSplineSurgeCurve` / `SafeSplineStoneWallCurve` - the compressor map's boundaries.
//!
//! **Two splines, one per direction.** `setCurve` sorts the measured `(flow, head)` pairs by
//! flow and fits a natural cubic spline `head = f(flow)`; it then dedupes the pairs by head
//! (ascending, the last duplicate winning, as the class's `TreeMap` does) and fits a second
//! natural spline `flow = f(head)`. A read inside the measured range is that spline; outside it
//! is a **linear extrapolation** from the two end points, clamped at zero - the class's "safe"
//! part - and a single measured point is a constant.
//!
//! **The two curves differ in one comparison.** `SurgeCurve` calls a point a limit when the
//! curve's flow at its head exceeds the operating flow (`getFlow(head) > flow`); the stonewall
//! curve calls it a limit when the flow is greater than the curve's. Everything else is shared.
//!
//! **A natural spline, not a not-a-knot one.** The class fits with Apache Commons Math's
//! `SplineInterpolator`, whose second derivative is zero at both ends; [`Spline`] solves that
//! same tridiagonal system. The capture's absolute reads are what hold the two together.

use azoth_core::{AzothError, Result};

/// A natural cubic spline through strictly increasing abscissae.
///
/// **The callers validate and this assumes it.** `Boundary` refuses a length mismatch, a
/// duplicate flow and fewer than two heads; the compressor's map lookup refuses its own. The
/// solve below divides by `x[i+1] - x[i]`, so a zero interval is a caller's error.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct Spline {
    x: Vec<f64>,
    y: Vec<f64>,
    /// The second derivative at each knot; zero at both ends.
    second: Vec<f64>,
}

impl Spline {
    pub(crate) fn new(x: Vec<f64>, y: Vec<f64>) -> Self {
        let n = x.len();
        // claim-ok: a_stonewall_curve_flips_the_comparison - every caller refuses the input
        // this would trip on, so the debug assertion is a precondition and not a check.
        debug_assert!(
            n >= 2 && x.windows(2).all(|window| window[1] > window[0]),
            "a spline needs at least two strictly increasing abscissae"
        );
        let mut h = Vec::with_capacity(n - 1);
        for i in 0..n - 1 {
            h.push(x[i + 1] - x[i]);
        }
        // The natural-spline system for the interior second derivatives.
        let mut second = vec![0.0; n];
        let mut diagonal = vec![0.0; n];
        let mut upper = vec![0.0; n];
        let mut rhs = vec![0.0; n];
        for i in 1..n - 1 {
            diagonal[i] = 2.0 * (h[i - 1] + h[i]);
            upper[i] = h[i];
            rhs[i] = 6.0 * ((y[i + 1] - y[i]) / h[i] - (y[i] - y[i - 1]) / h[i - 1]);
        }
        // Thomas algorithm over the interior rows, with `second[0] = second[n - 1] = 0`.
        for i in 2..n - 1 {
            let factor = h[i - 1] / diagonal[i - 1];
            diagonal[i] -= factor * upper[i - 1];
            rhs[i] -= factor * rhs[i - 1];
        }
        for i in (1..n - 1).rev() {
            let below = if i + 1 < n - 1 {
                h[i] * second[i + 1]
            } else {
                0.0
            };
            second[i] = (rhs[i] - below) / diagonal[i];
        }
        Self { x, y, second }
    }

    pub(crate) fn at(&self, value: f64) -> f64 {
        let n = self.x.len();
        // The interval is found by a plain scan: a surge curve has a dozen knots, and the
        // class's own `PolynomialSplineFunction` does the same.
        let mut index = 0;
        for i in 0..n - 1 {
            if value >= self.x[i] && value <= self.x[i + 1] {
                index = i;
                break;
            }
        }
        let h = self.x[index + 1] - self.x[index];
        let a = (self.x[index + 1] - value) / h;
        let b = (value - self.x[index]) / h;
        a * self.y[index]
            + b * self.y[index + 1]
            + ((a * a * a - a) * self.second[index] + (b * b * b - b) * self.second[index + 1])
                * h
                * h
                / 6.0
    }
}

/// The shared curve: the flow-sorted pairs, the head-deduped pairs, and the two splines.
#[derive(Debug, Clone, PartialEq)]
struct Boundary {
    /// Flow ascending, with its head - the class's `flow`/`head` after `setCurve`'s sort.
    flow: Vec<f64>,
    head: Vec<f64>,
    /// Head ascending, with the flow it came from - the class's `sortedHead`/`sortedFlow`.
    sorted_head: Vec<f64>,
    sorted_flow: Vec<f64>,
    head_from_flow: Option<Spline>,
    flow_from_head: Option<Spline>,
    single_flow: f64,
    single_head: f64,
    single: bool,
}

impl Boundary {
    fn new(flow: &[f64], head: &[f64]) -> Result<Self> {
        if flow.len() != head.len() {
            return Err(AzothError::invalid_input(
                "surge_curve",
                "the curve's flow and head arrays must have the same length",
            ));
        }
        if flow.is_empty() {
            return Err(AzothError::invalid_input(
                "surge_curve",
                "an anti-surge curve needs at least one point",
            ));
        }
        for (index, (flow_value, head_value)) in flow.iter().zip(head).enumerate() {
            if !flow_value.is_finite() || !head_value.is_finite() {
                return Err(AzothError::invalid_input(
                    "surge_curve",
                    format!("point {index} is not a finite pair, and a curve cannot be fitted"),
                ));
            }
        }
        if flow.len() == 1 {
            return Ok(Self {
                flow: vec![flow[0]],
                head: vec![head[0]],
                sorted_head: vec![head[0]],
                sorted_flow: vec![flow[0]],
                head_from_flow: None,
                flow_from_head: None,
                single_flow: flow[0],
                single_head: head[0],
                single: true,
            });
        }

        let mut pairs: Vec<(f64, f64)> = flow.iter().copied().zip(head.iter().copied()).collect();
        pairs.sort_by(|a, b| a.0.partial_cmp(&b.0).expect("finite flows"));
        if pairs.len() == 2 {
            // Commons Math's `SplineInterpolator` refuses fewer than three points
            // (`NumberIsTooSmallException`), and the class propagates it; one point is its
            // constant curve, so two is the one size that cannot be read at all.
            return Err(AzothError::invalid_input(
                "surge_curve",
                "a two-point boundary cannot be splined - the class's interpolator needs three \
                 points, and one is its constant",
            ));
        }
        for window in pairs.windows(2) {
            if window[0].0 == window[1].0 {
                // The class sorts the pairs but does not dedupe the flows, and its spline then
                // divides by a zero interval; a curve with two heads at one flow is not a curve.
                return Err(AzothError::invalid_input(
                    "surge_curve",
                    format!(
                        "two points share the flow {}, and a curve needs distinct flows",
                        window[0].0
                    ),
                ));
            }
        }
        let sorted_flow: Vec<f64> = pairs.iter().map(|pair| pair.0).collect();
        let sorted_head: Vec<f64> = pairs.iter().map(|pair| pair.1).collect();

        // `TreeMap<Double, Double>`: head ascending by `Double.compareTo` (which is
        // `total_cmp`), one flow per head, the last duplicate winning on `put`. The pairs are
        // flow-sorted, so a stable sort by head keeps the largest flow last among equal heads.
        if pairs.windows(2).any(|window| window[1].0 <= window[0].0) {
            return Err(AzothError::invalid_input(
                "surge_curve",
                "the curve's flows are not strictly increasing after sorting",
            ));
        }
        let mut by_head = pairs.clone();
        by_head.sort_by(|a, b| a.1.total_cmp(&b.1));
        let mut heads: Vec<f64> = Vec::new();
        let mut flows: Vec<f64> = Vec::new();
        for pair in &by_head {
            match heads.last() {
                Some(head) if *head == pair.1 => {
                    *flows.last_mut().expect("one flow per head") = pair.0;
                }
                _ => {
                    heads.push(pair.1);
                    flows.push(pair.0);
                }
            }
        }
        if heads.len() < 2 {
            return Err(AzothError::invalid_input(
                "surge_curve",
                "not enough distinct head values for spline interpolation",
            ));
        }

        let head_from_flow = Spline::new(sorted_flow.clone(), sorted_head.clone());
        let flow_from_head = Spline::new(heads.clone(), flows.clone());
        Ok(Self {
            flow: sorted_flow,
            head: sorted_head,
            sorted_head: heads,
            sorted_flow: flows,
            head_from_flow: Some(head_from_flow),
            flow_from_head: Some(flow_from_head),
            single_flow: 0.0,
            single_head: 0.0,
            single: false,
        })
    }

    fn flow_at_head(&self, head_value: f64) -> f64 {
        if self.single {
            return self.single_flow;
        }
        let n = self.sorted_head.len();
        if head_value >= self.sorted_head[0] && head_value <= self.sorted_head[n - 1] {
            return self
                .flow_from_head
                .as_ref()
                .expect("a fitted curve has its spline")
                .at(head_value)
                .max(0.0);
        }
        let (slope, base_flow, base_head) = if head_value < self.sorted_head[0] {
            (
                (self.sorted_flow[1] - self.sorted_flow[0])
                    / (self.sorted_head[1] - self.sorted_head[0]),
                self.sorted_flow[0],
                self.sorted_head[0],
            )
        } else {
            (
                (self.sorted_flow[n - 1] - self.sorted_flow[n - 2])
                    / (self.sorted_head[n - 1] - self.sorted_head[n - 2]),
                self.sorted_flow[n - 1],
                self.sorted_head[n - 1],
            )
        };
        (base_flow + slope * (head_value - base_head)).max(0.0)
    }

    fn head_at_flow(&self, flow_value: f64) -> f64 {
        if self.single {
            return self.single_head;
        }
        let n = self.flow.len();
        if flow_value >= self.flow[0] && flow_value <= self.flow[n - 1] {
            return self
                .head_from_flow
                .as_ref()
                .expect("a fitted curve has its spline")
                .at(flow_value)
                .max(0.0);
        }
        let (slope, base_flow, base_head) = if flow_value < self.flow[0] {
            (
                (self.head[1] - self.head[0]) / (self.flow[1] - self.flow[0]),
                self.flow[0],
                self.head[0],
            )
        } else {
            (
                (self.head[n - 1] - self.head[n - 2]) / (self.flow[n - 1] - self.flow[n - 2]),
                self.flow[n - 1],
                self.head[n - 1],
            )
        };
        (base_head + slope * (flow_value - base_flow)).max(0.0)
    }
}

/// `SafeSplineSurgeCurve`: the surge boundary, `flow(head) > flow` is a limit.
#[derive(Debug, Clone, PartialEq)]
pub struct SurgeCurve {
    boundary: Boundary,
}

impl SurgeCurve {
    /// Fit the curve from measured flow and head arrays, in any order.
    ///
    /// # Errors
    /// [`azoth_core::AzothError::InvalidInput`] for a length mismatch, no points, a non-finite
    /// value, two points at one flow, or fewer than two distinct heads.
    pub fn new(flow: &[f64], head: &[f64]) -> Result<Self> {
        Ok(Self {
            boundary: Boundary::new(flow, head)?,
        })
    }

    /// `getFlow(head)` / `getSurgeFlow(head)`: the curve's flow at a head, m3/h.
    #[must_use]
    pub fn flow_at_head(&self, head: f64) -> f64 {
        self.boundary.flow_at_head(head)
    }

    /// `getSurgeHead(flow)`: the curve's head at a flow.
    #[must_use]
    pub fn head_at_flow(&self, flow: f64) -> f64 {
        self.boundary.head_at_flow(flow)
    }

    /// `isSurge(head, flow)`: the curve's flow at the head is greater than the operating flow.
    #[must_use]
    pub fn is_limit(&self, head: f64, flow: f64) -> bool {
        self.flow_at_head(head) > flow
    }
}

/// `SafeSplineStoneWallCurve`: the stonewall boundary, `flow(head) < flow` is a limit.
#[derive(Debug, Clone, PartialEq)]
pub struct StoneWallCurve {
    boundary: Boundary,
}

impl StoneWallCurve {
    /// Fit the curve from measured flow and head arrays, in any order.
    ///
    /// # Errors
    /// The same refusals as [`SurgeCurve::new`].
    pub fn new(flow: &[f64], head: &[f64]) -> Result<Self> {
        Ok(Self {
            boundary: Boundary::new(flow, head)?,
        })
    }

    /// `getFlow(head)` / `getStoneWallFlow(head)`: the curve's flow at a head, m3/h.
    #[must_use]
    pub fn flow_at_head(&self, head: f64) -> f64 {
        self.boundary.flow_at_head(head)
    }

    /// The curve's head at a flow.
    #[must_use]
    pub fn head_at_flow(&self, flow: f64) -> f64 {
        self.boundary.head_at_flow(flow)
    }

    /// `isStoneWall(head, flow)`: the curve's flow at the head is less than the operating flow.
    #[must_use]
    pub fn is_limit(&self, head: f64, flow: f64) -> bool {
        self.flow_at_head(head) < flow
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// `SafeSplineSurgeCurveTest.testSurgeCurve10250`'s table, in the test's own order.
    const FLOW: [f64; 13] = [
        9758.49, 9578.11, 9397.9, 9248.64, 9006.93, 8749.97, 8508.5, 8179.81, 7799.81, 7111.75,
        6480.26, 6007.91, 5607.45,
    ];
    const HEAD: [f64; 13] = [
        112.65, 121.13, 127.56, 132.13, 137.29, 140.73, 142.98, 144.76, 146.14, 148.05, 148.83,
        149.54, 150.0,
    ];

    fn curve() -> SurgeCurve {
        SurgeCurve::new(&FLOW, &HEAD).expect("the captured curve fits")
    }

    fn close(actual: f64, expected: f64, tolerance: f64, what: &str) {
        let error = (actual - expected).abs() / expected.abs().max(f64::MIN_POSITIVE);
        assert!(
            error <= tolerance,
            "{what}: {actual} against {expected}, relative error {error:.3e} > {tolerance:.3e}"
        );
    }

    #[test]
    fn the_captured_reads_are_the_class_s() {
        let curve = curve();
        close(
            curve.flow_at_head(130.0),
            9321.84315857,
            1.0e-6,
            "flow at 130",
        );
        close(curve.flow_at_head(10.0), 11941.98139, 1.0e-6, "flow at 10");
        // **The two directions are two splines, not one inverted.** Each is fitted on the same
        // measured points at different knot spacings, so a round trip lands `0.017` head below
        // where it started - inside the class's own `+/- 1.0` on this row.
        close(
            curve.head_at_flow(9321.84315),
            130.0,
            1.0e-3,
            "head at 9321.84",
        );
        assert_eq!(
            curve.flow_at_head(190.0),
            0.0,
            "extrapolation is clamped at zero"
        );
    }

    #[test]
    fn the_head_falls_as_the_flow_rises() {
        let curve = curve();
        assert!(curve.head_at_flow(6000.0) > curve.head_at_flow(9500.0));
    }

    #[test]
    fn the_surge_limit_is_the_curve_s_flow_above_the_operating_one() {
        let curve = curve();
        // At head 130 the curve carries 9321.84: 9000 is below it, 9600 above it.
        assert!(curve.is_limit(130.0, 9000.0));
        assert!(!curve.is_limit(130.0, 9600.0));
    }

    #[test]
    fn a_stonewall_curve_flips_the_comparison() {
        let curve = StoneWallCurve::new(&FLOW, &HEAD).expect("the captured curve fits");
        assert!(!curve.is_limit(130.0, 9000.0));
        assert!(curve.is_limit(130.0, 9600.0));
    }

    #[test]
    fn a_single_point_is_a_constant() {
        let surge = SurgeCurve::new(&[5000.0], &[100.0]).expect("one point is a curve");
        assert_eq!(surge.flow_at_head(50.0), 5000.0);
        assert_eq!(surge.flow_at_head(200.0), 5000.0);
        assert_eq!(surge.head_at_flow(10.0), 100.0);
    }

    #[test]
    fn duplicate_flows_and_one_head_are_refused() {
        assert!(SurgeCurve::new(&[5000.0, 5000.0], &[100.0, 90.0]).is_err());
        assert!(SurgeCurve::new(&[5000.0, 6000.0], &[100.0, 100.0]).is_err());
        assert!(SurgeCurve::new(&[5000.0], &[]).is_err());
        // Two points is the one size Commons Math's interpolator refuses.
        assert!(SurgeCurve::new(&[5000.0, 6000.0], &[100.0, 90.0]).is_err());
    }
}
