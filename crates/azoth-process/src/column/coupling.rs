//! `DistillationColumn`'s pressure-drop coupling, as a solve.
//!
//! **This is a solver change and not a report.** `hydraulicPressureDropCouplingEnabled` makes
//! `updatePressureProfileFromHydraulics` run inside `updateColumnTearVariables`, and
//! `hasActiveColumnTearVariables` - `!sideDrawSpecifications.isEmpty() ||
//! !pumparounds.isEmpty() || hydraulicPressureDropCouplingEnabled` - is what puts the column on
//! `solveWithColumnTearVariables`. So enabling it changes *which loop* the column takes.
//!
//! **And that loop converges here, which the plan doubted.** The coordinated loop's side-draw
//! row does not converge (`process_column_tear.tsv`: 18 of 30 candidates rejected, a `4.4e-4`
//! residual, and an empty candidate history with a pumparound beside it), and the plan reasoned
//! that this "is a new tear of the same kind". `process_hydraulic_coupling.tsv` measures the
//! opposite: **2 tear iterations, 0 rejected candidates, 0 rollbacks**, on every row. The
//! difference is the tear variable's shape - that row searches a candidate list over several
//! coupled variables, while this one is a single scalar updated by a contraction.
//!
//! # What each pass does
//!
//! `updatePressureProfileFromHydraulics` builds a designer at the stated internals type - **its
//! own constructor defaults and nothing else**, `calcColumnInternals` sets no geometry - takes
//! `getTotalPressureDrop`, and calls `applyHydraulicPressureDrop`, which rewrites **one end** so
//! the difference between the two equals that drop: `bottom = top + drop` where the top is
//! positive, `top = max(1e-6 bara, bottom - drop)` otherwise, at `drop` in bar. The residual it
//! answers is the relative change it made to the end it rewrote, and a residual above `1e-12`
//! is what marks the tear variables as changed.

use azoth_core::units::{Pressure, pascals};
use azoth_core::{AzothError, Result};

use crate::column::designer::{
    DEFAULT_DESIGNER_FLOOD_FRACTION, DEFAULT_DOWNCOMMER_AREA_FRACTION, DEFAULT_HOLE_AREA_FRACTION,
    DEFAULT_HOLE_DIAMETER_MM, DEFAULT_TRAY_SPACING_M, DEFAULT_WEIR_HEIGHT_M, DesignerGeometry,
    UNSIZED_COLUMN_DIAMETER_M, designer_report,
};
use crate::kernels::distillation_column::{ColumnOutcome, ColumnSetup, distillation_column};

/// `DistillationColumn.columnTearTolerance`'s own initialiser.
pub const COLUMN_TEAR_TOLERANCE: f64 = 1.0e-4;
/// `DistillationColumn.maxColumnTearIterations`' own initialiser. **`getColumnTearIterationLimit`
/// takes the larger of this and `maxPumparoundIterations`**, whose initialiser is `8` - so the
/// limit this port uses is `12` and not `8`.
pub const MAX_COLUMN_TEAR_ITERATIONS: usize = 12;
/// `DistillationColumn.maxPumparoundIterations`' initialiser, which the limit is taken against.
pub const MAX_PUMPAROUND_ITERATIONS: usize = 8;
/// The floor `applyHydraulicPressureDrop` clamps a rewritten top pressure to, in bara.
const MIN_TOP_PRESSURE_BARA: f64 = 1.0e-6;
/// The residual above which the tear's variables count as changed.
const CHANGED_RESIDUAL: f64 = 1.0e-12;

/// The coupling's own settings: the internals type and the tear's two knobs.
#[derive(Debug, Clone)]
pub struct HydraulicCoupling {
    /// `hydraulicPressureDropInternalsType`, which `enableHydraulicPressureDropCoupling` sets
    /// beside the flag. **The class's own default is `sieve`**, and it is the designer's *type*
    /// alone: `calcColumnInternals` sets no geometry.
    pub internals_type: String,
    /// `columnTearTolerance`.
    pub tolerance: f64,
    /// `maxColumnTearIterations`, which the limit is taken against the pumparound cap.
    pub max_iterations: usize,
}

impl Default for HydraulicCoupling {
    fn default() -> Self {
        Self {
            internals_type: "sieve".to_string(),
            tolerance: COLUMN_TEAR_TOLERANCE,
            max_iterations: MAX_COLUMN_TEAR_ITERATIONS,
        }
    }
}

/// What the coupled solve did, which is the class's own tear bookkeeping.
#[derive(Debug, Clone, PartialEq)]
pub struct ColumnTear {
    /// `getLastColumnTearIterationCount`: how many passes the outer loop ran.
    pub iterations: u32,
    /// `getLastColumnTearResidual`: the last relative change the coupling made.
    pub residual: f64,
    /// `isLastColumnTearConverged`.
    pub converged: bool,
    /// `getLastColumnTearInnerIterationCount`: **the sweeps of every solve the loop ran, summed**
    /// - which is how a two-pass loop shows `59` where the last solve alone took `22`.
    pub inner_iterations: u32,
    /// `getLastHydraulicPressureDropPa`: the drop the last designer summed, Pa.
    pub hydraulic_pressure_drop: Pressure,
    /// The top pressure the solve ended on, which the coupling rewrote.
    pub top_pressure: Pressure,
    /// The bottom pressure the solve ended on.
    pub bottom_pressure: Pressure,
}

/// `DistillationColumn`'s coupled solve: the column's own solve inside the tear the coupling puts
/// it on.
///
/// # Errors
/// * [`AzothError::InvalidInput`] where a side-draw flow specification or a pumparound return is
///   stated beside the coupling. **The class does take that combination** - with the coupling on,
///   `hasSingleSideDrawFlowSpecificationOnly` and `hasPumparoundTearVariablesOnly` are both false,
///   so the general loop runs with the pressure *and* the draw among its tear variables - and
///   this port carries the pressure variable alone, so it refuses rather than solving half of it.
///   The pair's own measurement is `process_column_tear.tsv`'s.
/// * Whatever the column's own solve raises, at any pass.
pub fn distillation_column_coupled(
    setup: &ColumnSetup,
    coupling: &HydraulicCoupling,
) -> Result<(ColumnOutcome, ColumnTear)> {
    if !setup.side_draw_flows.is_empty() || !setup.pumparound_returns.is_empty() {
        return Err(AzothError::invalid_input(
            "hydraulic_pressure_drop_coupling",
            "the pressure-drop coupling is stated beside a side-draw flow specification or a \
             pumparound return: the class runs all of them as coordinated tear variables of one \
             loop, and this port carries the pressure alone",
        ));
    }

    let limit = coupling
        .max_iterations
        .max(MAX_PUMPAROUND_ITERATIONS)
        .max(1);
    // **The state the loop mutates**, which is the setup with its two ends rewritten.
    let mut state = setup.clone();
    let mut inner_iterations = 0_u32;
    let mut drop = pascals(0.0);
    let mut residual = f64::INFINITY;

    for pass in 0..limit {
        let solved = distillation_column(&state)?;
        inner_iterations += solved.iterations;
        let (pass_residual, pass_drop) = update_tear_variables(&mut state, &solved, coupling)?;
        residual = pass_residual;
        drop = pass_drop;
        let iterations = u32::try_from(pass + 1).unwrap_or(u32::MAX);

        if residual <= coupling.tolerance {
            // **A last solve where the update moved something**, which is the class's own
            // re-solve inside its convergence branch.
            if residual > CHANGED_RESIDUAL {
                state.initial_state = None;
                let again = distillation_column(&state)?;
                inner_iterations += again.iterations;
                return Ok((
                    again,
                    ColumnTear {
                        iterations,
                        residual,
                        converged: true,
                        inner_iterations,
                        hydraulic_pressure_drop: drop,
                        top_pressure: state.top_pressure,
                        bottom_pressure: state.bottom_pressure,
                    },
                ));
            }
            return Ok((
                solved,
                ColumnTear {
                    iterations,
                    residual,
                    converged: true,
                    inner_iterations,
                    hydraulic_pressure_drop: drop,
                    top_pressure: state.top_pressure,
                    bottom_pressure: state.bottom_pressure,
                },
            ));
        }
        if residual <= CHANGED_RESIDUAL {
            // Nothing moved, so another pass would repeat this one: the class stops here.
            return Ok((
                solved,
                ColumnTear {
                    iterations,
                    residual,
                    converged: false,
                    inner_iterations,
                    hydraulic_pressure_drop: drop,
                    top_pressure: state.top_pressure,
                    bottom_pressure: state.bottom_pressure,
                },
            ));
        }
        // `setDoInitializion(true)`: the next pass re-seeds rather than warm-starting.
        state.initial_state = None;
    }

    // The cap, and the class's own last solve before it reports where it stopped.
    state.initial_state = None;
    let again = distillation_column(&state)?;
    inner_iterations += again.iterations;
    let tear = ColumnTear {
        iterations: u32::try_from(limit).unwrap_or(u32::MAX),
        residual,
        converged: residual <= coupling.tolerance,
        inner_iterations,
        hydraulic_pressure_drop: drop,
        top_pressure: state.top_pressure,
        bottom_pressure: state.bottom_pressure,
    };
    Ok((again, tear))
}

/// `updateColumnTearVariables`' pressure half: the designer's drop, and the end it rewrites.
fn update_tear_variables(
    state: &mut ColumnSetup,
    solved: &ColumnOutcome,
    coupling: &HydraulicCoupling,
) -> Result<(f64, Pressure)> {
    // **`calcColumnInternals(type)`: the designer's own constructor defaults and nothing else.**
    // The class sets the internals type and calls `calculate`; the tray spacing, the weir height,
    // the hole geometry and the flood fraction are all the designer's initialisers.
    let report = designer_report(
        &solved.trays,
        &solved.distillate.components,
        &DesignerGeometry {
            internals_type: coupling.internals_type.clone(),
            tray_spacing: azoth_core::units::meters(DEFAULT_TRAY_SPACING_M),
            weir_height: azoth_core::units::meters(DEFAULT_WEIR_HEIGHT_M),
            hole_diameter: azoth_core::units::millimeters(DEFAULT_HOLE_DIAMETER_MM),
            hole_area_fraction: DEFAULT_HOLE_AREA_FRACTION,
            downcommer_area_fraction: DEFAULT_DOWNCOMMER_AREA_FRACTION,
            design_flood_fraction: DEFAULT_DESIGNER_FLOOD_FRACTION,
            column_diameter_override: azoth_core::units::meters(UNSIZED_COLUMN_DIAMETER_M),
        },
    )?;
    let drop = pascals(report.total_pressure_drop.value.max(0.0));
    let residual = apply_hydraulic_pressure_drop(state, drop);
    Ok((residual, drop))
}

/// `applyHydraulicPressureDrop`: rewrite one end, and answer the relative change.
///
/// **The top branch wins where the top is positive**, which is the class's own order, and the
/// other end is left where it was. The relative change is taken against the value the rewritten
/// end had, so a first pass that moves the bottom from `20` to `19.0178` answers `0.0491` and a
/// second that moves it by `7.6e-7` answers under the tolerance.
fn apply_hydraulic_pressure_drop(state: &mut ColumnSetup, drop: Pressure) -> f64 {
    let previous_top = state.top_pressure.value;
    let previous_bottom = state.bottom_pressure.value;

    if previous_top.is_finite() && previous_top > 0.0 {
        let rewritten = previous_top + drop.value;
        state.bottom_pressure = pascals(rewritten);
        if !(previous_bottom.is_finite() && previous_bottom > 0.0) {
            return 1.0;
        }
        return (rewritten - previous_bottom).abs() / previous_bottom.abs().max(CHANGED_RESIDUAL);
    }
    if previous_bottom.is_finite() && previous_bottom > 0.0 {
        let rewritten = (previous_bottom - drop.value).max(MIN_TOP_PRESSURE_BARA * 1.0e5);
        state.top_pressure = pascals(rewritten);
        if !(previous_top.is_finite() && previous_top > 0.0) {
            return 1.0;
        }
        return (rewritten - previous_top).abs() / previous_top.abs().max(CHANGED_RESIDUAL);
    }
    0.0
}

#[cfg(test)]
mod test {
    use azoth_core::units::{kelvins, pascals};

    use crate::kernels::distillation_column::{
        ColumnSetup, SideDrawFlow, SideDrawPhase, SolverType, Specification,
    };
    use crate::stream::Stream;

    use super::{ColumnTear, HydraulicCoupling, distillation_column_coupled};

    /// The capture's own binary column, which every row of this probe starts from.
    fn binary() -> ColumnSetup {
        ColumnSetup {
            feed: Stream::from_pt(
                vec!["methane".into(), "n-butane".into()],
                vec![0.5, 0.5],
                7.490704036290964,
                pascals(20.0e5),
                kelvins(300.0),
            )
            .expect("the fluid resolves"),
            feed_stage: 2,
            number_of_stages: 4,
            has_reboiler: true,
            has_condenser: true,
            top_pressure: pascals(19.0e5),
            bottom_pressure: pascals(20.0e5),
            condenser_temperature: Some(kelvins(253.15)),
            reboiler_temperature: Some(kelvins(373.15)),
            temperature_tolerance: 1.0e-6,
            murphree_efficiency: None,
            absorber_murphree: None,
            initial_state: None,
            max_iterations: 200,
            top_specification: None::<Specification>,
            bottom_specification: None,
            top_feed: None,
            tray_temperatures: None,
            solver_type: SolverType::DirectSubstitution,
            reactive: crate::kernels::ReactiveSection::None,
            gas_side_draw_fractions: None,
            liquid_side_draw_fractions: None,
            pumparound_fractions: None,
            side_draw_flows: Vec::new(),
            pumparound_returns: Vec::new(),
            pumparound_inlets: Vec::new(),
            pumparound_tolerance: None,
            pumparound_max_iterations: None,
        }
    }

    /// The fields the capture's two coupled rows hold, held here.
    fn check_tear(tear: &ColumnTear, outcome: &crate::kernels::distillation_column::ColumnOutcome) {
        let close = |got: f64, want: f64| (got - want).abs() / want.abs() < 1e-5;
        // **The tear's residual is a difference of two nearly-equal pressures**, so it is not a
        // relative-tolerance quantity: the two libraries' ends agree to `1e-11` of `1.9e6` Pa
        // while the residual itself is `7.6e-7`, and the ratio of the two differences is what a
        // tight gate would be measuring. It is held to five significant figures instead.
        let close_residual = |got: f64, want: f64| (got - want).abs() / want.abs() < 1e-4;
        assert!(tear.converged, "the coupling converges on this state");
        assert_eq!(tear.iterations, 2, "the tear's own pass count");
        assert!(
            close_residual(tear.residual, 7.564_364_502_033_816e-7),
            "tear residual: {}",
            tear.residual
        );
        // **The sweep counts are not oracle quantities**, and the capture's are NeqSim's: it
        // sums `15 + 22 + 22 = 59` over the loop's three solves where this port's own solve takes
        // `12 + 14 + 14 = 40` on the same states - the same standing difference the case file
        // records, where `iterations` is deliberately not an expected value. What *is* the
        // oracle is the pass count, the drop, the rewritten end and the residual above.
        assert_eq!(tear.inner_iterations, 40, "this port's own sweeps, summed");
        assert!(
            close(tear.hydraulic_pressure_drop.value, 1_781.574_745_225_591),
            "the summed drop: {}",
            tear.hydraulic_pressure_drop.value
        );
        assert!(
            close(tear.bottom_pressure.value, 19.017_815_747_452_257e5),
            "the rewritten bottom: {}",
            tear.bottom_pressure.value
        );
        // **And the profile carries it**: the rewrite is not a scalar the tear keeps, it is the
        // bottom stage's own pressure, so the solve that follows runs at it.
        assert!(
            close(outcome.trays[0].pressure.value, 19.017_815_747_452_257e5),
            "the reboiler's pressure: {}",
            outcome.trays[0].pressure.value
        );
        assert!(
            close(
                outcome.trays.last().expect("the condenser").pressure.value,
                19.0e5
            ),
            "the condenser keeps the stated top"
        );
    }

    /// **`process_hydraulic_coupling.tsv`'s `binary_coupling_on` row, end to end.**
    #[test]
    fn the_coupling_converges_on_the_captured_binary_column() {
        let (outcome, tear) = distillation_column_coupled(&binary(), &HydraulicCoupling::default())
            .expect("the coupled solve converges");
        check_tear(&tear, &outcome);
    }

    /// **The valve type moves the drop and nothing else**, which the capture's third row shows.
    #[test]
    fn the_internals_type_moves_the_drop_alone() {
        let (_, sieve) = distillation_column_coupled(&binary(), &HydraulicCoupling::default())
            .expect("the coupled solve converges");
        let (_, valve) = distillation_column_coupled(
            &binary(),
            &HydraulicCoupling {
                internals_type: "valve".to_string(),
                ..HydraulicCoupling::default()
            },
        )
        .expect("the coupled solve converges");
        assert_eq!(
            valve.iterations, 2,
            "the pass count is the tear's, not the type's"
        );
        assert!(
            valve.hydraulic_pressure_drop.value < sieve.hydraulic_pressure_drop.value,
            "the valve tray drops less: {} against {}",
            valve.hydraulic_pressure_drop.value,
            sieve.hydraulic_pressure_drop.value
        );
        // The bottom the capture prints for the valve row, to five figures.
        let close = |got: f64, want: f64| (got - want).abs() / want.abs() < 1e-5;
        assert!(
            close(valve.bottom_pressure.value, 19.017_766_632_564_35e5),
            "the valve row's bottom: {}",
            valve.bottom_pressure.value
        );
    }

    /// **A draw beside the coupling is refused rather than half-solved.**
    #[test]
    fn a_side_draw_flow_beside_the_coupling_is_refused() {
        let mut setup = binary();
        setup.side_draw_flows = vec![SideDrawFlow {
            tray: 3,
            phase: SideDrawPhase::Liquid,
            target: 0.1,
            tolerance: 1.0e-4,
            max_iterations: 12,
        }];
        let error = distillation_column_coupled(&setup, &HydraulicCoupling::default())
            .expect_err("the pair is refused");
        assert!(
            format!("{error}").contains("coordinated tear variables"),
            "{error}"
        );
    }
}
