//! The column solve, against `validation/neqsim/captures/process_column.tsv`.
//!
//! Four rows, and they carry three different statements between them. A **converged** binary
//! column, which is the port's oracle; the same column at a **looser tolerance**, which is a
//! different measurement rather than a sloppier version of the first; and NeqSim's own
//! **deethanizer regression state**, which on a Peng-Robinson fluid does not converge at all -
//! so it is captured as evidence and declared uncased, and the only thing asserted about it is
//! that this port's answer to the same column is self-consistent.

use azoth_core::units::{kelvins, pascals};
use azoth_process::Stream;
use azoth_process::kernels::distillation_column::{ColumnSetup, SolverType, distillation_column};

fn binary_feed() -> Stream {
    Stream::from_pt(
        vec!["methane".into(), "n-butane".into()],
        vec![0.5, 0.5],
        7.490704036290964,
        pascals(20.0e5),
        kelvins(300.0),
    )
    .expect("the fluid resolves")
}

/// The captured binary column: four stages, the feed on stage 2, a condenser at -20 C and a
/// reboiler at 100 C, 19 bara at the top and 20 at the bottom.
fn binary_column(tolerance: f64) -> ColumnSetup {
    ColumnSetup {
        feed: binary_feed(),
        feed_stage: 2,
        number_of_stages: 4,
        has_reboiler: true,
        has_condenser: true,
        top_pressure: pascals(19.0e5),
        bottom_pressure: pascals(20.0e5),
        condenser_temperature: Some(kelvins(253.15)),
        reboiler_temperature: Some(kelvins(373.15)),
        temperature_tolerance: tolerance,
        max_iterations: 200,
        top_specification: None,
        bottom_specification: None,
        top_feed: None,
        tray_temperatures: None,
        solver_type: SolverType::DirectSubstitution,
        reactive: azoth_process::kernels::ReactiveSection::None,
        gas_side_draw_fractions: None,
        liquid_side_draw_fractions: None,
        pumparound_fractions: None,
    }
}

fn relative(actual: f64, expected: f64, tolerance: f64, what: &str) {
    let scale = 1.0 + actual.abs().max(expected.abs());
    assert!(
        (actual - expected).abs() <= tolerance * scale,
        "{what}: {actual} vs {expected}, {} relative",
        (actual - expected).abs() / scale
    );
}

fn absolute(actual: f64, expected: f64, tolerance: f64, what: &str) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{what}: {actual} vs {expected}, {} absolute",
        (actual - expected).abs()
    );
}

/// **The whole profile, which is what a column's answer is.**
///
/// A converged column is a fixed point, and NeqSim's and azoth's reach it **one iteration
/// apart** - 15 against 14 - so the agreement below is two solvers finding the same state by
/// different paths rather than one solver replayed. Every tray temperature is within `1.4e-4`
/// K (a relative `4e-7`), every traffic rate within `5e-6` relative, and both products'
/// compositions within `3e-15`.
///
/// The residuals are the two libraries' *own*, and they are not comparable: the class's
/// `getLastMassResidual` is a norm over its MESH equations, and this port's is the products'
/// worst component imbalance against the feed. They are reported side by side in the capture
/// for that reason, and asserted separately here.
#[test]
fn a_column_solves_the_captured_binary_profile() {
    let out = distillation_column(&binary_column(1.0e-6)).expect("the column converges");

    // NeqSim: 15 iterations, a mean tray-temperature change of 6.374e-7 K.
    assert!(
        (out.iterations as i32 - 15).abs() <= 1,
        "iterations: {} against NeqSim's 15",
        out.iterations
    );
    assert!(
        out.temperature_residual <= 1.0e-6,
        "the gate is what the solve was held to: {} K",
        out.temperature_residual
    );

    let temperatures = [
        373.15,
        336.15382381733946,
        303.0710539116486,
        302.1225394582587,
        296.84208939635556,
        253.15,
    ];
    let gas = [
        1.373278181743518,
        0.5259031551756898,
        4.5761815164886315,
        4.558039582740156,
        4.451759798394696,
        3.7914997294607584,
    ];
    let liquid = [
        3.6992043068302056,
        5.072482507910388,
        4.225107543947319,
        0.784681868969297,
        0.7665399361896209,
        0.6602601517057779,
    ];
    assert_eq!(out.trays.len(), 6, "the two ends are trays");
    for (i, tray) in out.trays.iter().enumerate() {
        absolute(
            tray.temperature.value,
            temperatures[i],
            1.0e-3,
            "tray temperature",
        );
        relative(tray.gas_n, gas[i], 1.0e-5, "tray vapour");
        relative(tray.liquid_n, liquid[i], 1.0e-5, "tray liquid");
        // The pressure profile is the class's own linear rule, and it is exact.
        absolute(
            tray.pressure.value,
            (20.0 - i as f64 * 0.2) * 1.0e5,
            1.0,
            "tray pressure",
        );
    }

    // The products, which are the ends' own outlets.
    relative(out.distillate.n, 3.7914997294607584, 1.0e-5, "distillate");
    relative(
        out.distillate.z[0],
        0.9659099052324753,
        1.0e-5,
        "the distillate's methane",
    );
    relative(out.bottoms.n, 3.6992043068302056, 1.0e-5, "bottoms");
    relative(
        out.bottoms.z[1],
        0.9775343702231742,
        1.0e-5,
        "the bottoms' n-butane",
    );

    // **The closure the class's own tests assert**: the products reconcile against the feed.
    relative(
        out.distillate.n + out.bottoms.n,
        binary_feed().n,
        1.0e-5,
        "the product mole balance",
    );
    for c in 0..2 {
        relative(
            out.distillate.n * out.distillate.z[c] + out.bottoms.n * out.bottoms.z[c],
            binary_feed().n * binary_feed().z[c],
            1.0e-5,
            "the component balance",
        );
    }

    // The duties, to the libraries' enthalpy offset rather than to their digits: the two are
    // enthalpy *differences*, so the offset mostly cancels and what is left is its temperature
    // dependence - measured at `0.39` W of `21323` and `1.06` W of `47786`.
    absolute(
        out.condenser_duty.value,
        -21323.042789303567,
        5.0,
        "condenser duty",
    );
    absolute(
        out.reboiler_duty.value,
        47786.58389381015,
        5.0,
        "reboiler duty",
    );

    // And this port's own closure, which is the two residuals it reports.
    assert!(
        out.mass_residual < 1.0e-6,
        "mass residual {}",
        out.mass_residual
    );
    assert!(
        out.energy_residual < 1.0e-6,
        "energy residual {}",
        out.energy_residual
    );
}

/// **A looser tolerance is a different measurement, not a sloppier one.** With the gate at
/// `1e-2` NeqSim stops after 7 iterations at `3.7913503170929044` mol/s of distillate against
/// the rigorous row's `3.7914997294607584` - the same separation, four figures of it, and the
/// row is in the capture to say so.
///
/// **The agreement here is `1.7e-4` relative and not `5e-6`**, which is the point: where a
/// solve stops is where its answer is, so a looser gate buys a different state rather than the
/// same one computed faster. The port stops at its own iteration, so the two states are two
/// points on the same approach and the band is the gap between where each stopped.
#[test]
fn a_looser_gate_stops_sooner_at_the_same_answer() {
    let out = distillation_column(&binary_column(1.0e-2)).expect("the column converges");

    assert!(out.iterations < 15, "a looser gate takes fewer iterations");
    relative(out.distillate.n, 3.7913503170929044, 1.0e-3, "distillate");
    relative(out.bottoms.n, 3.6993537191980588, 1.0e-3, "bottoms");
    // The same separation to three figures, which is what makes the tolerance a tolerance.
    absolute(
        out.distillate.z[0],
        0.9659,
        1.0e-3,
        "the distillate's methane",
    );
}

/// **NeqSim's own deethanizer does not converge on a Peng-Robinson fluid, and this port
/// does.** That divergence is the row's finding, and it makes the row uncased: the capture
/// holds the class's *partial* state as evidence, and the port's state is a different one.
///
/// `NaphtaliSandholmPublishedStateTest` writes the column on `SystemSrkEos`; azoth's
/// `Stream::mixture()` resolves PR and has no other route, and on PR the class's hard-capped
/// 80 iterations end at a mean tray-temperature change of `7.3e-3` K against the `1e-5` it asks
/// for. Raising the cap does not rescue it: 200 iterations reach `2.9e-2` K, which is *worse* -
/// an iteration wandering rather than crawling. This port reaches `6.7e-6` K in 32.
///
/// **The named difference is the relaxation controller's input.** The class decides whether to
/// damp from a combined residual carrying the scaled temperature, mass *and* energy terms, and
/// its energy term is a MESH norm it weights up to ten; this port carries the class's
/// temperature term alone, because those two norms are D5's machinery. So where the class
/// damps and crawls, this port keeps a relaxation of one and converges - which is a solver-path
/// difference and not a physics one, and it is recorded rather than papered over: the row is
/// declared uncased in the layer diff, and what is asserted here is that the port's answer is
/// **self-consistent**, not that it is NeqSim's.
#[test]
fn the_deethanizer_converges_here_where_neqsim_does_not() {
    let feed = Stream::from_pt(
        vec![
            "methane".into(),
            "ethane".into(),
            "propane".into(),
            "i-butane".into(),
            "n-butane".into(),
            "i-pentane".into(),
            "n-pentane".into(),
            "n-hexane".into(),
        ],
        vec![0.22, 0.34, 0.20, 0.08, 0.08, 0.03, 0.03, 0.02],
        73.24414949904373,
        pascals(25.0e5),
        kelvins(283.15),
    )
    .expect("the fluid resolves");
    let out = distillation_column(&ColumnSetup {
        feed,
        feed_stage: 6,
        number_of_stages: 8,
        has_reboiler: true,
        has_condenser: true,
        top_pressure: pascals(24.0e5),
        bottom_pressure: pascals(25.0e5),
        condenser_temperature: Some(kelvins(273.15)),
        reboiler_temperature: Some(kelvins(353.15)),
        temperature_tolerance: 1.0e-5,
        max_iterations: 80,
        top_specification: None,
        bottom_specification: None,
        top_feed: None,
        tray_temperatures: None,
        solver_type: SolverType::DirectSubstitution,
        reactive: azoth_process::kernels::ReactiveSection::None,
        gas_side_draw_fractions: None,
        liquid_side_draw_fractions: None,
        pumparound_fractions: None,
    });

    let out = out.expect(
        "this port converges on the deethanizer, where the class does not - see this test's \
         documentation for the measured reason",
    );

    // The class stops at `7.3e-3` K after 80 and at `2.9e-2` K after 200; the port is inside
    // the gate it was given.
    assert!(
        out.temperature_residual <= 1.0e-5,
        "{} K",
        out.temperature_residual
    );
    assert!(out.iterations < 80, "in {} iterations", out.iterations);

    // **Self-consistency, which is the assertion this row can carry**: whatever state it
    // reached, the products reconcile against the feed and the energy closes.
    assert!(
        out.mass_residual < 1.0e-5,
        "mass residual {}",
        out.mass_residual
    );
    assert!(
        out.energy_residual < 1.0e-5,
        "energy residual {}",
        out.energy_residual
    );
    relative(
        out.distillate.n + out.bottoms.n,
        setup_feed_moles(),
        1.0e-6,
        "the product mole balance",
    );
    // And the separation is a separation: the overhead is richer in methane than the feed.
    assert!(
        out.distillate.z[0] > 0.22 && out.bottoms.z[0] < 0.22,
        "a column that does not split its light key is not solving anything: {} and {}",
        out.distillate.z[0],
        out.bottoms.z[0]
    );
}

/// The deethanizer's feed, in mol/s, as the capture's `feed_mol_per_sec`.
fn setup_feed_moles() -> f64 {
    73.24414949904373
}

/// A stage count and a feed stage outside the column are refused, with the range named.
#[test]
fn a_column_refuses_a_stage_it_does_not_have() {
    let mut setup = binary_column(1.0e-6);
    setup.feed_stage = 9;
    assert!(distillation_column(&setup).is_err());

    let mut setup = binary_column(1.0e-6);
    setup.number_of_stages = 0;
    assert!(distillation_column(&setup).is_err());

    let mut setup = binary_column(1.0e-6);
    setup.temperature_tolerance = 0.0;
    assert!(distillation_column(&setup).is_err());
}

/// The three per-tray draw vectors a side-draw test states, in the order a setup carries them.
type DrawVectors = (Option<Vec<f64>>, Option<Vec<f64>>, Option<Vec<f64>>);

/// The captured binary column with one tray's draw stated, which is the same setup three of the
/// tests below vary.
fn drawn_column(vectors: DrawVectors) -> ColumnSetup {
    let mut setup = binary_column(1.0e-6);
    setup.gas_side_draw_fractions = vectors.0;
    setup.liquid_side_draw_fractions = vectors.1;
    setup.pumparound_fractions = vectors.2;
    setup
}

/// A per-tray vector of one kind: a zero everywhere but the named tray.
fn on_tray(index: usize, fraction: f64) -> Vec<f64> {
    let mut vector = vec![0.0; 6];
    vector[index] = fraction;
    vector
}

/// **The column's draw is the tray's own phase, and the tray's own traffic is what is left.**
///
/// Against `validation/neqsim/captures/process_side_draw.tsv` and the same block in
/// `process_column.tsv`: tray 3 of the four-stage column drawing a quarter of its vapour. The
/// capture's `columnReportsSideDrawAsOutletStream` and `columnEnergyBalanceIncludesSideDraw-
/// Streams` are the class's own tests of this, and their numbers are what the case is held to.
/// What is asserted here is the *identity* rather than the numbers, which is the port's own
/// statement of the mechanism: the draw plus what the tray keeps is the vapour the tray formed,
/// and the draw is the stated fraction of it.
///
/// **The draw carries the tray's own state**, so it is at the tray's temperature and pressure
/// and not the product's: measured in the capture at 302.176 K and 19.4 bara against a
/// distillate of 253.15 K at 19.0.
#[test]
fn a_column_draw_is_the_fraction_of_the_trays_own_phase() {
    let drawn = distillation_column(&drawn_column((Some(on_tray(3, 0.25)), None, None)))
        .expect("the column converges with a draw open");

    assert_eq!(drawn.gas_side_draws.len(), 6);
    let kept = drawn.trays[3].gas_n;
    let side = drawn.gas_side_draws[3]
        .as_ref()
        .expect("tray 3 withdraws vapour")
        .n;
    // NeqSim's own row: the tray's vapour is `4.5419` mol/s, of which `3.4064` goes up and
    // `1.1355` is withdrawn.
    relative(
        side + kept,
        4.541891923509715,
        1.0e-5,
        "the tray's own vapour",
    );
    relative(
        side,
        0.25 * (side + kept),
        1.0e-12,
        "the quarter that is withdrawn",
    );
    relative(
        kept,
        3.4064189426322864,
        1.0e-5,
        "the traffic above the tray",
    );
    relative(side, 1.1354729808774289, 1.0e-5, "the draw itself");
    absolute(
        drawn.gas_side_draws[3].as_ref().expect("the draw").t.value,
        drawn.trays[3].temperature.value,
        1.0e-9,
        "the draw is at the tray's own temperature",
    );
    // Every other tray drew nothing, and the vectors say so rather than repeating a flow.
    assert_eq!(
        drawn
            .liquid_side_draws
            .iter()
            .filter(|d| d.is_some())
            .count(),
        0
    );
    assert_eq!(drawn.pumparounds.iter().filter(|d| d.is_some()).count(), 0);
    assert_eq!(
        drawn.gas_side_draws.iter().filter(|d| d.is_some()).count(),
        1
    );
}

/// **The closure counts the draws, which is what the class's own mass-balance test asserts.**
///
/// `columnReportsSideDrawAsOutletStream` asserts `getMassBalance("kg/hr")` is zero to the feed
/// flow x `1e-6`; the port's equivalent is the products' component imbalance against the feed.
/// **A closure that ignored what the trays withdrew would report a phantom imbalance rather
/// than a gate miss**: measured, the two states below read `0.248` and `0.070` where they read
/// `4.1e-8` and `9.7e-9`.
///
/// The gate is `1e-6` because the *solve's own* closure is not tighter than that: this column
/// with no draw at all leaves `3.9e-8`, and the draw states sit in the same band
/// (`1.0e-8` to `4.2e-8`). Holding the draw states to a band the undrawn state cannot meet
/// would be asserting the draw's arithmetic through the solver's stopping rule.
#[test]
fn the_columns_closures_count_what_the_trays_withdrew() {
    for vectors in [
        (Some(on_tray(3, 0.25)), None, None),
        (None, Some(on_tray(3, 0.10)), Some(on_tray(3, 0.05))),
    ] {
        let out = distillation_column(&drawn_column(vectors))
            .expect("the column converges with a draw open");
        let drawn: f64 = out
            .gas_side_draws
            .iter()
            .chain(out.liquid_side_draws.iter())
            .chain(out.pumparounds.iter())
            .flatten()
            .map(|draw| draw.n)
            .sum();
        assert!(drawn > 0.0, "the state under test draws something");
        assert!(
            out.mass_residual < 1.0e-6,
            "mass {} with {drawn} mol/s withdrawn",
            out.mass_residual
        );
        assert!(
            out.energy_residual < 1.0e-6,
            "energy {} with {drawn} mol/s withdrawn",
            out.energy_residual
        );
    }
}

/// **A fraction on an end is refused by name**, because this port's ends are
/// `column::reboiler` and `column::condenser` rather than stages.
#[test]
fn a_fraction_on_an_end_is_refused_by_name() {
    for index in [0, 5] {
        for which in 0..3 {
            let mut vectors = (None, None, None);
            match which {
                0 => vectors.0 = Some(on_tray(index, 0.1)),
                1 => vectors.1 = Some(on_tray(index, 0.1)),
                _ => vectors.2 = Some(on_tray(index, 0.1)),
            }
            let error = distillation_column(&drawn_column(vectors))
                .expect_err("an end has no tray outlet to split");
            let message = error.to_string();
            assert!(
                message.contains("reboiler") || message.contains("condenser"),
                "tray {index} names the end it is not: {message}"
            );
        }
    }
    // And a vector that is not one entry per tray, which is the shape the fraction is stated in.
    let error = distillation_column(&drawn_column((Some(vec![0.1, 0.2]), None, None)))
        .expect_err("six trays need six fractions");
    assert!(error.to_string().contains("2 fraction(s)"), "{error}");
}
