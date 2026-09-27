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
use azoth_process::column::murphree::Murphree;
use azoth_process::column::murphree::{correct_vapour, corrects};
use azoth_process::kernels::distillation_column::ColumnState;
use azoth_process::kernels::distillation_column::{
    ColumnSetup, SideDrawFlow, SideDrawPhase, SolverType, distillation_column,
};

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
        murphree_efficiency: None,
        initial_state: None,
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
        side_draw_flows: Vec::new(),
        pumparound_returns: Vec::new(),
        pumparound_inlets: Vec::new(),
        pumparound_tolerance: None,
        pumparound_max_iterations: None,
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
/// A converged column is a fixed point, and NeqSim's and azoth's reach **the same one** - the
/// profile below agrees to `1.4e-4` K, every traffic rate to `5e-6` relative and both products'
/// compositions to `3e-15` - so this is two solvers finding one state by different paths rather
/// than one solver replayed. The *iteration counts* differ (15 against 12) because the two stop
/// on different gates, and where a solve stops is its own measurement.
///
/// The residuals are the two libraries' *own*, and they are not comparable: the class's
/// `getLastMassResidual` is a norm over its MESH equations, and this port's is the products'
/// worst component imbalance against the feed. They are reported side by side in the capture
/// for that reason, and asserted separately here.
#[test]
fn a_column_solves_the_captured_binary_profile() {
    let out = distillation_column(&binary_column(1.0e-6)).expect("the column converges");

    // NeqSim: 15 iterations, a mean tray-temperature change of 6.374e-7 K. **The count is not
    // the port's to match**: the class stops on three residuals - temperature, a MESH mass norm
    // and the energy closure - where this port stops on the temperature alone, so where each
    // stops is a different measurement rather than a different answer. Measured, the port takes
    // 12 from the class's own seed and reached the same state in 14 before that seed was fixed.
    assert!(
        (out.iterations as i32 - 15).abs() <= 4,
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

/// **The Murphree correction, tray by tray, on the column the ideal row above already pins.**
///
/// The row is `binary_methane_butane_4_stages` with an efficiency of `0.6` stated - one of the
/// two `DistillationSolverBenchmarkTest` uses - so the *difference* from the ideal profile is
/// the correction and nothing else. It is a large difference rather than a perturbation: the
/// distillate goes from `0.96591` to `0.97459` methane and the condenser duty from `-21323` W to
/// `-18651`, and **the correction is reproduced to about `1e-6` on every quantity**: tray
/// temperatures within `1.3e-4` K, traffic rates within `2.3e-6` relative, both products within
/// `4e-8`, and both duties within `1.1` W of `18651` and `61658`.
///
/// **The two ends do not move at all**, which is the class's own two guards made visible: stage
/// 0 is the reboiler and stage 5 the condenser, and `applyMurphreeCorrection` returns on both
/// before it reads anything. Their temperatures are the pins.
///
/// **The ends' *traffic* is the product and not the flash**, because `finalizeSolve` writes the
/// reconciled products back onto the end trays before anything can read them - so `tray0_liquid`
/// is the bottoms and `tray5_gas` the distillate, which is what the capture prints.
#[test]
fn a_murphree_efficiency_reproduces_the_captured_corrected_profile() {
    let mut setup = binary_column(1.0e-6);
    setup.murphree_efficiency = Some(Murphree::from_column_wide(0.6));
    let out = distillation_column(&setup).expect("the column converges");

    let temperatures = [
        373.15,
        336.17610149619395,
        303.1308245367222,
        302.1846100734558,
        296.56967983176895,
        253.15,
    ];
    let gas = [
        1.5935901648726694,
        0.6104268433574532,
        4.497779636973569,
        3.904812548909508,
        3.80082289353683,
        3.7264760720867867,
    ];
    let liquid = [
        3.764227964204178,
        5.889255034047075,
        4.906090702189604,
        1.3027391363431713,
        0.709771702362906,
        0.6057817557197408,
    ];
    assert_eq!(out.trays.len(), 6);
    for (i, tray) in out.trays.iter().enumerate() {
        relative(
            tray.temperature.value,
            temperatures[i],
            1.0e-5,
            "tray temperature",
        );
        relative(tray.gas_n, gas[i], 1.0e-5, "tray vapour");
        relative(tray.liquid_n, liquid[i], 1.0e-5, "tray liquid");
    }

    relative(out.distillate.n, 3.7264760720867867, 1.0e-6, "distillate");
    relative(
        out.distillate.z[0],
        0.9745893826201694,
        1.0e-6,
        "the distillate's methane",
    );
    relative(out.bottoms.n, 3.764227964204178, 1.0e-6, "bottoms");
    relative(
        out.bottoms.z[1],
        0.969829668983505,
        1.0e-6,
        "the bottoms' n-butane",
    );

    // The duties are enthalpy *differences* between two libraries whose absolute enthalpies
    // differ, so they are held to a watt rather than to a relative band - measured at `0.37` W
    // of `18651` and `1.02` W of `61658`.
    absolute(
        out.condenser_duty.value,
        -18650.975463115934,
        5.0,
        "condenser duty",
    );
    absolute(
        out.reboiler_duty.value,
        61658.439796998275,
        5.0,
        "reboiler duty",
    );
}

/// **A corrected column's products carry the feed exactly and its solve does not.**
///
/// `applyMurphreeCorrection` moves components between a stage's vapour and its liquid without
/// moving its moles, so what a stage hands up is not what its own flash produced and the interior
/// imbalances no longer cancel: on this row the tray terminals miss the feed's methane by
/// `5.6e-2` relative. `updateProductsFromExternalComponentBalance` is the class's own answer -
/// rescale the two products per component so they carry the feed - and it is on the run path, so
/// this port publishes the same record.
///
/// **`mass_residual` is the solve's own number and is not the class's `getLastMassResidual`**,
/// which compares a stage's inlets with the stage's own system and is therefore identically zero
/// on every column. The port's is a measurement of the closure, and a corrected column does not
/// close it; that is why it is reported here and gated only on a conserving column.
#[test]
fn a_corrected_columns_products_are_reconciled_against_its_feed() {
    let mut setup = binary_column(1.0e-6);
    setup.murphree_efficiency = Some(Murphree::from_column_wide(0.6));
    let out = distillation_column(&setup).expect("the column converges");

    let feed = binary_feed();
    assert!(
        out.mass_residual > 1.0e-2,
        "the solve's own closure is the measurement: {}",
        out.mass_residual
    );
    // The published pair carries the feed to the library's own digits.
    for c in 0..2 {
        let delivered = out.distillate.n * out.distillate.z[c] + out.bottoms.n * out.bottoms.z[c];
        relative(delivered, feed.n * feed.z[c], 1.0e-9, "the reconciled feed");
    }
    relative(
        out.distillate.n + out.bottoms.n,
        feed.n,
        1.0e-9,
        "the reconciled mole balance",
    );
}

/// **A high efficiency reaches a different fixed point here than in NeqSim, and that is
/// measured rather than excused.**
///
/// The capture holds `binary_murphree_0_85` beside the `0.6` row, and beside *it* the same state
/// at a gate of `1e-9` instead of `1e-6`: `binary_murphree_0_85_tight` takes 22 iterations
/// against 15 and moves every tray by less than `8e-10` K, so **NeqSim's `0.85` row is a
/// converged state and not a partial one**. This port finds a different one - nearer the ideal
/// profile, and self-consistent to `1e-13` K - and the two differ by up to `7` K on tray 3.
///
/// **The two implementations reach different fixed points at `0.85`, and the cause is NOT
/// established.** At `0.6` they agree to `1e-6` (the row above); at `0.85`, with both solvers
/// reporting convergence, they differ by up to `7` K on tray 3.
///
/// **`applyRelaxationFast` is not the reason, and that is measured rather than argued.** The
/// class relaxes the *streams* a stage is handed by its own relaxation factor, which this port
/// does not - so the obvious guess is that the path is what selects the fixed point. It was
/// implemented and measured on both rows: **`0.6` still matches, and `0.85` lands on the same
/// state to every printed digit** (`302.156771` on tray 3, with and without). A relaxation
/// cannot move a fixed point, and these two rows say it does not.
///
/// **The phase guards are refuted too.** The class's correction skips a stage whose system holds
/// fewer than two phases - and so does this port, whose `(Some, Some)` is the same statement -
/// so a stage skipped by one and corrected by the other would be the divergence. The capture's
/// `trayN_phases` rows say **every stage holds two phases on both rows**, with `phase0type=GAS`,
/// so neither implementation skips anything and the guard is not the difference either.
///
/// So the difference is in the *maps* and not in the paths, and it is not localised. What this
/// test holds is the *measurement* rather than a tolerance: the shape of the port's answer, and
/// the size of the disagreement - so that a later tranche which closes it moves this test
/// deliberately rather than by accident.
#[test]
fn a_high_murphree_efficiency_reaches_this_ports_own_fixed_point() {
    let mut setup = binary_column(1.0e-6);
    setup.murphree_efficiency = Some(Murphree::from_column_wide(0.85));
    let out = distillation_column(&setup).expect("the column converges");

    // The port's own answer: a profile between the ideal one and the `0.6` one.
    let temperatures = [
        373.15,
        336.1669805048,
        303.1060154799,
        302.1567710365,
        296.7533017195,
        253.15,
    ];
    for (i, tray) in out.trays.iter().enumerate() {
        absolute(
            tray.temperature.value,
            temperatures[i],
            1.0e-6,
            "tray temperature",
        );
    }
    // And the ends are still the pins, still untouched by the correction.
    absolute(out.trays[0].temperature.value, 373.15, 1.0e-9, "reboiler");
    absolute(out.trays[5].temperature.value, 253.15, 1.0e-9, "condenser");

    // **The disagreement is the measurement, so it is asserted as one**: NeqSim's tray 3 is
    // `313.37` K against this port's `302.16`, and a port that quietly drifted onto NeqSim's
    // state without the relaxation path would fail here rather than pass unnoticed.
    let disagreement = (out.trays[3].temperature.value - 313.36769810181676).abs();
    assert!(
        disagreement > 5.0,
        "tray 3 is {disagreement} K from the capture, and this row exists to hold that gap"
    );
}

/// **A gas side draw on a corrected stage is refused, and the class's own arithmetic is why.**
/// `getGasOutStream` returns the corrected vapour from the cache without applying
/// `gasSideDrawFraction`, while the draw's own stream is built from the uncorrected phase - so
/// the two outlets together carry more vapour than the stage made. The port refuses the pair
/// rather than reproducing a stage whose balance does not close.
#[test]
fn a_gas_draw_on_a_corrected_stage_is_refused() {
    let mut setup = binary_column(1.0e-6);
    setup.murphree_efficiency = Some(Murphree::from_column_wide(0.7));
    setup.gas_side_draw_fractions = Some(vec![0.0, 0.0, 0.0, 0.2, 0.0, 0.0]);
    let error = distillation_column(&setup).expect_err("the pair is refused");
    assert!(
        error.to_string().contains("getGasOutStream"),
        "the refusal names the class: {error}"
    );
}

/// The two ends draw nothing, so a fraction stated there is refused before the correction is
/// ever consulted - which is what keeps the two refusals from overlapping.
#[test]
fn an_ideal_stage_still_draws_its_gas() {
    let mut setup = binary_column(1.0e-6);
    setup.gas_side_draw_fractions = Some(vec![0.0, 0.0, 0.0, 0.2, 0.0, 0.0]);
    distillation_column(&setup).expect("an ideal column draws freely");
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
        murphree_efficiency: None,
        initial_state: None,
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
        side_draw_flows: Vec::new(),
        pumparound_returns: Vec::new(),
        pumparound_inlets: Vec::new(),
        pumparound_tolerance: None,
        pumparound_max_iterations: None,
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

/// **The side-draw *flow* specification: a draw's mass flow as a tear over the whole solve.**
///
/// Every other draw this port carries is a fraction of a stage's own phase, decided once.
/// `addSideDrawFlowSpecification` states the mass flow a draw must deliver and moves its
/// fraction until it does - so the fraction is a *controlled variable* and the column is solved
/// repeatedly.
///
/// **The oracle is a substitution, and the row it substitutes for is in the capture beside it.**
/// The class's own state is a one-tray column fed pure methane, and on it the vapour *is* the
/// whole feed: NeqSim publishes a bottom product of `0.0` mol/s whose enthalpy is
/// `-Infinity`, and this port declines to fabricate that - a `Stream`'s molar enthalpy is
/// finite by construction, which is `TrayOutcome`'s own documented judgement. The row is kept
/// as the evidence for that refusal. What is reproduced here is the same one-tray, no-ends,
/// gas-draw shape on the binary feed this model is already held to, which does have a liquid.
///
/// **The search is the class's own and it is legible.** The class seeds an uncontrolled draw at
/// `0.05`, so the first candidate is `0.05` of the tray's vapour; the second is the
/// multiplicative proposal `0.05 x target/actual`; and on this state that is the answer - **two
/// candidates, the flow on target to the digit, and no rejection at all**. The fraction the
/// search lands on is `0.07362063713783379` and the flow is `25.0` kg/hr exactly.
#[test]
fn a_specified_side_draw_flow_is_torn_until_the_draw_delivers_it() {
    let mut setup = ColumnSetup {
        feed: binary_feed(),
        feed_stage: 0,
        number_of_stages: 1,
        has_reboiler: false,
        has_condenser: false,
        top_pressure: pascals(20.0e5),
        bottom_pressure: pascals(20.0e5),
        condenser_temperature: None,
        reboiler_temperature: None,
        temperature_tolerance: 1.0e-6,
        max_iterations: 200,
        murphree_efficiency: None,
        initial_state: None,
        top_specification: None,
        bottom_specification: None,
        top_feed: None,
        tray_temperatures: None,
        solver_type: SolverType::DirectSubstitution,
        reactive: azoth_process::kernels::ReactiveSection::None,
        gas_side_draw_fractions: None,
        liquid_side_draw_fractions: None,
        pumparound_fractions: None,
        side_draw_flows: Vec::new(),
        pumparound_returns: Vec::new(),
        pumparound_inlets: Vec::new(),
        pumparound_tolerance: None,
        pumparound_max_iterations: None,
    };
    setup.side_draw_flows = vec![SideDrawFlow {
        tray: 0,
        phase: SideDrawPhase::Gas,
        target: 25.0 / 3600.0,
        tolerance: 1.0e-5,
        max_iterations: 12,
    }];
    let out = distillation_column(&setup).expect("the tear finds the fraction");

    let tear = out
        .tear
        .as_ref()
        .expect("a flow specification leaves a trace");
    assert_eq!(
        tear.iterations, 2,
        "the seed and the multiplicative proposal"
    );
    assert!(tear.converged, "the residual is {} ", tear.residual);
    absolute(tear.residual, 0.0, 1.0e-15, "the residual");
    assert_eq!(tear.rejected_candidates, 0);
    assert_eq!(tear.rollbacks, 0);
    assert_eq!(tear.inner_iterations, 2, "one inner solve per candidate");
    relative(
        tear.fraction,
        0.07362063713783379,
        1.0e-12,
        "the fraction the search lands on",
    );
    relative(
        tear.actual_flow * 3600.0,
        25.0,
        1.0e-15,
        "the draw's own flow, kg/hr",
    );
    assert!(
        tear.history.starts_with("#1 fraction=0.05"),
        "the class seeds an uncontrolled draw at five per cent: {}",
        tear.history
    );

    // And the column the tear published is a two-phase one-tray column, as the capture's is.
    relative(
        out.trays[0].gas_n,
        3.7940693256162255,
        1.0e-6,
        "the tray's vapour",
    );
    relative(
        out.trays[0].liquid_n,
        3.395114823364233,
        1.0e-6,
        "the tray's liquid",
    );
    let draw = out.gas_side_draws[0]
        .as_ref()
        .expect("the tray withdraws vapour");
    relative(
        draw.mass_flow().expect("a mass flow") * 3600.0,
        25.0,
        1.0e-6,
        "the draw",
    );
}

/// **Two of the tear's three neighbours are refused rather than reduced**, each by name.
///
/// `solveWithColumnTearVariables` is the class's *coordinated* problem - several tear variables
/// updated together under one convergence test - and this port carries the independent
/// single-variable search. So a second specification, and a specification beside a pumparound,
/// are refused with that reason; and a specification on an end is refused because the ends here
/// are `column::reboiler` and `column::condenser` rather than stages.
#[test]
fn a_coordinated_side_draw_tear_is_refused_by_name() {
    let single = SideDrawFlow {
        tray: 0,
        phase: SideDrawPhase::Gas,
        target: 25.0 / 3600.0,
        tolerance: 1.0e-5,
        max_iterations: 12,
    };
    let mut setup = ColumnSetup {
        feed: binary_feed(),
        feed_stage: 0,
        number_of_stages: 1,
        has_reboiler: false,
        has_condenser: false,
        top_pressure: pascals(20.0e5),
        bottom_pressure: pascals(20.0e5),
        condenser_temperature: None,
        reboiler_temperature: None,
        temperature_tolerance: 1.0e-6,
        max_iterations: 200,
        murphree_efficiency: None,
        initial_state: None,
        top_specification: None,
        bottom_specification: None,
        top_feed: None,
        tray_temperatures: None,
        solver_type: SolverType::DirectSubstitution,
        reactive: azoth_process::kernels::ReactiveSection::None,
        gas_side_draw_fractions: None,
        liquid_side_draw_fractions: None,
        pumparound_fractions: None,
        side_draw_flows: vec![single.clone(), single.clone()],
        pumparound_returns: Vec::new(),
        pumparound_inlets: Vec::new(),
        pumparound_tolerance: None,
        pumparound_max_iterations: None,
    };
    let error = distillation_column(&setup).expect_err("several tear variables are coordinated");
    assert!(
        error.to_string().contains("solveWithColumnTearVariables"),
        "{error}"
    );

    setup.side_draw_flows = vec![single];
    setup.pumparound_fractions = Some(vec![0.05]);
    let error = distillation_column(&setup).expect_err("a pumparound is coordinated too");
    assert!(error.to_string().contains("pumparound"), "{error}");
}

/// **The pumparound's return: a liquid draw that comes back to another tray.**
///
/// The draw alone this port already carried - `pumparound_fractions` is one of the three
/// fraction vectors - and what it did not carry is the *return*: `addLiquidPumparound("PA", 1,
/// 3, 0.10, 5.0)` withdraws a tenth of tray 1's liquid, cools it by five kelvin, and feeds it in
/// at tray 3. The class converges it as an ordinary fixed point on the return's own flow, and
/// NeqSim's captured row is this port's own binary column with that one mechanism added, so the
/// difference from `a_column_solves_the_captured_binary_profile` is the recycle and nothing else.
///
/// **The mechanism is legible in what it produces.** The return is the draw at exactly
/// `T - 5 K`, its flow settles to `2.8e-5` relative change against the class's own `1e-4`
/// tolerance, the cooler takes `-438.78` W, and the recycle moves the profile: tray 1 lands at
/// `334.58` K against the ideal column's `336.15`, and tray 3 at `304.82` against `302.12`.
#[test]
fn a_pumparound_return_is_recycled_until_its_flow_settles() {
    use azoth_process::column::pumparound::PumparoundReturn;

    let mut setup = binary_column(1.0e-6);
    setup.pumparound_returns = vec![PumparoundReturn {
        draw_tray: 1,
        return_tray: 3,
        fraction: 0.10,
        temperature_drop: 5.0,
    }];
    setup.pumparound_tolerance = Some(1.0e-4);
    setup.pumparound_max_iterations = Some(12);
    let out = distillation_column(&setup).expect("the recycle converges");

    let pumparound = out
        .pumparound
        .as_ref()
        .expect("a pumparound with a return leaves a trace");
    assert!(
        pumparound.converged,
        "relative change {}",
        pumparound.relative_change
    );
    assert!(
        pumparound.relative_change <= 1.0e-4,
        "the class's own tolerance: {}",
        pumparound.relative_change
    );

    // **The return is the draw at `T - drop`**, which is `updateReturnStream`'s own clone,
    // restate and re-flash: the captured draw is `0.571186834474117` mol/s at
    // `334.5849285970497` K and the return is the same flow at `329.5849285970497`.
    relative(
        pumparound.return_n[0],
        0.571186834474117,
        1.0e-4,
        "the return's flow",
    );
    absolute(
        pumparound.duty[0],
        -438.78378433523903,
        1.0,
        "the cooler's duty",
    );

    // And the recycle moved the column, which is what makes it worth carrying: measured, NeqSim
    // reports `334.5849285970497` on tray 1 and `304.82237776677346` on tray 3.
    absolute(
        out.trays[1].temperature.value,
        334.5849285970497,
        1.0e-3,
        "the draw tray's temperature",
    );
    absolute(
        out.trays[3].temperature.value,
        304.82237776677346,
        1.0e-3,
        "the return tray's temperature",
    );
    // **The draw is the tenth of the tray's own liquid *phase*, not of its outlet**, which is
    // the half that already existed: `getLiquidOutStream` is the liquid scaled by
    // `1 - liquid - pumparound`, so the outlet is nine tenths of the phase and the draw is the
    // tenth that left it. NeqSim's own row is `0.571186834474117` mol/s.
    let draw = out.pumparounds[1]
        .as_ref()
        .expect("tray 1 draws a pumparound");
    relative(
        draw.n,
        0.10 * (draw.n + out.trays[1].liquid_n),
        1.0e-9,
        "the draw's share",
    );
    relative(draw.n, 0.571186834474117, 1.0e-4, "the draw's own flow");
}

/// **The per-stage override resolves where the class resolves it, and `NaN` is the
/// fall-through.** `validation/neqsim/captures/process_column_efficiency.tsv` runs this same
/// binary column three times: at `0.6` column-wide, at `0.6` with stage 3 overridden to `0.85`,
/// and at `0.6` with the same override spelled as an array of `NaN`s elsewhere. NeqSim's second
/// and third blocks are **identical apart from their label**, so the array *is* the override and
/// `NaN` *is* `getEffectiveMurphreeEfficiency`'s fall-through.
///
/// **The first half is exact and is the resolution.** An all-`NaN` vector over a column-wide
/// `0.6` has to give the column-wide row's own profile to the last digit - every stage falls
/// through, so the two states are one - and that is what separates a resolution from a
/// per-stage clamp or an off-by-one.
///
/// **The second half is a divergence, and it is the same one the `0.85` column-wide row has.**
/// NeqSim's overridden state sits at `321.44` K on tray 3; this port reaches `302.19`, which is
/// the *column-wide* state to `1.3e-3` K. The override is applied - the kernel resolves stage 3
/// to `0.85` and `corrects` admits it - so this is not a dropped parameter but a second fixed
/// point the port's map does not reach, which is item 5's open question rather than this
/// item's. It is asserted as a *size* here for that reason, the same way
/// [`a_high_murphree_efficiency_reaches_this_ports_own_fixed_point`] holds its 7 K.
#[test]
fn a_per_stage_override_resolves_where_the_class_resolves_it() {
    let mut plain = binary_column(1.0e-6);
    plain.murphree_efficiency = Some(Murphree::from_column_wide(0.6));
    let column_wide = distillation_column(&plain).expect("the column converges");

    // Every stage falls through, so this is the same state spelled the declared way.
    let mut all_nan = binary_column(1.0e-6);
    all_nan.murphree_efficiency = Some(Murphree {
        column_wide: 0.6,
        per_stage: Some(vec![f64::NAN; 6]),
    });
    let fallen_through = distillation_column(&all_nan).expect("the column converges");
    for (index, (a, b)) in column_wide
        .trays
        .iter()
        .zip(&fallen_through.trays)
        .enumerate()
    {
        absolute(
            a.temperature.value,
            b.temperature.value,
            0.0,
            &format!("tray {index} of an all-NaN vector and the column-wide value"),
        );
    }

    // The column-wide state, against the capture's own row.
    absolute(
        column_wide.trays[3].temperature.value,
        302.1846100734558,
        1.0e-3,
        "the column-wide tray 3",
    );

    // And the override, which the port applies and which lands near that same state rather than
    // on NeqSim's. The gap is the measurement.
    let mut overridden = binary_column(1.0e-6);
    overridden.murphree_efficiency = Some(Murphree {
        column_wide: 0.6,
        per_stage: Some(vec![f64::NAN, f64::NAN, f64::NAN, 0.85, f64::NAN, f64::NAN]),
    });
    let out = distillation_column(&overridden).expect("the column converges");
    let disagreement = (out.trays[3].temperature.value - 321.4369463103764).abs();
    assert!(
        disagreement > 19.0,
        "tray 3 is {disagreement} K from the capture, and this row exists to hold that gap"
    );
}

/// **A warm start from a solved column lands on that column's own fixed point**, which is what
/// a state *is*: `DistillationColumn.initializeTrayStateFromColumn`, installed instead of the
/// cold seed.
///
/// **One pass, because the state is the answer.** The gate is the mean tray-temperature change
/// between passes, and a column started at its own endpoint has that at zero - so a warm start
/// that took more than one pass would mean the install was lossy.
///
/// **And it is measured on a corrected column too**, which is the half that matters: the
/// Murphree correction reads the neighbour below's *raw flash*, and this tranche's plan asserted
/// a state had to carry that pair. `Network::run` rewrites `equilibrium[i]` from the flash
/// before the correction reads `equilibrium[i - 1]`, and the sweep walks upward - so the pair is
/// this pass's by the time it is read, and a state without one is faithful. These two rows are
/// that measurement: `E = 0.6` corrects four of the six trays, and the warm start still lands on
/// the cold path's temperatures to `1e-9` K.
#[test]
fn a_warm_start_reproduces_the_cold_paths_tray_temperatures() {
    for efficiency in [None, Some(Murphree::from_column_wide(0.6))] {
        let mut setup = binary_column(1.0e-6);
        setup.murphree_efficiency = efficiency.clone();
        let cold = distillation_column(&setup).expect("the cold solve converges");

        let state =
            ColumnState::of(&cold, &setup.feed.components).expect("the state rebuilds from it");
        let mut warm_setup = binary_column(1.0e-6);
        warm_setup.murphree_efficiency = efficiency;
        warm_setup.initial_state = Some(state);
        let warm = distillation_column(&warm_setup).expect("the warm solve converges");

        assert_eq!(
            warm.iterations, 1,
            "a state at the fixed point is one pass, and this took {}",
            warm.iterations
        );
        for (index, (before, after)) in cold.trays.iter().zip(&warm.trays).enumerate() {
            // **One pass moves a tray by at most the cold solve's own gate times the tray
            // count**, and the reason for the factor is the gate's own shape: the port stops on
            // `temperature_residual`, which is the *mean* tray change, so a single tray may move
            // by up to N times it. Measured, tray 1 moves `1.02e-6` K against a `1e-6` gate on
            // six trays. That is the instrument's floor - six orders under the `7` K it exists
            // to resolve - and it belongs in the divergence test's own doc.
            absolute(
                before.temperature.value,
                after.temperature.value,
                setup.temperature_tolerance * cold.trays.len() as f64,
                &format!("tray {index} of the warm start"),
            );
            // **The flows carry a rebuild's own floor and the temperatures do not.** A
            // `ColumnState` is `Stream`s rebuilt from each tray's `(T, P, z)`, and the reboiler's
            // vapour comes back `3.7e-8` relative from where the cold path left it - measured,
            // and the reason this one is not held to `1e-9`. It is five orders under the
            // `7` K this instrument exists to resolve.
            relative(
                before.gas_n,
                after.gas_n,
                1.0e-6,
                &format!("tray {index}'s vapour"),
            );
        }
        // Measured at `3.5e-9` relative, and held at the flows' floor for the same reason.
        relative(
            cold.distillate.n,
            warm.distillate.n,
            1.0e-6,
            "the warm start's distillate",
        );
    }
}

/// One row of a capture's `key=value` lines, selected by its label.
fn capture_row(capture: &str, label: &str) -> std::collections::BTreeMap<String, String> {
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../validation/neqsim/captures")
        .join(capture);
    let text = std::fs::read_to_string(&path).expect("the capture is committed");
    let block = text
        .split("\n\n")
        .find(|block| block.lines().next() == Some(label))
        .unwrap_or_else(|| panic!("{capture} has no `{label}` row"));
    block
        .lines()
        .filter_map(|line| line.split_once('='))
        .map(|(name, value)| (name.to_string(), value.to_string()))
        .collect()
}

fn row_number(row: &std::collections::BTreeMap<String, String>, key: &str) -> f64 {
    row.get(key)
        .unwrap_or_else(|| panic!("the row has no `{key}`"))
        .parse()
        .unwrap_or_else(|_| panic!("`{key}` is not a number"))
}

fn row_composition(row: &std::collections::BTreeMap<String, String>, key: &str) -> Vec<f64> {
    row.get(key)
        .unwrap_or_else(|| panic!("the row has no `{key}`"))
        .split_whitespace()
        .map(|pair| {
            pair.split_once(':')
                .expect("a composition is `name:value`")
                .1
                .parse()
                .expect("a mole fraction is a number")
        })
        .collect()
}

/// **The divergence instrument.** NeqSim's `0.85` endpoint, installed as this port's
/// `ColumnState`, run for **one** pass of this port's map, and compared tray by tray against the
/// state it came from.
///
/// **Why one pass decides.** NeqSim's endpoint is a *converged* state - the `_tight` row takes
/// 22 passes at a `1e-9` gate moving no tray by more than `8e-10` K - so it is a fixed point of
/// NeqSim's map to `1e-7` K, and this port's answer differs from it by `7` K on tray 3. Seven
/// orders between the instrument's floor (about `1e-6` K, measured in
/// [`a_warm_start_reproduces_the_cold_paths_tray_temperatures`]) and the signal is what makes a
/// single pass a measurement rather than a guess.
///
/// **The state is NeqSim's reads fed through this port's own arithmetic**, never a formula
/// retyped here: each tray's `y_eq` is the capture's flashed vapour, its neighbour below's is
/// `y_in`, and the corrected stream the tray *hands up* is
/// [`azoth_process::column::murphree::correct_vapour`] applied to that pair. The class's caches
/// are invalidated by `finalizeTrayProperties` at the end of its solve, so the handed-up vapour
/// is not in the capture for the interior trays - recomputing it is the only faithful reading,
/// and it is the port's arithmetic rather than a second implementation.
///
/// **`tray0` and `tray5` are the two special cases, and both are the capture's own shape.** The
/// reboiler is never corrected, so its gas is the flash's; and `tray5_gas_n` is the *published
/// distillate* rather than the condenser's flash, because `finalizeSolve` writes the products
/// back onto the ends. The capture shows the second: `tray5_out_gas_z` is `0.96882`, the
/// distillate, where `tray5_eq_gas_z` is `0.96591`.
///
/// **Gauss-Seidel localises the first tray *affected*, not the tray at fault.** A tray reading a
/// neighbour that is already wrong moves because of that, so a reading of "tray 2 moved" names
/// the lowest tray whose neighbourhood is inconsistent, and the tray whose own equation differs
/// is at or below it.
fn divergence_row(label: &str) -> azoth_process::kernels::distillation_column::ColumnOutcome {
    use azoth_process::kernels::ReactiveSection;
    use azoth_process::kernels::distillation_column::SolverType;

    let row = capture_row("process_column_divergence.tsv", label);
    let components: Vec<String> = vec!["methane".into(), "n-butane".into()];
    let trays = row_number(&row, "tray_count") as usize;
    let efficiency = row_number(&row, "murphree_efficiency");

    let mut equilibrium_gas = Vec::with_capacity(trays);
    let mut equilibrium_liquid = Vec::with_capacity(trays);
    let mut temperatures = Vec::with_capacity(trays);
    for index in 0..trays {
        let at = |key: &str| row_number(&row, &format!("tray{index}_{key}"));
        let p = pascals(at("pressure_bara") * 1.0e5);
        let t = kelvins(at("temperature_K"));
        equilibrium_gas.push(
            Stream::from_pt(
                components.clone(),
                row_composition(&row, &format!("tray{index}_eq_gas_z")),
                at("gas_n"),
                p,
                t,
            )
            .expect("the tray's vapour resolves"),
        );
        equilibrium_liquid.push(
            Stream::from_pt(
                components.clone(),
                row_composition(&row, &format!("tray{index}_eq_liquid_z")),
                at("liquid_n"),
                p,
                t,
            )
            .expect("the tray's liquid resolves"),
        );
        temperatures.push(at("temperature_K"));
    }

    // **The vapour each tray hands up**, which is the corrected blend where the class corrects
    // and the flash where it does not - `corrects` is the port's own guard, so this cannot
    // disagree with the solve about which trays are corrected.
    let mut gas = Vec::with_capacity(trays);
    for index in 0..trays {
        // `corrects` is the port's own guard, so this cannot disagree with the solve about which
        // trays are corrected.
        if corrects(index, trays, true, efficiency) {
            gas.push(
                correct_vapour(
                    (
                        Some(&equilibrium_gas[index]),
                        Some(&equilibrium_liquid[index]),
                    ),
                    (
                        Some(&equilibrium_gas[index - 1]),
                        Some(&equilibrium_liquid[index - 1]),
                    ),
                    efficiency,
                )
                .expect("the corrected vapour resolves")
                .expect("a corrected stage hands up a stream"),
            );
        } else {
            gas.push(equilibrium_gas[index].clone());
        }
    }
    // And the condenser hands up the published distillate, which is what the class wrote there.
    let top = trays - 1;
    gas[top] = Stream::from_pt(
        components,
        row_composition(&row, "distillate_z"),
        row_number(&row, "distillate_n"),
        pascals(row_number(&row, "distillate_P") * 1.0e5),
        kelvins(row_number(&row, "distillate_T")),
    )
    .expect("the distillate resolves");

    let mut setup = binary_column(1.0e9);
    setup.murphree_efficiency = Some(Murphree::from_column_wide(efficiency));
    setup.initial_state = Some(ColumnState {
        gas,
        liquid: equilibrium_liquid,
        temperatures,
    });
    setup.solver_type = SolverType::DirectSubstitution;
    setup.reactive = ReactiveSection::None;
    setup.max_iterations = 5;
    distillation_column(&setup).expect("one pass of the map")
}

/// **The control first, and it is what makes the reading a measurement.**
///
/// `0.6` is the row this port reproduces to `1e-4` K, so one pass from NeqSim's `0.6` endpoint
/// must move no tray by more than the instrument's floor. If it moves by `7` K the instrument
/// is broken - a lossy install, a wrong inlet rule, a wrong flash - and the `0.85` reading below
/// is void rather than interesting.
#[test]
fn the_instruments_control_settles_on_the_row_this_port_reproduces() {
    let before = capture_row("process_column_divergence.tsv", "divergence_0_6");
    let after = divergence_row("divergence_0_6");
    let mut worst = 0.0_f64;
    for (index, tray) in after.trays.iter().enumerate() {
        let disagreement = (tray.temperature.value
            - row_number(&before, &format!("tray{index}_temperature_K")))
        .abs();
        worst = worst.max(disagreement);
    }
    assert!(
        worst < 1.0e-4,
        "the control moved a tray by {worst} K, so the instrument is not faithful and the \
         `0.85` reading below says nothing"
    );
}

/// **The `0.85` reading, pre-registered.** Either NeqSim's endpoint is a fixed point of this
/// port's map - in which case the `7` K is a question of which of two fixed points each path
/// reached, and not of a different arithmetic - or the first tray to move is where the two maps
/// part.
#[test]
fn neqsims_0_85_endpoint_under_one_pass_of_this_ports_map() {
    let before = capture_row("process_column_divergence.tsv", "divergence_0_85");
    let trays = row_number(&before, "tray_count") as usize;
    let after = divergence_row("divergence_0_85");

    let mut moved: Vec<(usize, f64, f64)> = Vec::new();
    for (index, tray) in after.trays.iter().enumerate() {
        let was = row_number(&before, &format!("tray{index}_temperature_K"));
        let delta = tray.temperature.value - was;
        if delta.abs() > 1.0e-4 {
            moved.push((index, was, delta));
        }
    }
    // **The localisation is the result, so it is asserted and not printed.** Measured on the
    // `0.85` row: `{2: -2.729881334332674, 3: -9.834999320712939, 4: -10.189198031777039}` and
    // nothing on trays 0, 1 or 5. Trays 2 and 3 are held at `1e-3` K so a port that moved them
    // differently fails rather than reaching this point; the *sign* matters as much as the size,
    // because the port's own `0.85` state is **above** NeqSim's on tray 3 by `11.2` K and this
    // pass moves NeqSim's *down* by `9.83` - the two maps part in the same direction.
    let expected = [
        (2_usize, -2.729881334332674_f64),
        (3, -9.834999320712939),
        (4, -10.189198031777039),
    ];
    let seen: Vec<usize> = moved.iter().map(|(index, _, _)| *index).collect();
    assert_eq!(
        seen,
        vec![2, 3, 4],
        "the trays that move, and the three that do not: the ends are uncorrected and tray 1 \
         reads a reboiler that has not moved"
    );
    for (index, delta) in &expected {
        let (_, _, actual) = moved
            .iter()
            .find(|(moved, _, _)| moved == index)
            .expect("the tray moved");
        assert!(
            (actual - delta).abs() < 1.0e-3,
            "tray {index} moves {actual} K where the measurement says {delta}"
        );
    }
    // **And the two ends do not move at all**, which is the correction's own guard showing:
    // `corrects` is false on the reboiler and the condenser, so whatever the interior does, the
    // ends are the same two states on both sides.
    for index in [0, trays - 1] {
        let was = row_number(&before, &format!("tray{index}_temperature_K"));
        assert!(
            (after.trays[index].temperature.value - was).abs() < 1.0e-6,
            "tray {index} is an end and must be untouched"
        );
    }
    // **The reading, printed rather than asserted away.** Which trays move, and by how much, is
    // the experiment's whole output; a bare pass/fail here would throw it out.
    for (index, was, delta) in &moved {
        println!("tray {index}: {was} K moves {delta} K");
        // **The classification the plan pre-registered**, on the three quantities the outcome
        // publishes: the flashed liquid (which the correction never touches), the vapour the
        // tray hands up, and the flash's two traffic rates.
        let tray = &after.trays[*index];
        let eq_liquid = row_composition(&before, &format!("tray{index}_eq_liquid_z"));
        let eq_gas = row_composition(&before, &format!("tray{index}_eq_gas_z"));
        println!(
            "  liquid_z  was {eq_liquid:?}\n            is {:?}",
            tray.liquid_z
        );
        println!(
            "  handed-up gas_z was {eq_gas:?}\n            is {:?}",
            tray.gas_z
        );
        println!(
            "  gas_n {:.6} -> {:.6}, liquid_n {:.6} -> {:.6}",
            row_number(&before, &format!("tray{index}_gas_n")),
            tray.gas_n,
            row_number(&before, &format!("tray{index}_liquid_n")),
            tray.liquid_n
        );
    }
    println!("{} of 6 trays moved by more than 1e-4 K", moved.len());
}
