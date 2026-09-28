//! The rate-based packed column's segment model, against
//! `validation/neqsim/captures/process_rate_based_packed_column.tsv`.
//!
//! **The capture is SRK and PR, and both are on purpose.** `RateBasedPackedColumnTest` runs
//! `SystemSrkEos`, and this library's process layer resolves PR - so the probe runs every state
//! on both cubics and the pair measures what the cubic moved. The PR rows are the ones this
//! port is held to; the `_srk` row and the two TEG rows are evidence the probe prints and this
//! file does not read.

use azoth_core::units::{kelvins, pascals};
use azoth_process::Stream;
use azoth_process::segment::{ProfileSettings, SnapshotSettings, solve_fixed_point_profile};

/// `feedMass`: a feed stated as a mass flow, which is how the class's own tests state theirs.
fn feed(
    names: &[&str],
    z: &[f64],
    temperature: f64,
    pressure_bara: f64,
    kg_per_hour: f64,
) -> Stream {
    let components: Vec<String> = names.iter().map(|name| (*name).to_string()).collect();
    let mut probe = Stream::from_pt(
        components.clone(),
        z.to_vec(),
        1.0,
        pascals(pressure_bara * 1.0e5),
        kelvins(temperature),
    )
    .expect("the feed resolves");
    let molar_mass = probe.molar_mass().expect("a molar mass");
    probe.n = (kg_per_hour / 3600.0) / molar_mass.value;
    Stream::from_pt(
        components,
        z.to_vec(),
        probe.n,
        pascals(pressure_bara * 1.0e5),
        kelvins(temperature),
    )
    .expect("the feed resolves at its own flow")
}

/// `configuredColumn`'s own state: the class's settings, on PR.
fn configured(height: f64) -> Stream {
    let _ = height;
    feed(&["methane", "CO2"], &[0.9, 0.1], 313.15, 50.0, 1000.0)
}

fn lean(height: f64) -> Stream {
    let _ = height;
    feed(&["water", "CO2"], &[1.0, 0.0], 303.15, 50.0, 2000.0)
}

fn settings(
    height: f64,
    heat_none: bool,
    billet_schultes: bool,
) -> (ProfileSettings, SnapshotSettings) {
    (
        ProfileSettings {
            packed_height: height,
            number_of_segments: 4,
            tolerance: 1.0e-9,
            max_iterations: 20,
        },
        SnapshotSettings {
            column_diameter: 1.0,
            packed_height: height,
            packing: "Pall-Ring-50".to_string(),
            heat_transfer_none: heat_none,
            mass_transfer_correction: 3.0,
            heat_transfer_correction: 1.0,
            billet_schultes,
        },
    )
}

/// **The class's own claim, and it is a direction rather than a digit.** `testAbsorbsCarbonDioxideFromGas`
/// asserts that the gas outlet's CO2 fraction falls, that the transfer total is positive, that
/// there are four segments and that the component balance closes to `1e-5` mol/s.
#[test]
fn the_absorber_absorbs_carbon_dioxide() {
    let (profile, snapshot) = settings(6.0, false, false);
    let gas_in = configured(6.0);
    let liquid_in = lean(6.0);
    let outcome = solve_fixed_point_profile(
        &gas_in,
        &liquid_in,
        &profile,
        &snapshot,
        &["CO2".to_string()],
        true,
    )
    .expect("the profile solves");

    let co2 = outcome
        .gas_outlet
        .components
        .iter()
        .position(|name| name == "CO2")
        .expect("the gas carries CO2");
    assert!(
        outcome.gas_outlet.z[co2] < gas_in.z[co2],
        "CO2 falls in the gas outlet: {} against {}",
        outcome.gas_outlet.z[co2],
        gas_in.z[co2]
    );
    let total = outcome
        .component_transfer_totals
        .iter()
        .find(|(name, _)| name == "CO2")
        .expect("CO2 transfers")
        .1;
    assert!(total > 0.0, "positive CO2 transfer is absorption: {total}");
    assert_eq!(outcome.segments.len(), 4);

    // The component balance, which is the class's own strongest assertion on this state.
    let before = co2_moles(&gas_in) + co2_moles(&liquid_in);
    let after = co2_moles(&outcome.gas_outlet) + co2_moles(&outcome.liquid_outlet);
    assert!(
        (before - after).abs() < 1.0e-5,
        "the CO2 balance closes: {before} in, {after} out"
    );
}

fn co2_moles(stream: &Stream) -> f64 {
    stream
        .components
        .iter()
        .position(|name| name == "CO2")
        .map_or(0.0, |index| stream.z[index] * stream.n)
}

/// `testZeroPackedHeightGivesNoTransfer`, and both halves are **exact**: a bed of no height
/// short-circuits the gate on the first pass and the transfer block is behind the class's own
/// `if (segmentHeight > 0.0)`.
#[test]
fn a_bed_of_no_height_transfers_nothing() {
    let (profile, snapshot) = settings(0.0, false, false);
    let outcome = solve_fixed_point_profile(
        &configured(0.0),
        &lean(0.0),
        &profile,
        &snapshot,
        &["CO2".to_string()],
        true,
    )
    .expect("the profile solves");
    assert_eq!(outcome.total_absolute_molar_transfer, 0.0);
    assert!(outcome.component_transfer_totals.is_empty());
    assert!(outcome.converged, "a bed of no height is a converged solve");
    assert_eq!(outcome.iterations, 1);
    // **And the profile is still four segments**, at a height of zero: `numberOfSegments` is
    // untouched by the short circuit, so a port that expected none would be wrong.
    assert_eq!(outcome.segments.len(), 4);
    assert!(
        outcome
            .segments
            .iter()
            .all(|segment| segment.height_from_bottom == 0.0)
    );
}

/// `testCanDisableExplicitHeatTransfer`, and the zeros are exact.
#[test]
fn disabling_heat_transfer_is_exactly_zero() {
    let (profile, snapshot) = settings(6.0, true, false);
    let outcome = solve_fixed_point_profile(
        &configured(6.0),
        &lean(6.0),
        &profile,
        &snapshot,
        &["CO2".to_string()],
        true,
    )
    .expect("the profile solves");
    for segment in &outcome.segments {
        assert_eq!(segment.overall_heat_transfer_coefficient, 0.0);
        assert_eq!(segment.heat_transfer_rate, 0.0);
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

/// **The PR row of the class's own absorber state, reproduced.** The capture is the probe's
/// own output on PR, and every number here is read from its `co2_water_absorber` block.
///
/// The convergence pair is the sharpest of them: the port stops after **fifteen** passes at a
/// residual of `6.02e-10` mol/s, and NeqSim's own row says `15` and `6.009670053264138e-10`.
/// That is the loop, the segment step, the interface flash, the hydraulics and the two-film
/// transfer agreeing at once - a wrong index anywhere in that chain moves the state by orders
/// of magnitude rather than by digits.
#[test]
fn the_classs_absorber_state_is_reproduced_on_pr() {
    let (profile, snapshot) = settings(6.0, false, false);
    let outcome = solve_fixed_point_profile(
        &configured(6.0),
        &lean(6.0),
        &profile,
        &snapshot,
        &["CO2".to_string()],
        true,
    )
    .expect("the profile solves");

    assert!(outcome.converged, "the gate is met");
    assert_eq!(outcome.iterations, 15, "NeqSim takes fifteen passes");
    relative(
        outcome.convergence_residual,
        6.009670053264138e-10,
        1.0e-4,
        "residual",
    );
    relative(
        outcome.total_absolute_molar_transfer,
        0.021771811808519004,
        1.0e-6,
        "total absolute transfer",
    );
    let co2_total = outcome
        .component_transfer_totals
        .iter()
        .find(|(name, _)| name == "CO2")
        .expect("CO2 transfers")
        .1;
    relative(co2_total, 0.0034653825611094258, 1.0e-5, "CO2 total");

    relative(
        outcome.gas_outlet.n,
        14.74081280540249,
        1.0e-8,
        "gas outlet flow",
    );
    relative(
        outcome.gas_outlet.z[1],
        0.09978842114433098,
        1.0e-8,
        "gas outlet CO2",
    );
    relative(
        outcome.liquid_outlet.n,
        30.841964163798856,
        1.0e-8,
        "liquid outlet flow",
    );
    relative(
        outcome.liquid_outlet.z[1],
        0.00011235931347143501,
        1.0e-4,
        "liquid outlet CO2",
    );

    // The bottom segment, which is the one the two-film balance is won or lost in.
    let first = &outcome.segments[0];
    relative(
        first.gas_temperature,
        308.7440274651043,
        1.0e-6,
        "segment 1 gas T",
    );
    relative(
        first.liquid_temperature,
        305.43212846307154,
        1.0e-5,
        "segment 1 liquid T",
    );
    relative(first.k_ga, 135.79988356254145, 1.0e-6, "segment 1 kGa");
    relative(first.k_la, 0.009432860559780317, 1.0e-6, "segment 1 kLa");
    relative(
        first.wetted_area,
        57.01800105586562,
        1.0e-5,
        "segment 1 wetted area",
    );
    relative(
        first.percent_flood,
        3.404657951671574,
        1.0e-4,
        "segment 1 flood",
    );
    relative(
        first.pressure_drop_per_meter,
        0.007530567871836917,
        1.0e-6,
        "segment 1 pressure drop",
    );
    relative(
        first.heat_transfer_rate,
        2731.7582800165587,
        1.0e-3,
        "segment 1 heat rate",
    );
    relative(
        first.interface_temperature,
        309.3469076998667,
        1.0e-5,
        "segment 1 interface T",
    );
    relative(
        first.net_molar_transfer,
        -0.0018659748797604063,
        1.0e-5,
        "segment 1 net transfer",
    );

    // **The profile is not monotone, and that is the measurement**: the three lower segments
    // strip CO2 back out of the liquid and the top one absorbs it, so a port that summed a
    // single direction would land somewhere else entirely.
    let nets: Vec<f64> = outcome
        .segments
        .iter()
        .map(|s| s.net_molar_transfer)
        .collect();
    assert!(
        nets[0] < 0.0 && nets[1] < 0.0 && nets[2] < 0.0 && nets[3] > 0.0,
        "{nets:?}"
    );
}

/// **The snapshot's reference diffusivity is the class's constant, and its pair matrix is
/// real.** `averageDiffusivity` averages a vector NeqSim never populates, so the reference is
/// `1.5e-5` m²/s for the gas and `1.5e-9` for the liquid - while the *pair* coefficient the
/// film model scales against it is `9.457085848059925e-7`, which is the value
/// `eos.phase_transport`'s own capture holds for methane/CO2 at 313 K and 50 bar.
#[test]
fn the_reference_diffusivity_is_the_classs_constant_and_the_matrix_is_real() {
    use azoth_process::segment::calculate_transport_snapshot;
    let (_, snapshot) = settings(6.0, false, false);
    let (snap, gas_phase, _liquid_phase) =
        calculate_transport_snapshot(&configured(6.0), &lean(6.0), 1.5, &snapshot)
            .expect("the snapshot");

    assert_eq!(snap.gas_diffusivity, 1.5e-5);
    assert_eq!(snap.liquid_diffusivity, 1.5e-9);
    relative(
        gas_phase.transport.d_binary[0][1],
        9.457085848059925e-7,
        1.0e-9,
        "the gas pair diffusivity",
    );
    assert!(
        snap.fallbacks.diffusivity,
        "the class's constant stood in, and the record says so"
    );
    // The two properties the class *does* answer, against the same capture.
    relative(snap.k_ga, 135.79988356254145, 1.0e-6, "kGa");
    relative(
        snap.overall_heat_transfer_coefficient,
        483220.7509944843,
        1.0e-4,
        "the overall heat coefficient",
    );
}

/// **The union transfer list, which is the model's default, refuses - and it refused before
/// this tranche too, for a different reason.**
///
/// The union of the two inlets' components always names the solvent, which the gas does not
/// carry, so this is the one state that reaches the class's `componentIndex(phase, component)
/// < 0` branch. Measured on both sides of that fix: without it the profile runs and does
/// **not converge** (`union.converged` is false at the iteration cap, which is what the wrong
/// gas-side film does to the solvent's transfer); with it, the solvent's film is the base
/// coefficient and the first pass lands the liquid's composition at `1 + 2.9e-9`, which
/// `crate::stream::Stream` refuses - "the composition sums to 1.0000000028839944, not to
/// one", which is the guard that keeps a *caller's* error visible and is now being asked about
/// the port's own renormalisation.
///
/// Both are findings about the union list rather than about the film index, and neither is
/// this step's to fix. What is pinned here is that the state does not quietly work: the
/// branch is exercised, and the refusal names the composition rather than a fabricated
/// coefficient. The film index itself is pinned in
/// `segment::film::tests::a_component_the_phase_does_not_carry_is_the_base_coefficient`.
#[test]
fn the_union_transfer_list_names_a_component_the_gas_does_not_carry() {
    let (profile, snapshot) = settings(6.0, false, false);
    let error =
        solve_fixed_point_profile(&configured(6.0), &lean(6.0), &profile, &snapshot, &[], true)
            .expect_err("the union list does not solve in this port");
    let message = format!("{error}");
    assert!(
        message.contains("composition sums to"),
        "and it refuses for a reason of its own: {message}"
    );
}

/// **`BILLET_SCHULTES_1999` is two constants of the packing and not a correlation, and the
/// capture's own pair measures it.** `process_rate_based_billet.tsv` runs the class's
/// configured absorber state twice on PR, once at each correlation; every number here is read
/// from those two blocks.
///
/// `Pall-Ring-50` resolves to the file's plastic row, `cp = 0.698` and `ch = 2.725`, so the
/// multipliers are `ch/6/0.4 = 1.1354166...` on `kGa` and `cp = 0.698` on `kLa` - **neither
/// floored at the class's `0.1`**, which is why this state does not exercise the clamp. What it
/// does pin is that the multiplier is applied *after* `finiteNonNegative` and *before*
/// `massTransferCorrectionFactor`, and that both film coefficients carry it into the heat
/// coefficients rather than only into the transfer.
#[test]
fn the_billet_schultes_multiplier_is_the_packings_own_two_constants() {
    let (profile, onda) = settings(6.0, false, false);
    let (_, billet) = settings(6.0, false, true);
    let solve = |snapshot: &SnapshotSettings| {
        solve_fixed_point_profile(
            &configured(6.0),
            &lean(6.0),
            &profile,
            snapshot,
            &["CO2".to_string()],
            true,
        )
        .expect("the profile solves")
    };
    let onda = solve(&onda);
    let billet = solve(&billet);

    // The bottom segment, against the capture's own two readings.
    relative(
        onda.segments[0].k_ga,
        135.79988356254145,
        1.0e-6,
        "ONDA kGa",
    );
    relative(
        billet.segments[0].k_ga,
        154.18945112830227,
        1.0e-6,
        "billet kGa",
    );
    relative(
        onda.segments[0].k_la,
        0.009432860559780317,
        1.0e-6,
        "ONDA kLa",
    );
    relative(
        billet.segments[0].k_la,
        0.006583152243472645,
        1.0e-6,
        "billet kLa",
    );

    // **And the pair's ratio is the multiplier, exactly on the gas side and only nearly on the
    // liquid one** - which is a property of the state and not of the arithmetic. The class
    // scales a *base* that both film coefficients are computed from, so the ratio of two runs
    // is the multiplier times the ratio of the two bases; the state moved between them, and the
    // two sides do not feel it equally. Measured on the capture's own pair: the gas ratio is
    // `1.13541666...` to the last digit, and the liquid ratio is `0.6978956` against `0.698`,
    // because the base `kLa` itself moved by `1.495e-4`.
    let (gas_multiplier, liquid_multiplier) = (
        billet.segments[0].k_ga / onda.segments[0].k_ga,
        billet.segments[0].k_la / onda.segments[0].k_la,
    );
    relative(
        gas_multiplier,
        2.725 / 6.0 / 0.4,
        1.0e-9,
        "the gas multiplier",
    );
    relative(
        liquid_multiplier,
        0.6978956385236719,
        1.0e-5,
        "the liquid multiplier, as the capture itself measures it",
    );
    assert!(
        (liquid_multiplier - 0.698).abs() < 2.0e-4,
        "and it is the packing's own `cp`: {liquid_multiplier} against 0.698"
    );

    // The state moves with it, and in the direction the capture does: less liquid-film capacity
    // retains less CO2, so the gas keeps more of it.
    assert!(
        billet.gas_outlet.z[1] > onda.gas_outlet.z[1],
        "CO2 in the gas outlet: {} billet against {} onda",
        billet.gas_outlet.z[1],
        onda.gas_outlet.z[1]
    );
    relative(
        billet.component_transfer_totals[0].1,
        0.0024291784706160805,
        1.0e-6,
        "billet CO2 total",
    );
}

/// One block of a capture, by its first line.
fn capture_row(capture: &str, label: &str) -> std::collections::BTreeMap<String, String> {
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../validation/neqsim/captures")
        .join(capture);
    let text = std::fs::read_to_string(&path).expect("the capture is committed");
    text.split("\n\n")
        .find(|block| block.lines().next() == Some(label))
        .unwrap_or_else(|| panic!("{capture} has no `{label}` row"))
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

/// **The two solvers this port refuses do not converge, and the capture is the measurement
/// rather than the class's own comment about itself.**
///
/// `process_rate_based_solvers.tsv` runs `configuredColumn`'s own state - the PR row
/// [`the_classs_absorber_state_is_reproduced_on_pr`] holds - five ways. The first is the pair
/// this port carries, and it is the control: it lands on the port's own fifteen passes and
/// `6.009670053264138e-10`, so what the other four show is a statement about the two solvers and
/// not about a mis-built state.
///
/// What they show is not that the branches are slow. **They are refused because neither
/// publishes a state anything can be held to**: `SIMULTANEOUS_RESIDUAL` takes the whole twenty
/// pass cap and its vapour leaves at `2.2e-29` K, and `EQUATION_ORIENTED` stalls after two Newton
/// passes with a residual six orders above its own gate and freezes the gas at `153` K beside a
/// `454` K liquid. In both equation-oriented rows the component balances close **exactly** and
/// the energy balance does not - the stall is the energy equation alone.
///
/// **The class's own test of that Newton is the reason this is a recorded non-port rather than a
/// tranche.** It caps the Newton at two iterations and the homotopy at one, and asserts
/// `Double.isFinite` on the residual rather than that the norm met the tolerance; its companion
/// caps them at one. A branch no upstream test holds to a residual has no oracle to port against.
#[test]
fn neither_refused_solver_converges_on_the_classs_own_state() {
    const CAPTURE: &str = "process_rate_based_solvers.tsv";

    // The control: this port's own pair, and the numbers the port's own test asserts.
    let control = capture_row(CAPTURE, "sequential_fixed_point");
    assert_eq!(
        control.get("segment_solver").map(String::as_str),
        Some("SEQUENTIAL_EXPLICIT"),
        "the control runs the pair this port carries"
    );
    assert_eq!(
        control.get("column_solver").map(String::as_str),
        Some("FIXED_POINT_PROFILE")
    );
    let (profile, snapshot) = settings(6.0, false, false);
    let outcome = solve_fixed_point_profile(
        &configured(6.0),
        &lean(6.0),
        &profile,
        &snapshot,
        &["CO2".to_string()],
        true,
    )
    .expect("the profile solves");
    assert_eq!(
        row_number(&control, "iterations") as usize,
        outcome.iterations,
        "the control row is this port's own pass count"
    );
    relative(
        outcome.convergence_residual,
        row_number(&control, "convergence_residual"),
        1.0e-4,
        "the control row's residual, against this port's",
    );

    // `SIMULTANEOUS_RESIDUAL`: the cap, nine orders above the gate, and a vapour that is not a
    // temperature.
    let simultaneous = capture_row(CAPTURE, "simultaneous_fixed_point");
    assert_eq!(
        row_number(&simultaneous, "convergence_tolerance_mol_per_s"),
        1.0e-9
    );
    assert_eq!(
        row_number(&simultaneous, "iterations"),
        row_number(&simultaneous, "max_iterations"),
        "the segment solve takes the whole cap"
    );
    assert!(
        row_number(&simultaneous, "convergence_residual") > 1.0,
        "and ends nine orders above its gate: {} mol/s",
        row_number(&simultaneous, "convergence_residual")
    );
    assert!(
        row_number(&simultaneous, "gas_out_T") < 1.0,
        "publishing a vapour at {} K",
        row_number(&simultaneous, "gas_out_T")
    );

    // `EQUATION_ORIENTED` alone: stalls in two passes, publishes a `153` K gas beside a `454` K
    // liquid, closes both component balances exactly and leaves the energy balance at `1.1e8` J.
    let equation_oriented = capture_row(CAPTURE, "sequential_equation_oriented");
    assert_eq!(
        row_number(&equation_oriented, "column_residual_tolerance"),
        1.0e-6
    );
    assert_eq!(row_number(&equation_oriented, "iterations"), 2.0);
    assert!(
        row_number(&equation_oriented, "column_residual_norm") > 1.0,
        "the norm is {} against a 1e-6 gate",
        row_number(&equation_oriented, "column_residual_norm")
    );
    assert_eq!(
        row_number(&equation_oriented, "gas_component_balance_residual"),
        0.0
    );
    assert_eq!(
        row_number(&equation_oriented, "liquid_component_balance_residual"),
        0.0
    );
    assert!(
        row_number(&equation_oriented, "column_energy_balance_residual") > 1.0e7,
        "and the energy balance is what stalled: {} J",
        row_number(&equation_oriented, "column_energy_balance_residual")
    );
    assert!(row_number(&equation_oriented, "gas_out_T") < 200.0);
    assert!(row_number(&equation_oriented, "liquid_out_T") > 400.0);

    // Both together, which is worse on both counts: a larger norm and a transfer two orders below
    // the control's.
    let both = capture_row(CAPTURE, "simultaneous_equation_oriented");
    assert!(
        row_number(&both, "column_residual_norm") > 1.0,
        "the pair is no better: {}",
        row_number(&both, "column_residual_norm")
    );
    relative(
        row_number(&both, "total_absolute_molar_transfer_mol_per_s"),
        1.6818984914629715e-4,
        1.0e-6,
        "the pair's transfer",
    );
    assert!(
        row_number(&both, "total_absolute_molar_transfer_mol_per_s")
            < row_number(&control, "total_absolute_molar_transfer_mol_per_s"),
        "the pair transfers less than the control, not more"
    );

    // **And the class's own equation-oriented test state, at its own four settings.** The two
    // settings that matter are the cap and the homotopy: the test supplies both, which is what
    // makes its assertions reachable without the Newton ever converging.
    let own = capture_row(CAPTURE, "equation_oriented_test_settings");
    assert_eq!(row_number(&own, "max_column_residual_iterations"), 2.0);
    assert_eq!(row_number(&own, "column_homotopy_steps"), 1.0);
    assert_eq!(row_number(&own, "column_residual_tolerance"), 1.0e-5);
    assert!(
        row_number(&own, "column_residual_norm") > 1.0,
        "and it does not converge there either: {} against 1e-5",
        row_number(&own, "column_residual_norm")
    );
}

/// **A refusal reads the declaration rather than carrying its own sentence**, and what it
/// reads is the row the spec's `[[unported]]` array declares.
///
/// The numbers behind the refusal - `1.76` mol/s against a `1e-9` gate, `65.9` against
/// `1e-6` - are the spec's own `description` of the two inputs, and the capture measures
/// them. What travels through the code is the key, the class and the capture, which is
/// what makes the two implementations comparable at all.
#[test]
fn both_refusals_read_the_row_their_spec_declares() {
    use azoth_process::kernels::rate_based_packed_column::{ColumnSolver, SegmentSolver};
    use azoth_process::unported;

    for (key, class) in [
        (
            "segment_solver=simultaneous_residual",
            "RateBasedPackedColumn.SegmentSolver.SIMULTANEOUS_RESIDUAL",
        ),
        (
            "column_solver=equation_oriented",
            "RateBasedPackedColumn.ColumnSolver.EQUATION_ORIENTED",
        ),
    ] {
        let row = unported::row(key).unwrap_or_else(|| panic!("`{key}` has no declared row"));
        assert_eq!(row.class, class, "the row for `{key}` names its class");
        assert_eq!(row.model, "process.rate_based_packed_column");
        assert!(
            row.capture.ends_with("process_rate_based_solvers.tsv"),
            "the row for `{key}` points at the capture"
        );
    }

    let segment = SegmentSolver::parse("simultaneous_residual")
        .expect_err("refused")
        .to_string();
    assert!(
        segment.contains("SIMULTANEOUS_RESIDUAL")
            && segment.contains("process_rate_based_solvers.tsv"),
        "the segment refusal names the class and the capture: {segment}"
    );

    let column = ColumnSolver::parse("equation_oriented")
        .expect_err("refused")
        .to_string();
    assert!(
        column.contains("EQUATION_ORIENTED") && column.contains("process_rate_based_solvers.tsv"),
        "the column refusal names the class and the capture: {column}"
    );
}
