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
