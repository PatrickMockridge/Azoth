//! The absorber and the stripper, against `validation/neqsim/captures/process_absorber.tsv`.

use azoth_core::units::{kelvins, meters, pascals};
use azoth_process::Stream;
use azoth_process::column::absorber_murphree::AbsorberMurphree;
use azoth_process::column::capacity::{fs_limits, gas_load_limits};
use azoth_process::column::mechanical::{
    DEFAULT_CONTACTOR_INTERNALS_TYPE, DEFAULT_MATERIAL_GRADE, DEFAULT_MAX_FLOODING_FACTOR,
    DEFAULT_MAX_OPERATION_PRESSURE_BARA, DEFAULT_TRAY_EFFICIENCY, DEFAULT_TRAY_TYPE,
    MechanicalGeometry, mechanical_design,
};
use azoth_process::column::murphree::Murphree;
use azoth_process::kernels::distillation_column::{
    ColumnSetup, SolverType, distillation_column, tray_streams,
};

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

/// **The lean-oil absorber of `AbsorptionColumnTest`, unpinned.**
///
/// The class's own case pins every tray to one temperature, and that is captured as evidence
/// rather than as an oracle: `SimpleTray.setOutletTemperature` on every stage makes
/// `solveSequential`'s own gate - the mean tray-temperature change - zero, so the solve stops
/// after one sweep with `iterations=1`, a temperature residual of `0.0`, and **no liquid
/// traffic at all**. Left unpinned, the same column converges in eight iterations to a real
/// counter-current state, and that is this row.
#[test]
fn the_lean_oil_absorber_solves_unpinned() {
    let out = distillation_column(&lean_oil()).expect("the absorber converges");

    // NeqSim takes 17 at this gate and the port 19: where a solve stops is its own.
    assert!(out.iterations <= 25, "in {} iterations", out.iterations);
    let temperatures = [
        299.10907966626354,
        299.7940316636134,
        300.49322271941304,
        301.26748639388506,
        301.7126561904438,
    ];
    let gas = [
        31.02657698982344,
        30.994044668451405,
        30.961128515891744,
        30.91441340721759,
        30.587796968266012,
    ];
    let liquid = [
        1.9280715291956925,
        2.1020464244421477,
        2.069511412099089,
        2.0365923915807302,
        1.9898748368771202,
    ];
    assert_eq!(out.trays.len(), 5, "no condenser and no reboiler");
    for (i, tray) in out.trays.iter().enumerate() {
        absolute(
            tray.temperature.value,
            temperatures[i],
            2.0e-3,
            "tray temperature",
        );
        relative(tray.gas_n, gas[i], 1.0e-5, "tray vapour");
        relative(tray.liquid_n, liquid[i], 1.0e-5, "tray liquid");
    }

    // The products are the class's own getters: the sweet gas overhead and the rich oil.
    relative(out.distillate.n, 30.58779316625974, 1.0e-5, "gas out");
    relative(out.bottoms.n, 1.928065918154764, 1.0e-5, "liquid out");
    relative(
        out.distillate.z[3],
        0.0070019486397466715,
        1.0e-4,
        "the gas out's n-butane",
    );
    relative(
        out.bottoms.z[5],
        0.756165997919113,
        1.0e-4,
        "the rich oil's n-heptane",
    );
    // **The separation the class's own test asserts**: the lean oil takes n-butane out of the
    // gas, and the heavier the component the more strongly it is taken.
    assert!(
        out.distillate.z[3] * out.distillate.n < 0.010 * 30.852602094576984,
        "the gas must lose n-butane"
    );
    assert!(
        out.bottoms.z[3] * out.bottoms.n > 0.0,
        "the rich oil must gain n-butane"
    );

    // **The two capacity families on the same state**, which the case file holds as numbers and
    // cannot hold as flags: `TestCase` carries a boolean for an *input* and not for an
    // expectation, so the two verdicts are asserted here, as the packed column's own flag is.
    // The K value is the capture's `lean_oil_absorber_default_diameter` row, and it is the
    // **fallback's**: `getPhase(0)` of the liquid outlet's own system is the vapour the bottom
    // tray carries, so `rho_l - rho_g` is negative and `1000.0` is substituted. Without the
    // fallback the same state reads `0.00826846105398939`, a factor of `1.235` - which is what
    // makes this assertion a measurement of the branch rather than of the arithmetic alone.
    let (bottom_vapour, _) = tray_streams(&out.trays[0], &out.distillate.components)
        .expect("the bottom tray's vapour rebuilds");
    let fs = fs_limits(&out.distillate, meters(1.0), 3.0).expect("the Fs family evaluates");
    let gas_load = gas_load_limits(&out.distillate, &bottom_vapour, meters(1.0), 0.15)
        .expect("the K evaluates");
    assert!(
        gas_load.gas_load_factor_within_design_limit,
        "the K of {:.9} is inside the class's own 0.15",
        gas_load.gas_load_factor
    );
    assert!(
        fs.fs_factor_within_design_limit,
        "the Fs of {:.9} is inside the absorber's own 3.0",
        fs.fs_factor
    );
    relative(
        gas_load.gas_load_factor,
        0.006694885232955643,
        1.0e-4,
        "the gas-load factor",
    );
    relative(fs.fs_factor, 0.21051436885536753, 1.0e-4, "the Fs factor");
}

/// **The hydrocarbon stripper of `StrippingColumnTest`, unpinned**, whose own case states
/// `MESH_RESIDUAL` - a rung this port refuses, so the row's own strategy is captured beside
/// the answer it produces and the port is held to the answer.
#[test]
fn the_hydrocarbon_stripper_solves_unpinned() {
    let out = distillation_column(&hydrocarbon_stripper()).expect("the stripper converges");

    assert!(out.iterations <= 20, "in {} iterations", out.iterations);
    let temperatures = [
        321.706608455782,
        325.7629838093825,
        329.5399755743931,
        333.18535921380345,
        337.4756767231461,
    ];
    for (i, tray) in out.trays.iter().enumerate() {
        absolute(
            tray.temperature.value,
            temperatures[i],
            2.0e-3,
            "tray temperature",
        );
    }
    relative(out.distillate.n, 3.086389912429422, 1.0e-5, "overhead gas");
    relative(out.bottoms.n, 2.4096710196150735, 1.0e-5, "lean liquid");
    // **The separation `StrippingColumnTest` asserts**, on the state the port reaches rather
    // than on the class's: the stripping gas gains both lights and the rich liquid loses them.
    assert!(
        out.distillate.z[1] > 0.0080,
        "the overhead gas must gain propane: {}",
        out.distillate.z[1]
    );
    assert!(
        out.distillate.z[3] > 0.0005,
        "the overhead gas must gain n-pentane: {}",
        out.distillate.z[3]
    );
    assert!(
        out.bottoms.z[1] < 0.08 && out.bottoms.z[3] < 0.15,
        "the lean liquid must lose both: propane {}, pentane {}",
        out.bottoms.z[1],
        out.bottoms.z[3]
    );
}

/// The lean-oil absorber: gas at stage 0, solvent at the top stage, and no ends.
fn lean_oil() -> ColumnSetup {
    ColumnSetup {
        feed: Stream::from_pt(
            vec![
                "methane".into(),
                "ethane".into(),
                "propane".into(),
                "n-butane".into(),
                "n-pentane".into(),
                "n-heptane".into(),
            ],
            vec![0.920, 0.040, 0.025, 0.010, 0.005, 0.0],
            30.852602094576984,
            pascals(15.0e5),
            kelvins(303.15),
        )
        .expect("the fluid resolves"),
        feed_stage: 0,
        number_of_stages: 5,
        has_reboiler: false,
        has_condenser: false,
        top_pressure: pascals(15.0e5),
        bottom_pressure: pascals(15.0e5),
        condenser_temperature: None,
        reboiler_temperature: None,
        temperature_tolerance: 1.0e-4,
        murphree_efficiency: None,
        absorber_murphree: None,
        initial_state: None,
        max_iterations: 80,
        top_specification: None,
        bottom_specification: None,
        top_feed: Some(
            Stream::from_pt(
                vec![
                    "methane".into(),
                    "ethane".into(),
                    "propane".into(),
                    "n-butane".into(),
                    "n-pentane".into(),
                    "n-heptane".into(),
                ],
                vec![0.0, 0.0, 0.0, 0.0, 0.0, 1.0],
                1.6632569898375,
                pascals(15.0e5),
                kelvins(293.15),
            )
            .expect("the fluid resolves"),
        ),
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

/// The hydrocarbon stripper: stripping gas at stage 0, rich liquid at the top, no ends.
fn hydrocarbon_stripper() -> ColumnSetup {
    ColumnSetup {
        feed: Stream::from_pt(
            vec![
                "methane".into(),
                "propane".into(),
                "n-butane".into(),
                "n-pentane".into(),
                "n-heptane".into(),
            ],
            vec![0.9900, 0.0080, 0.0015, 0.0005, 0.0],
            2.547079374624363,
            pascals(12.0e5),
            kelvins(343.15),
        )
        .expect("the fluid resolves"),
        feed_stage: 0,
        number_of_stages: 5,
        has_reboiler: false,
        has_condenser: false,
        top_pressure: pascals(12.0e5),
        bottom_pressure: pascals(12.0e5),
        condenser_temperature: None,
        reboiler_temperature: None,
        temperature_tolerance: 1.0e-4,
        murphree_efficiency: None,
        absorber_murphree: None,
        initial_state: None,
        max_iterations: 80,
        top_specification: None,
        bottom_specification: None,
        top_feed: Some(
            Stream::from_pt(
                vec![
                    "methane".into(),
                    "propane".into(),
                    "n-butane".into(),
                    "n-pentane".into(),
                    "n-heptane".into(),
                ],
                vec![0.02, 0.08, 0.12, 0.15, 0.63],
                2.9489815574232177,
                pascals(12.0e5),
                kelvins(343.15),
            )
            .expect("the fluid resolves"),
        ),
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

/// **The class's own isothermal case is refused by the port's own closure gate.**
///
/// Pinning every tray makes the base's gate - the mean tray-temperature change - zero, so the
/// solve stops after its first sweep. The state that sweep lands on has an energy closure of
/// `0.28` against the class's default gate of `1.6e-2`, so the port refuses it: the class's own
/// test accepts the same state only because it loosens that gate to `5e-2`, which is a settable
/// this port has not carried. The capture keeps the row as evidence.
#[test]
fn a_pinned_column_stops_after_one_sweep_and_is_refused() {
    let mut setup = lean_oil();
    setup.tray_temperatures = Some(vec![298.15; 5]);
    setup.temperature_tolerance = 1.0e-2;

    let error = distillation_column(&setup).expect_err("the closure gate refuses the state");
    assert!(
        matches!(error, azoth_core::AzothError::SolverNotConverged { .. }),
        "{error}"
    );
}

/// **The stripper's id, on the same state its case states**, which is what verifies the model
/// layer while the extension cannot be rebuilt: `process.stripping_column` is the absorber's
/// model under the class's own names, so this asserts the names *and* the numbers.
#[test]
fn the_stripping_column_id_reaches_the_absorber_model() {
    let setup = hydrocarbon_stripper();
    let rich = setup.top_feed.as_ref().expect("the stripper's rich liquid");
    let out = azoth_process::stripping_column(
        &setup.feed.components,
        &rich.components,
        setup.feed.n,
        &setup.feed.z,
        setup.feed.p,
        setup.feed.t,
        rich.n,
        &rich.z,
        rich.p,
        rich.t,
        setup.number_of_stages,
        setup.top_pressure,
        setup.bottom_pressure,
        setup.temperature_tolerance,
        setup.max_iterations,
        None,
        None,
        None,
        None,
        None,
        None,
        None,
        None,
        None,
        None,
        None,
        None,
        None,
        None,
        None,
        None,
        None,
        None,
        None,
        None,
        None,
        None,
        None,
        // The two capacity inputs, unstated: this test is about the product names and numbers.
        None,
        None,
    )
    .expect("the stripper converges");

    relative(
        out.overhead_gas_n,
        3.086389912429422,
        1.0e-5,
        "the stripped gas",
    );
    relative(
        out.lean_liquid_n,
        2.4096710196150735,
        1.0e-5,
        "the lean liquid",
    );
    relative(
        out.overhead_gas_z[1],
        0.07619406979290001,
        1.0e-3,
        "the stripped gas's propane",
    );
    absolute(
        out.tray_temperature[0].value,
        321.706608455782,
        2.0e-3,
        "the bottom tray's temperature",
    );
}

/// **The reactive section on an absorber, which is inherited, routed, and refused where the
/// class does not converge it either.**
///
/// `AbsorptionColumn` inherits `setReactive` and overrides no `run`, so the section is the
/// base's own - and an absorber has no ends, so `setReactive(true)` flags **every** stage:
/// the capture's `reactive_tray_flags=11111`. No `AbsorptionColumnTest`, `StrippingColumnTest`
/// or `PackedColumnTest` calls `setReactive`, so the two rows in the capture are this port's
/// own measurement of what the class does with it.
///
/// **On a fluid that reacts the class does not reach a column result** - the hydrocarbon
/// absorber at this port's own gate ends `FALLBACK_PRODUCTS` with a mass residual of `6.4e5`
/// after 16 iterations, and its own log calls the published products not a rigorous column
/// result. This port refuses it one stage earlier and names the reason the reactive tray
/// tranche already codified: the reactive flash's two phases converge to one composition, so
/// nothing says which of them leaves by the vapour port.
#[test]
fn a_reactive_absorber_is_refused_where_the_class_does_not_converge_it_either() {
    let mut setup = lean_oil();
    setup.reactive = azoth_process::kernels::ReactiveSection::All;
    let error = distillation_column(&setup).expect_err("neither library reaches a column result");
    assert!(
        error.to_string().contains("without a vapour/liquid label"),
        "the refusal names the reactive flash's own reason: {error}"
    );
}

/// **And on a fluid with no independent reaction the section is the equilibrium flash**, so the
/// two columns are the same machine - the property `testReactiveColumnMassBalanceNR0` asserts
/// and the reason the packed column's own reactive case is an identity.
#[test]
fn a_reactive_absorber_on_a_fluid_that_does_not_react_is_the_same_column() {
    // **Two species and no more.** The carbon-and-hydrogen formula matrix has rank two, so a
    // third species would give the fluid an independent reaction and the section would stop
    // being the identity; methane against n-heptane is two, and the heavy one is a liquid at
    // 15 bar and 293 K, which an absorber's solvent has to be.
    let light = |reactive| ColumnSetup {
        feed: Stream::from_pt(
            vec!["methane".into(), "n-heptane".into()],
            vec![0.98, 0.02],
            30.852602094576984,
            pascals(15.0e5),
            kelvins(303.15),
        )
        .expect("the fluid resolves"),
        feed_stage: 0,
        number_of_stages: 5,
        has_reboiler: false,
        has_condenser: false,
        top_pressure: pascals(15.0e5),
        bottom_pressure: pascals(15.0e5),
        condenser_temperature: None,
        reboiler_temperature: None,
        temperature_tolerance: 1.0e-4,
        murphree_efficiency: None,
        absorber_murphree: None,
        initial_state: None,
        max_iterations: 80,
        top_specification: None,
        bottom_specification: None,
        top_feed: Some(
            Stream::from_pt(
                vec!["methane".into(), "n-heptane".into()],
                vec![0.0, 1.0],
                1.6632569898375,
                pascals(15.0e5),
                kelvins(293.15),
            )
            .expect("the fluid resolves"),
        ),
        tray_temperatures: None,
        solver_type: SolverType::DirectSubstitution,
        gas_side_draw_fractions: None,
        liquid_side_draw_fractions: None,
        pumparound_fractions: None,
        side_draw_flows: Vec::new(),
        pumparound_returns: Vec::new(),
        pumparound_inlets: Vec::new(),
        pumparound_tolerance: None,
        pumparound_max_iterations: None,
        reactive,
    };
    let plain = distillation_column(&light(azoth_process::kernels::ReactiveSection::None))
        .expect("the plain column converges");
    let reactive = distillation_column(&light(azoth_process::kernels::ReactiveSection::All))
        .expect("and the reactive one is the same machine");

    for (a, b) in plain.trays.iter().zip(&reactive.trays) {
        relative(
            b.temperature.value,
            a.temperature.value,
            1.0e-9,
            "tray temperature",
        );
        relative(b.gas_n, a.gas_n, 1.0e-9, "tray vapour");
    }
    relative(
        reactive.distillate.n,
        plain.distillate.n,
        1.0e-10,
        "distillate",
    );
    relative(reactive.bottoms.n, plain.bottoms.n, 1.0e-10, "bottoms");
    assert_eq!(
        plain.iterations, reactive.iterations,
        "the two routes take the same number of passes"
    );
}

/// **The absorber's Murphree override, and NeqSim's own solve does not converge with it.**
///
/// `validation/neqsim/captures/process_absorber_efficiency.tsv` runs the uncorrected row's own
/// lean-oil state twice: at `0.6` column-wide, and at `1.0` column-wide with **methane** overridden
/// to `0.6` through `setComponentMurphreeEfficiency(String, double)`. NeqSim takes the
/// 80-iteration cap on both, reporting `FAILED` and `FALLBACK_PRODUCTS` with temperature residuals
/// of `3.25` and `22.66` K - where the same state with no correction converges in 17 passes at
/// `9.4e-5` and `RIGOROUS_CONVERGED`. **The correction is what breaks it.**
///
/// This port converges on the same state in seven passes at `6.5e-3`, so the row is evidence of a
/// divergence rather than an oracle - the same family as the `0_85` column row, and asserted as
/// one: what is pinned is that the port *reaches* a state at the row's own gate, that the state
/// **moves** when the correction is stated (so a port that dropped the parameter would fail here),
/// and the resolution the override reads.
///
/// **The correction's arithmetic is oracled at the kernel instead**, where it can be: the
/// allocator's limiting-component loop is `column::absorber_murphree`'s own unit test, because
/// the class invalidates both caches at the end of a solve and a capture can never show them.
#[test]
fn the_absorbers_murphree_override_converges_where_the_class_does_not() {
    // **Both sides at the row's own gate**, so the comparison is between two states and not
    // between a converged one and a partial one.
    let mut uncorrected = lean_oil();
    uncorrected.temperature_tolerance = 1.0e-2;
    let plain = distillation_column(&uncorrected).expect("the uncorrected row converges");

    let mut setup = lean_oil();
    setup.temperature_tolerance = 1.0e-2;
    setup.absorber_murphree = Some(AbsorberMurphree {
        base: Murphree::from_column_wide(0.6),
        per_component: None,
    });
    let corrected = distillation_column(&setup)
        .expect("this port reaches the row's gate where the class takes its iteration cap");

    // **The correction bites**, which is what a dropped parameter would not do - and it bites
    // *weakly*, which is the measurement: this absorber is liquid-film limited, so an efficiency
    // on the vapour blend moves the top tray's methane by `5.3e-4` and no more. A threshold
    // above that would pass a port that had dropped the parameter.
    assert!(
        (corrected.trays[4].gas_z[0] - plain.trays[4].gas_z[0]).abs() > 1.0e-4,
        "the override moves the profile: {} against {}",
        corrected.trays[4].gas_z[0],
        plain.trays[4].gas_z[0]
    );
    assert!(
        corrected.temperature_residual <= 1.0e-2,
        "the row's own gate is met: {}",
        corrected.temperature_residual
    );

    // **And the resolution, which is what the capture pins to its own numbers.**
    // `getComponentMurphreeEfficiency` reads the component's own value, then the base's two
    // steps - so a component at `0.6` under a column-wide `1.0` resolves to `0.6` for that
    // component and `1.0` for the others.
    let resolved = AbsorberMurphree {
        base: Murphree::from_column_wide(1.0),
        per_component: Some(vec![0.6, f64::NAN, f64::NAN, f64::NAN, f64::NAN, f64::NAN]),
    };
    assert_eq!(resolved.resolve(2, 0), 0.6, "methane's own value");
    assert_eq!(
        resolved.resolve(2, 1),
        1.0,
        "and every other component falls through"
    );
    assert_eq!(
        AbsorberMurphree {
            base: Murphree::from_column_wide(0.6),
            per_component: None,
        }
        .resolve(4, 3),
        0.6,
        "and with no component map it is the base's own resolution"
    );
}

/// The class's own constructor defaults, which every capture row of the mechanical design states.
fn mechanical_geometry() -> MechanicalGeometry {
    MechanicalGeometry {
        tray_type: DEFAULT_TRAY_TYPE.to_string(),
        contactor_internals_type: DEFAULT_CONTACTOR_INTERNALS_TYPE.to_string(),
        tray_efficiency: DEFAULT_TRAY_EFFICIENCY,
        max_flooding_factor: DEFAULT_MAX_FLOODING_FACTOR,
        material_grade: DEFAULT_MATERIAL_GRADE.to_string(),
        max_operation_pressure_bara: DEFAULT_MAX_OPERATION_PRESSURE_BARA,
        tray_spacing: meters(0.6),
        weir_height: meters(0.05),
        hole_diameter: meters(12.7e-3),
        hole_area_fraction: 0.1,
        downcommer_area_fraction: 0.1,
        column_diameter_override: meters(-1.0),
    }
}

/// **`process_column_mechanical_design.tsv`'s `absorber_mechanical` row — the one claim the
/// mechanical module's doc comment makes that nothing else holds.**
///
/// `DistillationColumnMechanicalDesign.calcDesign` reads tray 0's *liquid-outlet* density from a
/// fluid it never asks to initialise, so on this column it reads `0.0`, publishes an `Infinity`
/// weir loading and collapses the tray pressure drop to its own bare `5.0` mbar
/// (equinor/neqsim#4140). The port publishes the numbers the class *would* read after one
/// `initProperties()` — `645.1756172849634`, hence `2.664675513235188` and `8.164586402782746` —
/// and this is where that is held. **If azoth's density for this outlet is not that number the
/// module's claim is false, and the deviation is a finding rather than a tolerance to widen.**
#[test]
fn the_absorber_row_publishes_the_post_init_density() {
    let out = distillation_column(&lean_oil()).expect("the absorber converges");
    let components = [
        "methane",
        "ethane",
        "propane",
        "n-butane",
        "n-pentane",
        "n-heptane",
    ]
    .map(str::to_string);
    let report = mechanical_design(
        &out.trays,
        &components,
        &out.distillate,
        &out.bottoms,
        None,
        None,
        &mechanical_geometry(),
    )
    .expect("the mechanical design computes");

    let close = |got: f64, want: f64, what: &str| {
        assert!(
            (got - want).abs() / want.abs() < 1e-5,
            "{what}: {got} against {want}"
        );
    };
    assert_eq!(report.actual_trays, 8, "the actual trays");
    assert_eq!(report.material_grade, "SA-516-70", "the material grade");
    close(report.vessel_height.value, 8.8, "vessel_height");
    // **The capture's block is two designs, and the port publishes the second.** The probe calls
    // `calcDesign` once *before* it asks the liquid outlet's fluid to initialise and once after,
    // so the un-suffixed keys are the **first** design's: `column_diameter_m = 1.0`,
    // `column_wall_thickness_mm = 432.39`, `flooding_factor = 0.6988`, its `total_pressure_drop`
    // and its four `internals_*`. Only the `_after_init` keys are the repaired one, and only those
    // (with the block's path-independent keys) are held here.
    close(report.vessel_diameter.value, 0.5, "vessel_diameter");
    close(report.weir_loading, 2.664_675_513_235_188, "weir_loading");
    close(
        report.tray_pressure_drop_mbar,
        8.164_586_402_782_746,
        "tray_pressure_drop_mbar",
    );
    // The class's thickness is linear in the diameter and the binary row pins it at `0.5` m, so at
    // the repaired diameter it is half the `432.38993710691824` the first design took at `1.0` m.
    close(
        report.vessel_wall_thickness_mm,
        216.194_968_553_459_12,
        "vessel_wall_thickness_mm",
    );
    // Not a key of the block, but the block's own tray-0 readings give it: `0.0494` m3/s of vapour
    // (`111695.677` mol/hr at `0.01843` kg/mol over `11.5763` kg/m3), `u_flood = 0.7398` at the
    // repaired `645.1756` kg/m3, over the `0.1767` m2 net area of the `0.5` m shell.
    close(
        report.flooding_factor,
        0.377_840_006_003_444_86,
        "flooding_factor",
    );
    // The internals, which the block only carries at the defect path's `1.0` m, are held by their
    // own oracle instead: `process_internals_designer.tsv`'s `absorber_designer_defaults`.
    close(
        report.total_pressure_drop_bar,
        0.023_629_991_163_397_854,
        "total_pressure_drop_bar",
    );
    assert_eq!(report.reboiler_duty_kw, 0.0, "no reboiler");
    assert_eq!(report.condenser_duty_kw, 0.0, "no condenser");
}
