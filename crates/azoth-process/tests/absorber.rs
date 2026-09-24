//! The absorber and the stripper, against `validation/neqsim/captures/process_absorber.tsv`.

use azoth_core::units::{kelvins, pascals};
use azoth_process::Stream;
use azoth_process::kernels::distillation_column::{ColumnSetup, SolverType, distillation_column};

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
    }
}
