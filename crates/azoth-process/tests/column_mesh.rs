//! Naphtali-Sandholm, against `validation/neqsim/captures/process_column_solvers.tsv`.
//!
//! The capture runs both columns under every strategy `ColumnSolverFactory` carries, and the
//! rows that matter here are the two `naphtali_sandholm` ones: the mesh solve converges the
//! **binary** in four iterations and the **deethanizer** in six, where every other strategy
//! falls back to `DAMPED_SUBSTITUTION` and stops at the iteration cap.
//!
//! **The port's mesh solve starts from the substitution core's answer** - the class's own
//! `initializeTrayStateFromColumn` warm start, which is why it converges in five iterations to
//! NeqSim's four. The substitution core's fixed point is *not* a root of the MESH equations
//! (its worst residual is `6e-3` on the binary), so the fifth iteration is the one that closes
//! the equations the class's own solve closes - and the state it lands on is NeqSim's.

use azoth_core::units::{kelvins, pascals};
use azoth_process::Stream;
use azoth_process::kernels::distillation_column::{
    ColumnSetup, SolverType, Specification, SpecificationKind, distillation_column,
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

/// The captured binary column, under a stated solve.
fn binary(solver_type: SolverType) -> ColumnSetup {
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
        max_iterations: 200,
        top_specification: None,
        bottom_specification: None,
        solver_type,
    }
}

/// NeqSim's own deethanizer, under a stated solve.
fn deethanizer(solver_type: SolverType) -> ColumnSetup {
    ColumnSetup {
        feed: Stream::from_pt(
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
        .expect("the fluid resolves"),
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
        solver_type,
    }
}

/// **The mesh solve on the captured binary column**, against NeqSim's own `NAPHTALI_SANDHOLM`
/// row - which reaches the same fixed point as its substitution row to seven figures, and this
/// port's to five: every tray temperature within `1.6e-4` K, every traffic rate within `5e-6`
/// relative, and the distillate within `3e-13`.
///
/// **The declared divergence budget is the libraries' enthalpy offset**, as it is for the
/// substitution row: the two duties are enthalpy differences and agree to `0.4` W of `21323`
/// and `1.1` W of `47786`.
#[test]
fn the_binary_column_converges_under_the_mesh_solve() {
    let out = distillation_column(&binary(SolverType::NaphtaliSandholm))
        .expect("the mesh solve converges");

    // NeqSim: 4 iterations, status RIGOROUS_CONVERGED. The port takes 5, because its seed is
    // the substitution core's answer rather than the class's cold BP/Boston-Sullivan profile.
    assert!(out.iterations <= 8, "in {} iterations", out.iterations);

    let temperatures = [
        373.15,
        336.1538231834974,
        303.07105352147533,
        302.12253908101036,
        296.84208915683996,
        253.14999999999998,
    ];
    let gas = [
        1.3732781432407117,
        0.525903134833042,
        4.576181492897451,
        4.55803955854207,
        4.451759778195404,
        3.7914996439673083,
    ];
    let liquid = [
        3.6992043923236553,
        5.072482535564368,
        4.225107527156697,
        0.7846818489301433,
        0.7665399145747628,
        0.6602601342280967,
    ];
    assert_eq!(out.trays.len(), 6, "the two ends are trays");
    for (i, tray) in out.trays.iter().enumerate() {
        absolute(
            tray.temperature.value,
            temperatures[i],
            2.0e-4,
            "tray temperature",
        );
        relative(tray.gas_n, gas[i], 1.0e-5, "tray vapour");
        relative(tray.liquid_n, liquid[i], 1.0e-5, "tray liquid");
    }

    relative(out.distillate.n, 3.7914996439673083, 1.0e-6, "distillate");
    relative(out.bottoms.n, 3.6992043923236553, 1.0e-6, "bottoms");
    relative(
        out.distillate.z[0],
        0.9659099052335057,
        1.0e-5,
        "the distillate's methane",
    );
    relative(
        out.bottoms.z[1],
        0.9775343702217995,
        1.0e-6,
        "the bottoms' n-butane",
    );
    absolute(
        out.condenser_duty.value,
        -21323.042177706004,
        5.0,
        "condenser duty",
    );
    absolute(
        out.reboiler_duty.value,
        47786.58295288164,
        5.0,
        "reboiler duty",
    );

    // The port's own closure, which is what a converged MESH solve must show.
    assert!(out.mass_residual < 1.0e-8, "mass {}", out.mass_residual);
    assert!(
        out.energy_residual < 1.0e-8,
        "energy {}",
        out.energy_residual
    );
    // The mesh residual is the class's own convergence measure, and this is the field the
    // column reports it in for this solve.
    assert!(
        out.temperature_residual < 1.0e-8,
        "mesh residual {}",
        out.temperature_residual
    );
}

/// **The deethanizer, which no other strategy converges.** NeqSim's `DIRECT_SUBSTITUTION`
/// stops at the iteration cap with a residual of `7.3e-3` K after 80 iterations; its mesh
/// solve reaches `RIGOROUS_CONVERGED` in six. Every tray temperature here is within `3.5e-5` K
/// of that row and both products within `5e-8` relative.
#[test]
fn the_deethanizer_converges_under_the_mesh_solve() {
    let out = distillation_column(&deethanizer(SolverType::NaphtaliSandholm))
        .expect("the mesh solve converges");

    assert!(out.iterations <= 10, "in {} iterations", out.iterations);

    let temperatures = [
        353.15,
        333.7864384739617,
        322.64122258578095,
        315.7792368087982,
        310.55162130108437,
        303.844339682725,
        290.9017562675362,
        287.7398599203877,
        283.9719467063039,
        273.15,
    ];
    assert_eq!(out.trays.len(), 10);
    for (i, tray) in out.trays.iter().enumerate() {
        absolute(
            tray.temperature.value,
            temperatures[i],
            5.0e-4,
            "tray temperature",
        );
    }

    relative(out.distillate.n, 40.286439082548355, 1.0e-6, "distillate");
    relative(out.bottoms.n, 32.95771041675448, 1.0e-6, "bottoms");
    relative(
        out.distillate.z[0],
        0.3998917170028344,
        1.0e-5,
        "the distillate's methane",
    );
    relative(
        out.bottoms.z[2],
        0.33148810298445563,
        1.0e-5,
        "the bottoms' propane",
    );
    absolute(
        out.condenser_duty.value,
        -104063.02333018146,
        20.0,
        "condenser duty",
    );
    absolute(
        out.reboiler_duty.value,
        547533.2856140495,
        60.0,
        "reboiler duty",
    );

    assert!(out.mass_residual < 1.0e-8, "mass {}", out.mass_residual);
    assert!(
        out.energy_residual < 1.0e-8,
        "energy {}",
        out.energy_residual
    );
}

/// **The pair, pinned against each other on both columns.**
///
/// This is the deliverable's own check: one column, two solves, one compare. The substitution
/// core's fixed point and the mesh solve's agree to `1.6e-4` K on the binary and `3.5e-5` K on
/// the deethanizer - which is the same band each of them agrees with NeqSim's own solve of the
/// same kind, so the two libraries' two methods land on one state.
///
/// **The mesh solve closes the equations the substitution core only approximates.** Its own
/// residual at the substitution core's answer is `6e-3` on the binary - the substitution core
/// converges on its *own* closure, which its trays' flashes decide, and the MESH equations are
/// a slightly different system.
#[test]
fn the_two_solves_reach_the_same_state() {
    for (name, setup, band) in [
        ("binary", binary(SolverType::DirectSubstitution), 2.0e-4),
        (
            "deethanizer",
            deethanizer(SolverType::DirectSubstitution),
            5.0e-4,
        ),
    ] {
        let mut mesh_setup = setup.clone();
        mesh_setup.solver_type = SolverType::NaphtaliSandholm;
        let substitution = distillation_column(&setup).expect("the substitution core converges");
        let mesh = distillation_column(&mesh_setup).expect("the mesh solve converges");

        assert_eq!(substitution.trays.len(), mesh.trays.len());
        for (i, (a, b)) in substitution.trays.iter().zip(mesh.trays.iter()).enumerate() {
            absolute(
                b.temperature.value,
                a.temperature.value,
                band,
                &format!("{name} tray {i} temperature"),
            );
            relative(b.gas_n, a.gas_n, 1.0e-4, &format!("{name} tray {i} vapour"));
            relative(
                b.liquid_n,
                a.liquid_n,
                1.0e-4,
                &format!("{name} tray {i} liquid"),
            );
        }
        relative(mesh.distillate.n, substitution.distillate.n, 1.0e-5, name);
        relative(mesh.bottoms.n, substitution.bottoms.n, 1.0e-5, name);
        absolute(
            mesh.condenser_duty.value,
            substitution.condenser_duty.value,
            20.0,
            &format!("{name} condenser duty"),
        );
        absolute(
            mesh.reboiler_duty.value,
            substitution.reboiler_duty.value,
            60.0,
            &format!("{name} reboiler duty"),
        );
    }
}

/// A mesh solve needs a pinned reboiler, and says which class closes the alternative.
#[test]
fn a_reboiler_with_no_pin_is_refused_by_name() {
    let mut setup = binary(SolverType::NaphtaliSandholm);
    setup.reboiler_temperature = None;
    let error = distillation_column(&setup).expect_err("the mesh solve refuses it");
    let message = error.to_string();
    assert!(
        message.contains("solveBubblePointMethod"),
        "the refusal names the class that would close it: {message}"
    );
}

/// A specification drives the substitution core's outer loop, which the mesh solve has no
/// counterpart for.
#[test]
fn a_specification_is_refused_by_name() {
    let mut setup = binary(SolverType::NaphtaliSandholm);
    setup.top_specification = Some(Specification {
        kind: SpecificationKind::ProductPurity,
        target: 0.98,
        component: Some("methane".into()),
    });
    let error = distillation_column(&setup).expect_err("the mesh solve refuses it");
    let message = error.to_string();
    assert!(
        message.contains("solveWithSpecificationTargets"),
        "the refusal names the class that would close it: {message}"
    );
}
