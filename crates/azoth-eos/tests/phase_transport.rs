//! Spec-driven tests for the `eos.phase_transport` model.
//!
//! The oracle is `validation/neqsim/PhaseTransportProbe.java` and its capture
//! `captures/phase_transport_probe.tsv`: four fluids' six phases, each with the viscosity and
//! conductivity NeqSim answers, the binary matrix behind them and **the model class each
//! property came from**. The cases hold the two scalars; the matrices and the model classes are
//! held here, because a matrix is a shape the case generator does not carry.

use azoth_core::units::{kelvins, pascals};
use azoth_eos::databank::mixture_of;
use azoth_eos::phase_transport::{PhaseKind, phase_transport};
use azoth_test_support as common;

const MODEL_ID: &str = "eos.phase_transport";

/// The fluid and the phase a test runs: what the capture's row says.
struct State {
    label: &'static str,
    components: &'static [&'static str],
    phase: PhaseKind,
    t: f64,
    p: f64,
    z: &'static [f64],
    /// NeqSim's own viscosity and conductivity for that phase.
    mu: f64,
    k: f64,
    /// Its binary matrix, row-major, zeros on the diagonal.
    d: &'static [&'static [f64]],
}

/// **The six phases the capture holds**, with NeqSim's numbers for each.
///
/// The two CO2/water rows are the state `RateBasedPackedColumnTest` runs its absorber at, and
/// they are the pair the segment model needs: a gas phase whose diffusivity is Chapman-Enskog
/// and an **aqueous** phase whose viscosity, conductivity and diffusivity all come from
/// *different* classes again.
const STATES: &[State] = &[
    State {
        label: "methane_nitrogen_298k_1atm_gas",
        components: &["methane", "nitrogen"],
        phase: PhaseKind::Gas,
        t: 298.15,
        p: 101325.0,
        z: &[0.5, 0.5],
        mu: 1.4506127676908967e-5,
        k: 0.029958794899910723,
        d: &[&[0.0, 3.555070057344887e-5], &[3.555070057344887e-5, 0.0]],
    },
    State {
        label: "methane_co2_313k_50bar_gas",
        components: &["methane", "CO2"],
        phase: PhaseKind::Gas,
        t: 313.15,
        p: 5.0e6,
        z: &[0.5, 0.5],
        mu: 1.6159655797683304e-5,
        k: 0.030071678065251736,
        d: &[&[0.0, 9.457085848059925e-7], &[9.457085848059925e-7, 0.0]],
    },
    State {
        label: "methane_nbutane_300k_20bar_gas",
        components: &["methane", "n-butane"],
        phase: PhaseKind::Gas,
        t: 300.0,
        p: 2.0e6,
        z: &[0.8356249351028912, 0.1643750648971088],
        mu: 1.1173254776369514e-5,
        k: 0.032654729853716064,
        d: &[&[0.0, 1.6228425436789486e-6], &[1.6228425436789486e-6, 0.0]],
    },
    State {
        label: "methane_nbutane_300k_20bar_oil",
        components: &["methane", "n-butane"],
        phase: PhaseKind::Oil,
        t: 300.0,
        p: 2.0e6,
        z: &[0.09542082642782022, 0.9045791735721797],
        mu: 1.3466360688772563e-4,
        k: 0.09534463941597507,
        // **Asymmetric**, and by a factor of six: the correlation divides by the *solvent*
        // component's pure viscosity, and the common-phase ladder answers zero for n-butane's
        // LIQVISC model 2 - which the clamp turns into 0.01 cP.
        d: &[&[0.0, 9.552345101712192e-9], &[1.5074311608978982e-9, 0.0]],
    },
    State {
        label: "co2_water_313k_50bar_gas",
        components: &["CO2", "water"],
        phase: PhaseKind::Gas,
        t: 313.15,
        p: 5.0e6,
        z: &[0.9978562046113857, 0.002143795388614308],
        mu: 1.8462596972222117e-5,
        k: 0.024804749452710006,
        d: &[&[0.0, 6.307823319063167e-7], &[6.307823319063167e-7, 0.0]],
    },
    State {
        label: "co2_water_313k_50bar_aqueous",
        components: &["CO2", "water"],
        phase: PhaseKind::Aqueous,
        t: 313.15,
        p: 5.0e6,
        z: &[0.0006691762234084198, 0.9993308237765915],
        mu: 6.526078179332641e-4,
        k: 0.6344057039895419,
        d: &[&[0.0, 1.919021499279848e-9], &[3.9028244166575445e-9, 0.0]],
    },
];

fn run(state: &State) -> azoth_eos::PhaseTransportResult {
    let (mixture, ideal_gas) = mixture_of(state.components, azoth_eos::Cubic::Pr, None)
        .unwrap_or_else(|e| panic!("{}: the fluid resolves: {e}", state.label));
    phase_transport(
        &mixture,
        &ideal_gas,
        state.phase,
        kelvins(state.t),
        pascals(state.p),
        state.z,
    )
    .unwrap_or_else(|e| panic!("{}: the phase transports: {e}", state.label))
}

/// **Every phase of every captured fluid, against NeqSim's own numbers.**
#[test]
fn the_captured_phases_are_reproduced() {
    for state in STATES {
        let out = run(state);
        relative(
            out.mu.value,
            state.mu,
            1.0e-9,
            &format!("{} mu", state.label),
        );
        // **NeqSim's PFCT conductivity agrees to its own band and not to the digits**: the
        // port's `eos.thermal_conductivity` case is stated at `1e-5` and this state is measured
        // at `5.7e-5`, because the correlation's reference methane flash is a solve of its own.
        // The polynom conductivities an aqueous phase takes agree to `1e-15`.
        relative(out.k.value, state.k, 1.0e-4, &format!("{} k", state.label));
        for (i, row) in state.d.iter().enumerate() {
            for (j, expected) in row.iter().enumerate() {
                relative(
                    out.d_binary[i][j],
                    *expected,
                    1.0e-9,
                    &format!("{} d[{i}][{j}]", state.label),
                );
            }
        }
    }
}

/// **The effective vector is assembled from the matrix the model answered**, which is the
/// assembly `eos.effective_diffusion` carries and the one NeqSim does *not* do for a liquid.
///
/// Measured: the class's aqueous and oil phases report an effective vector of zero, because
/// their `SiddiqiLucasMethod` never computes one - so this test is the port's own consistency
/// rather than an oracle, and it says so.
#[test]
fn the_effective_vector_is_the_assembly_over_the_matrix() {
    for state in STATES {
        let out = run(state);
        let assembly = azoth_eos::effective_diffusion::effective_diffusion(&out.d_binary, state.z)
            .expect("the assembly runs");
        for (i, (got, expected)) in out
            .d_effective
            .iter()
            .zip(&assembly.effective_diffusion)
            .enumerate()
        {
            assert_eq!(
                got, expected,
                "{}: d_effective[{i}] is the assembly over the same matrix",
                state.label
            );
        }
        // And it is a real number, not the class's zero.
        assert!(
            out.d_effective.iter().all(|value| *value > 0.0),
            "{}: the port assembles an effective coefficient",
            state.label
        );
    }
}

/// **The dispatch is not a relabelling.** The same fluid's gas and aqueous phases answer
/// through different correlations, so a dispatch that ran one family for both would fail on at
/// least one of them.
#[test]
fn the_phase_kind_decides_the_correlations() {
    let (mixture, ideal_gas) =
        mixture_of(&["CO2", "water"], azoth_eos::Cubic::Pr, None).expect("the fluid");
    let gas = phase_transport(
        &mixture,
        &ideal_gas,
        PhaseKind::Gas,
        kelvins(313.15),
        pascals(5.0e6),
        &[0.9978562046113857, 0.002143795388614308],
    )
    .expect("the gas phase transports");
    let aqueous = phase_transport(
        &mixture,
        &ideal_gas,
        PhaseKind::Aqueous,
        kelvins(313.15),
        pascals(5.0e6),
        &[0.0006691762234084198, 0.9993308237765915],
    )
    .expect("the aqueous phase transports");
    // PFCT for the gas, the polynom models for the aqueous: their conductivities differ by a
    // factor of twenty-five, and their diffusivities by six.
    assert!(aqueous.k.value > 20.0 * gas.k.value);
    assert!(gas.d_binary[0][1] > 100.0 * aqueous.d_binary[0][1]);
}

/// A pure phase is refused by name: NeqSim answers a `NaN` there, and a `NaN` is not a state.
#[test]
fn a_single_component_phase_is_refused() {
    let (mixture, ideal_gas) =
        mixture_of(&["methane"], azoth_eos::Cubic::Pr, None).expect("the fluid");
    let error = phase_transport(
        &mixture,
        &ideal_gas,
        PhaseKind::Gas,
        kelvins(300.0),
        pascals(1.0e5),
        &[1.0],
    )
    .expect_err("one component has no pairwise diffusivity");
    assert!(format!("{error}").contains("NaN"), "{error}");
}

/// A relative comparison against the capture, with the band the two libraries' own arithmetic
/// allows: every number here is a quotient of the same correlations, so the agreement is to
/// `1e-9` rather than to a band - which is what makes a wrong dispatch fail loudly.
fn relative(actual: f64, expected: f64, tolerance: f64, what: &str) {
    let scale = actual.abs().max(expected.abs());
    assert!(
        (actual - expected).abs() <= tolerance * scale,
        "{what}: {actual} vs {expected}, {} relative",
        (actual - expected).abs() / scale
    );
}

/// The model's own cases, which hold the two scalars through the spec rather than here.
#[test]
fn every_case_in_the_spec() {
    let spec = azoth_eos::model_gen::model(MODEL_ID).expect("the model should be in its table");
    assert!(!spec.cases.is_empty());
    for case in spec.cases {
        let state = STATES
            .iter()
            .find(|state| state.label == case.id)
            .unwrap_or_else(|| panic!("case `{}` has no captured state", case.id));
        let out = run(state);
        common::assert_close(
            out.mu.value,
            case.expected_value("mu").expect("the case states mu"),
            case.tolerance,
            &format!("{}::{} (mu)", spec.id, case.id),
        );
        common::assert_close(
            out.k.value,
            case.expected_value("k").expect("the case states k"),
            // The case states `1e-12`; the PFCT conductivity's own band is wider, as the
            // model's test says.
            case.tolerance.max(1.0e-4),
            &format!("{}::{} (k)", spec.id, case.id),
        );
        common::assert_consistent(&out, &format!("{}::{}", spec.id, case.id));
    }
}
