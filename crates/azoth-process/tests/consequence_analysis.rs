//! `process.consequence_analysis` against
//! `validation/neqsim/captures/process_consequence_analysis.tsv`.

use azoth_core::units::{
    kelvins, kilograms, kilograms_per_mole, kilograms_per_second,
    kilograms_per_square_meter_second, meters, meters_per_second, pascals, seconds,
};
use azoth_process::models::consequence_analysis::{
    ConsequenceAnalysisResult, ConsequenceOutcomeSpec, consequence_analysis,
};

fn close(actual: f64, expected: f64, what: &str) {
    let error = (actual - expected).abs() / expected.abs().max(f64::MIN_POSITIVE);
    assert!(
        error <= 1.0e-9,
        "{what}: {actual} against {expected}, relative error {error:.3e}"
    );
}

fn spec(kind: &str, probability: f64, a: f64, b: f64, n: f64) -> ConsequenceOutcomeSpec {
    ConsequenceOutcomeSpec {
        outcome_type: kind.to_string(),
        conditional_probability: probability,
        probit_a: a,
        probit_b: b,
        probit_n: n,
        exposure_seconds: None,
        mass_flow: None,
        heat_of_combustion: None,
        radiative_fraction: None,
        relative_humidity: None,
        pool_diameter: None,
        burn_rate: None,
        air_density: None,
        cloud_fuel_mass: None,
        strength_class: None,
        ambient_pressure: None,
        emission_rate: None,
        effective_stack_height: None,
        wind_speed: None,
        stability: None,
        terrain: None,
        molar_mass: None,
        temperature: None,
        pressure: None,
    }
}

fn thermal_jet(probability: f64, exposure: f64) -> ConsequenceOutcomeSpec {
    let mut row = spec("jet_fire", probability, -14.9, 2.56, 1.0);
    row.exposure_seconds = Some(seconds(exposure));
    row.mass_flow = Some(kilograms_per_second(30.0));
    row.heat_of_combustion = Some(50.0e6);
    row.radiative_fraction = Some(0.25);
    row
}

fn thermal_pool(probability: f64) -> ConsequenceOutcomeSpec {
    let mut row = spec("pool_fire", probability, -14.9, 2.56, 1.0);
    row.exposure_seconds = Some(seconds(60.0));
    row.pool_diameter = Some(meters(10.0));
    row.burn_rate = Some(kilograms_per_square_meter_second(0.06));
    row
}

fn blast_vce(probability: f64) -> ConsequenceOutcomeSpec {
    let mut row = spec("vce", probability, -77.1, 6.91, 1.0);
    row.cloud_fuel_mass = Some(kilograms(1000.0));
    row.heat_of_combustion = Some(46.0e6);
    row.strength_class = Some(4.0);
    row
}

fn h2s_plume(probability: f64, emission: f64) -> ConsequenceOutcomeSpec {
    let mut row = spec("toxic_dispersion", probability, -31.42, 3.008, 1.43);
    row.emission_rate = Some(kilograms_per_second(emission));
    row.effective_stack_height = Some(meters(10.0));
    row.wind_speed = Some(meters_per_second(5.0));
    row.stability = Some("D".to_string());
    row.terrain = Some("rural".to_string());
    row.molar_mass = Some(kilograms_per_mole(0.0340819));
    row.temperature = Some(kelvins(288.15));
    row.pressure = Some(pascals(101_325.0));
    row.exposure_seconds = Some(seconds(600.0));
    row
}

fn run(
    release: f64,
    distance: f64,
    outcomes: &[ConsequenceOutcomeSpec],
) -> ConsequenceAnalysisResult {
    consequence_analysis(release, meters(distance), outcomes).expect("the captured state is valid")
}

#[test]
fn the_class_test_state_reproduces_the_capture() {
    let result = run(1.0e-4, 10.0, &[thermal_jet(0.05, 60.0)]);
    assert_eq!(result.outcomes.len(), 1);
    let row = &result.outcomes[0];
    assert_eq!(row.outcome, "Jet fire");
    assert_eq!(row.outcome_type, "jet_fire");
    close(row.frequency_per_year, 5.0e-6, "the branch frequency");
    close(row.fatality_probability, 1.0, "the saturated probability");
    close(
        row.fatality_frequency_per_year,
        5.0e-6,
        "the fatality frequency",
    );
    close(
        result.individual_fatality_risk_per_year,
        5.0e-6,
        "the individual risk",
    );
}

#[test]
fn the_pool_vce_and_toxic_branches_reproduce_the_capture() {
    let pool = run(1.0e-4, 20.0, &[thermal_pool(0.10)]);
    close(
        pool.outcomes[0].fatality_probability,
        1.055_358_503_372_905_2e-4,
        "the pool probability",
    );
    close(
        pool.individual_fatality_risk_per_year,
        1.055_358_503_372_905_2e-9,
        "the pool risk",
    );

    let vce = run(1.0e-4, 20.0, &[blast_vce(0.30)]);
    close(
        vce.outcomes[0].frequency_per_year,
        3.0e-5,
        "the VCE frequency",
    );
    close(
        vce.outcomes[0].fatality_probability,
        0.0,
        "the VCE probability below the probit onset",
    );

    let toxic = run(1.0e-4, 50.0, &[h2s_plume(0.20, 1.0)]);
    close(
        toxic.outcomes[0].frequency_per_year,
        2.0e-5,
        "the toxic frequency",
    );
    close(
        toxic.outcomes[0].fatality_probability,
        0.0,
        "the toxic probability",
    );
}

#[test]
fn the_mixed_scenario_keeps_caller_order_and_sums_its_risk() {
    let result = run(
        1.0e-3,
        20.0,
        &[
            thermal_jet(0.05, 60.0),
            thermal_pool(0.10),
            blast_vce(0.30),
            h2s_plume(0.20, 1.0),
        ],
    );
    assert_eq!(result.outcomes.len(), 4);
    let labels: Vec<&str> = result
        .outcomes
        .iter()
        .map(|row| row.outcome.as_str())
        .collect();
    assert_eq!(
        labels,
        vec!["Jet fire", "Pool fire", "VCE", "Toxic dispersion"]
    );
    close(
        result.outcomes[0].fatality_probability,
        0.999_185_140_471_441,
        "the jet probability",
    );
    close(
        result.outcomes[1].fatality_probability,
        1.055_358_503_372_905_2e-4,
        "the pool probability",
    );
    close(
        result.individual_fatality_risk_per_year,
        4.996_981_060_860_579e-5,
        "the summed risk",
    );
}

#[test]
fn a_zero_dose_is_a_zero_probability() {
    let exposure = run(1.0e-4, 10.0, &[thermal_jet(0.05, 0.0)]);
    close(
        exposure.outcomes[0].frequency_per_year,
        5.0e-6,
        "the frequency",
    );
    close(
        exposure.outcomes[0].fatality_probability,
        0.0,
        "a zero exposure's probability",
    );

    let emission = run(1.0e-4, 50.0, &[h2s_plume(0.20, 0.0)]);
    close(
        emission.outcomes[0].fatality_probability,
        0.0,
        "a zero emission's probability",
    );
}

#[test]
fn the_guards_are_the_score_contract() {
    assert!(
        consequence_analysis(-1.0e-4, meters(10.0), &[]).is_err(),
        "a negative release frequency is refused"
    );
    assert!(
        consequence_analysis(1.0e-4, meters(0.0), &[thermal_jet(0.05, 60.0)]).is_err(),
        "a zero receiver distance is refused"
    );
    let mut bad_probability = thermal_jet(0.05, 60.0);
    bad_probability.conditional_probability = 1.5;
    assert!(
        consequence_analysis(1.0e-4, meters(10.0), &[bad_probability]).is_err(),
        "a probability above one is refused"
    );
    let mut nan_probability = thermal_jet(0.05, 60.0);
    nan_probability.conditional_probability = f64::NAN;
    assert!(
        consequence_analysis(1.0e-4, meters(10.0), &[nan_probability]).is_err(),
        "a NaN probability is refused"
    );
    let mut missing = thermal_jet(0.05, 60.0);
    missing.mass_flow = None;
    assert!(
        consequence_analysis(1.0e-4, meters(10.0), &[missing]).is_err(),
        "a jet branch without its mass flow is refused"
    );
    let mut foreign = thermal_jet(0.05, 60.0);
    foreign.pool_diameter = Some(meters(10.0));
    assert!(
        consequence_analysis(1.0e-4, meters(10.0), &[foreign]).is_err(),
        "a jet branch carrying a pool field is refused"
    );
    let mut unknown = thermal_jet(0.05, 60.0);
    unknown.outcome_type = "flash_fire".to_string();
    assert!(
        consequence_analysis(1.0e-4, meters(10.0), &[unknown]).is_err(),
        "an unknown outcome type is refused"
    );
}
