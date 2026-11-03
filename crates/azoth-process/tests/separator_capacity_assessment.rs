//! `process.separator_capacity_assessment` against
//! `validation/neqsim/captures/process_separator_capacity_assessment.tsv`.
//!
//! The recorded cases run through `test_model_cases` and `test_cross_impl`; what is here is
//! the branch order, the absences, the fold and the two-id agreement.

use azoth_core::units::{cubic_meters_per_second, kelvins, meters, meters_per_second, pascals};
use azoth_process::models::separator::{SeparatorResult, separator};
use azoth_process::models::separator_capacity_assessment::{
    SeparatorCapacityAssessmentResult, separator_capacity_assessment,
};

const COMPONENTS: [&str; 2] = ["methane", "n-butane"];

fn components() -> Vec<String> {
    COMPONENTS.iter().map(|name| (*name).to_string()).collect()
}

#[allow(clippy::too_many_arguments)]
fn assess(
    components: &[String],
    feed_n: f64,
    z: &[f64],
    p: f64,
    t: f64,
    diameter: f64,
    load_factor: f64,
    orientation: &str,
    level: f64,
    provenance: &str,
    mechanical: Option<f64>,
) -> SeparatorCapacityAssessmentResult {
    separator_capacity_assessment(
        components,
        feed_n,
        z,
        pascals(p),
        kelvins(t),
        pascals(0.0),
        0.0,
        meters(diameter),
        meters_per_second(load_factor),
        orientation,
        level,
        provenance,
        None,
        mechanical.map(cubic_meters_per_second),
    )
    .expect("the captured state separates")
}

fn close(actual: f64, expected: f64, tolerance: f64, what: &str) {
    let error = (actual - expected).abs() / expected.abs().max(f64::MIN_POSITIVE);
    assert!(
        error <= tolerance,
        "{what}: {actual} against {expected}, relative error {error:.3e} > {tolerance:.3e}"
    );
}

#[test]
fn the_vertical_snapshot_matches_the_capture() {
    let result = assess(
        &components(),
        1.0,
        &[0.7, 0.3],
        2_000_000.0,
        300.0,
        1.0,
        0.11,
        "vertical",
        0.8,
        "design basis A",
        None,
    );
    assert_eq!(result.status, "feasible");
    assert_eq!(result.design_basis_source, "equipment");
    assert_eq!(result.governing_capacity_basis, "soudersBrown");
    close(
        result.maximum_gas_velocity.expect("velocity").value,
        0.559_861_417_115_024,
        2.0e-4,
        "maximum_gas_velocity",
    );
    close(
        result.souders_brown_maximum_gas_flow.expect("flow").value,
        0.439_714_128_759_232_6,
        2.0e-4,
        "souders_brown_maximum_gas_flow",
    );
    close(
        result.utilization.expect("utilisation"),
        0.002_126_127_241_504_300_3,
        2.0e-4,
        "utilization",
    );
    assert_eq!(
        result.mechanical_design_maximum_gas_flow.map(|f| f.value),
        Some(0.0)
    );
}

#[test]
fn a_mechanical_limit_at_or_below_the_souders_brown_flow_governs() {
    let result = assess(
        &components(),
        1.0,
        &[0.7, 0.3],
        2_000_000.0,
        300.0,
        1.0,
        0.11,
        "vertical",
        0.8,
        "mechanical design basis",
        Some(0.351_771_303_007_386_1),
    );
    assert_eq!(result.design_basis_source, "mechanicalDesign");
    assert_eq!(result.governing_capacity_basis, "mechanicalDesignGasFlow");
    close(
        result.maximum_gas_flow.expect("maximum flow").value,
        0.351_771_303_007_386_1,
        1.0e-12,
        "maximum_gas_flow",
    );
}

#[test]
fn a_blank_provenance_is_unassessed_and_carries_no_numbers() {
    let result = assess(
        &components(),
        1.0,
        &[0.7, 0.3],
        2_000_000.0,
        300.0,
        1.0,
        0.11,
        "vertical",
        0.8,
        " ",
        None,
    );
    assert_eq!(result.status, "unassessed");
    assert_eq!(result.diagnostic, "Design-basis provenance is required");
    assert_eq!(result.design_basis_source, "");
    assert!(result.gas_area.is_none());
    assert!(result.gas_density.is_none());
    assert!(result.liquid_density.is_none());
    assert!(result.maximum_gas_velocity.is_none());
    assert!(result.current_gas_flow.is_none());
    assert!(result.mechanical_design_maximum_gas_flow.is_none());
    assert!(result.souders_brown_maximum_gas_flow.is_none());
    assert!(result.maximum_gas_flow.is_none());
    assert!(result.utilization.is_none());
}

/// **The class checks the gas phase before the level**, so a liquid-only feed with a level of
/// one carries the geometry sentence.
#[test]
fn the_gas_phase_is_checked_before_the_liquid_level() {
    let result = assess(
        &["n-butane".to_string()],
        1.0,
        &[1.0],
        500_000.0,
        250.0,
        1.0,
        0.11,
        "vertical",
        1.0,
        "liquid only",
        None,
    );
    assert_eq!(result.status, "unassessed");
    assert_eq!(
        result.diagnostic,
        "Finite geometry, design K-factor, and a calculated gas phase are required"
    );
}

#[test]
fn an_overflowing_design_k_is_the_evidence_sentence() {
    let result = assess(
        &components(),
        1.0,
        &[0.7, 0.3],
        2_000_000.0,
        300.0,
        1.0,
        1.0e308,
        "vertical",
        0.8,
        "design basis A",
        None,
    );
    assert_eq!(result.status, "unassessed");
    assert_eq!(
        result.diagnostic,
        "Separator capacity calculation returned non-finite evidence"
    );
}

#[test]
fn the_assessment_and_the_capacity_surface_share_their_numbers() {
    let assessment = assess(
        &components(),
        1.0,
        &[0.7, 0.3],
        2_000_000.0,
        300.0,
        1.0,
        0.11,
        "vertical",
        0.8,
        "design basis A",
        None,
    );
    let capacity: SeparatorResult = separator(
        &components(),
        1.0,
        &[0.7, 0.3],
        pascals(2_000_000.0),
        kelvins(300.0),
        pascals(0.0),
        0.0,
        None,
        Some(meters(1.0)),
        Some(meters_per_second(0.11)),
        Some("vertical"),
        Some(0.8),
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
    )
    .expect("the same vessel separates");
    assert_eq!(
        assessment.maximum_gas_velocity.map(|value| value.value),
        capacity.max_allowable_gas_velocity.map(|value| value.value)
    );
    assert_eq!(assessment.utilization, capacity.capacity_utilization);
}

#[test]
fn a_non_finite_mechanical_limit_is_refused() {
    let result = separator_capacity_assessment(
        &components(),
        1.0,
        &[0.7, 0.3],
        pascals(2_000_000.0),
        kelvins(300.0),
        pascals(0.0),
        0.0,
        meters(1.0),
        meters_per_second(0.11),
        "vertical",
        0.8,
        "design basis A",
        None,
        Some(cubic_meters_per_second(f64::INFINITY)),
    );
    assert!(result.is_err());
}
