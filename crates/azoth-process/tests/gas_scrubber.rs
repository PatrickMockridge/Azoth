//! `process.gas_scrubber`'s mesh-pad configuration against
//! `validation/neqsim/captures/process_gas_scrubber_mechanical_design.tsv`.
//!
//! The booleans have no slot in a `TestCase`, so this is where they are held: the capture's
//! `port_fixture_pr_mesh_pad` row configures a pad and leaves both cyclone kinds unset.

use azoth_core::units::{kelvins, meters, meters_per_second, pascals};
use azoth_process::models::gas_scrubber::gas_scrubber;

fn names() -> Vec<String> {
    ["methane", "ethane", "propane", "n-pentane"]
        .iter()
        .map(|name| (*name).to_string())
        .collect()
}

#[test]
fn a_configured_mesh_pad_turns_its_flag_on() {
    let result = gas_scrubber(
        &names(),
        107.408_524_456_061_74,
        &[0.7, 0.1, 0.1, 0.1],
        pascals(2_000_000.0),
        kelvins(298.15),
        pascals(0.0),
        0.0,
        Some(meters(0.515_375_309_329_634_3)),
        Some(meters_per_second(0.11)),
        None,
        Some(true),
        Some("vertical"),
        None,
        None,
        None,
        None,
        None,
        None,
        Some(1.5),
        Some(0.208_610_948_791_937),
        Some(0.15),
        Some(meters(0.75)),
    )
    .expect("the captured state separates");

    assert_eq!(result.has_mesh_pad, Some(true));
    assert_eq!(result.mesh_pad_area, Some(0.208_610_948_791_937));
    assert_eq!(result.mesh_pad_thickness, Some(0.150));
    assert_eq!(result.has_inlet_cyclones, Some(false));
    assert_eq!(result.has_demisting_cyclones, Some(false));
    let diameter = result.required_diameter.expect("a stated safety factor").value;
    let error = (diameter - 0.603_241_350_383_609_7).abs() / 0.603_241_350_383_609_7;
    assert!(error <= 2.0e-4, "sizeFromFlow: {diameter} ({error:.3e})");
}
