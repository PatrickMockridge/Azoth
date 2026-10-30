//! `process.hem_release` against the states its cases state.

use azoth_core::units::{kelvins, meters, pascals};
use azoth_process::models::hem_release::ReleaseStation;
use azoth_process::{HemReleaseResult, hem_release};

fn pure(component: &str, t: f64, p: f64, d: f64, cd: f64, back: f64) -> HemReleaseResult {
    let components = vec![component.to_string()];
    hem_release(
        &components,
        &[1.0],
        kelvins(t),
        pascals(p),
        meters(d),
        cd,
        pascals(back),
    )
    .expect("the state runs")
}

fn station(
    result: &HemReleaseResult,
    which: ReleaseStation,
) -> &azoth_process::models::hem_release::HemReleaseStation {
    result
        .stations
        .iter()
        .find(|row| row.station == which)
        .expect("the station is present")
}

#[test]
fn no_forward_flow_is_a_valid_zero_with_no_ambient_station() {
    let result = pure("methane", 300.0, 500_000.0, 0.01, 1.0, 500_000.0);
    assert_eq!(result.status, "valid");
    assert_eq!(result.mass_flow_rate.expect("a true zero").value, 0.0);
    assert!(!result.choked);
    assert!(
        result
            .diagnostic_codes
            .iter()
            .any(|code| code == "NO_FORWARD_FLOW")
    );
    assert_eq!(result.stations.len(), 3);
    assert!(
        !result
            .stations
            .iter()
            .any(|row| row.station == ReleaseStation::AmbientExpanded)
    );
}

#[test]
fn a_receiving_pressure_above_the_upstream_is_the_same_zero() {
    let result = pure("methane", 300.0, 500_000.0, 0.01, 1.0, 600_000.0);
    assert_eq!(result.status, "valid");
    assert_eq!(result.mass_flow_rate.expect("a true zero").value, 0.0);
    assert!(
        result
            .diagnostic_codes
            .iter()
            .any(|code| code == "NO_FORWARD_FLOW")
    );
}

#[test]
fn cold_co2_is_unsupported_and_carries_no_payload() {
    let result = pure("CO2", 210.0, 2_000_000.0, 0.01, 0.62, 100_000.0);
    assert_eq!(result.status, "unsupported");
    assert!(result.mass_flow_rate.is_none());
    assert!(result.throat_sound_speed.is_none());
    assert!(result.stations.is_empty());
    assert_eq!(
        result.diagnostic_codes,
        vec!["UNSUPPORTED_PHYSICS".to_string()]
    );
}

#[test]
fn a_choked_methane_release_has_four_stations_and_a_conserved_entropy() {
    let result = pure("methane", 300.0, 5_000_000.0, 0.01, 0.62, 101_325.0);
    assert_eq!(result.status, "valid");
    assert!(result.choked);
    let rate = result.mass_flow_rate.expect("a usable rate").value;
    assert!(rate > 0.0);
    assert_eq!(result.stations.len(), 4);

    let upstream = station(&result, ReleaseStation::UpstreamStagnation);
    let throat = station(&result, ReleaseStation::ThroatCritical);
    let exit = station(&result, ReleaseStation::OrificeExit);
    assert_eq!(upstream.velocity.value, 0.0);
    assert!((throat.pressure.value - exit.pressure.value).abs() <= 1.0e-9 * throat.pressure.value);
    // The expansion holds the specific entropy; the station carries the molar one.
    let entropy_error = (throat.molar_entropy.value - upstream.molar_entropy.value).abs()
        / upstream.molar_entropy.value.abs();
    assert!(
        entropy_error <= 1.0e-9,
        "entropy moved by {entropy_error:.3e}"
    );
    assert!(throat.velocity.value > 0.0);
    assert!(
        result
            .throat_sound_speed
            .expect("a choked throat has a speed")
            .value
            > 0.0
    );
}

#[test]
fn an_unchoked_release_stops_at_the_receiving_pressure() {
    let result = pure("methane", 300.0, 500_000.0, 0.01, 1.0, 490_000.0);
    assert_eq!(result.status, "valid");
    assert!(!result.choked);
    assert!(result.throat_sound_speed.is_none());
    let throat = station(&result, ReleaseStation::ThroatCritical);
    assert!((throat.pressure.value - 490_000.0).abs() <= 1.0e-9 * 490_000.0);
}

#[test]
fn a_flashing_pure_liquid_keeps_its_mass_fraction_and_warns_at_the_boundary() {
    let result = pure("propane", 300.0, 2_000_000.0, 0.01, 0.62, 300_000.0);
    assert_eq!(result.status, "valid_with_warnings");
    assert!(result.choked);
    assert!(
        result
            .diagnostic_codes
            .iter()
            .any(|code| code == "THROAT_MACH_MISMATCH")
    );
    let ambient = station(&result, ReleaseStation::AmbientExpanded);
    assert_eq!(ambient.phase.as_str(), "two_phase");
    assert!(ambient.gas_mass_fraction > 0.0 && ambient.gas_mass_fraction < 1.0);
}

#[test]
fn the_area_and_the_coefficient_scale_the_rate_and_nothing_else() {
    let small = pure("methane", 300.0, 500_000.0, 0.01, 0.5, 101_325.0);
    let large = pure("methane", 300.0, 500_000.0, 0.02, 1.0, 101_325.0);
    let small_rate = small.mass_flow_rate.expect("a usable rate").value;
    let large_rate = large.mass_flow_rate.expect("a usable rate").value;
    let ratio = large_rate / small_rate;
    assert!((ratio - 8.0).abs() <= 1.0e-12, "the ratio is {ratio}");
    let small_throat = station(&small, ReleaseStation::ThroatCritical);
    let large_throat = station(&large, ReleaseStation::ThroatCritical);
    assert!(
        (small_throat.pressure.value - large_throat.pressure.value).abs()
            <= 1.0e-9 * small_throat.pressure.value
    );
}

#[test]
fn a_mixture_carries_both_compositions() {
    let components = vec!["methane".to_string(), "ethane".to_string()];
    let result = hem_release(
        &components,
        &[0.9, 0.1],
        kelvins(300.0),
        pascals(500_000.0),
        meters(0.01),
        0.62,
        pascals(101_325.0),
    )
    .expect("the state runs");
    assert_eq!(result.status, "valid");
    assert_eq!(result.component_mole_fractions, vec![0.9, 0.1]);
    let total: f64 = result.component_mass_fractions.iter().sum();
    assert!((total - 1.0).abs() <= 1.0e-12);
    assert_eq!(result.component_mass_fractions.len(), 2);
}
