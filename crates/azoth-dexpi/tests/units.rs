//! The reference-unit table converts what it names and refuses what it does not.

use azoth_dexpi::units::{QuantityKind, UnitError, from_si, kind, to_si};

#[test]
fn the_process_reference_units_convert() {
    let kg_per_h = "Core/PhysicalQuantities.MassFlowRateUnit.KilogramPerHour";
    assert!((to_si(kg_per_h, 3600.0).unwrap() - 1.0).abs() < 1e-12);
    assert!((from_si(kg_per_h, 1.0).unwrap() - 3600.0).abs() < 1e-9);

    let bar = "Core/PhysicalQuantities.PressureAbsoluteUnit.Bar";
    assert!((to_si(bar, 40.0).unwrap() - 4.0e6).abs() < 1e-6);
    assert!((from_si(bar, 4.0e6).unwrap() - 40.0).abs() < 1e-9);

    let celsius = "Core/PhysicalQuantities.TemperatureUnit.DegreeCelsius";
    assert!((to_si(celsius, 25.0).unwrap() - 298.15).abs() < 1e-12);
    assert!((from_si(celsius, 298.15).unwrap() - 25.0).abs() < 1e-12);

    assert_eq!(
        kind("Core/PhysicalQuantities.LengthUnit.Millimetre"),
        Some(QuantityKind::Length)
    );
}

#[test]
fn an_unreviewed_unit_is_refused_by_name() {
    let error = to_si(
        "Core/PhysicalQuantities.MassFlowRateUnit.SlugPerFortnight",
        1.0,
    )
    .expect_err("an unreviewed unit refuses");
    assert!(matches!(error, UnitError::Unsupported { .. }));
    assert!(error.message().contains("SlugPerFortnight"));
}

#[test]
fn gauge_pressure_is_refused_rather_than_given_an_atmosphere() {
    let error = to_si("Core/PhysicalQuantities.PressureGaugeUnit.Bar", 1.0)
        .expect_err("gauge pressure refuses");
    assert!(matches!(error, UnitError::Gauge { .. }));
}

#[test]
fn a_non_finite_value_is_refused() {
    let bar = "Core/PhysicalQuantities.PressureAbsoluteUnit.Bar";
    let error = to_si(bar, f64::NAN).expect_err("NaN refuses");
    assert!(matches!(error, UnitError::NonFinite { .. }));
}
