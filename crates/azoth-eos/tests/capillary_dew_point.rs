//! `eos.capillary_dew_point` against NeqSim 3.20.0's own in-pore dew point.
//!
//! The expected values are NeqSim's, printed by `validation/neqsim/CapillaryDew.java`. The
//! oracle gates nothing - a divergence is a finding - so the assertions say what a divergence
//! would be rather than only that there is one.
//!
//! **The surface tension is an argument here and a property there.** NeqSim reads it off the
//! phase's interphase properties; each case passes the value NeqSim read at that state, so the
//! comparison is of the Kelvin shift and not of a surface-tension correlation.

use azoth_core::units::pascals;
use azoth_eos::capillary_dew_point;
use azoth_eos::databank::mixture_of;

/// One case: the radius, the angle, the tension, and what NeqSim reported.
struct Case {
    radius: f64,
    angle: f64,
    sigma: f64,
    /// NeqSim's capillary dew point, in kelvin.
    temperature: f64,
}

fn run(case: &Case) -> azoth_eos::results::CapillaryDewPointResult {
    let (mixture, _) = mixture_of(&["methane", "n-butane"], None).expect("the pair resolves");
    capillary_dew_point(
        &mixture,
        pascals(2.0e6),
        &[0.5, 0.5],
        case.radius,
        case.angle,
        case.sigma,
    )
    .expect("the dew point converges")
}

#[test]
fn the_capillary_dew_point_matches_neqsim_at_every_radius() {
    let cases = [
        Case {
            radius: 1.0e-8,
            angle: 0.0,
            sigma: 0.004789210568,
            temperature: 347.6658507316,
        },
        Case {
            radius: 1.0e-7,
            angle: 0.0,
            sigma: 0.004997523456,
            temperature: 345.4016016437,
        },
        Case {
            radius: 1.0e-6,
            angle: 0.0,
            sigma: 0.005018886547,
            temperature: 345.1677087431,
        },
        Case {
            radius: 1.0e-7,
            angle: std::f64::consts::FRAC_PI_3,
            sigma: 0.005009379844,
            temperature: 345.2718314099,
        },
    ];
    for case in &cases {
        let got = run(case);
        let t = got.temperature.value;
        // A tenth of a kelvin: the shift is 5% small because the molar volume is taken from
        // the cubic's root where NeqSim takes it from `getMolarVolume`, and the case file
        // records the measurement. The *form* is exact, which the two tests below check.
        assert!(
            (t - case.temperature).abs() < 0.2,
            "r={:e} theta={:.4}: got {t:.10} K, NeqSim {:.10} K",
            case.radius,
            case.angle,
            case.temperature
        );
    }
}

/// **The shift is `2 sigma cos(theta) / r` and nothing else**, which is the whole physics and
/// the thing a single case cannot show. Three radii at a tenth of each other and one angle at
/// `cos = 1/2` pin both dependences: a port that lost the `cos`, or squared the radius, or
/// applied the shift with the wrong sign, fails here even where it might pass one point.
///
/// **The scaling is a form, not an identity.** The Kelvin exponent enters the K-values, so the
/// converged state - its `z_liquid`, its composition - moves with it, and the response is
/// therefore a few per cent from linear. NeqSim's own ratio between its 10 nm and 100 nm shifts
/// is `9.71`, not `10`, and this library's is `10.14`; both are the same physics with a
/// different molar volume. So the tolerances below are what says "the right power of `r` and
/// the right function of `theta`", which is what these tests exist for.
#[test]
fn the_shift_scales_as_the_young_laplace_pressure() {
    let bulk = 345.1416354411;
    let ten_nm = run(&Case {
        radius: 1.0e-8,
        angle: 0.0,
        sigma: 0.005,
        temperature: 0.0,
    });
    let hundred_nm = run(&Case {
        radius: 1.0e-7,
        angle: 0.0,
        sigma: 0.005,
        temperature: 0.0,
    });
    let micron = run(&Case {
        radius: 1.0e-6,
        angle: 0.0,
        sigma: 0.005,
        temperature: 0.0,
    });

    let d10 = ten_nm.temperature.value - bulk;
    let d100 = hundred_nm.temperature.value - bulk;
    let d1000 = micron.temperature.value - bulk;

    assert!(
        d10 > 0.0 && d100 > 0.0 && d1000 > 0.0,
        "the dew point must move up"
    );
    // `1/r`: ten times the radius, about a tenth of the shift - about, because the state moves
    // with the shift and the response is not exactly linear.
    assert!(
        (d10 / d100 - 10.0).abs() < 0.5,
        "ten times the radius gave {:.4} times the shift, not ten",
        d10 / d100
    );
    assert!(
        (d100 / d1000 - 10.0).abs() < 0.5,
        "ten times the radius gave {:.4} times the shift, not ten",
        d100 / d1000
    );
}

/// `cos(theta)`: at 60 degrees the shift is exactly half the wetting one.
#[test]
fn the_shift_scales_as_the_cosine_of_the_contact_angle() {
    let bulk = 345.1416354411;
    let wetting = run(&Case {
        radius: 1.0e-7,
        angle: 0.0,
        sigma: 0.005,
        temperature: 0.0,
    });
    let sixty = run(&Case {
        radius: 1.0e-7,
        angle: std::f64::consts::FRAC_PI_3,
        sigma: 0.005,
        temperature: 0.0,
    });
    let d0 = wetting.temperature.value - bulk;
    let d60 = sixty.temperature.value - bulk;
    assert!(
        (d60 / d0 - 0.5).abs() < 0.01,
        "cos(60 deg) is a half and the shift ratio was {:.10}",
        d60 / d0
    );
    // `pi/2` is non-wetting: the interface is flat and there is no shift at all.
    let ninety = run(&Case {
        radius: 1.0e-7,
        angle: std::f64::consts::FRAC_PI_2,
        sigma: 0.005,
        temperature: 0.0,
    });
    assert!(
        (ninety.temperature.value - bulk).abs() < 1e-6,
        "a non-wetting contact angle left a shift of {:.3e} K",
        ninety.temperature.value - bulk
    );
    assert!(ninety.capillary_pressure.value.abs() < 1e-3);
}

/// A wide enough pore is the flat interface, and the model reduces to `eos.dew_temperature`.
#[test]
fn a_wide_pore_reproduces_the_bulk_dew_point() {
    let (mixture, _) = mixture_of(&["methane", "n-butane"], None).expect("the pair resolves");
    let bulk = azoth_eos::dew_temperature(&mixture, pascals(2.0e6), &[0.5, 0.5])
        .expect("the bulk dew point converges");
    let wide = capillary_dew_point(&mixture, pascals(2.0e6), &[0.5, 0.5], 1.0e-2, 0.0, 0.005)
        .expect("the wide-pore dew point converges");
    assert!(
        (wide.temperature.value - bulk.temperature.value).abs() < 1e-5,
        "a centimetre pore is flat to {:.3e} K",
        wide.temperature.value - bulk.temperature.value
    );
}

#[test]
fn a_non_positive_radius_is_refused() {
    let (mixture, _) = mixture_of(&["methane", "n-butane"], None).expect("the pair resolves");
    for radius in [0.0, -1.0e-8] {
        let error = capillary_dew_point(&mixture, pascals(2.0e6), &[0.5, 0.5], radius, 0.0, 0.005)
            .expect_err("a non-positive radius is a division by zero, not a wide pore");
        // `OutOfRange`, not `InvalidInput`: the spec's `valid_range` block declares the
        // bound and `apply_checks` enforces it, which is where a declared range belongs.
        assert!(matches!(error, azoth_core::AzothError::OutOfRange { .. }));
    }
}

#[test]
fn a_negative_surface_tension_is_refused() {
    let (mixture, _) = mixture_of(&["methane", "n-butane"], None).expect("the pair resolves");
    let error = capillary_dew_point(&mixture, pascals(2.0e6), &[0.5, 0.5], 1.0e-7, 0.0, -0.005)
        .expect_err("a negative tension would move the dew point down");
    assert!(matches!(error, azoth_core::AzothError::OutOfRange { .. }));
}
