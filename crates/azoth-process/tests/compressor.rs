//! `process.compressor`'s power and its chart-type refusals.
//!
//! The chart's own rows and the chart branch are held by `compressor_chart.rs`'s tests and
//! the model case; this is the isentropic power the record now publishes and the named
//! refusals the spec declares.

use azoth_core::units::{cubic_meters_per_hour, kelvins, meters, pascals};
use azoth_process::models::compressor::compressor;

#[test]
fn the_isentropic_power_is_the_enthalpy_step() {
    let result = compressor(
        &["methane".to_string(), "n-butane".to_string()],
        1.0,
        &[0.9, 0.1],
        pascals(3_000_000.0),
        kelvins(320.0),
        Some(pascals(6_000_000.0)),
        0.75,
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
    .expect("the captured state compresses");

    // `process_compressor.tsv`'s `power_kW=2.4188668614585818`.
    let expected = 2418.8668614585818;
    let error = (result.power.value - expected).abs() / expected;
    assert!(
        error < 1.0e-5,
        "power: {} against {expected}, relative error {error:.3e}",
        result.power.value
    );
    assert_eq!(result.polytropic_head.value, 0.0);
    assert_eq!(result.polytropic_efficiency, 0.0);
}

#[test]
fn the_khader_chart_needs_its_conditions() {
    // Every chart type is carried; the Khader one is the only map that cannot be built
    // without the reference fluid's chart conditions.
    let curve = [
        azoth_process::models::compressor::CompressorCurvePoint {
            speed: 10_000.0,
            flow: cubic_meters_per_hour(2000.0),
            head: meters(4000.0),
            efficiency: 0.8,
        },
        azoth_process::models::compressor::CompressorCurvePoint {
            speed: 10_000.0,
            flow: cubic_meters_per_hour(2500.0),
            head: meters(3800.0),
            efficiency: 0.79,
        },
        azoth_process::models::compressor::CompressorCurvePoint {
            speed: 11_000.0,
            flow: cubic_meters_per_hour(2400.0),
            head: meters(4600.0),
            efficiency: 0.81,
        },
    ];
    let error = compressor(
        &["methane".to_string()],
        1.0,
        &[1.0],
        pascals(3_000_000.0),
        kelvins(320.0),
        None,
        0.75,
        Some(&curve),
        Some(10_000.0),
        Some("khader_2015"),
        None,
        None,
        None,
        None,
        None,
        None,
    )
    .expect_err("the Khader chart needs chart_conditions");
    assert_eq!(error.field(), Some("chart_conditions"));
    for chart_type in ["interpolate", "interpolate_and_extrapolate"] {
        // No curve on this state, so the isentropic branch runs; the type itself is accepted.
        assert!(
            compressor(
                &["methane".to_string()],
                1.0,
                &[1.0],
                pascals(3_000_000.0),
                kelvins(320.0),
                Some(pascals(6_000_000.0)),
                0.75,
                None,
                None,
                Some(chart_type),
                None,
                None,
                None,
                None,
                None,
                None,
            )
            .is_ok(),
            "{chart_type} is carried"
        );
    }
}

/// `process_anti_surge.tsv`'s machine rows, held to the model directly.
///
/// The cross-implementation test compares the two kernels with each other and the model
/// case compares one backend against the recorded value; these pin the Rust model's own
/// anti-surge surface to the capture, including the recycle block the case cannot fully
/// state.
mod anti_surge {
    use super::*;
    use azoth_core::units::{cubic_meters_per_hour, meters};
    use azoth_process::models::compressor::{
        CompressorCurvePoint, CompressorResult, CompressorSurgePoint,
    };

    fn surge_curve() -> Vec<CompressorSurgePoint> {
        vec![
            CompressorSurgePoint {
                flow: cubic_meters_per_hour(2789.0),
                head: meters(8154.943934760448),
            },
            CompressorSurgePoint {
                flow: cubic_meters_per_hour(2550.0),
                head: meters(7339.4495412844035),
            },
            CompressorSurgePoint {
                flow: cubic_meters_per_hour(2500.0),
                head: meters(7135.575942915392),
            },
            CompressorSurgePoint {
                flow: cubic_meters_per_hour(2200.0),
                head: meters(6625.891946992864),
            },
        ]
    }

    fn machine(n: f64, margin: Option<f64>, active: bool, factor: Option<f64>) -> CompressorResult {
        compressor(
            &[
                "nitrogen".to_string(),
                "CO2".to_string(),
                "methane".to_string(),
                "ethane".to_string(),
                "propane".to_string(),
                "i-butane".to_string(),
                "n-butane".to_string(),
                "i-pentane".to_string(),
                "n-pentane".to_string(),
            ],
            n,
            &[
                0.01208165393330526,
                0.013435200224588424,
                0.8820509735507028,
                0.05271812147827308,
                0.03291624055024163,
                0.0008221540435942168,
                0.004882792917443703,
                0.0005614710541619043,
                0.000531392247688945,
            ],
            pascals(4_800_000.0),
            kelvins(297.15),
            None,
            1.0,
            Some(&curve()),
            Some(11918.0),
            None,
            None,
            None,
            Some(&surge_curve()),
            margin,
            Some(active),
            factor,
        )
        .expect("the captured machine runs")
    }

    /// The eight-speed `CompressorChartTest` table, in metres, as the case states it.
    fn curve() -> Vec<azoth_process::models::compressor::CompressorCurvePoint> {
        vec![
            CompressorCurvePoint {
                speed: 12913.0,
                flow: cubic_meters_per_hour(2789.1285),
                head: meters(8158.766564729867),
                efficiency: 0.7724522384095729,
            },
            CompressorCurvePoint {
                speed: 12913.0,
                flow: cubic_meters_per_hour(3174.0375),
                head: meters(8042.1406727828735),
                efficiency: 0.794154186459363,
            },
            CompressorCurvePoint {
                speed: 12913.0,
                flow: cubic_meters_per_hour(3689.2288),
                head: meters(7769.031600407748),
                efficiency: 0.80737960012489,
            },
            CompressorCurvePoint {
                speed: 12913.0,
                flow: cubic_meters_per_hour(4179.4503),
                head: meters(7325.9734964322115),
                efficiency: 0.805229826589649,
            },
            CompressorCurvePoint {
                speed: 12913.0,
                flow: cubic_meters_per_hour(4570.2768),
                head: meters(6830.39755351682),
                efficiency: 0.792210931638144,
            },
            CompressorCurvePoint {
                speed: 12913.0,
                flow: cubic_meters_per_hour(4954.7728),
                head: meters(6177.991845056064),
                efficiency: 0.754719133864634,
            },
            CompressorCurvePoint {
                speed: 12913.0,
                flow: cubic_meters_per_hour(5246.0329),
                head: meters(5407.737003058104),
                efficiency: 0.696034181197298,
            },
            CompressorCurvePoint {
                speed: 12913.0,
                flow: cubic_meters_per_hour(5661.0331),
                head: meters(4049.7451580020384),
                efficiency: 0.587322388482707,
            },
            CompressorCurvePoint {
                speed: 12298.0,
                flow: cubic_meters_per_hour(2571.1753),
                head: meters(7361.0805300713555),
                efficiency: 0.770107837113504,
            },
            CompressorCurvePoint {
                speed: 12298.0,
                flow: cubic_meters_per_hour(2943.7254),
                head: meters(7322.82364933741),
                efficiency: 0.7930699741363889,
            },
            CompressorCurvePoint {
                speed: 12298.0,
                flow: cubic_meters_per_hour(3440.2675),
                head: meters(7023.537206931702),
                efficiency: 0.8089411890211351,
            },
            CompressorCurvePoint {
                speed: 12298.0,
                flow: cubic_meters_per_hour(3837.4448),
                head: meters(6710.917431192661),
                efficiency: 0.807190194665918,
            },
            CompressorCurvePoint {
                speed: 12298.0,
                flow: cubic_meters_per_hour(4253.0898),
                head: meters(6189.266055045871),
                efficiency: 0.795313242980328,
            },
            CompressorCurvePoint {
                speed: 12298.0,
                flow: cubic_meters_per_hour(4668.6643),
                head: meters(5576.146788990825),
                efficiency: 0.7559126228963671,
            },
            CompressorCurvePoint {
                speed: 12298.0,
                flow: cubic_meters_per_hour(4997.1926),
                head: meters(4819.051987767584),
                efficiency: 0.696846136362097,
            },
            CompressorCurvePoint {
                speed: 12298.0,
                flow: cubic_meters_per_hour(5387.4952),
                head: meters(3643.9449541284407),
                efficiency: 0.600043057990909,
            },
            CompressorCurvePoint {
                speed: 11683.0,
                flow: cubic_meters_per_hour(2415.3793),
                head: meters(6641.957186544343),
                efficiency: 0.770043065299874,
            },
            CompressorCurvePoint {
                speed: 11683.0,
                flow: cubic_meters_per_hour(2763.0706),
                head: meters(6577.50254841998),
                efficiency: 0.791690958847856,
            },
            CompressorCurvePoint {
                speed: 11683.0,
                flow: cubic_meters_per_hour(3141.7095),
                head: meters(6382.446483180428),
                efficiency: 0.8080381699756749,
            },
            CompressorCurvePoint {
                speed: 11683.0,
                flow: cubic_meters_per_hour(3594.7436),
                head: meters(6030.774719673802),
                efficiency: 0.806543975614197,
            },
            CompressorCurvePoint {
                speed: 11683.0,
                flow: cubic_meters_per_hour(4047.6467),
                head: meters(5509.225280326197),
                efficiency: 0.7885323891027051,
            },
            CompressorCurvePoint {
                speed: 11683.0,
                flow: cubic_meters_per_hour(4494.1889),
                head: meters(4791.630988786951),
                efficiency: 0.7366647742706129,
            },
            CompressorCurvePoint {
                speed: 11683.0,
                flow: cubic_meters_per_hour(4853.7353),
                head: meters(3995.412844036697),
                efficiency: 0.662735600426727,
            },
            CompressorCurvePoint {
                speed: 11683.0,
                flow: cubic_meters_per_hour(5138.7858),
                head: meters(3225.1478083588177),
                efficiency: 0.57671664571658,
            },
            CompressorCurvePoint {
                speed: 11098.0,
                flow: cubic_meters_per_hour(2247.2043),
                head: meters(5975.06625891947),
                efficiency: 0.770716623789093,
            },
            CompressorCurvePoint {
                speed: 11098.0,
                flow: cubic_meters_per_hour(2799.7342),
                head: meters(5806.595310907237),
                efficiency: 0.804629750233093,
            },
            CompressorCurvePoint {
                speed: 11098.0,
                flow: cubic_meters_per_hour(3178.3428),
                head: meters(5572.344546381244),
                efficiency: 0.8113908111690721,
            },
            CompressorCurvePoint {
                speed: 11098.0,
                flow: cubic_meters_per_hour(3656.1551),
                head: meters(5142.324159021407),
                efficiency: 0.796374242667478,
            },
            CompressorCurvePoint {
                speed: 11098.0,
                flow: cubic_meters_per_hour(4102.778),
                head: meters(4529.2762487257905),
                efficiency: 0.75380928428817,
            },
            CompressorCurvePoint {
                speed: 11098.0,
                flow: cubic_meters_per_hour(4394.1591),
                head: meters(3915.840978593272),
                efficiency: 0.695332969549779,
            },
            CompressorCurvePoint {
                speed: 11098.0,
                flow: cubic_meters_per_hour(4648.3224),
                head: meters(3354.576962283384),
                efficiency: 0.637997587622339,
            },
            CompressorCurvePoint {
                speed: 11098.0,
                flow: cubic_meters_per_hour(4840.4998),
                head: meters(2936.890927624873),
                efficiency: 0.588120614497758,
            },
            CompressorCurvePoint {
                speed: 10453.0,
                flow: cubic_meters_per_hour(2072.8397),
                head: meters(5334.301732925586),
                efficiency: 0.769705872525642,
            },
            CompressorCurvePoint {
                speed: 10453.0,
                flow: cubic_meters_per_hour(2463.9483),
                head: meters(5204.617737003057),
                efficiency: 0.798335492585324,
            },
            CompressorCurvePoint {
                speed: 10453.0,
                flow: cubic_meters_per_hour(2836.4078),
                head: meters(5048.756371049949),
                efficiency: 0.809468133671171,
            },
            CompressorCurvePoint {
                speed: 10453.0,
                flow: cubic_meters_per_hour(3202.5266),
                head: meters(4722.9969418960245),
                efficiency: 0.805806471927835,
            },
            CompressorCurvePoint {
                speed: 10453.0,
                flow: cubic_meters_per_hour(3599.6333),
                head: meters(4318.909276248725),
                efficiency: 0.7804621582254261,
            },
            CompressorCurvePoint {
                speed: 10453.0,
                flow: cubic_meters_per_hour(3978.0203),
                head: meters(3797.16615698267),
                efficiency: 0.730403707523258,
            },
            CompressorCurvePoint {
                speed: 10453.0,
                flow: cubic_meters_per_hour(4257.0022),
                head: meters(3209.826707441386),
                efficiency: 0.665572286338589,
            },
            CompressorCurvePoint {
                speed: 10453.0,
                flow: cubic_meters_per_hour(4517.345),
                head: meters(2609.3781855249745),
                efficiency: 0.598624822515064,
            },
            CompressorCurvePoint {
                speed: 9224.0,
                flow: cubic_meters_per_hour(1835.9552),
                head: meters(4144.525993883793),
                efficiency: 0.775063036680357,
            },
            CompressorCurvePoint {
                speed: 9224.0,
                flow: cubic_meters_per_hour(2208.455),
                head: meters(4040.937818552497),
                efficiency: 0.802056198362559,
            },
            CompressorCurvePoint {
                speed: 9224.0,
                flow: cubic_meters_per_hour(2618.1322),
                head: meters(3832.9051987767584),
                efficiency: 0.810339108025933,
            },
            CompressorCurvePoint {
                speed: 9224.0,
                flow: cubic_meters_per_hour(2940.8034),
                head: meters(3533.1600407747196),
                efficiency: 0.7960859626879391,
            },
            CompressorCurvePoint {
                speed: 9224.0,
                flow: cubic_meters_per_hour(3244.7852),
                head: meters(3154.974515800204),
                efficiency: 0.7638145344044049,
            },
            CompressorCurvePoint {
                speed: 9224.0,
                flow: cubic_meters_per_hour(3530.1279),
                head: meters(2763.669724770642),
                efficiency: 0.708027503005902,
            },
            CompressorCurvePoint {
                speed: 9224.0,
                flow: cubic_meters_per_hour(3753.3738),
                head: meters(2372.201834862385),
                efficiency: 0.6464373671605711,
            },
            CompressorCurvePoint {
                speed: 9224.0,
                flow: cubic_meters_per_hour(3895.9746),
                head: meters(2085.076452599388),
                efficiency: 0.6052993499823419,
            },
            CompressorCurvePoint {
                speed: 8609.0,
                flow: cubic_meters_per_hour(1711.3386),
                head: meters(3595.3618756371047),
                efficiency: 0.778175271586685,
            },
            CompressorCurvePoint {
                speed: 8609.0,
                flow: cubic_meters_per_hour(1965.8848),
                head: meters(3530.6727828746175),
                efficiency: 0.80065165942218,
            },
            CompressorCurvePoint {
                speed: 8609.0,
                flow: cubic_meters_per_hour(2356.9431),
                head: meters(3335.65749235474),
                efficiency: 0.810631362122632,
            },
            CompressorCurvePoint {
                speed: 8609.0,
                flow: cubic_meters_per_hour(2685.9247),
                head: meters(3166.615698267074),
                efficiency: 0.7989550517712991,
            },
            CompressorCurvePoint {
                speed: 8609.0,
                flow: cubic_meters_per_hour(3008.5154),
                head: meters(2762.334352701325),
                efficiency: 0.761983240929369,
            },
            CompressorCurvePoint {
                speed: 8609.0,
                flow: cubic_meters_per_hour(3337.2855),
                head: meters(2318.8786952089704),
                efficiency: 0.69289982774309,
            },
            CompressorCurvePoint {
                speed: 8609.0,
                flow: cubic_meters_per_hour(3591.5092),
                head: meters(1836.0142711518856),
                efficiency: 0.608567149372229,
            },
            CompressorCurvePoint {
                speed: 8200.0,
                flow: cubic_meters_per_hour(1636.5807),
                head: meters(3281.549439347604),
                efficiency: 0.780924334304045,
            },
            CompressorCurvePoint {
                speed: 8200.0,
                flow: cubic_meters_per_hour(2002.8708),
                head: meters(3177.9408766564725),
                efficiency: 0.8093535515686671,
            },
            CompressorCurvePoint {
                speed: 8200.0,
                flow: cubic_meters_per_hour(2338.0319),
                head: meters(2969.714576962283),
                efficiency: 0.8079044377662341,
            },
            CompressorCurvePoint {
                speed: 8200.0,
                flow: cubic_meters_per_hour(2642.1245),
                head: meters(2735.2701325178386),
                efficiency: 0.788639325223295,
            },
            CompressorCurvePoint {
                speed: 8200.0,
                flow: cubic_meters_per_hour(2896.4894),
                head: meters(2435.3618756371047),
                efficiency: 0.752170936751143,
            },
            CompressorCurvePoint {
                speed: 8200.0,
                flow: cubic_meters_per_hour(3113.6264),
                head: meters(2174.5565749235475),
                efficiency: 0.703105081673411,
            },
            CompressorCurvePoint {
                speed: 8200.0,
                flow: cubic_meters_per_hour(3274.8764),
                head: meters(1913.61875637105),
                efficiency: 0.6555075685335691,
            },
            CompressorCurvePoint {
                speed: 8200.0,
                flow: cubic_meters_per_hour(3411.2977),
                head: meters(1665.6778797145769),
                efficiency: 0.610391468300337,
            },
        ]
    }

    fn close(actual: f64, expected: f64, tolerance: f64, what: &str) {
        let error = (actual - expected).abs() / expected.abs().max(f64::MIN_POSITIVE);
        assert!(
            error <= tolerance,
            "{what}: {actual} against {expected}, relative error {error:.3e} > {tolerance:.3e}"
        );
    }

    #[test]
    fn the_read_outs_match_the_capture() {
        let result = machine(489.4970603565765, Some(0.10), false, None);
        close(
            result.surge_flow_rate.expect("a curve is set"),
            1_615.282_214_403_193 / 3600.0,
            0.03,
            "surge_flow_rate",
        );
        close(
            result.distance_to_surge.expect("a curve is set"),
            -0.504_410_127_751_431_7,
            0.03,
            "distance_to_surge",
        );
        close(
            result.control_line_flow.expect("a curve is set"),
            1_776.810_435_843_512_3 / 3600.0,
            0.03,
            "control_line_flow",
        );
        close(
            result.required_recycle_fraction.expect("a curve is set"),
            0.549_463_752_501_301_5,
            0.03,
            "required_recycle_fraction",
        );
        assert_eq!(result.recycle_fraction, None, "the block was not active");
        assert_eq!(result.anti_surge_recycle_power, Some(0.0));
    }

    #[test]
    fn the_recycle_block_holds_the_machine_off_surge() {
        let result = machine(146.84911810697295, None, true, Some(1.05));
        close(
            result.recycle_fraction.expect("the block was active"),
            9.065_402_690_878_43,
            0.07,
            "recycle_fraction",
        );
        close(result.power.value, 2_411_274.262_873_462_4, 0.03, "power");
        close(
            result.surge_flow_rate.expect("a curve is set"),
            20_729.726_580_828_03 / 3600.0,
            0.07,
            "surge_flow_rate",
        );
    }

    #[test]
    fn the_machine_well_above_surge_recycles_nothing() {
        let result = machine(1957.9882414263066, None, true, Some(1.05));
        assert_eq!(result.recycle_fraction, Some(0.0));
        close(
            result.distance_to_surge.expect("a curve is set"),
            0.455_804_891_686_615_53,
            0.16,
            "distance_to_surge",
        );
        close(result.power.value, 2_891_032.446_677_318_3, 0.03, "power");
    }

    #[test]
    fn the_control_factor_is_the_recycle_target_s_multiple() {
        let lower = machine(146.84911810697295, None, true, Some(1.0));
        let nominal = machine(146.84911810697295, None, true, Some(1.05));
        let lower = lower.recycle_fraction.expect("the block was active");
        let nominal = nominal.recycle_fraction.expect("the block was active");
        assert!(lower > 0.0 && nominal > 0.0);
        assert!(
            nominal > lower,
            "a larger factor targets a larger flow: {nominal} against {lower}"
        );
    }

    #[test]
    fn a_negative_control_margin_is_refused() {
        let error = compressor(
            &[
                "nitrogen".to_string(),
                "CO2".to_string(),
                "methane".to_string(),
                "ethane".to_string(),
                "propane".to_string(),
                "i-butane".to_string(),
                "n-butane".to_string(),
                "i-pentane".to_string(),
                "n-pentane".to_string(),
            ],
            1.0,
            &[
                0.01208165393330526,
                0.013435200224588424,
                0.8820509735507028,
                0.05271812147827308,
                0.03291624055024163,
                0.0008221540435942168,
                0.004882792917443703,
                0.0005614710541619043,
                0.000531392247688945,
            ],
            pascals(4_800_000.0),
            kelvins(297.15),
            None,
            1.0,
            Some(&curve()),
            Some(11918.0),
            None,
            None,
            None,
            Some(&surge_curve()),
            Some(-0.1),
            None,
            None,
        )
        .expect_err("the spec refuses a negative margin");
        assert_eq!(error.field(), Some("surge_control_margin"));
    }
}
