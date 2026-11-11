//! `process.map_turbo_expander_compressor` against `validation/neqsim/captures/process_map_turbo_expander.tsv`.

use azoth_core::units::{cubic_meters_per_hour, kelvins, meters, pascals};
use azoth_process::models::map_turbo_expander_compressor::{
    MapCurvePoint, MapTurboExpanderResult, map_turbo_expander_compressor,
};

fn curve() -> Vec<MapCurvePoint> {
    vec![
        MapCurvePoint {
            speed: 12913.0,
            flow: cubic_meters_per_hour(2789.1285),
            head: meters(8158.766564729867),
            efficiency: 0.7724522384095729,
        },
        MapCurvePoint {
            speed: 12913.0,
            flow: cubic_meters_per_hour(3174.0375),
            head: meters(8042.1406727828735),
            efficiency: 0.794154186459363,
        },
        MapCurvePoint {
            speed: 12913.0,
            flow: cubic_meters_per_hour(3689.2288),
            head: meters(7769.031600407748),
            efficiency: 0.80737960012489,
        },
        MapCurvePoint {
            speed: 12913.0,
            flow: cubic_meters_per_hour(4179.4503),
            head: meters(7325.9734964322115),
            efficiency: 0.805229826589649,
        },
        MapCurvePoint {
            speed: 12913.0,
            flow: cubic_meters_per_hour(4570.2768),
            head: meters(6830.39755351682),
            efficiency: 0.792210931638144,
        },
        MapCurvePoint {
            speed: 12913.0,
            flow: cubic_meters_per_hour(4954.7728),
            head: meters(6177.991845056064),
            efficiency: 0.754719133864634,
        },
        MapCurvePoint {
            speed: 12913.0,
            flow: cubic_meters_per_hour(5246.0329),
            head: meters(5407.737003058104),
            efficiency: 0.696034181197298,
        },
        MapCurvePoint {
            speed: 12913.0,
            flow: cubic_meters_per_hour(5661.0331),
            head: meters(4049.7451580020384),
            efficiency: 0.587322388482707,
        },
        MapCurvePoint {
            speed: 12298.0,
            flow: cubic_meters_per_hour(2571.1753),
            head: meters(7361.0805300713555),
            efficiency: 0.770107837113504,
        },
        MapCurvePoint {
            speed: 12298.0,
            flow: cubic_meters_per_hour(2943.7254),
            head: meters(7322.82364933741),
            efficiency: 0.7930699741363889,
        },
        MapCurvePoint {
            speed: 12298.0,
            flow: cubic_meters_per_hour(3440.2675),
            head: meters(7023.537206931702),
            efficiency: 0.8089411890211351,
        },
        MapCurvePoint {
            speed: 12298.0,
            flow: cubic_meters_per_hour(3837.4448),
            head: meters(6710.917431192661),
            efficiency: 0.807190194665918,
        },
        MapCurvePoint {
            speed: 12298.0,
            flow: cubic_meters_per_hour(4253.0898),
            head: meters(6189.266055045871),
            efficiency: 0.795313242980328,
        },
        MapCurvePoint {
            speed: 12298.0,
            flow: cubic_meters_per_hour(4668.6643),
            head: meters(5576.146788990825),
            efficiency: 0.7559126228963671,
        },
        MapCurvePoint {
            speed: 12298.0,
            flow: cubic_meters_per_hour(4997.1926),
            head: meters(4819.051987767584),
            efficiency: 0.696846136362097,
        },
        MapCurvePoint {
            speed: 12298.0,
            flow: cubic_meters_per_hour(5387.4952),
            head: meters(3643.9449541284407),
            efficiency: 0.600043057990909,
        },
        MapCurvePoint {
            speed: 11683.0,
            flow: cubic_meters_per_hour(2415.3793),
            head: meters(6641.957186544343),
            efficiency: 0.770043065299874,
        },
        MapCurvePoint {
            speed: 11683.0,
            flow: cubic_meters_per_hour(2763.0706),
            head: meters(6577.50254841998),
            efficiency: 0.791690958847856,
        },
        MapCurvePoint {
            speed: 11683.0,
            flow: cubic_meters_per_hour(3141.7095),
            head: meters(6382.446483180428),
            efficiency: 0.8080381699756749,
        },
        MapCurvePoint {
            speed: 11683.0,
            flow: cubic_meters_per_hour(3594.7436),
            head: meters(6030.774719673802),
            efficiency: 0.806543975614197,
        },
        MapCurvePoint {
            speed: 11683.0,
            flow: cubic_meters_per_hour(4047.6467),
            head: meters(5509.225280326197),
            efficiency: 0.7885323891027051,
        },
        MapCurvePoint {
            speed: 11683.0,
            flow: cubic_meters_per_hour(4494.1889),
            head: meters(4791.630988786951),
            efficiency: 0.7366647742706129,
        },
        MapCurvePoint {
            speed: 11683.0,
            flow: cubic_meters_per_hour(4853.7353),
            head: meters(3995.412844036697),
            efficiency: 0.662735600426727,
        },
        MapCurvePoint {
            speed: 11683.0,
            flow: cubic_meters_per_hour(5138.7858),
            head: meters(3225.1478083588177),
            efficiency: 0.57671664571658,
        },
        MapCurvePoint {
            speed: 11098.0,
            flow: cubic_meters_per_hour(2247.2043),
            head: meters(5975.06625891947),
            efficiency: 0.770716623789093,
        },
        MapCurvePoint {
            speed: 11098.0,
            flow: cubic_meters_per_hour(2799.7342),
            head: meters(5806.595310907237),
            efficiency: 0.804629750233093,
        },
        MapCurvePoint {
            speed: 11098.0,
            flow: cubic_meters_per_hour(3178.3428),
            head: meters(5572.344546381244),
            efficiency: 0.8113908111690721,
        },
        MapCurvePoint {
            speed: 11098.0,
            flow: cubic_meters_per_hour(3656.1551),
            head: meters(5142.324159021407),
            efficiency: 0.796374242667478,
        },
        MapCurvePoint {
            speed: 11098.0,
            flow: cubic_meters_per_hour(4102.778),
            head: meters(4529.2762487257905),
            efficiency: 0.75380928428817,
        },
        MapCurvePoint {
            speed: 11098.0,
            flow: cubic_meters_per_hour(4394.1591),
            head: meters(3915.840978593272),
            efficiency: 0.695332969549779,
        },
        MapCurvePoint {
            speed: 11098.0,
            flow: cubic_meters_per_hour(4648.3224),
            head: meters(3354.576962283384),
            efficiency: 0.637997587622339,
        },
        MapCurvePoint {
            speed: 11098.0,
            flow: cubic_meters_per_hour(4840.4998),
            head: meters(2936.890927624873),
            efficiency: 0.588120614497758,
        },
        MapCurvePoint {
            speed: 10453.0,
            flow: cubic_meters_per_hour(2072.8397),
            head: meters(5334.301732925586),
            efficiency: 0.769705872525642,
        },
        MapCurvePoint {
            speed: 10453.0,
            flow: cubic_meters_per_hour(2463.9483),
            head: meters(5204.617737003057),
            efficiency: 0.798335492585324,
        },
        MapCurvePoint {
            speed: 10453.0,
            flow: cubic_meters_per_hour(2836.4078),
            head: meters(5048.756371049949),
            efficiency: 0.809468133671171,
        },
        MapCurvePoint {
            speed: 10453.0,
            flow: cubic_meters_per_hour(3202.5266),
            head: meters(4722.9969418960245),
            efficiency: 0.805806471927835,
        },
        MapCurvePoint {
            speed: 10453.0,
            flow: cubic_meters_per_hour(3599.6333),
            head: meters(4318.909276248725),
            efficiency: 0.7804621582254261,
        },
        MapCurvePoint {
            speed: 10453.0,
            flow: cubic_meters_per_hour(3978.0203),
            head: meters(3797.16615698267),
            efficiency: 0.730403707523258,
        },
        MapCurvePoint {
            speed: 10453.0,
            flow: cubic_meters_per_hour(4257.0022),
            head: meters(3209.826707441386),
            efficiency: 0.665572286338589,
        },
        MapCurvePoint {
            speed: 10453.0,
            flow: cubic_meters_per_hour(4517.345),
            head: meters(2609.3781855249745),
            efficiency: 0.598624822515064,
        },
        MapCurvePoint {
            speed: 9224.0,
            flow: cubic_meters_per_hour(1835.9552),
            head: meters(4144.525993883793),
            efficiency: 0.775063036680357,
        },
        MapCurvePoint {
            speed: 9224.0,
            flow: cubic_meters_per_hour(2208.455),
            head: meters(4040.937818552497),
            efficiency: 0.802056198362559,
        },
        MapCurvePoint {
            speed: 9224.0,
            flow: cubic_meters_per_hour(2618.1322),
            head: meters(3832.9051987767584),
            efficiency: 0.810339108025933,
        },
        MapCurvePoint {
            speed: 9224.0,
            flow: cubic_meters_per_hour(2940.8034),
            head: meters(3533.1600407747196),
            efficiency: 0.7960859626879391,
        },
        MapCurvePoint {
            speed: 9224.0,
            flow: cubic_meters_per_hour(3244.7852),
            head: meters(3154.974515800204),
            efficiency: 0.7638145344044049,
        },
        MapCurvePoint {
            speed: 9224.0,
            flow: cubic_meters_per_hour(3530.1279),
            head: meters(2763.669724770642),
            efficiency: 0.708027503005902,
        },
        MapCurvePoint {
            speed: 9224.0,
            flow: cubic_meters_per_hour(3753.3738),
            head: meters(2372.201834862385),
            efficiency: 0.6464373671605711,
        },
        MapCurvePoint {
            speed: 9224.0,
            flow: cubic_meters_per_hour(3895.9746),
            head: meters(2085.076452599388),
            efficiency: 0.6052993499823419,
        },
        MapCurvePoint {
            speed: 8609.0,
            flow: cubic_meters_per_hour(1711.3386),
            head: meters(3595.3618756371047),
            efficiency: 0.778175271586685,
        },
        MapCurvePoint {
            speed: 8609.0,
            flow: cubic_meters_per_hour(1965.8848),
            head: meters(3530.6727828746175),
            efficiency: 0.80065165942218,
        },
        MapCurvePoint {
            speed: 8609.0,
            flow: cubic_meters_per_hour(2356.9431),
            head: meters(3335.65749235474),
            efficiency: 0.810631362122632,
        },
        MapCurvePoint {
            speed: 8609.0,
            flow: cubic_meters_per_hour(2685.9247),
            head: meters(3166.615698267074),
            efficiency: 0.7989550517712991,
        },
        MapCurvePoint {
            speed: 8609.0,
            flow: cubic_meters_per_hour(3008.5154),
            head: meters(2762.334352701325),
            efficiency: 0.761983240929369,
        },
        MapCurvePoint {
            speed: 8609.0,
            flow: cubic_meters_per_hour(3337.2855),
            head: meters(2318.8786952089704),
            efficiency: 0.69289982774309,
        },
        MapCurvePoint {
            speed: 8609.0,
            flow: cubic_meters_per_hour(3591.5092),
            head: meters(1836.0142711518856),
            efficiency: 0.608567149372229,
        },
        MapCurvePoint {
            speed: 8200.0,
            flow: cubic_meters_per_hour(1636.5807),
            head: meters(3281.549439347604),
            efficiency: 0.780924334304045,
        },
        MapCurvePoint {
            speed: 8200.0,
            flow: cubic_meters_per_hour(2002.8708),
            head: meters(3177.9408766564725),
            efficiency: 0.8093535515686671,
        },
        MapCurvePoint {
            speed: 8200.0,
            flow: cubic_meters_per_hour(2338.0319),
            head: meters(2969.714576962283),
            efficiency: 0.8079044377662341,
        },
        MapCurvePoint {
            speed: 8200.0,
            flow: cubic_meters_per_hour(2642.1245),
            head: meters(2735.2701325178386),
            efficiency: 0.788639325223295,
        },
        MapCurvePoint {
            speed: 8200.0,
            flow: cubic_meters_per_hour(2896.4894),
            head: meters(2435.3618756371047),
            efficiency: 0.752170936751143,
        },
        MapCurvePoint {
            speed: 8200.0,
            flow: cubic_meters_per_hour(3113.6264),
            head: meters(2174.5565749235475),
            efficiency: 0.703105081673411,
        },
        MapCurvePoint {
            speed: 8200.0,
            flow: cubic_meters_per_hour(3274.8764),
            head: meters(1913.61875637105),
            efficiency: 0.6555075685335691,
        },
        MapCurvePoint {
            speed: 8200.0,
            flow: cubic_meters_per_hour(3411.2977),
            head: meters(1665.6778797145769),
            efficiency: 0.610391468300337,
        },
    ]
}

#[allow(clippy::too_many_arguments)]
fn machine(expander_bara: f64, mode: &str, chart: bool) -> MapTurboExpanderResult {
    let points = curve();
    let conditions = [25.0, 50.0, 50.0, 20.0];
    let (curve_arg, chart_type, chart_conditions, diameter) = if chart {
        (
            Some(&points[..]),
            Some("khader_2015"),
            Some(&conditions[..]),
            Some(meters(0.9)),
        )
    } else {
        (None, None, None, None)
    };
    map_turbo_expander_compressor(
        &[
            "nitrogen".to_string(),
            "CO2".to_string(),
            "methane".to_string(),
            "ethane".to_string(),
            "propane".to_string(),
            "i-butane".to_string(),
            "n-butane".to_string(),
        ],
        1468.4911810697297,
        &[
            0.012094871975027354,
            0.01344989912575656,
            0.8830159893203787,
            0.05277579821136417,
            0.03295225285810357,
            0.000823053528591073,
            0.004888134980778689,
        ],
        pascals(expander_bara * 1.0e5),
        kelvins(298.15),
        &[
            "nitrogen".to_string(),
            "CO2".to_string(),
            "methane".to_string(),
            "ethane".to_string(),
            "propane".to_string(),
            "i-butane".to_string(),
            "n-butane".to_string(),
        ],
        1468.4911810697297,
        &[
            0.012094871975027354,
            0.01344989912575656,
            0.8830159893203787,
            0.05277579821136417,
            0.03295225285810357,
            0.000823053528591073,
            0.004888134980778689,
        ],
        pascals(4_800_000.0),
        kelvins(293.15),
        pascals(3_000_000.0),
        0.85,
        Some(mode),
        Some(0.98),
        None,
        Some(pascals(5_800_000.0)),
        Some(0.78),
        None,
        None,
        None,
        curve_arg,
        chart_type,
        chart_conditions,
        diameter,
    )
    .expect("the captured machine runs")
}

fn close(actual: f64, expected: f64, tolerance: f64, what: &str) {
    let error = (actual - expected).abs() / expected.abs().max(f64::MIN_POSITIVE);
    assert!(
        error <= tolerance,
        "{what}: {actual} against {expected}, relative error {error:.3e} > {tolerance:.3e}"
    );
}

#[test]
fn the_balanced_rows_are_the_capture_s() {
    let balanced = machine(60.0, "balanced_speed", true);
    close(
        balanced.shaft_speed,
        10418.063973844051,
        0.03,
        "shaft_speed",
    );
    close(
        balanced.available_shaft_power.value,
        1_638_745.031_561_063_1,
        0.03,
        "available_shaft_power",
    );
    close(
        balanced.consumed_compressor_power.value,
        1_638_745.132_115_545_7,
        0.03,
        "consumed_compressor_power",
    );
    assert!(!balanced.speed_limited);
    assert_eq!(balanced.operating_status, "balanced");
    // **The identity is exact by construction**: the residual is what the two report.
    close(
        balanced.power_balance_residual.value,
        balanced.available_shaft_power.value - balanced.consumed_compressor_power.value,
        1.0e-12,
        "the residual identity",
    );

    let turndown = machine(31.0, "balanced_speed", true);
    assert!(turndown.speed_limited);
    assert_eq!(turndown.operating_status, "under_power_surge");
    close(turndown.shaft_speed, 8200.0, 1.0e-9, "the pegged speed");
    close(
        turndown.available_shaft_power.value,
        90_817.129_694_174_33,
        0.03,
        "the turndown's available power",
    );
    assert!(turndown.power_balance_residual.value < 0.0);
}

#[test]
fn the_specified_row_reports_its_residual() {
    let specified = machine(60.0, "specified_pressures", false);
    close(specified.shaft_speed, 3000.0, 1.0e-12, "the reported speed");
    close(
        specified.consumed_compressor_power.value,
        770_477.316_556_879_1,
        0.03,
        "consumed_compressor_power",
    );
    close(
        specified.power_balance_residual.value,
        868_267.715_004_184,
        0.03,
        "power_balance_residual",
    );
    close(
        specified.compressor_outlet_p.value,
        5_800_000.0,
        1.0e-9,
        "the stated pressure",
    );
    assert_eq!(specified.operating_status, "balanced");
    assert!(!specified.speed_limited);
}
