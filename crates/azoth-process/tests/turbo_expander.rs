//! `process.turbo_expander_compressor` against `validation/neqsim/captures/process_turbo_expander.tsv`.

use azoth_core::units::{kelvins, meters, pascals, square_meters, watts};
use azoth_process::models::turbo_expander_compressor::{
    TurboExpanderResult, TurboQnPoint, TurboUcPoint, turbo_expander_compressor,
};

const COMPONENTS: [&str; 8] = [
    "nitrogen", "CO2", "methane", "ethane", "propane", "i-butane", "n-butane", "n-hexane",
];

const Z: [f64; 8] = [
    0.006_000_600_060_006_002,
    0.014_001_400_140_014_003,
    0.862_086_208_620_862_2,
    0.080_008_000_800_080_01,
    0.030_003_000_300_030_006,
    0.002_400_240_024_002_400_4,
    0.004_000_400_040_004_001,
    0.001_500_150_015_001_500_5,
];

fn uc_curve() -> Vec<TurboUcPoint> {
    [
        (0.996_475_135_962_444_9, 0.984_090_909_090_909),
        (0.759_083_511_321_354_1, 0.796_590_909_090_909),
        (0.984_295_619_176_559, 0.993_181_818_181_818_3),
        (0.882_779_980_339_782_1, 0.936_363_636_363_636_4),
        (0.955_246_026_988_092_2, 0.994_318_181_818_181_8),
        (1.0, 1.0),
    ]
    .into_iter()
    .map(|(uc, efficiency)| TurboUcPoint { uc, efficiency })
    .collect()
}

fn qn_efficiency_curve() -> Vec<TurboQnPoint> {
    [
        (0.5, 0.88),
        (0.7, 0.91),
        (0.85, 0.95),
        (1.0, 1.0),
        (1.2, 0.97),
        (1.4, 0.85),
        (1.6, 0.6),
    ]
    .into_iter()
    .map(|(qn, factor)| TurboQnPoint { qn, factor })
    .collect()
}

fn qn_head_curve() -> Vec<TurboQnPoint> {
    [
        (0.5, 1.1),
        (0.8, 1.05),
        (1.0, 1.0),
        (1.2, 0.9),
        (1.4, 0.7),
        (1.6, 0.4),
    ]
    .into_iter()
    .map(|(qn, factor)| TurboQnPoint { qn, factor })
    .collect()
}

fn with_bearing(bearing_loss: Option<f64>) -> TurboExpanderResult {
    let components: Vec<String> = COMPONENTS.iter().map(|name| (*name).to_string()).collect();
    turbo_expander_compressor(
        &components,
        6721.52342728258,
        &Z,
        pascals(6_095_000.0),
        kelvins(250.15),
        &components,
        6241.700991745513,
        &Z,
        pascals(4_200_000.0),
        kelvins(290.15),
        pascals(4_200_000.0),
        meters(0.424),
        6850.0,
        0.7,
        0.033_28,
        None,
        None,
        Some(square_meters(0.016_37)),
        None,
        bearing_loss.map(watts),
        Some(0.88),
        Some(0.81),
        Some(meters(20.47 * 1000.0 / 9.81)),
        None,
        Some(&uc_curve()),
        Some(&qn_efficiency_curve()),
        Some(&qn_head_curve()),
        None,
        None,
    )
    .expect("the captured state balances")
}

fn close(actual: f64, expected: f64, tolerance: f64, what: &str) {
    let error = (actual - expected).abs() / expected.abs().max(f64::MIN_POSITIVE);
    assert!(
        error <= tolerance,
        "{what}: {actual} against {expected}, relative error {error:.3e} > {tolerance:.3e}"
    );
}

#[test]
fn the_reference_state_balances_at_the_captured_speed() {
    let result = with_bearing(None);
    close(result.speed, 6494.405432697534, 0.03, "speed");
    close(
        result.power_expander.value,
        2_884_181.338_380_164_5,
        0.03,
        "power_expander",
    );
    close(
        result.power_compressor.value,
        2_884_172.349_666_985,
        0.03,
        "power_compressor",
    );
    close(
        result.compressor_outlet_p.value,
        49.662_749_105_318_27e5,
        0.03,
        "compressor_outlet_p",
    );
    close(
        result.expander_outlet_t.value,
        232.189_731_590_532,
        0.03,
        "expander_outlet_t",
    );
    close(
        result.expander_isentropic_efficiency,
        0.837_139_534_945_472_1,
        0.03,
        "efficiency",
    );
    close(
        result.compressor_polytropic_head.value,
        1941.075,
        0.03,
        "polytropic_head",
    );
}

#[test]
fn the_bearing_loss_is_the_only_imbalance() {
    let none = with_bearing(Some(0.0));
    close(
        none.power_expander.value,
        none.power_compressor.value,
        1.0e-9,
        "the shaft balance without a bearing loss",
    );
    let lossy = with_bearing(None);
    let bearing = 10.0 * (lossy.speed / 6850.0) * (lossy.speed / 6850.0);
    close(
        lossy.power_expander.value,
        lossy.power_compressor.value + bearing,
        1.0e-9,
        "the shaft balance with the class's 10 W bearing loss",
    );
}
