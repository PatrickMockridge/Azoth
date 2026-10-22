//! The manual PHREEQC CO2-chloride dataset the Pitzer hydrate flash runs on.

use azoth_eos::pitzer_phase::{PitzerDatasetSelection, pitzer_phase, pitzer_phase_with_dataset};

/// The reference coefficients at a stated composition, from NeqSim's own phase after
/// `applyPhreeqcCo2ChlorideParameters(zero zeta)` and a `TPflash` at 278 K and 30 bara:
/// `ln phi` for CO2, water, Na+ and Cl- are 3.6416399262585983, -8.158005468187374,
/// 24.004688973360274 and 23.778694992509898, and the water activity is 0.9454110770199476.
#[test]
fn the_manual_dataset_reproduces_the_reference_coefficients() {
    let names = ["CO2", "water", "Na+", "Cl-"];
    let x = [
        0.02143766156847338,
        0.9478045231457128,
        0.015378907642906867,
        0.015378907642906867,
    ];
    let manual = pitzer_phase_with_dataset(
        &names,
        278.0,
        3.0e6,
        &x,
        PitzerDatasetSelection::PhreeqcCo2Chloride,
    )
    .expect("the manual dataset covers this brine");
    for (index, expected) in [
        3.6416399262585983,
        -8.158005468187374,
        24.004688973360274,
        23.778694992509898,
    ]
    .into_iter()
    .enumerate()
    {
        assert!(
            (manual.ln_phi[index] - expected).abs() < 1.0e-12,
            "ln phi[{index}] is {} against the reference {expected}",
            manual.ln_phi[index]
        );
    }

    // The automatic rule falls back to the legacy CSV on this topology, which is why the
    // manual dataset is a selection rather than the model's default.
    let automatic = pitzer_phase(&names, 278.0, 3.0e6, &x).expect("the automatic rule answers");
    assert!(manual.water_activity < automatic.water_activity);
}
