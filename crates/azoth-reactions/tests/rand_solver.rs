//! The modified-RAND solve, against the reactive flash's own capture.
//!
//! The oracle is `validation/neqsim/ReactiveFlashProbe.java`, whose capture is
//! `captures/reactive_flash_probe.tsv`. The water-gas shift at 600 K and 1 bar is the state
//! its authors chose: `CO + H2O = CO2 + H2` has one independent reaction over three elements,
//! the equilibrium constant is about 14, so the products are strongly favoured and the
//! direction of the answer is not in doubt.
//!
//! **The answer is `getEquilibriumMoles()[0]`**, the moles the RAND solve leaves in the single
//! phase - not the phase composition, which the flash's own trial-phase removal leaves in two
//! near-identical copies. The element inventory `b` is the feed's own `A n`, which is what
//! `FormulaMatrix.computeElementVector` builds when the flash has no frozen target to hand.
//!
//! The fugacity coefficients come from `azoth-eos`' SRK, which is the cubic the probe's
//! `SystemSrkEos` carries; the standard potentials come from the component databank's three
//! formation columns and its heat-capacity polynomial, which is what `computeG0` reads.

use azoth_core::units::{kelvins, pascals};
use azoth_eos::databank::mixture_of;
use azoth_eos::{Cubic, RootSide};
use azoth_reactions::databank::formation_properties;
use azoth_reactions::formula_matrix::FormulaMatrix;
use azoth_reactions::rand_solver::{ThermoData, solve, standard_potentials};

/// The water-gas shift's fluid, in the order the capture prints it.
const NAMES: [&str; 4] = ["CO", "water", "CO2", "hydrogen"];
const FEED: [f64; 4] = [0.25, 0.25, 0.25, 0.25];
const TEMPERATURE: f64 = 600.0;
const PRESSURE_BARA: f64 = 1.0;

/// From `captures/reactive_flash_probe.tsv`, `fluid=wgs-600K`, `equilibrium_moles[0]`.
const CAPTURED_MOLES: [f64; 4] = [
    0.079_140_502_548_023_4,
    0.079_140_502_560_693_7,
    0.420_859_061_906_977_14,
    0.420_859_061_880_687_67,
];

fn thermo_data() -> Vec<ThermoData> {
    let (_, ideal) = mixture_of(&NAMES, Cubic::Srk, None).expect("the databank carries them");
    NAMES
        .iter()
        .enumerate()
        .map(|(index, name)| {
            let formation = formation_properties(name)
                .expect("the component table parses")
                .expect("the component databank carries it");
            ThermoData {
                enthalpy_of_formation: formation.enthalpy_of_formation,
                absolute_entropy: formation.absolute_entropy,
                gibbs_energy_of_formation: formation.gibbs_energy_of_formation,
                cp: [
                    ideal.cp_a[index],
                    ideal.cp_b[index],
                    ideal.cp_c[index],
                    ideal.cp_d[index],
                    ideal.cp_e[index],
                ],
            }
        })
        .collect()
}

#[test]
fn the_water_gas_shift_reproduces_the_captured_equilibrium() {
    let matrix =
        FormulaMatrix::build(&NAMES.map(String::from)).expect("the components are in the databank");
    assert_eq!(matrix.element_names, ["C", "O", "H"]);
    assert_eq!(matrix.independent_reactions(), 1, "one reaction runs");

    // The element inventory, `A n` of the feed.
    let b: Vec<f64> = matrix
        .matrix
        .iter()
        .map(|row| row.iter().zip(FEED).map(|(a, n)| a * n).sum())
        .collect();

    // The pressure is in **bara**, which is the unit NeqSim's `PP` is in against a `P_REF` of
    // one - so at one bar the `ln(P/P_ref)` term is zero and the potentials are the standard
    // ones.
    let g0 = standard_potentials(&thermo_data(), TEMPERATURE, PRESSURE_BARA);

    let (mixture, _) = mixture_of(&NAMES, Cubic::Srk, None).expect("the databank carries them");
    let reduced = mixture
        .reduced_parameters(kelvins(TEMPERATURE), pascals(PRESSURE_BARA * 1.0e5))
        .expect("a state");
    let mut ln_phi = |x: &[f64]| -> azoth_core::Result<Vec<f64>> {
        let state = mixture.phase_state(&reduced, x, RootSide::Vapour)?;
        Ok(state.ln_phi.clone())
    };

    let solution = solve(&matrix.matrix, &g0, &b, &FEED, &mut ln_phi).expect("the solve runs");

    assert!(solution.converged, "the captured state converges");
    // **The band is the measurement.** The two codes stop on the same `1e-9` residual but
    // reach it along different paths - the linear solve here is this crate's own elimination
    // and the fugacity coefficients are SRK's own rounding - and the worst component lands
    // `3.6e-7` from the capture. A tighter bound would be asserting bit equality between two
    // iterative codes, which is not a stricter test but an impossible one.
    for (index, (got, want)) in solution.moles.iter().zip(CAPTURED_MOLES).enumerate() {
        let relative = (got - want).abs() / want.abs();
        assert!(
            relative < 1.0e-5,
            "component {index} ({}): {got} against the capture's {want}, a relative {relative:.3e}",
            NAMES[index]
        );
    }

    // The carbon balance is exact by construction, which is what the class's own test checks.
    let carbon = solution.moles[0] + solution.moles[2];
    assert!((carbon - 0.5).abs() < 1.0e-6, "carbon came out at {carbon}");
    // And the products are favoured, so the answer is on the right side of the feed.
    assert!(
        solution.moles[2] > 0.25,
        "CO2 rose to {}",
        solution.moles[2]
    );
    assert!(solution.moles[0] < 0.25, "CO fell to {}", solution.moles[0]);
}
