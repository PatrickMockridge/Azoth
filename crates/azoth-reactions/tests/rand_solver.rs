//! The modified-RAND solve, against the reactive flash's own capture.
//!
//! The oracle is `validation/neqsim/ReactiveFlashProbe.java`, whose capture is
//! `captures/reactive_flash_probe.tsv`. The water-gas shift at 600 K and 1 bar is the state
//! its authors chose: `CO + H2O = CO2 + H2` has one independent reaction over three elements,
//! the equilibrium constant is about 14, so the products are strongly favoured and the
//! direction of the answer is not in doubt.
//!
//! **The answer is `getEquilibriumMoles()[0]`**, the moles the RAND solve leaves in the first
//! phase. The driver ran its outer loop on that state rather than its single-phase branch -
//! `SystemSrkEos` constructs with two phase objects, so `np = 2` and `maxPhases = 2` when the
//! flash is handed the system, and both phases end up holding the same answer. The element
//! inventory `b` is the feed's own `A n`, which is what `FormulaMatrix.computeElementVector`
//! builds when the flash has no frozen target to hand.
//!
//! The fugacity coefficients come from `azoth-eos`' SRK, which is the cubic the probe's
//! `SystemSrkEos` carries; the standard potentials come from the component databank's three
//! formation columns and its heat-capacity polynomial, which is what `computeG0` reads.

use azoth_core::units::{kelvins, pascals};
use azoth_eos::databank::mixture_of;
use azoth_eos::{Cubic, RootSide};
use azoth_reactions::databank::formation_properties;
use azoth_reactions::formula_matrix::FormulaMatrix;
use azoth_reactions::rand_solver::{ThermoData, solve_single_phase, standard_potentials};

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

    let solution =
        solve_single_phase(&matrix.matrix, &g0, &b, &FEED, &mut ln_phi).expect("the solve runs");

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

/// **The driver's own numbers, not only the composition.** `ReactiveMultiphaseTPflash` reports
/// three things for the path it takes on this fluid - the equilibrium moles, their total, and
/// the Gibbs measure `computeGibbsEnergy` - and all three are in the capture. The total is the
/// one that is easy to get wrong by assuming: it is `getTotalMoles` of the solve, and a
/// reaction that splits one species into two moves it away from the feed's.
#[test]
fn the_drivers_single_phase_answer_is_the_captured_one() {
    use azoth_reactions::reactive_flash::single_phase_equilibrium;

    let matrix =
        FormulaMatrix::build(&NAMES.map(String::from)).expect("the components are in the databank");
    let b: Vec<f64> = matrix
        .matrix
        .iter()
        .map(|row| row.iter().zip(FEED).map(|(a, n)| a * n).sum())
        .collect();
    let g0 = standard_potentials(&thermo_data(), TEMPERATURE, PRESSURE_BARA);

    let (mixture, _) = mixture_of(&NAMES, Cubic::Srk, None).expect("the databank carries them");
    let reduced = mixture
        .reduced_parameters(kelvins(TEMPERATURE), pascals(PRESSURE_BARA * 1.0e5))
        .expect("a state");
    let mut ln_phi = |x: &[f64]| -> azoth_core::Result<Vec<f64>> {
        Ok(mixture
            .phase_state(&reduced, x, RootSide::Vapour)?
            .ln_phi
            .clone())
    };

    let outcome =
        single_phase_equilibrium(&matrix.matrix, &g0, &b, &FEED, &mut ln_phi).expect("it solves");
    assert!(outcome.converged, "the captured state converges");

    // `equilibrium_total_moles` from the capture.
    let captured_total = 0.999_999_128_896_382_f64;
    assert!(
        (outcome.total_moles - captured_total).abs() / captured_total < 1.0e-5,
        "the total is {} against the capture's {captured_total}",
        outcome.total_moles
    );

    // `final_gibbs_energy` from the capture: `sum_i x_i (ln x_i + ln phi_i)`, dimensionless -
    // **and exactly twice this port's**, because the driver weighs each phase by
    // `phase.getBeta()` and the capture's two phases - the pair `SystemSrkEos` constructs, both
    // at `beta = 1.0` - sum to two. The port models the phase the solve leaves, so it reports
    // the thermodynamic value; `reactive_flash.rs` reproduces the capture's number from those
    // two phases and their weights. The ratio is asserted here rather than the number, which is
    // what makes the doubling a measurement instead of a rounding difference.
    let captured_gibbs = -2.259_535_542_715_054_7_f64;
    let ratio = captured_gibbs / outcome.gibbs_energy;
    assert!(
        (ratio - 2.0).abs() < 1.0e-4,
        "the capture's Gibbs measure is {captured_gibbs} against this port's {}, a ratio of \
         {ratio} rather than the two the duplicate phase accounts for",
        outcome.gibbs_energy
    );
}
