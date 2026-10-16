//! `eos.calc_ionic_composition` - molalities and activity coefficients of one phase.
//!
//! Spec: `specs/models/eos/calc_ionic_composition.toml`. A *direct* model: the class solves
//! nothing, it reads a phase. What makes it more than arithmetic is the third column.
//!
//! Ported from `CalcIonicComposition`, 90 lines, whose `run` formats three strings per
//! component - and from `Phase`'s two reference-state constructions, which the third column
//! calls and which are the whole of the work here.
//!
//! # The reference states, and why the coefficient is not a ratio to water
//!
//! `Phase.getActivityCoefficient(k, p)` is `exp(ln phi_k - ln phi_k^reference)`, where the
//! reference is chosen per component by the `REFERENCESTATETYPE` column:
//!
//! * `solvent` takes [`pure_log_fugacity`] - a phase of that one component, at the same `T`
//!   and `P`, so the coefficient is the symmetric (Raoult) one;
//! * anything else takes [`infinite_dilute_log_fugacity`] - a phase of `1e-10` mol of the
//!   component beside `10` mol of water, so the coefficient is the unsymmetric one.
//!
//! So a brine's ionic coefficient is against the ion at **infinite dilution in water**, not
//! against water in the same phase. The two are different numbers and the difference is the
//! model: at the 40 bara state `Na+` reads `0.687` unsymmetrically, and a ratio to water's
//! own coefficient would read about one.
//!
//! Both references are *phases of the same fluid at a different composition*, which is why
//! this model is a model and not a column of arithmetic.

use azoth_core::Result;
use azoth_core::units::{Pressure, ThermodynamicTemperature};

use crate::databank::Entry;
use crate::mixture::RootSide;
use crate::model_gen;
use crate::results::CalcIonicCompositionResult;

/// The solute's amount in the reference binary, the class's own.
const DILUTE_MOLES: f64 = 1.0e-10;

/// The solvent's amount beside it, the class's own.
const SOLVENT_MOLES: f64 = 10.0;

/// The three columns of `CalcIonicComposition`'s table, per component.
///
/// # Errors
/// * [`azoth_core::AzothError::InvalidInput`] if `x` is not one entry per component or is not
///   a composition, a name is not in the databank, or no component is named `water`.
/// * [`azoth_core::AzothError::PropertyUnavailable`] if a component carries no heat-capacity
///   coefficients, or the databank states no molar mass for one.
/// * [`azoth_core::AzothError::OutOfRange`] if `T` or `P` is not positive, or no volume root
///   exists at the phase or at one of its reference states.
#[allow(non_snake_case)] // `T` and `P` are the symbols in the chemistry
pub fn calc_ionic_composition(
    components: &[String],
    T: ThermodynamicTemperature,
    P: Pressure,
    x: &[f64],
) -> Result<CalcIonicCompositionResult> {
    let spec = &model_gen::CALC_IONIC_COMPOSITION_SPEC;
    let mut warnings = Vec::new();
    azoth_core::apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "T" => Some(T.value),
            "P" => Some(P.value),
            _ => None,
        },
        &mut warnings,
    )?;

    let n = components.len();
    if x.len() != n {
        return Err(azoth_core::AzothError::invalid_input(
            "x",
            format!("a mixture of {n} components has {} mole fractions", x.len()),
        ));
    }
    let sum: f64 = x.iter().sum();
    if (sum - 1.0).abs() > 1.0e-09 {
        return Err(azoth_core::AzothError::invalid_input(
            "x",
            format!(
                "the mole fractions sum to {sum}, not to one. Renormalising them here would \
                 make a composition error invisible in every number downstream, so it is \
                 refused instead"
            ),
        ));
    }

    // **The solvent is found by name**, as the class finds it: it asks the phase for its
    // `water` component rather than for one whose reference state is `solvent`.
    let water = components
        .iter()
        .position(|name| name.eq_ignore_ascii_case("water"))
        .ok_or_else(|| {
            azoth_core::AzothError::invalid_input(
                "components",
                "the two concentration columns are per kilogram of water, and no component is \
                 named `water`",
            )
        })?;

    let names: Vec<&str> = components.iter().map(String::as_str).collect();
    let entries: Vec<Entry> = names
        .iter()
        .map(|name| crate::databank::entry(name, None))
        .collect::<Result<Vec<_>>>()?;
    let water_mass = entries[water].molar_mass.ok_or_else(|| {
        azoth_core::AzothError::property_unavailable(
            names[water],
            "a molar mass",
            "the concentration columns are per kilogram of it",
        )
    })?;

    let (mixture, _) = crate::furst_electrolyte::electrolyte_cpa_statoil_mixture_of(&names, None)?;
    let reduced = mixture.reduced_parameters(T, P)?;
    warnings.extend(reduced.warnings.iter().cloned());
    let state = mixture.phase_state(&reduced, x, RootSide::Liquid)?;

    let mut molality = Vec::with_capacity(n);
    let mut mass_ratio = Vec::with_capacity(n);
    let mut activity_coefficient = Vec::with_capacity(n);
    for i in 0..n {
        // The class's `mmol/kgSolvent` column is this times 1000 and its `mg/kgSolvent` is
        // this times the molar mass times `1e6`; the vocabulary states the quantity itself,
        // so both are carried in their own unit rather than as the class's multiple.
        let per_kg = x[i] / x[water] / water_mass;
        molality.push(per_kg);
        let mass = entries[i].molar_mass.ok_or_else(|| {
            azoth_core::AzothError::property_unavailable(
                names[i],
                "a molar mass",
                "the mass column is the molar one times it",
            )
        })?;
        mass_ratio.push(per_kg * mass);

        let reference = if entries[i].reference_state == crate::databank::SOLVENT {
            pure_log_fugacity(names[i], T, P)?
        } else {
            infinite_dilute_log_fugacity(names[i], names[water], T, P)?
        };
        activity_coefficient.push((state.ln_phi[i] - reference).exp());
    }

    Ok(CalcIonicCompositionResult {
        molality,
        mass_ratio,
        activity_coefficient,
        warnings,
    })
}

/// `ln phi` of one component in a phase of itself, at the same `T` and `P`.
///
/// The class's `getLogPureComponentFugacity`, which builds a one-component phase and reads
/// its fugacity coefficient. The root is the liquid one because the phase being asked about
/// is an aqueous one - `initRefPhases` copies `this.getType()`.
#[allow(non_snake_case)] // `T` and `P` are the symbols in the chemistry
fn pure_log_fugacity(name: &str, T: ThermodynamicTemperature, P: Pressure) -> Result<f64> {
    let (mixture, _) = crate::furst_electrolyte::electrolyte_cpa_statoil_mixture_of(&[name], None)?;
    let reduced = mixture.reduced_parameters(T, P)?;
    let state = mixture.phase_state(&reduced, &[1.0], RootSide::Liquid)?;
    Ok(state.ln_phi[0])
}

/// `ln phi` of one component at infinite dilution in the solvent, at the same `T` and `P`.
///
/// The class's `getLogInfiniteDiluteFugacity(k, p)`, whose binary is `1e-10` mol of the
/// component beside `10` mol of the solvent. The composition is that ratio and not a limit:
/// the phase is solved there.
#[allow(non_snake_case)] // `T` and `P` are the symbols in the chemistry
fn infinite_dilute_log_fugacity(
    name: &str,
    solvent: &str,
    T: ThermodynamicTemperature,
    P: Pressure,
) -> Result<f64> {
    let total = DILUTE_MOLES + SOLVENT_MOLES;
    let x = [DILUTE_MOLES / total, SOLVENT_MOLES / total];
    let (mixture, _) =
        crate::furst_electrolyte::electrolyte_cpa_statoil_mixture_of(&[name, solvent], None)?;
    let reduced = mixture.reduced_parameters(T, P)?;
    let state = mixture.phase_state(&reduced, &x, RootSide::Liquid)?;
    Ok(state.ln_phi[0])
}
