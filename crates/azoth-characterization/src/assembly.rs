//! A `Mixture` rebuilt from a characterization table.
//!
//! This is the step the class does with `addTBPfraction` + `profile.applyTo`, and it is a
//! library helper rather than a calculation: its answer is a fluid, and a fluid is not a result.
//! A caller holds [`FluidComponent`] rows - from `characterization.combine_reservoir_fluids`,
//! from a table of its own, or from a re-cut - and gets the mixture and its ideal-gas model,
//! which are the two things a flash and an enthalpy both take.
//!
//! **Nothing is mutated.** The mixture is constructed from its components and its interaction
//! matrix in one call, because `Mixture` is frozen once built and NeqSim's own route - create a
//! system, remove every component, add a TBP fraction, then overwrite sixteen of its fields -
//! is a sequence with a half-built state in it.
//!
//! # What a base row and a pseudo row resolve to
//!
//! A base row is a name: its `Tc`, `Pc`, `omega`, the interaction pairs it has with its
//! neighbours and its heat-capacity polynomial are the databank's, and a row that carries
//! properties of its own has them ignored. A pseudo row is the opposite, and must carry every
//! number a component needs: `Tc`, `Pc`, `omega`, the molar mass and all four heat-capacity
//! coefficients.
//!
//! The Peneloux translation is not a row field, because it is a function of the Rackett
//! compressibility, the acentric factor and the critical constants - the same derivation the
//! databank path makes, run over the row's `racket_z`. Both cubics' forms are carried, and the
//! mixture's cubic chooses.
//!
//! `ideal_gas_enthalpy_of_formation` is read by the blend and has no home here: neither
//! `Component` nor `IdealGasModel` carries it, so it is dropped at assembly rather than stored
//! somewhere nothing reads.

use azoth_core::units::{kelvins, pascals};
use azoth_core::{AzothError, Result};
use azoth_eos::databank;
use azoth_eos::{
    Alpha, Component, Cubic, IdealGasModel, Mixture, pr_peneloux_shift, srk_kappa,
    srk_peneloux_shift,
};

use crate::combiner::CP_D_GETTER_OFFSET;
use crate::component_row::FluidComponent;

/// A mixture and its ideal-gas model, from a characterization table.
///
/// # Errors
/// * [`azoth_core::AzothError::InvalidInput`] if a row is a pseudo-component without the
///   critical constants, molar mass or heat-capacity coefficients a component needs, or if
///   `cubic` is one this helper does not build.
/// * [`azoth_core::AzothError::PropertyUnavailable`] if a base row's name is not in the
///   databank, or its substance has no heat-capacity polynomial.
/// * Propagates [`Component::new`]'s range checks and the two Peneloux shifts' checks.
pub fn mixture_from_characterization(
    components: &[FluidComponent],
    cubic: Cubic,
) -> Result<(Mixture, IdealGasModel)> {
    if !matches!(cubic, Cubic::Pr | Cubic::Srk) {
        return Err(AzothError::invalid_input(
            "cubic",
            "a TBP fluid is characterised for Peng-Robinson or Soave-Redlich-Kwong; the other \
             cubics have no Peneloux translation an assembled cut can carry",
        ));
    }

    let mut built = Vec::with_capacity(components.len());
    let mut names = Vec::with_capacity(components.len());
    let mut base_names: Vec<Option<String>> = Vec::with_capacity(components.len());
    // One row per component: its five heat-capacity coefficients, in K**0..K**4.
    let mut cp: Vec<[f64; 5]> = Vec::with_capacity(components.len());
    for (index, row) in components.iter().enumerate() {
        names.push(row.name.clone());
        if row.pseudo_component {
            built.push(pseudo_component(index, row, cubic)?);
            cp.push(pseudo_cp(index, row)?);
            base_names.push(None);
        } else {
            let entry = databank::entry(&row.name, None)?;
            if entry.class == databank::ION {
                return Err(AzothError::invalid_input(
                    "components",
                    format!(
                        "{} is an ion, and a cubic has no notion of one; an ion belongs to the \
                         electrolyte models",
                        row.name
                    ),
                ));
            }
            let coefficients = entry.cp.ok_or_else(|| {
                AzothError::property_unavailable(
                    row.name.clone(),
                    "heat-capacity coefficients".to_string(),
                    "a mixture's ideal-gas model needs the polynomial, and this substance has \
                     none - which is what a keycard-added substance carries"
                        .to_string(),
                )
            })?;
            let mut component = entry.component()?;
            if matches!(cubic, Cubic::Srk) {
                // Soave's coefficient from the acentric factor, unless the row already
                // carries a fitted one - every shipped databank row leaves it empty.
                let fitted = component
                    .alpha_params
                    .first()
                    .copied()
                    .unwrap_or(srk_kappa(entry.omega)?.kappa);
                component = component.with_alpha_params(vec![fitted]);
            }
            built.push(component);
            cp.push(coefficients);
            base_names.push(Some(entry.name));
        }
    }

    let n = built.len();
    let mut kij = vec![0.0; n * n];
    for first in 0..n {
        for second in (first + 1)..n {
            let (Some(left), Some(right)) = (&base_names[first], &base_names[second]) else {
                continue;
            };
            let value = databank::kij(left, right, cubic, None);
            kij[first * n + second] = value;
            kij[second * n + first] = value;
        }
    }

    let coefficient = |position: usize| cp.iter().map(|row| row[position]).collect();
    let ideal_gas = IdealGasModel {
        cp_a: coefficient(0),
        cp_b: coefficient(1),
        cp_c: coefficient(2),
        cp_d: coefficient(3),
        cp_e: coefficient(4),
    };

    let mut mixture = Mixture::new(built, kij)?
        .with_names(names)
        .with_cubic(cubic);
    if matches!(cubic, Cubic::Srk) {
        // Every component carries a coefficient: the base rows take Soave's from their acentric
        // factor, which is the same number `Alpha::Srk` would compute.
        mixture = mixture.with_alpha(Alpha::SrkFitted);
    }
    Ok((mixture, ideal_gas))
}

/// One pseudo-component as a `Component`, its fitted alpha coefficient and Peneloux shifts.
fn pseudo_component(index: usize, row: &FluidComponent, cubic: Cubic) -> Result<Component> {
    let missing = |field: &str| {
        AzothError::invalid_input(
            format!("components[{index}].{field}"),
            format!("a pseudo-component needs its own {field}, and this row states none"),
        )
    };
    let tc = kelvins(
        row.critical_temperature
            .ok_or_else(|| missing("critical_temperature"))?,
    );
    let pc = pascals(
        row.critical_pressure
            .ok_or_else(|| missing("critical_pressure"))?,
    );
    let omega = row
        .acentric_factor
        .ok_or_else(|| missing("acentric_factor"))?;
    let molar_mass = row.molar_mass.ok_or_else(|| missing("molar_mass"))?;

    let pr_shift = pr_peneloux_shift(omega, tc, pc, row.racket_z)?;
    let srk_shift = srk_peneloux_shift(omega, tc, pc, row.racket_z)?;
    let mut component = Component::new(tc, pc, omega)?
        .with_molar_mass(Some(molar_mass))
        .with_class("hc".to_owned())
        .with_parachor(row.parachor.unwrap_or(0.0))
        .with_wax_data(
            false,
            row.heat_of_fusion.unwrap_or(0.0),
            row.triple_point_temperature.unwrap_or(0.0),
        )
        .with_volume_shift(pr_shift.c.value)
        .with_volume_shift_srk(srk_shift.c.value);
    if matches!(cubic, Cubic::Srk) {
        // NeqSim's `addTBPfraction` installs the fitted `m`; a row that states none is a cut
        // whose alpha is Soave's own, which is the correlation's default rather than a zero.
        let m = match row.attractive_m {
            Some(m) => m,
            None => srk_kappa(omega)?.kappa,
        };
        component = component.with_alpha_params(vec![m]);
    }
    Ok(component)
}

/// A pseudo-component's heat-capacity polynomial, all four coefficients required.
fn pseudo_cp(index: usize, row: &FluidComponent) -> Result<[f64; 5]> {
    let missing = |field: &str| {
        AzothError::invalid_input(
            format!("components[{index}].{field}"),
            format!(
                "a pseudo-component needs its own {field}; a mixture's enthalpy needs the \
                 polynomial rather than a row of zeros"
            ),
        )
    };
    Ok([
        row.cp_a.ok_or_else(|| missing("cp_a"))?,
        row.cp_b.ok_or_else(|| missing("cp_b"))?,
        row.cp_c.ok_or_else(|| missing("cp_c"))?,
        // The row states what `getCpD` states; `setCpD` stores the number without the
        // getter's offset, exactly as `PseudoComponentCombiner.applyTo` does.
        row.cp_d.ok_or_else(|| missing("cp_d"))? - CP_D_GETTER_OFFSET,
        // NeqSim's TBP path never sets `CpE`, so the quartic and higher term is absent.
        0.0,
    ])
}
