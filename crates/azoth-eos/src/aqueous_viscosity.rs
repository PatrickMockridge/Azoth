//! `eos.aqueous_viscosity` - the liquid viscosity NeqSim gives an aqueous phase.
//!
//! Spec: `specs/models/eos/aqueous_viscosity.toml`. A port of the liquid `Viscosity`
//! class - the one `WaterPhysicalProperties` selects for `PhaseType.AQUEOUS` - through
//! `calcViscosity` and `calcPureComponentViscosity`.
//!
//! ```text
//! mu_i  = f(LIQVISCMODEL_i, LIQVISC_i, T) * (correction_i + 1)/2
//! mu    = exp(sum_i w_i ln mu_i)                    <- Grunberg-Nissan, G_ij = 0
//! ```
//!
//! **The mixing rule's interaction term is zero, and that is a measurement rather than
//! an omission.** `PhysicalPropertyMixingRule.initMixingRules` fills `Gij` with an inner
//! loop that starts at `k = l` and breaks on `k == l`, so its body never runs and every
//! phase reports a zero matrix - the probe prints `G[0][1] = 0.0` for a water/methanol
//! mixture whose `Gij` column in `INTER.csv` is not empty. Reproducing the class means
//! reproducing the zero.

use azoth_core::units::{Pressure, ThermodynamicTemperature, pascal_seconds};
use azoth_core::{AzothError, Result, apply_checks};

use crate::databank::ION;
use crate::mixture::Mixture;
use crate::model_gen;
use crate::results::AqueousViscosityResult;

/// The viscosity a component takes when its row names no model at all: NeqSim's `else` branch.
///
/// The above-critical sentinel is `eos.liquid_viscosity_pure`'s now, since that is where the
/// ladder lives.
const NO_MODEL_CP: f64 = 0.7;
/// `mPa*s` to `Pa*s`, applied once at the end: every correlation here is written in cP.
const CP_TO_PA_S: f64 = 1.0e-3;

/// The liquid viscosity of an aqueous phase, in Pa·s.
///
/// # Errors
/// * [`azoth_core::AzothError::InvalidInput`] if `z` is the wrong length, or **if any
///   component is an ion** - see the spec's own assumption: NeqSim's `LIQVISC` rows for
///   `na+` and `cl-` are methanol's four numbers, and a brine computed from them comes out
///   less viscous than pure water.
/// * [`azoth_core::AzothError::PropertyUnavailable`] if a component carries no molar mass,
///   which the mass-fraction weighting needs.
pub fn aqueous_viscosity(
    mixture: &Mixture,
    t: ThermodynamicTemperature,
    p: Pressure,
    z: &[f64],
) -> Result<AqueousViscosityResult> {
    let spec = &model_gen::AQUEOUS_VISCOSITY_SPEC;
    let mut warnings = Vec::new();

    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "T" => Some(t.value),
            "P" => Some(p.value),
            _ => None,
        },
        &mut warnings,
    )?;

    let components = mixture.components();
    let n = components.len();
    if z.len() != n {
        return Err(AzothError::invalid_input(
            "z",
            format!(
                "a mixture of {n} components needs {n} mole fractions, got {}",
                z.len()
            ),
        ));
    }

    // An ion's liquid-viscosity row is another substance's, so a brine would rest on it.
    // A `Component` carries critical constants and no name, so the ion is named from the
    // mixture's own name list - and by index where a mixture built from constants has none.
    let names = mixture.names();
    let ions: Vec<String> = components
        .iter()
        .enumerate()
        .filter(|(_, c)| c.class == ION)
        .map(|(i, _)| match names {
            Some(names) if i < names.len() => names[i].clone(),
            _ => format!("the component at index {i}"),
        })
        .collect();
    if !ions.is_empty() {
        return Err(AzothError::invalid_input(
            "components",
            format!(
                "{} is an ion, and NeqSim's `LIQVISC` row for an ion carries the same four \
                 numbers as methanol's - so a brine computed from it comes out less viscous \
                 than pure water, which is the wrong direction. The class that would close \
                 this is a salt-aware viscosity, not this correlation",
                ions.join(", ")
            ),
        ));
    }

    // Grunberg-Nissan weights by **mass** fraction, so the molar masses are load-bearing.
    let mut masses = Vec::with_capacity(n);
    for (component, zi) in components.iter().zip(z) {
        let m = component.molar_mass.ok_or_else(|| {
            AzothError::property_unavailable(
                "component",
                "molar mass",
                "a card-added component needs its own molar mass: the mixing rule weights by \
                 mass fraction",
            )
        })?;
        masses.push(zi * m);
    }
    let total_mass: f64 = masses.iter().sum();
    if total_mass <= 0.0 {
        return Err(AzothError::invalid_input(
            "z",
            "a phase whose mass is zero has no mass fractions to weight by",
        ));
    }

    let mut weighted = 0.0;
    for (index, component) in components.iter().enumerate() {
        let weight = masses[index] / total_mass;
        // A `Some` weight of zero contributes nothing, and `ln(0)` would poison it - the
        // class's own loop has the same structure but never sees an absent component
        // because a phase's mole fractions are positive.
        if weight <= 0.0 {
            continue;
        }
        weighted += weight * pure_viscosity(component, t.value, p.value).ln();
    }

    Ok(AqueousViscosityResult {
        viscosity: pascal_seconds(weighted.exp() * CP_TO_PA_S),
        warnings,
    })
}

/// One component's pure-liquid viscosity in cP, with the pressure correction applied.
///
/// **The ladder is `eos.liquid_viscosity_pure`'s, called and not copied.** This id is that
/// model's mixture rule; a second copy of the branch would be a second answer to the same
/// question, and the *liquid* variant is the one this phase's own viscosity class carries.
fn pure_viscosity(component: &crate::mixture::Component, temperature: f64, pressure: f64) -> f64 {
    let [l1, l2, l3, l4] = component.liqvisc;
    crate::liquid_viscosity_pure::liquid_viscosity_pure(
        crate::liquid_viscosity_pure::LiquidViscosityLadder::Liquid,
        component.liqvisc_model,
        l1,
        l2,
        l3,
        l4,
        component.tc,
        component.pc,
        component.omega,
        azoth_core::units::kelvins(temperature),
        azoth_core::units::pascals(pressure),
    )
    .map_or(NO_MODEL_CP, |r| r.mu.value / CP_TO_PA_S)
}
