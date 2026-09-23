//! The mass-transfer rate matrix, from `Kinetics.calcReacMatrix`.
//!
//! Krishna and Standart's film model as NeqSim writes it: a pseudo-first-order coefficient per
//! component, built from the fluid's reactions, the two phases of the film, and each component's
//! effective diffusion coefficient.
//!
//! ```text
//! coefficient = sum over reactions of  k * prod_j c_j^(-nu_j) * prod_{k: nu_k nu_j > 0} c_k^(nu_k/nu_j)
//! ```
//!
//! where `c_j = x_j * rho / M_j` is a molar concentration, `k` is the reaction's rate factor and
//! the second product runs only over the siblings that appear on the **same side** of the
//! reaction as the component being asked about (`nu_k nu_j > 0`), excluding water and excluding
//! the component itself.
//!
//! # Three things the class does that a formula does not say
//!
//! **The mixture of phases is deliberate in the original and is kept.** The concentrations are
//! built from `interPhase`'s composition with `phase`'s *density* and *molar masses* - a pairing
//! that only makes sense because the two phases of a film share a mixture; the port takes both
//! and reads each field from the one the class reads it from.
//!
//! **The irreversibility test is `1/K` scaled by the reactant concentrations**, and it is
//! computed on every reaction whether or not the component appears in it. It sets a flag rather
//! than changing the coefficient.
//!
//! **`phiInfinite` is assigned, not accumulated**, so where more than one reaction has an inner
//! product the *last* one wins - and where none does, the class leaves the field holding
//! whatever the previous call left, which on a fresh object is zero. The port returns
//! [`RateMatrix::phi_infinite`] as an `Option` instead: `None` is "no reaction produced one",
//! which is the same fact without the stale value.
//!
//! # The interface and the bulk phase
//!
//! NeqSim's caller is `ReactiveKrishnaStandartFilmModel`, which lives in `fluidmechanics/` and is
//! beyond this port. So the caller supplies the pair - `phase` is the bulk and `inter_phase` the
//! film's interface, which in the capture are the same object because that is what the probe
//! passed - and nothing here knows which is which.

use azoth_core::{AzothError, CalcResult, Result, Warning, apply_checks};

/// The floor on `|1/K|` scaled by the reactant concentrations below which the class calls a
/// reaction irreversible, from `if (Math.abs(irr) < 1e-3)`.
pub const IRREVERSIBLE_THRESHOLD: f64 = 1.0e-3;

/// One reaction as the matrix reads it: its species, its stoichiometry, its rate factor and its
/// equilibrium constant at the state.
#[derive(Debug, Clone, PartialEq)]
pub struct KineticReaction {
    /// The species, in the reaction's own order.
    pub names: Vec<String>,
    /// The stoichiometric coefficients, signed and in the same order.
    pub stoc_coefs: Vec<f64>,
    /// `getRateFactor` at the interface's temperature, which the caller computes -
    /// [`crate::kinetic_rate_law::rate_factor`] is that law.
    pub rate_factor: f64,
    /// `getK` at the bulk phase's state, which the caller computes -
    /// `reactions.equilibrium_constant` is that id.
    pub equilibrium_constant: f64,
}

/// One phase as the matrix reads it: its composition, its components' molar masses and its
/// density - in kilograms per mole and kilograms per cubic metre, which is what makes the
/// concentration `x rho / M` a mol/m³.
#[derive(Debug, Clone, PartialEq)]
pub struct KineticPhase {
    /// The species, in the phase's own order.
    pub names: Vec<String>,
    /// Each species' mole fraction, in the same order.
    pub fractions: Vec<f64>,
    /// Each species' molar mass, in kg/mol, in the same order.
    pub molar_masses: Vec<f64>,
    /// The phase's mass density, in kg/m³.
    pub density: f64,
}

impl KineticPhase {
    /// The position of a species, which is how the class finds it: by name.
    ///
    /// # Errors
    /// [`AzothError::InvalidInput`] where the species is not in the phase.
    fn index(&self, name: &str) -> Result<usize> {
        self.names
            .iter()
            .position(|candidate| candidate == name)
            .ok_or_else(|| AzothError::InvalidInput {
                field: "names".to_string(),
                reason: format!("`{name}` is not in the phase"),
            })
    }
}

/// **The class builds three concentrations from two phases and they do not all pair the same
/// way.** Kept as three functions rather than one so the pairing is visible:
///
/// * the rate law's own `c_j` takes the *interface's* mole fraction, the **bulk's** density and
///   the bulk's molar mass;
/// * the reactive concentration takes the interface's fraction and density with the **bulk's**
///   molar mass;
/// * the sibling concentration is the bulk's throughout.
///
/// In the capture the two phases are the same object - the probe passes the aqueous phase for
/// both, which is what NeqSim's own film caller does when the film is not resolved - so the three
/// agree there and the pairing is not pinned by a number. It is pinned by the source, and it is
/// written out here so a caller who passes two different phases gets what the class would.
fn reaction_concentration(
    inter_phase: &KineticPhase,
    phase: &KineticPhase,
    name: &str,
) -> Result<f64> {
    let here = inter_phase.index(name)?;
    let mass = phase.index(name)?;
    Ok(inter_phase.fractions[here] * phase.density / phase.molar_masses[mass])
}

/// The interface's fraction and density over the bulk phase's molar mass.
fn reactive_concentration(
    inter_phase: &KineticPhase,
    phase: &KineticPhase,
    name: &str,
) -> Result<f64> {
    let here = inter_phase.index(name)?;
    let mass = phase.index(name)?;
    Ok(inter_phase.fractions[here] * inter_phase.density / phase.molar_masses[mass])
}

/// The bulk phase's own `x rho / M`.
fn bulk_concentration(phase: &KineticPhase, name: &str) -> Result<f64> {
    let index = phase.index(name)?;
    Ok(phase.fractions[index] * phase.density / phase.molar_masses[index])
}

/// What one component's rate matrix row answers with.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RateMatrix {
    /// `reacCoef`: the pseudo-first-order coefficient, summing every reaction's own.
    pub coefficient: f64,
    /// `getPhiInfinite` - **`None` where no reaction produced one**, which the class reports as
    /// whatever the previous call left.
    pub phi_infinite: Option<f64>,
    /// `isIrreversible`: whether any reaction's scaled `1/K` came in under
    /// [`IRREVERSIBLE_THRESHOLD`] while the component's row was being built.
    ///
    /// **The class's field is never reset**, so read off a shared object it is the union
    /// over every row taken from it so far. This is the per-call answer, which is what a
    /// fresh object reports - the capture holds both.
    pub irreversible: bool,
}

/// `calcReacMatrix`: one component's row of the mass-transfer rate matrix.
///
/// `diffusion` is the phase's effective diffusion coefficients, one per species in `phase`'s
/// order - the class reads them by index, and so does this.
///
/// # Errors
/// [`AzothError::InvalidInput`] where a reaction names a species the phase does not carry, where
/// a reaction's names and coefficients are different lengths, or where `diffusion` is not one
/// value per species.
pub fn rate_matrix(
    reactions: &[KineticReaction],
    phase: &KineticPhase,
    inter_phase: &KineticPhase,
    component: &str,
    diffusion: &[f64],
) -> Result<RateMatrix> {
    if diffusion.len() != phase.names.len() {
        return Err(AzothError::InvalidInput {
            field: "diffusion".to_string(),
            reason: format!(
                "{} coefficient(s) against {} species",
                diffusion.len(),
                phase.names.len()
            ),
        });
    }
    let component_index = phase.index(component)?;

    let mut coefficient = 0.0_f64;
    let mut phi_infinite: Option<f64> = None;
    let mut irreversible = false;

    for reaction in reactions {
        if reaction.names.len() != reaction.stoc_coefs.len() {
            return Err(AzothError::InvalidInput {
                field: "reactions".to_string(),
                reason: format!(
                    "`{}` has {} name(s) and {} coefficient(s)",
                    reaction.names.join(" "),
                    reaction.names.len(),
                    reaction.stoc_coefs.len()
                ),
            });
        }
        let mut ktemp = reaction.rate_factor;

        // `irr = (1/K) * prod_j [c_j]^(-nu_j)`, in which the concentration is the *interface's*
        // composition over the *bulk* phase's density and molar masses. Nothing reads it but the
        // threshold below.
        let mut irr = 1.0 / reaction.equilibrium_constant;
        for (j, name) in reaction.names.iter().enumerate() {
            let concentration = reaction_concentration(inter_phase, phase, name)?;
            irr *= concentration.powf(-reaction.stoc_coefs[j]);
        }
        if irr.abs() < IRREVERSIBLE_THRESHOLD {
            irreversible = true;
        }

        for (j, name) in reaction.names.iter().enumerate() {
            if name != component {
                continue;
            }
            for (k, sibling) in reaction.names.iter().enumerate() {
                // The same side of the reaction, not the component itself, and never water.
                if reaction.stoc_coefs[k] * reaction.stoc_coefs[j] <= 0.0
                    || k == j
                    || sibling == "water"
                {
                    continue;
                }
                let exponent = reaction.stoc_coefs[k] / reaction.stoc_coefs[j];
                let reactive = reactive_concentration(inter_phase, phase, component)?;
                let sibling_concentration = bulk_concentration(phase, sibling)?;
                ktemp *= sibling_concentration.powf(exponent);

                let sibling_index = phase.index(sibling)?;
                let forward = (diffusion[component_index] / diffusion[sibling_index]).sqrt();
                let backward = (diffusion[sibling_index] / diffusion[component_index]).sqrt();
                phi_infinite =
                    Some(forward + backward * sibling_concentration / (exponent * reactive));
            }
        }
        coefficient += ktemp;
    }

    Ok(RateMatrix {
        coefficient,
        phi_infinite,
        irreversible,
    })
}

/// Result of `reactions.kinetics`: one entry per component of `components`, in that order.
#[derive(Debug, Clone, PartialEq)]
pub struct KineticsResult {
    /// `reacCoef` per component.
    pub coefficient: Vec<f64>,
    /// `getPhiInfinite` per component, **zero where no reaction produced one** - the value
    /// the class's field is constructed with, which is what a fresh object reads back.
    pub phi_infinite: Vec<f64>,
    /// A mask: 1 where any reaction's scaled `1/K` came in under
    /// [`IRREVERSIBLE_THRESHOLD`] while this component's row was built.
    pub irreversible: Vec<f64>,
    /// Caveats.
    pub warnings: Vec<Warning>,
}

impl CalcResult for KineticsResult {
    const CALC_ID: &'static str = "reactions.kinetics";
    const FIELDS: &'static [&'static str] =
        &["coefficient", "phi_infinite", "irreversible", "warnings"];

    fn warnings(&self) -> &[Warning] {
        &self.warnings
    }
}

/// `reactions.kinetics`: the whole phase's mass-transfer rate matrix.
///
/// The reactions cross as **concatenated name lists** rather than a matrix over
/// `components`, because a reaction's own name order decides which sibling's phi it keeps
/// and a rectangular matrix cannot state that order. `reaction_lengths[i]` is how many
/// names reaction `i` contributes, and the two flattened vectors are aligned one for one
/// with `reaction_components`.
///
/// # Errors
/// [`AzothError::InvalidInput`] where the lengths disagree, where `diffusion` is not one
/// value per component, or where a reaction names a species the phase does not carry.
#[allow(clippy::too_many_arguments)] // the class's own parameter list, one input per fact
pub fn kinetics(
    components: &[&str],
    reaction_components: &[&str],
    reaction_lengths: &[f64],
    reaction_coefficients: &[f64],
    rate_factors: &[f64],
    equilibrium_constants: &[f64],
    fractions: &[f64],
    molar_masses: &[f64],
    density: f64,
    inter_fractions: &[f64],
    inter_density: f64,
    diffusion: &[f64],
) -> Result<KineticsResult> {
    let spec = &crate::model_gen::KINETICS_SPEC;
    let mut warnings = Vec::new();
    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "density" => Some(density),
            "inter_density" => Some(inter_density),
            _ => None,
        },
        &mut warnings,
    )?;

    let count = components.len();
    for (field, length) in [
        ("fractions", fractions.len()),
        ("molar_masses", molar_masses.len()),
        ("diffusion", diffusion.len()),
        ("inter_fractions", inter_fractions.len()),
    ] {
        if length != count {
            return Err(AzothError::InvalidInput {
                field: field.to_string(),
                reason: format!("{length} value(s) against {count} component(s)"),
            });
        }
    }
    if reaction_components.len() != reaction_coefficients.len() {
        return Err(AzothError::InvalidInput {
            field: "reaction_coefficients".to_string(),
            reason: format!(
                "{} name(s) and {} coefficient(s)",
                reaction_components.len(),
                reaction_coefficients.len()
            ),
        });
    }
    if rate_factors.len() != reaction_lengths.len()
        || equilibrium_constants.len() != reaction_lengths.len()
    {
        return Err(AzothError::InvalidInput {
            field: "rate_factors".to_string(),
            reason: format!(
                "{} reaction(s) and {} rate factor(s) against {} equilibrium constant(s)",
                reaction_lengths.len(),
                rate_factors.len(),
                equilibrium_constants.len()
            ),
        });
    }

    // Cut the concatenated lists into the reactions the class iterates, each with its own
    // names in its own order.
    let mut reactions = Vec::with_capacity(reaction_lengths.len());
    let mut cursor = 0usize;
    for (index, length) in reaction_lengths.iter().enumerate() {
        if !length.is_finite() || *length < 1.0 || length.fract() != 0.0 {
            return Err(AzothError::InvalidInput {
                field: "reaction_lengths".to_string(),
                reason: format!("reaction {index} names {length} species, which is not a count"),
            });
        }
        let length = *length as usize;
        if cursor + length > reaction_components.len() {
            return Err(AzothError::InvalidInput {
                field: "reaction_lengths".to_string(),
                reason: format!(
                    "reaction {index} runs past the {count}-species reaction list at {cursor}",
                    count = reaction_components.len()
                ),
            });
        }
        reactions.push(KineticReaction {
            names: reaction_components[cursor..cursor + length]
                .iter()
                .map(|name| (*name).to_string())
                .collect(),
            stoc_coefs: reaction_coefficients[cursor..cursor + length].to_vec(),
            rate_factor: rate_factors[index],
            equilibrium_constant: equilibrium_constants[index],
        });
        cursor += length;
    }
    if cursor != reaction_components.len() {
        return Err(AzothError::InvalidInput {
            field: "reaction_lengths".to_string(),
            reason: format!(
                "the lengths sum to {cursor} and the reaction list has {} species",
                reaction_components.len()
            ),
        });
    }

    let names: Vec<String> = components.iter().map(|name| (*name).to_string()).collect();
    let phase = KineticPhase {
        names: names.clone(),
        fractions: fractions.to_vec(),
        molar_masses: molar_masses.to_vec(),
        density,
    };
    let inter_phase = KineticPhase {
        names,
        fractions: inter_fractions.to_vec(),
        molar_masses: molar_masses.to_vec(),
        density: inter_density,
    };

    let mut result = KineticsResult {
        coefficient: Vec::with_capacity(count),
        phi_infinite: Vec::with_capacity(count),
        irreversible: Vec::with_capacity(count),
        warnings,
    };
    for name in components {
        let row = rate_matrix(&reactions, &phase, &inter_phase, name, diffusion)?;
        result.coefficient.push(row.coefficient);
        result.phi_infinite.push(row.phi_infinite.unwrap_or(0.0));
        result
            .irreversible
            .push(if row.irreversible { 1.0 } else { 0.0 });
    }
    Ok(result)
}
