//! The van der Waals-Platteeuw hydrate fugacity, as NeqSim's `ComponentHydratePVTsim`
//! builds it.
//!
//! The route this carries is the fitted one: a guest's Langmuir constant is
//! `C = A/T exp(B/T)` from the eight `HydrateA*`/`B*` columns, and the water fugacity
//! coefficient follows from the cavity occupancies and a per-structure chemical-potential
//! constant. NeqSim's **base** `ComponentHydrate` integrates a Kihara potential over the
//! guest's Lennard-Jones parameters instead, and `ComponentHydratePVTsim` overrides that
//! away - which is why the three `*HYDRATE` LJ columns stay `vendored` in the ledger while
//! these eight are `used`.
//!
//! # The two structures
//!
//! Structure I has two small 5^12 cavities and six large 5^12 6^2 per forty-six waters;
//! structure II sixteen 5^12 and eight 5^12 6^4 per a hundred and thirty-six. NeqSim keeps
//! those as arrays in the component's constructor and takes the **lower** of the two
//! structures' water fugacity coefficients as the stable one, which is what
//! [`water_fugacity_coefficient`] reports.
//!
//! # Units
//!
//! Fugacities here are in **pascals**, and the coefficient is `f_w/P` and so dimensionless.
//! NeqSim works in bara; every place that matters is either a ratio or carries its own
//! conversion, and the chemical-potential constant's volume term is `dV P` with `P` the
//! pressure in pascals - which is what its own `pres * 1e5` is doing.

use azoth_core::AzothError;

use crate::databank::Entry;

/// NeqSim's `ThermodynamicConstantsInterface.R`, which `calcDeltaChemPot` divides by.
pub const R: f64 = 8.314_462_1;

/// The reference temperature of the chemical-potential constant, K.
const T_REFERENCE: f64 = 273.15;

/// How many cavities of each type a water molecule belongs to: structure I then II, small
/// then large. NeqSim's `ComponentHydrate` constructor, and hardcoded rather than read.
pub const CAVITIES_PER_WATER: [[f64; 2]; 2] = [[1.0 / 23.0, 3.0 / 23.0], [2.0 / 17.0, 1.0 / 17.0]];

/// The chemical-potential constants `(dGf, dHf, Cp, dV)` of `calcDeltaChemPot`, structure I
/// then II, in J/mol, J/mol, J/(mol K) and m3/mol.
const CHEMICAL_POTENTIAL: [[f64; 4]; 2] = [
    [1264.0, -4858.0, -39.16, 4.6e-6],
    [883.0, -5201.0, -39.16, 5.0e-6],
];

/// A guest's fitted pair at one structure and cavity, `C = A/T exp(B/T)`.
///
/// A cell left empty in the table is a zero, and a zero here means the guest does not occupy
/// that cavity at all: ethane's and propane's small-cavity columns are empty, which is why
/// their occupancies there are exactly zero rather than merely small.
#[must_use]
pub fn langmuir(entry: &Entry, structure: usize, cavity: usize, t: f64) -> f64 {
    entry.hydrate_langmuir_a[structure][cavity] / t
        * (entry.hydrate_langmuir_b[structure][cavity] / t).exp()
}

/// The occupancy of one cavity type by one guest, `C f_i / (1 + sum_j C_j f_j)`.
///
/// The sum runs over the guests the table marks as formers; water is not one, and NeqSim
/// excludes it by name in the same loop.
///
/// **The product `C f` is taken in bar**, whatever unit the caller states its fugacities in:
/// the Langmuir constant is fitted against neqsim's own bar-valued fugacities, and `C` carries
/// the reciprocal, so this is the one place a unit reaches the arithmetic rather than cancelling.
/// Feeding it pascals saturates every cage - the occupancies go to one and the answer the
/// solver finds is a different one - which is how this was found.
#[must_use]
pub fn occupancy(
    entries: &[&Entry],
    ref_fugacities: &[f64],
    structure: usize,
    cavity: usize,
    t: f64,
) -> Vec<f64> {
    const PA_PER_BAR: f64 = 1.0e5;
    let mut denominator = 1.0;
    for (entry, fugacity) in entries.iter().zip(ref_fugacities) {
        if entry.hydrate_former {
            denominator += langmuir(entry, structure, cavity, t) * fugacity / PA_PER_BAR;
        }
    }
    entries
        .iter()
        .zip(ref_fugacities)
        .map(|(entry, fugacity)| {
            if entry.hydrate_former {
                langmuir(entry, structure, cavity, t) * fugacity / PA_PER_BAR / denominator
            } else {
                0.0
            }
        })
        .collect()
}

/// `calcDeltaChemPot`: the change in chemical potential on forming the hydrate, dimensionless.
fn chemical_potential(structure: usize, t: f64, p: f64) -> f64 {
    let [dgf, dhf, cp, dvolume] = CHEMICAL_POTENTIAL[structure];
    dgf / R / T_REFERENCE
        - (-dhf * (1.0 / R / t - 1.0 / R / T_REFERENCE)
            + cp / R * (t / T_REFERENCE).ln()
            + cp * T_REFERENCE / R * (1.0 / t - 1.0 / T_REFERENCE))
        + dvolume / R / t * p
}

/// The hydrate's water fugacity coefficient on one structure.
///
/// ```text
/// f_w/P = (f_w^fluid/P) exp( sum_cav n_cav ln(1 - sum_j theta_j) + dMu + ln(f_w^ref/f_w^fluid) )
/// ```
///
/// **`f_w^ref` is the reference water phase's fugacity and is *not* the pressure.** NeqSim
/// builds that reference as a one-component phase of the **host's own class** -
/// `setSolidRefFluidPhase` clones the phase type the hydrate is attached to - so its water
/// component is the host's (`ComponentSrk` for an SRK fluid) and its fugacity is a real
/// number: at the probe's first state `0.0188312033887596` bar against the fluid's
/// `0.0188311981257476`, which is why the term nearly vanishes and why substituting the
/// pressure for it put the coefficient out by a factor of `P/f_w^ref`. The caller states it,
/// because only the caller knows which equation the fluid is.
///
/// Filled with the three terms, the coefficient at that state is `1.8831197e-4` - the
/// capture's own - with the cavity sum `-0.577685053` and the chemical-potential change
/// `+0.577684724` cancelling to a part in `1e7`.
///
/// # Errors
/// * [`AzothError::OutOfRange`] if a cavity type is fully occupied, where `ln(1 - sum)` has
///   no value. That is a real refusal rather than a clamp: a saturation of one means the
///   occupancy model has left the region it was fitted to.
pub fn water_fugacity_coefficient(
    entries: &[&Entry],
    ref_fugacities: &[f64],
    t: f64,
    p: f64,
    water_index: usize,
    structure: usize,
    reference_water_fugacity: f64,
) -> Result<f64, AzothError> {
    let mut val = 0.0;
    for (cavity, per_water) in CAVITIES_PER_WATER[structure].iter().enumerate() {
        let occupied: f64 = occupancy(entries, ref_fugacities, structure, cavity, t)
            .iter()
            .sum();
        if occupied >= 1.0 {
            return Err(AzothError::out_of_range(
                "occupancy",
                occupied,
                format!(
                    "cavity {cavity} of structure {} is fully occupied at {t} K, where the \
                     cavity sum's logarithm has no value",
                    structure + 1
                ),
            ));
        }
        val += per_water * (1.0 - occupied).ln();
    }

    let alpha_water = ref_fugacities[water_index];
    let water_alpha_ref = (reference_water_fugacity / alpha_water).ln();

    Ok(alpha_water * (val + chemical_potential(structure, t, p) + water_alpha_ref).exp() / p)
}

/// The stable structure and its water fugacity coefficient.
///
/// NeqSim evaluates both structures and keeps the **lower** coefficient, which is the more
/// stable hydrate at that state.
///
/// # Errors
/// * [`AzothError::OutOfRange`] from either structure's cavity sum.
pub fn stable_structure(
    entries: &[&Entry],
    ref_fugacities: &[f64],
    t: f64,
    p: f64,
    water_index: usize,
    reference_water_fugacity: f64,
) -> Result<(usize, f64), AzothError> {
    let first = water_fugacity_coefficient(
        entries,
        ref_fugacities,
        t,
        p,
        water_index,
        0,
        reference_water_fugacity,
    )?;
    let second = water_fugacity_coefficient(
        entries,
        ref_fugacities,
        t,
        p,
        water_index,
        1,
        reference_water_fugacity,
    )?;
    Ok(if first <= second {
        (0, first)
    } else {
        (1, second)
    })
}
