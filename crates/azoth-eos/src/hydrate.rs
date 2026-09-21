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
pub fn langmuir(guest: &HydrateGuest, structure: usize, cavity: usize, t: f64) -> f64 {
    guest.langmuir_a[structure][cavity] / t * (guest.langmuir_b[structure][cavity] / t).exp()
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
    guests: &[HydrateGuest],
    ref_fugacities: &[f64],
    structure: usize,
    cavity: usize,
    t: f64,
) -> Vec<f64> {
    const PA_PER_BAR: f64 = 1.0e5;
    let mut denominator = 1.0;
    for (guest, fugacity) in guests.iter().zip(ref_fugacities) {
        if guest.former {
            denominator += langmuir(guest, structure, cavity, t) * fugacity / PA_PER_BAR;
        }
    }
    guests
        .iter()
        .zip(ref_fugacities)
        .map(|(guest, fugacity)| {
            if guest.former {
                langmuir(guest, structure, cavity, t) * fugacity / PA_PER_BAR / denominator
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

/// The exponent of a structure's water fugacity coefficient: the part that is its own.
///
/// ```text
/// f_w^hydrate/P = (f_w^ref/P) exp( sum_cav n_cav ln(1 - sum_j theta_j) + dMu )
/// ```
///
/// **The fluid's own water fugacity cancels out of this**, which is why it is not an
/// argument. NeqSim writes the coefficient as
/// `(f_w^fluid/P) exp(sum + dMu + ln(f_w^ref/f_w^fluid))`, where the fluid's water fugacity is
/// multiplied in and divided out again inside the logarithm - it is not what the hydrate's
/// water fugacity depends on, and `f_w^ref` is. Written out that way the coefficient is
/// `0 * inf` on a fluid with no water, which is a state a hydrate fraction's *bound* is, so
/// the cancellation is taken here and the structure comparison below is a comparison of
/// finite exponents rather than of `NaN`s.
///
/// Filled with the three terms, the coefficient at the probe's first state is `1.8831197e-4`,
/// the capture's own, with the cavity sum `-0.577685053` and the chemical-potential change
/// `+0.577684724` cancelling to a part in `1e7`.
///
/// # Errors
/// * [`AzothError::OutOfRange`] if a cavity type is fully occupied, where `ln(1 - sum)` has
///   no value. That is a real refusal rather than a clamp: a saturation of one means the
///   occupancy model has left the region it was fitted to.
fn exponent(
    guests: &[HydrateGuest],
    ref_fugacities: &[f64],
    t: f64,
    p: f64,
    structure: usize,
) -> Result<f64, AzothError> {
    let mut val = 0.0;
    for (cavity, per_water) in CAVITIES_PER_WATER[structure].iter().enumerate() {
        let occupied: f64 = occupancy(guests, ref_fugacities, structure, cavity, t)
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

    Ok(val + chemical_potential(structure, t, p))
}

/// The hydrate's water fugacity coefficient on one structure, `f_w^hydrate/P`.
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
/// # Errors
/// * [`AzothError::OutOfRange`] from the structure's cavity sum.
pub fn water_fugacity_coefficient(
    guests: &[HydrateGuest],
    ref_fugacities: &[f64],
    t: f64,
    p: f64,
    structure: usize,
    reference_water_fugacity: f64,
) -> Result<f64, AzothError> {
    Ok(exponent(guests, ref_fugacities, t, p, structure)?.exp() * reference_water_fugacity / p)
}

/// The stable structure and its water fugacity coefficient.
///
/// NeqSim evaluates both structures and keeps the **lower** coefficient, which is the more
/// stable hydrate at that state. The reference term is the same for both - it is the
/// *reference fluid's* fugacity over the pressure, and neither depends on the structure - so
/// the comparison is the exponents'.
///
/// # Errors
/// * [`AzothError::OutOfRange`] from either structure's cavity sum.
pub fn stable_structure(
    guests: &[HydrateGuest],
    ref_fugacities: &[f64],
    t: f64,
    p: f64,
    reference_water_fugacity: f64,
) -> Result<(usize, f64), AzothError> {
    let first = exponent(guests, ref_fugacities, t, p, 0)?;
    let second = exponent(guests, ref_fugacities, t, p, 1)?;
    let (structure, chosen) = if first <= second {
        (0, first)
    } else {
        (1, second)
    };
    Ok((structure, chosen.exp() * reference_water_fugacity / p))
}

/// A gas the mixture's fluid carries, and whether the hydrate takes it into a cage.
#[derive(Debug, Clone, PartialEq)]
pub struct HydrateGuest {
    /// The substance's name, as the databank spells it.
    pub name: String,
    /// The Langmuir pair, `[structure][cavity]`, in K.
    pub langmuir_a: [[f64; 2]; 2],
    /// The same pair's `B`, in K.
    pub langmuir_b: [[f64; 2]; 2],
    /// Whether it occupies a cage at all.
    pub former: bool,
}

/// What a mixture needs to have a hydrate calculated for it, carried beside the fluid.
///
/// The `Mixture` holds critical constants and nothing else, so a substance's *name* - which is
/// what the hydrate's tables are keyed by, and which decides whether it is a guest - has to
/// travel here. The same shape [`crate::furst_electrolyte`] uses, and for the same reason.
#[derive(Debug, Clone, PartialEq)]
pub struct Hydration {
    /// One entry per component, in the fluid's own order.
    pub guests: Vec<HydrateGuest>,
    /// Where water is, or `None` if the fluid has none - and a hydrate needs one.
    pub water_index: Option<usize>,
}

/// A fluid with the hydrate's tables attached, built from names.
///
/// # Errors
/// * [`AzothError::InvalidInput`] if a name is not in the databank, if there is no water, or
///   if nothing in the mixture is a hydrate former.
pub fn hydrate_mixture_of(
    names: &[&str],
    cubic: crate::Cubic,
    overlay: Option<&crate::databank::Overlay>,
) -> azoth_core::Result<(
    crate::mixture::Mixture,
    crate::molar_enthalpy_entropy::IdealGasModel,
)> {
    let entries: Vec<Entry> = names
        .iter()
        .map(|name| crate::databank::entry(name, overlay))
        .collect::<azoth_core::Result<Vec<_>>>()?;

    let guests: Vec<HydrateGuest> = entries
        .iter()
        .map(|entry| HydrateGuest {
            name: entry.name.clone(),
            langmuir_a: entry.hydrate_langmuir_a,
            langmuir_b: entry.hydrate_langmuir_b,
            former: entry.hydrate_former,
        })
        .collect();

    let water_index = entries.iter().position(|entry| entry.name == "water");
    if water_index.is_none() {
        return Err(AzothError::invalid_input(
            "components",
            "no water: a hydrate is water's, so a fluid without it has no formation \
             temperature"
                .to_string(),
        ));
    }
    if !guests.iter().any(|guest| guest.former) {
        return Err(AzothError::invalid_input(
            "components",
            "nothing in this mixture is a hydrate former, so no cage would have a guest and \
             the hydrate's fugacity would be its empty one"
                .to_string(),
        ));
    }

    let (mixture, ideal_gas) = crate::databank::mixture_of(names, cubic, overlay)?;
    Ok((
        mixture.with_hydration(Hydration {
            guests,
            water_index,
        }),
        ideal_gas,
    ))
}

/// The cavities and the water molecules of one unit cell: structure I then II.
///
/// NeqSim's `ComponentHydrate` constructor, and the counts its `46/54` and `136/160` bounds
/// come from - which are the *fully occupied* limits, not the composition at a state.
pub const CAVITIES_PER_CELL: [[f64; 2]; 2] = [[2.0, 6.0], [16.0, 8.0]];

/// The water molecules in one unit cell.
pub const WATER_PER_CELL: [f64; 2] = [46.0, 136.0];

/// The hydrate's mole fractions at a state, from **both** cavity types of a structure.
///
/// A cell of structure I holds `46` waters with two small and six large cavities, and one of
/// structure II `136` with sixteen and eight. A cavity of type `cav` holds guest `i` with
/// probability `theta_icav`, so the cell's guest count is the cavities' own weighted sum of
/// the occupancies, and each guest's is the one weighted by its own.
///
/// **The second cavity is the whole point.** A guest occupies one cage type at a time and its
/// occupancy in each is a different function of the state, so a composition built from one
/// cavity type leaves the other's guests at zero. For structure I that is the large
/// `5^12 6^2` cage, where most of the ethane and propane sit: the phase's fractions would not
/// sum to one, and its bound would not move with the temperature.
#[must_use]
pub fn composition(
    guests: &[HydrateGuest],
    ref_fugacities: &[f64],
    structure: usize,
    t: f64,
    water_index: usize,
) -> Vec<f64> {
    let mut counts = vec![0.0; guests.len()];
    let mut total_guests = 0.0;
    for (cavity, &per_cell) in CAVITIES_PER_CELL[structure].iter().enumerate() {
        for (index, occupied) in occupancy(guests, ref_fugacities, structure, cavity, t)
            .iter()
            .enumerate()
        {
            counts[index] += per_cell * occupied;
            total_guests += per_cell * occupied;
        }
    }
    let water = WATER_PER_CELL[structure];
    let total = water + total_guests;
    let mut fractions: Vec<f64> = counts.iter().map(|count| count / total).collect();
    fractions[water_index] = water / total;
    fractions
}
