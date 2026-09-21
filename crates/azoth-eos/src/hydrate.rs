//! The van der Waals-Platteeuw hydrate fugacity, in the two fitted forms NeqSim carries.
//!
//! Both build a guest's Langmuir constant as `C = A/T exp(B/T)` and take the water fugacity
//! coefficient from the cavity occupancies; they differ in which fitted pair, whether a
//! chemical-potential constant is added, and what the reference term is. [`HydrateModel`] is
//! the choice, and NeqSim makes it from the *system's model name* rather than from the fluid.
//!
//! - `ComponentHydratePVTsim` (the default, [`HydrateModel::Pvtsim`]) reads the eight
//!   `HydrateA*`/`B*` columns, adds `calcDeltaChemPot`, and references the **reference fluid's**
//!   own water fugacity.
//! - `ComponentHydrateGF` ([`HydrateModel::GuoFinch`]) reads the eight `A*_GF`/`B*_GF` columns,
//!   adds nothing, and references the **empty lattice's** vapour pressure on that structure.
//!
//! NeqSim's **base** `ComponentHydrate` integrates a Kihara potential over the guest's
//! Lennard-Jones parameters instead, and both of the above override that away - which is why
//! the three `*HYDRATE` LJ columns stay `vendored` in the ledger while the other sixteen are
//! `used`.
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

/// The reference pressure the empty-lattice vapour pressure is stated against, bar. NeqSim's
/// `ThermodynamicConstantsInterface.referencePressure`.
const REFERENCE_PRESSURE_BAR: f64 = 1.01325;

/// The empty structure's vapour pressure constants, structure I then II: Sloan (1990), as
/// NeqSim's `ComponentHydrate.emptyHydrateVapourPressureConstant`.
const EMPTY_VAPOUR_PRESSURE: [[f64; 2]; 2] = [[17.44, -6003.9], [17.332, -6017.6]];

/// Avlonitis (1994)'s molar volume of the empty hydrate, `(v0, k1, k2, k3)` in cm3/mol and
/// per K, structure I then II. NeqSim's `ComponentHydrate.getMolarVolumeHydrate`.
const HYDRATE_VOLUME: [[f64; 4]; 2] = [
    [22.35, 3.1075e-4, 5.9537e-7, 1.3707e-10],
    [22.57, 1.9335e-4, 2.1768e-7, -1.4786e-10],
];

/// Which fitted hydrate model a calculation runs.
///
/// The two differ in three places and nowhere else: which Langmuir pair a guest's constant is
/// built from, whether the chemical-potential constant is added, and what multiplies the
/// exponent to make a coefficient. NeqSim picks between them by the *system's model name*
/// (`PhaseHydrate`'s constructor), so the name is a model-level input rather than a property of
/// the fluid.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum HydrateModel {
    /// NeqSim's `ComponentHydratePVTsim`, the default: the fitted pair plus `calcDeltaChemPot`,
    /// against the reference *fluid's* water fugacity.
    #[default]
    Pvtsim,
    /// NeqSim's `ComponentHydrateGF`: the Guo-Finch pair, no chemical-potential constant, and
    /// the empty lattice's own vapour pressure as the reference.
    GuoFinch,
}

impl std::str::FromStr for HydrateModel {
    type Err = std::convert::Infallible;

    /// The names a spec carries. An unknown string is `Pvtsim`, which is what the databank's
    /// `eos` column does with an unknown cubic - a spec's own vocabulary check is what should
    /// refuse, not this.
    fn from_str(text: &str) -> Result<Self, Self::Err> {
        Ok(match text.trim().to_lowercase().as_str() {
            "guo_finch" | "guofinch" | "gf" => Self::GuoFinch,
            _ => Self::Pvtsim,
        })
    }
}

/// The molar volume of the empty hydrate on one structure, m3/mol. Avlonitis (1994).
///
/// The third value NeqSim carries (`v0 = 19.6522`) is the *ice* branch, reached only where
/// `hydrateStructure == -1`, which nothing in NeqSim sets.
#[must_use]
pub fn molar_volume_hydrate(structure: usize, t: f64) -> f64 {
    let [v0, k1, k2, k3] = HYDRATE_VOLUME[structure];
    let dt = t - T_REFERENCE;
    v0 * (1.0 + k1 * dt + k2 * dt * dt + k3 * dt * dt * dt) / 1.0e6
}

/// The empty hydrate's vapour pressure on one structure, **in bar**.
///
/// This is the Guo-Finch route's reference term, and it is the reason that route breaks the
/// structure comparison: it is a function of the structure, so the term does not cancel out of
/// a comparison of the two coefficients the way the PVTsim route's reference does.
#[must_use]
pub fn empty_hydrate_vapour_pressure(structure: usize, t: f64) -> f64 {
    let [c0, c1] = EMPTY_VAPOUR_PRESSURE[structure];
    (c0 + c1 / t).exp() * REFERENCE_PRESSURE_BAR
}

/// A guest's fitted pair at one structure and cavity, `C = A/T exp(B/T)`.
///
/// A cell left empty in the table is a zero, and a zero here means the guest does not occupy
/// that cavity at all: ethane's and propane's small-cavity columns are empty, which is why
/// their occupancies there are exactly zero rather than merely small.
#[must_use]
pub fn langmuir(
    guest: &HydrateGuest,
    model: HydrateModel,
    structure: usize,
    cavity: usize,
    t: f64,
) -> f64 {
    let (a, b) = match model {
        HydrateModel::Pvtsim => (guest.langmuir_a, guest.langmuir_b),
        HydrateModel::GuoFinch => (guest.guo_finch_a, guest.guo_finch_b),
    };
    a[structure][cavity] / t * (b[structure][cavity] / t).exp()
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
    model: HydrateModel,
    structure: usize,
    cavity: usize,
    t: f64,
) -> Vec<f64> {
    const PA_PER_BAR: f64 = 1.0e5;
    let mut denominator = 1.0;
    for (guest, fugacity) in guests.iter().zip(ref_fugacities) {
        if guest.former {
            denominator += langmuir(guest, model, structure, cavity, t) * fugacity / PA_PER_BAR;
        }
    }
    guests
        .iter()
        .zip(ref_fugacities)
        .map(|(guest, fugacity)| {
            if guest.former {
                langmuir(guest, model, structure, cavity, t) * fugacity / PA_PER_BAR / denominator
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
    model: HydrateModel,
    t: f64,
    p: f64,
    structure: usize,
) -> Result<f64, AzothError> {
    let mut val = 0.0;
    for (cavity, per_water) in CAVITIES_PER_WATER[structure].iter().enumerate() {
        let occupied: f64 = occupancy(guests, ref_fugacities, model, structure, cavity, t)
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

    Ok(match model {
        HydrateModel::Pvtsim => val + chemical_potential(structure, t, p),
        // The Guo-Finch route's reference is the empty lattice's own vapour pressure, carried
        // by `reference_term` below, so there is no chemical-potential constant here.
        HydrateModel::GuoFinch => val,
    })
}

/// What multiplies `exp(exponent)` to make the coefficient, `[structure][cavity]` aside.
///
/// The two models state their reference differently and the difference is the whole of the
/// structure comparison's shape. `Pvtsim`'s reference is the *reference fluid's* water
/// fugacity over the pressure, which is the same number for both structures - so a comparison
/// of the coefficients is a comparison of the exponents. `GuoFinch`'s is the empty lattice's
/// vapour pressure **on that structure**, which is not, and the pressure term that goes with
/// it is a Poynting correction about that vapour pressure rather than about the system's.
fn reference_term(
    model: HydrateModel,
    structure: usize,
    t: f64,
    p: f64,
    reference_water_fugacity: f64,
) -> f64 {
    match model {
        HydrateModel::Pvtsim => reference_water_fugacity / p,
        HydrateModel::GuoFinch => {
            let p_empty = empty_hydrate_vapour_pressure(structure, t) * 1.0e5;
            p_empty * (molar_volume_hydrate(structure, t) / (R * t) * (p - p_empty)).exp() / p
        }
    }
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
    model: HydrateModel,
    t: f64,
    p: f64,
    structure: usize,
    reference_water_fugacity: f64,
) -> Result<f64, AzothError> {
    Ok(
        exponent(guests, ref_fugacities, model, t, p, structure)?.exp()
            * reference_term(model, structure, t, p, reference_water_fugacity),
    )
}

/// The stable structure and its water fugacity coefficient.
///
/// NeqSim evaluates both structures and keeps the **lower** coefficient, which is the more
/// stable hydrate at that state.
///
/// **The comparison is of the coefficients, not of the exponents**, and for the Guo-Finch
/// route that is not a distinction without a difference: its reference term is the empty
/// lattice's vapour pressure on that structure, so the two structures do not share a factor and
/// an exponent comparison would pick the wrong one wherever the vapour pressures differ by more
/// than the cavity sums do. The PVTsim route's reference is the same number for both, which is
/// why an exponent comparison happens to agree there.
///
/// # Errors
/// * [`AzothError::OutOfRange`] from either structure's cavity sum.
pub fn stable_structure(
    guests: &[HydrateGuest],
    ref_fugacities: &[f64],
    model: HydrateModel,
    t: f64,
    p: f64,
    reference_water_fugacity: f64,
) -> Result<(usize, f64), AzothError> {
    let first = water_fugacity_coefficient(
        guests,
        ref_fugacities,
        model,
        t,
        p,
        0,
        reference_water_fugacity,
    )?;
    let second = water_fugacity_coefficient(
        guests,
        ref_fugacities,
        model,
        t,
        p,
        1,
        reference_water_fugacity,
    )?;
    Ok(if first <= second {
        (0, first)
    } else {
        (1, second)
    })
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
    /// The Guo-Finch pair, `[structure][cavity]`, in K and in K. Same shape, different fit.
    pub guo_finch_a: [[f64; 2]; 2],
    /// The same pair's `B`, in K.
    pub guo_finch_b: [[f64; 2]; 2],
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
    /// Which fitted model the guests' constants are read with.
    pub model: HydrateModel,
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
    model: HydrateModel,
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
            guo_finch_a: entry.hydrate_guo_finch_a,
            guo_finch_b: entry.hydrate_guo_finch_b,
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
            model,
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
    model: HydrateModel,
    structure: usize,
    t: f64,
    water_index: usize,
) -> Vec<f64> {
    let mut counts = vec![0.0; guests.len()];
    let mut total_guests = 0.0;
    for (cavity, &per_cell) in CAVITIES_PER_CELL[structure].iter().enumerate() {
        for (index, occupied) in occupancy(guests, ref_fugacities, model, structure, cavity, t)
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
