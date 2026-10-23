//! `PseudoComponentCombiner`'s arithmetic, less the fluid it usually takes.
//!
//! The class is 1841 lines and every public method takes or returns a `SystemInterface`. What
//! this module holds is the part that is a function of a table: extraction, the equal-mass cut
//! points, the distribution of each fluid's own rows onto them, and the blend of the per-fluid
//! profiles.
//!
//! **The arithmetic is the class's, down to the summation order.** A boundary is the sorting key
//! of the row whose running mass first reaches a target, not an interpolation; a row joins a
//! group by a strict `key > boundary`; a group's molar mass is its accumulated mass over its
//! accumulated moles and its density the mass over the volume its rows occupy. Two of the
//! combined properties - the boiling point and the three critical ones - are weighted by fluid
//! *mass* fraction times *mole* fraction within that fluid, where the other fifteen are
//! mass-weighted; that split is Pedersen's, and it is the whole of `combineProfiles`.

use azoth_core::{AzothError, Result};

use crate::component_row::{CharacterizationComponent, FluidComponent};
use crate::delumping::sorting_key;

/// `PseudoComponentCombiner.MASS_TOLERANCE`, on an accumulated mass or mole amount.
pub(crate) const MASS_TOLERANCE: f64 = 1.0e-12;

/// `Component.getCpD`'s own offset: the getter returns the stored coefficient plus this.
///
/// **The class reads a source through the getter and assigns through the setter**, so a
/// blend's coefficient is the mean of the getter's values, and the rebuilt component's getter
/// states that mean plus the offset again. A row carries the getter's number - which is what
/// `PseudoComponentContribution` read - and `assembly` stores it back exactly as `applyTo`
/// would. Every other getter the class reads is plain, so this is the one place the
/// distinction exists.
pub(crate) const CP_D_GETTER_OFFSET: f64 = 1.0e-10;

/// `!(value > tolerance)`, NaN included: the class's own comparison drops an accumulation it
/// cannot order rather than keeping it, where `value <= tolerance` would keep the NaN.
pub(crate) fn below(value: f64, tolerance: f64) -> bool {
    value.is_nan() || value <= tolerance
}

/// The seventeen properties a contribution carries beyond its moles, molar mass and density.
///
/// **An absent property is a NaN**, which is how NeqSim spells one in its row: every
/// `Double.isFinite` guard in the class is a test for it.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub(crate) struct Properties {
    pub(crate) normal_boiling_point: f64,
    pub(crate) critical_temperature: f64,
    pub(crate) critical_pressure: f64,
    pub(crate) acentric_factor: f64,
    pub(crate) critical_volume: f64,
    pub(crate) racket_z: f64,
    pub(crate) racket_z_cpa: f64,
    pub(crate) parachor: f64,
    pub(crate) critical_viscosity: f64,
    pub(crate) triple_point_temperature: f64,
    pub(crate) heat_of_fusion: f64,
    pub(crate) ideal_gas_enthalpy_of_formation: f64,
    pub(crate) cp_a: f64,
    pub(crate) cp_b: f64,
    pub(crate) cp_c: f64,
    pub(crate) cp_d: f64,
    pub(crate) attractive_m: f64,
}

impl Properties {
    /// The row's properties, with an absent field as a NaN.
    fn from_row(row: &FluidComponent) -> Self {
        let value = |field: Option<f64>| field.unwrap_or(f64::NAN);
        Self {
            normal_boiling_point: value(row.normal_boiling_point),
            critical_temperature: value(row.critical_temperature),
            critical_pressure: value(row.critical_pressure),
            acentric_factor: value(row.acentric_factor),
            critical_volume: value(row.critical_volume),
            racket_z: value(row.racket_z),
            racket_z_cpa: value(row.racket_z_cpa),
            parachor: value(row.parachor),
            critical_viscosity: value(row.critical_viscosity),
            triple_point_temperature: value(row.triple_point_temperature),
            heat_of_fusion: value(row.heat_of_fusion),
            ideal_gas_enthalpy_of_formation: value(row.ideal_gas_enthalpy_of_formation),
            cp_a: value(row.cp_a),
            cp_b: value(row.cp_b),
            cp_c: value(row.cp_c),
            cp_d: value(row.cp_d),
            attractive_m: value(row.attractive_m),
        }
    }

    /// The properties an array in [`Self::row_values`]'s order states.
    fn from_values(values: [f64; 17]) -> Self {
        let [
            normal_boiling_point,
            critical_temperature,
            critical_pressure,
            acentric_factor,
            critical_volume,
            racket_z,
            racket_z_cpa,
            parachor,
            critical_viscosity,
            triple_point_temperature,
            heat_of_fusion,
            ideal_gas_enthalpy_of_formation,
            cp_a,
            cp_b,
            cp_c,
            cp_d,
            attractive_m,
        ] = values;
        Self {
            normal_boiling_point,
            critical_temperature,
            critical_pressure,
            acentric_factor,
            critical_volume,
            racket_z,
            racket_z_cpa,
            parachor,
            critical_viscosity,
            triple_point_temperature,
            heat_of_fusion,
            ideal_gas_enthalpy_of_formation,
            cp_a,
            cp_b,
            cp_c,
            cp_d,
            attractive_m,
        }
    }

    /// Each property as a row states it: finite or absent.
    fn row_values(self) -> [f64; 17] {
        [
            self.normal_boiling_point,
            self.critical_temperature,
            self.critical_pressure,
            self.acentric_factor,
            self.critical_volume,
            self.racket_z,
            self.racket_z_cpa,
            self.parachor,
            self.critical_viscosity,
            self.triple_point_temperature,
            self.heat_of_fusion,
            self.ideal_gas_enthalpy_of_formation,
            self.cp_a,
            self.cp_b,
            self.cp_c,
            self.cp_d,
            self.attractive_m,
        ]
    }
}

/// One pseudo-component, as `extractComponents` lifts it.
#[derive(Debug, Clone)]
pub(crate) struct Contribution {
    pub(crate) moles: f64,
    pub(crate) molar_mass: f64,
    pub(crate) density: f64,
    pub(crate) properties: Properties,
    pub(crate) mass: f64,
}

impl Contribution {
    /// The row as a contribution, or the refusal for a pseudo-component that cannot be one.
    ///
    /// A molar mass and a density are required: the mass is `moles * molar_mass` and a group's
    /// volume is the mass over its density, so a row without either is a row the class could
    /// not have held in the first place.
    fn from_row(index: usize, row: &FluidComponent) -> Result<Self> {
        let Some(molar_mass) = row.molar_mass else {
            return Err(AzothError::invalid_input(
                format!("components[{index}].molar_mass"),
                "a pseudo-component carries its own molar mass, and this row states none",
            ));
        };
        let Some(density) = row.density else {
            return Err(AzothError::invalid_input(
                format!("components[{index}].density"),
                "a pseudo-component carries its own normal liquid density, and this row states \
                 none",
            ));
        };
        Ok(Self {
            moles: row.moles,
            molar_mass,
            density,
            properties: Properties::from_row(row),
            mass: row.moles * molar_mass,
        })
    }

    /// The contribution with its mole amount - and hence its mass - scaled, its position on the
    /// grid unchanged. The common slate's weighted imaginary pool is built from these.
    pub(crate) fn scaled_by_mass(&self, factor: f64) -> Self {
        let moles = self.moles * factor;
        Self {
            moles,
            molar_mass: self.molar_mass,
            density: self.density,
            properties: self.properties,
            mass: moles * self.molar_mass,
        }
    }

    /// `PseudoComponentContribution.sortingKey`: the boiling point, or the molar mass without a
    /// finite positive one.
    pub(crate) fn sorting_key(&self) -> f64 {
        sorting_key(Some(self.properties.normal_boiling_point), self.molar_mass)
    }
}

/// One fluid's extraction: its base components, merged by name, and its pseudo-components in
/// ascending sorting key.
#[derive(Debug, Clone)]
pub(crate) struct Extraction {
    pub(crate) base: Vec<(String, f64)>,
    pub(crate) pseudo: Vec<Contribution>,
}

/// `extractComponents`, for every fluid the rows name, in first-appearance order.
///
/// A row whose mole amount is not positive is skipped - the class's own `!(moles > 0)` - so a
/// zero is a row that is not on the grid rather than a row that contributes nothing.
pub(crate) fn extract(rows: &[FluidComponent]) -> Result<Vec<Extraction>> {
    let mut fluids: Vec<f64> = Vec::new();
    for (index, row) in rows.iter().enumerate() {
        let valid = row.fluid.is_finite() && row.fluid >= 0.0 && row.fluid.fract() == 0.0;
        if !valid {
            return Err(AzothError::invalid_input(
                format!("components[{index}].fluid"),
                "a fluid index is a whole number from zero, and this row's is not",
            ));
        }
        if !fluids.contains(&row.fluid) {
            fluids.push(row.fluid);
        }
    }

    let mut extractions: Vec<Extraction> = fluids
        .iter()
        .map(|_| Extraction {
            base: Vec::new(),
            pseudo: Vec::new(),
        })
        .collect();
    for (index, row) in rows.iter().enumerate() {
        if below(row.moles, 0.0) {
            continue;
        }
        let fluid = fluids
            .iter()
            .position(|candidate| *candidate == row.fluid)
            .expect("every row's fluid was collected above");
        if row.pseudo_component {
            extractions[fluid]
                .pseudo
                .push(Contribution::from_row(index, row)?);
        } else {
            let base = &mut extractions[fluid].base;
            match base.iter_mut().find(|(name, _)| *name == row.name) {
                Some((_, moles)) => *moles += row.moles,
                None => base.push((row.name.clone(), row.moles)),
            }
        }
    }

    for extraction in &mut extractions {
        extraction
            .pseudo
            .sort_by(|left, right| left.sorting_key().total_cmp(&right.sorting_key()));
    }
    Ok(extractions)
}

/// The base components of every fluid, merged by name in the fluids' own order.
pub(crate) fn merge_bases(extractions: &[Extraction]) -> Vec<(String, f64)> {
    let mut merged: Vec<(String, f64)> = Vec::new();
    for extraction in extractions {
        for (name, moles) in &extraction.base {
            match merged.iter_mut().find(|(known, _)| known == name) {
                Some((_, total)) => *total += moles,
                None => merged.push((name.clone(), *moles)),
            }
        }
    }
    merged
}

/// `determineQuantileBoundaries`: the sorting key of the row each equal-mass target falls in.
///
/// The targets are equal **masses**, `total * (i + 1) / target`, and the walk places a boundary
/// at the row whose running mass first reaches one - so a boundary is always a key some row
/// carries. A target nothing reaches is padded with an infinity, which `distribute` then never
/// crosses.
pub(crate) fn quantile_boundaries(contributions: &[Contribution], target: usize) -> Vec<f64> {
    if target <= 1 || contributions.is_empty() {
        return Vec::new();
    }
    let mut sorted: Vec<&Contribution> = contributions.iter().collect();
    sorted.sort_by(|left, right| left.sorting_key().total_cmp(&right.sorting_key()));

    let total_mass: f64 = sorted.iter().map(|contribution| contribution.mass).sum();
    if below(total_mass, MASS_TOLERANCE) {
        return Vec::new();
    }

    let targets: Vec<f64> = (0..target - 1)
        .map(|index| total_mass * (index + 1) as f64 / target as f64)
        .collect();
    let mut boundaries = Vec::with_capacity(targets.len());
    let mut cumulative = 0.0;
    let mut target_index = 0;
    for contribution in &sorted {
        let next = cumulative + contribution.mass;
        while target_index < targets.len() && next >= targets[target_index] - MASS_TOLERANCE {
            boundaries.push(contribution.sorting_key());
            target_index += 1;
        }
        cumulative = next;
    }
    while boundaries.len() < targets.len() {
        boundaries.push(f64::INFINITY);
    }
    boundaries
}

/// A weighted mean that ignores a non-finite sample, as the class's accumulators do.
#[derive(Debug, Clone, Copy, Default)]
struct Weighted {
    sum: f64,
    weight: f64,
    any: bool,
}

impl Weighted {
    fn add(&mut self, weight: f64, value: f64) {
        if value.is_finite() {
            self.sum += weight * value;
            self.weight += weight;
            self.any = true;
        }
    }

    /// The mean, or a NaN where no sample was finite - `hasX ? sum / mass : NaN`.
    fn mean(&self, mass: f64) -> f64 {
        if self.any { self.sum / mass } else { f64::NAN }
    }

    /// The weighted mean of the samples themselves, `sum(w X) / sum(w)`, or a NaN below the
    /// tolerance - `WeightedMean.mean()`.
    fn weighted_mean(&self) -> f64 {
        if self.weight > MASS_TOLERANCE {
            self.sum / self.weight
        } else {
            f64::NAN
        }
    }
}

/// `PseudoComponentGroupBuilder`: the running totals one group is built from.
#[derive(Debug, Clone, Default)]
struct GroupBuilder {
    mass: f64,
    moles: f64,
    volume: f64,
    density_mass: f64,
    normal_boiling_point: Weighted,
    critical_temperature: Weighted,
    critical_pressure: Weighted,
    acentric_factor: Weighted,
    critical_volume: Weighted,
    racket_z: Weighted,
    racket_z_cpa: Weighted,
    parachor: Weighted,
    critical_viscosity: Weighted,
    triple_point_temperature: Weighted,
    heat_of_fusion: Weighted,
    ideal_gas_enthalpy_of_formation: Weighted,
    cp_a: Weighted,
    cp_b: Weighted,
    cp_c: Weighted,
    cp_d: Weighted,
    attractive_m: Weighted,
}

impl GroupBuilder {
    /// Add a row's mass and its mass-weighted properties.
    fn add_contribution(&mut self, source: &Contribution, mass_portion: f64) {
        if below(mass_portion, MASS_TOLERANCE) || below(source.mass, 0.0) {
            return;
        }
        let ratio = mass_portion / source.mass;
        self.mass += mass_portion;
        self.moles += source.moles * ratio;
        if source.density > 0.0 {
            self.volume += mass_portion / source.density;
            self.density_mass += mass_portion * source.density;
        }
        let properties = source.properties;
        self.normal_boiling_point
            .add(mass_portion, properties.normal_boiling_point);
        self.critical_temperature
            .add(mass_portion, properties.critical_temperature);
        self.critical_pressure
            .add(mass_portion, properties.critical_pressure);
        self.acentric_factor
            .add(mass_portion, properties.acentric_factor);
        self.critical_volume
            .add(mass_portion, properties.critical_volume);
        self.racket_z.add(mass_portion, properties.racket_z);
        self.racket_z_cpa.add(mass_portion, properties.racket_z_cpa);
        self.parachor.add(mass_portion, properties.parachor);
        self.critical_viscosity
            .add(mass_portion, properties.critical_viscosity);
        self.triple_point_temperature
            .add(mass_portion, properties.triple_point_temperature);
        self.heat_of_fusion
            .add(mass_portion, properties.heat_of_fusion);
        self.ideal_gas_enthalpy_of_formation
            .add(mass_portion, properties.ideal_gas_enthalpy_of_formation);
        self.cp_a.add(mass_portion, properties.cp_a);
        self.cp_b.add(mass_portion, properties.cp_b);
        self.cp_c.add(mass_portion, properties.cp_c);
        self.cp_d.add(mass_portion, properties.cp_d);
        self.attractive_m.add(mass_portion, properties.attractive_m);
    }

    /// Add a whole profile, which is what blending the fluids does.
    fn add_profile(&mut self, profile: &Profile) {
        if !profile.is_valid() {
            return;
        }
        let mass_portion = profile.mass();
        self.mass += mass_portion;
        self.moles += profile.moles;
        if profile.density > 0.0 {
            self.volume += mass_portion / profile.density;
            self.density_mass += mass_portion * profile.density;
        }
        let properties = profile.properties;
        self.normal_boiling_point
            .add(mass_portion, properties.normal_boiling_point);
        self.critical_temperature
            .add(mass_portion, properties.critical_temperature);
        self.critical_pressure
            .add(mass_portion, properties.critical_pressure);
        self.acentric_factor
            .add(mass_portion, properties.acentric_factor);
        self.critical_volume
            .add(mass_portion, properties.critical_volume);
        self.racket_z.add(mass_portion, properties.racket_z);
        self.racket_z_cpa.add(mass_portion, properties.racket_z_cpa);
        self.parachor.add(mass_portion, properties.parachor);
        self.critical_viscosity
            .add(mass_portion, properties.critical_viscosity);
        self.triple_point_temperature
            .add(mass_portion, properties.triple_point_temperature);
        self.heat_of_fusion
            .add(mass_portion, properties.heat_of_fusion);
        self.ideal_gas_enthalpy_of_formation
            .add(mass_portion, properties.ideal_gas_enthalpy_of_formation);
        self.cp_a.add(mass_portion, properties.cp_a);
        self.cp_b.add(mass_portion, properties.cp_b);
        self.cp_c.add(mass_portion, properties.cp_c);
        self.cp_d.add(mass_portion, properties.cp_d);
        self.attractive_m.add(mass_portion, properties.attractive_m);
    }

    /// `buildWithOverrides`: the profile, or an empty one where the group carries no mass.
    fn build_with_overrides(
        &self,
        boiling_point: Option<f64>,
        critical_temperature: Option<f64>,
        critical_pressure: Option<f64>,
        acentric_factor: Option<f64>,
    ) -> Profile {
        if below(self.mass, MASS_TOLERANCE) || below(self.moles, MASS_TOLERANCE) {
            return Profile::empty();
        }
        let density = if self.volume > MASS_TOLERANCE {
            self.mass / self.volume
        } else if self.density_mass > MASS_TOLERANCE {
            self.density_mass / self.mass
        } else {
            f64::NAN
        };
        let mass = self.mass;
        Profile {
            moles: self.moles,
            molar_mass: self.mass / self.moles,
            density,
            properties: Properties {
                normal_boiling_point: boiling_point
                    .unwrap_or_else(|| self.normal_boiling_point.mean(mass)),
                critical_temperature: critical_temperature
                    .unwrap_or_else(|| self.critical_temperature.mean(mass)),
                critical_pressure: critical_pressure
                    .unwrap_or_else(|| self.critical_pressure.mean(mass)),
                acentric_factor: acentric_factor.unwrap_or_else(|| self.acentric_factor.mean(mass)),
                critical_volume: self.critical_volume.mean(mass),
                racket_z: self.racket_z.mean(mass),
                racket_z_cpa: self.racket_z_cpa.mean(mass),
                parachor: self.parachor.mean(mass),
                critical_viscosity: self.critical_viscosity.mean(mass),
                triple_point_temperature: self.triple_point_temperature.mean(mass),
                heat_of_fusion: self.heat_of_fusion.mean(mass),
                ideal_gas_enthalpy_of_formation: self.ideal_gas_enthalpy_of_formation.mean(mass),
                cp_a: self.cp_a.mean(mass),
                cp_b: self.cp_b.mean(mass),
                cp_c: self.cp_c.mean(mass),
                cp_d: self.cp_d.mean(mass),
                attractive_m: self.attractive_m.mean(mass),
            },
        }
    }
}

/// `PseudoComponentProfile`: one cut of a fluid, or of the blend.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct Profile {
    pub(crate) moles: f64,
    pub(crate) molar_mass: f64,
    pub(crate) density: f64,
    pub(crate) properties: Properties,
}

impl Profile {
    /// The class's empty profile: no moles, and every property absent.
    fn empty() -> Self {
        Self {
            moles: 0.0,
            molar_mass: 0.0,
            density: f64::NAN,
            properties: Properties {
                normal_boiling_point: f64::NAN,
                critical_temperature: f64::NAN,
                critical_pressure: f64::NAN,
                acentric_factor: f64::NAN,
                critical_volume: f64::NAN,
                racket_z: f64::NAN,
                racket_z_cpa: f64::NAN,
                parachor: f64::NAN,
                critical_viscosity: f64::NAN,
                triple_point_temperature: f64::NAN,
                heat_of_fusion: f64::NAN,
                ideal_gas_enthalpy_of_formation: f64::NAN,
                cp_a: f64::NAN,
                cp_b: f64::NAN,
                cp_c: f64::NAN,
                cp_d: f64::NAN,
                attractive_m: f64::NAN,
            },
        }
    }

    pub(crate) fn mass(&self) -> f64 {
        self.moles * self.molar_mass
    }

    /// `moles > 0 && molarMass > 0 && density > 0`: what makes a cut a component rather than a
    /// gap in the slate.
    pub(crate) fn is_valid(&self) -> bool {
        self.moles > 0.0 && self.molar_mass > 0.0 && self.density > 0.0
    }
}

/// `distributeToProfiles`: each row into the group its key falls in.
///
/// The walk is forward and the comparison is **strict**, so a key exactly on a boundary stays
/// in the lower group. `contributions` must already be in ascending sorting key.
pub(crate) fn distribute(
    contributions: &[Contribution],
    boundaries: &[f64],
    target: usize,
) -> Vec<Profile> {
    if target == 0 {
        return Vec::new();
    }
    if contributions.is_empty() {
        return vec![Profile::empty(); target];
    }
    let mut builders: Vec<GroupBuilder> = vec![GroupBuilder::default(); target];
    let mut group = 0;
    let mut current_boundary = boundaries.first().copied().unwrap_or(f64::INFINITY);
    for contribution in contributions {
        let key = contribution.sorting_key();
        while group < target - 1 && key > current_boundary {
            group += 1;
            current_boundary = boundaries.get(group).copied().unwrap_or(f64::INFINITY);
        }
        builders[group].add_contribution(contribution, contribution.mass);
    }
    builders
        .iter()
        .map(|builder| builder.build_with_overrides(None, None, None, None))
        .collect()
}

/// `combineProfiles`: the per-fluid cuts blended into one slate.
///
/// A fluid's weight is its share of the pooled pseudo-component **mass**, and within it a cut's
/// share is its mole fraction of that fluid. The boiling point and the three critical
/// properties are weighted by the product of the two; the remaining thirteen are mass-weighted
/// by the builder. A denominator at or below the tolerance leaves the property to the builder
/// rather than dividing by it.
pub(crate) fn combine_profiles(
    per_fluid_profiles: &[Vec<Profile>],
    fluid_mass_totals: &[f64],
    fluid_mole_totals: &[f64],
) -> Vec<Profile> {
    let Some(first) = per_fluid_profiles.first() else {
        return Vec::new();
    };
    let group_count = first.len();
    let total_fluid_mass: f64 = fluid_mass_totals.iter().sum();
    let mut combined = Vec::with_capacity(group_count);

    for group in 0..group_count {
        let mut builder = GroupBuilder::default();
        let mut boiling_numerator = 0.0;
        let mut boiling_denominator = 0.0;
        let mut temperature_numerator = 0.0;
        let mut temperature_denominator = 0.0;
        let mut pressure_numerator = 0.0;
        let mut pressure_denominator = 0.0;
        let mut acentric_numerator = 0.0;
        let mut acentric_denominator = 0.0;

        for (fluid, profiles) in per_fluid_profiles.iter().enumerate() {
            let profile = profiles[group];
            builder.add_profile(&profile);
            let fluid_mass = fluid_mass_totals[fluid];
            let fluid_moles = fluid_mole_totals[fluid];
            let fluid_fraction = if total_fluid_mass > MASS_TOLERANCE {
                fluid_mass / total_fluid_mass
            } else {
                0.0
            };
            let mole_fraction = if fluid_moles > MASS_TOLERANCE {
                profile.moles / fluid_moles
            } else {
                0.0
            };
            let weight = fluid_fraction * mole_fraction;
            let properties = profile.properties;
            if properties.normal_boiling_point.is_finite() {
                boiling_numerator += weight * properties.normal_boiling_point;
                boiling_denominator += weight;
            }
            if properties.critical_temperature.is_finite() {
                temperature_numerator += weight * properties.critical_temperature;
                temperature_denominator += weight;
            }
            if properties.critical_pressure.is_finite() {
                pressure_numerator += weight * properties.critical_pressure;
                pressure_denominator += weight;
            }
            if properties.acentric_factor.is_finite() {
                acentric_numerator += weight * properties.acentric_factor;
                acentric_denominator += weight;
            }
        }

        let over = |numerator: f64, denominator: f64| {
            if denominator > MASS_TOLERANCE {
                Some(numerator / denominator)
            } else {
                None
            }
        };
        combined.push(builder.build_with_overrides(
            over(boiling_numerator, boiling_denominator),
            over(temperature_numerator, temperature_denominator),
            over(pressure_numerator, pressure_denominator),
            over(acentric_numerator, acentric_denominator),
        ));
    }
    combined
}

/// `AttractiveTermSrk.setm`'s side effect: the acentric factor the fitted `m` implies.
///
/// **The setter rewrites the component's acentric factor**, solving
/// `-0.176 w^2 + 1.574 w + (0.48 - m) = 0` with the class's own Newton-Halley iteration from a
/// 0.25 seed. `applyTo` assigns the averaged acentric factor first and `m` last, so where the row
/// carries an `m` the answer's omega is this root and not the average - and it is the solver's
/// number rather than a closed-form root.
fn implied_acentric_factor(m: f64) -> f64 {
    const EPS: f64 = 1.0e-10;
    let value = |x: f64| -0.176 * x * x + 1.574 * x + (0.48 - m);
    let derivative = |x: f64| -0.352 * x + 1.574;
    let second = |_x: f64| -0.352;
    let mut x = 0.25;
    let mut next = x;
    let mut alternate = x;
    let mut iterations = 0_u32;
    loop {
        iterations += 1;
        x = if value(next).abs() < value(alternate).abs() {
            next
        } else {
            alternate
        };
        let (slope, curvature) = (derivative(x), second(x));
        let radical = slope * slope - 2.0 * value(x) * curvature;
        if radical > 0.0 {
            next = x - slope / curvature + radical.sqrt() / curvature;
            alternate = x - slope / curvature - radical.sqrt() / curvature;
        } else {
            next = x - value(x) / slope;
            alternate = next;
        }
        if next < -0.99 {
            next = -0.99;
        }
        if next > 1.5 {
            next = 1.0;
        }
        if below(value(x).abs(), EPS) || iterations > 500 {
            break;
        }
    }
    next
}

/// `computeCommonSlate`: the shared slate, by Pedersen's Chapter 5.6 weighting.
///
/// Each fluid's cut contributes weighted by the fluid's slate weight times the cut's mole
/// fraction **of that fluid's total** - base components included, which is what makes a dilute
/// fluid's pseudo-components count for less than its weight alone would say. The molar mass is
/// the weighted mean of the cuts' own; the density is that mean over the weighted mean of their
/// molar volumes, so a cut that states no positive density contributes to neither.
pub(crate) fn compute_common_slate(
    per_fluid_profiles: &[Vec<Profile>],
    fluid_mole_totals: &[f64],
    weights: &[f64],
    target: usize,
) -> Vec<Profile> {
    let mut slate = Vec::with_capacity(target);
    for index in 0..target {
        let mut molar_mass = Weighted::default();
        let mut molar_volume = Weighted::default();
        let mut properties = [Weighted::default(); 17];
        let mut slate_moles = 0.0;
        for (fluid, profiles) in per_fluid_profiles.iter().enumerate() {
            let profile = profiles[index];
            if !profile.is_valid() {
                continue;
            }
            let total = fluid_mole_totals[fluid];
            let mole_fraction = if total > MASS_TOLERANCE {
                profile.moles / total
            } else {
                0.0
            };
            let weight = weights[fluid] * mole_fraction;
            if below(weight, 0.0) {
                continue;
            }
            slate_moles += profile.moles;
            molar_mass.add(weight, profile.molar_mass);
            if profile.density > 0.0 {
                molar_volume.add(weight, profile.molar_mass / profile.density);
            }
            for (slot, value) in properties.iter_mut().zip(profile.properties.row_values()) {
                slot.add(weight, value);
            }
        }
        let mean_molar_mass = molar_mass.weighted_mean();
        let mean_molar_volume = molar_volume.weighted_mean();
        let density = if mean_molar_mass.is_finite()
            && mean_molar_volume.is_finite()
            && mean_molar_volume > 0.0
        {
            mean_molar_mass / mean_molar_volume
        } else {
            f64::NAN
        };
        slate.push(Profile {
            moles: slate_moles,
            molar_mass: if mean_molar_mass.is_finite() {
                mean_molar_mass
            } else {
                0.0
            },
            density,
            properties: Properties::from_values(properties.map(|slot| slot.weighted_mean())),
        });
    }
    slate
}

/// A profile as a result row, a finite property present and an absent one left out.
pub(crate) fn profile_row(name: String, profile: &Profile) -> CharacterizationComponent {
    let finite = |value: f64| value.is_finite().then_some(value);
    let properties = profile.properties;
    let [
        normal_boiling_point,
        critical_temperature,
        critical_pressure,
        acentric_factor,
        critical_volume,
        racket_z,
        racket_z_cpa,
        parachor,
        critical_viscosity,
        triple_point_temperature,
        heat_of_fusion,
        ideal_gas_enthalpy_of_formation,
        cp_a,
        cp_b,
        cp_c,
        cp_d,
        attractive_m,
    ] = properties.row_values();
    CharacterizationComponent {
        name,
        pseudo_component: true,
        moles: profile.moles,
        molar_mass: finite(profile.molar_mass),
        density: finite(profile.density),
        normal_boiling_point: finite(normal_boiling_point),
        critical_temperature: finite(critical_temperature),
        critical_pressure: finite(critical_pressure),
        // `setm` runs after `setAcentricFactor` in `applyTo`, so the fitted coefficient's root
        // is what the rebuilt component's getter states.
        acentric_factor: finite(if attractive_m.is_finite() {
            implied_acentric_factor(attractive_m)
        } else {
            acentric_factor
        }),
        critical_volume: finite(critical_volume),
        racket_z: finite(racket_z),
        racket_z_cpa: finite(racket_z_cpa),
        parachor: finite(parachor),
        critical_viscosity: finite(critical_viscosity),
        triple_point_temperature: finite(triple_point_temperature),
        heat_of_fusion: finite(heat_of_fusion),
        ideal_gas_enthalpy_of_formation: finite(ideal_gas_enthalpy_of_formation),
        cp_a: finite(cp_a),
        cp_b: finite(cp_b),
        cp_c: finite(cp_c),
        cp_d: finite(cp_d + CP_D_GETTER_OFFSET),
        attractive_m: finite(attractive_m),
    }
}

/// A base component as a result row: its name and amount, every property absent.
pub(crate) fn base_row(name: String, moles: f64) -> CharacterizationComponent {
    CharacterizationComponent {
        name,
        pseudo_component: false,
        moles,
        molar_mass: None,
        density: None,
        normal_boiling_point: None,
        critical_temperature: None,
        critical_pressure: None,
        acentric_factor: None,
        critical_volume: None,
        racket_z: None,
        racket_z_cpa: None,
        parachor: None,
        critical_viscosity: None,
        triple_point_temperature: None,
        heat_of_fusion: None,
        ideal_gas_enthalpy_of_formation: None,
        cp_a: None,
        cp_b: None,
        cp_c: None,
        cp_d: None,
        attractive_m: None,
    }
}
