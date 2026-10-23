//! `characterization.tbp_fraction_component` - what `addTBPfraction` makes of a cut.
//!
//! Spec: `specs/models/characterization/tbp_fraction_component.toml`. Oracle:
//! `validation/neqsim/captures/tbp_fraction_component_probe.tsv`.
//!
//! # A cut is two numbers and nineteen outputs
//!
//! [`crate::tbp_cut_properties`] answers what the TBP model computes from a molar mass and a
//! gravity. `SystemInterface.addTBPfraction` goes further: it hands those to a component and then
//! adds the numbers a cubic's density, viscosity, wax and heat-capacity models read - the
//! critical volume, the Peneloux Rackett compressibility, the parachor, the critical viscosity,
//! the fusion and formation data, the ideal-gas heat-capacity polynomial, and the Lennard-Jones
//! pair. This id is that composition, for the class's **default** model.
//!
//! # Two acentric factors, and they are not the same number
//!
//! The model's Edmister form and the fitted `m` do not imply each other. `addTBPfraction` sets
//! the model's value, then `AttractiveTermSrk.setm` overwrites it with the root of
//! `-0.176 w^2 + 1.574 w + (0.48 - m)` - so the component carries the *root*, and the critical
//! volume, which recomputes the model's own, is built from the Edmister value. Both are used
//! here, each where the class uses it.
//!
//! # The Rackett compressibility is a flash
//!
//! It is not a function of the cut's two numbers. `calcRacketZ` reads the molar volume of a
//! **one-component SRK phase** at 15 C and one atmosphere, built from the cut's own constants
//! and fitted `m`, and the class turns that volume and the Peneloux coefficient into
//! `Z_RA = 0.29441 - c / (0.40768 R Tc/Pc)`. The binary here is reproduced with
//! [`azoth_eos::pt_flash`], whose liquid root is the same cubic root the class solves.
//! NeqSim's gas constant is `8.3144621`, and the formula is a difference of two volumes, so the
//! port carries that constant rather than the library's CODATA one.

use azoth_core::units::{
    MassDensity, MolarMass, ThermodynamicTemperature, angstroms, cubic_meters_per_mole,
    joules_per_mole, kelvins, pascal_seconds, pascals,
};
use azoth_core::{AzothError, Result, WarningCode, apply_checks};
use azoth_eos::{Alpha, Component, Cubic, Mixture, pt_flash, srk_kappa};

use crate::combiner::implied_acentric_factor;
use crate::model_gen;
use crate::results::TbpFractionComponentResult;
use crate::tbp_cut_properties::{TbpModel, tbp_cut_properties};

/// NeqSim's `ThermodynamicConstantsInterface.R`, in J/(mol K). **Not the library's CODATA
/// constant**: the Rackett formula divides a difference of molar volumes by `R Tc/Pc`, and at
/// `1e-7` relative a CODATA constant is still wrong in the twelfth digit.
const NEQSIM_R: f64 = 8.314_462_1;

/// The reference system's temperature, `273.15 + 15`, and its pressure in the class's bar.
const REFERENCE_TEMPERATURE: f64 = 288.15;
/// See [`REFERENCE_TEMPERATURE`].
const REFERENCE_PRESSURE_BAR: f64 = 1.01325;

/// `TPSRKModel`'s `TPBracketcoefs`: the Rackett intercept and the Peneloux coefficient.
const RACKETT_INTERCEPT: f64 = 0.29441;
/// See [`RACKETT_INTERCEPT`].
const PENELOUX_COEFFICIENT: f64 = 0.40768;

/// The component `addTBPfraction` builds from a cut, by the class's default TBP model.
///
/// # Errors
/// * [`azoth_core::AzothError::OutOfRange`] if the molar mass or density is not positive.
/// * [`azoth_core::AzothError::SolverNotConverged`] or `OutOfRange` from the reference flash if
///   the cut has no liquid root at 15 C and one atmosphere, which a cut below its own critical
///   temperature always has.
///
/// # Example
/// ```
/// use azoth_core::units::{kilograms_per_cubic_meter, kilograms_per_mole};
/// use azoth_characterization::tbp_fraction_component;
///
/// let cut = tbp_fraction_component(
///     kilograms_per_mole(0.095),
///     kilograms_per_cubic_meter(730.0),
/// )?;
/// assert!((cut.boiling_temperature.value - 369.89575).abs() < 1e-6);
/// assert!((cut.critical_temperature.value - 532.4869590128519).abs() < 1e-9);
/// # Ok::<(), azoth_core::AzothError>(())
/// ```
pub fn tbp_fraction_component(
    molar_mass: MolarMass,
    density: MassDensity,
) -> Result<TbpFractionComponentResult> {
    let spec = &model_gen::TBP_FRACTION_COMPONENT_SPEC;
    let mut warnings = Vec::new();
    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "molar_mass" => Some(molar_mass.value),
            "density" => Some(density.value),
            _ => None,
        },
        &mut warnings,
    )?;
    let cut = tbp_cut_properties(Some(TbpModel::PedersenSrk), molar_mass, density, None)?;
    // The inner model's own warnings are the cut's state; a skipped optional bound is its
    // bookkeeping and not this id's.
    warnings.extend(
        cut.warnings
            .iter()
            .filter(|warning| warning.code != WarningCode::RangeCheckSkipped)
            .cloned(),
    );
    let Some(attraction_exponent) = cut.attraction_exponent else {
        return Err(AzothError::property_unavailable(
            "m".to_string(),
            "attraction exponent".to_string(),
            "the default TBP model always states one, and this path is only that model".to_string(),
        ));
    };
    let molar_mass_grams = molar_mass.value * 1000.0;
    let gravity = density.value / 1000.0;
    let tc = cut.tc;
    let pc = cut.pc;
    let pc_bar = pc.value / 1.0e5;
    // The component's own acentric factor is the root the fitted `m` implies, clamped by the
    // solver; `init` then re-derives the attractive term's coefficient from that factor, so the
    // component carries Soave's form on the root - not the fitted `m`, which for a heavy cut is
    // nowhere near it.
    let acentric = implied_acentric_factor(attraction_exponent);
    let attraction_exponent = srk_kappa(acentric)?.kappa;
    // The critical volume is built from the model's Edmister value, not from the root.
    let model_acentric = cut.acentric_factor;

    let racket_z = racket_z(
        molar_mass_grams,
        gravity,
        tc,
        pc_bar,
        attraction_exponent,
        acentric,
    )?;

    let mut critical_volume =
        (0.2918 - 0.0928 * model_acentric) * NEQSIM_R * tc.value / pc_bar * 10.0;
    if critical_volume < 0.0 {
        critical_volume = (0.2918 - 0.0928) * NEQSIM_R * tc.value / pc_bar * 10.0;
    }
    let watson = (cut.boiling_temperature.value * 1.8).cbrt() / gravity;
    let correction = ((12.8 - watson) * (10.0 - watson) / (10.0 * acentric)).powi(2);
    let cp_a = (-0.33886 + 0.02827 * watson - 0.26105 * correction
        + 0.59332 * acentric * correction)
        * 4.18682
        * molar_mass.value
        * 1.0e3;
    let cp_b = (-(0.9291 - 1.1543 * watson + 0.0368 * watson * watson) * 1.0e-4
        + correction * (4.56 - 9.48 * acentric) * 1.0e-4)
        * 4.18682
        * molar_mass.value
        * 1.8
        * 1.0e3;
    let cp_c = (-1.6658e-7 + correction * (0.536 - 0.6828 * acentric) * 1.0e-7)
        * 4.18682
        * molar_mass.value
        * 1.8
        * 1.8
        * 1.0e3;
    let triple_point = 374.5 + 0.02617 * molar_mass_grams - 20172.0 / molar_mass_grams;
    let pc_atmosphere = pc_bar / REFERENCE_PRESSURE_BAR;

    Ok(TbpFractionComponentResult {
        critical_temperature: tc,
        critical_pressure: pc,
        boiling_temperature: kelvins(cut.boiling_temperature.value),
        acentric_factor: acentric,
        attraction_exponent,
        critical_volume: cubic_meters_per_mole(critical_volume / 1.0e6),
        racket_z,
        racket_z_cpa: racket_z,
        parachor: 59.3 + 2.34 * molar_mass_grams,
        critical_viscosity: pascal_seconds(
            7.94830 * molar_mass_grams.sqrt() * pc_bar.powf(2.0 / 3.0) / tc.value.powf(1.0 / 6.0)
                * 1.0e-7,
        ),
        triple_point_temperature: kelvins(triple_point),
        heat_of_fusion: joules_per_mole(0.1426 / 0.238845 * molar_mass_grams * triple_point),
        ideal_gas_enthalpy_of_formation: joules_per_mole(-1462600.0 * molar_mass.value - 47566.0),
        cp_a,
        cp_b,
        cp_c,
        // NeqSim's TBP path never sets `CpE`; `cp_d`'s getter adds `1e-10`, which is the number a
        // component states and the number this row states.
        cp_d: 1.0e-10,
        lennard_jones_diameter: angstroms(
            (2.3551 - 0.087 * acentric) * (tc.value / pc_atmosphere).cbrt(),
        ),
        lennard_jones_energy: kelvins(tc.value * (0.7915 + 0.1693 * acentric * acentric)),
        warnings,
    })
}

/// `calcRacketZ`: the reference SRK phase's molar volume, turned into the Peneloux constant.
fn racket_z(
    molar_mass_grams: f64,
    gravity: f64,
    tc: ThermodynamicTemperature,
    pc_bar: f64,
    attraction_exponent: f64,
    acentric: f64,
) -> Result<f64> {
    let component = Component::new(tc, pascals(pc_bar * 1.0e5), acentric)?
        .with_alpha_params(vec![attraction_exponent]);
    let reference = Mixture::new(vec![component], vec![0.0])?
        .with_cubic(Cubic::Srk)
        .with_alpha(Alpha::SrkFitted);
    let flash = pt_flash(
        &reference,
        kelvins(REFERENCE_TEMPERATURE),
        pascals(REFERENCE_PRESSURE_BAR * 1.0e5),
        &[1.0],
    )?;
    // NeqSim's `getMolarVolume` is `Z R T / P` with `P` in bar and this `R`; the value is in the
    // same mixed unit as `M / (10 rho)` below, which is what makes the difference a real volume.
    let volume = flash.z_liquid * NEQSIM_R * REFERENCE_TEMPERATURE / REFERENCE_PRESSURE_BAR;
    let peneloux = volume - molar_mass_grams / (gravity * 10.0);
    Ok(RACKETT_INTERCEPT - peneloux / (PENELOUX_COEFFICIENT * NEQSIM_R * tc.value / pc_bar))
}
