//! `eos.tbp_fraction_properties` - a TBP pseudo-component's properties from two numbers.
//!
//! ```text
//! tc, pc, tb, m = Pedersen's TBP correlations(molar_mass, density)
//! omega         = 3/7 log10(pc/1.01325) / (tc/tb - 1) - 1
//! ```
//!
//! Spec: `specs/calcs/eos/tbp_fraction_properties.toml`, which carries the coefficient sets,
//! where the two branches switch, and the measurement of what the heavy set does above
//! `mw = 1120` g/mol.
//!
//! NeqSim's `addTBPfraction` is this and nothing else: a plus fraction is split into cuts
//! described by a molar mass and a normal liquid density, and these correlations are what
//! turn a cut into a component a cubic can take. The units are the correlation's own inside -
//! **g/mol and g/cm3** - and the conversion is at the boundary, because that is where NeqSim
//! puts it (`molarMass = 1000 * molarMass`).

use azoth_core::units::{kelvins, pascals};
use azoth_core::{Result, apply_checks};

use crate::results::TbpFractionPropertiesResult;
use crate::spec_gen;

/// The molar mass above which the critical-property coefficients change set, in g/mol.
const HEAVY_CUT_MOLAR_MASS: f64 = 1120.0;

/// The molar mass above which the boiling point takes the power law, in g/mol.
const POWER_LAW_MOLAR_MASS: f64 = 540.0;

/// `1.01325` bar, the reference pressure the acentric factor is reduced against.
const REFERENCE_PRESSURE_BAR: f64 = 1.01325;

/// The coefficients, in NeqSim's own units: `[tc, pc, m]`, each `[c0, c1, c2, c3, c4]`.
///
/// `PedersenTBPModelSRK.TBPfractionCoefOil` and `TBPfractionCoefsHeavyOil`, in that order.
const OIL: [[f64; 5]; 3] = [
    [163.12, 86.052, 0.43475, -1877.4, 0.0],
    [-0.13408, 2.5019, 208.46, -3987.2, 1.0],
    [0.7431, 0.0048122, 0.0096707, -3.7184e-6, 0.0],
];

/// The heavy-oil set, in force at and above `mw = 1120` g/mol.
const HEAVY: [[f64; 5]; 3] = [
    [8.3063e2, 1.75228e1, 4.55911e-2, -1.13484e4, 0.0],
    [8.02988e-1, 1.78396, 1.56740e2, -6.96559e3, 0.25],
    [-4.7268e-2, 6.02931e-2, 1.21051, -5.76676e-3, 0.0],
];

/// A TBP cut's critical properties, boiling point, acentric factor and alpha exponent.
///
/// # Errors
/// * [`azoth_core::AzothError::OutOfRange`] if the molar mass or the density is not positive,
///   each of which the correlation divides by or raises to a fractional power.
pub fn tbp_fraction_properties(
    molar_mass: f64,
    density: f64,
) -> Result<TbpFractionPropertiesResult> {
    let spec = &spec_gen::TBP_FRACTION_PROPERTIES_SPEC;
    let mut warnings = Vec::new();
    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "molar_mass" => Some(molar_mass),
            "density" => Some(density),
            _ => None,
        },
        &mut warnings,
    )?;

    // NeqSim's own units, which is what the coefficients are fitted in.
    let mw = molar_mass * 1000.0;
    let d = density / 1000.0;
    let set = if mw < HEAVY_CUT_MOLAR_MASS {
        &OIL
    } else {
        &HEAVY
    };

    let tc = set[0][0] * d + set[0][1] * mw.ln() + set[0][2] * mw + set[0][3] / mw;
    let pc_bar = (0.01325
        + set[1][0]
        + set[1][1] * d.powf(set[1][4])
        + set[1][2] / mw
        + set[1][3] / mw.powi(2))
    .exp();
    let tb = if mw < POWER_LAW_MOLAR_MASS {
        2.0e-6 * mw.powi(3) - 0.0035 * mw.powi(2) + 2.4003 * mw + 171.74
    } else {
        97.58 * mw.powf(0.3323) * d.powf(0.04609)
    };
    let exponent = set[2][0] + set[2][1] * mw + set[2][2] * d + set[2][3] * mw.powi(2);
    let acentric = 3.0 / 7.0 * (pc_bar / REFERENCE_PRESSURE_BAR).log10() / (tc / tb - 1.0) - 1.0;

    apply_checks(
        spec.derived_checks(),
        |name| (name == "acentric_factor").then_some(acentric),
        &mut warnings,
    )?;

    Ok(TbpFractionPropertiesResult {
        tc: kelvins(tc),
        pc: pascals(pc_bar * 1.0e5),
        boiling_temperature: kelvins(tb),
        acentric_factor: acentric,
        attraction_exponent: exponent,
        warnings,
    })
}
