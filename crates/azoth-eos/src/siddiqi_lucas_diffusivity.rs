//! `eos.siddiqi_lucas_diffusivity` - the liquid binary diffusivity from the
//! Siddiqi-Lucas correlation.
//!
//! Spec: `specs/calcs/eos/siddiqi_lucas_diffusivity.toml`, which records the two
//! solvent correlations and the viscosity floor.

use azoth_core::units::{
    DynamicViscosity, MolarVolume, ThermodynamicTemperature, square_meters_per_second,
};
use azoth_core::{Result, apply_checks};

use crate::results::SiddiqiLucasDiffusivityResult;
use crate::spec_gen;

/// Which Siddiqi-Lucas solvent correlation to apply.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SiddiqiLucasForm {
    /// The water-rich-solvent correlation.
    Aqueous,
    /// The organic-solvent correlation.
    Organic,
}

impl SiddiqiLucasForm {
    /// The spec's spelling.
    #[must_use]
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::Aqueous => "aqueous",
            Self::Organic => "organic",
        }
    }
}

impl std::str::FromStr for SiddiqiLucasForm {
    type Err = String;

    fn from_str(s: &str) -> std::result::Result<Self, Self::Err> {
        match s {
            "aqueous" => Ok(Self::Aqueous),
            "organic" => Ok(Self::Organic),
            other => Err(format!(
                "unknown Siddiqi-Lucas form `{other}`; expected `aqueous` or `organic`"
            )),
        }
    }
}

/// The binary diffusion coefficient at infinite dilution, from the Siddiqi-Lucas
/// correlation.
///
/// `form` selects the solvent correlation; `VA` and `VB` are the solute and solvent
/// molar volumes and `eta` the solvent viscosity, floored at 0.01 cP.
///
/// # Errors
/// * [`azoth_core::AzothError::OutOfRange`] if `T` is not positive.
///
/// # Example
/// ```
/// use azoth_core::units::{cubic_meters_per_mole, kelvins, pascal_seconds};
/// use azoth_eos::{SiddiqiLucasForm, siddiqi_lucas_diffusivity};
///
/// let r = siddiqi_lucas_diffusivity(SiddiqiLucasForm::Organic,
///     cubic_meters_per_mole(4.0203262233375156e-5),
///     cubic_meters_per_mole(8.816478555304741e-5), kelvins(298.15),
///     pascal_seconds(9.163064908813372e-4))?;
/// assert!((r.d.value - 1.984475615508378e-9).abs() < 1e-15);
/// # Ok::<(), azoth_core::AzothError>(())
/// ```
#[allow(non_snake_case)] // `VA`, `VB` and `T` are the symbols in the published equation
pub fn siddiqi_lucas_diffusivity(
    form: SiddiqiLucasForm,
    VA: MolarVolume,
    VB: MolarVolume,
    T: ThermodynamicTemperature,
    eta: DynamicViscosity,
) -> Result<SiddiqiLucasDiffusivityResult> {
    let spec = &spec_gen::SIDDIQI_LUCAS_DIFFUSIVITY_SPEC;
    let mut warnings = Vec::new();

    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "VA" => Some(VA.value),
            "VB" => Some(VB.value),
            "T" => Some(T.value),
            "eta" => Some(eta.value),
            _ => None,
        },
        &mut warnings,
    )?;

    let va_cm3 = VA.value * 1.0e6;
    let vb_cm3 = VB.value * 1.0e6;
    let eta_cp = (eta.value * 1000.0).max(0.01);

    let d_cm2s = match form {
        SiddiqiLucasForm::Aqueous => 2.98e-7 * eta_cp.powf(-1.026) * va_cm3.powf(-0.5473) * T.value,
        SiddiqiLucasForm::Organic => {
            9.89e-8 * eta_cp.powf(-0.907) * va_cm3.powf(-0.45) * vb_cm3.powf(0.265) * T.value
        }
    };
    let d = d_cm2s * 1.0e-4;

    apply_checks(
        spec.derived_checks(),
        |name| (name == "d").then_some(d),
        &mut warnings,
    )?;

    Ok(SiddiqiLucasDiffusivityResult {
        d: square_meters_per_second(d),
        warnings,
    })
}
