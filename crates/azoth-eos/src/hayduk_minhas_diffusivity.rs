//! `eos.hayduk_minhas_diffusivity` - the liquid binary diffusivity from the
//! Hayduk-Minhas correlation.
//!
//! Spec: `specs/calcs/eos/hayduk_minhas_diffusivity.toml`, which records the two
//! solvent correlations and the clamps.

use azoth_core::units::{
    DynamicViscosity, MolarVolume, ThermodynamicTemperature, square_meters_per_second,
};
use azoth_core::{Result, apply_checks};

use crate::results::HaydukMinhasDiffusivityResult;
use crate::spec_gen;

/// Which Hayduk-Minhas solvent correlation to apply.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum HaydukMinhasForm {
    /// The hydrocarbon-solvent correlation.
    Paraffin,
    /// The water-rich-solvent correlation.
    Aqueous,
}

impl HaydukMinhasForm {
    /// The spec's spelling.
    #[must_use]
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::Paraffin => "paraffin",
            Self::Aqueous => "aqueous",
        }
    }
}

impl std::str::FromStr for HaydukMinhasForm {
    type Err = String;

    fn from_str(s: &str) -> std::result::Result<Self, Self::Err> {
        match s {
            "paraffin" => Ok(Self::Paraffin),
            "aqueous" => Ok(Self::Aqueous),
            other => Err(format!(
                "unknown Hayduk-Minhas form `{other}`; expected `paraffin` or `aqueous`"
            )),
        }
    }
}

/// The binary diffusion coefficient at infinite dilution, from the Hayduk-Minhas
/// correlation.
///
/// `form` selects the solvent correlation; `VA` is the solute molar volume at its
/// normal boiling point and `eta` the solvent viscosity. `VA` is clamped to `[20,
/// 500]` cm**3/mol and `eta` to `[0.1, 100]` cP before the correlation is applied.
///
/// # Errors
/// * [`azoth_core::AzothError::OutOfRange`] if `T` is not positive.
///
/// # Example
/// ```
/// use azoth_core::units::{cubic_meters_per_mole, kelvins, pascal_seconds};
/// use azoth_eos::{HaydukMinhasForm, hayduk_minhas_diffusivity};
///
/// let r = hayduk_minhas_diffusivity(HaydukMinhasForm::Paraffin,
///     cubic_meters_per_mole(4.0203262233375156e-5), kelvins(298.15),
///     pascal_seconds(9.163064908813372e-4))?;
/// assert!((r.d.value - 4.391778089044218e-9).abs() < 1e-15);
/// # Ok::<(), azoth_core::AzothError>(())
/// ```
#[allow(non_snake_case)] // `VA` and `T` are the symbols in the published equation
pub fn hayduk_minhas_diffusivity(
    form: HaydukMinhasForm,
    VA: MolarVolume,
    T: ThermodynamicTemperature,
    eta: DynamicViscosity,
) -> Result<HaydukMinhasDiffusivityResult> {
    let spec = &spec_gen::HAYDUK_MINHAS_DIFFUSIVITY_SPEC;
    let mut warnings = Vec::new();

    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "VA" => Some(VA.value),
            "T" => Some(T.value),
            "eta" => Some(eta.value),
            _ => None,
        },
        &mut warnings,
    )?;

    let va_cm3 = (VA.value * 1.0e6).clamp(20.0, 500.0);
    let eta_cp = (eta.value * 1000.0).clamp(0.1, 100.0);

    let d_cm2s = match form {
        HaydukMinhasForm::Paraffin => {
            13.3e-8 * T.value.powf(1.47) * eta_cp.powf(10.2 / va_cm3 - 0.791) / va_cm3.powf(0.71)
        }
        HaydukMinhasForm::Aqueous => {
            let vi_term = (va_cm3.powf(-0.19) - 0.292).max(0.01);
            1.25e-8 * vi_term * T.value.powf(1.52) * eta_cp.powf(9.58 / va_cm3 - 1.12)
        }
    };
    let d = d_cm2s * 1.0e-4;

    apply_checks(
        spec.derived_checks(),
        |name| (name == "d").then_some(d),
        &mut warnings,
    )?;

    Ok(HaydukMinhasDiffusivityResult {
        d: square_meters_per_second(d),
        warnings,
    })
}
