//! `eos.pr_peneloux_shift` - the Peng-Robinson Peneloux volume-translation parameter.
//!
//! ```text
//! c = 0.50033*(0.25969 - Z_RA)*R*Tc/Pc,  Z_RA = 0.29056 - 0.08775*omega
//! ```
//!
//! Spec: `specs/calcs/eos/pr_peneloux_shift.toml`, which records the dead NeqSim
//! expression this is *not* and why the Rackett compressibility is the correlation
//! rather than a databank value.

use azoth_core::units::{Pressure, ThermodynamicTemperature, cubic_meters_per_mole};
use azoth_core::{Result, apply_checks};

use crate::pr_molar_volume::MOLAR_GAS_CONSTANT;
use crate::results::PrPenelouxShiftResult;
use crate::spec_gen;

/// The Peng-Robinson Peneloux volume-translation parameter for a pure component.
///
/// `omega`, `Tc` and `Pc` are the caller's. The shift is subtracted from the
/// untranslated molar volume: `v_corr = v - c`.
///
/// # Errors
/// * [`azoth_core::AzothError::OutOfRange`] if `Tc` or `Pc` is not positive.
///
/// # Example
/// ```
/// use azoth_core::units::{kelvins, pascals};
/// use azoth_eos::pr_peneloux_shift;
///
/// let r = pr_peneloux_shift(0.152, kelvins(369.83), pascals(4_248_000.0), None)?;
/// assert!((r.c.value - (-6.349504285100194e-06)).abs() < 1e-18);
///
/// // Water, whose databank row carries a Rackett compressibility - and a shift of the
/// // opposite sign to the one the fallback correlation would give it.
/// let r = pr_peneloux_shift(0.3443, kelvins(647.096), pascals(22_064_000.0), Some(0.235662374))?;
/// assert!((r.c.value - 2.929079216846814e-06).abs() < 1e-18);
/// # Ok::<(), azoth_core::AzothError>(())
/// ```
#[allow(non_snake_case)] // `Tc` and `Pc` are the symbols in the published equation
pub fn pr_peneloux_shift(
    omega: f64,
    Tc: ThermodynamicTemperature,
    Pc: Pressure,
    z_ra: Option<f64>,
) -> Result<PrPenelouxShiftResult> {
    let spec = &spec_gen::PR_PENELOUX_SHIFT_SPEC;
    let mut warnings = Vec::new();

    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "omega" => Some(omega),
            "Tc" => Some(Tc.value),
            "Pc" => Some(Pc.value),
            "z_ra" => z_ra,
            _ => None,
        },
        &mut warnings,
    )?;

    // NeqSim's own rule: the field when it carries one, the correlation where it does
    // not. `Some(0.0)` and `None` are the same statement, because zero *is* the table's
    // spelling of absence - `ComponentPR.getVolumeCorrection` tests the value, not the
    // presence of a column.
    let z_ra = match z_ra {
        Some(value) if value.abs() >= 1e-10 => value,
        _ => 0.29056 - 0.08775 * omega,
    };
    let c = 0.50033 * (0.25969 - z_ra) * MOLAR_GAS_CONSTANT * Tc.value / Pc.value;

    apply_checks(
        spec.derived_checks(),
        |name| (name == "c").then_some(c),
        &mut warnings,
    )?;

    Ok(PrPenelouxShiftResult {
        c: cubic_meters_per_mole(c),
        warnings,
    })
}
