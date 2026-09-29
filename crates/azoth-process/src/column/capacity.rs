//! `DistillationColumn`'s two capacity-limit families, as the class computes them.
//!
//! **Two families, and they are not one.** `getFsFactor` is `u·sqrt(rho_g)` over the **total**
//! cross-section `pi·D²/4`; `getGasLoadFactor` is the Souders-Brown `Ks = u·sqrt(rho_g/(rho_l -
//! rho_g))` over the same area. Both were measured before either was ported, in
//! `validation/neqsim/captures/process_column_capacity.tsv`, and three of that capture's readings
//! shape this module:
//!
//! **The two read different densities.** The Fs family takes the *system's* `getDensity("kg/m3")`
//! and the K family `getPhase(0)`'s, so an over-cooled overhead with two phases has them a per
//! cent apart - `15.34351809893286` against `14.4273646532561` on the capture's two-phase row.
//!
//! **The K family's liquid density is usually not the liquid's.** `getPhase(0)` of NeqSim's liquid
//! outlet is the *gas* that outlet carries, gas-first ordering being the class's own, so
//! `rho_l - rho_g` falls under its `10.0` floor and `1000.0` is substituted. On the capture's
//! solved absorber that is the ordinary path: the fallback's K is `0.006694885232955643` against
//! the real liquid's `0.00826846105398939`.
//!
//! **Both velocities are superficial over the total area**, which is what separates them from
//! `TrayHydraulicsCalculator.getFsFactor()`, whose area is the net one less a downcomer. The two
//! ids are not interchangeable and the tray model refuses to be wired in as if they were.

use azoth_core::Result;
use azoth_core::units::{Length, meters};

use crate::segment::phase::{Pick, phase_view, system_mass_density, system_volumetric_flow};
use crate::stream::Stream;

/// `DistillationColumn`'s own default for `maxAllowableFsFactor`, from its constructor.
pub const DEFAULT_MAX_ALLOWABLE_FS_FACTOR: f64 = 2.5;
/// `AbsorptionColumn`'s, set in its constructor - so the absorber pair reads `3.0` where the
/// distillation column and `PackedColumn` read `2.5`.
pub const DEFAULT_MAX_ALLOWABLE_FS_FACTOR_ABSORBER: f64 = 3.0;
/// `AbsorptionColumn.DEFAULT_MAX_ALLOWABLE_GAS_LOAD_FACTOR`.
pub const DEFAULT_MAX_ALLOWABLE_GAS_LOAD_FACTOR: f64 = 0.15;
/// `AbsorptionColumn.MIN_LIQUID_GAS_DENSITY_DIFFERENCE`: the floor a liquid density has to clear
/// above the gas's before it is used at all.
pub const MIN_LIQUID_GAS_DENSITY_DIFFERENCE: f64 = 10.0;
/// `AbsorptionColumn.DEFAULT_LIQUID_DENSITY`: what that floor substitutes.
pub const DEFAULT_LIQUID_DENSITY: f64 = 1000.0;

/// The two families at one solved state, as the class's four getters per family answer them.
///
/// **The gas-load half is `Option` because only the absorber pair has it.** `getGasLoadFactor` is
/// `AbsorptionColumn`'s own method, not the base's, so a distillation or packed column publishes
/// the Fs family alone - and a `None` here is that absence rather than a refused reading.
#[derive(Debug, Clone, PartialEq)]
pub struct CapacityLimits {
    /// `getFsFactor`: `u·sqrt(rho_g)` over the total cross-section.
    pub fs_factor: f64,
    /// `getFsFactorUtilization`: the factor over the limit, or zero where the limit is not
    /// positive.
    pub fs_factor_utilization: f64,
    /// `isFsFactorWithinDesignLimit`: `Fs <= limit`.
    pub fs_factor_within_design_limit: bool,
    /// `getMinimumDiameterForFsLimit`: the diameter at which the factor would equal the limit.
    ///
    /// **A function of the flow, the density and the limit alone**, so it is invariant under the
    /// diameter the factor is read at - the capture's three binary rows share one value.
    pub minimum_diameter_for_fs_limit: Length,
    /// `getGasLoadFactor`, on the absorber pair only.
    pub gas_load_factor: Option<f64>,
    /// `getGasLoadFactorUtilization`, on the absorber pair only.
    pub gas_load_factor_utilization: Option<f64>,
    /// `isGasLoadFactorWithinDesignLimit`, on the absorber pair only.
    pub gas_load_factor_within_design_limit: Option<bool>,
    /// `getMinimumDiameterForGasLoadLimit`, on the absorber pair only.
    pub minimum_diameter_for_gas_load_limit: Option<Length>,
}

/// Both families at the products of a solved column.
///
/// `gas_out` and `liquid_out` are the class's `getGasOutStream()` and `getLiquidOutStream()`, and
/// `gas_load_limit` is `Some` exactly where the class has the method - the absorber pair.
///
/// # Errors
/// Whatever the two outlets' flashes, their label rules or their densities raise. **No guard
/// here is a refusal the class makes**: its own zero answers (a non-positive area, a non-positive
/// limit, a missing outlet) are reproduced as zeros.
pub fn capacity_limits(
    gas_out: &Stream,
    liquid_out: &Stream,
    internal_diameter: Length,
    max_allowable_fs_factor: f64,
    gas_load_limit: Option<f64>,
) -> Result<CapacityLimits> {
    let diameter = internal_diameter.value;
    let area = std::f64::consts::PI * diameter * diameter / 4.0;

    // `getFsFactor`: the area guard comes first and answers zero.
    let fs_factor = if area <= 0.0 {
        0.0
    } else {
        system_volumetric_flow(gas_out)? / area * system_mass_density(gas_out)?.sqrt()
    };
    // `getMinimumDiameterForFsLimit` **does not consult the area at all** - its own guard is the
    // limit's - so a column with no diameter still owes this number.
    let minimum_diameter_for_fs_limit = if max_allowable_fs_factor > 0.0 {
        let flow = system_volumetric_flow(gas_out)?;
        let density = system_mass_density(gas_out)?;
        meters(
            (4.0 * flow * density.sqrt() / (std::f64::consts::PI * max_allowable_fs_factor)).sqrt(),
        )
    } else {
        meters(0.0)
    };

    let mut limits = CapacityLimits {
        fs_factor,
        fs_factor_utilization: utilization(fs_factor, max_allowable_fs_factor),
        fs_factor_within_design_limit: fs_factor <= max_allowable_fs_factor,
        minimum_diameter_for_fs_limit,
        gas_load_factor: None,
        gas_load_factor_utilization: None,
        gas_load_factor_within_design_limit: None,
        minimum_diameter_for_gas_load_limit: None,
    };

    if let Some(limit) = gas_load_limit {
        // `getGasSuperficialVelocity`: the same total area, and the same guard.
        let velocity = if area <= 0.0 {
            0.0
        } else {
            system_volumetric_flow(gas_out)? / area
        };
        let gas_load_factor = if not_positive(velocity) {
            0.0
        } else {
            // **`getPhase(0)` of both outlets**, which is what `Pick::Gas` spells: the gas phase
            // where there is one and the single phase otherwise.
            let gas_density = phase_view(gas_out, Pick::Gas)?.density;
            let liquid_density =
                resolve_liquid_density(gas_density, phase_view(liquid_out, Pick::Gas)?.density);
            if not_positive(gas_density) || liquid_density.is_nan() {
                0.0
            } else {
                velocity * (gas_density / (liquid_density - gas_density)).sqrt()
            }
        };
        let minimum = if limit > 0.0 {
            let gas_density = phase_view(gas_out, Pick::Gas)?.density;
            let liquid_density =
                resolve_liquid_density(gas_density, phase_view(liquid_out, Pick::Gas)?.density);
            if not_positive(gas_density) || liquid_density.is_nan() {
                meters(0.0)
            } else {
                let permissible = limit * ((liquid_density - gas_density) / gas_density).sqrt();
                if permissible > 0.0 {
                    meters(
                        (4.0 * system_volumetric_flow(gas_out)?
                            / (std::f64::consts::PI * permissible))
                            .sqrt(),
                    )
                } else {
                    meters(0.0)
                }
            }
        } else {
            meters(0.0)
        };
        limits.gas_load_factor = Some(gas_load_factor);
        limits.gas_load_factor_utilization = Some(utilization(gas_load_factor, limit));
        limits.gas_load_factor_within_design_limit = Some(gas_load_factor <= limit);
        limits.minimum_diameter_for_gas_load_limit = Some(minimum);
    }

    Ok(limits)
}

/// `getFsFactorUtilization` and `getGasLoadFactorUtilization`, which share one shape: a limit that
/// is not positive answers zero rather than an infinity.
fn utilization(value: f64, limit: f64) -> f64 {
    if limit > 0.0 { value / limit } else { 0.0 }
}

/// The class's own `!(x > 0.0)` guard, which a `NaN` also fails.
///
/// **Spelled out rather than negated** because the negation is the point: `NaN > 0.0` is false,
/// so a `NaN` takes the zero branch exactly as a negative does - which is what makes the
/// class's guards total rather than partial.
fn not_positive(value: f64) -> bool {
    value.is_nan() || value <= 0.0
}

/// `resolveLiquidDensityForGasLoad`: the near-dry fallback, and its `NaN`.
///
/// **`NaN` propagates rather than substituting.** A field the difference cannot be taken against
/// leaves the comparison false in both branches, so the caller's `isNaN` test is what refuses it -
/// which is why this returns `NaN` rather than a zero the caller would have to distinguish.
#[must_use]
pub fn resolve_liquid_density(gas_density: f64, liquid_density_field: f64) -> f64 {
    let density = if liquid_density_field - gas_density < MIN_LIQUID_GAS_DENSITY_DIFFERENCE {
        DEFAULT_LIQUID_DENSITY
    } else {
        liquid_density_field
    };
    if density <= gas_density {
        f64::NAN
    } else {
        density
    }
}
