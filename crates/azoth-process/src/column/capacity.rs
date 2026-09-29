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

/// `DistillationColumn`'s own constructor default for `internalDiameter`, which the whole family
/// inherits and none of the three overrides.
pub const DEFAULT_INTERNAL_DIAMETER_M: f64 = 1.0;
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

/// The Fs family at one solved state, as the class's four getters answer it.
#[derive(Debug, Clone, PartialEq)]
pub struct FsLimits {
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
}

/// The Souders-Brown family, which `AbsorptionColumn` has and the base class does not.
#[derive(Debug, Clone, PartialEq)]
pub struct GasLoadLimits {
    /// `getGasLoadFactor`.
    pub gas_load_factor: f64,
    /// `getGasLoadFactorUtilization`.
    pub gas_load_factor_utilization: f64,
    /// `isGasLoadFactorWithinDesignLimit`.
    pub gas_load_factor_within_design_limit: bool,
    /// `getMinimumDiameterForGasLoadLimit`.
    pub minimum_diameter_for_gas_load_limit: Length,
}

/// `getFsFactor` and its three siblings, at the products of a solved column.
///
/// **The Fs family is the base class's**, so every column has it - and `PackedColumn`'s diameter
/// is the one its own sizing resolved rather than the one a caller stated.
///
/// # Errors
/// Whatever the gas outlet's flash, its label rule or its density raises. **No guard here is a
/// refusal the class makes**: its own zero answers (a non-positive area, a non-positive limit) are
/// reproduced as zeros.
pub fn fs_limits(
    gas_out: &Stream,
    internal_diameter: Length,
    max_allowable_fs_factor: f64,
) -> Result<FsLimits> {
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

    Ok(FsLimits {
        fs_factor,
        fs_factor_utilization: utilization(fs_factor, max_allowable_fs_factor),
        fs_factor_within_design_limit: fs_factor <= max_allowable_fs_factor,
        minimum_diameter_for_fs_limit,
    })
}

/// `getGasLoadFactor` and its three siblings, on the absorber pair.
///
/// **`getPhase(0)` of both outlets, and the second one is not the liquid.** `getLiquidOutStream`
/// answers a stream whose thermo system is the *tray's*, so its phase 0 - the array is gas-first -
/// is the **vapour that tray carries**, and that is what `rho_l` is taken from. So the class's
/// `MIN_LIQUID_GAS_DENSITY_DIFFERENCE` floor is met on every solved state and
/// `DEFAULT_LIQUID_DENSITY = 1000.0` is substituted: **the fallback is the ordinary path rather
/// than an edge case.**
///
/// **The difference is measurable rather than a formality.** On the capture's stripper the field
/// density is `9.221459471348782` against a gas of `11.267630874486859`, so `rho_l - rho_g` is
/// negative and the `10.0` floor is what answers - while the *liquid's* own density on that state
/// would give `0.0011870273688000244` against the class's `0.0009402581113639056`, a factor of
/// `1.26`. `liquid_out_vapour` is therefore the tray's vapour and not the liquid product: a port
/// that read the product would be 26 per cent out on this row and exact on the absorber's, which
/// is why one row cannot settle it and two can.
///
/// # Errors
/// Whatever the two flashes, their label rules or their densities raise.
pub fn gas_load_limits(
    gas_out: &Stream,
    liquid_out_vapour: &Stream,
    internal_diameter: Length,
    max_allowable_gas_load_factor: f64,
) -> Result<GasLoadLimits> {
    let diameter = internal_diameter.value;
    let area = std::f64::consts::PI * diameter * diameter / 4.0;

    // `getGasSuperficialVelocity`: the same total area, and the same guard as the Fs family's.
    let velocity = if area <= 0.0 {
        0.0
    } else {
        system_volumetric_flow(gas_out)? / area
    };

    let gas_load_factor = if not_positive(velocity) {
        0.0
    } else {
        let (gas_density, liquid_density) = densities(gas_out, liquid_out_vapour)?;
        if not_positive(gas_density) || liquid_density.is_nan() {
            0.0
        } else {
            velocity * (gas_density / (liquid_density - gas_density)).sqrt()
        }
    };

    let minimum_diameter_for_gas_load_limit = if max_allowable_gas_load_factor > 0.0 {
        let (gas_density, liquid_density) = densities(gas_out, liquid_out_vapour)?;
        if not_positive(gas_density) || liquid_density.is_nan() {
            meters(0.0)
        } else {
            let permissible = max_allowable_gas_load_factor
                * ((liquid_density - gas_density) / gas_density).sqrt();
            if permissible > 0.0 {
                meters(
                    (4.0 * system_volumetric_flow(gas_out)? / (std::f64::consts::PI * permissible))
                        .sqrt(),
                )
            } else {
                meters(0.0)
            }
        }
    } else {
        meters(0.0)
    };

    Ok(GasLoadLimits {
        gas_load_factor,
        gas_load_factor_utilization: utilization(gas_load_factor, max_allowable_gas_load_factor),
        gas_load_factor_within_design_limit: gas_load_factor <= max_allowable_gas_load_factor,
        minimum_diameter_for_gas_load_limit,
    })
}

/// The two densities both gas-load getters read: `getPhase(0)` of the gas outlet and of the
/// liquid outlet's own system, with the near-dry fallback applied to the second.
fn densities(gas_out: &Stream, liquid_out_vapour: &Stream) -> Result<(f64, f64)> {
    let gas_density = phase_view(gas_out, Pick::Gas)?.density;
    let liquid_density = resolve_liquid_density(
        gas_density,
        phase_view(liquid_out_vapour, Pick::Gas)?.density,
    );
    Ok((gas_density, liquid_density))
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
