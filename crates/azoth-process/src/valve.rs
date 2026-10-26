//! `ControlValveSizing`, the `ValveMechanicalDesign` default: the IEC 60534 simplified sizing.
//!
//! **The default sizing standard and not the IEC class.** `ValveMechanicalDesign`'s default
//! `valveSizingMethod` is `ControlValveSizing` (`setValveSizingStandard` selects the IEC
//! variants), and this module is that class: `calcKv`, its inverse flow and its outlet-pressure
//! solve, on the liquid and gas branches, with the opening characteristic and the fouling
//! fraction. **The other standards are refused by name**, in the model spec's `[[unported]]`
//! rows - each is a second sizing machine with its own corrections.
//!
//! **The gas formula is the class's, unit mixing included.** It reads `Q` at standard conditions
//! (273.15 K, 101.325 kPa) from the *actual* volumetric flow through `P1`, `T` and `Z`, then
//! `Kv = Q_std/(N9 P1 Y) sqrt(M T Z / x)`, with `N9 = 24.6`, `M` in g/mol, `T` in kelvin and
//! `P1` in kPa. The liquid form is `Kv = Q/(sqrt(dP/SG))` with `SG = rho/1000`. Both divide by
//! the opening factor afterwards, so the reported `Kv` is the rated value at 100% opening.
//!
//! **The outlet-pressure solve is a bisection and it is the class's own.** The gas branch is not
//! invertible in closed form, so the class bisects `P2` in `[0.1, P1 - 0.001]` bara against a
//! `Kv` tolerance of `1e-6`; the liquid branch is the analytic `P2 = P1 - (Q/Kv)^2 SG`.

use azoth_core::units::{Pressure, pascals};
use azoth_core::{AzothError, Result};
use azoth_eos::{Phase, pt_flash};

use crate::segment::phase::{
    liquid_view, system_mass_density, system_volumetric_flow, vapour_view,
};
use crate::stream::Stream;

/// `ControlValveSizing.KV_TO_CV_FACTOR`: `Cv = Kv * 1.156`.
pub const KV_TO_CV_FACTOR: f64 = 1.156;
/// `ControlValveSizing.N9`, the IEC 60534 constant for standard m3/h, kPa, K and g/mol.
pub const N9: f64 = 24.6;
/// `ControlValveSizing.T_STD`: the standard temperature the gas flow is quoted at, K.
pub const STANDARD_TEMPERATURE_K: f64 = 273.15;
/// `ControlValveSizing.P_STD_KPA`: the standard pressure the gas flow is quoted at, kPa.
pub const STANDARD_PRESSURE_KPA: f64 = 101.325;
/// `ControlValveSizing.xT`, the pressure-drop ratio factor.
pub const DEFAULT_XT: f64 = 0.137;
/// The ratio `Cg = Cv * 1360`.
pub const CG_PER_CV: f64 = 1360.0;
/// The bisection's `Kv` tolerance and cap, the class's own.
pub const BISECTION_TOLERANCE: f64 = 1.0e-6;
/// See [`BISECTION_TOLERANCE`].
pub const MAX_BISECTION_ITERATIONS: usize = 100;
/// `ThermodynamicConstantsInterface.R`, J/(mol*K), the constant `getGamma2` subtracts.
pub const R_J_PER_MOL_K: f64 = 8.314;
/// `ThrottlingValve.percentValveOpening`'s initialiser.
pub const DEFAULT_OPENING_PERCENT: f64 = 100.0;
/// `ValveMechanicalDesign`'s default characteristic.
pub const DEFAULT_CHARACTERISTIC: ValveCharacteristic = ValveCharacteristic::Linear;
/// `EqualPercentageCharacteristic`'s default rangeability.
pub const DEFAULT_RANGEABILITY: f64 = 50.0;
/// `ModifiedParabolicCharacteristic`'s default exponent.
pub const DEFAULT_PARABOLIC_EXPONENT: f64 = 2.0;

/// The class's four `ValveCharacteristic` implementations, by the names its setter takes.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ValveCharacteristic {
    /// `LinearCharacteristic`: the opening itself.
    Linear,
    /// `EqualPercentageCharacteristic`: `R^(x-1)`, floored at `1/R` and capped at `1`.
    EqualPercentage,
    /// `QuickOpeningCharacteristic`: `sqrt(x)`.
    QuickOpening,
    /// `ModifiedParabolicCharacteristic`: `x^n`.
    ModifiedParabolic,
}

impl ValveCharacteristic {
    /// The characteristic a declaration's name selects.
    ///
    /// # Errors
    /// [`AzothError::InvalidInput`] on a name the class does not carry.
    pub fn named(name: &str) -> Result<Self> {
        match name {
            "linear" => Ok(Self::Linear),
            "equal_percentage" => Ok(Self::EqualPercentage),
            "quick_opening" => Ok(Self::QuickOpening),
            "modified_parabolic" => Ok(Self::ModifiedParabolic),
            other => Err(AzothError::invalid_input(
                "valve_characterization",
                format!(
                    "'{other}' is not one of linear, equal_percentage, quick_opening, \
                     modified_parabolic"
                ),
            )),
        }
    }

    /// `getOpeningFactor(percentOpening)`, the fraction of the rated `Kv` the opening gives.
    #[must_use]
    pub fn opening_factor(self, percent_opening: f64) -> f64 {
        let opening = percent_opening.clamp(0.0, 100.0);
        let x = opening / 100.0;
        match self {
            Self::Linear => x,
            Self::EqualPercentage => {
                if opening <= 0.0 {
                    1.0 / DEFAULT_RANGEABILITY
                } else if opening >= 100.0 {
                    1.0
                } else {
                    DEFAULT_RANGEABILITY.powf(x - 1.0)
                }
            }
            Self::QuickOpening => x.sqrt(),
            Self::ModifiedParabolic => x.powf(DEFAULT_PARABOLIC_EXPONENT),
        }
    }
}

/// The sizing's own configuration, as `ValveMechanicalDesign` and `ThrottlingValve` hold it.
#[derive(Debug, Clone, PartialEq)]
pub struct ValveSizingGeometry {
    /// The rated `Kv` at 100% opening, or `None` for the class's own sizing (`calcKv`).
    pub kv: Option<f64>,
    /// `percentValveOpening`, in percent.
    pub opening_percent: f64,
    /// `foulingFraction`: the effective `Kv` is `Kv * (1 - foulingFraction)`.
    pub fouling_fraction: f64,
    /// The opening characteristic.
    pub characteristic: ValveCharacteristic,
    /// `xT`.
    pub xt: f64,
    /// `allowChoked`: whether the sizing caps the flow at the choked limit.
    pub allow_choked: bool,
}

impl Default for ValveSizingGeometry {
    fn default() -> Self {
        Self {
            kv: None,
            opening_percent: DEFAULT_OPENING_PERCENT,
            fouling_fraction: 0.0,
            characteristic: DEFAULT_CHARACTERISTIC,
            xt: DEFAULT_XT,
            allow_choked: false,
        }
    }
}

/// What `calcValveSize` and its inverse solve report.
#[derive(Debug, Clone, PartialEq)]
pub struct ValveSizingReport {
    /// The rated `Kv` at 100% opening - the class's `getKv()`.
    pub kv: f64,
    /// `Cv = Kv * 1.156`.
    pub cv: f64,
    /// `Cg = Cv * 1360`.
    pub cg: f64,
    /// `x >= Fgamma * xT` on a gas valve.
    pub choked: bool,
    /// `isGasValve()`: a gas phase carrying more than half the volume.
    pub gas_valve: bool,
    /// `calculateMolarFlow()`: the mole flow the adjusted `Kv` and the drop give, mol/s.
    pub calculated_molar_flow: f64,
    /// `findOutletPressureForFixedKv`: the outlet pressure the trim gives, in bara.
    pub calculated_outlet_pressure: f64,
}

/// `isGasValve()`: a gas phase carrying more than half the volume.
///
/// # Errors
/// Whatever the inlet's flash and its two phase views raise.
pub fn gas_valve(inlet: &Stream) -> Result<bool> {
    let (mixture, _) = inlet.mixture()?;
    let flash = pt_flash(&mixture, inlet.t, inlet.p, &inlet.z)?;
    match flash.phase {
        Phase::AllVapour => Ok(true),
        Phase::AllLiquid | Phase::Trivial => Ok(false),
        _ => {
            let y = flash.vapour_fraction.unwrap_or(0.0);
            let gas = vapour_view(inlet)?;
            let liquid = liquid_view(inlet)?;
            let gas_volume = y * gas.molar_mass / gas.density;
            let liquid_volume = (1.0 - y) * liquid.molar_mass / liquid.density;
            Ok(gas_volume / (gas_volume + liquid_volume) > 0.5)
        }
    }
}

/// The inlet reads the sizing takes, in the class's own units.
struct Reads {
    p1_bara: f64,
    t_k: f64,
    vapor_flow_m3hr: f64,
    density_kg_m3: f64,
    z: f64,
    gamma: f64,
    molar_mass_gr_mol: f64,
}

fn reads(inlet: &Stream) -> Result<Reads> {
    let density = system_mass_density(inlet)?;
    let vapor_flow_m3s = system_volumetric_flow(inlet)?;
    let molar_mass = inlet.molar_mass()?.value;
    // The cubic's own Z, not the translated volume's: the class's `getZ()` is the EOS root.
    let z = {
        let cubic_density = inlet.density()?;
        inlet.p.value * molar_mass / (cubic_density * R_J_PER_MOL_K * inlet.t.value)
    };
    let cp = inlet.molar_heat_capacity()?;
    Ok(Reads {
        p1_bara: inlet.p.value / 1.0e5,
        t_k: inlet.t.value,
        vapor_flow_m3hr: vapor_flow_m3s * 3600.0,
        density_kg_m3: density,
        z,
        gamma: cp / (cp - R_J_PER_MOL_K),
        molar_mass_gr_mol: molar_mass * 1000.0,
    })
}

/// `getActualKv` then `adjustKv`: the characteristic's fraction of the rated `Kv`, less fouling.
#[must_use]
pub fn adjusted_kv(rated_kv: f64, geometry: &ValveSizingGeometry) -> f64 {
    rated_kv
        * geometry
            .characteristic
            .opening_factor(geometry.opening_percent)
        * (1.0 - geometry.fouling_fraction)
}

/// `calcKv`: the rated `Kv` the current flow and drop require.
///
/// # Errors
/// [`AzothError::InvalidInput`] on a non-positive drop or opening, and whatever the reads raise.
pub fn calc_kv(inlet: &Stream, outlet_bara: f64, geometry: &ValveSizingGeometry) -> Result<f64> {
    let read = reads(inlet)?;
    let flow = GasOrLiquid::classify(inlet)?;
    let delta_p = read.p1_bara - outlet_bara;
    if delta_p <= 0.0 {
        return Err(AzothError::invalid_input(
            "outlet_pressure",
            format!(
                "sizing needs a positive drop, and {outlet_bara} bara against {} bara is \
                 {delta_p} bar",
                read.p1_bara
            ),
        ));
    }
    let opening_factor = geometry
        .characteristic
        .opening_factor(geometry.opening_percent);
    if opening_factor <= 0.0 {
        return Err(AzothError::invalid_input(
            "valve_opening",
            format!(
                "an opening of {}% gives an opening factor of {opening_factor}, and `calcKv` \
                 divides by it",
                geometry.opening_percent
            ),
        ));
    }
    let kv = match flow {
        GasOrLiquid::Liquid => {
            let sg = read.density_kg_m3 / 1000.0;
            read.vapor_flow_m3hr / (delta_p / sg).sqrt()
        }
        GasOrLiquid::Gas => {
            let p1_kpa = read.p1_bara * 100.0;
            let x = delta_p / read.p1_bara;
            let f_gamma = read.gamma / 1.40;
            let choked = x >= f_gamma * geometry.xt;
            let x_effective = if choked && geometry.allow_choked {
                f_gamma * geometry.xt
            } else {
                x
            };
            let y = (1.0 - x_effective / (3.0 * f_gamma * geometry.xt)).max(2.0 / 3.0);
            let standard_flow = read.vapor_flow_m3hr
                * (p1_kpa / STANDARD_PRESSURE_KPA)
                * (STANDARD_TEMPERATURE_K / read.t_k)
                / read.z;
            let x_for_scale = if choked && geometry.allow_choked {
                geometry.xt * f_gamma
            } else {
                x
            };
            standard_flow / (N9 * p1_kpa * y)
                * (read.molar_mass_gr_mol * read.t_k * read.z / x_for_scale).sqrt()
        }
    };
    Ok(kv / opening_factor)
}

/// `calculateMolarFlow` (through `calculateFlowRateFromValveOpening`): the mole flow an adjusted
/// `Kv` and the drop give, mol/s.
///
/// # Errors
/// Whatever the reads raise. A non-positive drop answers zero, which is the class's own branch.
pub fn calculated_molar_flow(
    inlet: &Stream,
    outlet_bara: f64,
    adjusted: f64,
    geometry: &ValveSizingGeometry,
) -> Result<f64> {
    let read = reads(inlet)?;
    let delta_p = read.p1_bara - outlet_bara;
    if delta_p <= 0.0 {
        return Ok(0.0);
    }
    let flow_m3s = match GasOrLiquid::classify(inlet)? {
        GasOrLiquid::Liquid => {
            let sg = read.density_kg_m3 / 1000.0;
            adjusted * (delta_p / sg).sqrt() / 3600.0
        }
        GasOrLiquid::Gas => {
            let p1_kpa = read.p1_bara * 100.0;
            let x = delta_p / read.p1_bara;
            let f_gamma = read.gamma / 1.40;
            let choked = x >= f_gamma * geometry.xt;
            let x_effective = if choked && geometry.allow_choked {
                f_gamma * geometry.xt
            } else {
                x
            };
            let y = (1.0 - x_effective / (3.0 * f_gamma * geometry.xt)).max(2.0 / 3.0);
            let x_for_scale = if choked && geometry.allow_choked {
                geometry.xt * f_gamma
            } else {
                x
            };
            let standard_flow = adjusted * N9 * p1_kpa * y
                / (read.molar_mass_gr_mol * read.t_k * read.z / x_for_scale).sqrt();
            let actual = standard_flow
                * (STANDARD_PRESSURE_KPA / p1_kpa)
                * (read.t_k / STANDARD_TEMPERATURE_K)
                * read.z;
            actual / 3600.0
        }
    };
    Ok(flow_m3s * read.density_kg_m3 / (read.molar_mass_gr_mol / 1000.0))
}

/// `findOutletPressureForFixedKv`: the outlet pressure an adjusted `Kv` gives, in bara.
///
/// # Errors
/// [`AzothError::InvalidInput`] on a non-positive adjusted `Kv`, and whatever the reads raise.
pub fn calculated_outlet_pressure(
    inlet: &Stream,
    adjusted: f64,
    geometry: &ValveSizingGeometry,
) -> Result<f64> {
    if adjusted <= 0.0 {
        return Err(AzothError::invalid_input(
            "kv",
            format!("an adjusted Kv of {adjusted} gives no outlet pressure"),
        ));
    }
    let read = reads(inlet)?;
    if !matches!(GasOrLiquid::classify(inlet)?, GasOrLiquid::Gas) {
        let sg = read.density_kg_m3 / 1000.0;
        let q_m3hr = read.vapor_flow_m3hr;
        let delta_p_bar = (q_m3hr / adjusted).powi(2) * sg;
        return Ok(read.p1_bara - delta_p_bar);
    }
    let p1_kpa = read.p1_bara * 100.0;
    let standard_flow = read.vapor_flow_m3hr
        * (p1_kpa / STANDARD_PRESSURE_KPA)
        * (STANDARD_TEMPERATURE_K / read.t_k)
        / read.z;
    let f_gamma = read.gamma / 1.40;
    let mut low = 0.1_f64;
    let mut high = read.p1_bara - 0.001;
    let mut mid = (low + high) / 2.0;
    for _ in 0..MAX_BISECTION_ITERATIONS {
        mid = (low + high) / 2.0;
        let p2_kpa = mid * 100.0;
        let x = (p1_kpa - p2_kpa) / p1_kpa;
        let choked = x >= f_gamma * geometry.xt;
        let x_effective = if choked && geometry.allow_choked {
            f_gamma * geometry.xt
        } else {
            x
        };
        let y = (1.0 - x_effective / (3.0 * f_gamma * geometry.xt)).max(2.0 / 3.0);
        let x_for_scale = if choked && geometry.allow_choked {
            geometry.xt * f_gamma
        } else {
            x
        };
        let kv = standard_flow / (N9 * p1_kpa * y)
            * (read.molar_mass_gr_mol * read.t_k * read.z / x_for_scale).sqrt();
        if (kv - adjusted).abs() / adjusted < BISECTION_TOLERANCE {
            break;
        }
        if kv > adjusted {
            high = mid;
        } else {
            low = mid;
        }
    }
    Ok(mid)
}

enum GasOrLiquid {
    Gas,
    Liquid,
}

impl GasOrLiquid {
    fn classify(inlet: &Stream) -> Result<Self> {
        Ok(if gas_valve(inlet)? {
            Self::Gas
        } else {
            Self::Liquid
        })
    }
}

/// The whole report: the sized (or stated) `Kv`, the two derived coefficients, and the inverse
/// solves the model publishes.
///
/// # Errors
/// Whatever [`calc_kv`], [`calculated_molar_flow`] or [`calculated_outlet_pressure`] raise.
pub fn sizing_report(
    inlet: &Stream,
    outlet_bara: f64,
    geometry: &ValveSizingGeometry,
) -> Result<ValveSizingReport> {
    let gas = gas_valve(inlet)?;
    let read = reads(inlet)?;
    let kv = match geometry.kv {
        Some(kv) => kv,
        // **A non-positive drop has no sizing**, and the class's own `calcKv` answers `NaN`
        // there: a stated outlet above the inlet is the `acceptNegativeDP` row, whose stream
        // state is still a state. So the report carries the NaNs rather than refusing.
        None if read.p1_bara - outlet_bara > 0.0 => calc_kv(inlet, outlet_bara, geometry)?,
        None => f64::NAN,
    };
    if !kv.is_finite() {
        return Ok(ValveSizingReport {
            kv,
            cv: f64::NAN,
            cg: f64::NAN,
            choked: false,
            gas_valve: gas,
            calculated_molar_flow: 0.0,
            calculated_outlet_pressure: f64::NAN,
        });
    }
    let adjusted = adjusted_kv(kv, geometry);
    let choked = if gas {
        let x = (read.p1_bara - outlet_bara) / read.p1_bara;
        x >= (read.gamma / 1.40) * geometry.xt
    } else {
        false
    };
    Ok(ValveSizingReport {
        kv,
        cv: kv * KV_TO_CV_FACTOR,
        cg: kv * KV_TO_CV_FACTOR * CG_PER_CV,
        choked,
        gas_valve: gas,
        calculated_molar_flow: calculated_molar_flow(inlet, outlet_bara, adjusted, geometry)?,
        calculated_outlet_pressure: calculated_outlet_pressure(inlet, adjusted, geometry)?,
    })
}

/// The report's solved pressure as a [`Pressure`], for the model boundary.
#[must_use]
pub fn outlet_pressure(report: &ValveSizingReport) -> Pressure {
    pascals(report.calculated_outlet_pressure * 1.0e5)
}

#[cfg(test)]
mod tests {
    use super::*;
    use azoth_core::units::{kelvins, pascals};

    /// The capture's `sizing_gas_rich` state: 0.1/0.2/0.4/0.3 methane/ethane/propane/n-butane
    /// at 296.25 K and 1.82925 bara, 14942.5 kg/hr - 97.21961972278174 mol/s.
    fn rich_gas() -> Stream {
        Stream::from_pt(
            ["methane", "ethane", "propane", "n-butane"]
                .iter()
                .map(|n| (*n).to_string())
                .collect(),
            vec![0.1, 0.2, 0.4, 0.3],
            97.219_619_722_781_74,
            pascals(1.82925e5),
            kelvins(296.25),
        )
        .expect("the rich gas resolves")
    }

    fn methane() -> Stream {
        Stream::from_pt(
            vec!["methane".to_string()],
            vec![1.0],
            82.228_663_245_306_04,
            pascals(10.0e5),
            kelvins(293.15),
        )
        .expect("the methane resolves")
    }

    fn water() -> Stream {
        Stream::from_pt(
            vec!["water".to_string()],
            vec![1.0],
            26100.471_829_031_365,
            pascals(14.8e5),
            kelvins(289.15),
        )
        .expect("the water resolves")
    }

    /// **`1e-4` on the gas rows and `1e-9` on the liquid one.** The gas sizing reads `gamma`,
    /// and the two libraries' heat capacities differ by the divergence the plug-flow reactor's
    /// own test measures; the choked branch scales `Kv` as `1/sqrt(Fgamma)`, so that difference
    /// reaches the coefficient. The liquid form reads only the density, and it agrees to
    /// rounding.
    fn close(got: f64, want: f64, what: &str) {
        assert!(
            (got - want).abs() / want.abs() < 1e-9,
            "{what}: {got} against {want}"
        );
    }

    fn close_gas(got: f64, want: f64, what: &str) {
        assert!(
            (got - want).abs() / want.abs() < 1e-4,
            "{what}: {got} against {want}"
        );
    }

    fn close_loose(got: f64, want: f64, what: &str) {
        assert!(
            (got - want).abs() / want.abs() < 1e-5,
            "{what}: {got} against {want}"
        );
    }

    /// **`sizing_gas_rich`**: the default sizing on a gas valve, choked.
    #[test]
    fn the_rich_gas_row_is_reproduced() {
        let report = sizing_report(&rich_gas(), 1.41325, &ValveSizingGeometry::default())
            .expect("the sizing computes");
        assert!(report.gas_valve);
        assert!(report.choked);
        close_gas(report.kv, 607.954_033_813_071_4, "kv");
        close_gas(report.cv, 702.794_863_087_910_5, "cv");
        close_gas(report.cg, 955_801.013_799_558_3, "cg");
        close(report.calculated_molar_flow, 97.219_619_722_781_74, "flow");
        close_loose(
            report.calculated_outlet_pressure,
            1.413_249_538_421_630_6,
            "solved outlet",
        );
    }

    /// **`sizing_methane`**: a gas valve that is not choked.
    #[test]
    fn the_methane_row_is_reproduced() {
        let report = sizing_report(&methane(), 9.0, &ValveSizingGeometry::default())
            .expect("the sizing computes");
        assert!(report.gas_valve);
        assert!(!report.choked);
        close_gas(report.kv, 78.646_163_426_829_79, "kv");
        close_gas(report.cv, 90.914_964_921_415_23, "cv");
        close(report.calculated_molar_flow, 82.228_663_245_306_04, "flow");
        close_loose(
            report.calculated_outlet_pressure,
            8.999_995_656_967_162,
            "solved outlet",
        );
    }

    /// **`sizing_water_liquid`**: the liquid branch, analytic in both directions.
    #[test]
    fn the_water_row_is_reproduced() {
        let report = sizing_report(&water(), 4.46, &ValveSizingGeometry::default())
            .expect("the sizing computes");
        assert!(!report.gas_valve);
        assert!(!report.choked);
        // The liquid `Kv` reads the density, and the translated one carries the two
        // libraries' rounding at `3e-8`; the flow and the solved pressure are exact.
        assert!(
            (report.kv - 528.496_159_725_506_8).abs() / 528.496_159_725_506_8 < 1e-6,
            "kv: {}",
            report.kv
        );
        assert!(
            (report.cv - 610.941_560_642_685_8).abs() / 610.941_560_642_685_8 < 1e-6,
            "cv: {}",
            report.cv
        );
        close(report.calculated_molar_flow, 26100.471_829_031_365, "flow");
        close_loose(
            report.calculated_outlet_pressure,
            4.460_000_000_000_003,
            "solved outlet",
        );
    }

    /// **The opening and the fouling scale the adjusted `Kv`**, and the flow with it: a linear
    /// characteristic at half opening is half the flow, and a ten-per-cent fouling is nine
    /// tenths of it.
    #[test]
    fn the_opening_and_fouling_scale_the_flow() {
        let stated = ValveSizingGeometry {
            kv: Some(607.954_033_813_071_4),
            ..ValveSizingGeometry::default()
        };
        let full = sizing_report(&rich_gas(), 1.41325, &stated)
            .expect("the sizing computes")
            .calculated_molar_flow;
        let half = sizing_report(
            &rich_gas(),
            1.41325,
            &ValveSizingGeometry {
                opening_percent: 50.0,
                ..stated.clone()
            },
        )
        .expect("the sizing computes")
        .calculated_molar_flow;
        close(half / full, 0.5, "half opening");
        let fouled = sizing_report(
            &rich_gas(),
            1.41325,
            &ValveSizingGeometry {
                fouling_fraction: 0.1,
                ..stated
            },
        )
        .expect("the sizing computes")
        .calculated_molar_flow;
        close(fouled / full, 0.9, "ten-per-cent fouling");
    }

    /// **A non-positive drop answers `NaN`**, which is the class's own `calcKv` (it takes the
    /// square root of a negative ratio), and the `acceptNegativeDP` row's stream is still a
    /// state the model publishes.
    #[test]
    fn a_zero_drop_answers_no_sizing() {
        let report = sizing_report(&rich_gas(), 1.82925, &ValveSizingGeometry::default())
            .expect("the report computes");
        assert!(report.kv.is_nan());
        assert!(report.calculated_outlet_pressure.is_nan());
        assert_eq!(report.calculated_molar_flow, 0.0);
    }

    /// The four `getOpeningFactor` forms.
    #[test]
    fn the_characteristics_are_the_class_own() {
        assert!((ValveCharacteristic::Linear.opening_factor(50.0) - 0.5).abs() < 1e-15);
        assert!(
            (ValveCharacteristic::EqualPercentage.opening_factor(50.0) - 50.0_f64.powf(-0.5)).abs()
                < 1e-15
        );
        assert!(
            (ValveCharacteristic::QuickOpening.opening_factor(50.0) - 0.5_f64.sqrt()).abs() < 1e-15
        );
        assert!((ValveCharacteristic::ModifiedParabolic.opening_factor(50.0) - 0.25).abs() < 1e-15);
    }
}
