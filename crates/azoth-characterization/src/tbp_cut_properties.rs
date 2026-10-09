//! `characterization.tbp_cut_properties` - a TBP cut's critical properties.
//!
//! Spec: `specs/calcs/characterization/tbp_cut_properties.toml`, which carries the citations.
//! The oracle is `validation/neqsim/captures/characterization_probe.tsv`, ten models over
//! twelve cuts.
//!
//! # The unit boundary
//!
//! NeqSim's correlations take a molar mass in **g/mol** and a specific gravity in **g/cm3**.
//! The declared inputs are kg/mol and kg/m3 - the SI a caller has - and this converts once, at
//! the top. `density` happens to be the same *number* as the specific gravity, so only the
//! molar mass actually moves; the density division is written out anyway because it is a unit
//! conversion and not a no-op a reader should have to prove.
//!
//! # Two upstream behaviours reproduced rather than repaired
//!
//! **`PedersenSRKHeavyOil` is a no-op.** `TBPfractionModel$PedersenTBPModelSRKHeavyOil`
//! re-declares `TBPfractionCoefOil` and `TBPfractionCoefsHeavyOil` as *its own* fields,
//! shadowing the parent's, and the parent's `calcTC` reads the parent's - which the subclass
//! never touches. So the model returns `PedersenSRK`'s numbers. Measured: at `M = 500` g/mol,
//! `d = 0.88`, `tc = 891.945256085426` from both. Where this port encodes it, it says so.
//!
//! **`TwuModel` computes a molar mass and discards it.** `calcTC`, `calcPC` and
//! `calcCriticalVolume` each open with `double MW = solveMW(TB);` and never reference `MW`
//! again. The call is kept here for that reason and its result dropped.
//!
//! # What is not here
//!
//! `calcRacketZ` (needs a flashed reference system) and `calcParachorParameter`,
//! `calcCriticalViscosity`, `calcCriticalVolume` (no unit NeqSim states anywhere; and
//! `addTBPfraction` never calls `calcCriticalVolume` - it sets `critVol`, computed elsewhere).
//! The spec's `assumptions` carry the same statement.

use azoth_core::units::{MassDensity, MolarMass, ThermodynamicTemperature, kelvins, pascals};
use azoth_core::{Result, apply_checks};
use std::str::FromStr;

use crate::results::TbpCutPropertiesResult;
use crate::spec_gen;

/// The reference pressure NeqSim's acentric-factor correlations divide by, **in bar** - the unit
/// `pc` is in when they run. `ThermodynamicConstantsInterface.referencePressure`.
const REFERENCE_PRESSURE_BAR: f64 = 1.01325;

/// Which of NeqSim's ten TBP models evaluates a cut.
///
/// `TBPfractionModel.getAvailableModels()` names them; the spelling each variant carries is the
/// spec's, which is what `FromStr` accepts.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum TbpModel {
    /// `PedersenTBPModelSRK` - the class's default, and the model the seam id
    /// `eos.tbp_fraction_properties` carries.
    #[default]
    PedersenSrk,
    /// `PedersenTBPModelSRKHeavyOil`. **A no-op upstream** - see the module documentation.
    PedersenSrkHeavyOil,
    /// `PedersenTBPModelPR`.
    PedersenPr,
    /// `PedersenTBPModelPR2`, which differs only in its boiling point.
    PedersenPr2,
    /// `PedersenTBPModelPRHeavyOil`, whose coefficient set is the heavy one on both branches.
    PedersenPrHeavyOil,
    /// `RiaziDaubert`, which falls back to `PedersenSRK` above 300 g/mol.
    RiaziDaubert,
    /// `LeeKesler`.
    LeeKesler,
    /// `TwuModel`.
    Twu,
    /// `CavettModel`.
    Cavett,
    /// `StandingModel`, which is `RiaziDaubert` without the switch.
    Standing,
}

impl FromStr for TbpModel {
    type Err = std::convert::Infallible;

    /// NeqSim's `getModel` defaults an unrecognised name to `PedersenSRK` rather than refusing,
    /// and this does the same: the spec's own vocabulary check is what should refuse.
    fn from_str(text: &str) -> std::result::Result<Self, Self::Err> {
        Ok(match text.trim() {
            "pedersen_srk_heavy_oil" => Self::PedersenSrkHeavyOil,
            "pedersen_pr" => Self::PedersenPr,
            "pedersen_pr2" => Self::PedersenPr2,
            "pedersen_pr_heavy_oil" => Self::PedersenPrHeavyOil,
            "riazi_daubert" => Self::RiaziDaubert,
            "lee_kesler" => Self::LeeKesler,
            "twu" => Self::Twu,
            "cavett" => Self::Cavett,
            "standing" => Self::Standing,
            _ => Self::PedersenSrk,
        })
    }
}

/// NeqSim's coefficient sets, in its own scale: `[c0, c1, c2, c3, c4]` per row, row 0 for `tc`,
/// row 1 for `pc` (whose `c4` is the density exponent) and row 2 for the alpha exponent `m`.
const SRK_OIL: [[f64; 5]; 3] = [
    [163.12, 86.052, 0.43475, -1877.4, 0.0],
    [-0.13408, 2.5019, 208.46, -3987.2, 1.0],
    [0.7431, 0.0048122, 0.0096707, -3.7184e-6, 0.0],
];
const SRK_HEAVY_OIL: [[f64; 5]; 3] = [
    [8.3063e2, 1.75228e1, 4.55911e-2, -1.13484e4, 0.0],
    [8.02988e-1, 1.78396, 1.56740e2, -6.96559e3, 0.25],
    [-4.7268e-2, 6.02931e-2, 1.21051, -5.76676e-3, 0.0],
];
const PR_OIL: [[f64; 5]; 3] = [
    [73.4043, 97.3562, 0.618744, -2059.32, 0.0],
    [0.0728462, 2.18811, 163.91, -4043.23, 0.25],
    [0.373765, 0.00549269, 0.0117934, -4.93049e-6, 0.0],
];
const PR_HEAVY_OIL: [[f64; 5]; 3] = [
    [9.13222e2, 1.01134e1, 4.54194e-2, -1.3587e4, 0.0],
    [1.28155, 1.26838, 1.67106e2, -8.10164e3, 0.25],
    [-2.3838e-1, 6.10147e-2, 1.32349, -6.52067e-3, 0.0],
];

/// The switch both Pedersen correlations make, at this molar mass in g/mol.
const PEDERSEN_HEAVY_SWITCH: f64 = 1120.0;
/// The boiling-point correlation's own switch, in g/mol.
const BOILING_POINT_SWITCH: f64 = 540.0;
/// `RiaziDaubert`'s fallback to its parent, in g/mol. The test is `>`, so 300 is not over it.
const RIAZI_DAUBERT_SWITCH: f64 = 300.0;

/// The Pedersen coefficient set a model uses at this molar mass.
fn pedersen_coefs(model: TbpModel, m: f64) -> &'static [[f64; 5]; 3] {
    let heavy = m >= PEDERSEN_HEAVY_SWITCH;
    match model {
        // `PedersenSRKHeavyOil` shadows its parent's fields and changes nothing, so it reads the
        // parent's set by the same switch.
        TbpModel::PedersenSrk | TbpModel::PedersenSrkHeavyOil => {
            if heavy { &SRK_HEAVY_OIL } else { &SRK_OIL }
        }
        // Its constructor assigns the heavy set to *both* fields, so both branches are heavy.
        TbpModel::PedersenPrHeavyOil => &PR_HEAVY_OIL,
        _ => {
            if heavy { &PR_HEAVY_OIL } else { &PR_OIL }
        }
    }
}

/// The Pedersen critical temperature, `[K]`.
fn pedersen_tc(coefs: &[[f64; 5]; 3], m: f64, d: f64) -> f64 {
    let c = coefs[0];
    c[0] * d + c[1] * m.ln() + c[2] * m + c[3] / m
}

/// The Pedersen critical pressure, **in bar**.
fn pedersen_pc(coefs: &[[f64; 5]; 3], m: f64, d: f64) -> f64 {
    let c = coefs[1];
    (0.01325 + c[0] + c[1] * d.powf(c[4]) + c[2] / m + c[3] / m.powi(2)).exp()
}

/// The Pedersen alpha exponent `m`.
fn pedersen_m(coefs: &[[f64; 5]; 3], m: f64, d: f64) -> f64 {
    let c = coefs[2];
    c[0] + c[1] * m + c[2] * d + c[3] * m.powi(2)
}

/// The boiling point the SRK family correlates, `[K]`.
fn srk_tb(m: f64, d: f64) -> f64 {
    if m < BOILING_POINT_SWITCH {
        2.0e-6 * m.powi(3) - 0.0035 * m.powi(2) + 2.4003 * m + 171.74
    } else {
        97.58 * m.powf(0.3323) * d.powf(0.04609)
    }
}

/// `TBPBaseModel.calcTB`'s fallback, `[K]`, which the three models that do not override it reach.
fn base_tb(m: f64, d: f64) -> f64 {
    (m / 5.805e-5 * d.powf(0.9371)).powf(1.0 / 2.3776)
}

/// `PedersenTBPModelPR2`'s boiling point: Soreide's correlation, returned in R and divided back
/// to K.
fn soreide_tb(m: f64, d: f64) -> f64 {
    (1928.3
        - 1.695e5
            * m.powf(-0.03522)
            * d.powf(3.266)
            * (-4.922e-3 * m - 4.7685 * d + 3.462e-3 * m * d).exp())
        / 1.8
}

/// Edmister's acentric factor from the critical pair.
fn acentric_edmister(tc: f64, tb: f64, pc_bar: f64) -> f64 {
    3.0 / 7.0 * (pc_bar / REFERENCE_PRESSURE_BAR).log10() / (tc / tb - 1.0) - 1.0
}

/// Kesler-Lee's acentric factor, two-branched on the reduced boiling point.
fn acentric_kesler_lee(tc: f64, tb: f64, pc_bar: f64, d: f64) -> f64 {
    let tbr = tb / tc;
    let pbr = REFERENCE_PRESSURE_BAR / pc_bar;
    if tbr < 0.8 {
        (pbr.ln() - 5.92714 + 6.09649 / tbr + 1.28862 * tbr.ln() - 0.169347 * tbr.powi(6))
            / (15.2518 - 15.6875 / tbr - 13.4721 * tbr.ln() + 0.43577 * tbr.powi(6))
    } else {
        let kw = tb.cbrt() / d;
        -7.904 + 0.1352 * kw - 0.007465 * kw * kw + 8.359 * tbr + (1.408 - 0.01063 * kw) / tbr
    }
}

/// `TwuModel.calculateTfunc`, the residual the n-alkane solve drives to zero.
fn twu_tfunc(mw: f64, tb: f64) -> f64 {
    let phi = mw.ln();
    (5.1264 + 2.71579 * phi - 0.28659 * phi * phi - 39.8544 / phi - 0.122488 / (phi * phi)).exp()
        - 13.7512 * phi
        + 19.6197 * phi * phi
        - tb
}

/// `TwuModel.solveMW`: a damped Newton with a central-difference gradient.
///
/// The loop condition is NeqSim's own, `abs(error) > 1e-6 && iter < 1000 || iter < 3`, which
/// binds three iterations at least and a thousand at most.
fn twu_solve_mw(tb: f64) -> f64 {
    let mut mw = tb / (5.8 - 0.0052 * tb);
    let mut iter = 0_u32;
    loop {
        iter += 1;
        let prev = mw;
        let delta = 1.0;
        let gradient = (twu_tfunc(mw + delta, tb) - twu_tfunc(mw - delta, tb)) / (2.0 * delta);
        mw -= 0.5 * twu_tfunc(mw, tb) / gradient;
        let error = (mw - prev).abs();
        if !((error > 1e-6 && iter < 1000) || iter < 3) {
            break;
        }
    }
    mw
}

/// `TwuModel`'s critical temperature and pressure, `(tc [K], pc [bar])`.
///
/// The n-alkane reference and the two shift functions are NeqSim's; `fV`/`VC` are computed for
/// `pc` and `tc` needs only `fT`.
fn twu_tc_pc(d: f64, tb: f64) -> (f64, f64) {
    // Computed by the class and never read - kept, with its result dropped, because that is what
    // the class does.
    let _discarded = twu_solve_mw(tb);
    let tc_n_alkane =
        tb / (0.533272 + 0.343831e-3 * tb + 2.526167e-7 * tb.powi(2) - 1.65848e-10 * tb.powi(3)
            + 4.60774e24 * tb.powf(-13.0));
    let phi = 1.0 - tb / tc_n_alkane;
    let sg_n_alkane =
        0.843593 - 0.128624 * phi - 3.36159 * phi.powi(3) - 13749.0 * phi.powf(12.0);
    let pc_n_alkane = (0.318317
        + 0.099334 * phi.sqrt()
        + 2.89698 * phi
        + 3.0054 * phi * phi
        + 8.65163 * phi.powi(4))
    .powi(2);
    let vc_n_alkane = (0.82055 + 0.715468 * phi + 2.21266 * phi.powi(3) + 13411.1 * phi.powf(14.0))
        .powf(-8.0);

    let delta_st = (5.0 * (sg_n_alkane - d)).exp() - 1.0;
    let f_t = delta_st
        * (-0.270159 * tb.powf(-0.5)
            + (0.0398285 - 0.706691 * tb.powf(-0.5)) * delta_st);
    let tc = tc_n_alkane * ((1.0 + 2.0 * f_t) / (1.0 - 2.0 * f_t)).powi(2);

    let delta_sp = (0.5 * (sg_n_alkane - d)).exp() - 1.0;
    let delta_sv = (4.0 * (sg_n_alkane * sg_n_alkane - d * d)).exp() - 1.0;
    let f_v = delta_sv
        * (0.347776 * tb.powf(-0.5)
            + (-0.182421 + 2.24890 * tb.powf(-0.5)) * delta_sv);
    let vc = vc_n_alkane * ((1.0 + 2.0 * f_v) / (1.0 - 2.0 * f_v)).powi(2);
    let f_p = delta_sp
        * ((2.53262 - 34.4321 * tb.powf(-0.5) - 0.00230193 * tb)
            + (-11.4277 + 187.934 * tb.powf(-0.5) + 0.00414963 * tb) * delta_sp);
    let pc = pc_n_alkane * (tc / tc_n_alkane) * (vc_n_alkane / vc)
        * ((1.0 + 2.0 * f_p) / (1.0 - 2.0 * f_p)).powi(2);
    // MPa to bar, the class's own conversion.
    (tc, pc * 10.0)
}

/// `LeeKesler`'s critical pair, `(tc [K], pc [bar])`, from its boiling point.
fn lee_kesler_tc_pc(tb: f64, d: f64) -> (f64, f64) {
    let tc = 189.8 + 450.6 * d + (0.4244 + 0.1174 * d) * tb + (0.1441 - 1.0069 * d) * 1e5 / tb;
    let log_pc = 3.3864
        - 0.0566 / d
        - ((0.43639 + 4.1216 / d + 0.21343 / (d * d)) * 1e-3 * tb)
        + ((0.47579 + 1.182 / d + 0.15302 / (d * d)) * 1e-6 * tb * tb)
        - ((2.4505 + 9.9099 / (d * d)) * 1e-10 * tb.powi(3));
    (tc, log_pc.exp() * 10.0)
}

/// `RiaziDaubert`'s own pair, `(tc [K], pc [bar])`, which `StandingModel` uses at every molar
/// mass and `RiaziDaubert` only up to 300 g/mol.
fn riazi_daubert_tc_pc(m: f64, d: f64) -> (f64, f64) {
    let tc = 5.0 / 9.0
        * 554.4
        * (-1.3478e-4 * m - 0.61641 * d).exp()
        * m.powf(0.2998)
        * d.powf(1.0555);
    let pc = 0.068947
        * 4.5203e4
        * (-1.8078e-3 * m + -0.3084 * d).exp()
        * m.powf(-0.8063)
        * d.powf(1.6015);
    (tc, pc)
}

/// `CavettModel`'s API gravity, `141.5/SG - 131.5`.
fn cavett_api(d: f64) -> f64 {
    141.5 / d - 131.5
}

/// A cut's six published properties, in NeqSim's own working units: `tc` and `tb` in K, `pc` in
/// bar, the exponent absent where the model has none.
struct CutProperties {
    tc: f64,
    pc_bar: f64,
    tb: f64,
    acentric_factor: f64,
    attraction_exponent: Option<f64>,
}

/// Evaluate one model over one cut, in g/mol and g/cm3.
fn evaluate(model: TbpModel, m: f64, d: f64, supplied_tb: Option<f64>) -> CutProperties {
    // Every model short-circuits on a supplied boiling point; only the correlation beneath it
    // differs.
    let correlated_tb =
        |f: fn(f64, f64) -> f64| supplied_tb.filter(|t| *t > 0.0).unwrap_or_else(|| f(m, d));

    match model {
        TbpModel::PedersenSrk | TbpModel::PedersenSrkHeavyOil => {
            let c = pedersen_coefs(model, m);
            let tc = pedersen_tc(c, m, d);
            let pc_bar = pedersen_pc(c, m, d);
            let tb = correlated_tb(srk_tb);
            CutProperties {
                tc,
                pc_bar,
                tb,
                acentric_factor: acentric_edmister(tc, tb, pc_bar),
                attraction_exponent: Some(pedersen_m(c, m, d)),
            }
        }
        TbpModel::PedersenPr | TbpModel::PedersenPr2 | TbpModel::PedersenPrHeavyOil => {
            let c = pedersen_coefs(model, m);
            let tc = pedersen_tc(c, m, d);
            let pc_bar = pedersen_pc(c, m, d);
            let tb = if model == TbpModel::PedersenPr2 {
                correlated_tb(soreide_tb)
            } else {
                correlated_tb(srk_tb)
            };
            CutProperties {
                tc,
                pc_bar,
                tb,
                acentric_factor: acentric_edmister(tc, tb, pc_bar),
                attraction_exponent: Some(pedersen_m(c, m, d)),
            }
        }
        TbpModel::RiaziDaubert => {
            let (tc, pc_bar) = if m > RIAZI_DAUBERT_SWITCH {
                let c = &SRK_OIL;
                (pedersen_tc(c, m, d), pedersen_pc(c, m, d))
            } else {
                riazi_daubert_tc_pc(m, d)
            };
            let tb = correlated_tb(|m, d| 97.58 * m.powf(0.3323) * d.powf(0.04609));
            CutProperties {
                tc,
                pc_bar,
                tb,
                acentric_factor: acentric_kesler_lee(tc, tb, pc_bar, d),
                attraction_exponent: None,
            }
        }
        TbpModel::LeeKesler => {
            let tb = correlated_tb(base_tb);
            let (tc, pc_bar) = lee_kesler_tc_pc(tb, d);
            CutProperties {
                tc,
                pc_bar,
                tb,
                acentric_factor: acentric_kesler_lee(tc, tb, pc_bar, d),
                attraction_exponent: None,
            }
        }
        TbpModel::Twu => {
            let tb = correlated_tb(base_tb);
            let (tc, pc_bar) = twu_tc_pc(d, tb);
            CutProperties {
                tc,
                pc_bar,
                tb,
                acentric_factor: acentric_edmister(tc, tb, pc_bar),
                attraction_exponent: None,
            }
        }
        TbpModel::Cavett => {
            let tb = correlated_tb(|m, d| 97.58 * m.powf(0.3323) * d.powf(0.04609));
            let (tc_base, pc_base) = lee_kesler_tc_pc(tb, d);
            let api = cavett_api(d);
            let (tc, pc_bar) = if api < 30.0 {
                (
                    tc_base * (1.0 + 0.002 * (30.0 - api)),
                    pc_base * (1.0 + 0.001 * (30.0 - api)),
                )
            } else {
                (tc_base, pc_base)
            };
            let acentric_factor = if tb / tc >= 1.0 {
                acentric_kesler_lee(tc, tb, pc_bar, d)
            } else {
                acentric_edmister(tc, tb, pc_bar).clamp(0.0, 1.5)
            };
            CutProperties {
                tc,
                pc_bar,
                tb,
                acentric_factor,
                attraction_exponent: None,
            }
        }
        TbpModel::Standing => {
            let (tc, pc_bar) = riazi_daubert_tc_pc(m, d);
            let tb = correlated_tb(base_tb);
            CutProperties {
                tc,
                pc_bar,
                tb,
                acentric_factor: acentric_kesler_lee(tc, tb, pc_bar, d),
                attraction_exponent: None,
            }
        }
    }
}

/// A TBP cut's critical properties, by any of NeqSim's ten models.
///
/// # Errors
/// [`azoth_core::AzothError::OutOfRange`] if `molar_mass` or `density` is not positive, or if a
/// supplied `boiling_point` is not positive.
///
/// # Example
/// ```
/// use azoth_core::units::{kilograms_per_cubic_meter, kilograms_per_mole};
/// use azoth_characterization::{TbpModel, tbp_cut_properties};
///
/// let r = tbp_cut_properties(
///     Some(TbpModel::LeeKesler),
///     kilograms_per_mole(0.5),
///     kilograms_per_cubic_meter(880.0),
///     None,
/// )?;
/// assert!((r.tc.value - 906.211218915845).abs() < 1e-9);
/// assert!(r.attraction_exponent.is_none());
/// # Ok::<(), azoth_core::AzothError>(())
/// ```
pub fn tbp_cut_properties(
    model: Option<TbpModel>,
    molar_mass: MolarMass,
    density: MassDensity,
    boiling_point: Option<ThermodynamicTemperature>,
) -> Result<TbpCutPropertiesResult> {
    let spec = &spec_gen::TBP_CUT_PROPERTIES_SPEC;
    let mut warnings = Vec::new();

    let supplied_tb = boiling_point.map(|t| t.value);
    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "molar_mass" => Some(molar_mass.value),
            "density" => Some(density.value),
            "boiling_point" => supplied_tb,
            // `model` is an enum and carries no bound; `None` is a skipped check by the same
            // rule as the absent boiling point.
            _ => None,
        },
        &mut warnings,
    )?;

    let model = model.unwrap_or_default();
    // The boundary: g/mol and g/cm3, which is what every correlation above is written in.
    let m = molar_mass.value * 1.0e3;
    let d = density.value / 1.0e3;

    let cut = evaluate(model, m, d, supplied_tb);

    let watson_k = (1.8 * cut.tb).cbrt() / d;

    apply_checks(
        spec.derived_checks(),
        |name| match name {
            "acentric_factor" => Some(cut.acentric_factor),
            "attraction_exponent" => cut.attraction_exponent,
            _ => None,
        },
        &mut warnings,
    )?;

    Ok(TbpCutPropertiesResult {
        tc: kelvins(cut.tc),
        pc: pascals(cut.pc_bar * 1.0e5),
        boiling_temperature: kelvins(cut.tb),
        acentric_factor: cut.acentric_factor,
        attraction_exponent: cut.attraction_exponent,
        watson_k,
        warnings,
    })
}
