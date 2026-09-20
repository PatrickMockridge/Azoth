//! The three Helmholtz terms a Furst electrolyte phase adds to a cubic.
//!
//! `PhaseModifiedFurstElectrolyteEos.getF()` is
//!
//! ```text
//! F = F_SRK + FSR2 * sr2On + FLR * lrOn + FBorn * bornOn
//! ```
//!
//! and this module is the three terms with the derivatives the phase's `dFdT`, `dFdV`,
//! `dFdVdV` and `dFdTdT` read. The values are checked against
//! `validation/neqsim/captures/furst_probe.tsv` layer by layer; the state they are checked
//! at is `SystemFurstElectrolyteEosTest`'s own mixture.
//!
//! # The two state variables, and which volume
//!
//! Every term is a function of `(T, V, composition)`, and `V` is NeqSim's `getMolarVolume()`
//! **times `1e-5`** - the source writes that product everywhere, and it is the molar volume
//! in m³/mol that this module takes.
//!
//! **`dFdV` is the derivative with respect to the total volume, not the molar one.** The
//! distinction is invisible on a one-mole phase and decides whether the derivative is
//! intensive. Measured, by running the probe's mixture at ten times the moles: `FSR2`, `FLR`
//! and `FBorn` all scale by ten, while `dFdV` and the three terms' `dV` do not move at all.
//! `FSR2V` is written against `(V n)^2` and `epsdV` against `V n`, which is what makes that
//! come out right; reading either in isolation suggests the opposite, and the control is
//! what settles it.
//!
//! # NeqSim's own constants, again
//!
//! `electronCharge` is `1.6021917e-19` and `vacumPermittivity` is `8.85419e-12` - the values
//! of their era, not the current CODATA ones - and `R` is `8.3144621`. Each differs in the
//! seventh significant figure or later, and each lands in `alphaLR2`, which is the whole of
//! the long-range term's scale.

use azoth_core::{AzothError, Result};

use crate::furst_dielectric::{NEQSIM_AVOGADRO, NEQSIM_PI};

/// NeqSim's `electronCharge`, from `ThermodynamicConstantsInterface`.
pub const ELECTRON_CHARGE: f64 = 1.6021917e-19;

/// NeqSim's `vacumPermittivity`, from the same interface.
pub const VACUUM_PERMITTIVITY: f64 = 8.85419e-12;

/// NeqSim's `R`. The muzzle velocity of the constants above, and the reason every
/// electrostatic layer is checked against the probe rather than against a textbook.
pub const R: f64 = 8.3144621;

/// The state a term is evaluated at.
///
/// Everything the three terms read, assembled once per phase evaluation: the temperature,
/// the molar volume in m³/mol, the phase's total mole number, the packing fractions, the
/// solvent and phase dielectric constants, the shielding parameter and the two sums.
#[derive(Debug, Clone, PartialEq)]
pub struct FurstState {
    /// Absolute temperature, in K.
    pub temperature: f64,
    /// Molar volume, in m³/mol.
    pub molar_volume: f64,
    /// The phase's total mole number.
    pub moles: f64,
    /// The packing fraction over every component, `calcEps`.
    pub packing: f64,
    /// The packing fraction over the ions, `calcEpsIonic`.
    pub ionic_packing: f64,
    /// The solvent mixture's dielectric constant.
    pub solvent_dielectric: f64,
    /// Its temperature derivative.
    pub solvent_dielectric_dt: f64,
    /// The phase's dielectric constant, `1 + (eps_s - 1)(1 - eps_i)/(1 + eps_i/2)`.
    pub dielectric: f64,
    /// Its temperature derivative.
    pub dielectric_dt: f64,
    /// Its second temperature derivative.
    pub dielectric_dtdt: f64,
    /// Its volume derivative, `calcDiElectricConstantdV`.
    pub dielectric_dv: f64,
    /// Its second volume derivative.
    pub dielectric_dvdv: f64,
    /// Its mixed derivative.
    pub dielectric_dtdv: f64,
    /// The shielding parameter `gamma`, from the implicit solve.
    pub shielding: f64,
    /// Its temperature derivative, by implicit differentiation.
    pub shielding_dt: f64,
    /// `XLR = sum_i n_i z_i^2 gamma/(1 + gamma sigma_i)`, over the ions.
    pub xlr: f64,
    /// Its temperature derivative.
    pub xlr_dt: f64,
    /// `bornX = sum_i n_i z_i^2/sigma_i`, over the ions.
    pub born_x: f64,
    /// The short-range parameter `W = -sum_ij n_i n_j Wij(T)`.
    pub w: f64,
    /// `dW/dT`, NeqSim's `WT` - the **plain** derivative and not its `T`-scaled form,
    /// because that is what the capture prints and a second convention here would be one
    /// more thing to get wrong in a chain rule.
    pub w_dt: f64,
    /// `d^2W/dT^2`, NeqSim's `WTT`.
    pub w_dtdt: f64,
}

/// `sr2On`, `lrOn` and `bornOn`: which of the three terms the model runs.
///
/// NeqSim's `SystemFurstElectrolyteEos` leaves all three on, and nothing in `src/main`
/// turns one off, so this is a stated default rather than a caller's choice - the fields
/// exist because a model that could silently be a different model is worse than one that
/// says what it is.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct TermSwitches {
    /// The short-range `W` term.
    pub short_range: bool,
    /// The MSA long-range term.
    pub long_range: bool,
    /// The Born solvation term.
    pub born: bool,
}

impl TermSwitches {
    /// All three on, which is how both Fürst systems are constructed.
    #[must_use]
    pub fn all_on() -> Self {
        Self {
            short_range: true,
            long_range: true,
            born: true,
        }
    }
}

/// `alphaLR2 = e^2 N_A/(eps0 eps R T)`.
///
/// `eps` is the **phase's** dielectric constant and not the solvent's - the two differ by
/// the ionic correction, and the probe prints them as `78.1382247597158` against
/// `78.3148814553987`.
#[must_use]
pub fn alpha_lr2(dielectric: f64, temperature: f64) -> f64 {
    ELECTRON_CHARGE * ELECTRON_CHARGE * NEQSIM_AVOGADRO
        / (VACUUM_PERMITTIVITY * dielectric * R * temperature)
}

/// `dalphaLR2/dT`, NeqSim's `alphaLRdT`, plain.
///
/// `alphaLR2 = scale/(eps T)`, so the derivative is `-scale/(eps T^2) - scale eps'/(eps^2)`,
/// written with the `1/T^2` cancelled against the `T` in the scale to keep it in one place.
#[must_use]
pub fn alpha_lr2_dt(state: &FurstState) -> f64 {
    let (e, e0, na, t, eps) = (
        ELECTRON_CHARGE,
        VACUUM_PERMITTIVITY,
        NEQSIM_AVOGADRO,
        state.temperature,
        state.dielectric,
    );
    let scale = e * e * na / (e0 * R);
    // `alphaLR2 = scale/(eps T)`, so the derivative is
    // `-scale/(eps T^2) - scale eps'/(eps^2 T)` - the second term carries a `1/T` that is
    // easy to drop and that a central difference of `alphaLR2` catches immediately.
    -scale / (eps * t * t) - scale * state.dielectric_dt / (eps * eps * t)
}

/// `P dalphaLR2/dP` is not this; this is `dalphaLR2/dV`, NeqSim's `alphaLRdV`.
#[must_use]
pub fn alpha_lr2_dv(state: &FurstState) -> f64 {
    let (e, e0, na, t, eps) = (
        ELECTRON_CHARGE,
        VACUUM_PERMITTIVITY,
        NEQSIM_AVOGADRO,
        state.temperature,
        state.dielectric,
    );
    -e * e * na / (e0 * eps * eps * R * t) * state.dielectric_dv
}

/// `FBorn`'s prefactor `K(T) = N_A e^2/(4 pi eps0 R T)`.
///
/// The whole temperature dependence of the Born term's scale, and the reason its three
/// partials are all `K` times a power of `T`.
#[must_use]
pub fn born_scale(temperature: f64) -> f64 {
    NEQSIM_AVOGADRO * ELECTRON_CHARGE * ELECTRON_CHARGE
        / (4.0 * NEQSIM_PI * VACUUM_PERMITTIVITY * R * temperature)
}

/// The short-range term `FSR2 = W/(V n (1 - eps))`.
///
/// # Errors
/// * [`AzothError::OutOfRange`] if `1 - eps` is not positive. It is a volume fraction, so
///   `eps >= 1` is a phase whose components' excluded volumes exceed its own - not a state,
///   and NeqSim divides by it.
pub fn fsr2(state: &FurstState) -> Result<f64> {
    let denominator = state.molar_volume * state.moles * (1.0 - state.packing);
    if !denominator.is_finite() || denominator <= 0.0 {
        return Err(AzothError::out_of_range(
            "packing",
            state.packing,
            "the short-range term divides by `V n (1 - eps)` and the excluded volume is \
             larger than the phase's own, so the term has no value there",
        ));
    }
    Ok(state.w / denominator)
}

/// `T dFSR2/dT`, which is `FSR2W * WT` with `FSR2W = 1/(V n (1 - eps))`.
#[must_use]
pub fn t_d_fsr2_dt(state: &FurstState) -> f64 {
    state.w_dt * state.temperature / (state.molar_volume * state.moles * (1.0 - state.packing))
}

/// `dFSR2/dV`, against the total volume and in NeqSim's own scaling.
///
/// `(FSR2V + FSR2eps epsdV) * 1e-5`, each partial taken from the source rather than
/// re-derived - see the module note on which volume. The `1e-5` is NeqSim's molar-volume
/// scaling and is reproduced rather than converted away, because converting it would put
/// this quantity in different units from `dFdV`'s cubic part.
#[must_use]
pub fn fsr2_dv(state: &FurstState) -> f64 {
    let vn = state.molar_volume * state.moles;
    let one_minus = 1.0 - state.packing;
    let fsr2_v = -state.w / (vn * vn * one_minus);
    let fsr2_eps = state.w / (vn * one_minus * one_minus);
    let eps_dv = -state.packing / vn;
    (fsr2_v + fsr2_eps * eps_dv) * 1.0e-5
}

/// `d^2FSR2/dV^2`, against the total volume.
#[must_use]
pub fn fsr2_dvdv(state: &FurstState) -> f64 {
    let vn = state.molar_volume * state.moles;
    let one_minus = 1.0 - state.packing;
    let eps_dv = -state.packing / vn;
    let eps_dvdv = 2.0 * state.packing / (vn * vn);
    let fsr2_vv = 2.0 * state.w / (vn * vn * vn * one_minus);
    let fsr2_eps_v = -state.w / (vn * vn * one_minus * one_minus);
    let fsr2_eps_eps = 2.0 * state.w / (vn * one_minus * one_minus * one_minus);
    let fsr2_eps = state.w / (vn * one_minus * one_minus);
    (fsr2_vv + 2.0 * fsr2_eps_v * eps_dv + fsr2_eps_eps * eps_dv * eps_dv + fsr2_eps * eps_dvdv)
        * 1.0e-10
}

/// `d^2FSR2/dT dV`.
#[must_use]
pub fn fsr2_dtdv(state: &FurstState) -> f64 {
    let vn = state.molar_volume * state.moles;
    let one_minus = 1.0 - state.packing;
    let eps_dv = -state.packing / vn;
    let fsr2_vw = -1.0 / (vn * vn * one_minus);
    let fsr2_eps_w = 1.0 / (vn * one_minus * one_minus);
    (fsr2_vw * state.w_dt + fsr2_eps_w * eps_dv * state.w_dt) * 1.0e-5
}

/// `T^2 d^2FSR2/dT^2`, which is `FSR2W * WTT`.
#[must_use]
pub fn t2_d2_fsr2_dt2(state: &FurstState) -> f64 {
    state.w_dtdt * state.temperature.powi(2)
        / (state.molar_volume * state.moles * (1.0 - state.packing))
}

/// The MSA long-range term `FLR = -alphaLR2 XLR/(4 pi) + n V gamma^3/(3 pi N_A)`.
#[must_use]
pub fn flr(state: &FurstState) -> f64 {
    let first = -alpha_lr2(state.dielectric, state.temperature) * state.xlr / (4.0 * NEQSIM_PI);
    let second = state.moles * state.molar_volume * state.shielding.powi(3)
        / (3.0 * NEQSIM_PI * NEQSIM_AVOGADRO);
    first + second
}

/// `T dFLR/dT`, NeqSim's `dFLRdT` times `T`.
#[must_use]
pub fn t_d_flr_dt(state: &FurstState) -> f64 {
    let alpha = alpha_lr2(state.dielectric, state.temperature);
    let alpha_dt = alpha_lr2_dt(state);
    let term1 = -1.0 / (4.0 * NEQSIM_PI) * (alpha_dt * state.xlr + alpha * state.xlr_dt);
    let term2 = state.moles * state.molar_volume / (3.0 * NEQSIM_PI * NEQSIM_AVOGADRO)
        * 3.0
        * state.shielding.powi(2)
        * state.shielding_dt;
    (term1 + term2) * state.temperature
}

/// `dFLR/dV`, against the total volume.
#[must_use]
pub fn flr_dv(state: &FurstState) -> f64 {
    let d_f_d_alpha = -state.xlr / (4.0 * NEQSIM_PI);
    let flr_v = state.shielding.powi(3) / (3.0 * NEQSIM_PI * NEQSIM_AVOGADRO);
    (flr_v + d_f_d_alpha * alpha_lr2_dv(state)) * 1.0e-5
}

/// `d^2FLR/dV^2`, against the total volume.
#[must_use]
pub fn flr_dvdv(state: &FurstState) -> f64 {
    let d_f_d_alpha = -state.xlr / (4.0 * NEQSIM_PI);
    let e = ELECTRON_CHARGE;
    let alpha_dvdv = 2.0 * e * e * NEQSIM_AVOGADRO
        / (VACUUM_PERMITTIVITY * state.dielectric.powi(3) * R * state.temperature)
        * state.dielectric_dv
        * state.dielectric_dv
        - e * e * NEQSIM_AVOGADRO
            / (VACUUM_PERMITTIVITY * state.dielectric * state.dielectric * R * state.temperature)
            * state.dielectric_dvdv;
    d_f_d_alpha * alpha_dvdv * 1.0e-10
}

/// `d^2FLR/dT dV`.
#[must_use]
pub fn flr_dtdv(state: &FurstState) -> f64 {
    let d_f_d_alpha = -state.xlr / (4.0 * NEQSIM_PI);
    let e = ELECTRON_CHARGE;
    let eps = state.dielectric;
    let t = state.temperature;
    let alpha_dtdv = e * e * NEQSIM_AVOGADRO / (VACUUM_PERMITTIVITY * eps * eps * R * t * t)
        * state.dielectric_dv
        + 2.0 * e * e * NEQSIM_AVOGADRO / (VACUUM_PERMITTIVITY * eps * eps * eps * R * t)
            * state.dielectric_dt
            * state.dielectric_dv
        - e * e * NEQSIM_AVOGADRO / (VACUUM_PERMITTIVITY * eps * eps * R * t)
            * state.dielectric_dtdv;
    d_f_d_alpha * alpha_dtdv * 1.0e-5
}

/// `T^2 d^2FLR/dT^2`, NeqSim's `dFLRdTdT` times `T^2`.
///
/// **This is NeqSim's own simplification and not the second derivative.** The source says
/// so: "the full second derivative would require `d2gamma/dT2` and `d2XLR/dT2` which are not
/// yet implemented". What it keeps is `-1/(4 pi)(alphaLR'' XLR + 2 alphaLR' XLR')` and the
/// whole of the `gamma` term is missing. A port that completed the derivative would be a
/// different model from the one being ported, and the divergence would be invisible,
/// because both are plausible numbers.
#[must_use]
pub fn t2_d2_flr_dt2(state: &FurstState) -> f64 {
    let (e, e0, na, t, eps) = (
        ELECTRON_CHARGE,
        VACUUM_PERMITTIVITY,
        NEQSIM_AVOGADRO,
        state.temperature,
        state.dielectric,
    );
    // NeqSim's `alphaLRdTdT`. **The source writes its `2 scale eps'/(eps^2 T^2)` term as
    // two identical lines in a row**, and transcribing only one of them is a mistake that
    // stays invisible until a value is compared - which is what happened here. Written as
    // the analytic form instead, because `2 scale/(eps T^3) + 2 scale eps'/(eps^2 T^2)
    // - scale eps''/(eps^2 T) + 2 scale eps'^2/(eps^3 T)` is what differentiates
    // `alphaLR2 = scale/(eps T)` twice, and the two copies are the `2` on the second term.
    let scale = e * e * na / (e0 * R);
    let alpha_dtdt = 2.0 * scale / (eps * t.powi(3))
        + 2.0 * scale * state.dielectric_dt / (eps * eps * t * t)
        - scale * state.dielectric_dtdt / (eps * eps * t)
        + 2.0 * scale * state.dielectric_dt.powi(2) / (eps.powi(3) * t);
    let alpha_dt = alpha_lr2_dt(state);
    let term1 = -1.0 / (4.0 * NEQSIM_PI) * alpha_dtdt * state.xlr;
    let cross = -2.0 / (4.0 * NEQSIM_PI) * alpha_dt * state.xlr_dt;
    (term1 + cross) * t * t
}

/// The Born solvation term `FBorn = K(T)(1/eps_s - 1) bornX`.
#[must_use]
pub fn fborn(state: &FurstState) -> f64 {
    born_scale(state.temperature) * (1.0 / state.solvent_dielectric - 1.0) * state.born_x
}

/// `T dFBorn/dT`, NeqSim's `dFBorndT` times `T`.
///
/// `FBornT + FBornD * eps_s'`, with `FBornT = -(K/T)(1/eps_s - 1) bornX` and
/// `FBornD = -K bornX/eps_s^2`.
#[must_use]
pub fn t_d_fborn_dt(state: &FurstState) -> f64 {
    let k = born_scale(state.temperature);
    let t = state.temperature;
    let d = 1.0 / state.solvent_dielectric - 1.0;
    let f_born_t = -(k / t) * d * state.born_x;
    let f_born_d = -k / state.solvent_dielectric.powi(2) * state.born_x;
    (f_born_t + f_born_d * state.solvent_dielectric_dt) * t
}

/// `T^2 d^2FBorn/dT^2`, NeqSim's `dFBorndTdT` times `T^2`.
///
/// **Also simplified, and more so than the long-range one.** NeqSim keeps
/// `FBornTT + FBornTD eps_s'` with `FBornTT = 2(K/T^2)(1/eps_s - 1) bornX` and
/// `FBornTD = (K/T^2) bornX/eps_s^2`. The `FBornD eps_s''` and `FBornDD (eps_s')^2` terms
/// are absent, so this is not the second derivative of the expression above it - and the
/// two are checked against each other nowhere upstream.
#[must_use]
pub fn t2_d2_fborn_dt2(state: &FurstState) -> f64 {
    let k = born_scale(state.temperature);
    let t = state.temperature;
    let d = 1.0 / state.solvent_dielectric - 1.0;
    // `FBornTT` carries `T^3` and `FBornTD` carries `T^2` - one power apart, because the
    // first differentiates the `1/T` in the scale twice and the second only once. `k` is
    // `C/T`, so `FBornTT = 2k/T^2` and `FBornTD = k/T`. Giving them the same power is the
    // mistake, and it is 0.9% on this number.
    let f_born_tt = 2.0 * k / (t * t) * d * state.born_x;
    let f_born_td = k / t / state.solvent_dielectric.powi(2) * state.born_x;
    (f_born_tt + f_born_td * state.solvent_dielectric_dt) * t * t
}

/// One component's inputs to the composition derivatives.
///
/// Assembled by the phase rather than stored, because three of the four come from somewhere
/// else: the ion diameter is the one the phase resolves (derived from the fitted covolume),
/// the dielectric constant is the component's own at `T`, and `w_i` is a row sum of the
/// short-range table rather than a property of the component at all.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ComponentDerivatives {
    /// The ionic charge, in elementary charges.
    pub charge: f64,
    /// The ion diameter in **metres**, as the phase resolves it. Zero for a component with
    /// none, which makes its Born radius zero rather than infinite.
    pub diameter_m: f64,
    /// The component's own dielectric constant at the state's temperature.
    pub dielectric: f64,
    /// NeqSim's `calcWi`, which is **`-2 sum_j n_j Wij(i, j, T)`** and not the row sum.
    ///
    /// The `-2` is the handler's own: `calcW`, `calcWij` and `calcWi` each return their sum
    /// negated and doubled, so a caller reading them as the sums they are named for is wrong
    /// by that factor in `W`, in the pair table and here. Taken as the handler returns it.
    pub w_i: f64,
}

/// One component's electrolyte contribution to `dFdN_i`, by term.
///
/// Carried apart rather than summed because each is oracled apart: the capture prints
/// `dFSR2dN`, `dFLRdN` and `dFBorndN` per component, and a sum that matched while its parts
/// did not would be three compensating errors.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CompositionContribution {
    /// The short-range `W` term's contribution.
    pub short_range: f64,
    /// The MSA long-range term's.
    pub long_range: f64,
    /// The Born term's.
    pub born: f64,
}

impl CompositionContribution {
    /// The three added, which is what the fugacity coefficient takes.
    #[must_use]
    pub fn total(&self) -> f64 {
        self.short_range + self.long_range + self.born
    }
}

/// The three electrolyte contributions to `dFdN_i`, which is what `ln phi_i` is built from.
///
/// `ComponentEos.fugcoef` is `exp(dFdN - ln(PV/RT))`, and this library's cubic `ln phi` is
/// already `dFdN - ln Z` for its own part - so these are the *changes*, and adding them to
/// the cubic's logarithm gives the electrolyte phase's fugacity coefficient. The identity is
/// checked rather than assumed: the capture's `lnPhi[2]` is `-275.908826771314` and its
/// `dFdN[2]` `-280.551091323773`, which differ by `-ln Z` to the last digit.
///
/// # The gamma chain is missing, upstream
///
/// `dFLRdN` is `FLRXLR XLRi + dFdAlphaLR alphai`, and **the chain through the shielding
/// parameter is commented out in the source** - `FLRGammaLR * gammaLRdn` is written and
/// disabled, and `dFLRdNdN` has the same hole. A port that completed it would have a
/// different fugacity coefficient from the one being ported, so it is reproduced and the
/// omission stated where a reader meets it.
///
/// # Errors
/// * [`AzothError::InvalidInput`] if the two slices are not one per component.
/// * [`AzothError::OutOfRange`] if no neutral component carries any moles, which the solvent
///   dielectric constant's composition derivative divides by.
pub fn ln_phi_contributions(
    state: &FurstState,
    components: &[ComponentDerivatives],
    mole_numbers: &[f64],
) -> Result<Vec<CompositionContribution>> {
    if components.len() != mole_numbers.len() {
        return Err(AzothError::invalid_input(
            "components",
            format!(
                "{} components and {} mole numbers; the composition derivative is one per \
                 component",
                components.len(),
                mole_numbers.len()
            ),
        ));
    }
    let neutral_moles: f64 = components
        .iter()
        .zip(mole_numbers)
        .filter(|(c, _)| c.charge == 0.0)
        .map(|(_, n)| n)
        .sum();
    if !neutral_moles.is_finite() || neutral_moles <= 0.0 {
        return Err(AzothError::out_of_range(
            "mole_numbers",
            neutral_moles,
            "the solvent dielectric constant's composition derivative divides by the neutral \
             components' total moles, and a phase with none has no such derivative",
        ));
    }

    let vn = state.molar_volume * state.moles;
    let one_minus = 1.0 - state.packing;
    let fsr2_eps = state.w / (vn * one_minus * one_minus);
    let fsr2_w = 1.0 / (vn * one_minus);
    let flr_xlr = -alpha_lr2(state.dielectric, state.temperature) / (4.0 * NEQSIM_PI);
    let d_f_d_alpha = -state.xlr / (4.0 * NEQSIM_PI);
    let k = born_scale(state.temperature);
    let f_born_x = k * (1.0 / state.solvent_dielectric - 1.0);
    let f_born_d = -k / state.solvent_dielectric.powi(2) * state.born_x;
    let eps_ionic_half = 1.0 + state.ionic_packing / 2.0;

    let mut out = Vec::with_capacity(components.len());
    for component in components {
        // `dEpsdNi`, and `dEpsIonicdNi` which is the same expression zeroed for a neutral.
        let scale = NEQSIM_AVOGADRO * NEQSIM_PI / 6.0 * component.diameter_m.powi(3) / vn;
        let eps_ionic_i = if component.charge == 0.0 { 0.0 } else { scale };

        // `calcSolventdiElectricdn`: zero for an ion, else the component's own constant less
        // the mixture's, over the neutral moles.
        let solvent_dn = if component.charge != 0.0 {
            0.0
        } else {
            (component.dielectric - state.solvent_dielectric) / neutral_moles
        };

        // `calcdiElectricdn = dYdf X + Y dXdf`.
        let x = (1.0 - state.ionic_packing) / eps_ionic_half;
        let y = state.solvent_dielectric - 1.0;
        let d_x = eps_ionic_i * -1.5 / (eps_ionic_half * eps_ionic_half);
        let dielectric_dn = solvent_dn * x + y * d_x;

        // `alphai = -e^2 N_A/(eps0 eps^2 R T) deps/dn_i`.
        let alpha_i = -ELECTRON_CHARGE * ELECTRON_CHARGE * NEQSIM_AVOGADRO
            / (VACUUM_PERMITTIVITY * state.dielectric.powi(2) * R * state.temperature)
            * dielectric_dn;

        let xlr_i = component.charge.powi(2) * state.shielding
            / (1.0 + state.shielding * component.diameter_m);
        let born_i = if component.diameter_m > 0.0 {
            component.charge.powi(2) / component.diameter_m
        } else {
            0.0
        };

        let fsr2 = fsr2_eps * scale + fsr2_w * component.w_i;
        let flr = flr_xlr * xlr_i + d_f_d_alpha * alpha_i;
        let born = f_born_x * born_i + f_born_d * solvent_dn;
        out.push(CompositionContribution {
            short_range: fsr2,
            long_range: flr,
            born,
        });
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The aqueous phase of `SystemFurstElectrolyteEosTest`'s own mixture, from
    /// `validation/neqsim/captures/furst_probe.tsv`. Every number here is a line of that
    /// capture, so a change to any of them is a change against the oracle rather than
    /// against this crate's own arithmetic.
    fn aqueous() -> FurstState {
        FurstState {
            temperature: 298.15,
            molar_volume: 2.385_525_337_515_67e-5,
            moles: 1.001_901_653_436_76,
            packing: 0.212_659_929_859_501,
            ionic_packing: 0.001_524_427_057_062_24,
            solvent_dielectric: 78.314_881_455_398_7,
            solvent_dielectric_dt: -0.359_218_709_298_880,
            dielectric: 78.138_224_759_715_8,
            dielectric_dt: -0.358_397_930_827_548,
            dielectric_dtdt: 0.001_952_585_159_916_30,
            dielectric_dv: 7_385.673_022_877_76,
            dielectric_dvdv: -617_561_262.424_822,
            dielectric_dtdv: -34.315_152_279_088_9,
            shielding: 302_510_558.954_255,
            shielding_dt: 169_175.175_516_858,
            xlr: 542_989.558_795_128,
            xlr_dt: 272.589_064_172_091,
            born_x: 5_401_911.027_213_43,
            w: -3.254_442_418_358_84e-07,
            w_dt: 3.687_911_635_695_15e-10,
            w_dtdt: -3.188_926_435_637_36e-12,
        }
    }

    fn close(got: f64, want: f64, relative: f64, what: &str) {
        let scale = want.abs().max(f64::MIN_POSITIVE);
        assert!(
            (got - want).abs() / scale < relative,
            "{what}: got {got}, NeqSim gives {want} (relative {})",
            (got - want).abs() / scale
        );
    }

    /// **`alphaLR2` and its derivatives, which are the long-range term's whole scale.**
    ///
    /// `alphaLR2` is checked against the probe directly; the two derivatives are checked
    /// against a central difference of it, because the probe prints neither and a second
    /// transcription of the same algebra would agree with a mistake in the first.
    #[test]
    fn the_long_range_scale_is_the_oracle_and_its_derivatives_are_the_calculus() {
        let state = aqueous();
        close(
            alpha_lr2(state.dielectric, state.temperature),
            9.014_890_778_729_25e-09,
            1.0e-14,
            "alphaLR2",
        );

        // d/dT, against a central difference at fixed V and composition. The dielectric
        // constant has to move with the temperature for this to mean anything, so the
        // difference is taken through a rebuilt state.
        let h = 1.0e-4 * state.temperature;
        let mut ahead = state.clone();
        let mut behind = state.clone();
        for (s, sign) in [(&mut ahead, 1.0), (&mut behind, -1.0)] {
            s.temperature += sign * h;
            s.dielectric = state.dielectric
                + sign * h * state.dielectric_dt
                + 0.5 * h * h * state.dielectric_dtdt;
        }
        let quotient = (alpha_lr2(ahead.dielectric, ahead.temperature)
            - alpha_lr2(behind.dielectric, behind.temperature))
            / (2.0 * h);
        close(alpha_lr2_dt(&state), quotient, 1.0e-6, "dalphaLR2/dT");

        // And d/dV, the same way - but on the dielectric constant directly and not through
        // a temperature-sized step, which at this state would move `eps` from 78 to 298 and
        // measure the curvature instead of the slope.
        let h_eps = 1.0e-6 * state.dielectric;
        let quotient = (alpha_lr2(state.dielectric + h_eps, state.temperature)
            - alpha_lr2(state.dielectric - h_eps, state.temperature))
            / (2.0 * h_eps / state.dielectric_dv);
        close(alpha_lr2_dv(&state), quotient, 1.0e-6, "dalphaLR2/dV");
    }

    /// **The three terms, against the probe's own rows.**
    ///
    /// `FSR2 = -0.0172943844995666`, `FLR = -0.000272971800726073` and
    /// `FBorn = -0.298937482116719` - so the Born term is the largest of the three and the
    /// MSA term the smallest, and all three are live rather than one carrying the others.
    #[test]
    fn the_three_terms_match_the_oracle() {
        let state = aqueous();
        close(
            fsr2(&state).expect("computes"),
            -0.017_294_384_499_566_6,
            1.0e-13,
            "FSR2",
        );
        close(flr(&state), -0.000_272_971_800_726_073, 1.0e-13, "FLR");
        close(fborn(&state), -0.298_937_482_116_719, 1.0e-13, "FBorn");
    }

    /// **Every derivative the phase's own `dFdT`, `dFdV`, `dFdVdV` and `dFdTdT` read.**
    ///
    /// The `T`-scaled rows are the capture's plain derivatives times the power of `T` this
    /// library carries them at, so the comparison is of the whole expression rather than of
    /// a factor.
    #[test]
    fn the_derivatives_match_the_oracle() {
        let state = aqueous();
        let t = state.temperature;

        close(
            t_d_fsr2_dt(&state),
            1.959_787_688_002_81e-05 * t,
            1.0e-13,
            "T dFSR2/dT",
        );
        close(
            fsr2_dv(&state),
            0.009_190_383_382_775_51,
            1.0e-12,
            "dFSR2/dV",
        );
        close(
            fsr2_dvdv(&state),
            -0.009_767_696_181_903_80,
            1.0e-12,
            "d2FSR2/dV2",
        );
        close(
            fsr2_dtdv(&state),
            -1.041_447_887_436_55e-05,
            1.0e-11,
            "dFSR2/dTdV",
        );
        close(
            t2_d2_fsr2_dt2(&state),
            -1.694_622_698_119_71e-07 * t * t,
            1.0e-13,
            "T2 d2FSR2/dT2",
        );

        close(
            t_d_flr_dt(&state),
            -4.801_729_009_900_02e-07 * t,
            1.0e-13,
            "T dFLR/dT",
        );
        close(flr_dv(&state), 4.913_648_794_711_00e-05, 1.0e-12, "dFLR/dV");
        close(
            flr_dvdv(&state),
            -3.085_598_006_341_87e-07,
            1.0e-12,
            "d2FLR/dV2",
        );
        close(
            flr_dtdv(&state),
            4.319_717_531_384_74e-10,
            1.0e-11,
            "dFLR/dTdV",
        );
        close(
            t2_d2_flr_dt2(&state),
            -3.917_030_762_662_23e-09 * t * t,
            1.0e-12,
            "T2 d2FLR/dT2",
        );

        close(
            t_d_fborn_dt(&state),
            0.001_020_376_258_110_79 * t,
            1.0e-13,
            "T dFBorn/dT",
        );
        close(
            t2_d2_fborn_dt2(&state),
            -6.785_233_895_932_52e-06 * t * t,
            1.0e-12,
            "T2 d2FBorn/dT2",
        );
    }

    /// **A shut-off term contributes nothing, and the switch is what says so.**
    ///
    /// The three switches exist because NeqSim has them; a model that ignored them would be
    /// a different model on any phase that set one.
    #[test]
    fn a_switched_off_term_contributes_nothing() {
        let all = TermSwitches::all_on();
        assert!(all.short_range && all.long_range && all.born);
        // The terms themselves do not read the switches - the phase's `getF` multiplies by
        // them - so what this pins is that they are three separate flags and not one.
        let off = TermSwitches { born: false, ..all };
        assert_ne!(off, all);
        assert!(off.short_range && off.long_range && !off.born);
    }

    /// **A packing fraction at or above one is refused rather than divided by.**
    ///
    /// `FSR2`'s denominator is `V n (1 - eps)`, and `eps >= 1` is a phase whose components'
    /// excluded volumes exceed its own. NeqSim divides and returns a large number with the
    /// wrong sign; this refuses, which is the tranche's rule for a domain edge.
    #[test]
    fn a_packing_fraction_at_or_above_one_is_refused() {
        let mut state = aqueous();
        state.packing = 1.0;
        let error = fsr2(&state).expect_err("the excluded volume is the whole phase");
        assert_eq!(error.field(), Some("packing"), "{error:?}");
    }

    /// **The composition derivatives, against the capture's own per-component rows.**
    ///
    /// The three terms are checked apart, because a sum that matched while its parts
    /// cancelled would be three coincidences. The state and the components are the aqueous
    /// phase's: the diameters are the ones the *phase* resolves, which for an ion is the
    /// value derived from the fitted covolume - 4.3437 Å - and not the correlation's 5.68.
    #[test]
    fn the_composition_derivatives_match_the_oracle() {
        let state = FurstState {
            w: -3.254_442_418_358_84e-07,
            ..aqueous()
        };
        // The short-range table, which `w_i` is a row sum of. Built here rather than
        // restated, and it is the same table `tests` in `furst_mixing` checks pairwise.
        let mixture = |name: &str, charge: f64, angstrom: f64, dielectric: f64| {
            crate::furst_mixing::FurstComponent {
                name: name.to_string(),
                charge,
                diameter: angstrom,
                dielectric_at_reference: dielectric,
            }
        };
        let table = crate::furst_mixing::wij_table(
            &[
                mixture("methane", 0.0, 2.52, 2.0),
                mixture("water", 0.0, 2.52, 78.332_147_573_168_0),
                mixture("na+", 1.0, 5.68, 0.0),
                mixture("cl-", -1.0, 3.60, 0.0),
            ],
            crate::databank::furst_wij,
        )
        .expect("the shipped mixture builds a table");
        let moles: Vec<f64> = [
            0.000_225_745_660_581_355,
            0.997_778_050_427_449,
            0.000_998_101_955_985_164,
            0.000_998_101_955_985_164,
        ]
        .iter()
        .map(|x| x * 1.001_901_653_436_76)
        .collect();
        let components = [
            (0.0, 2.52, 2.0),
            (0.0, 2.52, 78.332_147_573_168_0),
            (1.0, 4.343_717_662_170_81, 0.0),
            (-1.0, 3.226_081_589_674_81, 0.0),
        ]
        .iter()
        .enumerate()
        .map(
            |(i, &(charge, angstrom, dielectric))| ComponentDerivatives {
                charge,
                diameter_m: angstrom * 1.0e-10,
                dielectric,
                // `calcWi`: the row sum negated and doubled, which is how the handler
                // returns it - the same `-2` that `calcW` and `calcWij` carry.
                w_i: -2.0
                    * (0..4)
                        .map(|j| moles[j] * table.wij(i, j, state.temperature))
                        .sum::<f64>(),
            },
        )
        .collect::<Vec<_>>();

        let got = ln_phi_contributions(&state, &components, &moles).expect("computes");
        let want: [(f64, f64, f64); 4] = [
            (
                -0.026_957_328_895_281_1,
                -0.000_379_609_469_205_052,
                0.003_768_121_783_018_29,
            ),
            (
                -0.021_959_491_597_569_0,
                8.588_602_480_471_29e-08,
                -8.525_314_228_889_12e-07,
            ),
            (
                -17.315_614_106_443_6,
                -0.192_435_547_491_505,
                -127.400_565_779_252,
            ),
            (
                0.014_107_883_057_089_9,
                -0.197_975_495_079_022,
                -171.536_916_337_467,
            ),
        ];
        for (i, &(short, long, born)) in want.iter().enumerate() {
            let scale = short.abs().max(long.abs()).max(born.abs());
            assert!(
                (got[i].short_range - short).abs() < 1.0e-12 * scale.max(1.0),
                "component {i}: dFSR2dN = {}, the probe prints {short}",
                got[i].short_range
            );
            assert!(
                (got[i].long_range - long).abs() < 1.0e-12 * scale.max(1.0),
                "component {i}: dFLRdN = {}, the probe prints {long}",
                got[i].long_range
            );
            assert!(
                (got[i].born - born).abs() < 1.0e-12 * scale.max(1.0),
                "component {i}: dFBorndN = {}, the probe prints {born}",
                got[i].born
            );
        }
        // And the identity the whole thing exists for: `ln phi_i = dFdN_i - ln Z`.
        let ln_z = 0.009_635_852_009_232_98_f64.ln();
        let d_f_d_n = [
            3.733_727_075_442_27,
            -10.391_048_596_690_5,
            -280.551_091_323_773,
            -171.220_609_647_619,
        ];
        let ln_phi = [
            8.375_991_627_901_69,
            -5.748_784_044_231_08,
            -275.908_826_771_314,
            -166.578_345_095_159,
        ];
        for i in 0..4 {
            assert!(
                (d_f_d_n[i] - ln_z - ln_phi[i]).abs() < 1.0e-9,
                "ln phi[{i}] is {} and dFdN[{i}] - ln Z is {}",
                ln_phi[i],
                d_f_d_n[i] - ln_z
            );
        }
    }

    /// A phase with no neutral component has no solvent-dielectric composition derivative,
    /// so it is refused rather than divided by.
    #[test]
    fn a_phase_with_no_solvent_is_refused() {
        let state = aqueous();
        let components = [ComponentDerivatives {
            charge: 1.0,
            diameter_m: 4.34e-10,
            dielectric: 0.0,
            w_i: 0.0,
        }];
        let error = ln_phi_contributions(&state, &components, &[1.0])
            .expect_err("a lone ion has no solvent to differentiate against");
        assert_eq!(error.field(), Some("mole_numbers"), "{error:?}");
    }
}
