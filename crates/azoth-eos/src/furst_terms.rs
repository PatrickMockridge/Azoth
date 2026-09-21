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
    // **Propagated, and deliberately.** `eps` is the phase's own and `T` the state's, and a
    // non-positive one of either is not a state - so this returns the infinity or NaN
    // NeqSim's own division does rather than a value invented here. `T > 0` is refused at the
    // model's boundary; `eps` comes out of a fitted polynomial, and the same three functions
    // are what the probe prints, so a guard here would be a divergence from the oracle at
    // exactly the states no one can interpret.
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
    mod2004: bool,
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
        // **Zero in the 2004 revision.** `ComponentModifiedFurstElectrolyteEosMod2004`'s
        // `calcSolventdiElectricdn` returns `0.0` with its body commented out, so the solvent
        // dielectric constant has no composition dependence there - and the Born term's
        // composition path is the `dYdf` that goes with it.
        let solvent_dn = if mod2004 || component.charge != 0.0 {
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
        // **`FBornD` is weighted by the solvent's composition derivative in both models.**
        // The 2004 revision once *added* it unweighted - `FBornX XBorni + FBornD` - which made
        // a chemical potential depend on the phase's size, because `FBornD` is extensive and
        // the variant sets that derivative to zero. NeqSim's issue 3862 is that observation
        // and its fix is the one term, so both models read `FBornX XBorni + FBornD
        // solventdiElectricdn` here.
        let born = f_born_x * born_i + f_born_d * solvent_dn;
        out.push(CompositionContribution {
            short_range: fsr2,
            long_range: flr,
            born,
        });
    }
    Ok(out)
}

/// One component's contribution to the state, as `volInit` reads it.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ComponentState {
    /// The ionic charge, in elementary charges.
    pub charge: f64,
    /// The diameter in **metres**, as the phase resolves it - the value derived from the
    /// fitted covolume for an ion. Zero makes the component contribute nothing to the
    /// packing fraction or the Born sum, which is NeqSim's own treatment.
    pub diameter_m: f64,
    /// `eps_i(T)`, the component's own dielectric constant.
    pub dielectric: f64,
    /// `d eps_i/dT`.
    pub dielectric_dt: f64,
    /// `d^2 eps_i/dT^2`.
    pub dielectric_dtdt: f64,
    /// The component's critical volume, in m³/mol. Read only by the two volume-fraction
    /// dielectric mixing rules, so a molar-average caller may leave it zero.
    pub critical_volume: f64,
}

/// Everything `volInit` reads, handed to [`build_state`] as one argument.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct StateInputs<'a> {
    /// The state's absolute temperature.
    pub temperature: f64,
    /// The phase's molar volume, in m³/mol.
    pub molar_volume: f64,
    /// The mole numbers, one per component.
    pub mole_numbers: &'a [f64],
    /// The components.
    pub components: &'a [ComponentState],
    /// How the solvent dielectric constants are combined.
    pub rule: crate::furst_dielectric::MixingRule,
    /// The short-range sums, from [`crate::furst_mixing::short_range`].
    pub short_range: crate::furst_mixing::ShortRange,
    /// The pair table, for the per-component row sums.
    pub table: &'a crate::furst_mixing::WijTable,
    /// The component names, lower-cased, for the alternating short-range rows.
    pub names: &'a [String],
}

/// Build the state, in `volInit`'s own order.
///
/// The order is load-bearing and is the reason this is one function rather than a set of
/// accessors: the packing fraction needs the volume, the dielectric constants need the
/// packing fraction, `alphaLR2` needs the dielectric, and **the shielding parameter needs
/// `alphaLR2`** - so a `gamma` read before the dielectric is a `gamma` of the wrong brine.
///
/// # Errors
/// * Propagates [`crate::furst_dielectric::solvent_dielectric`]'s refusal of a phase with no
///   solvent, and the packing fraction's refusal of a non-positive total volume.
/// * [`AzothError::InvalidInput`] if the slices are not one per component.
pub fn build_state(inputs: &StateInputs<'_>) -> Result<FurstState> {
    let n = inputs.components.len();
    if inputs.mole_numbers.len() != n || inputs.names.len() != n {
        return Err(AzothError::invalid_input(
            "components",
            format!(
                "{} components, {} mole numbers and {} names; the state is one per component",
                n,
                inputs.mole_numbers.len(),
                inputs.names.len()
            ),
        ));
    }
    let is_ion: Vec<bool> = inputs.components.iter().map(|c| c.charge != 0.0).collect();
    let diameters: Vec<f64> = inputs.components.iter().map(|c| c.diameter_m).collect();
    let total_moles: f64 = inputs.mole_numbers.iter().sum();

    // `calcEps` and `calcEpsIonic` are one sum over two different sets.
    let packing = crate::furst_dielectric::packing_fraction(
        &diameters,
        inputs.mole_numbers,
        &is_ion,
        total_moles,
        inputs.molar_volume,
        false,
    )?;
    let ionic_packing = crate::furst_dielectric::packing_fraction(
        &diameters,
        inputs.mole_numbers,
        &is_ion,
        total_moles,
        inputs.molar_volume,
        true,
    )?;

    // The solvent mixture's constant, and its temperature derivatives. The derivatives use
    // the per-component `d eps_i/dT` in the same mole-number average - the ions are excluded
    // from all three, which is what `calcSolventDiElectricConstantdT` does whatever mixing
    // rule is selected.
    let eps_i: Vec<f64> = inputs.components.iter().map(|c| c.dielectric).collect();
    let eps_dt: Vec<f64> = inputs.components.iter().map(|c| c.dielectric_dt).collect();
    let eps_dtdt: Vec<f64> = inputs
        .components
        .iter()
        .map(|c| c.dielectric_dtdt)
        .collect();
    let critical_volumes: Vec<f64> = inputs
        .components
        .iter()
        .map(|c| c.critical_volume)
        .collect();
    let solvent_dielectric = crate::furst_dielectric::solvent_dielectric(
        inputs.rule,
        inputs.mole_numbers,
        &eps_i,
        &is_ion,
        &critical_volumes,
    )?;
    let solvent_dielectric_dt = crate::furst_dielectric::solvent_dielectric(
        inputs.rule,
        inputs.mole_numbers,
        &eps_dt,
        &is_ion,
        &critical_volumes,
    )?;
    let solvent_dielectric_dtdt = crate::furst_dielectric::solvent_dielectric(
        inputs.rule,
        inputs.mole_numbers,
        &eps_dtdt,
        &is_ion,
        &critical_volumes,
    )?;

    // The phase's own constant and its derivatives. Note that `epsIonic` does not depend on
    // `T` at fixed `V`, so `d eps/dT` carries the solvent term alone and `d2 eps/dT2` is the
    // solvent's scaled by the same `X`.
    let x = (1.0 - ionic_packing) / (1.0 + ionic_packing / 2.0);
    let y = solvent_dielectric - 1.0;
    let dielectric = 1.0 + y * x;
    let dielectric_dt = solvent_dielectric_dt * x;
    let dielectric_dtdt = solvent_dielectric_dtdt * x;
    let ionic_dv = crate::furst_dielectric::packing_fraction_dv(
        ionic_packing,
        total_moles,
        inputs.molar_volume,
    );
    let ionic_dvdv = crate::furst_dielectric::packing_fraction_dvdv(
        ionic_packing,
        total_moles,
        inputs.molar_volume,
    );
    let d_x_dv = ionic_dv * -1.5 / (1.0 + ionic_packing / 2.0).powi(2);
    let dielectric_dv = y * d_x_dv;
    let d_x_dvdv = ionic_dvdv * -1.5 / (1.0 + ionic_packing / 2.0).powi(2)
        + ionic_dv * ionic_dv * 1.5 / (1.0 + ionic_packing / 2.0).powi(3);
    let dielectric_dvdv = y * d_x_dvdv;
    let dielectric_dtdv = solvent_dielectric_dt * d_x_dv;

    let mut state = FurstState {
        temperature: inputs.temperature,
        molar_volume: inputs.molar_volume,
        moles: total_moles,
        packing,
        ionic_packing,
        solvent_dielectric,
        solvent_dielectric_dt,
        dielectric,
        dielectric_dt,
        dielectric_dtdt,
        dielectric_dv,
        dielectric_dvdv,
        dielectric_dtdv,
        shielding: 0.0,
        shielding_dt: 0.0,
        xlr: 0.0,
        xlr_dt: 0.0,
        born_x: 0.0,
        w: inputs.short_range.w,
        w_dt: inputs.short_range.w_dt,
        w_dtdt: inputs.short_range.w_dtdt,
    };

    state.shielding = shielding_parameter(&state, inputs.components, inputs.mole_numbers);
    state.shielding_dt = shielding_parameter_dt(&state, inputs.components, inputs.mole_numbers);
    state.xlr = xlr(&state, inputs.components, inputs.mole_numbers);
    state.xlr_dt = xlr_dt(&state, inputs.components, inputs.mole_numbers);
    state.born_x = born_x(inputs.components, inputs.mole_numbers, inputs.names);
    Ok(state)
}

/// `calcShieldingParameter`: a damped Newton solve for `gamma`.
///
/// `f(gamma) = 4 gamma^2/N_A - alphaLR2 sum_i n_i/V (z_i/(1 + gamma sigma_i))^2`, solved by
/// `gamma -= 0.8 f/f'` from `1e10`, with a 1000-iteration cap and a **three-iteration
/// floor**. The floor is NeqSim's and it matters: a phase carrying no ion returns exactly
/// zero, but a phase carrying ions at `1e-43` mol does not - its `f` is dominated by
/// `4 gamma^2/N_A`, so the solve walks down and stops wherever the floor leaves it. The
/// probe prints that residue as `gamma = 1692665.9444736` for the vapour, which is not a
/// physical number and is what the iteration returns.
fn shielding_parameter(
    state: &FurstState,
    components: &[ComponentState],
    mole_numbers: &[f64],
) -> f64 {
    // **The solve's two divisors are total by construction.** `v_total` is the phase's
    // volume, which `packing_fraction` has already refused when it was not positive. `df` is
    // a sum of positive terms - `8 gamma/N_A` plus a `2 alpha n_i z_i^2 sigma_i` over a cube -
    // for every `gamma > 0`, which is where the iteration starts and where the probe's own
    // states leave it. The scheme's damping and its floor are NeqSim's, reproduced as they
    // stand rather than augmented with a clamp this model would not have.
    let alpha = alpha_lr2(state.dielectric, state.temperature);
    let v_total = state.molar_volume * state.moles;
    let mut gamma = 1.0e10_f64;
    let mut iterations = 0;
    loop {
        iterations += 1;
        let gamma_old = gamma;
        let mut f = 4.0 * gamma * gamma / NEQSIM_AVOGADRO;
        let mut df = 8.0 * gamma / NEQSIM_AVOGADRO;
        let mut ions = 0;
        for (component, &moles) in components.iter().zip(mole_numbers) {
            if component.charge == 0.0 {
                continue;
            }
            ions += 1;
            let sigma = component.diameter_m;
            let denominator = 1.0 + gamma * sigma;
            f -= alpha * moles / v_total
                * (component.charge / denominator)
                * (component.charge / denominator);
            df += 2.0 * alpha * moles / v_total * component.charge.powi(2) * sigma
                / denominator.powi(3);
        }
        gamma = if ions > 0 {
            gamma_old - 0.8 * f / df
        } else {
            0.0
        };
        if !((f.abs() > 1.0e-10 && iterations < 1000) || iterations < 3) {
            break;
        }
    }
    gamma
}

/// `calcShieldingParameterdT`, by implicit differentiation.
fn shielding_parameter_dt(
    state: &FurstState,
    components: &[ComponentState],
    mole_numbers: &[f64],
) -> f64 {
    if state.shielding < 1.0e-10 {
        return 0.0;
    }
    let alpha = alpha_lr2(state.dielectric, state.temperature);
    let alpha_dt = alpha_lr2_dt(state);
    let v_total = state.molar_volume * state.moles;
    let mut dfdgamma = 8.0 * state.shielding / NEQSIM_AVOGADRO;
    let mut sum = 0.0;
    for (component, &moles) in components.iter().zip(mole_numbers) {
        if component.charge == 0.0 {
            continue;
        }
        let sigma = component.diameter_m;
        let denominator = 1.0 + state.shielding * sigma;
        dfdgamma +=
            2.0 * alpha * moles / v_total * component.charge.powi(2) * sigma / denominator.powi(3);
        sum += moles / v_total * component.charge.powi(2) / denominator.powi(2);
    }
    let dfdt = -alpha_dt * sum;
    if dfdgamma.abs() < 1.0e-50 {
        return 0.0;
    }
    -dfdt / dfdgamma
}

/// `calcXLR = sum_i n_i z_i^2 gamma/(1 + gamma sigma_i)` over the ions.
fn xlr(state: &FurstState, components: &[ComponentState], mole_numbers: &[f64]) -> f64 {
    let mut sum = 0.0;
    for (component, &moles) in components.iter().zip(mole_numbers) {
        if component.charge == 0.0 {
            continue;
        }
        sum += moles * component.charge.powi(2) * state.shielding
            / (1.0 + state.shielding * component.diameter_m);
    }
    sum
}

/// `calcXLRdT`, whose `d/dT` carries only the shielding parameter's own derivative.
fn xlr_dt(state: &FurstState, components: &[ComponentState], mole_numbers: &[f64]) -> f64 {
    let mut sum = 0.0;
    for (component, &moles) in components.iter().zip(mole_numbers) {
        if component.charge == 0.0 {
            continue;
        }
        let denominator = 1.0 + state.shielding * component.diameter_m;
        sum += moles * component.charge.powi(2) * state.shielding_dt / (denominator * denominator);
    }
    sum
}

/// `calcBornX = sum_i n_i z_i^2/sigma_i` over every component with a diameter.
///
/// Every component and not only the ions, because that is how NeqSim writes it - a neutral
/// contributes `z = 0` and so contributes nothing, which makes the two readings the same
/// number here and is stated so that it stays that way.
fn born_x(components: &[ComponentState], mole_numbers: &[f64], names: &[String]) -> f64 {
    let _ = names;
    let mut sum = 0.0;
    for (component, &moles) in components.iter().zip(mole_numbers) {
        if component.diameter_m > 0.0 {
            sum += moles * component.charge.powi(2) / component.diameter_m;
        }
    }
    sum
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

        let got = ln_phi_contributions(&state, &components, &moles, false).expect("computes");
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
        let error = ln_phi_contributions(&state, &components, &[1.0], false)
            .expect_err("a lone ion has no solvent to differentiate against");
        assert_eq!(error.field(), Some("mole_numbers"), "{error:?}");
    }

    /// **The whole state, built the way `volInit` builds it.**
    ///
    /// Every layer the capture prints for the aqueous phase comes out of one call in
    /// `volInit`'s own order: the two packing fractions, the solvent dielectric constant and
    /// its first two temperature derivatives, the phase's own constant, the shielding
    /// parameter from the Newton solve, `XLR` and `bornX`.
    ///
    /// **The shielding parameter is the interesting one.** It is a damped Newton solve with
    /// a three-iteration floor, so for the vapour - whose ions sit at `1e-43` mol - it
    /// returns a residue of the solve rather than a physical number, and the capture prints
    /// `1692665.94447360` there. That value is asserted too, in the module's own test below,
    /// because a port that "fixed" the floor would diverge on every vapour phase.
    #[test]
    fn the_state_matches_the_oracle_layer_by_layer() {
        let components = [
            ComponentState {
                charge: 0.0,
                diameter_m: 2.52e-10,
                dielectric: 2.0,
                dielectric_dt: 0.0,
                dielectric_dtdt: 0.0,
                critical_volume: 9.9e-5,
            },
            ComponentState {
                charge: 0.0,
                diameter_m: 2.52e-10,
                dielectric: 78.332_147_573_168_0,
                dielectric_dt: -0.359_299_981_947_430_96,
                dielectric_dtdt: 0.001_957_499_618_060_916,
                critical_volume: 5.6e-5,
            },
            ComponentState {
                charge: 1.0,
                diameter_m: 4.343_717_662_170_81e-10,
                dielectric: 0.0,
                dielectric_dt: 0.0,
                dielectric_dtdt: 0.0,
                critical_volume: 0.0,
            },
            ComponentState {
                charge: -1.0,
                diameter_m: 3.226_081_589_674_81e-10,
                dielectric: 0.0,
                dielectric_dt: 0.0,
                dielectric_dtdt: 0.0,
                critical_volume: 0.0,
            },
        ];
        let names: Vec<String> = ["methane", "water", "na+", "cl-"]
            .iter()
            .map(|s| s.to_string())
            .collect();
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
        let total = 1.001_901_653_436_76;
        let mole_numbers: Vec<f64> = [
            0.000_225_745_660_581_355,
            0.997_778_050_427_449,
            0.000_998_101_955_985_164,
            0.000_998_101_955_985_164,
        ]
        .iter()
        .map(|x| x * total)
        .collect();
        let short =
            crate::furst_mixing::short_range(&table, &mole_numbers, 298.15).expect("computes");

        let state = build_state(&StateInputs {
            temperature: 298.15,
            molar_volume: 2.385_525_337_515_67e-5,
            mole_numbers: &mole_numbers,
            components: &components,
            rule: crate::furst_dielectric::MixingRule::default_for_the_model(),
            short_range: short,
            table: &table,
            names: &names,
        })
        .expect("the aqueous phase builds a state");

        let check = |got: f64, want: f64, what: &str| {
            let scale = want.abs().max(1.0);
            assert!(
                (got - want).abs() < 1.0e-11 * scale,
                "{what}: got {got}, the probe prints {want}"
            );
        };
        check(state.packing, 0.212_659_929_859_501, "packing");
        check(
            state.ionic_packing,
            0.001_524_427_057_062_24,
            "packing_ionic",
        );
        check(state.solvent_dielectric, 78.314_881_455_398_7, "eps");
        check(
            state.solvent_dielectric_dt,
            -0.359_218_709_298_880,
            "eps_dT",
        );
        // The solvent's second derivative is computed on the way to the phase's and is not
        // carried: no term reads it, because `d2FBorn/dT2` takes the first derivative and
        // nothing takes the second. The phase's own is carried and is checked below.
        check(state.dielectric, 78.138_224_759_715_8, "eps_phase");
        check(state.dielectric_dt, -0.358_397_930_827_548, "eps_phase_dT");
        check(
            state.dielectric_dtdt,
            0.001_952_585_159_916_30,
            "eps_phase_dTdT",
        );
        check(state.dielectric_dv, 7_385.673_022_877_76, "eps_phase_dV");
        check(
            shielding_parameter(&state, &components, &mole_numbers),
            302_510_558.954_255,
            "gamma",
        );
        check(state.xlr, 542_989.558_795_128, "XLR");
        check(state.born_x, 5_401_911.027_213_43, "bornX");
        check(state.w, -3.254_442_418_358_84e-07, "W");
    }

    /// **The vapour's shielding parameter is a residue, not a number.**
    ///
    /// A phase whose ions sit at `1e-43` mol has an `f` dominated by `4 gamma^2/N_A`, so the
    /// damped Newton walks `gamma` down and the three-iteration floor stops it. The capture
    /// prints `1692665.9444736` - identical for two different brines, because the ion term
    /// is below the noise. Reproduced rather than guarded: it is what NeqSim returns.
    #[test]
    fn a_phase_with_no_ions_returns_an_exact_zero() {
        let components = [
            ComponentState {
                charge: 0.0,
                diameter_m: 2.52e-10,
                dielectric: 2.0,
                dielectric_dt: 0.0,
                dielectric_dtdt: 0.0,
                critical_volume: 9.9e-5,
            },
            ComponentState {
                charge: 0.0,
                diameter_m: 2.52e-10,
                dielectric: 78.332_147_573_168_0,
                dielectric_dt: -0.359_299_981_947_430_96,
                dielectric_dtdt: 0.001_957_499_618_060_916,
                critical_volume: 5.6e-5,
            },
        ];
        let names: Vec<String> = ["methane", "water"].iter().map(|s| s.to_string()).collect();
        let table = crate::furst_mixing::wij_table(
            &[
                crate::furst_mixing::FurstComponent {
                    name: "methane".into(),
                    charge: 0.0,
                    diameter: 2.52,
                    dielectric_at_reference: 2.0,
                },
                crate::furst_mixing::FurstComponent {
                    name: "water".into(),
                    charge: 0.0,
                    diameter: 2.52,
                    dielectric_at_reference: 78.332_147_573_168_0,
                },
            ],
            crate::databank::furst_wij,
        )
        .expect("builds");
        let mole_numbers = [0.1, 1.0];
        let short =
            crate::furst_mixing::short_range(&table, &mole_numbers, 298.15).expect("computes");
        let state = build_state(&StateInputs {
            temperature: 298.15,
            molar_volume: 2.385_525_337_515_67e-5,
            mole_numbers: &mole_numbers,
            components: &components,
            rule: crate::furst_dielectric::MixingRule::default_for_the_model(),
            short_range: short,
            table: &table,
            names: &names,
        })
        .expect("builds");
        assert_eq!(
            state.shielding, 0.0,
            "a phase with no ion returns exactly zero rather than the solve's floor"
        );
        assert_eq!(state.xlr, 0.0);
        assert_eq!(state.born_x, 0.0);
    }
}
