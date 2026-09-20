//! The dielectric surface a Furst electrolyte phase is built on.
//!
//! `PhaseModifiedFurstElectrolyteEos` adds three Helmholtz terms to an SRK cubic, and all
//! three read the **solvent's dielectric constant** and the **packing fraction**. This
//! module is those two quantities, and it is here rather than in the model for the reason
//! `electrolyte.rs` is: the arithmetic is the same for every term that reads it.
//!
//! ```text
//! eps_i(T)  = d0 + d1/T + d2 T + d3 T^2 + d4 T^3          per component, DIELECTRICPARAMETER1..5
//! eps_solv  = sum_i n_i eps_i(T) / sum_i n_i              over the components with no charge
//! eps_ionic = (N_A pi/6) sum_i n_i sigma_i^3 / (n V)      over the ions
//! eps       = 1 + (eps_solv - 1)(1 - eps_ionic)/(1 + eps_ionic/2)
//! ```
//!
//! # Two of NeqSim's constants are not the standard ones, and both are load-bearing here
//!
//! `avagadroNumber` is `6.023e23` where CODATA says `6.02214076e23`, and `pi` is
//! `3.14159265` where `Math.PI` carries sixteen digits. The ion radius test in
//! `tests/databank.rs` is where that was found, and the same two constants appear in
//! `packing_fraction` - so a port that wrote `std::f64::consts::PI` would be wrong in
//! every packing fraction, every shielding parameter and every Born radius, by an amount
//! no plausibility check would flag.
//!
//! # What is deliberately not here
//!
//! **The three other mixing rules.** `DielectricMixingRule` declares `OSTER` and
//! `LICHTENECKER` and `calcSolventDiElectricConstant`'s switch routes neither, so they
//! are unreachable. `VOLUME_AVERAGE` and `LOOYENGA` are implemented and reachable, and
//! the class's own comment warns that neither is thermodynamically consistent with
//! complete composition derivatives - which is why [`MixingRule::MolarAverage`] is the
//! default and the only one whose derivatives this library assembles.

use azoth_core::{AzothError, Result};

/// NeqSim's `avagadroNumber`, from `ThermodynamicConstantsInterface`.
///
/// Not CODATA's `6.02214076e23`. See the module note.
pub const NEQSIM_AVOGADRO: f64 = 6.023e23;

/// NeqSim's `pi`, from `ThermodynamicConstantsInterface`.
///
/// Not `Math.PI`. See the module note. The lint that would rather this were
/// `std::f64::consts::PI` is the mistake it exists to prevent.
#[allow(clippy::approx_constant)]
pub const NEQSIM_PI: f64 = 3.14159265;

/// How the solvents' dielectric constants are combined.
///
/// `MOLAR_AVERAGE` is NeqSim's default and the only rule whose temperature derivative it
/// computes consistently - `calcSolventDiElectricConstantdT` is the molar-average formula
/// whatever the rule, so under the other two the value and its derivative disagree. A
/// model that selected one of those would be reproducing that disagreement rather than a
/// choice, so both are carried here and neither is default.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MixingRule {
    /// `eps = sum_i x_i eps_i`.
    ///
    /// Over the components **with no ionic charge**: the ions are excluded from the
    /// average, which matters only through the denominator, since an ion's own
    /// `eps_i` is zero in the table.
    MolarAverage,
    /// `eps = sum_i phi_i eps_i` with `phi_i = n_i Vc_i / sum_j n_j Vc_j`.
    ///
    /// `Vc` is the component's **critical volume**, used as a stand-in for a molar
    /// volume because the phase's own is not yet solved when `volInit` runs.
    VolumeAverage,
    /// `eps = (sum_i phi_i eps_i^(1/3))^3`, on the same volume fractions.
    Looyenga,
}

impl MixingRule {
    /// The rule `PhaseModifiedFurstElectrolyteEos` is constructed with.
    #[must_use]
    pub fn default_for_the_model() -> Self {
        Self::MolarAverage
    }

    /// The rule's own name, as NeqSim's enum spells it.
    #[must_use]
    pub fn as_str(self) -> &'static str {
        match self {
            Self::MolarAverage => "MOLAR_AVERAGE",
            Self::VolumeAverage => "VOLUME_AVERAGE",
            Self::Looyenga => "LOOYENGA",
        }
    }
}

/// `eps_i(T)`, one component's dielectric constant from its five fitted coefficients.
///
/// `DIELECTRICPARAMETER1..5`, whose units are `dimensionless`, `K`, `1/K`, `1/K**2` and
/// `1/K**3` - so every term is dimensionless and `T` enters as written.
///
/// Total in `T` except at `T = 0`, where the `d1/T` term is undefined; the models that
/// call this refuse a non-positive temperature, so the domain is the caller's and it is
/// stated on each of them.
#[must_use]
pub fn component_dielectric(coefficients: &[f64; 5], temperature: f64) -> f64 {
    let [d0, d1, d2, d3, d4] = *coefficients;
    d0 + d1 / temperature
        + d2 * temperature
        + d3 * temperature * temperature
        + d4 * temperature.powi(3)
}

/// `d eps_i / dT`.
///
/// The integer powers are `powi` rather than `powf`, per the numerical-safety rule;
/// NeqSim writes `Math.pow(temperature, 2.0)`, which agrees with `powi(2)` to the last
/// bit for every double.
#[must_use]
pub fn component_dielectric_dt(coefficients: &[f64; 5], temperature: f64) -> f64 {
    let [_, d1, d2, d3, d4] = *coefficients;
    -d1 / temperature.powi(2) + d2 + 2.0 * d3 * temperature + 3.0 * d4 * temperature.powi(2)
}

/// `d^2 eps_i / dT^2`.
#[must_use]
pub fn component_dielectric_dtdt(coefficients: &[f64; 5], temperature: f64) -> f64 {
    let [_, d1, _, d3, d4] = *coefficients;
    2.0 * d1 / temperature.powi(3) + 2.0 * d3 + 6.0 * d4 * temperature
}

/// The solvent mixture's dielectric constant under the rule named.
///
/// `mole_numbers` and `eps_i` are indexed by component and `is_ion` says which entries the
/// average skips. `critical_volumes` is read only by the two volume-fraction rules, so a
/// caller of [`MixingRule::MolarAverage`] may pass an empty slice.
///
/// # Errors
/// * [`AzothError::InvalidInput`] if the three slices are not one length.
/// * [`AzothError::OutOfRange`] if no solvent component carries any moles, which is the
///   divisor. NeqSim returns `ans1 / 1e-50` there - a number of order `1e50` that is not a
///   dielectric constant - rather than refusing.
pub fn solvent_dielectric(
    rule: MixingRule,
    mole_numbers: &[f64],
    eps_i: &[f64],
    is_ion: &[bool],
    critical_volumes: &[f64],
) -> Result<f64> {
    checked_lengths(mole_numbers, eps_i, is_ion)?;
    let solvent_total: f64 = mole_numbers
        .iter()
        .zip(is_ion)
        .filter(|(_, ion)| !**ion)
        .map(|(n, _)| n)
        .sum();
    if !solvent_total.is_finite() || solvent_total <= 0.0 {
        return Err(AzothError::out_of_range(
            "mole_numbers",
            solvent_total,
            "a dielectric constant is a mole-number average over the solvent components, \
             and none of them carries any moles. NeqSim divides by its `1e-50` floor here, \
             which returns a number of order `1e50` that is not a dielectric constant",
        ));
    }
    match rule {
        MixingRule::MolarAverage => Ok(mole_numbers
            .iter()
            .zip(eps_i)
            .zip(is_ion)
            .filter(|(_, ion)| !**ion)
            .map(|((n, eps), _)| n * eps)
            .sum::<f64>()
            / solvent_total),
        MixingRule::VolumeAverage | MixingRule::Looyenga => {
            if critical_volumes.len() != mole_numbers.len() {
                return Err(AzothError::invalid_input(
                    "critical_volumes",
                    format!(
                        "the rule {} weights by volume fraction and needs one critical \
                         volume per component, but got {} for {} components",
                        rule.as_str(),
                        critical_volumes.len(),
                        mole_numbers.len()
                    ),
                ));
            }
            let total_volume: f64 = mole_numbers
                .iter()
                .zip(critical_volumes)
                .zip(is_ion)
                .filter(|(_, ion)| !**ion)
                .map(|((n, vc), _)| n * vc)
                .sum();
            if !total_volume.is_finite() || total_volume <= 0.0 {
                return Err(AzothError::out_of_range(
                    "critical_volumes",
                    total_volume,
                    "the volume-weighted dielectric constant divides by the solvent's total \
                     critical volume and it came out zero",
                ));
            }
            let sum: f64 = mole_numbers
                .iter()
                .zip(eps_i)
                .zip(critical_volumes)
                .zip(is_ion)
                .filter(|(_, ion)| !**ion)
                .map(|(((n, eps), vc), _)| {
                    let fraction = n * vc / total_volume;
                    match rule {
                        MixingRule::Looyenga => fraction * eps.cbrt(),
                        _ => fraction * eps,
                    }
                })
                .sum();
            Ok(match rule {
                MixingRule::Looyenga => sum.powi(3),
                _ => sum,
            })
        }
    }
}

/// The packing fraction `(N_A pi/6) sum_i n_i sigma_i^3 / (n V)`.
///
/// `diameters` are in **metres** - the table carries ångströms and the caller converts -
/// `molar_volume` in m³/mol and `total_moles` the phase's mole number. `ions_only`
/// selects `calcEpsIonic` over `calcEps`, which is the same sum over a different set.
///
/// # Errors
/// * [`AzothError::InvalidInput`] if the two slices are not one length.
/// * [`AzothError::OutOfRange`] if `total_moles * molar_volume` is not positive. It is the
///   divisor and it is a volume, so zero or below is not a state.
pub fn packing_fraction(
    diameters_m: &[f64],
    mole_numbers: &[f64],
    is_ion: &[bool],
    total_moles: f64,
    molar_volume: f64,
    ions_only: bool,
) -> Result<f64> {
    if diameters_m.len() != mole_numbers.len() || is_ion.len() != mole_numbers.len() {
        return Err(AzothError::invalid_input(
            "diameters",
            format!(
                "{} diameters, {} mole numbers and {} ion flags; the three are one per \
                 component",
                diameters_m.len(),
                mole_numbers.len(),
                is_ion.len()
            ),
        ));
    }
    let volume = total_moles * molar_volume;
    if !volume.is_finite() || volume <= 0.0 {
        return Err(AzothError::out_of_range(
            "molar_volume",
            volume,
            "the packing fraction is a volume ratio and its denominator is the phase's \
             total volume, which must be positive",
        ));
    }
    let sum: f64 = diameters_m
        .iter()
        .zip(mole_numbers)
        .zip(is_ion)
        .filter(|(_, ion)| !ions_only || **ion)
        .map(|((diameter, n), _)| n * diameter.powi(3))
        .sum();
    Ok(NEQSIM_AVOGADRO * NEQSIM_PI / 6.0 * sum / volume)
}

/// `d eps / dV` for a packing fraction, at constant composition.
///
/// `eps` is `C / V`, so the derivative is `-eps / V` - and it is written here rather than
/// in the caller because `calcEpsV` and `calcEpsIonicdV` are the same expression and two
/// copies of a derivative are two chances to disagree.
#[must_use]
pub fn packing_fraction_dv(packing: f64, total_moles: f64, molar_volume: f64) -> f64 {
    -packing / (total_moles * molar_volume)
}

/// `d^2 eps / dV^2`, which for `C / V` is `2 eps / V^2`.
#[must_use]
pub fn packing_fraction_dvdv(packing: f64, total_moles: f64, molar_volume: f64) -> f64 {
    2.0 * packing / (total_moles * molar_volume).powi(2)
}

/// The phase's dielectric constant: `1 + (eps_solv - 1)(1 - eps_ionic)/(1 + eps_ionic/2)`.
///
/// `calcDiElectricConstant`, which is *not* the solvent average the MSA and Born terms
/// read. It is the phase's effective constant with the ions' own volume excluded, and
/// `ComponentModifiedFurstElectrolyteEos` is what reads this one.
///
/// `eps_solv` is the **mixture** constant, so if the solvent average is itself not finite
/// the result is not; the caller states that, this does not re-check it.
#[must_use]
pub fn phase_dielectric(solvent: f64, ionic_packing: f64) -> f64 {
    1.0 + (solvent - 1.0) * (1.0 - ionic_packing) / (1.0 + ionic_packing / 2.0)
}

/// Three slices that must be one per component.
fn checked_lengths(a: &[f64], b: &[f64], c: &[bool]) -> Result<()> {
    if a.len() == b.len() && b.len() == c.len() {
        return Ok(());
    }
    Err(AzothError::invalid_input(
        "components",
        format!(
            "{} mole numbers, {} dielectric constants and {} ion flags; the three are one \
             per component",
            a.len(),
            b.len(),
            c.len()
        ),
    ))
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Aqueous and methanol dielectric coefficients, retyped from the compiled table so
    /// this test does not agree with any table at all - the values `COMP.csv` carries.
    const WATER: [f64; 5] = [-19.2905, 29814.5, -0.019678, 0.000132, -3.11e-07];
    const METHANE: [f64; 5] = [2.0, 0.0, 0.0, 0.0, 0.0];
    const METHANOL: [f64; 5] = [-5.3743, 11351.85, 0.0, 0.0, 0.0];

    /// **The polynomial and the mixture, against the oracle.**
    ///
    /// `validation/neqsim/captures/furst_probe.tsv` prints both for the aqueous phase of
    /// `SystemFurstElectrolyteEosTest`'s own mixture, and these are those numbers: the
    /// per-component constants, then the mole-number average over the two solvents, then
    /// its first two temperature derivatives. The ions are excluded from the average and
    /// the methane beside water is what makes it a test rather than a restatement of
    /// water's own value.
    #[test]
    fn the_solvent_dielectric_matches_the_oracle() {
        let t = 298.15;
        assert!(
            (component_dielectric(&WATER, t) - 78.332_147_573_168).abs() < 1.0e-12,
            "water's eps_i = {}",
            component_dielectric(&WATER, t)
        );
        assert!((component_dielectric(&METHANE, t) - 2.0).abs() < 1.0e-12);
        assert!(
            (component_dielectric(&METHANOL, t) - 32.699_991_464_028_2).abs() < 1.0e-12,
            "methanol's eps_i = {}",
            component_dielectric(&METHANOL, t)
        );

        // The aqueous phase of the probe's first case: methane, water and two ions, with
        // the ion mole numbers what make the exclusion observable.
        let total = 1.001_901_653_436_76;
        let x = [
            0.000_225_745_660_581_355,
            0.997_778_050_427_449,
            0.000_998_101_955_985_164,
            0.000_998_101_955_985_164,
        ];
        let n: Vec<f64> = x.iter().map(|fraction| fraction * total).collect();
        let coefficients = [METHANE, WATER, METHANOL, WATER];
        let eps_i: Vec<f64> = coefficients
            .iter()
            .map(|c| component_dielectric(c, t))
            .collect();
        let ion = [false, false, true, true];

        let mixture = solvent_dielectric(MixingRule::MolarAverage, &n, &eps_i, &ion, &[])
            .expect("the solvents carry moles");
        assert!(
            (mixture - 78.314_881_455_398_7).abs() < 1.0e-12,
            "eps = {mixture}, the probe prints 78.3148814553987"
        );

        let dt: Vec<f64> = coefficients
            .iter()
            .map(|c| component_dielectric_dt(c, t))
            .collect();
        let mixture_dt = solvent_dielectric(MixingRule::MolarAverage, &n, &dt, &ion, &[])
            .expect("signed derivatives average the same way");
        assert!(
            (mixture_dt - (-0.359_218_709_298_880)).abs() < 1.0e-12,
            "deps/dT = {mixture_dt}"
        );

        let dtdt: Vec<f64> = coefficients
            .iter()
            .map(|c| component_dielectric_dtdt(c, t))
            .collect();
        let mixture_dtdt = solvent_dielectric(MixingRule::MolarAverage, &n, &dtdt, &ion, &[])
            .expect("and the second");
        assert!(
            (mixture_dtdt - 0.001_957_056_837_135_53).abs() < 1.0e-15,
            "d2eps/dT2 = {mixture_dtdt}"
        );
    }

    /// **The packing fractions, against the same capture.**
    ///
    /// The diameters are the probe's own `lj` column - which for an ion is the value
    /// *derived* from the fitted covolume, not the table's. The total is over every
    /// component and the ionic one over the ions, so this is where the two sums are told
    /// apart; the phase's volume is NeqSim's scaled `V` and the kernel takes m³/mol.
    #[test]
    fn the_packing_fractions_match_the_oracle() {
        let total = 1.001_901_653_436_76;
        // NeqSim returns `getMolarVolume()` in units of `1e5` m³/mol - `calcEps`'s own
        // `1e-5` is the conversion back - so the Si value is this over `1e5`.
        let molar_volume = 2.385_525_337_515_67e-5;
        let x = [
            0.000_225_745_660_581_355,
            0.997_778_050_427_449,
            0.000_998_101_955_985_164,
            0.000_998_101_955_985_164,
        ];
        let diameters: Vec<f64> = [2.52, 2.52, 4.343_717_662_170_81, 3.226_081_589_674_81]
            .iter()
            .map(|angstrom| angstrom * 1.0e-10)
            .collect();
        let n: Vec<f64> = x.iter().map(|fraction| fraction * total).collect();
        let ion = [false, false, true, true];

        let all =
            packing_fraction(&diameters, &n, &ion, total, molar_volume, false).expect("computes");
        assert!(
            (all - 0.212_659_929_859_501).abs() < 1.0e-13,
            "packing = {all}, the probe prints 0.212659929859501"
        );

        let ionic =
            packing_fraction(&diameters, &n, &ion, total, molar_volume, true).expect("computes");
        assert!(
            (ionic - 0.001_524_427_057_062_24).abs() < 1.0e-15,
            "packing_ionic = {ionic}, the probe prints 0.00152442705706224"
        );

        // `-eps/V` and `2 eps/V^2`, against the probe's `packing_V` and `packing_VV`.
        let volume = total * molar_volume;
        assert!(
            (packing_fraction_dv(all, total, molar_volume) + all / volume).abs() < 1.0e-15,
            "d eps/dV is -eps/V"
        );
        assert!(
            (packing_fraction_dvdv(all, total, molar_volume) - 2.0 * all / volume.powi(2)).abs()
                < 1.0e-9
        );
    }

    /// **The phase's dielectric constant, with the ions' volume taken out of it.**
    ///
    /// `calcDiElectricConstant` is not the solvent average the MSA and Born terms read,
    /// and the difference is the ionic correction - `78.1382247597158` against
    /// `78.3148814553987` at this state, which is 0.23% and would be invisible in a
    /// plausible-looking `ln phi`.
    #[test]
    fn the_phase_dielectric_takes_the_ions_volume_out() {
        let solvent = 78.314_881_455_398_7;
        let ionic = 0.001_524_427_057_062_24;
        let phase = phase_dielectric(solvent, ionic);
        assert!(
            (phase - 78.138_224_759_715_8).abs() < 1.0e-12,
            "eps_phase = {phase}, the probe prints 78.1382247597158"
        );
        assert!(
            phase < solvent,
            "a phase carrying ions has the smaller constant, and this says so"
        );
    }

    /// **A mixture with no solvent is refused rather than answered with a huge number.**
    ///
    /// NeqSim divides by its `1e-50` floor, which returns something of order `1e50`. That
    /// is not a dielectric constant, and a term built on it would be a plausible-looking
    /// number rather than a failure.
    #[test]
    fn a_mixture_with_no_solvent_is_refused() {
        let n = [0.5, 0.5];
        let eps = [0.0, 0.0];
        let ion = [true, true];
        let error = solvent_dielectric(MixingRule::MolarAverage, &n, &eps, &ion, &[])
            .expect_err("no solvent carries moles");
        assert_eq!(error.field(), Some("mole_numbers"), "{error:?}");
    }
}
