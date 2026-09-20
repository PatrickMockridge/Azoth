//! The Fürst short-range `Wij` table, and the `W` it sums to.
//!
//! `PhaseModifiedFurstElectrolyteEos`'s short-range term is `FSR2 = W/(V n (1 - eps))`, and
//! `W` is the double sum
//!
//! ```text
//! W(T) = -sum_ij n_i n_j Wij(T)
//! Wij(T) = w0 + w1 (1/T - 1/298.15) + w2 ((298.15 - T)/T + ln(T/298.15))
//! ```
//!
//! This module is the table `w0`, `w1`, `w2` per pair, built the way
//! `EosMixingRuleHandler.ElectrolyteMixRule.calcWij` builds it.
//!
//! # Three passes, and the order matters
//!
//! NeqSim fills the table in three sweeps over the components, each *overwriting* what the
//! last wrote:
//!
//! 1. **The table**, `W1`/`W2`/`W3`, per pair - read first and then mostly overwritten.
//! 2. **The cation sweep.** For every pair whose `CalcWij` is zero, a cation against an
//!    anion gets `furstParamsCPA[4] (d + d)^4 + [5]` and a cation against a neutral gets
//!    the neutral's own solvent correlation - `furstParamsCPA[2] d + [3]` for water, a
//!    per-solvent set for the glycols and alcohols, and a predictive form in the solvent's
//!    dielectric constant for anything else.
//! 3. **The gas sweep**, then the **oil sweep**. A non-polar gas against an ion takes
//!    `furstParamsGasIon`, and methanol/MEG/ethanol against an ion takes
//!    `furstParamsOIIon` - both *replacing* what the cation sweep computed.
//!
//! So a methane/Na+ pair is the gas-ion value and a water/Na+ pair is the correlation's,
//! and neither is visible from the other's table alone. The sweeps are also why `W` cannot
//! be read pairwise out of a file: it is a function of the whole composition.
//!
//! # What the diameters are
//!
//! The correlation's `d` is `getStokesCationicDiameter` for a cation and
//! `getPaulingAnionicDiameter` for an anion, and for an ion both are **the table's
//! `LJDIAMETER`** - 5.68 Å for Na+, 3.60 Å for Cl-. That is not the diameter the phase's
//! own terms use: `ComponentModifiedFurstElectrolyteEos`'s constructor overwrites
//! `lennardJonesMolecularDiameter` with a value derived back out of the fitted covolume, so
//! the MSA and Born terms see 4.3437 Å where this correlation sees 5.68 Å. Feeding this
//! module the derived diameter is a mistake that changes `W` by 44% and leaves no trace.

use azoth_core::{AzothError, Result};

use crate::databank::furst_parameters;

/// The reference temperature in `Wij(T)`, written into the formula rather than a constant.
const T_REFERENCE: f64 = 298.15;

/// One component, as the construction sees it.
#[derive(Debug, Clone, PartialEq)]
pub struct FurstComponent {
    /// The databank name, lower-cased. The sweeps test it - `water`, `meg`, `methanol` -
    /// so it is part of the arithmetic and not a label.
    pub name: String,
    /// The ionic charge, in elementary charges. Zero is neutral.
    pub charge: f64,
    /// The **table's** Lennard-Jones diameter in ångströms, which is the ionic diameter the
    /// correlation uses.
    pub diameter: f64,
    /// The component's dielectric constant at 298.15 K, which is what the predictive
    /// cation-solvent correlation reads for a solvent it does not name.
    pub dielectric_at_reference: f64,
}

impl FurstComponent {
    /// Whether this component is a cation for the construction's purposes.
    ///
    /// `charge > 0` exactly, as NeqSim writes it: the anion test beside it is
    /// `charge < -0.01`, so a component with a charge between them is treated as neutral by
    /// both branches - which is what a fractional charge would be.
    #[must_use]
    pub fn is_cation(&self) -> bool {
        self.charge > 0.0
    }

    /// Whether it is an anion, by NeqSim's own threshold.
    #[must_use]
    pub fn is_anion(&self) -> bool {
        self.charge < -0.01
    }
}

/// The pair table: `w0`, `w1` and `w2` for every ordered pair, flattened row-major.
#[derive(Debug, Clone, PartialEq)]
pub struct WijTable {
    /// The constant term, `N x N` row-major.
    pub w0: Vec<f64>,
    /// The `1/T - 1/298.15` coefficient.
    pub w1: Vec<f64>,
    /// The `(298.15 - T)/T + ln(T/298.15)` coefficient.
    pub w2: Vec<f64>,
    /// How many components.
    pub n: usize,
}

impl WijTable {
    /// `Wij(T)` for one pair.
    #[must_use]
    pub fn wij(&self, i: usize, j: usize, temperature: f64) -> f64 {
        let at = i * self.n + j;
        self.w0[at]
            + self.w1[at] * (1.0 / temperature - 1.0 / T_REFERENCE)
            + self.w2[at]
                * ((T_REFERENCE - temperature) / temperature + (temperature / T_REFERENCE).ln())
    }

    /// `dWij/dT` for one pair.
    #[must_use]
    pub fn wij_dt(&self, i: usize, j: usize, temperature: f64) -> f64 {
        let at = i * self.n + j;
        -self.w1[at] / (temperature * temperature)
            - self.w2[at] * (T_REFERENCE - temperature) / (temperature * temperature)
    }

    /// `d^2Wij/dT^2` for one pair.
    #[must_use]
    pub fn wij_dtdt(&self, i: usize, j: usize, temperature: f64) -> f64 {
        let at = i * self.n + j;
        2.0 * self.w1[at] / temperature.powi(3)
            + self.w2[at] / (temperature * temperature)
            + 2.0 * self.w2[at] * (T_REFERENCE - temperature) / temperature.powi(3)
    }
}

/// The three sums the short-range term reads, at a composition and temperature.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ShortRange {
    /// `W = -sum_ij n_i n_j Wij(T)`.
    pub w: f64,
    /// `dW/dT`.
    pub w_dt: f64,
    /// `d^2W/dT^2`.
    pub w_dtdt: f64,
}

/// Build the pair table for a phase, the way `ElectrolyteMixRule.calcWij` builds it.
///
/// `fitted` is the databank's `furst_wij` for a pair - `(W1, W2, W3, CalcWij)` - or `None`
/// for a pair with no row, which the handler reads as all three zero with the flag unset.
/// It is a closure rather than a slice because the caller holds the names.
///
/// # Errors
/// * [`AzothError::InvalidInput`] if `components` is empty.
pub fn wij_table(
    components: &[FurstComponent],
    fitted: impl Fn(&str, &str) -> Option<(f64, f64, f64, bool)>,
) -> Result<WijTable> {
    let n = components.len();
    if n == 0 {
        return Err(AzothError::invalid_input(
            "components",
            "the short-range table is a sum over pairs and there are no components",
        ));
    }
    let mut table = WijTable {
        w0: vec![0.0; n * n],
        w1: vec![0.0; n * n],
        w2: vec![0.0; n * n],
        n,
    };
    let mut calc_wij = vec![false; n * n];

    // Pass 1: the table's own values, symmetric.
    for i in 0..n {
        for j in 0..n {
            if let Some((w0, w1, w2, _is_fitted)) = fitted(&components[i].name, &components[j].name)
            {
                table.w0[i * n + j] = w0;
                table.w1[i * n + j] = w1;
                table.w2[i * n + j] = w2;
            }
        }
    }
    for i in 0..n {
        for j in 0..n {
            // The flag is read from the *pair*, and NeqSim's own read is symmetric because
            // the row is.
            calc_wij[i * n + j] = fitted(&components[i].name, &components[j].name)
                .map(|(_, _, _, is_fitted)| is_fitted)
                .unwrap_or(false);
        }
    }

    // Pass 2: the cation sweep.
    for i in 0..n {
        if !components[i].is_cation() {
            continue;
        }
        let divalent = components[i].charge.round() >= 2.0;
        for j in 0..n {
            if calc_wij[i * n + j] {
                continue;
            }
            if components[j].is_anion() {
                let (w0, w1, w2) =
                    cation_anion(components[i].diameter, components[j].diameter, divalent);
                set_symmetric(&mut table, i, j, w0, w1, w2);
            } else if components[j].charge == 0.0 {
                let w0 = cation_solvent(&components[i], &components[j], divalent);
                // **The temperature coefficients carry the diameter as well**, as the
                // `w0` ones do: `TDep[0] d + TDep[1]`. `getFurstParamTDep(0)` alone is
                // `0.005` where the pair's `w1` is `0.0164`, and the three-fold difference
                // is invisible in `W` at the reference temperature - it shows up only in
                // `dW/dT`, which is what found it.
                let d = components[i].diameter;
                let (w1, w2) = if components[i].name == "mdea+" {
                    (0.0, 0.0)
                } else if divalent {
                    (t_dep(8) * d + t_dep(9), t_dep(10) * d + t_dep(11))
                } else {
                    (t_dep(0) * d + t_dep(1), t_dep(2) * d + t_dep(3))
                };
                set_symmetric(&mut table, i, j, w0, w1, w2);
            }
        }
    }

    // Pass 3: the gas sweep, then the oil sweep.
    for i in 0..n {
        let gas = gas_indices(&components[i].name);
        if let Some((cation, anion)) = gas {
            for j in 0..n {
                if calc_wij[i * n + j] {
                    continue;
                }
                // **`w0` only.** NeqSim's gas sweep writes `wij[0]` and the symmetric
                // entry and touches neither `wij[1]` nor `wij[2]`, so the temperature
                // coefficients the cation sweep gave the pair **survive** it. Zeroing them
                // here leaves `dW/dT` short by the methane pairs' contribution - 8.3e-14
                // against 3.69e-10, which is 0.02% and reads as rounding.
                if components[j].is_cation() {
                    set_w0_symmetric(&mut table, i, j, gas_ion(cation));
                } else if components[j].is_anion() {
                    set_w0_symmetric(&mut table, i, j, gas_ion(anion));
                }
            }
        }
    }
    for i in 0..n {
        let Some((cation, anion)) = oil_indices(&components[i].name) else {
            continue;
        };
        for j in 0..n {
            if calc_wij[i * n + j] {
                continue;
            }
            // **Only a non-zero parameter is applied**, so a zero here leaves the cation
            // sweep's value in place rather than overwriting it with nothing.
            let index = if components[j].is_cation() {
                Some(cation)
            } else if components[j].is_anion() {
                Some(anion)
            } else {
                None
            };
            if let Some(index) = index {
                let value = named_from("furstParamsOIIon", index);
                if value.abs() > 1.0e-20 {
                    set_w0_symmetric(&mut table, i, j, value);
                }
            }
        }
    }

    Ok(table)
}

/// The three sums at a composition.
///
/// # Errors
/// * [`AzothError::InvalidInput`] if the mole numbers do not match the table.
pub fn short_range(table: &WijTable, mole_numbers: &[f64], temperature: f64) -> Result<ShortRange> {
    if mole_numbers.len() != table.n {
        return Err(AzothError::invalid_input(
            "mole_numbers",
            format!(
                "the table is {} x {} and there are {} mole numbers",
                table.n,
                table.n,
                mole_numbers.len()
            ),
        ));
    }
    let (mut w, mut w_dt, mut w_dtdt) = (0.0, 0.0, 0.0);
    for i in 0..table.n {
        for j in 0..table.n {
            let weight = mole_numbers[i] * mole_numbers[j];
            w += weight * table.wij(i, j, temperature);
            w_dt += weight * table.wij_dt(i, j, temperature);
            w_dtdt += weight * table.wij_dtdt(i, j, temperature);
        }
    }
    Ok(ShortRange {
        w: -w,
        w_dt: -w_dt,
        w_dtdt: -w_dtdt,
    })
}

/// `Wij(cation, anion) = p4 (d_c + d_a)^4 + p5`, on the monovalent or divalent set.
fn cation_anion(cation_diameter: f64, anion_diameter: f64, divalent: bool) -> (f64, f64, f64) {
    let sum4 = (cation_diameter + anion_diameter).powi(4);
    if divalent {
        (
            cpa(8) * sum4 + cpa(9),
            t_dep(12) * sum4 + t_dep(13),
            t_dep(14) * sum4 + t_dep(15),
        )
    } else {
        (
            cpa(4) * sum4 + cpa(5),
            t_dep(4) * sum4 + t_dep(5),
            t_dep(6) * sum4 + t_dep(7),
        )
    }
}

/// `Wij(cation, neutral)`, by the neutral's **name**, as `calcWij`'s chain tests it.
///
/// A solvent the chain does not name takes the predictive form, which is a linear
/// correction in `1/eps - 1/eps_water` to the water fit - so the chain is exhaustive and an
/// unrecognised solvent is a different number rather than no number.
fn cation_solvent(cation: &FurstComponent, solvent: &FurstComponent, divalent: bool) -> f64 {
    let d = cation.diameter;
    match solvent.name.as_str() {
        "water" => {
            if divalent {
                cpa(6) * d + cpa(7)
            } else {
                cpa(2) * d + cpa(3)
            }
        }
        // **`TEG` reads the `MEG` set**, in NeqSim and here: the chain tests `TEG` twice
        // and the first match wins, so the branch written for TEG is unreachable and
        // `furstParamsCPA_TEG` is never read by anything. Reproduced rather than corrected,
        // because a port that silently improved on its source would be a different model;
        // it is reported upstream and the fix is one array swap here.
        "teg" | "triethylene glycol" | "meg" | "ethylene glycol" => meg_fit(d, divalent),
        "mdea" => {
            if divalent {
                named(6)
            } else {
                named(2)
            }
        }
        "piperazine" => {
            if divalent {
                cpa(6) * d + cpa(7)
            } else {
                cpa(2) * d + cpa(3)
            }
        }
        "methanol" => {
            if divalent {
                named_from("furstParamsCPA_MeOH", 6)
            } else {
                named_from("furstParamsCPA_MeOH", 2)
            }
        }
        "ethanol" => {
            if divalent {
                named_from("furstParamsCPA_EtOH", 6)
            } else {
                named_from("furstParamsCPA_EtOH", 2)
            }
        }
        "mea" => {
            if divalent {
                named_from("furstParamsCPA_MEA", 6) * d + named_from("furstParamsCPA_MEA", 7)
            } else {
                named_from("furstParamsCPA_MEA", 2) * d + named_from("furstParamsCPA_MEA", 3)
            }
        }
        _ => predictive(solvent.dielectric_at_reference, d, divalent),
    }
}

/// The MEG fit, which is what a `TEG` component is given.
///
/// One function for two solvent names and not two functions, because that is the defect:
/// `furstParamsCPA_TEG` differs from this set at `[2]` (`4.98e-5` against `8.0e-5`), `[3]`,
/// `[6]` and `[7]`, and a `TEG` pair gets `2.11x` the `Wij` its own fit gives.
fn meg_fit(diameter: f64, divalent: bool) -> f64 {
    if divalent {
        named_from("furstParamsCPA_MEG", 6) * diameter + named_from("furstParamsCPA_MEG", 7)
    } else {
        named_from("furstParamsCPA_MEG", 2) * diameter + named_from("furstParamsCPA_MEG", 3)
    }
}

/// `getPredictiveWij`: a linear correction to the water fit in `1/eps - 1/eps_water`.
fn predictive(solvent_dielectric: f64, diameter: f64, divalent: bool) -> f64 {
    const EPSILON_WATER_REFERENCE: f64 = 78.4;
    let delta = 1.0 / solvent_dielectric - 1.0 / EPSILON_WATER_REFERENCE;
    let slope_water = if divalent { cpa(6) } else { cpa(2) };
    let intercept_water = if divalent { cpa(7) } else { cpa(3) };
    let slope_coefficient = if divalent { -2.5e-3 } else { -3.0e-3 };
    let intercept_coefficient = if divalent { 1.5e-2 } else { 1.2e-2 };
    (slope_water + slope_coefficient * delta) * diameter
        + (intercept_water + intercept_coefficient * delta)
}

/// The `furstParamsGasIon` index pair for a non-polar gas, by name.
fn gas_indices(name: &str) -> Option<(usize, usize)> {
    Some(match name {
        "co2" | "carbon dioxide" => (0, 1),
        "methane" | "ch4" => (2, 3),
        "ethane" | "c2h6" => (4, 5),
        "propane" | "c3h8" => (6, 7),
        "n-butane" | "i-butane" | "n-c4" | "i-c4" | "butane" => (8, 9),
        "n-pentane" | "i-pentane" | "n-c5" | "i-c5" | "pentane" | "n-hexane" | "hexane"
        | "n-c6" | "n-heptane" | "heptane" | "n-c7" | "n-octane" | "octane" | "n-c8"
        | "n-nonane" | "nonane" | "n-c9" | "n-decane" | "decane" | "n-c10" => (10, 11),
        "nitrogen" | "n2" => (12, 13),
        "h2s" | "hydrogen sulfide" | "hydrogensulfide" => (14, 15),
        "hydrogen" | "h2" => (16, 17),
        _ => return None,
    })
}

/// The `furstParamsOIIon` index pair for an organic inhibitor, by name.
fn oil_indices(name: &str) -> Option<(usize, usize)> {
    Some(match name {
        "methanol" | "meoh" => (0, 1),
        "meg" | "ethylene glycol" => (2, 3),
        "ethanol" | "etoh" => (4, 5),
        _ => return None,
    })
}

/// `furstParamsGasIon[i]`.
fn gas_ion(index: usize) -> f64 {
    named_from("furstParamsGasIon", index)
}

/// `furstParamsCPA[i]`, with a missing coefficient read as zero.
fn cpa(index: usize) -> f64 {
    named_from("furstParamsCPA", index)
}

/// `furstParamsCPA_TDep[i]`.
fn t_dep(index: usize) -> f64 {
    named_from("furstParamsCPA_TDep", index)
}

/// `furstParams[i]` - the base model's own set.
fn named(index: usize) -> f64 {
    named_from("furstParams", index)
}

/// One coefficient of one set, zero when the set or the index is missing.
///
/// A missing coefficient reads as zero rather than failing, which is NeqSim's own behaviour
/// for a table it cannot reach - and the sets this reads are vendored whole, so a missing
/// one is a defect in the table rather than a caller's mistake.
fn named_from(set: &str, index: usize) -> f64 {
    furst_parameters(set).get(index).copied().unwrap_or(0.0)
}

/// Write a pair's constant term both ways round and leave the temperature coefficients.
///
/// The gas and oil sweeps write `wij[0]` alone, so a pair the cation sweep has already given
/// a `w1` and `w2` keeps them.
fn set_w0_symmetric(table: &mut WijTable, i: usize, j: usize, w0: f64) {
    let n = table.n;
    table.w0[i * n + j] = w0;
    table.w0[j * n + i] = w0;
}

/// Write a pair both ways round, as the handler does after every assignment.
fn set_symmetric(table: &mut WijTable, i: usize, j: usize, w0: f64, w1: f64, w2: f64) {
    let n = table.n;
    table.w0[i * n + j] = w0;
    table.w0[j * n + i] = w0;
    table.w1[i * n + j] = w1;
    table.w1[j * n + i] = w1;
    table.w2[i * n + j] = w2;
    table.w2[j * n + i] = w2;
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::databank;

    /// The aqueous phase of `SystemFurstElectrolyteEosTest`'s own mixture: methane, water,
    /// Na+ and Cl-. The diameters are the **table's** 5.68 and 3.60 Å and not the 4.3437 and
    /// 3.2261 the phase's own terms use - see the module note.
    fn components() -> Vec<FurstComponent> {
        vec![
            FurstComponent {
                name: "methane".into(),
                charge: 0.0,
                diameter: 2.52,
                dielectric_at_reference: 2.0,
            },
            FurstComponent {
                name: "water".into(),
                charge: 0.0,
                diameter: 2.52,
                dielectric_at_reference: 78.332_147_573_168_0,
            },
            FurstComponent {
                name: "na+".into(),
                charge: 1.0,
                diameter: 5.68,
                dielectric_at_reference: 0.0,
            },
            FurstComponent {
                name: "cl-".into(),
                charge: -1.0,
                diameter: 3.60,
                dielectric_at_reference: 0.0,
            },
        ]
    }

    fn table() -> WijTable {
        wij_table(&components(), databank::furst_wij).expect("builds")
    }

    /// **The pair table, pair by pair.**
    ///
    /// Every value here comes from the measured decomposition of the probe's `W`, so the
    /// four contributions can be read from the table rather than inferred from the sum. The
    /// two ion pairs are the correlation's and the two methane pairs are the gas sweep's -
    /// the methane/Na+ entry is *not* what the cation sweep computed before it, which is what
    /// says the second pass runs.
    #[test]
    fn the_pair_table_matches_the_measured_decomposition() {
        let table = table();
        let t = 298.15;
        // At the reference temperature the `w1` and `w2` terms vanish, so `Wij` is `w0`.
        let at = |i: usize, j: usize| table.wij(i, j, t);

        assert!(
            (at(2, 1) - 0.000_162_975_563_535_845_85).abs() < 1.0e-18,
            "Na+/water = {}, and `cpa[2] d + cpa[3]` gives 0.00016297556353584585",
            at(2, 1)
        );
        assert!(
            (at(2, 3) - (-0.000_248_050_280_932_059_94)).abs() < 1.0e-18,
            "Na+/Cl- = {}, and `cpa[4] (d+d)^4 + cpa[5]` gives -0.00024805028093205994",
            at(2, 3)
        );
        assert!(
            (at(0, 2) - 1.05e-04).abs() < 1.0e-18,
            "methane/Na+ = {}, and furstParamsGasIon[2] is 1.05e-04 - not the -4.0616e-04 \
             the predictive correlation gives, which is what the pair held before the gas \
             sweep overwrote it",
            at(0, 2)
        );
        assert!((at(0, 3) - 1.05e-04).abs() < 1.0e-18, "methane/Cl-");
    }

    /// **The three sums, against the probe.**
    ///
    /// `W = -3.25444241835884e-07`, `WT = 3.68791163569515e-10` and
    /// `WTT = -3.18892643563736e-12`. The derivatives are checked at the reference
    /// temperature, where `w0` alone does not determine them - `w1` and `w2` do - so this is
    /// the test that says the temperature coefficients were read from `furstParamsCPA_TDep`
    /// rather than left at their table values of zero.
    #[test]
    fn the_sums_match_the_oracle() {
        let table = table();
        let moles = [
            0.000_225_745_660_581_355 * 1.001_901_653_436_76,
            0.997_778_050_427_449 * 1.001_901_653_436_76,
            0.000_998_101_955_985_164 * 1.001_901_653_436_76,
            0.000_998_101_955_985_164 * 1.001_901_653_436_76,
        ];
        let sums = short_range(&table, &moles, 298.15).expect("computes");
        assert!(
            (sums.w - (-3.254_442_418_358_84e-07)).abs() < 1.0e-21,
            "W = {}, the probe prints -3.25444241835884e-07",
            sums.w
        );
        assert!(
            (sums.w_dt - 3.687_911_635_695_15e-10).abs() < 1.0e-24,
            "dW/dT = {}, the probe prints 3.68791163569515e-10",
            sums.w_dt
        );
        assert!(
            (sums.w_dtdt - (-3.188_926_435_637_36e-12)).abs() < 1.0e-26,
            "d2W/dT2 = {}, the probe prints -3.18892643563736e-12",
            sums.w_dtdt
        );
        assert!(
            sums.w_dt != 0.0,
            "a zero here would mean the temperature coefficients were never read, and `w0` \
             alone reproduces `W` at the reference temperature"
        );
    }

    /// **The ion diameter this correlation uses is the table's, and using the other one
    /// moves `W` by 44%.**
    ///
    /// The phase's own MSA and Born terms see the diameter derived from the fitted covolume
    /// - 4.3437 Å for Na+ against the table's 5.68. Feeding that here gives
    ///   `-1.9224146500595774e-07` against the oracle's `-3.25444241835884e-07`, so this
    ///   is the sabotage check the model's cases are read against.
    #[test]
    fn the_derived_diameter_is_the_wrong_one_for_this_correlation() {
        let mut wrong = components();
        wrong[2].diameter = 4.343_717_662_170_81;
        wrong[3].diameter = 3.226_081_589_674_81;
        let table = wij_table(&wrong, databank::furst_wij).expect("builds");
        let moles = [
            0.000_225_745_660_581_355 * 1.001_901_653_436_76,
            0.997_778_050_427_449 * 1.001_901_653_436_76,
            0.000_998_101_955_985_164 * 1.001_901_653_436_76,
            0.000_998_101_955_985_164 * 1.001_901_653_436_76,
        ];
        let sums = short_range(&table, &moles, 298.15).expect("computes");
        assert!(
            (sums.w - (-1.922_414_650_059_577_4e-07)).abs() < 1.0e-21,
            "with the derived diameters W = {}, and the measured wrong value is \
             -1.9224146500595774e-07",
            sums.w
        );
        assert!(
            (sums.w - (-3.254_442_418_358_84e-07)).abs() > 1.0e-07,
            "the two diameters must not agree, which is the whole point"
        );
    }

    /// **A fitted pair is left alone by the correlation.**
    ///
    /// `CalcWij` is `1` on six amine rows, and there the table's `W1` is what the phase uses,
    /// so the sweep skips the pair. Asserted on `MDEA+`/`methane`, whose `W1` is `1.23e-4`
    /// and whose gas-sweep value would be `1.05e-4`: the two differ, so which one comes out
    /// says which path ran.
    #[test]
    fn a_fitted_pair_keeps_its_own_value() {
        let components = vec![
            FurstComponent {
                name: "methane".into(),
                charge: 0.0,
                diameter: 2.52,
                dielectric_at_reference: 2.0,
            },
            FurstComponent {
                name: "mdea+".into(),
                charge: 1.0,
                diameter: 5.68,
                dielectric_at_reference: 0.0,
            },
        ];
        let table = wij_table(&components, databank::furst_wij).expect("builds");
        assert!(
            (table.wij(0, 1, 298.15) - 1.23e-04).abs() < 1.0e-18,
            "the fitted W1 is 1.23e-04 and the pair holds {}",
            table.wij(0, 1, 298.15)
        );
    }

    /// **A mixture of components that are none of the named kinds still gets a number.**
    ///
    /// The cation sweep's last arm is a predictive correlation in the solvent's dielectric
    /// constant, so an unnamed solvent is a different value rather than a zero - and a
    /// correlation is not a refusal.
    #[test]
    fn an_unnamed_solvent_takes_the_predictive_form() {
        let components = vec![
            FurstComponent {
                name: "toluene".into(),
                charge: 0.0,
                diameter: 5.0,
                dielectric_at_reference: 2.38,
            },
            FurstComponent {
                name: "na+".into(),
                charge: 1.0,
                diameter: 5.68,
                dielectric_at_reference: 0.0,
            },
        ];
        let table = wij_table(&components, databank::furst_wij).expect("builds");
        let value = table.wij(1, 0, 298.15);
        assert!(
            value != 0.0 && value.is_finite(),
            "the predictive form gives {value}"
        );
        // And it is not the water value, which the same diameter would give if the solvent's
        // dielectric constant were ignored.
        assert!(
            (value - (cpa(2) * 5.68 + cpa(3))).abs() > 1.0e-8,
            "toluene's constant is 2.38 and water's reference is 78.4, so the two must differ"
        );
    }
}
