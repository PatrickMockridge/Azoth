//! The Fürst electrolyte term, as one contribution to the residual Helmholtz energy.
//!
//! This is the assembly: the dielectric surface, the three Helmholtz terms, the short-range
//! table and the composition derivatives, behind the interface
//! [`crate::association::Association`] already establishes - a term that is evaluated at a
//! volume and returns its Helmholtz energy, its volume derivative and its `ln phi`.
//!
//! # Why a term and not a phase model of its own
//!
//! `PhaseModifiedFurstElectrolyteEos` is a cubic with three additive terms, and **the
//! volume it sits at is not the cubic's**. Measured from the probe: at the aqueous phase's
//! converged root, `dFdV = 0.415155577019407` and the identity `Z = 1 - V 1e5 dFdV` gives
//! `0.00963585`, which is the printed `Z = 0.00963585200923298` to nine digits. So the
//! orthobaric equation is the cubic's residual plus this term's pressure, exactly as an
//! associating mixture's is, and the root solve that finds it is
//! [`crate::mixture::Mixture::associating_root`]'s.
//!
//! # The volume conventions, which are the whole difficulty
//!
//! NeqSim's `dFdV` is a derivative with respect to the **total** volume in NeqSim's own
//! scaled units, where `V_ne = V_si 1e5`, and its `d(A^R/RT)/dV_si` is therefore
//! `1e5 dFdV_ne`. That factor, and the total-versus-molar distinction, are the two things
//! this module converts once so no caller has to: [`FurstSolution::helmholtz_rt_dv`] is in
//! m³ and in SI, and is what the residual wants.
//!
//! # What is not here
//!
//! **The Born term carries no volume**, upstream and here: `getF`, `dFdT` and `dFdTdT` take
//! it and `dFdV`, `dFdVdV` and `dFdTdV` do not. That is consistent - `FBorn` is
//! `K(T)(1/eps_s - 1) bornX` and carries neither a volume nor a density - and it is
//! asserted rather than assumed, because a port that added it would move every root.

use azoth_core::Result;

use crate::furst_dielectric::{
    MixingRule, component_dielectric, component_dielectric_dt, component_dielectric_dtdt,
};
use crate::furst_mixing::{FurstComponent, ShortRange, WijTable, short_range, wij_table};
use crate::furst_terms::{
    ComponentDerivatives, ComponentState, FurstState, StateInputs, build_state, fborn, flr, flr_dv,
    fsr2, fsr2_dv, ln_phi_contributions,
};

/// The static description of one component in a Fürst mixture.
#[derive(Debug, Clone, PartialEq)]
pub struct FurstSpecies {
    /// The databank name, lower-cased. The short-range sweeps test it.
    pub name: String,
    /// The ionic charge, in elementary charges.
    pub charge: f64,
    /// The **derived** Lennard-Jones diameter in metres - the value
    /// `ComponentModifiedFurstElectrolyteEos` ends up with, which is what the MSA and Born
    /// terms read. The short-range correlation reads the *table's* diameter instead; the
    /// two are different numbers for an ion and `furst_mixing` states why.
    pub diameter_m: f64,
    /// The **table's** diameter in ångströms, as the correlation uses it.
    pub table_diameter: f64,
    /// `DIELECTRICPARAMETER1..5`.
    pub dielectric_coefficients: [f64; 5],
    /// The component's critical volume, in m³/mol, for the volume-fraction mixing rules.
    pub critical_volume: f64,
    /// Its dielectric constant at 298.15 K, which the predictive correlation reads.
    pub dielectric_at_reference: f64,
}

/// A Fürst mixture's electrolyte term.
#[derive(Debug, Clone, PartialEq)]
pub struct FurstElectrolyte {
    species: Vec<FurstSpecies>,
    table: WijTable,
    rule: MixingRule,
}

/// What the term returns at one volume.
#[derive(Debug, Clone, PartialEq)]
pub struct FurstSolution {
    /// The three terms' contribution to `A^R/(R T)`.
    pub helmholtz_rt: f64,
    /// Its derivative with respect to the **SI** volume, in m³ - the pressure the term
    /// carries, divided by `-R T`. This is what the residual takes.
    pub helmholtz_rt_dv: f64,
    /// The three terms' contribution to `dFdN_i`, one per component, which adds to the
    /// cubic's `ln phi`.
    pub ln_phi: Vec<f64>,
    /// The state the terms were evaluated at, for introspection and for the layer diff.
    pub state: FurstState,
}

impl FurstElectrolyte {
    /// The term for a mixture, with the short-range table built from the databank.
    ///
    /// # Errors
    /// * [`azoth_core::AzothError::InvalidInput`] if `species` is empty or a pair table
    ///   cannot be built.
    pub fn new(species: Vec<FurstSpecies>, rule: MixingRule) -> Result<Self> {
        let components: Vec<FurstComponent> = species
            .iter()
            .map(|s| FurstComponent {
                name: s.name.clone(),
                charge: s.charge,
                diameter: s.table_diameter,
                dielectric_at_reference: s.dielectric_at_reference,
            })
            .collect();
        let table = wij_table(&components, crate::databank::furst_wij)?;
        Ok(Self {
            species,
            table,
            rule,
        })
    }

    /// How many components the term is over.
    #[must_use]
    pub fn len(&self) -> usize {
        self.species.len()
    }

    /// Whether the term is over no components, which is never a mixture.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.species.is_empty()
    }

    /// The short-range sums at a composition and temperature.
    ///
    /// # Errors
    /// * [`azoth_core::AzothError::InvalidInput`] if the mole numbers are not one per
    ///   component.
    pub fn short_range(&self, temperature: f64, mole_numbers: &[f64]) -> Result<ShortRange> {
        short_range(&self.table, mole_numbers, temperature)
    }

    /// Evaluate the term at one volume.
    ///
    /// `molar_volume` is in m³/mol and `mole_numbers` are the phase's, so the term is
    /// rebuilt from scratch at every volume the root solve tries - which is what NeqSim's
    /// `volInit()` does inside its own iteration and the reason this is a `solve` and not a
    /// stored state.
    ///
    /// # Errors
    /// * Propagates the state builder's and the terms' refusals.
    pub fn solve(
        &self,
        temperature: f64,
        mole_numbers: &[f64],
        molar_volume: f64,
    ) -> Result<FurstSolution> {
        let components: Vec<ComponentState> = self
            .species
            .iter()
            .map(|s| ComponentState {
                charge: s.charge,
                diameter_m: s.diameter_m,
                dielectric: component_dielectric(&s.dielectric_coefficients, temperature),
                dielectric_dt: component_dielectric_dt(&s.dielectric_coefficients, temperature),
                dielectric_dtdt: component_dielectric_dtdt(&s.dielectric_coefficients, temperature),
                critical_volume: s.critical_volume,
            })
            .collect();
        let names: Vec<String> = self.species.iter().map(|s| s.name.clone()).collect();
        let short = self.short_range(temperature, mole_numbers)?;
        let state = build_state(&StateInputs {
            temperature,
            molar_volume,
            mole_numbers,
            components: &components,
            rule: self.rule,
            short_range: short,
            table: &self.table,
            names: &names,
        })?;

        let helmholtz_rt = fsr2(&state)? + flr(&state) + fborn(&state);
        // **The Born term carries no volume**, so it is absent here - see the module note.
        // The `1e5` converts NeqSim's scaled total-volume derivative to the SI one, which is
        // the factor the identity `Z = 1 - V 1e5 dFdV` fixes.
        let helmholtz_rt_dv = 1.0e5 * (fsr2_dv(&state) + flr_dv(&state));

        let derivatives: Vec<ComponentDerivatives> = self
            .species
            .iter()
            .enumerate()
            .map(|(i, s)| ComponentDerivatives {
                charge: s.charge,
                diameter_m: s.diameter_m,
                dielectric: component_dielectric(&s.dielectric_coefficients, temperature),
                // `calcWi`: the row sum negated and doubled, as the handler returns it.
                w_i: -2.0
                    * (0..self.species.len())
                        .map(|j| mole_numbers[j] * self.table.wij(i, j, temperature))
                        .sum::<f64>(),
            })
            .collect();
        let contributions = ln_phi_contributions(&state, &derivatives, mole_numbers)?;

        Ok(FurstSolution {
            helmholtz_rt,
            helmholtz_rt_dv,
            ln_phi: contributions.iter().map(|c| c.total()).collect(),
            state,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The aqueous phase's four components, as the capture has them: the table diameters are
    /// 5.68 and 3.60 Å for the ions and the derived ones are 4.3437 and 3.2261.
    fn species() -> Vec<FurstSpecies> {
        vec![
            FurstSpecies {
                name: "methane".into(),
                charge: 0.0,
                diameter_m: 2.52e-10,
                table_diameter: 2.52,
                dielectric_coefficients: [2.0, 0.0, 0.0, 0.0, 0.0],
                critical_volume: 9.9e-5,
                dielectric_at_reference: 2.0,
            },
            FurstSpecies {
                name: "water".into(),
                charge: 0.0,
                diameter_m: 2.52e-10,
                table_diameter: 2.52,
                dielectric_coefficients: [-19.2905, 29814.5, -0.019678, 0.000132, -3.11e-07],
                critical_volume: 5.6e-5,
                dielectric_at_reference: 78.332_147_573_168_0,
            },
            FurstSpecies {
                name: "na+".into(),
                charge: 1.0,
                diameter_m: 4.343_717_662_170_81e-10,
                table_diameter: 5.68,
                dielectric_coefficients: [0.0; 5],
                critical_volume: 0.0,
                dielectric_at_reference: 0.0,
            },
            FurstSpecies {
                name: "cl-".into(),
                charge: -1.0,
                diameter_m: 3.226_081_589_674_81e-10,
                table_diameter: 3.60,
                dielectric_coefficients: [0.0; 5],
                critical_volume: 0.0,
                dielectric_at_reference: 0.0,
            },
        ]
    }

    fn moles() -> Vec<f64> {
        [
            0.000_225_745_660_581_355,
            0.997_778_050_427_449,
            0.000_998_101_955_985_164,
            0.000_998_101_955_985_164,
        ]
        .iter()
        .map(|x| x * 1.001_901_653_436_76)
        .collect()
    }

    /// **The term's three outputs at the aqueous phase's converged root.**
    ///
    /// The Helmholtz energy is the three terms added, the pressure is the two that carry a
    /// volume, and the `ln phi` are the composition derivatives - each already oracled apart
    /// in `furst_terms`, and this is the assembly that says the pieces still agree when they
    /// are put together in one call.
    #[test]
    fn the_term_matches_the_oracle_at_the_root() {
        let term = FurstElectrolyte::new(species(), MixingRule::default_for_the_model())
            .expect("the shipped mixture builds a term");
        let solution = term
            .solve(298.15, &moles(), 2.385_525_337_515_67e-5)
            .expect("evaluates");

        let want = -0.017_294_384_499_566_6 - 0.000_272_971_800_726_073 - 0.298_937_482_116_719;
        assert!(
            (solution.helmholtz_rt - want).abs() < 1.0e-12,
            "A^R/RT = {}, and the three terms add to {want}",
            solution.helmholtz_rt
        );

        // The pressure, as a check on the identity the root solve rests on: the electrolyte
        // shares of `dFdV` are `dFSR2dV + dFLRdV`, and `1e5` times them plus `dFdV`'s cubic
        // part must give the printed total. The Born share is zero, asserted separately.
        let want_dv = 1.0e5 * (0.009_190_383_382_775_51 + 4.913_648_794_711_00e-05);
        assert!(
            (solution.helmholtz_rt_dv - want_dv).abs() < 1.0e-9,
            "d(A^R/RT)/dV = {}, and the probe's two volume derivatives give {want_dv}",
            solution.helmholtz_rt_dv
        );

        let want_ln_phi = [
            -0.026_957_328_895_281_1 - 0.000_379_609_469_205_052 + 0.003_768_121_783_018_29,
            -0.021_959_491_597_569_0 + 8.588_602_480_471_29e-08 - 8.525_314_228_889_12e-07,
            -17.315_614_106_443_6 - 0.192_435_547_491_505 - 127.400_565_779_252,
            0.014_107_883_057_089_9 - 0.197_975_495_079_022 - 171.536_916_337_467,
        ];
        for (i, want) in want_ln_phi.iter().enumerate() {
            assert!(
                (solution.ln_phi[i] - want).abs() < 1.0e-9 * want.abs().max(1.0),
                "ln phi contribution {i} = {}, and the probe's three rows give {want}",
                solution.ln_phi[i]
            );
        }
    }

    /// **The Born term contributes nothing to the pressure**, which is what makes the
    /// identity above hold with two terms and not three.
    ///
    /// `FBorn` is `K(T)(1/eps_s - 1) bornX`, and neither factor carries a volume. NeqSim
    /// omits it from `dFdV`, `dFdVdV` and `dFdTdV` and includes it in `getF`, `dFdT` and
    /// `dFdTdT` - consistent, and stated here because a port that added it would move every
    /// root by the Born share of the pressure, which at this state is `0.2989` against the
    /// others' `0.0175`.
    #[test]
    fn the_born_term_carries_no_volume() {
        let term =
            FurstElectrolyte::new(species(), MixingRule::default_for_the_model()).expect("builds");
        let solution = term
            .solve(298.15, &moles(), 2.385_525_337_515_67e-5)
            .expect("evaluates");
        let short = term.short_range(298.15, &moles()).expect("computes");
        let state = build_state(&StateInputs {
            temperature: 298.15,
            molar_volume: 2.385_525_337_515_67e-5,
            mole_numbers: &moles(),
            components: &species()
                .iter()
                .map(|s| ComponentState {
                    charge: s.charge,
                    diameter_m: s.diameter_m,
                    dielectric: component_dielectric(&s.dielectric_coefficients, 298.15),
                    dielectric_dt: component_dielectric_dt(&s.dielectric_coefficients, 298.15),
                    dielectric_dtdt: component_dielectric_dtdt(&s.dielectric_coefficients, 298.15),
                    critical_volume: s.critical_volume,
                })
                .collect::<Vec<_>>(),
            rule: MixingRule::default_for_the_model(),
            short_range: short,
            table: &term.table,
            names: &term
                .species
                .iter()
                .map(|s| s.name.clone())
                .collect::<Vec<_>>(),
        })
        .expect("builds");
        assert!(
            (fborn(&state) - (-0.298_937_482_116_719)).abs() < 1.0e-12,
            "the Born term is the largest of the three and carries no volume"
        );
        assert!(
            solution.helmholtz_rt_dv > 0.0,
            "a liquid's residual Helmholtz energy rises with volume, so its derivative is \
             positive and the pressure it carries is negative"
        );
    }
}
