//! `eos.molar_enthalpy_entropy` - the absolute enthalpy and entropy of a mixture.
//!
//! ```text
//! H = H_ig(T_ref) + integral Cp dT      + H_dep
//! S = S_ig(T_ref) + integral Cp/T dT    - R ln(P/P_ref) - R sum z_i ln z_i + S_dep
//! ```
//!
//! Spec: `specs/models/eos/molar_enthalpy_entropy.toml`, which carries the assembly -
//! which term belongs on which side - and what the caller is responsible for: the
//! datum, the coefficients and the compressibility factor.
//!
//! The departure functions are [`crate::mixture`]'s and the integrals are exact
//! integrals of the polynomial [`crate::ideal_gas_cp`] registers.

use azoth_core::units::{
    Pressure, ThermodynamicTemperature, joules_per_mole, joules_per_mole_kelvin,
};
use azoth_core::{AzothError, Result, apply_checks};

use crate::mixture::Mixture;
use crate::pr_molar_volume::MOLAR_GAS_CONSTANT;

use crate::model_gen;
use crate::results::MolarEnthalpyEntropyResult;

/// The ideal-gas heat-capacity coefficients of every component in the mixture.
///
/// Five vectors rather than eight arguments because they belong together: a coefficient
/// set without the other four is not a polynomial, and passing them separately invites a
/// call that supplies three of the five.
///
/// **There is no datum here, and that is NeqSim's design rather than an omission.** Its
/// `getHID` is `integral Cp dT` from a fixed `referenceTemperature` of 273.15 K, and it
/// multiplies the formation enthalpy by zero; its entropy is the same integral of
/// `Cp/T` from the same temperature, less `R ln(P / referencePressure)`. So the ideal-gas
/// state is fully determined by the five coefficients, and `h_ref`, `s_ref`, `T_ref` and
/// `P_ref` were azoth's inventions.
#[derive(Debug, Clone, PartialEq)]
pub struct IdealGasModel {
    /// The constant term of each component's `Cp`, in J/(mol*K).
    pub cp_a: Vec<f64>,
    /// The coefficient of `T`, in J/(mol*K**2).
    pub cp_b: Vec<f64>,
    /// The coefficient of `T**2`, in J/(mol*K**3).
    pub cp_c: Vec<f64>,
    /// The coefficient of `T**3`, in J/(mol*K**4).
    pub cp_d: Vec<f64>,
    /// The coefficient of `T**4`, in J/(mol*K**5).
    pub cp_e: Vec<f64>,
}

/// The temperature NeqSim's ideal-gas integrals are measured from, in kelvin.
///
/// `ThermodynamicConstantsInterface.referenceTemperature`. Fixed rather than a caller's
/// choice: it is where `getHID` and `getIdEntropy` both start, so it is part of the
/// correlation rather than a datum a caller supplies.
pub const REFERENCE_TEMPERATURE: f64 = 273.15;

/// The pressure NeqSim's ideal-gas entropy is measured from, in pascals.
///
/// `ThermodynamicConstantsInterface.referencePressure`, which is 1.01325 bar.
pub const REFERENCE_PRESSURE: f64 = 1.01325e5;

impl IdealGasModel {
    /// Check the six vectors agree in length and the reference state is a state.
    ///
    /// # Errors
    /// * [`AzothError::InvalidInput`] if any vector's length differs from `n`.
    /// * [`AzothError::OutOfRange`] if `t_ref` or `p_ref` is not positive.
    fn validate(&self, n: usize) -> Result<()> {
        for (name, vector) in [
            ("cp_a", &self.cp_a),
            ("cp_b", &self.cp_b),
            ("cp_c", &self.cp_c),
            ("cp_d", &self.cp_d),
            ("cp_e", &self.cp_e),
        ] {
            if vector.len() != n {
                return Err(AzothError::invalid_input(
                    name,
                    format!(
                        "a mixture of {n} components needs one value per component, but \
                         {name} has {}",
                        vector.len()
                    ),
                ));
            }
        }
        Ok(())
    }
}

/// The absolute molar enthalpy and entropy of a mixture at a state.
///
/// `compressibility` is the cubic's root for the phase wanted, and the model does not
/// check that it describes that phase - see the spec's assumptions.
///
/// # Errors
/// * [`AzothError::OutOfRange`] if `T`, `P` or either reference coordinate is not
///   positive, or if `compressibility` is not admissible at the mixture's `B`.
/// * [`AzothError::InvalidInput`] if the composition or any of the ideal-gas vectors
///   is the wrong length, or if `z` is not a composition.
pub fn molar_enthalpy_entropy(
    mixture: &Mixture,
    ideal_gas: &IdealGasModel,
    t: ThermodynamicTemperature,
    p: Pressure,
    z: &[f64],
    compressibility: f64,
) -> Result<MolarEnthalpyEntropyResult> {
    let spec = &model_gen::MOLAR_ENTHALPY_ENTROPY_SPEC;
    let mut warnings = Vec::new();

    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "T" => Some(t.value),
            "P" => Some(p.value),
            "compressibility" => Some(compressibility),
            _ => None,
        },
        &mut warnings,
    )?;

    let n = mixture.len();
    ideal_gas.validate(n)?;
    check_composition(z, n)?;

    let reduced = mixture.reduced_parameters(t, p)?;
    warnings.extend(reduced.warnings.iter().cloned());
    let state = mixture.phase_state_at(&reduced, z, compressibility)?;

    // The ideal-gas part, which is NeqSim's `getHID` and `getIdEntropy`: each
    // component's polynomial integrated exactly from the fixed reference temperature to
    // the state. Closed forms rather than quadrature, because the polynomial is a
    // polynomial - and a five-term one, so the enthalpy gains a fifth-power term and
    // the entropy a fourth.
    let (t_ref, p_ref) = (REFERENCE_TEMPERATURE, REFERENCE_PRESSURE);
    let t2 = t.value * t.value;
    let t3 = t2 * t.value;
    let t4 = t3 * t.value;
    let t5 = t4 * t.value;
    let r2 = t_ref * t_ref;
    let r3 = r2 * t_ref;
    let r4 = r3 * t_ref;
    let r5 = r4 * t_ref;

    let mut h_ideal = 0.0;
    let mut s_ideal = 0.0;
    let mut cp_ideal = 0.0;
    for (i, &z_i) in z.iter().enumerate() {
        let (a, b, c, d, e) = (
            ideal_gas.cp_a[i],
            ideal_gas.cp_b[i],
            ideal_gas.cp_c[i],
            ideal_gas.cp_d[i],
            ideal_gas.cp_e[i],
        );
        // `integral Cp dT` from `T_ref` to `T`.
        h_ideal += z_i
            * (a * (t.value - t_ref)
                + b * (t2 - r2) / 2.0
                + c * (t3 - r3) / 3.0
                + d * (t4 - r4) / 4.0
                + e * (t5 - r5) / 5.0);
        // `integral Cp / T dT`, the same limits.
        s_ideal += z_i
            * (a * (t.value / t_ref).ln()
                + b * (t.value - t_ref)
                + c * (t2 - r2) / 2.0
                + d * (t3 - r3) / 3.0
                + e * (t4 - r4) / 4.0);
        // The polynomial itself at the state. It is the integrand of the enthalpy, so
        // reporting it costs one evaluation and no new assumption - and it is the same
        // `eos.ideal_gas_cp` a caller may evaluate directly.
        cp_ideal += z_i * (a + b * t.value + c * t2 + d * t3 + e * t4);
    }
    // The two ideal-gas terms that no coefficient switches off. The pressure one is per
    // mole of mixture - an ideal gas's entropy falls by `R ln(P/P_ref)` however many
    // components it has - and the mixing one is `-R sum z_i ln z_i`.
    s_ideal -= MOLAR_GAS_CONSTANT * (p.value / p_ref).ln();
    // **`0 ln 0 = 0`**, the convention the ideal-mixing term is defined with: a component the
    // mixture does not contain contributes nothing to its entropy, and `0.0 * (-inf)` is a
    // NaN that would poison the whole sum. Exposed by the first model that routes a component
    // entirely away - `process.component_splitter`, whose outlets carry exact zeros.
    s_ideal -= MOLAR_GAS_CONSTANT
        * z.iter()
            .filter(|&&zi| zi > 0.0)
            .map(|&zi| zi * zi.ln())
            .sum::<f64>();

    let h_departure = MOLAR_GAS_CONSTANT * t.value * state.h_dep_rt;
    let s_departure = MOLAR_GAS_CONSTANT * state.s_dep_r;

    let cp_departure = MOLAR_GAS_CONSTANT * state.cp_dep_r;

    Ok(MolarEnthalpyEntropyResult {
        h: joules_per_mole(h_ideal + h_departure),
        s: joules_per_mole_kelvin(s_ideal + s_departure),
        h_ideal: joules_per_mole(h_ideal),
        s_ideal: joules_per_mole_kelvin(s_ideal),
        h_departure: joules_per_mole(h_departure),
        s_departure: joules_per_mole_kelvin(s_departure),
        psi_bar: state.psi_bar,
        cp: joules_per_mole_kelvin(cp_ideal + cp_departure),
        cp_ideal: joules_per_mole_kelvin(cp_ideal),
        cp_departure: joules_per_mole_kelvin(cp_departure),
        warnings,
    })
}

/// A composition must be one entry per component, in `[0, 1]`, summing to one.
fn check_composition(values: &[f64], n: usize) -> Result<()> {
    if values.len() != n {
        return Err(AzothError::invalid_input(
            "z",
            format!(
                "a mixture of {n} components needs {n} mole fractions, but z has {}",
                values.len()
            ),
        ));
    }
    if let Some(bad) = values.iter().position(|&value| value < 0.0) {
        return Err(AzothError::invalid_input(
            "z",
            format!(
                "z[{bad}] is {} but a mole fraction cannot be negative",
                values[bad]
            ),
        ));
    }
    let sum: f64 = values.iter().sum();
    if (sum - 1.0).abs() > 1.0e-09 {
        return Err(AzothError::invalid_input(
            "z",
            format!(
                "the composition sums to {sum}, not to one. Renormalising it here would \
                 make a caller's error invisible in every number downstream, so it is \
                 refused instead"
            ),
        ));
    }
    Ok(())
}
