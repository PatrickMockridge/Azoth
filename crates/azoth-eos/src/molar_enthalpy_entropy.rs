//! `eos.molar_enthalpy_entropy` - the absolute enthalpy and entropy of a mixture.
//!
//! ```text
//! H = H_ig(T_ref) + integral Cp dT      + H_dep
//! S = S_ig(T_ref) + integral Cp/T dT    - R ln(P/P_ref) - R sum z_i ln z_i + S_dep
//! ```
//!
//! Spec: `specs/models/eos/molar_enthalpy_entropy.yaml`, which carries the assembly -
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
use crate::results::MolarEnthalpyEntropyResult;
use crate::{ideal_gas_cp::REFERENCE_TEMPERATURE, model_gen};

/// The caller's ideal-gas model and the datum it is referenced to.
///
/// One struct rather than eight arguments because the eight belong together: a
/// coefficient set without a reference state is not a thermodynamic model, and passing
/// them separately invites a call that supplies four of the six vectors.
#[derive(Debug, Clone, PartialEq)]
pub struct IdealGasModel {
    /// The constant term of each component's `Cp/R` polynomial.
    pub cp_a: Vec<f64>,
    /// The coefficient of `theta` in each component's polynomial.
    pub cp_b: Vec<f64>,
    /// The coefficient of `theta**2`.
    pub cp_c: Vec<f64>,
    /// The coefficient of `theta**3`.
    pub cp_d: Vec<f64>,
    /// Each component's ideal-gas molar enthalpy at [`Self::t_ref`].
    pub h_ref: Vec<f64>,
    /// Each component's ideal-gas molar entropy at [`Self::t_ref`] and [`Self::p_ref`].
    pub s_ref: Vec<f64>,
    /// The temperature the reference values are given at.
    pub t_ref: ThermodynamicTemperature,
    /// The pressure `s_ref` is given at. It does not enter the enthalpy.
    pub p_ref: Pressure,
}

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
            ("h_ref", &self.h_ref),
            ("s_ref", &self.s_ref),
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
        for (field, value) in [("T_ref", self.t_ref.value), ("P_ref", self.p_ref.value)] {
            // Explicit finiteness-and-positivity rather than `!(value > 0.0)`: both
            // reject NaN, and the explicit form also rejects an infinity, which
            // would make `ln(T/T_ref)` infinite.
            if !value.is_finite() || value <= 0.0 {
                return Err(AzothError::OutOfRange {
                    field: field.to_string(),
                    value,
                    detail: "the reference state is a state, and both of its coordinates \
                             are positive; `ln(T/T_ref)` and `ln(P/P_ref)` make neither \
                             zero nor negative usable"
                        .to_string(),
                });
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
            "T_ref" => Some(ideal_gas.t_ref.value),
            "P_ref" => Some(ideal_gas.p_ref.value),
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

    // The ideal-gas part: each component's reference value plus the exact integral of
    // its polynomial between the reference temperature and the state. The integrals
    // are closed forms rather than quadrature, because the polynomial is a polynomial.
    let theta = t.value / REFERENCE_TEMPERATURE;
    let theta_ref = ideal_gas.t_ref.value / REFERENCE_TEMPERATURE;
    let mut h_ideal = 0.0;
    let mut s_ideal = 0.0;
    let mut cp_ideal_over_r = 0.0;
    for (i, &z_i) in z.iter().enumerate() {
        let (a, b, c, d) = (
            ideal_gas.cp_a[i],
            ideal_gas.cp_b[i],
            ideal_gas.cp_c[i],
            ideal_gas.cp_d[i],
        );
        // `T = REFERENCE_TEMPERATURE * theta`, so `dT = REFERENCE_TEMPERATURE dtheta`.
        let dh = MOLAR_GAS_CONSTANT
            * REFERENCE_TEMPERATURE
            * (a * (theta - theta_ref)
                + b * (theta.powi(2) - theta_ref.powi(2)) / 2.0
                + c * (theta.powi(3) - theta_ref.powi(3)) / 3.0
                + d * (theta.powi(4) - theta_ref.powi(4)) / 4.0);
        // `integral Cp/T dT` is `R * integral (a + b theta + ...)/theta dtheta`.
        let ds = MOLAR_GAS_CONSTANT
            * (a * (theta / theta_ref).ln()
                + b * (theta - theta_ref)
                + c * (theta.powi(2) - theta_ref.powi(2)) / 2.0
                + d * (theta.powi(3) - theta_ref.powi(3)) / 3.0);
        h_ideal += z_i * (ideal_gas.h_ref[i] + dh);
        s_ideal += z_i * (ideal_gas.s_ref[i] + ds);
        // The polynomial itself, at the state's `theta`. It is the integrand of `dh`,
        // so reporting it costs one evaluation and no new assumption.
        cp_ideal_over_r += z_i * (a + b * theta + c * theta.powi(2) + d * theta.powi(3));
    }
    // The two ideal-gas terms that no coefficient switches off. The pressure one is
    // per mole of mixture - an ideal gas's entropy falls by `R ln(P/P_ref)` however
    // many components it has - and the mixing one is `-R sum z_i ln z_i`.
    s_ideal -= MOLAR_GAS_CONSTANT * (p.value / ideal_gas.p_ref.value).ln();
    s_ideal -= MOLAR_GAS_CONSTANT * z.iter().map(|&zi| zi * zi.ln()).sum::<f64>();

    let h_departure = MOLAR_GAS_CONSTANT * t.value * state.h_dep_rt;
    let s_departure = MOLAR_GAS_CONSTANT * state.s_dep_r;

    let cp_ideal = MOLAR_GAS_CONSTANT * cp_ideal_over_r;
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
