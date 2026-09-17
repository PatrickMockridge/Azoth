//! The hcp phase-I solid para-hydrogen Helmholtz equation of Sannerhaugen (2026), as
//! NeqSim's `ParaHydrogenSolidHelmholtzEquation` carries it.
//!
//! A Vinet cold curve, Debye and three Einstein external modes, an internal-vibration
//! mode, and the reported anharmonic terms. The stated range is temperatures below 200 K
//! and pressures below 10 GPa.

use crate::dual::{Dual, SolidState, state_from_helmholtz};

/// The gas constant of the thesis equation, J/(mol.K).
pub const R: f64 = 8.314_462_1;

const TRIPLE_POINT_TEMPERATURE: f64 = 13.8033;

const MAXIMUM_TEMPERATURE: f64 = 200.0;
const MAXIMUM_PRESSURE: f64 = 100_000.0e5;
const MINIMUM_VOLUME: f64 = 1.0e-7;
const MAXIMUM_VOLUME: f64 = 1.0e-3;

const V0: f64 = 23.14e-6;
const B0: f64 = 180.5e6;
const B0_PRIME: f64 = 7.131;

const THETA_D0: f64 = 127.0;
const GAMMA_D0: f64 = 2.843;
const Q_D: f64 = 1.647;
const A_D: f64 = 0.000_873_4;
const B_D: f64 = 0.9295;
const C_D: f64 = 1.020;

const EINSTEIN_WEIGHTS: [f64; 3] = [0.096_16, 0.1597, 0.004_324];
const EINSTEIN_TEMPERATURES: [f64; 3] = [51.62, 247.0, 187.9];
const EINSTEIN_GRUNEISEN: [f64; 3] = [2.290, 0.9952, 1.367];

const THETA_INTERNAL0: f64 = 1086.0;
const GAMMA_INTERNAL: f64 = 2.775;
const EXTERNAL_MODES_PER_MOLECULE: f64 = 5.0;

const B1: f64 = 0.083_37;
const B2: f64 = 22.97;
const B3: f64 = -5.5412;
const C1: f64 = -1.704;
const C2: f64 = -0.081_74;
const C3: f64 = -0.039_86;
const C4: f64 = 3.828;

/// The pressure at a temperature and molar volume, from `-a_V`.
#[must_use]
pub fn pressure(temperature: f64, molar_volume: f64) -> f64 {
    -helmholtz(temperature, molar_volume).dv
}

/// The Helmholtz energy at a temperature and molar volume.
#[must_use]
pub fn helmholtz_energy(temperature: f64, molar_volume: f64) -> f64 {
    helmholtz(temperature, molar_volume).value
}

/// The solid state at a pressure (Pa) and temperature (K).
///
/// # Panics
/// If the state is out of the published range or the volume root cannot be bracketed.
#[must_use]
pub fn properties(temperature: f64, pressure_pa: f64) -> SolidState {
    let molar_volume = solve_volume(temperature, pressure_pa);
    state_from_helmholtz(
        temperature,
        pressure_pa,
        molar_volume,
        helmholtz(temperature, molar_volume),
    )
}

/// The solid state at a temperature and molar volume, at the implied pressure.
#[must_use]
pub fn properties_at_volume(temperature: f64, molar_volume: f64) -> SolidState {
    let pressure_pa = pressure(temperature, molar_volume);
    state_from_helmholtz(
        temperature,
        pressure_pa,
        molar_volume,
        helmholtz(temperature, molar_volume),
    )
}

/// The second-order Helmholtz derivatives at a temperature and molar volume.
#[must_use]
pub fn helmholtz(temperature: f64, molar_volume: f64) -> Dual {
    let temp = Dual::temperature(temperature);
    let volume = Dual::volume(molar_volume);
    let reduced_volume = volume / V0;

    let vinet = vinet_energy(reduced_volume);

    let theta_d_temperature =
        THETA_D0 + A_D * ((-B_D * temp.powf(2.0) - C_D * temp.powf(3.0)).exp() - 1.0);
    let theta_d = ((GAMMA_D0 / Q_D) * (1.0 - reduced_volume.powf(Q_D))).exp() * theta_d_temperature;
    let weight_sum: f64 = EINSTEIN_WEIGHTS.iter().sum();
    let debye = (theta_d / temp).debye_free_energy()
        * temp
        * (R * EXTERNAL_MODES_PER_MOLECULE * (1.0 - weight_sum));

    let mut einstein = Dual::constant(0.0);
    for i in 0..3 {
        let theta =
            (EINSTEIN_GRUNEISEN[i] * (1.0 - reduced_volume)).exp() * EINSTEIN_TEMPERATURES[i];
        einstein = einstein + (theta / temp).log_one_minus_exp_negative() * EINSTEIN_WEIGHTS[i];
    }
    let einstein = einstein * temp * (R * EXTERNAL_MODES_PER_MOLECULE);

    let theta_internal = (GAMMA_INTERNAL * (1.0 - reduced_volume)).exp() * THETA_INTERNAL0;
    let internal = (theta_internal / temp).log_one_minus_exp_negative() * temp * R;

    let temperature_ratio = temp / THETA_D0;
    let anharmonic_temperature = temperature_ratio.powf(4.0)
        / (temperature_ratio.powf(2.0) * B2 + 1.0)
        * theta_d_temperature
        * B1
        * (B3 * (reduced_volume - 1.0)).exp()
        * R;
    let anharmonic_cold = (C1 * (C2 * (1.0 - reduced_volume)).exp()
        + reduced_volume.reciprocal() * C3 * (C4 * (1.0 - reduced_volume)).exp())
        * R;

    vinet + debye + einstein + internal + anharmonic_temperature + anharmonic_cold
}

/// The canonical Vinet cold-curve Helmholtz energy.
fn vinet_energy(reduced_volume: Dual) -> Dual {
    let x = reduced_volume.powf(1.0 / 3.0);
    let exponent_coefficient = 1.5 * (B0_PRIME - 1.0);
    let vinet_exponent = (Dual::constant(1.0) - x) * exponent_coefficient;
    ((vinet_exponent - 1.0) * vinet_exponent.exp() + 1.0)
        * (4.0 * B0 * V0 / (B0_PRIME - 1.0).powi(2))
}

/// Solve `-a_V = pressure_pa` for the molar volume by bracketing and bisection/Newton in
/// logarithmic volume.
fn solve_volume(temperature: f64, pressure_pa: f64) -> f64 {
    assert!(
        pressure_pa > 0.0 && pressure_pa <= MAXIMUM_PRESSURE,
        "pressure out of range"
    );
    assert!(
        temperature > 0.0 && temperature <= MAXIMUM_TEMPERATURE,
        "temperature out of range"
    );

    let mut guessed_volume = V0;
    if temperature >= TRIPLE_POINT_TEMPERATURE {
        let melting_volume = (27.1788 - 0.283_044 * temperature) * 1.0e-6;
        if melting_volume > MINIMUM_VOLUME {
            guessed_volume = melting_volume;
        }
    }

    let mut lower = guessed_volume.ln();
    let mut upper = lower;
    let mut lower_residual = pressure_residual(temperature, lower, pressure_pa);
    let mut upper_residual = lower_residual;
    let expansion = 1.1_f64.ln();
    for _ in 0..160 {
        if lower_residual * upper_residual < 0.0 {
            break;
        }
        lower = MINIMUM_VOLUME.ln().max(lower - expansion);
        upper = MAXIMUM_VOLUME.ln().min(upper + expansion);
        lower_residual = pressure_residual(temperature, lower, pressure_pa);
        upper_residual = pressure_residual(temperature, upper, pressure_pa);
    }
    assert!(
        lower_residual * upper_residual < 0.0,
        "could not bracket the para-hydrogen volume root"
    );

    let mut current = upper.min(lower.max(guessed_volume.ln()));
    for _ in 0..100 {
        let volume = current.exp();
        let helmholtz = helmholtz(temperature, volume);
        let residual = -helmholtz.dv - pressure_pa;
        if residual.abs() / pressure_pa.max(1000.0) < 1.0e-10 {
            return volume;
        }

        let derivative = -helmholtz.dvv * volume;
        let mut candidate = current - residual / derivative;
        if derivative >= 0.0 || !candidate.is_finite() || candidate <= lower || candidate >= upper {
            candidate = 0.5 * (lower + upper);
        }
        let candidate_residual = pressure_residual(temperature, candidate, pressure_pa);
        if lower_residual * candidate_residual <= 0.0 {
            upper = candidate;
        } else {
            lower = candidate;
            lower_residual = candidate_residual;
        }
        current = candidate;
    }
    panic!("para-hydrogen volume solver did not converge");
}

/// The pressure residual at a logarithmic molar volume.
fn pressure_residual(temperature: f64, log_volume: f64, pressure_pa: f64) -> f64 {
    pressure(temperature, log_volume.exp()) - pressure_pa
}
