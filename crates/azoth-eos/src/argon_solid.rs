//! The face-centered-cubic solid argon Helmholtz equation of Maltby, Hammer and
//! Wilhelmsen (2024), as NeqSim's `ArgonSolidHelmholtzEquation` carries it.
//!
//! A three-body-corrected Buckingham static lattice evaluated by the coordination-sphere
//! method, a zero-point term, Debye and Einstein vibrational contributions, and the
//! reported anharmonic corrections. The stated range is temperatures to 300 K and
//! pressures to 16 GPa.

use std::sync::OnceLock;

use crate::dual::{Dual, SolidState, state_from_helmholtz};

/// The gas constant of the published equation, J/(mol.K).
pub const R: f64 = 8.314_462_618;

const AVOGADRO: f64 = 6.022_140_76e23;
const BOLTZMANN: f64 = R / AVOGADRO;

const EPSILON: f64 = BOLTZMANN * 134.7;
const REPULSIVE_EXPONENT: f64 = 14.19;
const MINIMUM_POTENTIAL_DISTANCE: f64 = 3.802e-10;
const THREE_BODY_COEFFICIENT: f64 = BOLTZMANN * 3.202e5 * 1.0e-90;

const A_D: f64 = 10.4;
const B_D: f64 = 0.025_03;
const C_D: f64 = 0.000_156_8;
const V0: f64 = 22.56e-6;
const THETA_D0: f64 = 92.0;
const GAMMA_D0: f64 = 2.563;
const Q_D: f64 = 0.2874;

const B1: f64 = -0.000_447_5;
const B2: f64 = 2.041e-6;
const B3: f64 = 5.75e-7;
const C1: f64 = 0.7204;
const C2: f64 = -1.614;
const C3: f64 = -0.019_43;
const C4: f64 = -27.64;

const EINSTEIN_WEIGHTS: [f64; 3] = [0.0261, 0.037_84, 0.045_12];
const EINSTEIN_TEMPERATURES: [f64; 3] = [77.81, 550.0, 45.36];
const EINSTEIN_GRUNEISEN: [f64; 3] = [6.221, 1.617e-6, 3.1278];

const Z1: f64 = 140.0;
const Z2: f64 = 2.34;
const Z3: f64 = 0.683;
const Z4: f64 = 19.7e-6;

const MAXIMUM_TEMPERATURE: f64 = 300.0;
const MAXIMUM_PRESSURE: f64 = 160_000.0e5;
const MINIMUM_VOLUME: f64 = 1.0e-7;
const MAXIMUM_VOLUME: f64 = 1.0e-3;

const COORDINATION_SHELL_COUNT: usize = 20;
const FCC_ENUMERATION_LIMIT: i32 = 24;

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

    let static_lattice = static_lattice_energy(volume);
    let zero_point = ((Z2 / Z3) * (1.0 - (reduced_volume * (V0 / Z4)).powf(Z3))).exp() * (Z1 * R);

    let theta_d_temperature =
        THETA_D0 + A_D * ((-B_D * temp.powf(2.0) - C_D * temp.powf(3.0)).exp() - 1.0);
    let theta_d = ((GAMMA_D0 / Q_D) * (1.0 - reduced_volume.powf(Q_D))).exp() * theta_d_temperature;
    let weight_sum: f64 = EINSTEIN_WEIGHTS.iter().sum();
    let debye = (theta_d / temp).debye_free_energy() * temp * (3.0 * R * (1.0 - weight_sum));

    let mut einstein = Dual::constant(0.0);
    for i in 0..3 {
        let theta =
            (EINSTEIN_GRUNEISEN[i] * (1.0 - reduced_volume)).exp() * EINSTEIN_TEMPERATURES[i];
        einstein = einstein + (theta / temp).log_one_minus_exp_negative() * EINSTEIN_WEIGHTS[i];
    }
    let einstein = einstein * temp * (3.0 * R);

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

    static_lattice + zero_point + debye + einstein + anharmonic_temperature + anharmonic_cold
}

/// The published finite FCC Buckingham static-lattice sum, per mole.
fn static_lattice_energy(volume: Dual) -> Dual {
    let lattice_constant = (volume * (4.0 / AVOGADRO)).powf(1.0 / 3.0);
    let three_body_factor = 1.0
        - volume.reciprocal()
            * (THREE_BODY_COEFFICIENT * AVOGADRO / (EPSILON * buckingham_zero_distance().powi(6)));

    let (distances, counts) = fcc_shells();
    let mut shell_energy = Dual::constant(0.0);
    for i in 0..COORDINATION_SHELL_COUNT {
        let distance = lattice_constant * (0.5 * (distances[i] as f64).sqrt());
        let effective_potential = buckingham_pair_potential(distance) * three_body_factor;
        shell_energy = shell_energy + effective_potential * (counts[i] as f64);
    }
    shell_energy * (0.5 * AVOGADRO)
}

/// The Buckingham exp-6 pair potential.
fn buckingham_pair_potential(distance: Dual) -> Dual {
    let repulsive =
        (distance * (-REPULSIVE_EXPONENT / MINIMUM_POTENTIAL_DISTANCE) + REPULSIVE_EXPONENT).exp()
            * (EPSILON * 6.0 / (REPULSIVE_EXPONENT - 6.0));
    let attractive = (distance.reciprocal() * MINIMUM_POTENTIAL_DISTANCE).powf(6.0)
        * (EPSILON * REPULSIVE_EXPONENT / (REPULSIVE_EXPONENT - 6.0));
    repulsive - attractive
}

/// The scalar Buckingham pair potential, for the zero-distance search.
fn buckingham_pair_potential_scalar(distance: f64) -> f64 {
    let repulsive = EPSILON * 6.0 / (REPULSIVE_EXPONENT - 6.0)
        * (REPULSIVE_EXPONENT * (1.0 - distance / MINIMUM_POTENTIAL_DISTANCE)).exp();
    let attractive = EPSILON * REPULSIVE_EXPONENT / (REPULSIVE_EXPONENT - 6.0)
        * (MINIMUM_POTENTIAL_DISTANCE / distance).powi(6);
    repulsive - attractive
}

/// The physically relevant zero of the Buckingham potential, by bisection.
fn buckingham_zero_distance() -> f64 {
    static ZERO: OnceLock<f64> = OnceLock::new();
    *ZERO.get_or_init(|| {
        let mut lower = 0.5 * MINIMUM_POTENTIAL_DISTANCE;
        let mut upper = MINIMUM_POTENTIAL_DISTANCE;
        for _ in 0..100 {
            let middle = 0.5 * (lower + upper);
            if buckingham_pair_potential_scalar(middle) > 0.0 {
                lower = middle;
            } else {
                upper = middle;
            }
        }
        0.5 * (lower + upper)
    })
}

/// The first 20 FCC coordination shells: `(squared distance, coordination number)`.
fn fcc_shells() -> &'static (
    [i32; COORDINATION_SHELL_COUNT],
    [i32; COORDINATION_SHELL_COUNT],
) {
    static SHELLS: OnceLock<(
        [i32; COORDINATION_SHELL_COUNT],
        [i32; COORDINATION_SHELL_COUNT],
    )> = OnceLock::new();
    SHELLS.get_or_init(|| {
        let mut counts: Vec<(i32, i32)> = Vec::new();
        let mut by_distance: std::collections::BTreeMap<i32, i32> =
            std::collections::BTreeMap::new();
        for h in -FCC_ENUMERATION_LIMIT..=FCC_ENUMERATION_LIMIT {
            for k in -FCC_ENUMERATION_LIMIT..=FCC_ENUMERATION_LIMIT {
                for l in -FCC_ENUMERATION_LIMIT..=FCC_ENUMERATION_LIMIT {
                    let squared = h * h + k * k + l * l;
                    if squared == 0 || (h + k + l) % 2 != 0 {
                        continue;
                    }
                    *by_distance.entry(squared).or_insert(0) += 1;
                }
            }
        }
        for (distance, count) in by_distance {
            if counts.len() == COORDINATION_SHELL_COUNT {
                break;
            }
            counts.push((distance, count));
        }
        let mut distances = [0; COORDINATION_SHELL_COUNT];
        let mut coordination = [0; COORDINATION_SHELL_COUNT];
        for (i, (d, c)) in counts.iter().enumerate() {
            distances[i] = *d;
            coordination[i] = *c;
        }
        (distances, coordination)
    })
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

    let mut lower = V0.ln();
    let mut upper = lower;
    let mut lower_residual = pressure_residual(temperature, lower, pressure_pa);
    let mut upper_residual = lower_residual;
    let expansion = 1.02_f64.ln();
    if lower_residual > 0.0 {
        let mut iteration = 0;
        while upper_residual > 0.0 && iteration < 600 {
            lower = upper;
            lower_residual = upper_residual;
            upper = MAXIMUM_VOLUME.ln().min(upper + expansion);
            upper_residual = pressure_residual(temperature, upper, pressure_pa);
            iteration += 1;
        }
    } else {
        let mut iteration = 0;
        while lower_residual < 0.0 && iteration < 600 {
            upper = lower;
            upper_residual = lower_residual;
            lower = MINIMUM_VOLUME.ln().max(lower - expansion);
            lower_residual = pressure_residual(temperature, lower, pressure_pa);
            iteration += 1;
        }
    }
    assert!(
        lower_residual * upper_residual < 0.0,
        "could not bracket the argon volume root"
    );

    let mut current = 0.5 * (lower + upper);
    for _ in 0..100 {
        let volume = current.exp();
        let helmholtz = helmholtz(temperature, volume);
        let residual = -helmholtz.dv - pressure_pa;
        if residual.abs() / pressure_pa.max(1000.0) < 1.0e-9 {
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
        if upper - lower < 1.0e-14 {
            return (0.5 * (lower + upper)).exp();
        }
    }
    panic!("argon volume solver did not converge");
}

/// The pressure residual at a logarithmic molar volume.
fn pressure_residual(temperature: f64, log_volume: f64, pressure_pa: f64) -> f64 {
    pressure(temperature, log_volume.exp()) - pressure_pa
}
