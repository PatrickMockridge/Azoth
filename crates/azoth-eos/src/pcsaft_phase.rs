//! PC-SAFT's volume at a temperature, a pressure and a composition.
//!
//! ```text
//! solve P_calc(v) = P   for   v,   Z = P v/(R T)
//! ```
//!
//! The layers the pressure is built from are [`crate::pcsaft`]'s; this module is the
//! solve that turns them into a state, so the layers stay checkable one at a time.
//!
//! # The two roots
//!
//! A PC-SAFT isotherm is not monotonic: below the critical temperature the pressure dips
//! and rises again, so a temperature and a pressure can admit three volumes. Which one
//! describes the phase wanted is the caller's statement, not the fluid's, and it is the
//! same `RootSide` the cubic roots are asked for by.

use azoth_core::units::{Pressure, ThermodynamicTemperature, cubic_meters_per_mole};
use azoth_core::{AzothError, Result};

use crate::association::R;
use crate::mixture::RootSide;
use crate::pcsaft::{self, PcsaftComponent};
use crate::results::PcsaftPhaseResult;

/// A solved PC-SAFT volume.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PcsaftVolume {
    /// The molar volume, in m^3/mol.
    pub v: f64,
    /// The compressibility factor at that volume.
    pub z: f64,
    /// The packing fraction at that volume, which is what the hard-sphere terms are
    /// functions of and is below one at any state.
    pub eta: f64,
    /// How many Newton steps the vapour branch took, or zero on the liquid branch, whose
    /// root is bracketed and bisected. Reported so a caller can see a solve that did not
    /// converge rather than only its output.
    pub iterations: usize,
}

/// Solve `P_calc(v) = p` for the molar volume, in m^3/mol.
///
/// NeqSim's `PhasePCSAFTRahmat.calcVolume` is a damped Newton on the volume: `h = P -
/// P_calc`, `dh = -dP_calc/dv`, and `v += 0.9 h/dh`, stopping when the step is a relative
/// `1e-10` or after a hundred of them. That is the loop here, in SI and in the volume
/// itself rather than in NeqSim's reduced density.
///
/// The packing fraction reaches one at `v = (pi/6) N_A md3`, where every hard-sphere term
/// diverges. That volume is the floor: below it there is no state, and a solve that would
/// step under it refuses instead of answering from outside the domain.
///
/// # Errors
/// * [`AzothError::OutOfRange`] if `t` or `p` is not positive, or if the wanted branch has
///   no zero at this state.
/// * [`AzothError::InvalidInput`] as [`crate::pcsaft::state`].
pub fn molar_volume(
    components: &[PcsaftComponent],
    kij: &[f64],
    x: &[f64],
    t: f64,
    p: f64,
    side: RootSide,
) -> Result<PcsaftVolume> {
    if !p.is_finite() || p <= 0.0 {
        return Err(AzothError::OutOfRange {
            field: "P".to_string(),
            value: p,
            detail: "the volume solve divides by the pressure".to_string(),
        });
    }
    let rt = R * t;
    let ideal = rt / p;
    let md3 = pcsaft::state(components, kij, x, t, ideal)?.md3;
    let floor = std::f64::consts::PI / 6.0 * pcsaft::AVOGADRO * md3;

    // `h(v) = P_calc(v) - p` in Pa. `state` refuses a packing fraction of one, so a
    // volume outside the domain fails here rather than returning a large number.
    let h = |v: f64| -> Result<f64> {
        Ok(pcsaft::pressure_over_rt(components, kij, x, t, v)? * rt - p)
    };
    let slope = |v: f64| -> Result<f64> {
        Ok(pcsaft::d_pressure_over_rt_dv(components, kij, x, t, v)? * rt)
    };

    // A bisection on a fixed count rather than on the residual: near the floor `h` is a
    // difference of large terms, and one that stopped at an absolute tolerance there
    // would stop early.
    fn bisect(h: &impl Fn(f64) -> Result<f64>, lo: f64, hi: f64) -> Result<f64> {
        let (mut lo, mut hi) = (lo, hi);
        let sign = h(lo)?.signum();
        for _ in 0..200 {
            let mid = 0.5 * (lo + hi);
            if h(mid)?.signum() == sign {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        Ok(0.5 * (lo + hi))
    }

    let (v, iterations) = match side {
        RootSide::Vapour => {
            let mut v = ideal;
            let mut converged = None;
            for step in 0..100 {
                // `v += 0.9 (P - P_calc)/(dP_calc/dv)`, NeqSim's step: `h` is
                // `P_calc - P`, so the sign is Newton's with the residual negated.
                let delta = -0.9 * h(v)? / slope(v)?;
                if !delta.is_finite() || v + delta <= floor {
                    break;
                }
                let relative = delta.abs() / v;
                v += delta;
                if relative < 1.0e-10 {
                    converged = Some(step + 1);
                    break;
                }
            }
            match converged {
                Some(steps) => (v, steps),
                // Newton left the domain, or a step was not a number, or it stopped
                // moving. The vapour root is the **last** zero in volume, so the walk is
                // taken for the last bracket rather than the first.
                None => {
                    let (lo, hi) = walk(&h, floor, ideal)?
                        .last()
                        .copied()
                        .ok_or_else(|| no_root("vapour", t, p))?;
                    (bisect(&h, lo, hi)?, 0)
                }
            }
        }
        RootSide::Liquid => {
            // **The lowest zero above the floor, and Newton from the ideal gas cannot find
            // it.** At a pressure where the isotherm has one root the vapour one is it, so
            // a solve seeded dilute converges there and stops - the liquid root is a
            // separate zero the seed never points at. The walk is geometric because the
            // root's *ratio* to the floor is what is bounded, not its distance.
            let (lo, hi) = walk(&h, floor, ideal)?
                .first()
                .copied()
                .ok_or_else(|| no_root("liquid", t, p))?;
            (bisect(&h, lo, hi)?, 0)
        }
    };

    let eta = std::f64::consts::PI / 6.0 * pcsaft::AVOGADRO * md3 / v;
    Ok(PcsaftVolume {
        v,
        z: p * v / rt,
        eta,
        iterations,
    })
}

/// The sign-change brackets of `h` between the packing floor and a volume far enough into
/// ideality that `h` can only be negative, in increasing volume.
///
/// `h` is `+inf` at the floor, so the walk starts a hair above it, and at the top
/// `P_calc v/(RT)` is within `1e-4` of one - below the pressure - so at least one change
/// is always present.
fn walk(h: &impl Fn(f64) -> Result<f64>, floor: f64, ideal: f64) -> Result<Vec<(f64, f64)>> {
    let start = floor * (1.0 + 1.0e-9);
    let top = 1.0e4 * ideal;
    let steps = 128;
    let mut previous = start;
    let mut previous_value = h(previous)?;
    let mut out = Vec::new();
    for step in 1..=steps {
        let v = start * (top / start).powf(step as f64 / steps as f64);
        let value = h(v)?;
        if previous_value.signum() != value.signum() {
            out.push((previous, v));
        }
        previous = v;
        previous_value = value;
    }
    Ok(out)
}

/// The refusal for a branch whose root does not exist at a state, which is a real answer:
/// this fluid has no such phase here.
fn no_root(branch: &str, t: f64, p: f64) -> AzothError {
    AzothError::OutOfRange {
        field: "P".to_string(),
        value: p,
        detail: format!(
            "the {branch} branch has no zero at {t} K and {p} Pa. The isotherm at this \
             temperature does not have a root on that side, so there is no such phase \
             rather than a volume to report"
        ),
    }
}

/// A name list resolved to PC-SAFT parameters and the `KIJPCSAFT` matrix.
///
/// # Errors
/// * [`AzothError::InvalidInput`] if a name is not in the databank, or the component has no
///   PC-SAFT set - which the table spells as zeros in all three columns, 47 of its 286
///   rows, and a model that took them would answer for a fluid with no segments.
pub fn parameters_of(names: &[&str]) -> Result<(Vec<PcsaftComponent>, Vec<f64>)> {
    let mut components = Vec::with_capacity(names.len());
    for name in names {
        let entry = crate::databank::entry(name, None)?;
        components.push(PcsaftComponent {
            m: entry.m_saft,
            sigma: entry.sigma_saft,
            epsik: entry.epsik_saft,
        });
    }
    Ok((components, crate::databank::pcsaft_kij(names)))
}

/// `eos.pcsaft_phase`: the phase state PC-SAFT gives at a temperature and a pressure.
///
/// The names cross unresolved and are looked up here, so the two-kernel comparison covers
/// the resolution as well as the arithmetic - which for this model means the `KIJPCSAFT`
/// column and the zero-means-absent convention on the parameters.
///
/// # Errors
/// * [`AzothError::InvalidInput`] as [`parameters_of`], or if `z` is not one entry per
///   component or does not sum to one.
/// * [`AzothError::OutOfRange`] if `T` or `P` is not positive, or the wanted branch has no
///   root at this state.
pub fn pcsaft_phase(
    components: &[String],
    t: ThermodynamicTemperature,
    p: Pressure,
    z: &[f64],
    compressed_phase: &str,
) -> Result<PcsaftPhaseResult> {
    let spec = &crate::model_gen::PCSAFT_PHASE_SPEC;
    let mut warnings = Vec::new();
    azoth_core::range::apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "T" => Some(t.value),
            "P" => Some(p.value),
            _ => None,
        },
        &mut warnings,
    )?;

    let names: Vec<&str> = components.iter().map(String::as_str).collect();
    let (parameters, kij) = parameters_of(&names)?;
    let side = side_of(compressed_phase)?;
    let mut result = pcsaft_phase_of(&parameters, &kij, t, p, z, side)?;
    warnings.extend(result.warnings);
    result.warnings = warnings;
    Ok(result)
}

/// The side a case's `compressed_phase` names.
///
/// # Errors
/// [`AzothError::InvalidInput`] for anything but the two names the spec declares, rather
/// than a default: a caller that mistyped a branch would otherwise be handed a root.
pub fn side_of(compressed_phase: &str) -> Result<RootSide> {
    match compressed_phase.trim().to_lowercase().as_str() {
        "liquid" => Ok(RootSide::Liquid),
        "vapour" | "vapor" => Ok(RootSide::Vapour),
        other => Err(AzothError::invalid_input(
            "compressed_phase",
            format!("{other:?} is not a side. The spec's values are \"liquid\" and \"vapour\""),
        )),
    }
}

/// The same state, for a caller that resolved the fluid itself.
///
/// # Errors
/// * [`AzothError::InvalidInput`] if `z` is not one entry per component or does not sum to
///   one.
/// * [`AzothError::OutOfRange`] as [`molar_volume`] and
///   [`crate::pcsaft::ln_fugacity_coefficients`].
pub fn pcsaft_phase_of(
    components: &[PcsaftComponent],
    kij: &[f64],
    t: ThermodynamicTemperature,
    p: Pressure,
    z: &[f64],
    side: RootSide,
) -> Result<PcsaftPhaseResult> {
    let n = components.len();
    if z.len() != n {
        return Err(AzothError::invalid_input(
            "z",
            format!(
                "a mixture of {n} components needs {n} mole fractions, but z has {}",
                z.len()
            ),
        ));
    }
    let sum: f64 = z.iter().sum();
    if (sum - 1.0).abs() > 1.0e-9 {
        return Err(AzothError::invalid_input(
            "z",
            format!(
                "the mole fractions sum to {sum}, not to one. Renormalising them here would \
                 make a composition error invisible in every number downstream, so it is \
                 refused instead"
            ),
        ));
    }

    let solved = molar_volume(components, kij, z, t.value, p.value, side)?;
    Ok(PcsaftPhaseResult {
        z_factor: solved.z,
        ln_phi: crate::pcsaft::ln_fugacity_coefficients(components, kij, z, t.value, solved.v)?,
        v: cubic_meters_per_mole(solved.v),
        warnings: Vec::new(),
    })
}
