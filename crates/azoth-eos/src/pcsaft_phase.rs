//! PC-SAFT's volume at a temperature, a pressure and a composition.
//!
//! ```text
//! solve P_calc(v) = P   for   v,   Z = P v/(R T)
//! ```
//!
//! The layers the pressure is built from are [`crate::pcsaft`]'s; this module is the
//! solve that turns them into a volume, so the layers stay checkable one at a time.
//!
//! The *models* are thin modules over it - `pcsaft_rahmat_phase` is the one NeqSim's
//! `SystemPCSAFT` runs - because the solve is the same arithmetic and two copies would
//! invite them to disagree.
//!
//! # The two roots
//!
//! A PC-SAFT isotherm is not monotonic: below the critical temperature the pressure dips
//! and rises again, so a temperature and a pressure can admit three volumes. Which one
//! describes the phase wanted is the caller's statement, not the fluid's, and it is the
//! same `RootSide` the cubic roots are asked for by.

use azoth_core::{AzothError, Result};

use crate::association::R;
use crate::mixture::RootSide;
use crate::pcsaft::{self, PcsaftComponent};

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

    let (v, iterations) = crate::volume_solve::molar_volume(
        |v| Ok((h(v)?, slope(v)?)),
        ideal,
        floor,
        side,
        |branch| no_root(branch, t, p),
    )?;
    let eta = std::f64::consts::PI / 6.0 * pcsaft::AVOGADRO * md3 / v;
    Ok(PcsaftVolume {
        v,
        z: p * v / rt,
        eta,
        iterations,
    })
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

/// A PC-SAFT phase's residual enthalpy and entropy, per mole.
///
/// **`Z - 1` is the volume's own contribution and `-T dF/dT` the temperature's.** NeqSim
/// assembles the same two from `AresTV + T SresTV + P V - n R T` with
/// `SresTV = (-T dFdT - F) R`, and the two `F`-terms cancel: `Hres/(R T) = Z - 1 - T dFdT`.
/// The entropy then carries the `P`-to-`V` conversion, `SresTP = SresTV + n R ln Z`.
///
/// **This is a finite difference of NeqSim's own `F` at fixed volume before it is
/// NeqSim's `dFdT`**, because for PC-SAFT the two disagree:
/// [`crate::pcsaft::t_d_helmholtz_rt_dt`] carries the note and the write-up is at
/// `~/Desktop/neqsim-pcsaft-hard-chain-temperature-derivative.md`.
///
/// # Errors
/// As [`crate::pcsaft::state`].
pub fn departure(
    components: &[PcsaftComponent],
    kij: &[f64],
    z: &[f64],
    t: f64,
    v: f64,
    compressibility: f64,
) -> Result<(f64, f64)> {
    let state = crate::pcsaft::state(components, kij, z, t, v)?;
    let t_d_f = crate::pcsaft::t_d_helmholtz_rt_dt(components, z, t, v, &state);
    let h_over_rt = compressibility - 1.0 - t_d_f;
    let s_over_r = compressibility.ln() - t_d_f - state.f();
    Ok((h_over_rt, s_over_r))
}
