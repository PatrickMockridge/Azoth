//! `eos.vu_flash_single_comp` - the volume/internal-energy state of a pure component.
//!
//! Spec: `specs/models/eos/vu_flash_single_comp.toml`. NeqSim's `VUflashSingleComp`.
//!
//! A pure component has no two-phase split for a flash to find. One K-value is either
//! above one or below it, so [`crate::pt_flash`] reports a single phase at every
//! temperature and pressure and [`crate::vu_flash`]'s iteration - which is a Newton on
//! `(P, T)` around a flash - has nothing to converge on. It refuses, and it refuses for
//! every pure feed in the two-phase region, which is a large part of what a
//! depressurisation calculation asks for.
//!
//! The state is not unknown there, it is *constrained*: a pure component at a pressure
//! is on its saturation line, so the temperature is the saturation temperature and the
//! split is a lever rule on the two saturated internal energies. That is this model,
//! and it is exact for a pure component where a flash would be an approximation of it.

use azoth_core::units::{
    MolarEnergy, MolarVolume, Pressure, ThermodynamicTemperature, cubic_meters_per_mole, kelvins,
};
use azoth_core::{AzothError, Result, Warning, apply_checks};

use crate::algorithm_of;
use crate::mixture::{Mixture, RootSide};
use crate::model_gen;
use crate::molar_enthalpy_entropy::{IdealGasModel, molar_enthalpy_entropy};
use crate::pr_molar_volume::MOLAR_GAS_CONSTANT;
use crate::pure_saturation::pure_saturation;
use crate::results::VuFlashSingleCompResult;

/// One saturated phase's internal energy and molar volume at a state.
struct Saturated {
    u: f64,
    v: f64,
}

/// The saturated liquid and vapour of a pure component at a temperature and pressure.
///
/// Both are the *feed's* own composition - there is only one - so the two differ in the
/// cubic root alone, which is what "saturated" means for a pure component.
fn saturated(
    mixture: &Mixture,
    ideal_gas: &IdealGasModel,
    t: ThermodynamicTemperature,
    p: Pressure,
) -> Result<(Saturated, Saturated)> {
    let reduced = mixture.reduced_parameters(t, p)?;
    let x = [1.0];
    let liquid = mixture.phase_state(&reduced, &x, RootSide::Liquid)?;
    let vapour = mixture.phase_state(&reduced, &x, RootSide::Vapour)?;
    let mut out = Vec::with_capacity(2);
    for z in [liquid.z, vapour.z] {
        let state = molar_enthalpy_entropy(mixture, ideal_gas, t, p, &x, z)?;
        let v = z * MOLAR_GAS_CONSTANT * t.value / p.value;
        out.push(Saturated {
            u: state.h.value - p.value * v,
            v,
        });
    }
    let vapour_phase = out.pop().expect("two phases");
    let liquid_phase = out.pop().expect("two phases");
    Ok((liquid_phase, vapour_phase))
}

/// The temperature at which a pure component's saturation pressure is `p`.
///
/// Bisection on the registered `eos.pure_saturation`, which solves the forward problem -
/// the pressure at a temperature - and is monotone in it. The bracket is the subcritical
/// span the spec declares, and a pressure outside what that span reaches is a state this
/// model has no answer for rather than one to iterate towards.
fn saturation_temperature(
    mixture: &Mixture,
    p: Pressure,
    algorithm: &azoth_core::ModelAlgorithm,
) -> Result<f64> {
    let component = &mixture.components()[0];
    let bracket = algorithm.bracket.ok_or_else(|| AzothError::InvalidInput {
        field: "algorithm.bracket".to_string(),
        reason: format!(
            "scheme `{}` bisects a bracket but the spec declares none",
            algorithm.scheme
        ),
    })?;
    let tc = component.tc.value;
    if p.value >= component.pc.value {
        return Err(AzothError::OutOfRange {
            field: "P".to_string(),
            value: p.value,
            detail: format!(
                "the critical pressure is {:.6e} Pa, so {:.6e} Pa is above the saturation \
                 line at every temperature and a pure component there has no two-phase \
                 state to find - NeqSim's own guard is the same comparison",
                component.pc.value, p.value
            ),
        });
    }
    let (mut lo, mut hi) = (bracket.lower * tc, bracket.upper * tc);

    // Whether the saturation pressure at `t` is at most `p`, which is the question the
    // bisection asks. `None` is a temperature where `eos.pure_saturation` cannot be
    // evaluated - only reachable at the top of the bracket, where the roots have
    // coalesced - and there the answer is *below*: the saturation pressure at such a
    // temperature is at or above the critical pressure, and the critical pressure has
    // already been excluded above. So `None` is read as "above", which is one-sided and
    // needs no temperature at which the calculation is assumed to work.
    let at_most = |t: f64, p: f64| -> Result<Option<bool>> {
        match pure_saturation(component.tc, component.pc, component.omega, kelvins(t)) {
            Ok(result) => Ok(Some(result.p_sat.value <= p)),
            Err(AzothError::SolverNotConverged { .. }) => Ok(None),
            Err(other) => Err(other),
        }
    };

    for _ in 0..algorithm.max_iterations {
        let mid = 0.5 * (lo + hi);
        if (hi - lo) / mid <= algorithm.tolerance {
            break;
        }
        if at_most(mid, p.value)? == Some(true) {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    Ok(0.5 * (lo + hi))
}

/// The temperature and split of a pure component at a pressure and internal energy.
///
/// # Errors
/// * [`AzothError::InvalidInput`] if the mixture is not a single component, or the feed
///   is not one mole fraction summing to one.
/// * [`AzothError::OutOfRange`] if `P` is at or above the critical pressure, or `U` lies
///   outside the two saturated internal energies at it - both are states this model has
///   no split for rather than ones it failed to find.
///
/// # Example
/// ```
/// use azoth_core::units::{cubic_meters_per_mole, joules_per_mole, pascals};
/// use azoth_eos::databank;
///
/// let (mixture, ideal_gas) = databank::mixture_of(&["propane"], None)?;
/// use azoth_eos::vu_flash_single_comp::vu_flash_single_comp;
///
/// let r = vu_flash_single_comp(
///     &mixture,
///     &ideal_gas,
///     pascals(1.0e6),
///     cubic_meters_per_mole(0.001_059_848_052_516),
///     joules_per_mole(-7797.318_485_008),
/// )?;
/// assert!((r.beta - 0.5).abs() < 1e-6);
/// # Ok::<(), azoth_core::AzothError>(())
/// ```
pub fn vu_flash_single_comp(
    mixture: &Mixture,
    ideal_gas: &IdealGasModel,
    p: Pressure,
    v: MolarVolume,
    u: MolarEnergy,
) -> Result<VuFlashSingleCompResult> {
    let spec = &model_gen::VU_FLASH_SINGLE_COMP_SPEC;
    let mut warnings = Vec::new();

    if mixture.len() != 1 {
        return Err(AzothError::invalid_input(
            "components",
            format!(
                "a pure component's saturation state needs one component and this \
                 mixture has {}. A mixture's two-phase split is `eos.vu_flash`",
                mixture.len()
            ),
        ));
    }
    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "P" => Some(p.value),
            _ => None,
        },
        &mut warnings,
    )?;

    let algorithm = algorithm_of(spec)?;
    let temperature = saturation_temperature(mixture, p, algorithm)?;
    let (liquid, vapour) = saturated(mixture, ideal_gas, kelvins(temperature), p)?;
    let span = vapour.u - liquid.u;
    if span <= 0.0 || u.value < liquid.u || u.value > vapour.u {
        return Err(AzothError::OutOfRange {
            field: "U".to_string(),
            value: u.value,
            detail: format!(
                "the saturated liquid's internal energy at {:.4} K and {:.6e} Pa is \
                 {:.6e} J/mol and the saturated vapour's is {:.6e}, so a feed of \
                 {:.6e} J/mol is outside the two-phase region and has no split. A \
                 subcooled or superheated pure component has no vapour fraction",
                temperature, p.value, liquid.u, vapour.u, u.value
            ),
        });
    }

    let beta = (u.value - liquid.u) / span;
    let implied = (1.0 - beta) * liquid.v + beta * vapour.v;
    let given = v.value;
    if (implied - given).abs() > 1.0e-6 * implied.abs() {
        warnings.push(Warning::new(
            azoth_core::WarningCode::OutOfValidRange,
            format!(
                "the volume asked for, {given:.9e} m**3/mol, is not the volume this split \
                 implies, {implied:.9e}. The pressure and the internal energy fix the \
                 volume for a pure component - the two phases are saturated - so the \
                 volume is reported from them rather than solved for, and the caller's \
                 is inconsistent with the state they asked about"
            ),
        ));
    }

    Ok(VuFlashSingleCompResult {
        t: kelvins(temperature),
        beta,
        v: cubic_meters_per_mole(implied),
        phase: crate::results::Phase::TwoPhase,
        warnings,
    })
}
