//! `eos.liquid_viscosity_pure` - one component's pure-liquid viscosity, from the LIQVISC
//! correlation.
//!
//! Spec: `specs/calcs/eos/liquid_viscosity_pure.toml`, which records the four polynomial forms,
//! the pressure correction and **the two NeqSim classes whose ladders disagree**.

use azoth_core::units::{Pressure, ThermodynamicTemperature, pascal_seconds};
use azoth_core::{Result, apply_checks};

use crate::results::LiquidViscosityPureResult;
use crate::spec_gen;

/// The value NeqSim answers above a component's critical temperature, in cP.
const ABOVE_CRITICAL_CP: f64 = 0.5;
/// The value it answers for a component with no LIQVISC model, in cP.
const NO_MODEL_CP: f64 = 0.7;
/// A very small viscosity would divide; NeqSim clamps the *diffusivity's* use of this to 0.01 cP,
/// and the value itself is left as it comes out.
const CP_TO_PA_S: f64 = 1.0e-3;

/// **Which of NeqSim's two ladders to use.** They differ in exactly one branch.
///
/// `commonphasephysicalproperties.viscosity.Viscosity` leaves LIQVISC model 2's branch **empty**,
/// so a component with that model comes out at `0.0` cP and only the pressure correction scales
/// it; `liquidphysicalproperties.viscosity.Viscosity` implements it as `exp(l1 + l2/T)`. A phase
/// gets whichever class its own viscosity model extends, and both are on the rate-based column's
/// path: an oil phase's `PFCTViscosityMethodHeavyOil` inherits the first, an aqueous phase's
/// "polynom" model *is* the second.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LiquidViscosityLadder {
    /// `Viscosity` in `commonphasephysicalproperties` - model 2 answers nothing.
    CommonPhase,
    /// `Viscosity` in `liquidphysicalproperties` - every model answers.
    Liquid,
}

impl LiquidViscosityLadder {
    /// The spec's spelling.
    #[must_use]
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::CommonPhase => "common_phase",
            Self::Liquid => "liquid",
        }
    }
}

impl std::str::FromStr for LiquidViscosityLadder {
    type Err = String;

    fn from_str(s: &str) -> std::result::Result<Self, Self::Err> {
        match s {
            "common_phase" => Ok(Self::CommonPhase),
            "liquid" => Ok(Self::Liquid),
            other => Err(format!(
                "`{other}` is not one of the two ladders: common_phase or liquid"
            )),
        }
    }
}

/// One component's pure-liquid viscosity, from its LIQVISC correlation.
///
/// The four models are NeqSim's own: `1` a power law, `2` an Arrhenius pair, `3` a four-term
/// exponential and `4` the Yaws form. `l1..l4` are the component's `LIQVISC1..4` columns, `tc` and
/// `pc` its critical point and `omega` its acentric factor, which the pressure correction reads.
///
/// The result is in Pa*s, where NeqSim's own value is in cP.
///
/// # Errors
/// * [`azoth_core::AzothError::OutOfRange`] if `T`, `Tc` or `Pc` is not positive.
///
/// # Example
/// ```
/// use azoth_core::units::{kelvins, pascals};
/// use azoth_eos::liquid_viscosity_pure::{LiquidViscosityLadder, liquid_viscosity_pure};
///
/// // Water at 313.15 K and 50 bar, the CO2 absorber's aqueous phase.
/// let r = liquid_viscosity_pure(
///     LiquidViscosityLadder::Liquid, 3,
///     -27.952757828, 4665.22592993, 0.052323342, -3.8356e-5,
///     kelvins(647.3), pascals(22089000.0), 0.344,
///     kelvins(313.15), pascals(50.0e5),
/// )?;
/// assert!((r.mu.value - 6.528922494381046e-4).abs() < 1e-18);
/// # Ok::<(), azoth_core::AzothError>(())
/// ```
#[allow(non_snake_case)] // `T`, `P`, `Tc`, `Pc` are the symbols in the correlation
#[allow(clippy::too_many_arguments)] // The signature is the spec's declared inputs.
pub fn liquid_viscosity_pure(
    form: LiquidViscosityLadder,
    model: u32,
    l1: f64,
    l2: f64,
    l3: f64,
    l4: f64,
    Tc: ThermodynamicTemperature,
    Pc: Pressure,
    omega: f64,
    T: ThermodynamicTemperature,
    P: Pressure,
) -> Result<LiquidViscosityPureResult> {
    let spec = &spec_gen::LIQUID_VISCOSITY_PURE_SPEC;
    let mut warnings = Vec::new();

    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "Tc" => Some(Tc.value),
            "Pc" => Some(Pc.value),
            "T" => Some(T.value),
            "P" => Some(P.value),
            _ => None,
        },
        &mut warnings,
    )?;

    let uncorrected_cp = if T.value > Tc.value {
        ABOVE_CRITICAL_CP
    } else {
        match (form, model) {
            (_, 1) => l1 * T.value.powf(l2),
            // **Model 2 is the one branch the two ladders disagree on.**
            (LiquidViscosityLadder::Liquid, 2) => (l1 + l2 / T.value).exp(),
            (LiquidViscosityLadder::CommonPhase, 2) => 0.0,
            (_, 3) => (l1 + l2 / T.value + l3 * T.value + l4 * T.value.powi(2)).exp(),
            (_, 4) => 10f64.powf(l1 * (1.0 / T.value - 1.0 / l2)),
            _ => NO_MODEL_CP,
        }
    };

    let mu_cp = uncorrected_cp
        * (pressure_correction(T.value, P.value, Tc.value, Pc.value, omega) + 1.0)
        / 2.0;
    let mu = mu_cp * CP_TO_PA_S;

    apply_checks(
        spec.derived_checks(),
        |name| (name == "mu").then_some(mu),
        &mut warnings,
    )?;

    Ok(LiquidViscosityPureResult {
        mu: pascal_seconds(mu),
        warnings,
    })
}

/// `getViscosityPressureCorrection`: NeqSim's own four-coefficient form.
///
/// The Lucas pressure term, with the reduced temperature in the numerator of each exponential's
/// coefficient and the acentric factor in the denominator.
#[must_use]
pub fn pressure_correction(temperature: f64, pressure: f64, tc: f64, pc: f64, omega: f64) -> f64 {
    let reduced_t = temperature / tc;
    if reduced_t > 1.0 {
        return 1.0;
    }
    let delta_pr = pressure / pc;
    let a = 0.9991 - (4.674e-4 / (1.0523 * reduced_t.powf(-0.03877) - 1.0513));
    let d = (0.3257 / (1.0039 - reduced_t.powf(2.573)).powf(0.2906)) - 0.2086;
    let c = -0.07921 + 2.1616 * reduced_t - 13.4040 * reduced_t.powi(2)
        + 44.1706 * reduced_t.powi(3)
        - 84.8291 * reduced_t.powi(4)
        + 96.1209 * reduced_t.powi(5)
        - 59.8127 * reduced_t.powi(6)
        + 15.6719 * reduced_t.powi(7);
    (1.0 + d * (delta_pr / 2.118).powf(a)) / (1.0 + c * omega * delta_pr)
}
