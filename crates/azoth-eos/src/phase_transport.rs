//! `eos.phase_transport` - a phase's transport properties, by the dispatch NeqSim's
//! `PhysicalProperties` makes.
//!
//! Spec: `specs/models/eos/phase_transport.toml`, which records the three phase kinds' model
//! selections and the two ladders the diffusivity chain reads.

use azoth_core::units::{
    Pressure, ThermodynamicTemperature, pascal_seconds, watts_per_meter_kelvin,
};
use azoth_core::{AzothError, Result, apply_checks};

use crate::IdealGasModel;
use crate::chapman_enskog_diffusivity::chapman_enskog_diffusivity;
use crate::databank::{lennard_jones_pair, normal_boiling_molar_volume};
use crate::effective_diffusion::effective_diffusion;
use crate::liquid_conductivity_polynom::liquid_conductivity_polynom;
use crate::liquid_viscosity_pure::{LiquidViscosityLadder, liquid_viscosity_pure};
use crate::mixture::Mixture;
use crate::model_gen;
use crate::results::PhaseTransportResult;
use crate::siddiqi_lucas_diffusivity::{SiddiqiLucasForm, siddiqi_lucas_diffusivity};
use crate::thermal_conductivity::thermal_conductivity;
use crate::viscosity::viscosity;

/// The three phase kinds `PhysicalProperties` dispatches on.
///
/// NeqSim's `PhaseType`, and the dispatch is by it: an aqueous phase takes `WaterPhysicalProperties`,
/// whose viscosity and conductivity are the *polynom* correlations, while a gas and a hydrocarbon
/// liquid take PFCT for both.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PhaseKind {
    /// `PhaseType.GAS`: PFCT viscosity and conductivity, Chapman-Enskog diffusivity.
    Gas,
    /// `PhaseType.OIL`: PFCT for both, Siddiqi-Lucas diffusivity on the common-phase ladder.
    Oil,
    /// `PhaseType.AQUEOUS`: the polynom viscosity and conductivity, Siddiqi-Lucas diffusivity on
    /// the liquid ladder.
    Aqueous,
}

impl PhaseKind {
    /// The spec's spelling.
    #[must_use]
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::Gas => "gas",
            Self::Oil => "oil",
            Self::Aqueous => "aqueous",
        }
    }
}

impl std::str::FromStr for PhaseKind {
    type Err = String;

    fn from_str(s: &str) -> std::result::Result<Self, Self::Err> {
        match s {
            "gas" => Ok(Self::Gas),
            "oil" => Ok(Self::Oil),
            "aqueous" => Ok(Self::Aqueous),
            other => Err(format!(
                "`{other}` is not one of NeqSim's three phase kinds: gas, oil or aqueous"
            )),
        }
    }
}

/// The lowest pure-component viscosity the liquid diffusivity correlations take, in cP.
///
/// `SiddiqiLucasMethod.calcBinaryDiffusionCoefficient` clamps `eta` to this, and the clamp is
/// reachable: the common-phase ladder answers **zero** for a component whose LIQVISC model is 2.
const MIN_ETA_CP: f64 = 0.01;

/// A phase's viscosity, thermal conductivity, binary diffusivity matrix and effective
/// diffusivities - the four the segment model's transport snapshot reads.
///
/// **The dispatch is the answer, not a detail.** NeqSim's `PhysicalProperties` is created by
/// phase type and each subclass assigns its own models, so the same fluid's gas and liquid answer
/// their properties from *different correlations*. The capture holds both the values and the model
/// classes per phase, which is what the cases are held to.
///
/// # Errors
/// * [`azoth_core::AzothError::InvalidInput`] if `z` is not one entry per component, or if the
///   phase has fewer than two components - the effective diffusivity divides by an empty sum
///   there and NeqSim answers a `NaN`.
/// * Whatever the composed correlations refuse.
pub fn phase_transport(
    mixture: &Mixture,
    ideal_gas: &IdealGasModel,
    phase: PhaseKind,
    t: ThermodynamicTemperature,
    p: Pressure,
    z: &[f64],
) -> Result<PhaseTransportResult> {
    let spec = &model_gen::PHASE_TRANSPORT_SPEC;
    let mut warnings = Vec::new();

    let components = mixture.components();
    let count = components.len();
    if z.len() != count {
        return Err(AzothError::invalid_input(
            "z",
            format!(
                "a mixture of {count} components needs {count} mole fractions, got {}",
                z.len()
            ),
        ));
    }
    if count < 2 {
        return Err(AzothError::invalid_input(
            "components",
            "one component: the effective diffusivity divides by the sum over the *other* \
             components, so NeqSim answers a NaN here and this refuses instead of reporting \
             one"
            .to_string(),
        ));
    }

    // Every correlation here is a function of the components' molar masses, and NeqSim
    // refuses a component that has none rather than inventing one.
    let molar_mass: Vec<f64> = components
        .iter()
        .map(|component| {
            component.molar_mass.ok_or_else(|| {
                AzothError::invalid_input(
                    "components",
                    "a component carries no molar mass, and every correlation here is a                      function of one",
                )
            })
        })
        .collect::<Result<Vec<f64>>>()?;

    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "T" => Some(t.value),
            "P" => Some(p.value),
            "z" => z.iter().copied().reduce(f64::min),
            _ => None,
        },
        &mut warnings,
    )?;

    // ---- The viscosity and the conductivity, by phase kind.
    let (mu, k) = match phase {
        PhaseKind::Aqueous => (
            crate::aqueous_viscosity::aqueous_viscosity(mixture, t, p, z)?
                .viscosity
                .value,
            liquid_conductivity_polynom(
                &components
                    .iter()
                    .map(|component| component.liquid_conductivity)
                    .collect::<Vec<[f64; 3]>>(),
                &molar_mass
                    .iter()
                    .map(|mass| azoth_core::units::kilograms_per_mole(*mass))
                    .collect::<Vec<_>>(),
                z,
                t,
            )?
            .k
            .value,
        ),
        PhaseKind::Gas | PhaseKind::Oil => (
            viscosity(mixture, t, p, z)?.mu.value,
            thermal_conductivity(mixture, ideal_gas, t, p, z)?.k.value,
        ),
    };

    // ---- The binary diffusivity matrix, by phase kind.
    let d_binary = match phase {
        PhaseKind::Gas => {
            let mut matrix = vec![vec![0.0; count]; count];
            for i in 0..count {
                for j in 0..count {
                    if i == j {
                        continue;
                    }
                    let (sigma, eps, _) = lennard_jones_pair(
                        components[i].lennard_jones_diameter,
                        components[i].lennard_jones_energy,
                        components[j].lennard_jones_diameter,
                        components[j].lennard_jones_energy,
                        molar_mass[i],
                        molar_mass[j],
                    );
                    matrix[i][j] = chapman_enskog_diffusivity(
                        azoth_core::units::kilograms_per_mole(molar_mass[i]),
                        azoth_core::units::kilograms_per_mole(molar_mass[j]),
                        sigma * 1.0e-10,
                        azoth_core::units::kelvins(eps),
                        t,
                        p,
                    )?
                    .d
                    .value;
                }
            }
            matrix
        }
        PhaseKind::Oil | PhaseKind::Aqueous => {
            // **The default is the aqueous form for every liquid.** `SiddiqiLucasMethod`'s
            // `autoSelectCorrelation` is false, and `setDiffusionCoefficientModel` is the only
            // thing that turns it on - so an oil takes the aqueous correlation here.
            //
            // **And the two phase kinds take different ladders** for the `eta` it divides by:
            // the oil phase's PFCT viscosity class inherits the common-phase one.
            let ladder = match phase {
                PhaseKind::Aqueous => LiquidViscosityLadder::Liquid,
                _ => LiquidViscosityLadder::CommonPhase,
            };
            let mut matrix = vec![vec![0.0; count]; count];
            for i in 0..count {
                for j in 0..count {
                    if i == j {
                        continue;
                    }
                    let va = normal_boiling_molar_volume(
                        components[i].normal_liquid_density,
                        molar_mass[i],
                        components[i].critical_volume,
                    );
                    let vb = normal_boiling_molar_volume(
                        components[j].normal_liquid_density,
                        molar_mass[j],
                        components[j].critical_volume,
                    );
                    let pure_cp = liquid_viscosity_pure(
                        ladder,
                        components[j].liqvisc_model,
                        components[j].liqvisc[0],
                        components[j].liqvisc[1],
                        components[j].liqvisc[2],
                        components[j].liqvisc[3],
                        components[j].tc,
                        components[j].pc,
                        components[j].omega,
                        t,
                        p,
                    )?
                    .mu
                    .value
                        * 1000.0;
                    // **The class's own fallback, and on this path it is the normal one**: a
                    // pure-component viscosity that is not positive - exactly what the
                    // common-phase ladder's empty model 2 branch answers - takes the *phase's*
                    // viscosity instead. Measured on the captured oil phase, where the hole
                    // fires and `eta` is the phase's own `0.1347` cP rather than the `0.01`
                    // the clamp would otherwise floor it to.
                    let eta_cp = if pure_cp.is_finite() && pure_cp > 0.0 {
                        pure_cp
                    } else {
                        mu * 1000.0
                    };
                    matrix[i][j] = siddiqi_lucas_diffusivity(
                        SiddiqiLucasForm::Aqueous,
                        azoth_core::units::cubic_meters_per_mole(va * 1.0e-6),
                        azoth_core::units::cubic_meters_per_mole(vb * 1.0e-6),
                        t,
                        pascal_seconds(eta_cp.max(MIN_ETA_CP) * 1.0e-3),
                    )?
                    .d
                    .value;
                }
            }
            matrix
        }
    };

    // **A phase that is one substance has no effective diffusivity, and that is not this
    // id's refusal.** The assembly divides by the *other* components' fractions, so on a
    // pure phase there is nothing for the one component to diffuse into and
    // `eos.effective_diffusion` refuses - rightly, for the question *it* asks.
    // `RateBasedPackedColumnTest`'s own lean solvent is pure water, and NeqSim answers a
    // **zero** vector there rather than refusing: measured, the capture's absorber row prints
    // `0.0` for every entry of both phases. So the vector is zeros, and the three properties
    // that *are* defined - the viscosity, the conductivity and the pair matrix - are reported.
    let d_effective = match effective_diffusion(&d_binary, z) {
        Ok(assembled) => {
            warnings.extend(assembled.warnings.iter().cloned());
            assembled.effective_diffusion
        }
        Err(AzothError::OutOfRange { ref field, .. }) if field == "x" => vec![0.0; count],
        Err(other) => return Err(other),
    };

    apply_checks(
        spec.derived_checks(),
        |name| match name {
            "mu" => Some(mu),
            "k" => Some(k),
            _ => None,
        },
        &mut warnings,
    )?;

    Ok(PhaseTransportResult {
        mu: pascal_seconds(mu),
        k: watts_per_meter_kelvin(k),
        d_binary,
        d_effective,
        warnings,
    })
}
