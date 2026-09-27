//! `AbsorptionColumn.applyMurphreeCorrection`, which **replaces** the base column's.
//!
//! **The base blends one phase and this blends two.** `DistillationColumn.applyMurphreeCorrection`
//! writes the vapour leaving a stage towards the vapour entering it at the flash's own moles; an
//! absorber's override instead corrects **both** the vapour and the liquid towards their inlets'
//! counterparts, renormalises the corrected vapour fraction, and then **re-allocates the flash's
//! vapour moles across the components** through a limiting-component loop - so what a stage hands
//! up carries a different *composition* at the same total moles, and what it hands down carries
//! the remainder. The two are different arithmetic rather than two spellings of one.
//!
//! **Four levels of resolution, and this is where the extra two live.**
//! `getComponentMurphreeEfficiency(tray, component)` reads
//! `trayComponentMurphreeEfficiencies[tray][component]`, then `componentMurphreeEfficiencies
//! [component]`, then falls through to the base's two-step `getMurphreeEfficiency(tray)`. Only
//! the last three are declared - `AbsorptionColumn.setComponentMurphreeEfficiency(int, String,
//! double)` writes the first, which has no palette spelling and is not carried.
//!
//! **`finalizeTrayProperties` invalidates both caches at the end of a solve**, so a capture of a
//! converged absorber shows the equilibrium rebuild on both phases and not this correction. That
//! is why the oracle pins the *profile* rather than the corrected streams - reconstructing them
//! in a test would mean reproducing [`allocate_vapour_moles`] here, which is testing the
//! implementation against itself.

use azoth_core::{AzothError, Result};

use crate::column::murphree::Murphree;
use crate::stream::Stream;

/// The class's own tolerance, `AbsorptionColumn.MOLE_TOLERANCE`.
const MOLE_TOLERANCE: f64 = 1.0e-15;

/// The class's own ideal-stage tolerance, `applyMurphreeCorrection`'s.
const IDEAL_TOLERANCE: f64 = 1.0e-10;

/// `AbsorptionColumn`'s two efficiency maps over the base's two levels.
#[derive(Debug, Clone, PartialEq)]
pub struct AbsorberMurphree {
    /// The base's `murphreeEfficiency` and `perStageMurphreeEfficiency`, which
    /// `getComponentMurphreeEfficiency` falls through to.
    pub base: Murphree,
    /// `componentMurphreeEfficiencies`, one entry per component **in the column's own order** -
    /// the class keys it by normalized name and this is that map resolved against one order.
    pub per_component: Option<Vec<f64>>,
}

impl AbsorberMurphree {
    /// `getComponentMurphreeEfficiency`: the component's own value, else the base's two steps.
    #[must_use]
    pub fn resolve(&self, tray: usize, component: usize) -> f64 {
        if let Some(per_component) = self.per_component.as_ref() {
            if let Some(value) = per_component.get(component) {
                if !value.is_nan() {
                    return *value;
                }
            }
        }
        self.base.resolve(tray)
    }

    /// `clampMurphreeEfficiency` on every level, and the component vector's length checked
    /// against the components the column carries.
    ///
    /// # Errors
    /// [`AzothError::InvalidInput`] for a component vector that is not one entry per component.
    pub fn checked(&self, components: usize) -> Result<&Self> {
        if let Some(per_component) = self.per_component.as_ref() {
            if per_component.len() != components {
                return Err(AzothError::invalid_input(
                    "component_murphree_efficiency",
                    format!(
                        "{} component efficiency(ies) against {components} component(s): the \
                         vector is one entry per component in the column's own order, and a \
                         component that states none is written `NaN` rather than left out",
                        per_component.len()
                    ),
                ));
            }
        }
        Ok(self)
    }
}

/// `AbsorptionColumn.applyMurphreeCorrection`: the two phases leaving one stage, corrected.
///
/// `here` is the stage's own flash: its vapour is the equilibrium phase and its liquid the other.
/// `below` is the stage beneath, read the same way. **`below` is the tray's *gas* inlet** - at
/// stage 0 the class reads `gasInStream`, and elsewhere `getTray(trayIndex - 1).getThermoSystem()`,
/// which is the stage below's flash and not its corrected outlet.
///
/// `None` where the class would return without touching anything: a stage whose flash found one
/// phase, one whose inlet holds none, or one where no component's efficiency bites.
///
/// # Errors
/// [`AzothError::InvalidInput`] when the stage and its inlet do not carry the same substances in
/// the same order, and whatever building the corrected streams refuses.
pub fn correct(
    here: (Option<&Stream>, Option<&Stream>),
    below: Option<&Stream>,
    efficiency: &AbsorberMurphree,
    tray: usize,
) -> Result<Option<(Stream, Stream)>> {
    let (Some(equilibrium_gas), Some(equilibrium_liquid)) = here else {
        return Ok(None);
    };
    let Some(inlet) = below else {
        return Ok(None);
    };
    if equilibrium_gas.components != inlet.components {
        return Err(AzothError::invalid_input(
            "component_murphree_efficiency",
            "the stage and the stage below it do not carry the same substances in the same \
             order, so the correction has no basis for pairing one's components with the \
             other's: `AbsorptionColumn.getComponentMurphreeEfficiency` normalizes a component \
             *name* against the stage's own system",
        ));
    }

    let count = equilibrium_gas.components.len();
    let mut corrected = vec![0.0; count];
    let mut total = vec![0.0; count];
    let mut required = false;
    let mut sum = 0.0;
    for component in 0..count {
        let density = efficiency.resolve(tray, component);
        required |= density < 1.0 - IDEAL_TOLERANCE;
        let equilibrium = equilibrium_gas.z[component];
        let inlet_fraction = inlet.z[component];
        let blended = (inlet_fraction + density * (equilibrium - inlet_fraction)).max(0.0);
        corrected[component] = blended;
        sum += blended;
        // **Both phases' moles of this component**, which is what the allocator draws on.
        total[component] = (equilibrium_gas.z[component] * equilibrium_gas.n
            + equilibrium_liquid.z[component] * equilibrium_liquid.n)
            .max(0.0);
    }
    if !required || sum <= MOLE_TOLERANCE {
        return Ok(None);
    }
    for fraction in &mut corrected {
        *fraction /= sum;
    }

    let vapour_moles = equilibrium_gas.n;
    let liquid_moles = equilibrium_liquid.n;
    let allocated = allocate_vapour_moles(&corrected, &total, vapour_moles);
    let gas_z: Vec<f64> = allocated
        .iter()
        .map(|moles| {
            if vapour_moles > MOLE_TOLERANCE {
                moles / vapour_moles
            } else {
                0.0
            }
        })
        .collect();
    let liquid_z: Vec<f64> = allocated
        .iter()
        .zip(&total)
        .map(|(gas, available)| {
            let moles = (available - gas).max(0.0);
            if liquid_moles > MOLE_TOLERANCE {
                moles / liquid_moles
            } else {
                0.0
            }
        })
        .collect();

    Ok(Some((
        Stream::from_pt(
            equilibrium_gas.components.clone(),
            gas_z,
            vapour_moles,
            equilibrium_gas.p,
            equilibrium_gas.t,
        )?,
        Stream::from_pt(
            equilibrium_liquid.components.clone(),
            liquid_z,
            liquid_moles,
            equilibrium_liquid.p,
            equilibrium_liquid.t,
        )?,
    )))
}

/// `AbsorptionColumn.allocateVaporMoles`: the flash's vapour moles split across the components by
/// the corrected fraction, **with a component whose share exceeds what it has pinned at what it
/// has**.
///
/// The loop is the class's, including its two fallbacks: while a component's trial share is more
/// than the moles available to it, that component is fixed at its availability and the remainder
/// is redistributed. Where no component is limited, the whole remainder is handed out at once.
/// The `remainingWeight <= MOLE_TOLERANCE` branch divides by the *available* moles instead, which
/// is reached only when every unfixed component's corrected fraction has collapsed to zero.
#[must_use]
pub fn allocate_vapour_moles(
    mole_fraction: &[f64],
    available_moles: &[f64],
    vapour_moles: f64,
) -> Vec<f64> {
    let count = mole_fraction.len();
    let mut allocated = vec![0.0; count];
    let mut fixed = vec![false; count];
    let mut remaining = vapour_moles;

    for _ in 0..count {
        if remaining <= MOLE_TOLERANCE {
            break;
        }
        let weight: f64 = mole_fraction
            .iter()
            .zip(&fixed)
            .filter(|(_, fixed)| !**fixed)
            .map(|(fraction, _)| *fraction)
            .sum();

        let mut limited = false;
        for component in 0..count {
            if fixed[component] {
                continue;
            }
            let trial = if weight > MOLE_TOLERANCE {
                remaining * mole_fraction[component] / weight
            } else {
                remaining * available_moles[component] / sum_available(available_moles, &fixed)
            };
            if trial > available_moles[component] + MOLE_TOLERANCE {
                allocated[component] = available_moles[component];
                remaining -= allocated[component];
                fixed[component] = true;
                limited = true;
            }
        }
        if !limited {
            for component in 0..count {
                if !fixed[component] {
                    allocated[component] = if weight > MOLE_TOLERANCE {
                        remaining * mole_fraction[component] / weight
                    } else {
                        remaining * available_moles[component]
                            / sum_available(available_moles, &fixed)
                    };
                }
            }
            remaining = 0.0;
        }
    }
    allocated
}

/// `sumAvailableMoles`: the unfixed components' availability, floored at the tolerance.
fn sum_available(available: &[f64], fixed: &[bool]) -> f64 {
    let sum: f64 = available
        .iter()
        .zip(fixed)
        .filter(|(_, fixed)| !**fixed)
        .map(|(moles, _)| *moles)
        .sum();
    sum.max(MOLE_TOLERANCE)
}
