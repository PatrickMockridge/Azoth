//! The Murphree tray efficiency: `SimpleTray.setMurphreeEfficiency`, and the correction
//! `DistillationColumn.applyMurphreeCorrection` applies to a stage's vapour after it has run.
//!
//! **The correction blends the vapour leaving a stage towards the one entering it**:
//! `y_out = y_in + E(y_eq - y_in)`, clamped at zero and renormalised, at the *same* total
//! moles the flash found. The liquid is untouched - it keeps the equilibrium flash - which is
//! the standard post-correction form for a bubble-point sequential method.
//!
//! **It does not conserve components, and a caller has to know that.** What a stage hands up
//! carries the flash's *moles* at the blend's composition, so the stage's outlets no longer
//! carry its inlets and the interior imbalances do not cancel - the column's products miss the
//! feed by `5.6e-2` relative on the captured row. `DistillationColumn` answers that with
//! `updateProductsFromExternalComponentBalance`, which rescales the two products per component
//! against the feed; the kernel's `reconcile_products` is the same arithmetic, and the residual
//! it rescales away is the one it reports.
//!
//! **`y_in` is the stage below's *equilibrium* vapour and not its corrected outlet.** The class
//! reads `trays.get(i - 1).getThermoSystem().getPhase(0)`, and `setCachedGasOutStream` writes a
//! *separate* stream object without touching that system - so the stage above *receives* the
//! corrected vapour while the correction *compares against* the uncorrected one. The two are
//! different numbers and a port that reads the corrected one converges to a different profile,
//! so [`super::Network`] carries both.
//!
//! **Three stages are never corrected**, in the class's own order of tests: an ideal stage
//! (`E >= 1 - 1e-10`), the reboiler at index zero, and the condenser when there is one. The
//! capture's profile shows the shape of it - on the binary column the four interior stages move
//! and the two ends do not.

use azoth_core::Result;

use crate::stream::Stream;

/// The class's own ideal-stage tolerance: an efficiency within `1e-10` of one is no correction.
const IDEAL_TOLERANCE: f64 = 1.0e-10;

/// The floor below which the class does not renormalise the corrected composition.
const RENORMALISE_FLOOR: f64 = 1.0e-15;

/// **`DistillationColumn`'s two efficiency fields and the resolution between them.**
///
/// The class carries `murphreeEfficiency`, a column-wide scalar defaulting to `1.0`, and
/// `perStageMurphreeEfficiency`, a nullable array. `getEffectiveMurphreeEfficiency(stage)` is the
/// whole of the rule: an in-range entry that is not `NaN` wins, and everything else falls through
/// to the scalar. So **two levels, not four** - the per-component ones are
/// `AbsorptionColumn`'s, and arrive with its own correction.
///
/// `NaN` is the fall-through sentinel rather than an absence, which is why a per-stage vector
/// cannot simply be shorter than the column: a stage is either overridden or it is not, and
/// `[0.6, NaN, 0.9]` says so. `specs/models/process/absorption_column.toml`'s `tray_temperatures`
/// already spells one property this way.
#[derive(Debug, Clone, PartialEq)]
pub struct Murphree {
    /// `murphreeEfficiency`, which `setMurphreeEfficiency(double)` sets. The class clamps a
    /// request into `[0, 1]` rather than refusing it.
    pub column_wide: f64,
    /// `perStageMurphreeEfficiency`, which `setMurphreeEfficiency(int, double)` allocates and
    /// `setMurphreeEfficiencies(double[])` fills. `None` is the class's own initialiser.
    pub per_stage: Option<Vec<f64>>,
}

impl Murphree {
    /// The class's own initialiser: an ideal stage, column-wide, with no overrides.
    #[must_use]
    pub fn ideal() -> Self {
        Self {
            column_wide: 1.0,
            per_stage: None,
        }
    }

    /// The column-wide value alone, with no overrides - what a palette form can state.
    #[must_use]
    pub fn from_column_wide(efficiency: f64) -> Self {
        Self {
            column_wide: Self::clamp(efficiency),
            per_stage: None,
        }
    }

    /// `clampMurphreeEfficiency`: `max(0.0, min(1.0, efficiency))`.
    ///
    /// **The class clamps a request and does not refuse it**, which is the opposite of this
    /// library's rule for a fraction - so the clamp is a transcription to name, and the one place
    /// it is applied is at the declaration's edge.
    #[must_use]
    pub fn clamp(efficiency: f64) -> f64 {
        efficiency.clamp(0.0, 1.0) // numerics-ok: DistillationColumn.java:14661
    }

    /// `getEffectiveMurphreeEfficiency(stage)`: the override when the stage has a finite one,
    /// otherwise the column-wide value.
    #[must_use]
    pub fn resolve(&self, index: usize) -> f64 {
        if let Some(per_stage) = self.per_stage.as_ref() {
            if let Some(value) = per_stage.get(index) {
                if !value.is_nan() {
                    return *value;
                }
            }
        }
        self.column_wide
    }

    /// The stage count every override is stated over, or `None` where none is.
    ///
    /// `setMurphreeEfficiencies(double[])` refuses an array whose length is not the column's
    /// stage count, and so does this - a caller who meant one stage has `NaN` for the rest.
    ///
    /// # Errors
    /// [`azoth_core::AzothError::InvalidInput`] for a length the stage count does not match.
    pub fn checked(&self, tray_count: usize) -> Result<&Self> {
        if let Some(per_stage) = self.per_stage.as_ref() {
            if per_stage.len() != tray_count {
                return Err(azoth_core::AzothError::invalid_input(
                    "tray_murphree_efficiency",
                    format!(
                        "{} override(s) for a column of {} stage(s): \
                         `DistillationColumn.setMurphreeEfficiencies` refuses an array whose \
                         length is not the stage count, and a stage that states no override is \
                         written `NaN` rather than left out",
                        per_stage.len(),
                        tray_count
                    ),
                ));
            }
        }
        Ok(self)
    }
}

/// Whether `DistillationColumn.applyMurphreeCorrection` would correct the stage at `index`.
///
/// The three tests are the class's, in its order: an ideal stage returns first, then the
/// reboiler, then the condenser. `tray_count` counts the ends, so the condenser is
/// `tray_count - 1`.
#[must_use]
pub fn corrects(index: usize, tray_count: usize, has_condenser: bool, efficiency: f64) -> bool {
    efficiency < 1.0 - IDEAL_TOLERANCE && index > 0 && !(has_condenser && index >= tray_count - 1)
}

/// **`DistillationColumn.applyMurphreeCorrection`**: the vapour leaving one stage, blended.
///
/// The four pieces are the class's own reads. `here` is the stage's flash - its vapour is
/// `y_eq` and its *phase count* is what the class tests with `getNumberOfPhases() < 2`; `below`
/// is the stage beneath, read the same way. A stage whose flash found one phase is left alone,
/// and so is one whose neighbour below did, which is why the caller passes both outlets rather
/// than just the vapour: **the guard is that the system holds two phases**, and a `Stream` that
/// is absent is the only way this port can say a phase is not there.
///
/// `None` where the class would return without touching anything.
///
/// # Errors
/// Whatever building the corrected vapour's state refuses.
pub fn correct_vapour(
    here: (Option<&Stream>, Option<&Stream>),
    below: (Option<&Stream>, Option<&Stream>),
    efficiency: f64,
) -> Result<Option<Stream>> {
    let (Some(equilibrium), Some(_)) = here else {
        return Ok(None);
    };
    let (Some(inlet), Some(_)) = below else {
        return Ok(None);
    };
    if equilibrium.components != inlet.components {
        return Err(azoth_core::AzothError::invalid_input(
            "murphree_efficiency",
            "the stage and the stage below it do not carry the same substances in the same \
             order, so the correction has no basis for pairing one's components with the \
             other's: `DistillationColumn.applyMurphreeCorrection` indexes both by position",
        ));
    }

    let mut actual: Vec<f64> = equilibrium
        .z
        .iter()
        .zip(&inlet.z)
        .map(|(equilibrium, inlet)| (inlet + efficiency * (equilibrium - inlet)).max(0.0))
        .collect();
    let sum: f64 = actual.iter().sum();
    if sum > RENORMALISE_FLOOR {
        for fraction in &mut actual {
            *fraction /= sum;
        }
    }

    // **The total moles are the flash's own and the composition is the blended one.** The class
    // builds the corrected stream from `fluid.phaseToSystem(0)`, replaces every component's
    // fraction and mole count, and re-initialises - so what the stage above *receives* is a
    // stream whose overall composition is the blend, at the stage's own temperature and
    // pressure. Its own phase 0 after that re-initialisation is the equilibrium vapour again,
    // which is why a reader must take the stream's overall and not one of its phases.
    Ok(Some(Stream::from_pt(
        equilibrium.components.clone(),
        actual,
        equilibrium.n,
        equilibrium.p,
        equilibrium.t,
    )?))
}

#[cfg(test)]
mod tests {
    use super::Murphree;

    #[test]
    fn the_override_wins_where_it_is_finite() {
        let efficiency = Murphree {
            column_wide: 0.6,
            per_stage: Some(vec![f64::NAN, f64::NAN, f64::NAN, 0.85, f64::NAN, f64::NAN]),
        };
        assert_eq!(efficiency.resolve(3), 0.85);
        assert_eq!(efficiency.resolve(1), 0.6);
        assert_eq!(efficiency.resolve(5), 0.6);
        // An index past the array falls through rather than panicking, which is the class's own
        // `stage < length` test.
        assert_eq!(efficiency.resolve(9), 0.6);
    }

    #[test]
    fn the_class_clamps_a_request_rather_than_refusing_it() {
        assert_eq!(Murphree::clamp(1.4), 1.0);
        assert_eq!(Murphree::clamp(-0.2), 0.0);
        assert_eq!(Murphree::clamp(0.6), 0.6);
    }
}
