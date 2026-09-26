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
