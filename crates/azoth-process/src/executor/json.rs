//! A session's result as JSON.
//!
//! **One codec, on the session's result, and not one per model.** The middleware's gap table names
//! the need as "a result as JSON" beside `toml::to_string` for the flowsheet - the two directions
//! of the same round trip - and a codec per unit operation would be twenty-six chances for the
//! twenty-six records to disagree about what a quantity looks like on the wire.
//!
//! **Every number crosses as a magnitude and a unit, in that order's own shape.** A stream's
//! record is `n`, `z`, `P`, `T`, `h`, which is the palette's declaration and not a spelling of
//! this module's, and each scalar is `{"magnitude_si": …, "unit": …}` - the shape the Python
//! bridge already transports, so a front-end that reads one reads the other.
//!
//! **A non-finite number is refused rather than written.** JSON has no `NaN` and no `Infinity`,
//! and an encoder that meets one silently writes `null` - which is the one failure this codec
//! could hide, and it is a number a caller would read as "absent" rather than "wrong".

use azoth_core::{AzothError, Result};
use serde::Serialize;
use std::collections::BTreeMap;

use crate::executor::session::Session;

/// One scalar as a magnitude in its SI unit and the unit's name.
///
/// **The magnitude is SI and the unit is what names it**, which is the contract the Python
/// bridge's `PyQty` carries. A caller that wants another unit converts; a codec that offered one
/// would be a second place for a conversion to disagree.
#[derive(Debug, Clone, Serialize)]
pub struct Quantity {
    /// The value in the unit below.
    pub magnitude_si: f64,
    /// The unit's name, as the vocabulary spells it.
    pub unit: &'static str,
}

impl Quantity {
    /// A magnitude and its unit, refusing a value JSON cannot carry.
    fn new(magnitude_si: f64, unit: &'static str) -> Result<Self> {
        if !magnitude_si.is_finite() {
            return Err(AzothError::invalid_input(
                unit,
                format!(
                    "{magnitude_si} is not a number JSON can carry - it has no NaN and no \
                     Infinity, and an encoder that met one would write `null`, which reads as \
                     absent rather than as wrong"
                ),
            ));
        }
        Ok(Self { magnitude_si, unit })
    }
}

/// A stream's record, on the fields the palette's ports declare.
#[derive(Debug, Clone, Serialize)]
pub struct StreamRecord {
    /// Molar flow, `mol/s`.
    pub n: Quantity,
    /// Composition, one entry per component.
    pub z: Vec<f64>,
    /// Pressure, `Pa`. **The declaration's own name**, which is why the field is renamed: a port
    /// writes `P`, and a reader that finds `p` here and `P` in the palette has two spellings for
    /// one field.
    #[serde(rename = "P")]
    pub p: Quantity,
    /// Temperature, `K`, under the declaration's `T`.
    #[serde(rename = "T")]
    pub temperature: Quantity,
    /// Molar enthalpy, `J/mol`.
    pub h: Quantity,
}

/// One tear's convergence, as the capture records it.
#[derive(Debug, Clone, Serialize)]
pub struct TearReport {
    /// The recycle's own name.
    pub stream: String,
    /// Passes the tear was evaluated on.
    ///
    /// **A tear switched off by the low-flow cutoff stops counting**, because the class skips an
    /// inactive unit rather than running it - so a loop the outer loop ran twice can report one
    /// pass here, and the shipped `demo.toml` does exactly that.
    pub iterations: u32,
    /// Whether `solved()` held at the last pass.
    pub solved: bool,
    /// Whether the tear was evaluated at all; `false` is `deactivateOnLowFlow`.
    ///
    /// **`solved` and `active` together are what a reader needs**, because a tear that is not
    /// active is solved by being *absent* rather than by closing - its residuals are declared zero
    /// rather than measured.
    pub active: bool,
    /// The four residuals, when a pass measured them.
    pub residuals: Option<ResidualReport>,
}

/// `Recycle`'s four residuals, in the units `solved()` compares them.
#[derive(Debug, Clone, Serialize)]
pub struct ResidualReport {
    /// `flowBalanceCheck`, in its own mixed unit - see [`crate::recycle`].
    pub flow: f64,
    /// `compositionBalanceCheck`, a sum of absolute mole-fraction differences.
    pub composition: f64,
    /// `temperatureBalanceCheck`, a sum of percentage changes.
    pub temperature: f64,
    /// `pressureBalanceCheck`, the same.
    pub pressure: f64,
}

/// A whole session.
#[derive(Debug, Clone, Serialize)]
pub struct SessionReport {
    /// The flowsheet's id.
    pub flowsheet: String,
    /// Whether every tear solved within the cap.
    pub converged: bool,
    /// Passes taken.
    pub iterations: u32,
    /// Every stream, by the endpoint that produced it.
    pub streams: BTreeMap<String, StreamRecord>,
    /// Every declared tear, in declaration order.
    pub tears: Vec<TearReport>,
}

/// A session's result, as the document a caller reads or embeds.
///
/// **The one stream encoder**, extracted from `to_json` so that the middleware's envelope can
/// carry the value rather than parse the string back. `to_json` is this value written, so the two
/// cannot disagree - and the tests that pin the published bytes are what says so.
///
/// # Errors
/// [`AzothError::InvalidInput`] if a magnitude is not finite - which a session that ran should not
/// be able to reach, and which is refused rather than written as `null`.
pub fn session_report(session: &Session) -> Result<SessionReport> {
    let mut streams = BTreeMap::new();
    for (endpoint, stream) in &session.report().streams {
        streams.insert(
            endpoint.clone(),
            StreamRecord {
                n: Quantity::new(stream.n, "mol/s")?,
                z: stream.z.clone(),
                p: Quantity::new(stream.p.value, "Pa")?,
                temperature: Quantity::new(stream.t.value, "K")?,
                h: Quantity::new(stream.h.value, "J/mol")?,
            },
        );
    }
    let tears = session
        .report()
        .tears
        .iter()
        .map(|tear| TearReport {
            stream: tear.stream.clone(),
            iterations: tear.iterations,
            solved: tear.solved,
            active: tear.active,
            residuals: tear.residuals.map(|residuals| ResidualReport {
                flow: residuals.flow,
                composition: residuals.composition,
                temperature: residuals.temperature,
                pressure: residuals.pressure,
            }),
        })
        .collect();
    Ok(SessionReport {
        flowsheet: session.flowsheet().id.clone(),
        converged: session.report().converged,
        iterations: session.report().iterations,
        streams,
        tears,
    })
}

/// Write a session's result as JSON.
///
/// # Errors
/// Whatever [`session_report`] refuses, or a write that fails - which for this shape it cannot.
pub fn to_json(session: &Session) -> Result<String> {
    serde_json::to_string(&session_report(session)?).map_err(|error| {
        AzothError::invalid_input("json", format!("the report could not be written: {error}"))
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    /// **A non-finite magnitude is refused**, which is the one failure this codec could hide.
    /// Reachable only from here: a session that ran cannot produce one, and the refusal exists so
    /// that if one ever did, it is an error rather than a `null` a caller reads as "absent".
    #[test]
    fn a_non_finite_magnitude_is_refused() {
        assert!(Quantity::new(1.0, "K").is_ok());
        assert!(Quantity::new(f64::NAN, "K").is_err());
        assert!(Quantity::new(f64::INFINITY, "K").is_err());
        assert!(Quantity::new(f64::NEG_INFINITY, "K").is_err());
        let error = Quantity::new(f64::NAN, "K").expect_err("it refuses");
        assert!(error.to_string().contains("NaN"), "{error}");
    }
}
