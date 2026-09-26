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
//!
//! **The refusal is the *declared record's*, and [`WireScalar`] is the exception that proves
//! it.** `n`, `z`, `P`, `T` and `h` are the document's own fields: a session that ran cannot
//! produce a NaN among them, so one that did is a defect and the whole answer is refused. A unit
//! op's own result is a *report* rather than a declaration - a skipped solve has no residual and a
//! divergent one has an infinity - so a non-finite magnitude there is written `null`, which is
//! what the Python half already does. Refusing there would let one bad number inside one unit op
//! take the whole envelope with it, which is the failure `run_error` exists to prevent.

use azoth_core::units::{
    MassRate, MolarEnergy, Power, Pressure, ThermalConductance, ThermodynamicTemperature, Velocity,
};
use azoth_core::{AzothError, Result, Warning};
use serde::ser::SerializeSeq;
use serde::{Serialize, Serializer};
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

/// One dimension, as the wire writes it.
///
/// **This trait is the only place a unit's name is written down.** A result's quantity field
/// names its own serialiser ([`scalar`] or [`scalars`]), which reads the name from the field's
/// *type* - so `"Pa"` appears once in the crate rather than forty times, and a field of the wrong
/// dimension cannot be published under the wrong unit.
pub trait WireScalar {
    /// The unit the vocabulary spells this dimension with.
    const UNIT: &'static str;

    /// The magnitude in that unit.
    fn magnitude(&self) -> f64;
}

impl WireScalar for Pressure {
    const UNIT: &'static str = "Pa";
    fn magnitude(&self) -> f64 {
        self.value
    }
}

impl WireScalar for ThermodynamicTemperature {
    const UNIT: &'static str = "K";
    fn magnitude(&self) -> f64 {
        self.value
    }
}

impl WireScalar for MolarEnergy {
    const UNIT: &'static str = "J/mol";
    fn magnitude(&self) -> f64 {
        self.value
    }
}

impl WireScalar for Power {
    const UNIT: &'static str = "W";
    fn magnitude(&self) -> f64 {
        self.value
    }
}

impl WireScalar for Velocity {
    const UNIT: &'static str = "m/s";
    fn magnitude(&self) -> f64 {
        self.value
    }
}

impl WireScalar for MassRate {
    const UNIT: &'static str = "kg/s";
    fn magnitude(&self) -> f64 {
        self.value
    }
}

impl WireScalar for ThermalConductance {
    const UNIT: &'static str = "W/K";
    fn magnitude(&self) -> f64 {
        self.value
    }
}

/// A `#[serde(serialize_with = "…")]` for one quantity field.
pub fn scalar<T: WireScalar, S: Serializer>(
    value: &T,
    serializer: S,
) -> std::result::Result<S::Ok, S::Error> {
    Quantity {
        magnitude_si: value.magnitude(),
        unit: T::UNIT,
    }
    .serialize(serializer)
}

/// The same for a quantity a run may not have reached, which writes `null`.
///
/// **The quantity counterpart of an `Option<f64>` result field**, and it exists because a
/// process record has the same absence a flash does: `SaftFlashResult::beta` is `None` when the
/// flash found one phase, and an exchanger's effectiveness is `None` when it was handed a pinned
/// outlet instead of a rating. A `null` says "this run did not reach it", which is a different
/// statement from a zero and the one a reader comparing two runs needs.
pub fn optional_scalar<T: WireScalar, S: Serializer>(
    value: &Option<T>,
    serializer: S,
) -> std::result::Result<S::Ok, S::Error> {
    match value {
        Some(quantity) => scalar(quantity, serializer),
        None => serializer.serialize_none(),
    }
}

/// The same for a vector of them, one entry per element.
pub fn scalars<T: WireScalar, S: Serializer>(
    values: &[T],
    serializer: S,
) -> std::result::Result<S::Ok, S::Error> {
    let mut sequence = serializer.serialize_seq(Some(values.len()))?;
    for value in values {
        sequence.serialize_element(&Quantity {
            magnitude_si: value.magnitude(),
            unit: T::UNIT,
        })?;
    }
    sequence.end()
}

/// One warning, as the wire writes it.
#[derive(Debug, Clone, Serialize)]
pub struct WireWarning<'a> {
    /// The code's stable string form, e.g. `SOLVER_NOT_CONVERGED`.
    pub code: &'static str,
    /// What happened, in the library's own sentence.
    pub message: &'a str,
    /// The input or output it is about, when it is about one.
    pub field: Option<&'a str>,
}

/// A `#[serde(serialize_with = "…")]` for a result's `warnings`.
pub fn warnings<S: Serializer>(
    values: &[Warning],
    serializer: S,
) -> std::result::Result<S::Ok, S::Error> {
    let mut sequence = serializer.serialize_seq(Some(values.len()))?;
    for warning in values {
        sequence.serialize_element(&WireWarning {
            code: warning.code.as_str(),
            message: &warning.message,
            field: warning.field.as_deref(),
        })?;
    }
    sequence.end()
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
    /// Mass flow, `kg/s`.
    ///
    /// **`null` where the library cannot weigh the fluid**, which is [`crate::stream::Stream::molar_mass`]'s
    /// own refusal: a component built from critical constants alone carries no molar mass, and a
    /// record that defaulted it to zero would report a mass flow of zero for a real stream.
    pub mass_flow: Option<Quantity>,
    /// The mixture's molar mass, `kg/mol`, or `null` for the same reason.
    pub molar_mass: Option<Quantity>,
    /// The vapour fraction, or `null` where nothing computed one.
    ///
    /// **This is the one field a stream's own record cannot derive.** `n`, `z`, `P`, `T` and `h`
    /// are the record; the vapour fraction is what the *flash* said, and only the unit operation
    /// that ran one knows it. [`crate::stream::Stream`] carries it so that a phase a kernel knew
    /// is not thrown away when the stream is built, and `null` is the honest answer for a stream
    /// nothing flashed.
    pub vapour_fraction: Option<f64>,
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
    /// Every unit op's own answer, by instance id - **the numbers on no outlet stream**.
    ///
    /// **Absent for a unit op whose whole answer is its streams.** A mixer's result *is* its
    /// outlet and a separator's *are* its two, so those carry no entry here rather than a
    /// restatement of `streams` under a second key; what appears is what a kernel computed
    /// *beside* its outlets - a duty, a tray profile, a conversion, a convergence.
    ///
    /// Each entry is the unit op's own registered result, so its keys are that model's
    /// `CalcResult::FIELDS` - the same names the Python half publishes for the same operation.
    pub results: BTreeMap<String, serde_json::Value>,
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
        // Once, because `mass_flow` is `n` times it and a second call would resolve the fluid a
        // second time for the same two numbers.
        let molar_mass = stream.molar_mass().ok().map(|mass| mass.value);
        streams.insert(
            endpoint.clone(),
            StreamRecord {
                n: Quantity::new(stream.n, "mol/s")?,
                z: stream.z.clone(),
                p: Quantity::new(stream.p.value, "Pa")?,
                temperature: Quantity::new(stream.t.value, "K")?,
                h: Quantity::new(stream.h.value, "J/mol")?,
                mass_flow: molar_mass.and_then(|mass| finite(stream.n * mass, "kg/s")),
                molar_mass: molar_mass.and_then(|mass| finite(mass, "kg/mol")),
                vapour_fraction: stream.vapour_fraction.filter(|beta| beta.is_finite()),
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
        results: session.report().results.clone(),
        tears,
    })
}

/// A magnitude and its unit, or `None` where the number is one JSON cannot carry.
///
/// The record's own five fields refuse such a number rather than dropping it; these three are
/// enrichments a run may not have, so an absent one is an absent column rather than a failed
/// envelope.
fn finite(magnitude_si: f64, unit: &'static str) -> Option<Quantity> {
    if magnitude_si.is_finite() {
        Some(Quantity { magnitude_si, unit })
    } else {
        None
    }
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
    use azoth_core::units::{pascals, watts};

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

    /// A quantity field crosses as a magnitude and the unit **its type** names, so a serialiser
    /// cannot print a pressure under the wrong name.
    #[test]
    fn a_quantity_crosses_under_its_own_dimensions_unit() {
        let written = serde_json::to_value(Quantity {
            magnitude_si: 1950000.0,
            unit: Pressure::UNIT,
        })
        .expect("it writes");
        assert_eq!(written["unit"], "Pa");
        assert_eq!(written["magnitude_si"], 1950000.0);
        assert_eq!(ThermodynamicTemperature::UNIT, "K");
        assert_eq!(MolarEnergy::UNIT, "J/mol");
        assert_eq!(Power::UNIT, "W");
        assert_eq!(watts(1.0).magnitude(), 1.0);
        assert_eq!(pascals(1.0).magnitude(), 1.0);
    }

    /// **A result's non-finite number is written `null` and the record is still refused**, which
    /// is the asymmetry the module doc states: a report may have a gap where a declaration may
    /// not.
    #[test]
    fn a_report_writes_a_gap_where_a_declaration_refuses_one() {
        let magnitude = f64::NAN;
        assert!(Quantity::new(magnitude, "W").is_err());
        assert!(finite(magnitude, "W").is_none());
        assert!(finite(1.5, "W").is_some());
    }
}
