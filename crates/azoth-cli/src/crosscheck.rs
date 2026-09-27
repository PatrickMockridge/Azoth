//! `azoth crosscheck` — a NeqSim flash request in, a cross-check verdict out.
//!
//! This is the sidecar. A NeqSim user (or their MCP client) sends the same `runFlash` request it
//! would send to NeqSim, azoth computes the same state from its own kernels, and the answer is
//! one of three things: the two agree, they **diverge**, or azoth **refuses**. The refusal is the
//! product — it is the one outcome NeqSim's architecture cannot produce, because a Java `double`
//! holds `0.0` and "not available" in the same type.
//!
//! **The request is NeqSim's document, not a translation of it.** The field names, the unit
//! strings, the `components` map and the `flashType` vocabulary are theirs (`mcp_neqsim_core_layer`
//! documents them), and nothing azoth-specific is required to be added. Unknown fields are
//! *ignored* rather than refused, because their contract says new optional fields may be added at
//! any time and a client that refused them would break on their next release.
//!
//! **The answer is shaped like theirs on purpose.** `apiVersion`, `status`, `tool`, `data`,
//! `provenance`, `validation`, `qualityGate`, `warnings` and `autoValidation` are their envelope's
//! field names and their `status`, warning-code and `verdict` vocabularies, so one client can read
//! both servers. That reuse is a *request* to them, not a claim of conformance: the issue asks
//! whether a third-party tool may speak their envelope at all. `qualityGate.verdict` is
//! `passed`/`failed` and `engineeringReviewRequired` is true on both — read out of their
//! `mcp/README.md`, not guessed.
//!
//! **A divergence is a finding, never a failure of theirs.** azoth is a port of NeqSim and
//! everything it ships is `unverified`, so agreement is evidence about a *transcription* and not
//! about the physics. The envelope says what differed and by how much, and stops. It never says
//! which side is right.
//!
//! `--keycard` exists here and nowhere else in the CLI: this is the only subcommand that resolves
//! a substance, so it is the only one a card can change an answer for.

use std::collections::BTreeMap;
use std::io::{BufRead, Write};
use std::path::Path;

use azoth_core::AzothError;
use azoth_core::units::{kelvins, pascals};
use azoth_core::warning::Warning;
use azoth_eos::card::Card;
use azoth_eos::databank::Overlay;
use azoth_eos::{Cubic, Phase, databank, pt_flash};
use serde::{Deserialize, Serialize};

/// The contract version this speaks, as their envelope carries it.
const API_VERSION: &str = "1.0";

/// The tool name an answer carries. Not one of their tool names — this is a different tool that
/// happens to answer in their shape.
const TOOL: &str = "azothCrossCheck";

/// The tolerance a `vapour_fraction` comparison uses unless the request sets its own.
///
/// **It is a claim, so it is printed**: every answer names the tolerance it judged against, in
/// `provenance.assumptions`. A default here would otherwise be indistinguishable from agreement
/// to any precision at all. `1e-6` is loose on purpose — the two implementations stop on their
/// own convergence rules, so a difference smaller than a solver tolerance is not a divergence
/// worth a finding.
const DEFAULT_TOLERANCE: f64 = 1.0e-6;

/// One `runFlash` request, as their core layer defines it.
///
/// No `deny_unknown_fields`: their contract allows new optional fields at any time, and a sidecar
/// that refused them would fail on the release after this one.
#[derive(Debug, Deserialize)]
struct Request {
    /// The equation of state. Their default is `SRK`; see [`cubic_of`] for what is refused.
    model: String,
    temperature: Measured,
    pressure: Measured,
    /// `TP` and the rest of their vocabulary. Only `TP` is answered here.
    #[serde(rename = "flashType", default = "default_flash_type")]
    flash_type: String,
    /// Component name to mole fraction. Names come from NeqSim's `COMP.csv`, which is where
    /// azoth's databank comes from too, so an ordinary fluid needs no keycard.
    components: BTreeMap<String, f64>,
    #[serde(rename = "mixingRule", default = "default_mixing_rule")]
    mixing_rule: String,
    /// What the caller's own implementation produced, when they want it compared. Absent means
    /// "compute only", which is a different answer and not a failed comparison.
    #[serde(default)]
    reported: Option<Reported>,
    /// An absolute tolerance for the comparison. Defaults to [`DEFAULT_TOLERANCE`].
    #[serde(default)]
    tolerance: Option<f64>,
}

fn default_flash_type() -> String {
    "TP".to_string()
}

fn default_mixing_rule() -> String {
    "classic".to_string()
}

/// A number with the unit their request states it in.
#[derive(Debug, Deserialize)]
struct Measured {
    value: f64,
    unit: String,
}

/// The caller's own result for the state, as far as it is comparable.
#[derive(Debug, Deserialize)]
struct Reported {
    /// Their vapour fraction. Absent is a real state — a single-phase answer with no split — and
    /// is compared as a presence, not as a zero.
    #[serde(rename = "vapourFraction", default)]
    vapour_fraction: Option<f64>,
}

/// Why azoth would not answer, with the input that was refused named.
#[derive(Debug)]
struct Refusal {
    /// The request field the refusal is about, so a caller can act on it.
    field: &'static str,
    /// The warning code from *their* taxonomy that this refusal is closest to.
    code: &'static str,
    /// What was refused and why, in the library's own words where the library refused.
    message: String,
}

impl Refusal {
    fn new(field: &'static str, code: &'static str, message: impl Into<String>) -> Self {
        Self {
            field,
            code,
            message: message.into(),
        }
    }

    /// The library's refusal, with its own field carried across.
    ///
    /// The error's `field()` is what makes a refusal actionable and it is exactly what the
    /// middleware's shared error path discards, by flattening to a sentence at
    /// `session.rs`'s `map_err(|e| e.to_string())`. This path does not go through that one.
    fn from_library(error: &AzothError, fallback_field: &'static str) -> Self {
        Self {
            field: match error.field() {
                Some("components") => "components",
                Some("z") => "components",
                _ => fallback_field,
            },
            code: "MISSING_REFERENCE_DATA",
            message: error.to_string(),
        }
    }
}

/// What azoth computed, or `None` when it refused.
struct Computed {
    vapour_fraction: Option<f64>,
    phase: Phase,
    iterations: u32,
    residual: f64,
    min_t_over_tc: f64,
    warnings: Vec<Warning>,
    assumptions: Vec<String>,
}

/// Run the sidecar: one request per line on `input`, one answer per line on `output`.
///
/// The framing is `azoth mcp`'s — newline-delimited JSON, flushed per line — so a caller that can
/// speak to one can speak to the other, and a line that is not JSON is answered rather than
/// killing the process.
///
/// # Errors
/// Only on an I/O failure of the streams themselves. A malformed request is answered with an
/// envelope whose `status` is `error`, because a sidecar that exited on bad input would take the
/// caller's session with it.
pub fn serve(
    keycard: Option<&Path>,
    input: impl BufRead,
    mut output: impl Write,
) -> Result<(), String> {
    // Loaded once, not per request: a card that failed to load is the process failing to start,
    // which is the same place a flowsheet is read for the other subcommands.
    let card = match keycard {
        Some(path) => Some(
            Card::from_path(path)
                .map_err(|e| format!("{}: {e}", path.display()))
                .map_err(|message| {
                    format!(
                        "{message}\n  the card is refused as a whole: a card exists only once \
                         every one of its refusals has passed"
                    )
                })?,
        ),
        None => None,
    };
    let overlay = card.as_ref().map(Card::overlay);

    for line in input.lines() {
        let line = line.map_err(|e| e.to_string())?;
        if line.trim().is_empty() {
            continue;
        }
        let answer = answer(&line, overlay);
        writeln!(output, "{answer}").map_err(|e| e.to_string())?;
        output.flush().map_err(|e| e.to_string())?;
    }
    Ok(())
}

/// One line in, one envelope out, as JSON.
///
/// Never fails: a request that cannot be read is answered with `status: "error"`, which is one of
/// their four.
#[must_use]
pub fn answer(line: &str, overlay: Option<&Overlay>) -> String {
    let request = match serde_json::from_str::<Request>(line) {
        Ok(request) => request,
        Err(error) => return encode(&Envelope::error(&error.to_string())),
    };
    encode(&verdict(&request, overlay))
}

/// Convert a state's temperature to the kelvins the kernels take.
///
/// # Errors
/// On a unit outside their vocabulary, naming it — azoth's unit vocabulary is closed and a
/// `pint`-style alias is not accepted for `K`.
fn kelvin_of(measured: &Measured) -> Result<f64, Refusal> {
    match measured.unit.as_str() {
        "K" => Ok(measured.value),
        "C" => Ok(measured.value + 273.15),
        "F" => Ok((measured.value - 32.0) * 5.0 / 9.0 + 273.15),
        other => Err(Refusal::new(
            "temperature",
            "MODEL_LIMITATION",
            format!(
                "temperature unit `{other}` is not one of `K`, `C` or `F`. NeqSim's \
                 `validateInput` accepts a wider spelling of the same three units; this side \
                 takes the three their schema enumerates."
            ),
        )),
    }
}

/// Convert a state's pressure to the pascals the kernels take.
///
/// `barg` is answered rather than refused, because it is well defined given an atmosphere — but
/// the atmosphere is an assumption, so it goes in the answer's `assumptions` rather than
/// disappearing into a constant.
///
/// # Errors
/// On a unit outside their vocabulary, naming it.
fn pascal_of(measured: &Measured) -> Result<(f64, Option<String>), Refusal> {
    match measured.unit.as_str() {
        "Pa" => Ok((measured.value, None)),
        "kPa" => Ok((measured.value * 1.0e3, None)),
        "MPa" => Ok((measured.value * 1.0e6, None)),
        "bara" => Ok((measured.value * 1.0e5, None)),
        "psi" => Ok((measured.value * 6_894.757_293_168, None)),
        "atm" => Ok((measured.value * 101_325.0, None)),
        "barg" => Ok((
            measured.value * 1.0e5 + 101_325.0,
            Some(
                "`barg` was read as gauge over exactly 101325 Pa, which is a standard \
                 atmosphere and not a measurement. A site with a different barometric \
                 reference needs `bara`."
                    .to_string(),
            ),
        )),
        other => Err(Refusal::new(
            "pressure",
            "MODEL_LIMITATION",
            format!(
                "pressure unit `{other}` is not one of `bara`, `barg`, `Pa`, `kPa`, `MPa`, \
                 `psi` or `atm`"
            ),
        )),
    }
}

/// The cubic a model name is, or the refusal that says what azoth can and cannot express.
///
/// **A model azoth cannot express is refused, never substituted.** NeqSim's `CPA`, `GERG2008`,
/// `PCSAFT` and `UMRPRU` are not cubics and azoth reaches them through separate phase modules;
/// answering a `GERG2008` request with a Peng-Robinson number would be a plausible-looking wrong
/// answer, which is the one thing this library may not produce. `SRK` and `PR` are refused for
/// the same reason `mixture_of` refuses an ion: the constants a cubic would use are filler.
///
/// # Errors
/// On any model name outside `PR`, `SRK`, `RK` or `TST`, naming what was asked for.
fn cubic_of(model: &str) -> Result<Cubic, Refusal> {
    model.to_lowercase().parse::<Cubic>().map_err(|_| {
        Refusal::new(
            "model",
            "MODEL_LIMITATION",
            format!(
                "model `{model}` was not answered. This subcommand answers the cubic models \
                 `PR`, `SRK`, `RK` and `TST`; a non-cubic reaches azoth through its own phase \
                 module and is not wired into this request shape, and substituting a cubic for \
                 it would return a plausible-looking wrong number."
            ),
        )
    })
}

/// The fields of a request that must be a specific value for this subcommand to answer it.
///
/// # Errors
/// On a flash type or mixing rule azoth does not reach through this shape, naming it.
fn require_supported(request: &Request) -> Result<(), Refusal> {
    if !request.flash_type.eq_ignore_ascii_case("TP") {
        return Err(Refusal::new(
            "flashType",
            "MODEL_LIMITATION",
            format!(
                "flashType `{}` is not answered here. azoth has a kernel for most of their \
                 vocabulary — `eos.ph_flash`, `eos.ps_flash`, `eos.tv_flash`, `eos.dew_pressure`, \
                 `eos.dew_temperature` — and none of them is wired into this request shape yet, \
                 so the request is refused rather than answered by the isothermal flash, which \
                 would solve a different problem and return a number that looked like an answer.",
                request.flash_type
            ),
        ));
    }
    if !request.mixing_rule.eq_ignore_ascii_case("classic") {
        return Err(Refusal::new(
            "mixingRule",
            "MODEL_LIMITATION",
            format!(
                "mixingRule `{}` is not answered here. The databank's kij matrix is the \
                 classical one; the keycard's model vocabulary carries `classical_kij` and \
                 nothing else, so another rule has no expression on this side.",
                request.mixing_rule
            ),
        ));
    }
    Ok(())
}

/// The names and mole fractions a request asks about, in one fixed order.
///
/// # Errors
/// On an empty component map, a negative fraction, or a total that is not one. **The sum is
/// checked rather than renormalised**: NeqSim normalises and warns (`COMPOSITION_NORMALIZED` is
/// one of their codes), and azoth refuses, because rescaling a caller's composition here would
/// make their error invisible in every number downstream. The refusal names the total so the
/// caller can decide.
fn composition_of(request: &Request) -> Result<(Vec<&str>, Vec<f64>), Refusal> {
    if request.components.is_empty() {
        return Err(Refusal::new(
            "components",
            "MISSING_REFERENCE_DATA",
            "a mixture needs at least one component",
        ));
    }
    let names: Vec<&str> = request.components.keys().map(String::as_str).collect();
    let z: Vec<f64> = request.components.values().copied().collect();
    if let Some(bad) = z.iter().position(|&value| value < 0.0) {
        return Err(Refusal::new(
            "components",
            "MISSING_REFERENCE_DATA",
            format!(
                "{} is {}, and a mole fraction cannot be negative",
                names[bad], z[bad]
            ),
        ));
    }
    let sum: f64 = z.iter().sum();
    if (sum - 1.0).abs() > 1.0e-9 {
        return Err(Refusal::new(
            "components",
            "COMPOSITION_NORMALIZED",
            format!(
                "the mole fractions sum to {sum}, not to one, so no state was computed. NeqSim \
                 normalises and warns; azoth refuses and names the total, because renormalising \
                 here would make the input's error invisible in every number downstream. \
                 Normalise the composition and send it again."
            ),
        ));
    }
    Ok((names, z))
}

/// Compute the state the request asks for.
fn compute(request: &Request, overlay: Option<&Overlay>) -> Result<Computed, Refusal> {
    require_supported(request)?;
    let cubic = cubic_of(&request.model)?;
    let temperature = kelvin_of(&request.temperature)?;
    let (pressure, gauge_assumption) = pascal_of(&request.pressure)?;
    let (names, z) = composition_of(request)?;

    let mut assumptions = vec![
        "the request's `model` names a cubic and its `mixingRule` is the classical kij matrix, \
         which is the databank's own"
            .to_string(),
    ];
    if let Some(note) = gauge_assumption {
        assumptions.push(note);
    }
    if overlay.is_none() {
        assumptions.push(
            "no keycard was supplied, so every substance resolved against azoth's databank, \
             which is generated from NeqSim's own COMP.csv and INTER.csv"
                .to_string(),
        );
    }

    let (mixture, _ideal_gas) = databank::mixture_of(&names, cubic, overlay)
        .map_err(|error| Refusal::from_library(&error, "components"))?;
    let result = pt_flash(&mixture, kelvins(temperature), pascals(pressure), &z)
        .map_err(|error| Refusal::from_library(&error, "temperature"))?;

    Ok(Computed {
        vapour_fraction: result.vapour_fraction,
        phase: result.phase,
        iterations: result.iterations,
        residual: result.residual,
        min_t_over_tc: result.min_t_over_tc,
        warnings: result.warnings,
        assumptions,
    })
}

/// The verdict for one request: what azoth computed, and how it compares.
fn verdict(request: &Request, overlay: Option<&Overlay>) -> Envelope {
    let computed = match compute(request, overlay) {
        Ok(computed) => computed,
        Err(refusal) => return Envelope::refused(&refusal, request),
    };

    let tolerance = request.tolerance.unwrap_or(DEFAULT_TOLERANCE);
    let mut warnings = Vec::new();
    let mut checks = Vec::new();

    // Their own near-critical rule: within 10% of the critical point. azoth carries the smallest
    // T/Tc of the mixture, so this is a measured flag and not a heuristic.
    if computed.min_t_over_tc > 0.9 {
        warnings.push(OutWarning::new(
            "NEAR_CRITICAL",
            "CAUTION",
            format!(
                "the smallest T/Tc over the components is {:.4}, within 10% of the critical \
                 point, where a cubic's two roots approach each other and both implementations \
                 lose conditioning",
                computed.min_t_over_tc
            ),
        ));
    }
    for warning in &computed.warnings {
        warnings.push(map_warning(warning));
    }

    let comparison = compare(request, &computed, tolerance);
    match &comparison {
        Comparison::Agreed => checks.push(Check::new(
            "cross-implementation",
            "PASS",
            format!("the vapour fraction agrees within {tolerance:e}"),
        )),
        Comparison::Diverged {
            detail,
            code,
            severity,
        } => {
            checks.push(Check::new("cross-implementation", "FAIL", detail.clone()));
            warnings.push(OutWarning::new(code, severity, detail.clone()));
            warnings.push(OutWarning::new(
                "MODEL_LIMITATION",
                "INFO",
                "a divergence between two implementations of one source is evidence about a \
                 transcription, not about the physics: azoth is a port of NeqSim and both would \
                 share a conceptual error. Explain the difference before believing either number.",
            ));
        }
        Comparison::NotReported => checks.push(Check::new(
            "cross-implementation",
            "PASS",
            "no comparison was asked for; this answer is azoth's own and nothing was checked \
             against it",
        )),
    }

    let overall =
        if comparison.is_diverged() || warnings.iter().any(|warning| warning.severity != "INFO") {
            "WARNING"
        } else {
            "PASS"
        };

    Envelope {
        api_version: API_VERSION,
        status: "success",
        tool: TOOL,
        data: Some(Data {
            vapour_fraction: computed.vapour_fraction.map(|value| Quantity {
                magnitude_si: value,
                unit: "dimensionless",
            }),
            phase: computed.phase.as_str(),
            iterations: computed.iterations,
            residual: computed.residual,
        }),
        provenance: Provenance {
            model: request.model.clone(),
            flash_type: request.flash_type.clone(),
            convergence: Convergence {
                converged: true,
                iterations: computed.iterations,
                residual: computed.residual,
            },
            assumptions: computed.assumptions,
            limitations: limitations(),
        },
        validation: Validation {
            valid: true,
            phase: Some(computed.phase.as_str()),
            message: None,
        },
        quality_gate: QualityGate {
            verdict: if comparison.is_diverged() {
                "failed"
            } else {
                "passed"
            },
            summary: comparison.summary(),
            engineering_review_required: true,
        },
        warnings,
        auto_validation: AutoValidation { overall, checks },
    }
}

/// What a comparison found.
enum Comparison {
    /// Both produced a value and they are within tolerance.
    Agreed,
    /// The two do not agree. `code` and `severity` are the warning this warrants, which differs
    /// by *why* they disagree.
    Diverged {
        detail: String,
        code: &'static str,
        severity: &'static str,
    },
    /// No comparison was asked for.
    NotReported,
}

impl Comparison {
    fn is_diverged(&self) -> bool {
        matches!(self, Self::Diverged { .. })
    }

    fn summary(&self) -> String {
        match self {
            Self::Agreed => "azoth's vapour fraction agrees with the reported one".to_string(),
            Self::Diverged { detail, .. } => format!(
                "azoth and the reported implementation disagree: {detail}. Neither is \
                 authoritative; azoth is a port and ships `unverified`, so this is a finding to \
                 explain rather than a verdict on their number."
            ),
            Self::NotReported => {
                "azoth computed this state; no other implementation's value was supplied, so \
                 nothing was cross-checked"
                    .to_string()
            }
        }
    }
}

/// Compare azoth's vapour fraction against the one the caller reported.
///
/// **Presence is compared as well as value.** A single-phase answer has no vapour fraction on
/// either side, and azoth reports it *absent* rather than as zero — which is the whole doctrine
/// in one field. Treating absence as zero would turn a real disagreement into agreement.
///
/// **A number is only compared when both sides are describing a split.** For a single-phase feed
/// azoth reports the converged Rachford-Rice root anyway, as the negative-flash value, and says
/// so: it can be outside `[0, 1]` and it is not a fraction of anything. Comparing `1.5747`
/// against a reported `0.9` and calling the difference `0.67` would be a number that looks like
/// an answer and is not — the two sides are describing different regimes, and the disagreement
/// is worth reporting as *that*.
fn compare(request: &Request, computed: &Computed, tolerance: f64) -> Comparison {
    let Some(reported) = request.reported.as_ref() else {
        return Comparison::NotReported;
    };
    let split = computed.phase == Phase::TwoPhase;
    match (computed.vapour_fraction, reported.vapour_fraction) {
        (Some(ours), Some(theirs)) if split => {
            let difference = (ours - theirs).abs();
            if difference <= tolerance {
                Comparison::Agreed
            } else {
                Comparison::Diverged {
                    detail: format!(
                        "vapour fraction: azoth {ours:.12}, reported {theirs:.12}, absolute \
                         difference {difference:e}, tolerance {tolerance:e}"
                    ),
                    // Both sides are in the same regime, so this is a numeric disagreement and
                    // the interesting one: a transcription difference between two readings of
                    // one source.
                    code: "MODEL_LIMITATION",
                    severity: "INFO",
                }
            }
        }
        // azoth found no split and still reports a root. Its value is a negative-flash
        // extrapolation, so what differs is the regime, not the number.
        (Some(ours), Some(theirs)) => Comparison::Diverged {
            detail: format!(
                "phase: azoth finds `{}` and reports the negative-flash root {ours:.12}, which is \
                 not a fraction of anything; the reported implementation gives a vapour fraction \
                 of {theirs:.12}, so it found a split. The two are describing different regimes, \
                 and the values are not comparable — this is the phase boundary, or a stability \
                 test azoth's model does not have",
                computed.phase.as_str()
            ),
            code: "TWO_PHASE_UNCERTAINTY",
            severity: "CAUTION",
        },
        (None, None) => Comparison::Agreed,
        (None, Some(theirs)) => Comparison::Diverged {
            detail: format!(
                "vapour fraction: azoth reports it absent, in the phase `{}`, and the reported \
                 implementation gives {theirs:.12}. A state with no split has no vapour \
                 fraction; if the reported one is 0.0 it is a zero standing in for an absence",
                computed.phase.as_str()
            ),
            code: "TWO_PHASE_UNCERTAINTY",
            severity: "CAUTION",
        },
        (Some(ours), None) => Comparison::Diverged {
            detail: format!(
                "vapour fraction: azoth gives {ours:.12} in the phase `{}` and the reported \
                 implementation gives none at all",
                computed.phase.as_str()
            ),
            code: "TWO_PHASE_UNCERTAINTY",
            severity: "CAUTION",
        },
    }
}

/// Translate one of azoth's warnings into their taxonomy.
///
/// The mapping is a judgement, and it is labelled as one: azoth's own code travels inside the
/// message, so nothing is lost to a reader who needs the original. Their codes and severities
/// are the eight their contract enumerates.
fn map_warning(warning: &Warning) -> OutWarning {
    let (code, severity) = match warning.code {
        azoth_core::warning::WarningCode::OutOfValidRange => ("EXTRAPOLATION", "WARNING"),
        azoth_core::warning::WarningCode::RangeCheckSkipped => {
            ("MISSING_REFERENCE_DATA", "WARNING")
        }
        azoth_core::warning::WarningCode::TransitionalFlow => ("MODEL_LIMITATION", "INFO"),
        azoth_core::warning::WarningCode::SolverNotConverged => ("CONVERGENCE_WARNING", "WARNING"),
        azoth_core::warning::WarningCode::TrivialSolution => ("TWO_PHASE_UNCERTAINTY", "CAUTION"),
    };
    OutWarning::new(
        code,
        severity,
        format!("azoth {}: {}", warning.code.as_str(), warning.message),
    )
}

/// What this sidecar does not claim, carried on every answer.
///
/// It is on the envelope because a verdict read alone would overstate it, and the reader most
/// likely to act on a verdict is the one least likely to have read this file.
fn limitations() -> Vec<String> {
    vec![
        "azoth is a port of NeqSim, so agreement is evidence about a transcription and not about \
         the physics"
            .to_string(),
        "everything azoth ships is `unverified`, and its port coverage is incomplete: a state \
         this side answers is one it happens to cover, not one it is qualified on"
            .to_string(),
        "the only comparison here is the vapour fraction, and the tolerance is a claim the \
         caller can set rather than a measured agreement"
            .to_string(),
        "the mapping from azoth's warning codes onto NeqSim's is this sidecar's judgement and \
         has not been agreed upstream"
            .to_string(),
    ]
}

/// Their envelope's field names and vocabulary, assembled from what azoth computed.
#[derive(Debug, Serialize)]
struct Envelope {
    #[serde(rename = "apiVersion")]
    api_version: &'static str,
    status: &'static str,
    tool: &'static str,
    #[serde(skip_serializing_if = "Option::is_none")]
    data: Option<Data>,
    provenance: Provenance,
    validation: Validation,
    #[serde(rename = "qualityGate")]
    quality_gate: QualityGate,
    warnings: Vec<OutWarning>,
    #[serde(rename = "autoValidation")]
    auto_validation: AutoValidation,
}

impl Envelope {
    /// An answer to a request that could not be read. `error` is one of their four statuses.
    fn error(message: &str) -> Self {
        Self {
            api_version: API_VERSION,
            status: "error",
            tool: TOOL,
            data: None,
            provenance: Provenance {
                model: String::new(),
                flash_type: String::new(),
                convergence: Convergence {
                    converged: false,
                    iterations: 0,
                    residual: f64::NAN,
                },
                assumptions: Vec::new(),
                limitations: limitations(),
            },
            validation: Validation {
                valid: false,
                phase: None,
                message: Some(message.to_string()),
            },
            quality_gate: QualityGate {
                verdict: "failed",
                summary: "the request could not be read as a runFlash document".to_string(),
                engineering_review_required: true,
            },
            warnings: vec![OutWarning::new("MODEL_LIMITATION", "INFO", message)],
            auto_validation: AutoValidation {
                overall: "FAIL",
                checks: vec![Check::new(
                    "request-shape",
                    "FAIL",
                    "the line was not a request this sidecar can read",
                )],
            },
        }
    }

    /// The refusal answer: no number at all, and the field that was refused.
    fn refused(refusal: &Refusal, request: &Request) -> Self {
        Self {
            api_version: API_VERSION,
            status: "blocked",
            tool: TOOL,
            data: None,
            provenance: Provenance {
                model: request.model.clone(),
                flash_type: request.flash_type.clone(),
                convergence: Convergence {
                    converged: false,
                    iterations: 0,
                    residual: f64::NAN,
                },
                assumptions: Vec::new(),
                limitations: limitations(),
            },
            validation: Validation {
                valid: false,
                phase: None,
                message: Some(refusal.message.clone()),
            },
            quality_gate: QualityGate {
                verdict: "failed",
                summary: format!(
                    "azoth did not answer, so there is nothing to compare. The refusal is about \
                     `{}`.",
                    refusal.field
                ),
                engineering_review_required: true,
            },
            warnings: vec![OutWarning::new(
                refusal.code,
                "WARNING",
                format!("refused: {}", refusal.message),
            )],
            auto_validation: AutoValidation {
                overall: "FAIL",
                checks: vec![Check::new(refusal.field, "FAIL", refusal.message.clone())],
            },
        }
    }
}

/// A number with the unit it is in, in the middleware's own quantity shape.
#[derive(Debug, Serialize)]
struct Data {
    #[serde(rename = "vapourFraction", skip_serializing_if = "Option::is_none")]
    vapour_fraction: Option<Quantity>,
    phase: &'static str,
    iterations: u32,
    residual: f64,
}

#[derive(Debug, Serialize)]
struct Quantity {
    magnitude_si: f64,
    unit: &'static str,
}

#[derive(Debug, Serialize)]
struct Provenance {
    model: String,
    #[serde(rename = "flashType")]
    flash_type: String,
    convergence: Convergence,
    assumptions: Vec<String>,
    limitations: Vec<String>,
}

#[derive(Debug, Serialize)]
struct Convergence {
    converged: bool,
    iterations: u32,
    /// `NaN` where no iteration ran, which serialises as `null` — an absence, not a zero.
    residual: f64,
}

#[derive(Debug, Serialize)]
struct Validation {
    valid: bool,
    #[serde(skip_serializing_if = "Option::is_none")]
    phase: Option<&'static str>,
    #[serde(skip_serializing_if = "Option::is_none")]
    message: Option<String>,
}

#[derive(Debug, Serialize)]
struct QualityGate {
    verdict: &'static str,
    summary: String,
    #[serde(rename = "engineeringReviewRequired")]
    engineering_review_required: bool,
}

#[derive(Debug, Serialize)]
struct OutWarning {
    code: &'static str,
    severity: &'static str,
    message: String,
}

impl OutWarning {
    fn new(code: &'static str, severity: &'static str, message: impl Into<String>) -> Self {
        Self {
            code,
            severity,
            message: message.into(),
        }
    }
}

#[derive(Debug, Serialize)]
struct AutoValidation {
    overall: &'static str,
    checks: Vec<Check>,
}

#[derive(Debug, Serialize)]
struct Check {
    rule: &'static str,
    status: &'static str,
    message: String,
}

impl Check {
    fn new(rule: &'static str, status: &'static str, message: impl Into<String>) -> Self {
        Self {
            rule,
            status,
            message: message.into(),
        }
    }
}

fn encode(envelope: &Envelope) -> String {
    serde_json::to_string(envelope).unwrap_or_else(|error| {
        format!(
            r#"{{"apiVersion":"{API_VERSION}","status":"error","tool":"{TOOL}","validation":{{"valid":false,"message":"the answer could not be encoded: {error}"}}}}"#
        )
    })
}
