//! Running a flowsheet: the outer loop, the tear, and the streams by name.
//!
//! **This is `ProcessSystem.runSequential`'s loop and not a better one.** The shape is the
//! class's:
//!
//! ```text
//! do {
//!   iter++;
//!   for each unit in execution order: run it
//!   for each Recycle: run it, take its four residuals
//!   isConverged = every Recycle solved
//! } while ((!isConverged || (iter < 2 && hasRecycle && (requireRecycleConfirmation()
//!                                        || hasAutoDeactivatedRecycle()))) && iter < 100)
//! ```
//!
//! with three things from it worth naming. **`iter < 100` is a hard cap and there is no
//! failure.** A flowsheet that has not converged after a hundred passes returns what it has and
//! says so - the class's own loop simply falls out and publishes a `WARNING` severity event, and
//! the caller reads `isConverged`. **`iter < 2 && hasRecycle`** is the class's insistence that a
//! flowsheet with a recycle run at least twice, which `Recycle.solved()`'s own `iterations > 1`
//! restates from the other side; both disjuncts of its guard are true on an ordinary run, because
//! `RecycleController.init()` resets every recycle to zero iterations before the loop starts.
//! And **a unit that is not active is skipped rather than run** (`runUnitProfiled`), which is what
//! makes a tear the low-flow cutoff switched off keep its own counter while the loop around it
//! runs on.
//!
//! **A stream is named by the endpoint that produces it.** A connection is `from -> to`, and
//! `from` is a feed name or `instance.port`; the map is keyed by that, which is what makes a
//! name addressable by both a widget and a report.
//!
//! **The tear is read one pass behind and written at the end of the pass.** That is what a tear
//! *is*: `mix1` consumes `recycle_1` before `sep1` has produced this pass's `sep1.liquid`, and
//! the value it gets is the previous pass's. `Recycle`'s `lastIterationStream` is the same
//! device - it starts as the current stream, which is why its convergence test needs
//! `iterations > 1`.

use std::collections::BTreeMap;

use azoth_core::{AzothError, Result};

use crate::channel::{Direction, Multiplicity};
use crate::executor::dispatch::{Parameters, dispatch};
use crate::flowsheet::{Flowsheet, stream_path};
use crate::order::{ExecutionOrder, execution_order};
use crate::recycle::{
    Acceleration, BroydenAccelerator, RecycleSettings, Residuals, WEGSTEIN_DELAY_ITERATIONS,
    WEGSTEIN_Q_MAX, WEGSTEIN_Q_MIN, apply_composition, extract, mass_flow_kg_per_hr, residuals,
    solved, wegstein,
};
use crate::stream::Stream;
use crate::unit_op::UnitOpSpec;

/// The class's `iter < 100`, the outer loop's cap.
///
/// **There is no failure at the cap.** The class falls out of the loop and reports what it
/// reached; so does this, with `converged` false.
pub const MAX_PASSES: u32 = 100;

/// What one tear did across the run.
#[derive(Debug, Clone, PartialEq)]
pub struct TearRecord {
    /// The recycle's own name, the fresh name the tear is bound to.
    pub stream: String,
    /// Passes the tear was evaluated on.
    pub iterations: u32,
    /// Its four residuals at the last pass, `None` when no pass could compute them.
    pub residuals: Option<Residuals>,
    /// Whether the last pass's `solved()` held.
    pub solved: bool,
    /// Whether the tear was evaluated at all.
    ///
    /// **False means `deactivateOnLowFlow` switched it off**, its residuals are declared zero
    /// rather than measured, and `solved()` answered true for it — so a tear reading `solved` and
    /// `!active` converged by being *absent* rather than by closing. A front-end has to tell those
    /// apart, which is why this is not folded into `solved`.
    pub active: bool,
}

/// What a run reached.
///
/// **Not `PartialEq`, because a `Stream` is not** - a stream carries a computed enthalpy, and the
/// two implementations' are compared by a case's tolerance rather than by equality.
#[derive(Debug, Clone)]
pub struct RunReport {
    /// Passes taken.
    pub iterations: u32,
    /// Whether every tear solved within the cap.
    pub converged: bool,
    /// One record per declared recycle, in declaration order.
    pub tears: Vec<TearRecord>,
    /// Every stream the run produced, keyed by the endpoint that produced it.
    pub streams: BTreeMap<String, Stream>,
    /// Every unit op's own result, keyed by instance id, where it computed one.
    ///
    /// **The numbers a unit operation reached that are on no outlet stream** - a duty, a tray
    /// profile, a conversion, a convergence - published by the kernel that computed them rather
    /// than recomputed by anything downstream. A unit op whose whole answer is its outlets has no
    /// entry, which is a statement about the arithmetic and not a gap: a mixer's result *is* its
    /// mixed stream.
    pub results: BTreeMap<String, serde_json::Value>,
}

/// The acceleration's own state, carried across a run's passes.
///
/// **Wegstein stores a pair and Broyden stores a matrix**, and both are per tear: two tears in one
/// flowsheet accelerate independently, which is what the class gives each `Recycle` object of its
/// own. Which of the three methods runs is the declaration's, and `DirectSubstitution` uses
/// neither field.
#[derive(Debug, Clone)]
struct AccelerationState {
    broyden: BroydenAccelerator,
    previous_input: Option<Vec<f64>>,
    previous_output: Option<Vec<f64>>,
}

/// `Recycle.run`'s acceleration branch: the class's own methods, order and delays.
///
/// **It runs after the flash and before the residuals, and it writes onto the stream being
/// published** - so the residuals are measured against the accelerated state and not the flashed
/// one. That is the class's order, and it is why the step is visible at all.
///
/// **The two delays are in different places, and both are the class's.** Wegstein's is its own
/// field (`wegsteinDelayIterations = 2`), so the branch is taken from the third pass - and
/// `applyWegsteinAcceleration`'s `previousInputValues == null` then makes that third pass the
/// *identity* that seeds the secant, which puts the first real step on the fourth. Broyden has no
/// delay here at all; its delay is inside the accelerator, where the same two calls substitute
/// directly. Same first step, reached by a different route.
fn accelerate(
    state: &mut AccelerationState,
    settings: &RecycleSettings,
    iterations: u32,
    previous: &Stream,
    current: &Stream,
) -> Stream {
    match settings.acceleration {
        Acceleration::DirectSubstitution => current.clone(),
        Acceleration::Wegstein if iterations > WEGSTEIN_DELAY_ITERATIONS => {
            let input = extract(previous);
            let output = extract(current);
            let step = wegstein(
                &input,
                &output,
                state.previous_input.as_deref(),
                state.previous_output.as_deref(),
                WEGSTEIN_Q_MIN,
                WEGSTEIN_Q_MAX,
            );
            state.previous_input = Some(input);
            state.previous_output = Some(output);
            apply_composition(current, &step.values)
        }
        // A pass inside Wegstein's delay.
        Acceleration::Wegstein => current.clone(),
        Acceleration::Broyden => {
            let accelerated = state
                .broyden
                .accelerate(&extract(previous), &extract(current));
            apply_composition(current, &accelerated)
        }
    }
}

/// The boundary inlets: the document's own record, with the caller's overriding by name.
///
/// **A document is self-contained, and this is what that means mechanically.** Its `[[inputs]]`
/// build their own streams, so `run` takes an empty map for the ordinary case; a caller that
/// supplies one is stating a different value for a boundary the document already declares, which
/// is how a sweep or a front-end fills a form without editing the file.
///
/// A supplied name the document does not declare is **refused rather than added**: a boundary no
/// input declares is a stream nothing consumes, and the run-time twin of `UnknownFeed` is a
/// refusal rather than a stream that reaches no unit.
///
/// # Errors
/// [`AzothError::InvalidInput`] for either of those two, and whatever `input_streams` returns for
/// an input whose fluid does not resolve.
fn boundary(
    flowsheet: &Flowsheet,
    feeds: &BTreeMap<String, Stream>,
) -> Result<BTreeMap<String, Stream>> {
    let mut boundary = flowsheet.input_streams()?;
    for (name, stream) in feeds {
        if !boundary.contains_key(name) {
            return Err(AzothError::invalid_input(
                "inputs",
                format!(
                    "`{name}` is not an input this flowsheet declares, so nothing consumes it - \
                     a boundary the document does not state is a wiring error and not an extra \
                     feed"
                ),
            ));
        }
        boundary.insert(name.clone(), stream.clone());
    }
    Ok(boundary)
}

/// Run a flowsheet to its steady state.
///
/// `feeds` overrides the boundary inlets the document's `[[inputs]]` declare, keyed by their
/// names; an empty map is the ordinary call and runs the document as written.
///
/// # Errors
/// * [`AzothError::InvalidInput`] if an input's fluid does not resolve, if a supplied feed is not
///   one the document declares, if a connection names an instance or a port the palette does not
///   declare, if a single-multiplicity inlet receives more than one stream, or if a unit refuses
///   its arguments.
pub fn run(
    flowsheet: &Flowsheet,
    palette: &[UnitOpSpec],
    feeds: &BTreeMap<String, Stream>,
    order: ExecutionOrder,
) -> Result<RunReport> {
    let sequence = execution_order(flowsheet, order)?;
    let boundary = boundary(flowsheet, feeds)?;
    // **A declaration the port does not carry stops the run rather than being ignored.** That is
    // the whole of this branch's history: `acceleration_method` was accepted and did nothing, and
    // the refusal below is what makes the accepted names mean something.
    for recycle in &flowsheet.recycles {
        if let Some(error) = recycle.unsupported_acceleration() {
            return Err(error);
        }
    }
    let settings: Vec<RecycleSettings> = flowsheet
        .recycles
        .iter()
        .map(|recycle| recycle.settings())
        .collect::<Result<_>>()?;
    let mut tears: BTreeMap<String, Stream> = BTreeMap::new();
    let mut records: Vec<TearRecord> = flowsheet
        .recycles
        .iter()
        .map(|recycle| TearRecord {
            stream: recycle.stream.clone(),
            iterations: 0,
            residuals: None,
            solved: false,
            active: true,
        })
        .collect();

    // One state per declared tear, in declaration order - the class gives each `Recycle` object
    // its own, so two loops in one flowsheet do not share an accelerator.
    let mut acceleration: Vec<AccelerationState> = flowsheet
        .recycles
        .iter()
        .map(|_| AccelerationState {
            broyden: BroydenAccelerator::new(),
            previous_input: None,
            previous_output: None,
        })
        .collect();

    let mut iterations = 0u32;
    let mut converged = false;
    let mut streams: BTreeMap<String, Stream> = BTreeMap::new();
    // **Overwritten every pass, exactly as `streams` is.** A pass re-runs every unit op, so the
    // result this key held is the *previous* pass's - and a tear's converged answer is the last
    // pass's, not the second-to-last one's.
    let mut results: BTreeMap<String, serde_json::Value> = BTreeMap::new();

    for pass in 1..=MAX_PASSES {
        iterations = pass;
        streams = boundary.clone();
        results = BTreeMap::new();

        for &index in &sequence {
            let instance = &flowsheet.instances[index];
            let spec = palette
                .iter()
                .find(|spec| spec.id == instance.unit)
                .ok_or_else(|| {
                    AzothError::invalid_input(
                        "unit",
                        format!(
                            "`{}` is not a unit operation the palette declares",
                            instance.unit
                        ),
                    )
                })?;
            let inlets = inlets_of(flowsheet, spec, instance, &streams, &tears)?;
            let parameters = Parameters::new(spec, &instance.parameters);
            let outcome = dispatch(&instance.unit, &inlets, &parameters)?;
            if let Some(result) = outcome.result {
                results.insert(instance.id.clone(), result);
            }
            bind_products(flowsheet, spec, instance, outcome.streams, &mut streams)?;
        }

        // The tears, once every unit has run - which is what makes the value they publish this
        // pass the one the *next* pass reads.
        let mut all_solved = true;
        for (slot, recycle) in flowsheet.recycles.iter().enumerate() {
            let record = &mut records[slot];

            // **A tear the low-flow cutoff switched off is not evaluated again, so its own counter
            // stops where the cutoff left it.** The class's `runUnitProfiled` returns early for
            // equipment that is not active, so `Recycle.run` is never entered a second time - and
            // the extra pass the outer loop still insists on does not advance `getIterations()`.
            // That is why the capture beside this reads `recycle_iterations=1` for a loop the
            // class ran twice. The tear's published value is left as it was, which is what the
            // skipped unit's outlet stream is.
            if !record.active {
                all_solved &= record.solved;
                continue;
            }

            let current = streams.get(&recycle.from).ok_or_else(|| {
                AzothError::invalid_input(
                    "recycles",
                    format!(
                        "`{}` is torn from `{}`, which the run produced no stream for",
                        recycle.stream, recycle.from
                    ),
                )
            })?;
            let previous = tears.get(&recycle.stream);
            let current_kg = mass_flow_kg_per_hr(current)?;
            record.iterations = pass;

            // **`Recycle.run`'s low-flow cutoff, and it comes before any comparison.** A tear
            // whose inlet carries less than `minimumFlow` kg/hr is deactivated: its residuals are
            // *declared* zero rather than measured, it is marked inactive, and `solved()` answers
            // true for it - so a loop carrying nothing stops here rather than running until the
            // zero-flow floor closes it a pass later. The write below is the class's
            // `lastIterationStream = mixedStream.clone()`, which is what stops a deactivated tear
            // comparing itself against a stale pre-collapse snapshot forever.
            if current_kg < settings[slot].minimum_flow_kg_per_hr {
                let deactivated = Residuals::deactivated();
                let is_solved = solved(
                    &deactivated,
                    &settings[slot],
                    current_kg,
                    current_kg,
                    pass,
                    false,
                );
                record.residuals = Some(deactivated);
                record.active = false;
                record.solved = is_solved;
                all_solved &= is_solved;
                tears.insert(recycle.stream.clone(), current.clone());
                continue;
            }

            // What this pass publishes: the flashed stream, or the accelerated one where the
            // declaration asks for it. `None` means the first pass, which publishes the flash.
            let mut published = None;
            match previous {
                None => {
                    // The first pass has nothing to compare against, which is the class's own
                    // `lastIterationStream` before it is set.
                    record.solved = false;
                }
                Some(previous) => {
                    // **The class accelerates the stream it is about to publish**, between the
                    // flash and the residuals - so both the answer and the convergence test see
                    // the accelerated state rather than the flashed one.
                    let accelerated = accelerate(
                        &mut acceleration[slot],
                        &settings[slot],
                        pass,
                        previous,
                        current,
                    );
                    let measured = residuals(&accelerated, previous)?;
                    let previous_kg = mass_flow_kg_per_hr(previous)?;
                    let is_solved = solved(
                        &measured,
                        &settings[slot],
                        current_kg,
                        previous_kg,
                        pass,
                        true,
                    );
                    record.residuals = Some(measured);
                    record.solved = is_solved;
                    all_solved &= is_solved;
                    published = Some(accelerated);
                }
            }
            tears.insert(
                recycle.stream.clone(),
                published.unwrap_or_else(|| current.clone()),
            );
        }

        // **The class's own exit condition.** Every tear solved, and a flowsheet with a recycle
        // has run at least twice - which the tear's own `iterations > 1` enforces from the
        // residual side and which is restated here because it is the class's *other* clause:
        // `iter < 2 && hasRecycle && (requireRecycleConfirmation() || hasAutoDeactivatedRecycle())`.
        // The first of those two disjuncts is true whenever a recycle starts at zero iterations,
        // which `RecycleController.init()` has just reset it to, and the second is true exactly
        // when a tear was switched off above - so the clause holds in either case and the loop
        // runs a second pass even for a loop that has nothing to converge.
        converged = all_solved && (flowsheet.recycles.is_empty() || pass > 1);
        if converged || pass == MAX_PASSES {
            break;
        }
    }

    Ok(RunReport {
        iterations,
        converged,
        tears: records,
        streams,
        results,
    })
}

/// A unit's inlets, in the order its declaration gives its ports.
///
/// **`multiplicity = "many"` collects every connection that names the port** - a mixer's inlets -
/// and `"one"` refuses a second, because a single port carrying two streams is a flowsheet that
/// does not say which is which.
fn inlets_of(
    flowsheet: &Flowsheet,
    spec: &UnitOpSpec,
    instance: &crate::flowsheet::Instance,
    streams: &BTreeMap<String, Stream>,
    tears: &BTreeMap<String, Stream>,
) -> Result<Vec<Stream>> {
    let mut inlets = Vec::new();
    for port in spec
        .ports
        .iter()
        .filter(|port| port.direction == Direction::In)
    {
        let endpoint = format!("{}.{}", instance.id, port.name);
        let mut gathered = Vec::new();
        for connection in &flowsheet.connections {
            if connection.to != endpoint {
                continue;
            }
            gathered.push(streams.get(&connection.from).cloned().ok_or_else(|| {
                AzothError::invalid_input(
                    "connections",
                    format!(
                        "`{endpoint}` is fed by `{}`, which no unit has produced - a connection \
                         from a stream nothing writes is a wiring error",
                        connection.from
                    ),
                )
            })?);
        }
        // **A `[[recycles]]` entry is an edge, and it is the *only* edge a tear has.** The
        // shipped flowsheet declares `sep1.liquid -> mix1.feed` as a recycle and writes no
        // `[[connections]]` line for it, so a walk that read only the connections would never
        // feed the mixer from the loop at all - the tear would be switched off by the low-flow
        // cutoff and converge by carrying nothing. Its stream is the tear, which the *first* pass
        // has not produced,
        // so it contributes nothing then and the loop closes on the second.
        for recycle in &flowsheet.recycles {
            if recycle.to != endpoint {
                continue;
            }
            if let Some(stream) = tears.get(&recycle.stream) {
                gathered.push(stream.clone());
            }
        }
        if port.multiplicity == Multiplicity::One && gathered.len() > 1 {
            return Err(AzothError::invalid_input(
                "ports",
                format!(
                    "`{endpoint}` declares a single inlet and {} connections name it",
                    gathered.len()
                ),
            ));
        }
        inlets.extend(gathered);
    }
    if inlets.is_empty() {
        return Err(AzothError::invalid_input(
            "connections",
            format!(
                "`{}` has no inlet, so there is nothing to run it on",
                instance.id
            ),
        ));
    }
    Ok(inlets)
}

/// The streams a unit's outlets become, keyed by the endpoint that produces them.
fn bind_products(
    flowsheet: &Flowsheet,
    spec: &UnitOpSpec,
    instance: &crate::flowsheet::Instance,
    outlets: Vec<Stream>,
    streams: &mut BTreeMap<String, Stream>,
) -> Result<()> {
    let ports: Vec<&crate::channel::Port> = spec
        .ports
        .iter()
        .filter(|port| port.direction == Direction::Out)
        .collect();

    // **The ports are walked in declaration order, and a `many` one takes the remainder.** A
    // kernel returns the streams flat, in that same order, so the question is only where one
    // port's streams end - and the answer is the ports after it, each of which takes exactly one.
    // Refusing a count that does not fit is what the count check always did; what is new is that
    // a `many` port's count is not one, so the check is now against what the walk consumes.
    let mut distribution: Vec<(String, Option<usize>, Stream)> = Vec::with_capacity(outlets.len());
    let mut cursor = 0usize;
    for (position, port) in ports.iter().enumerate() {
        let after = ports.len() - position - 1;
        let take = if port.multiplicity == Multiplicity::Many {
            // Whatever is left once every later port has taken its one.
            outlets.len().saturating_sub(cursor).saturating_sub(after)
        } else {
            1
        };
        for index in 0..take {
            let Some(stream) = outlets.get(cursor).cloned() else {
                return Err(AzothError::invalid_input(
                    "ports",
                    format!(
                        "`{}` returned {} outlets for {} declared outlet ports, the last of which \
                         takes {}",
                        instance.id,
                        outlets.len(),
                        ports.len(),
                        take
                    ),
                ));
            };
            cursor += 1;
            let position = (port.multiplicity == Multiplicity::Many).then_some(index);
            distribution.push((port.name.clone(), position, stream));
        }
    }
    if cursor != outlets.len() {
        return Err(AzothError::invalid_input(
            "ports",
            format!(
                "`{}` returned {} outlets for {} declared outlet ports, which take {}",
                instance.id,
                outlets.len(),
                ports.len(),
                cursor
            ),
        ));
    }

    for (port, index, stream) in distribution {
        streams.insert(stream_path(&instance.id, &port, index), stream);
    }
    // **A boundary product is a stream too.** A connection whose `to` is a bare name is what the
    // environment consumes, and the run reports it under that name.
    for connection in &flowsheet.connections {
        if connection.to.contains('.') {
            continue;
        }
        if let Some(stream) = streams.get(&connection.from).cloned() {
            streams.insert(connection.to.clone(), stream);
        }
    }
    Ok(())
}

/// A live flowsheet and its named results.
///
/// **This is the object a widget and an agent both address**, and what makes that possible is that
/// every value has a **stable path**. The paths are the endpoints the run produced plus the
/// fields the palette declares, so `<endpoint>.<field>` — `p1.outlet.P`, `feed_1.n` — is derivable
/// from the declaration rather than invented here. `middleware.md`'s own example spells one
/// `p1.outlet_pressure`; the field names are the declaration's (`n`, `z`, `P`, `T`, `h`), so the
/// path is `p1.outlet.P`, and the page is corrected to match rather than the other way round.
///
/// **A runtime `Stream` has no name.** That is the gap this closes: the stream is addressed
/// *through* the session by the endpoint that produced it, and the session is what holds the two
/// together.
#[derive(Debug, Clone)]
pub struct Session {
    flowsheet: Flowsheet,
    report: RunReport,
}

impl Session {
    /// Run a flowsheet and hold what it reached.
    ///
    /// # Errors
    /// Whatever [`run`] returns.
    pub fn run(
        flowsheet: Flowsheet,
        palette: &[UnitOpSpec],
        feeds: &BTreeMap<String, Stream>,
        order: ExecutionOrder,
    ) -> Result<Self> {
        let report = run(&flowsheet, palette, feeds, order)?;
        Ok(Self { flowsheet, report })
    }

    /// The flowsheet this session is a state of.
    #[must_use]
    pub fn flowsheet(&self) -> &Flowsheet {
        &self.flowsheet
    }

    /// What the run reached, including the tear records.
    #[must_use]
    pub fn report(&self) -> &RunReport {
        &self.report
    }

    /// Whether every tear solved within the cap.
    #[must_use]
    pub fn converged(&self) -> bool {
        self.report.converged
    }

    /// Passes the run took.
    #[must_use]
    pub fn iterations(&self) -> u32 {
        self.report.iterations
    }

    /// One stream, by the endpoint that produced it.
    ///
    /// # Errors
    /// [`AzothError::InvalidInput`] if no stream has that name, listing what does - a path a
    /// caller guessed is worth a list rather than a `None`.
    pub fn stream(&self, endpoint: &str) -> Result<&Stream> {
        self.report.streams.get(endpoint).ok_or_else(|| {
            let mut known: Vec<&str> = self.report.streams.keys().map(String::as_str).collect();
            known.sort_unstable();
            AzothError::invalid_input(
                "path",
                format!("no stream is named `{endpoint}`; this run produced {known:?}"),
            )
        })
    }

    /// Every addressable path, sorted.
    ///
    /// **Both forms are listed**: `<endpoint>` for the stream itself and `<endpoint>.<field>` for
    /// each scalar the declaration gives it, so a front-end can enumerate what it may point at
    /// without guessing.
    #[must_use]
    pub fn paths(&self) -> Vec<String> {
        let mut paths = Vec::new();
        for endpoint in self.report.streams.keys() {
            paths.push(endpoint.clone());
            for field in STREAM_FIELDS {
                paths.push(format!("{endpoint}.{field}"));
            }
        }
        for tear in &self.report.tears {
            paths.push(tear.stream.clone());
        }
        paths.sort();
        paths
    }

    /// One named scalar, by path.
    ///
    /// # Errors
    /// [`AzothError::InvalidInput`] if the path has no field, if no stream has its endpoint, or
    /// if the field is `z` - a composition is a vector and not a scalar, and returning its first
    /// entry would be the kind of guess this session exists to avoid.
    pub fn value(&self, path: &str) -> Result<f64> {
        let (endpoint, field) = path.rsplit_once('.').ok_or_else(|| {
            AzothError::invalid_input(
                "path",
                format!("`{path}` names no field; a scalar is `<endpoint>.<field>` with field one of {STREAM_FIELDS:?}"),
            )
        })?;
        if field == "z" {
            return Err(AzothError::invalid_input(
                "path",
                format!(
                    "`{path}` is a composition: `z` is a vector and has no single value - read \
                     the stream with `stream({endpoint:?})`"
                ),
            ));
        }
        let stream = self.stream(endpoint)?;
        match field {
            "n" => Ok(stream.n),
            "P" => Ok(stream.p.value),
            "T" => Ok(stream.t.value),
            "h" => Ok(stream.h.value),
            other => Err(AzothError::invalid_input(
                "path",
                format!("`{other}` is not a field of a stream record; they are {STREAM_FIELDS:?}"),
            )),
        }
    }

    /// One tear's record.
    ///
    /// # Errors
    /// [`AzothError::InvalidInput`] if no recycle carries that name.
    pub fn tear(&self, name: &str) -> Result<&TearRecord> {
        self.report
            .tears
            .iter()
            .find(|tear| tear.stream == name)
            .ok_or_else(|| {
                AzothError::invalid_input(
                    "recycles",
                    format!("no recycle is named `{name}` in `{}`", self.flowsheet.id),
                )
            })
    }
}

/// The fields a port's record declares, and therefore the ones a path may name.
///
/// **The declaration's own names, not a spelling of this module's**: a palette port writes `n`,
/// `z`, `P`, `T`, `h`, so a path is `p1.outlet.P` and not `p1.outlet.pressure`.
pub const STREAM_FIELDS: [&str; 5] = ["n", "z", "P", "T", "h"];

#[cfg(test)]
mod tests {
    use super::*;
    use azoth_core::units::{kelvins, pascals};

    fn stream(methane: f64) -> Stream {
        Stream::from_pt(
            vec!["methane".to_string(), "n-butane".to_string()],
            vec![methane, 1.0 - methane],
            1.0,
            pascals(5.0e5),
            kelvins(300.0),
        )
        .expect("the two resolve")
    }

    fn settings(acceleration: Acceleration) -> RecycleSettings {
        RecycleSettings {
            acceleration,
            ..RecycleSettings::default()
        }
    }

    /// **The delays are the class's, and they are in two different places.**
    ///
    /// Wegstein's is a `Recycle` field, so its branch is not taken until the third pass - and the
    /// method's own missing previous pair then makes that third pass the identity that seeds the
    /// secant. Broyden's is inside the accelerator, so its branch is taken from the *second* pass
    /// but the accelerator substitutes directly for its first two calls. **Both put the first real
    /// step on the fourth pass**, which is the number this asserts.
    ///
    /// The loop starts at the second pass because that is where the session's does: a first pass
    /// has no published tear to compare against, so `accelerate` is not reached at all - and a
    /// test that called it there would be measuring a call the run never makes.
    #[test]
    fn both_accelerations_step_for_the_first_time_on_the_fourth_pass() {
        for acceleration in [Acceleration::Wegstein, Acceleration::Broyden] {
            let mut state = AccelerationState {
                broyden: BroydenAccelerator::new(),
                previous_input: None,
                previous_output: None,
            };
            let settings = settings(acceleration);
            let mut published = stream(0.9);
            for pass in 2..=5 {
                let current = stream(0.9 - 0.1 * f64::from(pass));
                let accelerated = accelerate(&mut state, &settings, pass, &published, &current);
                if pass < 4 {
                    assert_eq!(
                        accelerated.z, current.z,
                        "{acceleration:?} stepped on pass {pass}, which its delay forbids"
                    );
                } else {
                    assert_ne!(
                        accelerated.z, current.z,
                        "{acceleration:?} did not step on pass {pass}"
                    );
                }
                published = accelerated;
            }
        }
    }

    /// **The composition is the only field that moves**, which is `applyStreamValues`' own choice
    /// and the reason an acceleration is invisible on a loop that accumulates in its flow.
    #[test]
    fn an_accelerated_step_leaves_everything_but_the_composition_alone() {
        let mut state = AccelerationState {
            broyden: BroydenAccelerator::new(),
            previous_input: None,
            previous_output: None,
        };
        let settings = settings(Acceleration::Broyden);
        let mut published = stream(0.9);
        let mut current = stream(0.9);
        for pass in 2..=4 {
            current = stream(0.9 - 0.1 * f64::from(pass));
            published = accelerate(&mut state, &settings, pass, &published, &current);
        }
        let accelerated = accelerate(&mut state, &settings, 5, &published, &current);

        assert_eq!(accelerated.n, current.n, "the flow is the flash's");
        assert_eq!(
            accelerated.p.value, current.p.value,
            "the pressure is the flash's"
        );
        assert_eq!(
            accelerated.t.value, current.t.value,
            "the temperature is the flash's"
        );
        assert_eq!(
            accelerated.h.value, current.h.value,
            "the enthalpy is the flash's"
        );
        assert_ne!(accelerated.z, current.z, "and the composition is not");
    }

    /// Direct substitution does nothing at all, on any pass - which is the class's default and
    /// what every shipped flowsheet declares by saying nothing.
    #[test]
    fn direct_substitution_leaves_the_stream_exactly_as_the_flash_left_it() {
        let mut state = AccelerationState {
            broyden: BroydenAccelerator::new(),
            previous_input: None,
            previous_output: None,
        };
        let settings = settings(Acceleration::DirectSubstitution);
        let mut published = stream(0.9);
        for pass in 2..=6 {
            let current = stream(0.9 - 0.1 * f64::from(pass));
            let accelerated = accelerate(&mut state, &settings, pass, &published, &current);
            assert_eq!(accelerated.z, current.z);
            published = accelerated;
        }
    }
}
