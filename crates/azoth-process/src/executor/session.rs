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
//! } while ((!isConverged || (iter < 2 && hasRecycle && ...)) && iter < 100)
//! ```
//!
//! with two things from it worth naming. **`iter < 100` is a hard cap and there is no
//! failure.** A flowsheet that has not converged after a hundred passes returns what it has and
//! says so - the class's own loop simply falls out and publishes a `WARNING` severity event, and
//! the caller reads `isConverged`. And **`iter < 2 && hasRecycle`** is the class's insistence
//! that a flowsheet with a recycle run at least twice, which `Recycle.solved()`'s own
//! `iterations > 1` restates from the other side.
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
use crate::flowsheet::Flowsheet;
use crate::order::{ExecutionOrder, execution_order};
use crate::recycle::{RecycleSettings, Residuals, mass_flow_kg_per_hr, residuals, solved};
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
}

/// Run a flowsheet to its steady state.
///
/// `feeds` supplies the boundary inlets, keyed by the names the flowsheet's `feeds` declares.
///
/// # Errors
/// * [`AzothError::InvalidInput`] if a feed is missing, if a connection names an instance or a
///   port the palette does not declare, if a single-multiplicity inlet receives more than one
///   stream, or if a unit refuses its arguments.
pub fn run(
    flowsheet: &Flowsheet,
    palette: &[UnitOpSpec],
    feeds: &BTreeMap<String, Stream>,
    order: ExecutionOrder,
) -> Result<RunReport> {
    let sequence = execution_order(flowsheet, order)?;
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
        })
        .collect();

    let mut iterations = 0u32;
    let mut converged = false;
    let mut streams: BTreeMap<String, Stream> = BTreeMap::new();

    for pass in 1..=MAX_PASSES {
        iterations = pass;
        streams = feeds.clone();

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
            let outlets = dispatch(&instance.unit, &inlets, &parameters)?;
            bind_products(flowsheet, spec, instance, outlets, &mut streams)?;
        }

        // The tears, once every unit has run - which is what makes the value they publish this
        // pass the one the *next* pass reads.
        let mut all_solved = true;
        for (slot, recycle) in flowsheet.recycles.iter().enumerate() {
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
            let record = &mut records[slot];
            record.iterations = pass;
            match previous {
                None => {
                    // The first pass has nothing to compare against, which is the class's own
                    // `lastIterationStream` before it is set.
                    record.solved = false;
                }
                Some(previous) => {
                    let measured = residuals(current, previous)?;
                    let current_kg = mass_flow_kg_per_hr(current)?;
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
                }
            }
            tears.insert(recycle.stream.clone(), current.clone());
        }

        // **The class's own exit condition.** Every tear solved, and a flowsheet with a recycle
        // has run at least twice - which the tear's own `iterations > 1` already enforces, and
        // which is restated here because a second pass is what the class's clause is for.
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
        // feed the mixer from the loop at all - the tear would close on the zero-flow floor and
        // never on a comparison. Its stream is the tear, which the *first* pass has not produced,
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
    let names: Vec<&str> = spec
        .ports
        .iter()
        .filter(|port| port.direction == Direction::Out)
        .map(|port| port.name.as_str())
        .collect();
    if outlets.len() != names.len() {
        // **A splitter is the case, and it is a gap rather than a mismatch.** `unit_ops.splitter`
        // declares its outlets as one `many` port named `products`, so a two-way split returns
        // two streams for one port name and there is no way to say which is which: a connection
        // writes `split1.products`, and both outlets would answer to it. Refusing here is the
        // honest disposition - the alternative is to bind the last one and lose the first - and
        // the class that would close it is the port declaration, which needs a way to name an
        // outlet by index.
        let many = spec.ports.iter().any(|port| {
            port.direction == Direction::Out && port.multiplicity == Multiplicity::Many
        });
        let detail = if many && outlets.len() > names.len() {
            " - the entry declares a `many` outlet port, and a connection cannot name one stream \
             of it"
        } else {
            ""
        };
        return Err(AzothError::invalid_input(
            "ports",
            format!(
                "`{}` returned {} outlets for {} declared outlet ports{detail}",
                instance.id,
                outlets.len(),
                names.len()
            ),
        ));
    }
    for (name, stream) in names.into_iter().zip(outlets) {
        streams.insert(format!("{}.{}", instance.id, name), stream);
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
