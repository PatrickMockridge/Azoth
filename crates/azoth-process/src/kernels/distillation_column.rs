//! `unit_ops.distillation_column` - the tray-by-tray sequential-substitution solve.

use azoth_core::units::{Power, Pressure, ThermodynamicTemperature, kelvins, pascals, watts};
use azoth_core::{AzothError, Result};

use crate::column::tray::{self as stage, TrayOutcome};
use crate::stream::Stream;

/// The class's own adaptive-relaxation constants, from `DistillationColumn`'s field
/// initialisers.
const MIN_SEQUENTIAL_RELAXATION: f64 = 0.5;
const MAX_ADAPTIVE_RELAXATION: f64 = 1.2;
const RELAXATION_INCREASE_FACTOR: f64 = 1.2;
const RELAXATION_DECREASE_FACTOR: f64 = 0.5;
/// The floor on the *temperature* update's step, which is a different clamp from the streams'.
const MIN_TEMPERATURE_RELAXATION: f64 = 0.2;
/// `DEFAULT_MASS_BALANCE_TOLERANCE`, the class's own.
const MASS_BALANCE_TOLERANCE: f64 = 1.6e-2;
/// `DEFAULT_ENTHALPY_BALANCE_TOLERANCE`, the class's own.
const ENTHALPY_BALANCE_TOLERANCE: f64 = 1.6e-2;

/// One tray's solved state, as the capture prints it.
#[derive(Debug, Clone)]
pub struct TrayProfile {
    /// The tray's temperature.
    pub temperature: ThermodynamicTemperature,
    /// The tray's pressure.
    pub pressure: Pressure,
    /// The vapour leaving the tray, mol/s - zero where the flash found none.
    pub gas_n: f64,
    /// The liquid leaving the tray, mol/s - zero where the flash found none.
    pub liquid_n: f64,
    /// The vapour's composition.
    pub gas_z: Vec<f64>,
    /// The liquid's composition.
    pub liquid_z: Vec<f64>,
}

/// What a column solve hands back.
///
/// **Three residuals and no `converged` flag.** Where a solve stops is a measurement, and a
/// flag a caller can read is a number nothing gated. A solve that misses the gate is a
/// refusal naming the residuals - the precedent `e07` set - and this column needs it: on a PR
/// fluid the deethanizer's captured state is a *partial* solve after 80 hard-capped
/// iterations, and the capture keeps it as evidence rather than as an answer.
#[derive(Debug, Clone)]
pub struct ColumnOutcome {
    /// The trays, from the reboiler up to the condenser.
    pub trays: Vec<TrayProfile>,
    /// The overhead product.
    pub distillate: Stream,
    /// The bottom product.
    pub bottoms: Stream,
    /// The condenser's duty, W: its products' enthalpy less the vapour's.
    pub condenser_duty: Power,
    /// The reboiler's duty, W: its boilup and bottoms less the downcomer's.
    pub reboiler_duty: Power,
    /// Iterations taken.
    pub iterations: u32,
    /// The mean tray-temperature change at the last iteration, K.
    pub temperature_residual: f64,
    /// The products' worst component imbalance against the feed, relative.
    pub mass_residual: f64,
    /// `|H_feed + duties - H_products| / |H_feed|`.
    pub energy_residual: f64,
    /// **What the stages had to fall back on, once per kind.** A reactive tray falls back the way
    /// the class does - a reactive PH flash to a reactive TP one to a plain flash - and each kind
    /// is reported here rather than taken silently. Deduplicated by code, because a column runs
    /// its trays hundreds of times and a caller needs the kind rather than the count.
    pub warnings: Vec<azoth_core::warning::Warning>,
}

/// Which of the class's ten solving strategies the column runs.
///
/// **Two of the ten are ported and eight are refused by name.** `DIRECT_SUBSTITUTION` is the
/// class's own default and the sequential-substitution core; `NAPHTALI_SANDHOLM` is the
/// simultaneous MESH correction, which is a *different* answer only where the substitution
/// core's own path matters - and it is the only strategy under which NeqSim converges the
/// deethanizer at all, which is why it is the second solve rather than a later one. The
/// remaining eight are `ColumnSolverFactory`'s:
/// `DAMPED_SUBSTITUTION`, `WEGSTEIN`, `SUM_RATES`, `NEWTON`, `INSIDE_OUT`,
/// `MATRIX_INSIDE_OUT`, `MESH_RESIDUAL` and `AUTO` - and `AUTO` is a *ladder* over the
/// others rather than a method, which is what `ColumnSolverFactory`'s `candidateSolvers`
/// makes it: the ladder is stated and refused rather than silently resolved.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SolverType {
    /// Sequential substitution with the class's own adaptive relaxation.
    DirectSubstitution,
    /// Full MESH simultaneous correction, by `NaphtaliSandholmSolver`.
    NaphtaliSandholm,
}

impl SolverType {
    /// The name the model's enum declares, in the class's own spelling.
    #[must_use]
    pub fn name(self) -> &'static str {
        match self {
            Self::DirectSubstitution => "direct_substitution",
            Self::NaphtaliSandholm => "naphtali_sandholm",
        }
    }
}

/// Everything a column solve takes.
///
/// **The stage count excludes the two ends**, which is `DistillationColumn`'s own rule:
/// `numberOfTrays` is the middle trays and a reboiler or a condenser adds one each on top. So
/// `feed_stage` is 0-based over the trays *including* the ends - stage 0 is the reboiler when
/// there is one - which is what `addFeedStream(stream, tray)` is indexed by.
#[derive(Debug, Clone)]
pub struct ColumnSetup {
    /// The feed.
    pub feed: Stream,
    /// The stage the feed enters.
    pub feed_stage: usize,
    /// The theoretical trays between the ends.
    pub number_of_stages: usize,
    /// Whether the column has a reboiler at stage 0.
    pub has_reboiler: bool,
    /// Whether the column has a condenser at the top stage.
    pub has_condenser: bool,
    /// The pressure at the top stage.
    pub top_pressure: Pressure,
    /// The pressure at the bottom stage.
    pub bottom_pressure: Pressure,
    /// The condenser's temperature, which pins the top tray. **Absent means no pin**, and the
    /// end's flash is then at its own enthalpy - which is what `setCondenserTemperature` not
    /// having been called does, and what makes a duty specification reachable at all.
    pub condenser_temperature: Option<ThermodynamicTemperature>,
    /// The reboiler's temperature, which pins the bottom tray. Absent as the condenser's is.
    pub reboiler_temperature: Option<ThermodynamicTemperature>,
    /// The convergence tolerance on the mean tray-temperature change.
    pub temperature_tolerance: f64,
    /// The iteration cap.
    pub max_iterations: usize,
    /// The top product's specification, or `None` where the top is pinned by temperature.
    pub top_specification: Option<Specification>,
    /// The bottom product's specification, or `None`.
    pub bottom_specification: Option<Specification>,
    /// **The second inlet, which enters the top stage.** `AbsorptionColumn.addSolventInStream`
    /// is `addFeedStream(stream, getNumberOfTrays() - 1)` and `StrippingColumn.addRichLiquidStream`
    /// is the same call: the gas and the stripping gas enter stage 0 through `feed`, and the
    /// solvent and the rich liquid enter the top stage through this. The class refuses a second
    /// assignment of either inlet, so two is the whole of the arrangement.
    pub top_feed: Option<Stream>,
    /// **One outlet-temperature pin per tray**, `NaN` where a tray has none - which is what
    /// `SimpleTray.setOutletTemperature` states and what an absorber's own tests use: they pin
    /// every stage to one temperature and solve an isothermal column. A finite entry here wins
    /// over the two ends' own fields, which are the same mechanism reached by the setters the
    /// base class carries.
    pub tray_temperatures: Option<Vec<f64>>,
    /// Which solving strategy to run.
    pub solver_type: SolverType,
    /// **Which trays flash reactively**, over middle-tray indices, as the class states it.
    pub reactive: ReactiveSection,
}

/// Which trays run their flash reactively: `DistillationColumn.setReactive`.
///
/// **The section is stated over *middle*-tray indices, which is the class's own convention**:
/// `setReactive(true, 3, 7)` makes trays 4-8 of the middle section reactive, and a bound the
/// class leaves at minus one means every middle tray. The two ends are **never** reactive -
/// `replaceMiddleTrays` walks only the range between them.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ReactiveSection {
    /// Every tray is an equilibrium stage, which is the class's default.
    None,
    /// Every middle tray.
    All,
    /// A run of middle-tray indices, inclusive, 0-based among the middle trays.
    Section {
        /// The first reactive middle tray.
        start: usize,
        /// The last reactive middle tray.
        end: usize,
    },
}

impl ReactiveSection {
    /// Whether the middle tray at `index` is reactive.
    #[must_use]
    pub fn covers(self, index: usize) -> bool {
        match self {
            Self::None => false,
            Self::All => true,
            Self::Section { start, end } => index >= start && index <= end,
        }
    }
}

/// Which of `ColumnSpecification`'s five types a specification is.
///
/// The two the outer loop does not drive are the two the class applies **directly** to the end
/// itself: a reflux ratio and a duty are settings rather than targets, so
/// `applyDirectSpecification` writes them onto the condenser or the reboiler and no temperature
/// is searched for.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SpecificationKind {
    /// Mole fraction of a named component in the product.
    ProductPurity,
    /// The fraction of a named component's feed that leaves in the product.
    ComponentRecovery,
    /// The product's molar flow, in **mol/hr** - `ColumnSpecification.defaultTargetUnit`'s own
    /// unit for this type, and the one the class's `evaluateSpecError` reads it in.
    ProductFlowRate,
    /// The end's reflux ratio or boilup ratio, applied directly.
    RefluxRatio,
    /// The end's duty, in W, applied directly.
    Duty,
}

impl SpecificationKind {
    /// Whether the outer loop drives it, which is `needsAdjustment`'s own split.
    #[must_use]
    pub fn is_adjusted(self) -> bool {
        !matches!(self, Self::RefluxRatio | Self::Duty)
    }
}

/// One product specification: the type, the target and the component.
///
/// **A specification's location is the slot it is handed to, and not a field.** NeqSim's
/// `ColumnSpecification` carries a `ProductLocation`, but `validateColumnSpecification` refuses
/// a `TOP` specification given to `setBottomSpecification` - "Use setTopSpecification() for TOP
/// specs and setBottomSpecification() for BOTTOM specs" - so the field can only ever agree with
/// the slot it sits in. That is why the model's inputs name the end (`top_specification_*`,
/// `bottom_specification_*`) rather than carrying a location, and why there are two of them
/// rather than a list: the class holds exactly `topSpecification` and `bottomSpecification`.
#[derive(Debug, Clone)]
pub struct Specification {
    /// Which of the five types.
    pub kind: SpecificationKind,
    /// The target value, in the unit its type implies: dimensionless for a purity, a recovery
    /// or a ratio, mol/hr for a flow rate, W for a duty.
    pub target: f64,
    /// The component a purity or a recovery constrains. Unread by the other three types.
    pub component: Option<String>,
}

impl Specification {
    /// The value the product currently has, which is what the error is measured against.
    fn value(&self, product: &Stream, feed: &Stream) -> Result<f64> {
        let index_of = |name: &str| -> Result<usize> {
            product
                .components
                .iter()
                .position(|c| c == name)
                .ok_or_else(|| {
                    AzothError::invalid_input(
                        "component",
                        format!("{name} is not one of this fluid's components"),
                    )
                })
        };
        match self.kind {
            SpecificationKind::ProductPurity => {
                let name = self.component.as_deref().ok_or_else(|| {
                    AzothError::invalid_input(
                        "component",
                        "a purity specification constrains a component, and none was stated",
                    )
                })?;
                Ok(product.z[index_of(name)?])
            }
            SpecificationKind::ComponentRecovery => {
                let name = self.component.as_deref().ok_or_else(|| {
                    AzothError::invalid_input(
                        "component",
                        "a recovery specification constrains a component, and none was stated",
                    )
                })?;
                let index = index_of(name)?;
                let supplied = feed.n * feed.z[index];
                if supplied.abs() <= 1.0e-12 {
                    return Ok(0.0);
                }
                Ok(product.n * product.z[index] / supplied)
            }
            // `getFlowRate("mol/hr")`: the class's own target unit for this type.
            SpecificationKind::ProductFlowRate => Ok(product.n * 3600.0),
            SpecificationKind::RefluxRatio | SpecificationKind::Duty => Ok(0.0),
        }
    }
}

/// How one end of the column is run.
///
/// A **temperature** is the default and the outer loop's knob: the end's tray flashes at it.
/// The other two are what a *direct* specification writes - a reflux ratio replaces the tray
/// with the end's own reflux flash, and a duty replaces the pin with a heat input.
#[derive(Debug, Clone, Copy)]
enum EndMode {
    /// A pinned outlet temperature, which the outer loop moves while it searches.
    Temperature(f64),
    /// The end's reflux or boilup ratio, from a `RefluxRatio` specification.
    RefluxRatio(f64),
    /// A heat input, from a `Duty` specification.
    Duty(f64),
}

/// The network a solve carries between its steps.
struct Network {
    /// Each tray's vapour outlet.
    gas: Vec<Option<Stream>>,
    /// Each tray's liquid outlet.
    liquid: Vec<Option<Stream>>,
    /// The pressure of each tray.
    pressures: Vec<Pressure>,
    /// The tray the feed enters.
    feed_stage: usize,
    /// The tray count, which the ends are part of.
    tray_count: usize,
    /// The feed, which enters `feed_stage`.
    feed: Stream,
    /// The second inlet, which enters the top stage.
    top_feed: Option<Stream>,
    /// How each end is run.
    modes: Vec<EndMode>,
    /// **Which trays flash reactively, one flag per tray and resolved to absolute indices** -
    /// the class states the section over middle trays, and the ends are never in it.
    reactive: Vec<bool>,
    /// What the stages have had to fall back on, one entry per kind.
    warnings: Vec<azoth_core::warning::Warning>,
}

impl Network {
    /// The inlets a tray sees, from its neighbours' current outlets.
    ///
    /// **`seed` re-states each inlet at the tray's own temperature, and only `init` passes
    /// it.** `SimpleTray.init()` sets every inlet stream's temperature to the tray's before
    /// its first run, which is the mechanism by which `DistillationColumn.init`'s linear
    /// temperature profile reaches the trays at all; the sweeps call `run` without it, so
    /// there the inlets keep their own states. The external feed is never re-stated, which is
    /// what the class's own `refreshInternalExternalFeedSystems` restores.
    fn inlets_of(&self, i: usize, seed: Option<&[f64]>) -> Vec<Stream> {
        let mut inlets: Vec<Stream> = Vec::new();
        let at = |stream: &Stream| match seed {
            Some(temperatures) if temperatures[i].is_finite() => {
                at_temperature(stream, temperatures[i])
            }
            _ => stream.clone(),
        };
        if i > 0 {
            if let Some(g) = self.gas[i - 1].as_ref() {
                inlets.push(at(g));
            }
        }
        if i + 1 < self.tray_count {
            if let Some(l) = self.liquid[i + 1].as_ref() {
                inlets.push(at(l));
            }
        }
        inlets.extend(self.feeds_at(i).into_iter().cloned());
        inlets
    }

    /// The feed a tray carries, if any: the main one at its stage, the second at the top.
    ///
    /// **The two are the class's two inlets**, and a tray can carry both only if the feed
    /// stage *is* the top stage - which is a column whose feed and solvent arrive together.
    fn feeds_at(&self, i: usize) -> Vec<&Stream> {
        let mut feeds = Vec::new();
        if i == self.feed_stage {
            feeds.push(&self.feed);
        }
        if i + 1 == self.tray_count {
            if let Some(top) = self.top_feed.as_ref() {
                feeds.push(top);
            }
        }
        feeds
    }

    /// Run one tray from its neighbours' current state.
    fn run(&mut self, i: usize, seed: Option<&[f64]>) -> Result<TrayOutcome> {
        let inlets = self.inlets_of(i, seed);
        self.run_with(i, &inlets)
    }

    /// Run one tray from an explicit inlet list.
    ///
    /// **`init` needs this and the sweeps do not.** During `init` each link is an `addStream`
    /// on a specific neighbour - the reboiler against the *feed stage's* liquid, and not the
    /// first tray's, which is what the sweeps use - so the inlet list is the class's sequence
    /// of calls rather than a rule about the geometry.
    fn run_with(&mut self, i: usize, inlets: &[Stream]) -> Result<TrayOutcome> {
        let pressure = Some(self.pressures[i]);
        let out = match self.modes[i] {
            // The default: the tray's flash at the end's temperature, which the outer loop
            // moves while it searches. **A middle tray carries a NaN here and no pin**, so its
            // flash is at its own enthalpy - the same distinction the endpoints' temperatures
            // make everywhere else in this tier.
            EndMode::Temperature(t) => {
                let pin = if t.is_finite() {
                    Some(kelvins(t))
                } else {
                    None
                };
                stage::tray(inlets, pressure, pin, watts(0.0), self.reactive[i])?
            }
            // A heat input and no pin, which is what a DUTY specification writes.
            EndMode::Duty(duty) => {
                stage::tray(inlets, pressure, None, watts(duty), self.reactive[i])?
            }
            // The end's own reflux flash, which is what a REFLUX_RATIO specification writes -
            // `unit_ops.distillation_column`'s ends rather than its stage.
            EndMode::RefluxRatio(ratio) => {
                if i == 0 {
                    let out = crate::column::reboiler(
                        inlets,
                        pressure,
                        None,
                        watts(0.0),
                        crate::column::ReboilerMode::VaporBoilupRatio(ratio),
                    )?;
                    TrayOutcome {
                        temperature: out.temperature,
                        pressure: out.pressure,
                        gas: out.vapour,
                        liquid: out.liquid,
                        warnings: Vec::new(),
                    }
                } else {
                    let out = crate::column::condenser(
                        inlets,
                        pressure,
                        crate::column::CondenserMode::RefluxRatio(ratio),
                    )?;
                    TrayOutcome {
                        temperature: out.temperature,
                        pressure: out.pressure,
                        gas: out.distillate,
                        liquid: out.reflux,
                        warnings: Vec::new(),
                    }
                }
            }
        };
        for warning in &out.warnings {
            if !self.warnings.iter().any(|held| held.code == warning.code) {
                self.warnings.push(warning.clone());
            }
        }
        self.gas[i] = out.gas.clone();
        self.liquid[i] = out.liquid.clone();
        Ok(out)
    }

    /// A tray's temperature: its vapour's if it has one, its liquid's otherwise.
    fn temperature(&self, i: usize) -> f64 {
        self.gas[i]
            .as_ref()
            .or(self.liquid[i].as_ref())
            .map_or(f64::NAN, |s| s.t.value)
    }

    /// A tray's inlets' total enthalpy, W.
    fn inlet_enthalpy(&self, i: usize) -> f64 {
        self.inlets_of(i, None)
            .iter()
            .map(|s| s.n * s.h.value)
            .sum()
    }

    /// A tray's outlets' total enthalpy, W.
    fn outlet_enthalpy(&self, i: usize) -> f64 {
        self.gas[i].as_ref().map_or(0.0, |s| s.n * s.h.value)
            + self.liquid[i].as_ref().map_or(0.0, |s| s.n * s.h.value)
    }
}

/// Solve a distillation column by sequential substitution.
///
/// The arithmetic is `DistillationColumn`'s: `init` flashes the feed stage, seeds a **linear
/// temperature profile** towards both ends and links the trays upward and downward, and then
/// `solveSequential` sweeps - liquid down, vapour up, liquid down again - running each tray
/// from its neighbours' current outlets until the mean tray-temperature change falls below
/// the tolerance.
///
/// **The ends are pinned by temperature**, which is `setCondenserTemperature`'s mechanism: it
/// reaches the tray's own `outTemperature`, so its flash is a `TPflash` there rather than at
/// the mixed enthalpy. A middle tray has no pin and flashes at its enthalpy.
///
/// **The adaptive relaxation is the class's, and its controller reads the temperature
/// residual alone.** `DistillationColumn` combines it with the mass and the energy residuals
/// before deciding whether to damp; those two are MESH norms this port does not carry - they
/// are D5's subject - so the controller here sees one third of what the class's does. It
/// matters only when a solve struggles, and the rows this port is oracled on converge at a
/// relaxation of one, where the controller is the identity.
///
/// # Errors
/// [`azoth_core::AzothError::InvalidInput`] for a stage count or a feed stage outside the
/// column, and [`azoth_core::AzothError::SolverNotConverged`] when the gate is missed - with
/// the residuals in the message, because a column that stopped is not a column that solved.
pub fn distillation_column(setup: &ColumnSetup) -> Result<ColumnOutcome> {
    let tray_count =
        setup.number_of_stages + usize::from(setup.has_reboiler) + usize::from(setup.has_condenser);
    if setup.number_of_stages == 0 {
        return Err(AzothError::invalid_input(
            "number_of_stages",
            "a column with no stages between its ends is not a column",
        ));
    }
    if setup.feed_stage >= tray_count {
        return Err(AzothError::invalid_input(
            "feed_stage",
            format!(
                "stage {} is outside a column of {tray_count} trays, whose stages run 0 to {}",
                setup.feed_stage,
                tray_count - 1
            ),
        ));
    }
    if !setup.temperature_tolerance.is_finite() || setup.temperature_tolerance <= 0.0 {
        return Err(AzothError::invalid_input(
            "temperature_tolerance",
            format!(
                "a convergence tolerance is a positive number, and {} is not one",
                setup.temperature_tolerance
            ),
        ));
    }

    // Pressures run linearly from the bottom to the top, which is the class's own rule.
    let pressures: Vec<Pressure> = (0..tray_count)
        .map(|i| {
            pascals(
                setup.bottom_pressure.value
                    + (setup.top_pressure.value - setup.bottom_pressure.value) * i as f64
                        / (tray_count - 1) as f64,
            )
        })
        .collect();

    if setup.solver_type == SolverType::NaphtaliSandholm {
        return crate::column::naphtali_sandholm::solve(setup, &pressures);
    }
    // **An end's mode, and a direct specification is what changes it.** `applyDirectSpecification`
    // writes a reflux ratio or a duty onto the end itself, so those two types never reach the
    // outer loop; everything else is driven by moving the end's temperature.
    // **A pin wins over a directly-applied duty, and that is a measurement.** NeqSim's tray
    // `run` tests its stated outlet temperature first and takes a `TPflash` there, so a
    // condenser with both a pin and a `DUTY` specification never reads the heat input: the
    // capture's `spec_top_duty_under_a_pin_is_inert` row reports the pinned `-21323.04` W for a
    // target of `-20000`. A reflux ratio *does* win, because the end's own `run` tests
    // `refluxIsSet` before anything else. Reproducing the first is the port rule; quietly
    // honouring the duty would be an improvement NeqSim does not make.
    let end_mode = |specification: Option<&Specification>, pin: Option<f64>| -> EndMode {
        match specification.map(|s| (s.kind, s.target)) {
            Some((SpecificationKind::RefluxRatio, ratio)) => EndMode::RefluxRatio(ratio),
            Some((SpecificationKind::Duty, duty)) if pin.is_none() => EndMode::Duty(duty),
            _ => EndMode::Temperature(pin.unwrap_or(f64::NAN)),
        }
    };
    // **A per-tray pin is the tray's own stated outlet temperature**, and it is applied where
    // it is finite: `SimpleTray.setOutletTemperature` is one mechanism, whether a caller
    // reaches it tray by tray or through the ends' setters.
    let tray_pin = |i: usize| -> Option<f64> {
        setup
            .tray_temperatures
            .as_ref()
            .and_then(|pins| pins.get(i).copied())
            .filter(|temperature| temperature.is_finite())
    };
    if let Some(pins) = setup.tray_temperatures.as_ref() {
        if pins.len() != tray_count {
            return Err(AzothError::invalid_input(
                "tray_temperatures",
                format!(
                    "a column of {tray_count} trays needs one pin each, and {} were stated",
                    pins.len()
                ),
            ));
        }
    }
    let modes: Vec<EndMode> = (0..tray_count)
        .map(|i| {
            if let Some(temperature) = tray_pin(i) {
                return EndMode::Temperature(temperature);
            }
            if i == 0 && setup.has_reboiler {
                end_mode(
                    setup.bottom_specification.as_ref(),
                    setup.reboiler_temperature.map(|t| t.value),
                )
            } else if i == tray_count - 1 && setup.has_condenser {
                end_mode(
                    setup.top_specification.as_ref(),
                    setup.condenser_temperature.map(|t| t.value),
                )
            } else {
                EndMode::Temperature(f64::NAN)
            }
        })
        .collect();

    // **The section is stated over middle trays and the ends are never in it.** `replaceMiddleTrays`
    // walks from `hasReboiler ? 1 : 0` to `hasCondenser ? size - 1 : size`, so a middle index is
    // an absolute index less the reboiler's one - which is the whole of the conversion.
    let middle_offset = usize::from(setup.has_reboiler);
    let reactive: Vec<bool> = (0..tray_count)
        .map(|i| {
            let is_end =
                (i == 0 && setup.has_reboiler) || (i + 1 == tray_count && setup.has_condenser);
            if is_end {
                return false;
            }
            setup.reactive.covers(i - middle_offset)
        })
        .collect();

    let mut net = Network {
        reactive,
        warnings: Vec::new(),
        gas: vec![None; tray_count],
        liquid: vec![None; tray_count],
        pressures,
        feed_stage: setup.feed_stage,
        tray_count,
        feed: setup.feed.clone(),
        top_feed: setup.top_feed.clone(),
        modes,
    };

    let mut temperatures = seed_network(&mut net, setup)?;
    let (iterations, temperature_residual) = if adjustable_specification(setup) {
        specification_loop(&mut net, setup, &mut temperatures)?
    } else {
        sweep(&mut net, setup, &mut temperatures)?
    };

    // ---- The products, the duties and the closure ----
    let distillate = end_product(&net, tray_count - 1, "condenser")?;
    let bottoms = end_product(&net, 0, "reboiler")?;
    let reboiler_duty = if setup.has_reboiler {
        net.outlet_enthalpy(0) - net.inlet_enthalpy(0)
    } else {
        0.0
    };
    let condenser_duty = if setup.has_condenser {
        net.outlet_enthalpy(tray_count - 1) - net.inlet_enthalpy(tray_count - 1)
    } else {
        0.0
    };

    let feeds: Vec<&Stream> = std::iter::once(&setup.feed)
        .chain(setup.top_feed.iter())
        .collect();
    let feed_enthalpy: f64 = feeds.iter().map(|feed| feed.n * feed.h.value).sum();
    let products_enthalpy = distillate.n * distillate.h.value + bottoms.n * bottoms.h.value;
    let energy_residual = if feed_enthalpy.abs() > 0.0 {
        (feed_enthalpy + reboiler_duty + condenser_duty - products_enthalpy).abs()
            / feed_enthalpy.abs()
    } else {
        0.0
    };

    // The component closure: the worst component's imbalance against the feed, relative.
    let mut mass_residual = 0.0_f64;
    for (c, ci) in setup.feed.z.iter().enumerate() {
        let supplied: f64 = feeds
            .iter()
            .map(|feed| feed.n * feed.z.get(c).copied().unwrap_or(*ci))
            .sum();
        let delivered =
            distillate.n * fraction_of(&distillate, c) + bottoms.n * fraction_of(&bottoms, c);
        if supplied.abs() > 1.0e-12 {
            mass_residual = mass_residual.max((supplied - delivered).abs() / supplied.abs());
        }
    }

    let trays: Vec<TrayProfile> = (0..tray_count)
        .map(|i| TrayProfile {
            temperature: kelvins(net.temperature(i)),
            pressure: net.pressures[i],
            gas_n: net.gas[i].as_ref().map_or(0.0, |s| s.n),
            liquid_n: net.liquid[i].as_ref().map_or(0.0, |s| s.n),
            gas_z: net.gas[i].as_ref().map_or_else(Vec::new, |s| s.z.clone()),
            liquid_z: net.liquid[i]
                .as_ref()
                .map_or_else(Vec::new, |s| s.z.clone()),
        })
        .collect();

    if temperature_residual > setup.temperature_tolerance
        || mass_residual > MASS_BALANCE_TOLERANCE
        || energy_residual > ENTHALPY_BALANCE_TOLERANCE
    {
        // A refusal, and the three residuals travel with it: a solve that stopped is a state
        // a caller has to be able to diagnose, and the mass and the energy closure are what
        // say whether it stopped *near* a column.
        return Err(AzothError::SolverNotConverged {
            iterations,
            residual: temperature_residual.max(mass_residual).max(energy_residual),
            tolerance: setup
                .temperature_tolerance
                .min(MASS_BALANCE_TOLERANCE)
                .min(ENTHALPY_BALANCE_TOLERANCE),
        });
    }

    Ok(ColumnOutcome {
        trays,
        distillate,
        bottoms,
        condenser_duty: watts(condenser_duty),
        reboiler_duty: watts(reboiler_duty),
        iterations,
        temperature_residual,
        mass_residual,
        energy_residual,
        warnings: net.warnings.clone(),
    })
}

/// The product leaving an end: the vapour for a condenser, the liquid for a reboiler.
fn end_product(net: &Network, i: usize, which: &str) -> Result<Stream> {
    let stream = if i == 0 {
        net.liquid[i].clone()
    } else {
        net.gas[i].clone().or_else(|| net.liquid[i].clone())
    };
    stream.ok_or_else(|| {
        AzothError::invalid_input(
            "column",
            format!("the {which} has no product, which is a column whose ends did not solve"),
        )
    })
}

/// A stream re-stated at another temperature, at its own pressure and composition.
///
/// The tray `init` seeding: it changes an inlet's *state* without changing what it carries.
fn at_temperature(stream: &Stream, temperature: f64) -> Stream {
    Stream::from_pt(
        stream.components.clone(),
        stream.z.clone(),
        stream.n,
        stream.p,
        kelvins(temperature),
    )
    .unwrap_or_else(|_| stream.clone())
}

/// A stream's mole fraction of one component.
fn fraction_of(stream: &Stream, component: usize) -> f64 {
    stream.z.get(component).copied().unwrap_or(0.0)
}

/// `DistillationColumn.init`: the feed-stage flash, the linear temperature seed and the two
/// links.
///
/// Returns the seeded profile, which the sweeps then evolve. **The seed is the class's own
/// arithmetic**: a step towards each end, the condenser's divided by the trays above the feed
/// and the reboiler's by the trays below it, which is why the two differ.
fn seed_network(net: &mut Network, setup: &ColumnSetup) -> Result<Vec<f64>> {
    let first_feed = setup.feed_stage;
    let tray_count = net.tray_count;
    // The feed tray, from its feed alone.
    net.run_with(first_feed, std::slice::from_ref(&setup.feed))?;

    // The reboiler, against the *feed stage's* liquid and not the first tray's - which is
    // `init`'s own link (`trays.get(0).addStream(getTray(firstFeedTrayNumber).getLiquidOutStream())`)
    // and a different inlet from the one the sweeps give it.
    if setup.has_reboiler {
        let from_feed = net.liquid[first_feed]
            .clone()
            .expect("the feed stage has a liquid after its flash");
        net.run_with(0, &[from_feed])?;
    }

    // The linear temperature profile: the feed tray's temperature with a step towards each
    // end. The condenser's step divides by the trays above the feed and the reboiler's by the
    // trays below it, which is the class's own arithmetic and is why the two steps differ.
    let feed_temperature = net.temperature(first_feed);
    let pin_of = |mode: EndMode| match mode {
        EndMode::Temperature(t) if t.is_finite() => Some(t),
        _ => None,
    };
    let condenser_temperature = pin_of(net.modes[tray_count - 1]).unwrap_or(feed_temperature - 1.0);
    let reboiler_temperature = net
        .liquid
        .first()
        .and_then(Option::as_ref)
        .map_or(feed_temperature, |s| s.t.value);
    let mut temperatures = vec![f64::NAN; tray_count];
    temperatures[first_feed] = feed_temperature;
    let delta_up =
        (feed_temperature - condenser_temperature) / (tray_count as f64 - first_feed as f64 - 1.0);
    let delta_down = (reboiler_temperature - feed_temperature) / first_feed as f64;
    let mut delta = 0.0;
    for temperature in temperatures.iter_mut().skip(first_feed + 1) {
        delta += delta_up;
        *temperature = feed_temperature - delta;
    }
    let mut delta = 0.0;
    for temperature in temperatures[..first_feed].iter_mut().rev() {
        delta += delta_down;
        *temperature = feed_temperature + delta;
    }
    for (i, mode) in net.modes.iter().enumerate() {
        if let Some(t) = pin_of(*mode) {
            temperatures[i] = t;
        }
    }

    // Link upward, then downward, then the reboiler against the first tray's liquid. Each
    // link is seeded to the tray's own temperature, which is `SimpleTray.init()`'s mechanism
    // - it re-states every inlet the tray holds, the external feed included, and the class
    // restores the feed afterwards, which is why the sweeps below see the feed's own state.
    for (i, &temperature) in temperatures.iter().enumerate().skip(1) {
        let mut inlets = vec![at_temperature(
            net.gas[i - 1].as_ref().expect("the tray below has run"),
            temperature,
        )];
        // **The external feeds are never re-stated**, which is what the class's own
        // `refreshInternalExternalFeedSystems` restores - and a tray can carry two of them.
        inlets.extend(net.feeds_at(i).into_iter().cloned());
        net.run_with(i, &inlets)?;
    }
    for (i, &temperature) in temperatures
        .iter()
        .enumerate()
        .take(tray_count - 1)
        .skip(1)
        .rev()
    {
        let mut inlets = vec![at_temperature(
            net.gas[i - 1].as_ref().expect("the tray below has run"),
            temperature,
        )];
        inlets.push(at_temperature(
            net.liquid[i + 1].as_ref().expect("the tray above has run"),
            temperature,
        ));
        inlets.extend(net.feeds_at(i).into_iter().cloned());
        net.run_with(i, &inlets)?;
    }
    if setup.has_reboiler {
        let inlets = vec![at_temperature(
            net.liquid[1].as_ref().expect("the first tray has run"),
            temperatures[0],
        )];
        net.run_with(0, &inlets)?;
    }
    Ok(temperatures)
}

/// `DistillationColumn.solveSequential`: the sweeps, once, to the temperature gate.
///
/// Returns the iterations taken and the mean tray-temperature change at the last one. A call
/// after an earlier one starts from where the last stopped, which is what the class's outer
/// specification loop relies on - it passes `setDoInitializion(false)` after its first pass and
/// keeps the profile.
fn sweep(net: &mut Network, setup: &ColumnSetup, temperatures: &mut [f64]) -> Result<(u32, f64)> {
    let first_feed = setup.feed_stage;
    let tray_count = net.tray_count;
    let mut relaxation = 1.0_f64;
    let mut previous_combined = f64::INFINITY;
    let mut temperature_residual = f64::INFINITY;
    let mut iterations = 0_u32;

    for iter in 1..=setup.max_iterations {
        iterations = iter as u32;
        let old: Vec<f64> = (0..tray_count).map(|i| net.temperature(i)).collect();

        // Down the column: each tray takes the liquid from the one above it.
        for i in (2..=first_feed).rev() {
            net.run(i - 1, None)?;
        }
        // The reboiler takes the first stage's liquid.
        if setup.has_reboiler {
            net.run(0, None)?;
        }
        // Up the column: each tray takes the vapour from the one below it.
        for i in 1..tray_count {
            net.run(i, None)?;
        }
        // And down again, from below the condenser to the feed stage.
        for i in (first_feed..tray_count.saturating_sub(1)).rev() {
            net.run(i, None)?;
        }

        // The temperature update, with the class's own floor on the step.
        let effective = relaxation.clamp(MIN_TEMPERATURE_RELAXATION, 1.0);
        let mut sum = 0.0;
        for i in 0..tray_count {
            let updated = net.temperature(i);
            let updated = if updated.is_finite() { updated } else { old[i] };
            sum += (updated - old[i]).abs();
            temperatures[i] = old[i] + effective * (updated - old[i]);
        }
        temperature_residual = sum / tray_count as f64;

        // The adaptive controller, on the combined residual the class scales each term of.
        let combined = temperature_residual / setup.temperature_tolerance;
        if combined > previous_combined * 1.05 {
            relaxation = (relaxation * RELAXATION_DECREASE_FACTOR).max(MIN_SEQUENTIAL_RELAXATION);
        } else if combined < previous_combined * 0.98 {
            relaxation = (relaxation * RELAXATION_INCREASE_FACTOR).min(MAX_ADAPTIVE_RELAXATION);
        }
        previous_combined = combined;

        if temperature_residual <= setup.temperature_tolerance {
            break;
        }
    }
    Ok((iterations, temperature_residual))
}

/// Whether any specification needs the outer loop, which is `hasAdjustableSpecifications`.
fn adjustable_specification(setup: &ColumnSetup) -> bool {
    [&setup.top_specification, &setup.bottom_specification]
        .into_iter()
        .flatten()
        .any(|spec| spec.kind.is_adjusted())
}

/// The class's own defaults for a specification's convergence, from `ColumnSpecification`.
const SPECIFICATION_TOLERANCE: f64 = 1.0e-4;
const SPECIFICATION_MAX_ITERATIONS: usize = 20;
/// `secantStep`'s step cap and its physically reasonable window.
const SECANT_MAX_STEP: f64 = 50.0;
const SECANT_MIN_TEMPERATURE: f64 = 100.0;
const SECANT_MAX_TEMPERATURE: f64 = 1000.0;
/// The first iteration's second guess, an offset from the first - `topTemp1 = topTemp0 - 5.0`.
const SECANT_FIRST_OFFSET: f64 = 5.0;

/// `solveWithSpecificationTargets`: the outer secant on the ends' temperatures.
///
/// **One specification at each end, each driven by its own secant**, and both ends solved
/// together at every outer iteration: the class changes the condenser's and the reboiler's
/// temperatures and re-runs the whole column, then reads both errors off its products.
///
/// **There is no homotopy here, and that is a measurement.** `getEffectiveSpecificationHomotopySteps`
/// returns three stages only when `solverType == AUTO`, and the field it otherwise reads is
/// one - so the staged continuation the class does carry belongs to the solver this port
/// refuses, and the secant loop is the whole of the specification machinery it needs.
fn specification_loop(
    net: &mut Network,
    setup: &ColumnSetup,
    temperatures: &mut [f64],
) -> Result<(u32, f64)> {
    let tray_count = net.tray_count;
    let top = setup.top_specification.clone();
    let bottom = setup.bottom_specification.clone();
    let adjust_top = top.as_ref().is_some_and(|s| s.kind.is_adjusted());
    let adjust_bottom = bottom.as_ref().is_some_and(|s| s.kind.is_adjusted());

    // The guesses, from the stated temperatures - the class's own starting point.
    // `estimateFeedTemperature` plus the class's own 20 K offsets where an end is unpinned.
    let feed_temperature = setup.feed.t.value;
    let top_start = setup
        .condenser_temperature
        .map_or(feed_temperature - 20.0, |t| t.value);
    let bottom_start = setup
        .reboiler_temperature
        .map_or(feed_temperature + 20.0, |t| t.value);
    let mut top_pair = (top_start, top_start - SECANT_FIRST_OFFSET);
    let mut bottom_pair = (bottom_start, bottom_start + SECANT_FIRST_OFFSET);
    let mut top_errors = (f64::NAN, f64::NAN);
    let mut bottom_errors = (f64::NAN, f64::NAN);

    let mut iterations = 0_u32;
    let mut temperature_residual = f64::INFINITY;
    for outer in 0..SPECIFICATION_MAX_ITERATIONS {
        if adjust_top {
            let guess = if outer == 0 { top_pair.0 } else { top_pair.1 };
            net.modes[tray_count - 1] = EndMode::Temperature(guess);
        }
        if adjust_bottom {
            let guess = if outer == 0 {
                bottom_pair.0
            } else {
                bottom_pair.1
            };
            net.modes[0] = EndMode::Temperature(guess);
        }
        let (taken, residual) = sweep(net, setup, temperatures)?;
        // **The *last* inner solve's count, which is the class's own `lastIterationCount`** -
        // it is overwritten by every `solveSequential` call and the outer loop makes twenty of
        // them. Summing them here would report a number the capture has no counterpart for.
        iterations = taken;
        temperature_residual = residual;

        let top_error = match (&top, adjust_top) {
            (Some(spec), true) => {
                let product = net.gas[tray_count - 1].clone().ok_or_else(no_product)?;
                spec.value(&product, &setup.feed)? - spec.target
            }
            _ => 0.0,
        };
        let bottom_error = match (&bottom, adjust_bottom) {
            (Some(spec), true) => {
                let product = net.liquid[0].clone().ok_or_else(no_product)?;
                spec.value(&product, &setup.feed)? - spec.target
            }
            _ => 0.0,
        };

        let converged = (!adjust_top || top_error.abs() < SPECIFICATION_TOLERANCE)
            && (!adjust_bottom || bottom_error.abs() < SPECIFICATION_TOLERANCE);
        if converged {
            return Ok((iterations, temperature_residual));
        }

        if adjust_top {
            if outer == 0 {
                top_errors.0 = top_error;
            } else {
                top_errors.1 = top_error;
                let next = secant_step(top_pair, top_errors);
                top_pair = (top_pair.1, next);
                top_errors = (top_errors.1, f64::NAN);
            }
        }
        if adjust_bottom {
            if outer == 0 {
                bottom_errors.0 = bottom_error;
            } else {
                bottom_errors.1 = bottom_error;
                let next = secant_step(bottom_pair, bottom_errors);
                bottom_pair = (bottom_pair.1, next);
                bottom_errors = (bottom_errors.1, f64::NAN);
            }
        }
    }
    Err(AzothError::SolverNotConverged {
        iterations,
        residual: temperature_residual,
        tolerance: SPECIFICATION_TOLERANCE,
    })
}

/// `secantStep`, with its two guards: a step cap and a physically reasonable window.
fn secant_step(guess: (f64, f64), error: (f64, f64)) -> f64 {
    let denominator = error.1 - error.0;
    let mut next = if denominator.abs() < 1.0e-15 {
        guess.1 + 2.0
    } else {
        guess.1 - error.1 * (guess.1 - guess.0) / denominator
    };
    if next - guess.1 > SECANT_MAX_STEP {
        next = guess.1 + SECANT_MAX_STEP;
    } else if next - guess.1 < -SECANT_MAX_STEP {
        next = guess.1 - SECANT_MAX_STEP;
    }
    next.clamp(SECANT_MIN_TEMPERATURE, SECANT_MAX_TEMPERATURE)
}

/// An end with no product is an end whose solve did not produce one.
fn no_product() -> AzothError {
    AzothError::invalid_input(
        "column",
        "an end has no product, which is a column whose ends did not solve",
    )
}
