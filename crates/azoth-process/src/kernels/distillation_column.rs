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
    /// The condenser's temperature, which pins the top tray.
    pub condenser_temperature: ThermodynamicTemperature,
    /// The reboiler's temperature, which pins the bottom tray.
    pub reboiler_temperature: ThermodynamicTemperature,
    /// The convergence tolerance on the mean tray-temperature change.
    pub temperature_tolerance: f64,
    /// The iteration cap.
    pub max_iterations: usize,
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
    /// The feed.
    feed: Stream,
    /// The pinned temperature of each end, if it has one.
    pins: Vec<Option<ThermodynamicTemperature>>,
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
        if i == self.feed_stage {
            inlets.push(self.feed.clone());
        }
        inlets
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
        let out = stage::tray(inlets, Some(self.pressures[i]), self.pins[i], watts(0.0))?;
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
    let pins: Vec<Option<ThermodynamicTemperature>> = (0..tray_count)
        .map(|i| {
            if i == 0 && setup.has_reboiler {
                Some(setup.reboiler_temperature)
            } else if i == tray_count - 1 && setup.has_condenser {
                Some(setup.condenser_temperature)
            } else {
                None
            }
        })
        .collect();

    let mut net = Network {
        gas: vec![None; tray_count],
        liquid: vec![None; tray_count],
        pressures,
        feed_stage: setup.feed_stage,
        tray_count,
        feed: setup.feed.clone(),
        pins,
    };
    let first_feed = setup.feed_stage;

    // ---- `init` ----
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
    let condenser_temperature =
        net.pins[tray_count - 1].map_or(feed_temperature - 1.0, |t| t.value);
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
    for (i, pin) in net.pins.iter().enumerate() {
        if let Some(t) = pin {
            temperatures[i] = t.value;
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
        if i == first_feed {
            inlets.push(setup.feed.clone());
        }
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
        if i == first_feed {
            inlets.push(setup.feed.clone());
        }
        net.run_with(i, &inlets)?;
    }
    if setup.has_reboiler {
        let inlets = vec![at_temperature(
            net.liquid[1].as_ref().expect("the first tray has run"),
            temperatures[0],
        )];
        net.run_with(0, &inlets)?;
    }

    // ---- `solveSequential` ----
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

    let feed_enthalpy = setup.feed.n * setup.feed.h.value;
    let products_enthalpy = distillate.n * distillate.h.value + bottoms.n * bottoms.h.value;
    let energy_residual = if feed_enthalpy.abs() > 0.0 {
        (feed_enthalpy + reboiler_duty + condenser_duty - products_enthalpy).abs()
            / feed_enthalpy.abs()
    } else {
        0.0
    };

    // The component closure: the worst component's imbalance against the feed, relative.
    let mut mass_residual = 0.0_f64;
    for (c, &zi) in setup.feed.z.iter().enumerate() {
        let supplied = setup.feed.n * zi;
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
