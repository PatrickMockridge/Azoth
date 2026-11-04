//! What a run publishes beside its streams.
//!
//! **A kernel that ran a flash, a column solve or a reactor integration reached numbers that are
//! on no outlet stream**, and until the executor carried them a front end could draw a flowsheet
//! and not read one - which is what the *editor* tranche found: the duty was computed and dropped
//! at the call site. `KernelOutcome` is what carries it, and this file is the gate over it.
//!
//! **Two things are checked, and they are different in kind.** The first is a table against the
//! dispatch table in both directions, so an entry that starts publishing without a row here - or a
//! row whose entry stops - fails rather than passing quietly. The second runs the *shipped*
//! flowsheet and reads the wire, because a table of field names agrees with itself whether or not
//! a number ever crosses.
//!
//! **What is deliberately not here**: one row per unit operation. A unit operation's own numbers
//! are pinned by its case (`specs/models/process/*.toml`) and by the kernel tests beside this
//! file; what this file holds is that the *run* publishes them and names them what the model does.

use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

use azoth_core::CalcResult;
use azoth_process::executor::json::{SessionReport, session_report};
use azoth_process::executor::{DISPATCH, Session};
use azoth_process::models;
use azoth_process::{ExecutionOrder, Flowsheet, load_palette};

fn root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../..")
}

fn run_flowsheet(text: &str) -> SessionReport {
    let flowsheet: Flowsheet = toml::from_str(text).expect("it parses");
    let palette = load_palette(&root().join("specs/unit_ops")).expect("the palette loads");
    let session = Session::run(
        flowsheet,
        &palette,
        &BTreeMap::new(),
        ExecutionOrder::Insertion,
    )
    .expect("the flowsheet runs");
    session_report(&session).expect("it writes")
}

/// **What each dispatched entry publishes**, as the model's own `FIELDS`.
///
/// `None` is a real answer and not a gap: a mixer's result *is* its mixed stream and a separator's
/// *are* its two, so a run that published them here would be the same numbers under a second key.
/// `python/src/azoth/process/layers.py`'s `NO_INTERIOR` table names the models that way - that
/// their kernel forms nothing their records do not carry - and six of them are entries that
/// publish no record here. **`process.separator` left that table when its capacity surface
/// arrived**: utilisation and the three velocities are numbers neither outlet record carries.
fn published() -> Vec<(&'static str, Option<&'static [&'static str]>)> {
    vec![
        (
            "unit_ops.absorption_column",
            Some(models::AbsorptionColumnResult::FIELDS),
        ),
        ("unit_ops.component_splitter", None),
        ("unit_ops.compressor", None),
        ("unit_ops.cooler", Some(models::CoolerResult::FIELDS)),
        (
            "unit_ops.distillation_column",
            Some(models::DistillationColumnResult::FIELDS),
        ),
        // **The mixing pressure and the four velocities cross here**, none of them an input and
        // none on the discharge state: `Option`s for the one case where no fluid flowed.
        ("unit_ops.ejector", Some(models::EjectorResult::FIELDS)),
        ("unit_ops.expander", None),
        ("unit_ops.filter", Some(models::FilterResult::FIELDS)),
        ("unit_ops.flare", Some(models::FlareResult::FIELDS)),
        ("unit_ops.gas_scrubber", None),
        (
            "unit_ops.gibbs_reactor",
            Some(models::GibbsReactorResult::FIELDS),
        ),
        // **The duty and the rating cross here**, and the four rating numbers are `Option`s
        // because a pinned-outlet run reaches none of them: the record writes `null` where a
        // rating did not happen, which is the same absence a single-phase flash is.
        (
            "unit_ops.heat_exchanger",
            Some(models::HeatExchangerResult::FIELDS),
        ),
        ("unit_ops.heater", Some(models::HeaterResult::FIELDS)),
        ("unit_ops.manifold", None),
        // The map machine's executor arm runs the specified-pressure mode with no chart, so it
        // publishes its two outlets and nothing else; the speed and the powers are the model's
        // record, not the wire's.
        ("unit_ops.map_turbo_expander_compressor", None),
        ("unit_ops.mixer", None),
        ("unit_ops.pipe", Some(models::PipeResult::FIELDS)),
        (
            "unit_ops.plug_flow_reactor",
            Some(models::PlugFlowReactorResult::FIELDS),
        ),
        // **The packed column's entry became runnable**, so it appears here for the first time:
        // it dispatches the base column's solve, and its record is not on the wire either - the
        // same statement as its sibling's below.
        ("unit_ops.packed_column", None),
        ("unit_ops.pump", None),
        // **The rate-based packed column reaches an interior and publishes all of it**: its
        // executor arm calls the model's own `RateBasedPackedColumnResult::of`, so the wire
        // carries the segment profile rather than only the two outlets.
        (
            "unit_ops.rate_based_packed_column",
            Some(models::RateBasedPackedColumnResult::FIELDS),
        ),
        ("unit_ops.separator", None),
        (
            "unit_ops.shortcut_distillation_column",
            Some(models::ShortcutDistillationColumnResult::FIELDS),
        ),
        // **The capacity map's six optional fields cross here**, none of them on an outlet
        // stream; the two summaries and the products are the model's own record as well.
        ("unit_ops.splitter", Some(models::SplitterResult::FIELDS)),
        (
            "unit_ops.stirred_tank_reactor",
            Some(models::StirredTankReactorResult::FIELDS),
        ),
        (
            "unit_ops.stripping_column",
            Some(models::StrippingColumnResult::FIELDS),
        ),
        ("unit_ops.tank", None),
        ("unit_ops.three_phase_separator", None),
        ("unit_ops.throttling_valve", None),
        // The machine's executor arm publishes its two outlets and nothing else, like its
        // neighbours' - the speed and the powers are the model's record, not the wire's.
        ("unit_ops.turbo_expander_compressor", None),
    ]
}

/// **The table above and the dispatch table are the same table.** Both directions, because either
/// one alone leaves a hole: an entry added to `DISPATCH` with no row here, and a row here for an
/// entry that no longer dispatches.
#[test]
fn the_published_table_covers_the_dispatch_table_both_ways() {
    let table = published();
    let named: Vec<&str> = table.iter().map(|(id, _)| *id).collect();
    let dispatched: Vec<&str> = DISPATCH.iter().map(|(id, _)| *id).collect();
    for id in &dispatched {
        assert!(
            named.contains(id),
            "`{id}` dispatches and this file does not say what it publishes"
        );
    }
    for id in &named {
        assert!(
            dispatched.contains(id),
            "this file says what `{id}` publishes and it does not dispatch"
        );
    }
    assert_eq!(table.len(), dispatched.len());
}

/// **The shipped flowsheet, on the wire.**
///
/// The heater's duty is the number no outlet stream carries, and it is the one that made this
/// work necessary: `dispatch.rs` read `Ok(vec![out.outlet])` and the duty went out of scope. The
/// mixer beside it is the counter-case - it publishes nothing - so an implementation that
/// published a restatement of every stream record would fail here.
#[test]
fn the_shipped_demo_publishes_the_duty_no_stream_carries() {
    let text = std::fs::read_to_string(root().join("specs/flowsheets/demo.toml"))
        .expect("the shipped flowsheet is there");
    let flowsheet: Flowsheet = toml::from_str(&text).expect("it parses");
    let palette = load_palette(&root().join("specs/unit_ops")).expect("the palette loads");
    let session = Session::run(
        flowsheet,
        &palette,
        &BTreeMap::new(),
        ExecutionOrder::Insertion,
    )
    .expect("the flowsheet runs");
    let report: SessionReport = session_report(&session).expect("it writes");

    let heater = report
        .results
        .get("hx1")
        .expect("the heater reached a duty and it crosses");
    assert_eq!(heater["outlet_duty"]["unit"], "W");
    let duty = heater["outlet_duty"]["magnitude_si"]
        .as_f64()
        .expect("a magnitude");
    assert!(duty.is_finite() && duty != 0.0, "the duty is {duty}");

    // **Every key is one the model declares**, so the wire cannot grow a spelling of its own.
    let keys: Vec<&str> = heater
        .as_object()
        .expect("an object")
        .keys()
        .map(String::as_str)
        .collect();
    let mut declared: Vec<&str> = models::HeaterResult::FIELDS.to_vec();
    let mut found: Vec<&str> = keys.clone();
    declared.sort_unstable();
    found.sort_unstable();
    assert_eq!(found, declared);

    // The demo states no outlet-temperature limits, so the map is empty and its two
    // per-constraint fields are `null`.
    assert!(heater["max_outlet_temperature_utilization"].is_null());
    assert!(heater["min_outlet_temperature_utilization"].is_null());
    assert!(heater["bottleneck_constraint"].is_null());
    assert_eq!(heater["max_utilization"].as_f64(), Some(0.0));
    assert_eq!(heater["capacity_exceeded"].as_bool(), Some(false));
    assert_eq!(heater["hard_limit_exceeded"].as_bool(), Some(false));

    // A mixer's whole answer is its stream, so it publishes nothing rather than a restatement.
    assert!(
        !report.results.contains_key("mix1"),
        "the mixer published {:?}",
        report.results.get("mix1")
    );

    // **And the record carries what the run knew and the five declared fields do not.** The
    // separator's vapour fraction is the flash's own beta, so a null here would mean the
    // constructor threw the flash's answer away again.
    let vapour = &report.streams["sep1.vapour"];
    assert!(vapour.vapour_fraction.is_some(), "the flash answered");
    assert!(vapour.mass_flow.is_some(), "two resolved components weigh");
    assert_eq!(
        vapour.mass_flow.as_ref().expect("weighed").unit,
        "kg/s",
        "the mass flow crosses in the vocabulary's own unit"
    );
    assert!(vapour.molar_mass.is_some());
}

/// **The heater's bounded outlet temperatures cross the wire**, and the hard flag is the
/// minimum's: the class reads `maxValue` for the maximum, which its setter never sets.
#[test]
fn the_heater_publishes_its_outlet_temperature_constraints() {
    const BOUNDED: &str = r#"
id = "flowsheets.heater_capacity"
name = "One heater with outlet-temperature limits"

products = ["out"]

[[inputs]]
name = "feed_1"
components = ["methane", "n-butane"]
n = 1.0
z = [0.9, 0.1]
P = 3.0e6
T = 320.0

[[instances]]
id = "h1"
unit = "unit_ops.heater"
[instances.parameters]
outlet_temperature = 380.0
max_outlet_temperature = 400.0
min_outlet_temperature = 300.0

[[connections]]
from = "feed_1"
to = "h1.inlet"

[[connections]]
from = "h1.outlet"
to = "out"
"#;
    let report = run_flowsheet(BOUNDED);
    let heater = &report.results["h1"];
    assert_eq!(
        heater["max_outlet_temperature_utilization"].as_f64(),
        Some(0.95)
    );
    assert!(
        (heater["min_outlet_temperature_utilization"]
            .as_f64()
            .expect("a utilization")
            - 300.0 / 380.0)
            .abs()
            < 1e-12
    );
    assert_eq!(heater["max_utilization"].as_f64(), Some(0.95));
    assert_eq!(heater["capacity_exceeded"].as_bool(), Some(false));
    assert_eq!(heater["hard_limit_exceeded"].as_bool(), Some(false));
    assert_eq!(heater["bottleneck_constraint"], "maxOutletTemperature");

    // The minimum above the outlet: exceeded, and this is the constraint that trips the
    // class's hard flag.
    let hard = BOUNDED.replace(
        "min_outlet_temperature = 300.0",
        "min_outlet_temperature = 400.0",
    );
    let report = run_flowsheet(&hard);
    let heater = &report.results["h1"];
    assert_eq!(heater["capacity_exceeded"].as_bool(), Some(true));
    assert_eq!(heater["hard_limit_exceeded"].as_bool(), Some(true));
    assert_eq!(heater["bottleneck_constraint"], "minOutletTemperature");
}

/// **The splitter's capacity map crosses the wire, and its absence is a `null`.**
///
/// The map's six per-constraint fields are `Option`s because the class's constraints are lazy:
/// a fresh `Splitter` has none. One flowsheet sizes the header and enables analysis, so the
/// velocity crosses; the other leaves the four inputs out, so the six are `null` and the two
/// summaries are the class's empty-map answers. The keys are `FIELDS` exactly, as for the
/// heater.
#[test]
fn the_splitter_publishes_its_capacity_map() {
    const SIZED: &str = r#"
id = "flowsheets.splitter_capacity"
name = "One splitter with a sized header"

products = ["light", "heavy"]

[[inputs]]
name = "feed_1"
components = ["n-butane", "n-pentane"]
n = 1.0
z = [0.6, 0.4]
P = 1.0e6
T = 300.0

[[instances]]
id = "split1"
unit = "unit_ops.splitter"
[instances.parameters]
split_factors = [0.3, 0.7]
header_diameter = 0.02
capacity_analysis_enabled = true

[[connections]]
from = "feed_1"
to = "split1.feed"

[[connections]]
from = "split1.products[0]"
to = "light"

[[connections]]
from = "split1.products[1]"
to = "heavy"
"#;
    let report = run_flowsheet(SIZED);
    let split = &report.results["split1"];

    let mut found: Vec<&str> = split
        .as_object()
        .expect("an object")
        .keys()
        .map(String::as_str)
        .collect();
    let mut declared: Vec<&str> = models::SplitterResult::FIELDS.to_vec();
    found.sort_unstable();
    declared.sort_unstable();
    assert_eq!(found, declared);

    assert_eq!(split["pressure_drop_value"]["unit"], "Pa");
    assert_eq!(
        split["pressure_drop_value"]["magnitude_si"].as_f64(),
        Some(0.0)
    );
    assert_eq!(split["velocity_value"]["unit"], "m/s");
    assert_eq!(
        split["velocity_design"]["magnitude_si"].as_f64(),
        Some(30.0)
    );
    let velocity = split["velocity_value"]["magnitude_si"]
        .as_f64()
        .expect("a magnitude");
    assert!(
        (velocity - 0.34597673984807986).abs() < 1e-12,
        "the velocity is {velocity}"
    );
    assert_eq!(split["capacity_exceeded"].as_bool(), Some(false));
    assert_eq!(split["bottleneck_constraint"], "velocity");
    let products: Vec<f64> = split["products_n"]
        .as_array()
        .expect("a vector")
        .iter()
        .map(|n| n.as_f64().expect("a flow"))
        .collect();
    assert_eq!(products.len(), 2);
    assert!((products[0] - 0.3).abs() < 1e-12 && (products[1] - 0.7).abs() < 1e-12);

    // The same flowsheet with the four capacity inputs left out: the class's fresh Splitter.
    let empty = SIZED.replace(
        "header_diameter = 0.02\ncapacity_analysis_enabled = true\n",
        "",
    );
    let report = run_flowsheet(&empty);
    let split = &report.results["split1"];
    for field in [
        "pressure_drop_value",
        "pressure_drop_design",
        "pressure_drop_utilization",
        "velocity_value",
        "velocity_design",
        "velocity_utilization",
        "bottleneck_constraint",
        "design_header_diameter",
        "design_branch_diameter",
        "design_header_length",
        "design_wall_thickness",
        "manifold_volume",
        "design_total_pressure_drop",
        "header_weight",
        "branch_weight",
        "total_weight",
    ] {
        assert!(split[field].is_null(), "`{field}` is {:?}", split[field]);
    }
    assert_eq!(split["max_utilization"].as_f64(), Some(0.0));
    assert_eq!(split["capacity_exceeded"].as_bool(), Some(false));
}

/// **The mechanical design crosses the wire, and the capacity map reads the header it sized.**
///
/// `calcDesign` is the caller's explicit call, so `mechanical_design` is the flowsheet's flag.
/// The design's nine fields are the record's, not the stream's, and the velocity beside them is
/// read against the NPS-rounded header rather than a stated one.
#[test]
fn the_splitter_publishes_its_mechanical_design() {
    const DESIGNED: &str = r#"
id = "flowsheets.splitter_design"
name = "One splitter with a sized header"

products = ["light", "heavy"]

[[inputs]]
name = "feed_1"
components = ["n-butane", "n-pentane"]
n = 1.0
z = [0.6, 0.4]
P = 1.0e6
T = 300.0

[[instances]]
id = "split1"
unit = "unit_ops.splitter"
[instances.parameters]
split_factors = [0.3, 0.7]
mechanical_design = true
capacity_analysis_enabled = true

[[connections]]
from = "feed_1"
to = "split1.feed"

[[connections]]
from = "split1.products[0]"
to = "light"

[[connections]]
from = "split1.products[1]"
to = "heavy"
"#;
    let report = run_flowsheet(DESIGNED);
    let split = &report.results["split1"];

    assert_eq!(split["design_header_diameter"]["unit"], "m");
    assert_eq!(
        split["design_header_diameter"]["magnitude_si"].as_f64(),
        Some(0.0127)
    );
    assert_eq!(
        split["design_wall_thickness"]["magnitude_si"].as_f64(),
        Some(0.003)
    );
    assert_eq!(
        split["design_header_length"]["magnitude_si"].as_f64(),
        Some(1.5)
    );
    assert_eq!(split["header_weight"]["unit"], "kg");
    assert_eq!(split["design_total_pressure_drop"]["unit"], "Pa");
    assert_eq!(split["manifold_volume"]["unit"], "m**3");
    let velocity = split["velocity_value"]["magnitude_si"]
        .as_f64()
        .expect("a magnitude");
    assert!(
        (velocity - 0.8580239773996624).abs() < 1e-6,
        "the velocity is {velocity}"
    );
    assert_eq!(split["bottleneck_constraint"], "velocity");
}

/// **A result is the last pass's, and it agrees with the stream the same pass published.**
///
/// A result accumulated once and not overwritten would be the *first* pass's, and this flowsheet
/// has a recycle: the mixer is fed the tear, so the first pass runs the heater on a feed the loop
/// has not closed yet. Two numbers from two different passes is the defect, and the only way to
/// see it is to hold the published result against the published stream.
#[test]
fn a_result_is_the_same_passs_as_the_stream_beside_it() {
    let text = std::fs::read_to_string(root().join("specs/flowsheets/demo.toml"))
        .expect("the shipped flowsheet is there");
    let flowsheet: Flowsheet = toml::from_str(&text).expect("it parses");
    let palette = load_palette(&root().join("specs/unit_ops")).expect("the palette loads");
    let session = Session::run(
        flowsheet,
        &palette,
        &BTreeMap::new(),
        ExecutionOrder::Insertion,
    )
    .expect("the flowsheet runs");
    let report: SessionReport = session_report(&session).expect("it writes");
    assert!(report.iterations > 1, "a recycle takes at least two passes");

    let heater = &report.results["hx1"];
    let outlet = &report.streams["hx1.outlet"];
    // **The two shapes differ and both are deliberate.** A *record* writes every scalar as a
    // magnitude and a unit; a model *result* writes the fields its Python dataclass writes, and a
    // flow there is a bare number - which `FIELDS` and the cross-implementation test are what hold
    // the two halves to. So the comparison meets at the magnitude.
    assert_eq!(
        heater["outlet_n"].as_f64(),
        Some(outlet.n.magnitude_si),
        "the heater's flow is not the same pass's stream value"
    );
    for (field, magnitude) in [
        ("outlet_p", outlet.p.magnitude_si),
        ("outlet_t", outlet.temperature.magnitude_si),
    ] {
        assert_eq!(
            heater[field]["magnitude_si"].as_f64(),
            Some(magnitude),
            "the heater's `{field}` is not the same pass's stream value"
        );
    }
}
