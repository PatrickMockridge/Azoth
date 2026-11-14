//! The instrumentation layer's checker rules: tags, endpoints, ranges, safety declarations and
//! signals.

use std::path::Path;

use azoth_process::{load_palette, parse_flowsheet, validate};

fn palette() -> Vec<azoth_process::UnitOpSpec> {
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("../..");
    load_palette(&root.join("specs/unit_ops")).expect("the palette loads")
}

const BASE: &str = r#"
id = "flowsheets.instrumented"
name = "Instrumented"

products = ["OUT"]

[[inputs]]
name = "FEED"
components = ["methane"]
n = 1.0
z = [1.0]
P = 5.0e5
T = 300.0

[[instances]]
id = "p1"
unit = "unit_ops.pump"
[instances.parameters]
outlet_pressure = 6.0e5
isentropic_efficiency = 0.75

[[connections]]
from = "FEED"
to = "p1.inlet"

[[connections]]
from = "p1.outlet"
to = "OUT"

[[instruments]]
id = "PIC-101"
function = "PIC"
measures = "p1.outlet"
location = "field"
range = { min = 0.0, max = 1.0e7, unit = "Pa" }

[[instruments]]
id = "FV-101"
function = "FV"
actuates = "p1"
location = "field"

[[signals]]
from = "PIC-101"
to = "FV-101"
kind = "pneumatic"
"#;

fn codes(text: &str) -> Vec<&'static str> {
    let flowsheet = parse_flowsheet(text).expect("the fixture parses");
    let palette = palette();
    validate(&flowsheet, &palette)
        .iter()
        .map(|diagnostic| diagnostic.code())
        .collect()
}

#[test]
fn a_well_formed_layer_is_clean() {
    assert!(codes(BASE).is_empty(), "{:?}", codes(BASE));
}

#[test]
fn a_tag_declared_twice_is_refused() {
    let text = BASE.replace("id = \"FV-101\"", "id = \"PIC-101\"");
    assert!(codes(&text).contains(&"duplicate_instrument"));
}

#[test]
fn a_measurement_of_nothing_is_refused() {
    let text = BASE.replace("measures = \"p1.outlet\"", "measures = \"p9.outlet\"");
    assert!(codes(&text).contains(&"unknown_instrument_target"));
}

#[test]
fn a_function_of_lower_case_letters_is_refused() {
    let text = BASE.replace("function = \"PIC\"", "function = \"pic\"");
    assert!(codes(&text).contains(&"bad_instrument_function"));
}

#[test]
fn a_range_whose_minimum_is_not_below_its_maximum_is_refused() {
    let text = BASE.replace("min = 0.0, max = 1.0e7", "min = 2.0e7, max = 1.0e7");
    assert!(codes(&text).contains(&"bad_instrument_range"));
}

#[test]
fn an_unknown_range_unit_is_refused() {
    let text = BASE.replace("unit = \"Pa\" }", "unit = \"furlong\" }");
    assert!(codes(&text).contains(&"bad_instrument_range"));
}

#[test]
fn an_unknown_location_is_refused() {
    let text = BASE.replace("location = \"field\"", "location = \"attic\"");
    assert!(codes(&text).contains(&"bad_instrument_location"));
}

#[test]
fn a_sil_outside_one_to_four_is_refused() {
    let text = format!(
        "{BASE}\n[[instruments]]\nid = \"PSHH-101\"\nfunction = \"PSHH\"\n\
         safety = {{ function = \"SIF-101\", sil = 5 }}\n"
    );
    assert!(codes(&text).contains(&"bad_safety_function"));
}

#[test]
fn a_voting_architecture_that_is_not_n_oo_m_is_refused() {
    let text = format!(
        "{BASE}\n[[instruments]]\nid = \"PSHH-101\"\nfunction = \"PSHH\"\n\
         safety = {{ function = \"SIF-101\", sil = 2, voting = \"two out of three\" }}\n"
    );
    assert!(codes(&text).contains(&"bad_safety_function"));
}

#[test]
fn a_signal_to_nowhere_is_refused() {
    let text = BASE.replace("to = \"FV-101\"", "to = \"FV-999\"");
    assert!(codes(&text).contains(&"unknown_signal_end"));
}

#[test]
fn a_signal_declared_twice_is_refused() {
    let text =
        format!("{BASE}\n[[signals]]\nfrom = \"PIC-101\"\nto = \"FV-101\"\nkind = \"electric\"\n");
    assert!(codes(&text).contains(&"duplicate_signal"));
}

#[test]
fn an_unknown_signal_carrier_is_refused() {
    let text = BASE.replace("kind = \"pneumatic\"", "kind = \"smoke signal\"");
    assert!(codes(&text).contains(&"bad_signal_kind"));
}
