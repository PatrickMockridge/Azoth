use azoth_sheets::workspace::{Command, SheetWorkspace};

const DEMO: &str = include_str!("../../../specs/calcbooks/demo.toml");
const BLANK: &str = include_str!("../../../specs/calcbooks/blank.toml");

fn envelope(workspace: &SheetWorkspace) -> serde_json::Value {
    serde_json::from_str(&workspace.to_json().expect("the envelope writes")).expect("json")
}

#[test]
fn the_shipped_books_read_and_round_trip() {
    for document in [BLANK, DEMO] {
        let workspace = SheetWorkspace::open(document).expect("the book reads");
        let written = workspace.document().expect("the book writes");
        let again = SheetWorkspace::open(&written).expect("the written book reads");
        assert_eq!(again.document().expect("writes"), written);
    }
}

#[test]
fn the_demo_is_two_sheets_and_neither_is_ready_until_its_inputs_are_there() {
    let workspace = SheetWorkspace::open(DEMO).expect("the book reads");
    let value = envelope(&workspace);
    assert_eq!(value["book"]["sheets"].as_array().expect("sheets").len(), 2);
    assert_eq!(value["ok"], true);
    assert_eq!(value["dirty"], true, "nothing has been calculated yet");
    for sheet in value["book"]["sheets"].as_array().expect("sheets") {
        assert_eq!(sheet["ready"], true);
        assert_eq!(sheet["stale"], true);
        assert!(
            sheet["diagnostics"]
                .as_array()
                .expect("diagnostics")
                .is_empty()
        );
    }
}

#[test]
fn a_new_sheet_starts_incomplete_and_becomes_ready_when_its_required_inputs_are_set() {
    let mut workspace = SheetWorkspace::open(BLANK).expect("the book reads");
    workspace
        .apply(&Command::AddSheet {
            id: "pump_1".to_string(),
            calculation: "hydraulics.pump_power".to_string(),
            name: None,
        })
        .expect("the sheet is added");
    let value = envelope(&workspace);
    let sheet = &value["book"]["sheets"][0];
    assert_eq!(
        sheet["name"],
        "Pump shaft power from flow, head and efficiency"
    );
    assert_eq!(sheet["ready"], false);
    assert_eq!(
        sheet["missing"],
        serde_json::json!(["rho", "q", "H", "eta"])
    );
    assert_eq!(value["ok"], false);

    for (name, number) in [("rho", 998.0), ("q", 0.01), ("H", 30.0), ("eta", 0.75)] {
        workspace
            .apply(&Command::SetInput {
                sheet: "pump_1".to_string(),
                input: name.to_string(),
                value: serde_json::json!(number),
            })
            .expect("the input is set");
    }
    let value = envelope(&workspace);
    assert_eq!(value["book"]["sheets"][0]["ready"], true);
    assert_eq!(value["ok"], true);
    assert_eq!(value["book"]["sheets"][0]["stale"], true);
}

#[test]
fn a_command_that_names_an_unknown_sheet_or_input_is_refused() {
    let mut workspace = SheetWorkspace::open(DEMO).expect("the book reads");
    assert!(
        workspace
            .apply(&Command::SetInput {
                sheet: "nope".to_string(),
                input: "rho".to_string(),
                value: serde_json::json!(1.0),
            })
            .is_err()
    );
    assert!(
        workspace
            .apply(&Command::SetInput {
                sheet: "pump_power_1".to_string(),
                input: "nope".to_string(),
                value: serde_json::json!(1.0),
            })
            .is_err()
    );
    assert!(
        workspace
            .apply(&Command::SetInput {
                sheet: "pump_power_1".to_string(),
                input: "rho".to_string(),
                value: serde_json::json!("dense"),
            })
            .is_err()
    );
}

#[test]
fn an_unavailable_calculation_cannot_be_added_and_an_unknown_one_is_diagnosed_on_open() {
    let mut workspace = SheetWorkspace::open(BLANK).expect("the book reads");
    assert!(
        workspace
            .apply(&Command::AddSheet {
                id: "crane_1".to_string(),
                calculation: "hydraulics.crane_k_factors".to_string(),
                name: None,
            })
            .is_err()
    );

    let unknown = r#"
id = "calcbook_1"
name = "From the future"

[[sheets]]
id = "future_1"
name = "Future"
calculation = "hydraulics.not_a_calculation"
"#;
    let workspace =
        SheetWorkspace::open(unknown).expect("a book naming an unknown calc still opens");
    let value = envelope(&workspace);
    assert_eq!(value["ok"], false);
    assert_eq!(
        value["book"]["sheets"][0]["diagnostics"][0]["code"],
        "unknown_calculation"
    );
}

#[test]
fn removing_and_renaming_a_sheet_changes_the_document() {
    let mut workspace = SheetWorkspace::open(DEMO).expect("the book reads");
    workspace
        .apply(&Command::RenameSheet {
            id: "pump_power_1".to_string(),
            name: "Pump check".to_string(),
        })
        .expect("rename");
    workspace
        .apply(&Command::RemoveSheet {
            id: "darcy_weisbach_1".to_string(),
        })
        .expect("remove");
    let value = envelope(&workspace);
    assert_eq!(value["book"]["sheets"].as_array().expect("sheets").len(), 1);
    assert_eq!(value["book"]["sheets"][0]["name"], "Pump check");
}

#[test]
fn two_sheets_with_one_id_are_refused_on_open() {
    let twice = r#"
id = "calcbook_1"
name = "Twice"

[[sheets]]
id = "same"
name = "One"
calculation = "hydraulics.pump_power"

[[sheets]]
id = "same"
name = "Two"
calculation = "hydraulics.pump_power"
"#;
    assert!(SheetWorkspace::open(twice).is_err());
}

#[test]
fn calculate_is_refused_for_a_sheet_the_book_does_not_carry() {
    let mut workspace = SheetWorkspace::open(DEMO).expect("the book reads");
    assert!(workspace.calculate(Some("nope")).is_err());
}

#[test]
fn a_hand_written_toml_integer_is_a_quantity() {
    let book = r#"
id = "calcbook_1"
name = "Integers"

[[sheets]]
id = "pump_1"
name = "Pump"
calculation = "hydraulics.pump_power"

[sheets.inputs]
rho = 998
q = 1
H = 30
eta = 1
"#;
    let workspace = SheetWorkspace::open(book).expect("the book reads");
    let value = envelope(&workspace);
    assert_eq!(value["book"]["sheets"][0]["ready"], true);
    assert!(
        value["book"]["sheets"][0]["diagnostics"]
            .as_array()
            .expect("diagnostics")
            .is_empty()
    );
}
