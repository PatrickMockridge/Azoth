//! The command model: one typed edit at a time, and the checker is the only rule set.
//!
//! `docs/src/architecture/middleware.md` states it: every edit is a command that re-runs the
//! checker, which is what makes a red arrow a structured `OverfedPort` rather than a string. What
//! that costs is that a command must be exactly as strict as the document it edits - a typo'd
//! field read past in silence would be the same defect as a flowsheet key read past in silence,
//! one layer up.

use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

use azoth_process::middleware::command::{Command, apply};
use azoth_process::{Flowsheet, UnitOpSpec, load_palette, validate};

fn root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../..")
}

fn palette() -> Vec<UnitOpSpec> {
    load_palette(&root().join("specs/unit_ops")).expect("the shipped palette loads")
}

fn demo() -> Flowsheet {
    let text = std::fs::read_to_string(root().join("specs/flowsheets/demo.toml"))
        .expect("the shipped flowsheet is there");
    Flowsheet::from_toml(&text).expect("it parses")
}

/// One of every variant, so a new one cannot arrive without being carried here.
fn every_command() -> Vec<Command> {
    vec![
        Command::AddInstance {
            id: "p2".into(),
            unit: "unit_ops.pump".into(),
            parameters: BTreeMap::new(),
        },
        Command::RemoveInstance { id: "hx1".into() },
        Command::Connect {
            from: "p1.outlet".into(),
            to: "sep1.feed".into(),
        },
        Command::Disconnect {
            from: "p1.outlet".into(),
            to: "hx1.inlet".into(),
        },
        Command::SetParameter {
            instance: "p1".into(),
            name: "outlet_pressure".into(),
            value: serde_json::json!(3.0e6),
        },
        Command::UnsetParameter {
            instance: "p1".into(),
            name: "outlet_pressure".into(),
        },
        Command::AddInput {
            name: "feed_2".into(),
            components: vec!["methane".into()],
            n: 1.0,
            z: vec![1.0],
            p: 5.0e5,
            t: 300.0,
        },
        Command::RemoveInput {
            name: "feed_1".into(),
        },
        Command::AddProduct {
            name: "purge".into(),
        },
        Command::RemoveProduct {
            name: "vapour_product".into(),
        },
        Command::AddRecycle {
            stream: "r2".into(),
            from: "sep1.liquid".into(),
            to: "mix1.feed".into(),
        },
        Command::RemoveRecycle {
            stream: "recycle_1".into(),
        },
        Command::SetRecycle {
            stream: "recycle_1".into(),
            field: "flow_tolerance".into(),
            value: serde_json::json!(1e-6),
        },
        Command::UnsetRecycleField {
            stream: "recycle_1".into(),
            field: "max_iterations".into(),
        },
        Command::SetPosition {
            node: "instance:sep1".into(),
            x: 1234.0,
            y: 56.0,
        },
    ]
}

/// **The measurement this module's strictness rests on.** `deny_unknown_fields` on an internally
/// tagged enum's variants: the tag is consumed by the enum's content deserializer, and the
/// variant must still refuse a key it does not declare. Serde does not promise this combination,
/// so it is measured rather than assumed - and if it were to stop holding, the fallback is an
/// adjacently tagged `{"command": …, "data": {…}}` with the same deny on `data`.
#[test]
fn a_command_carrying_an_unknown_field_is_refused() {
    let refused = [
        r#"{"command": "remove_instance", "id": "hx1", "nope": 1}"#,
        r#"{"command": "connect", "from": "a", "to": "b", "nope": 1}"#,
        r#"{"command": "set_parameter", "instance": "p1", "name": "x", "value": 1, "nope": 1}"#,
        r#"{"command": "set_position", "node": "instance:p1", "x": 0, "y": 0, "nope": 1}"#,
    ];
    for text in refused {
        let error = serde_json::from_str::<Command>(text)
            .err()
            .unwrap_or_else(|| {
                panic!("{text} was accepted with a key the command does not declare")
            });
        assert!(
            error.to_string().contains("nope"),
            "the refusal does not name the key: {error}"
        );
    }
    // And the tag itself is not an unknown field.
    assert!(
        serde_json::from_str::<Command>(r#"{"command": "remove_instance", "id": "hx1"}"#).is_ok()
    );
}

#[test]
fn every_command_round_trips_through_json() {
    let commands = every_command();
    assert_eq!(commands.len(), 15, "the enum and this list have drifted");
    for command in commands {
        let text = serde_json::to_string(&command).expect("a command writes");
        let read: Command = serde_json::from_str(&text).expect("and reads back");
        assert_eq!(read, command, "{text}");
    }
}

/// **A command naming something that does not exist is a no-op plus the checker's own word for
/// it.** A refusal here would be a second rule set: `UnknownInstance` is already the answer.
#[test]
fn a_command_that_names_nothing_is_a_no_op_and_the_checker_says_so() {
    let palette = palette();
    let before = demo();
    for command in [
        Command::RemoveInstance { id: "nope".into() },
        Command::UnsetParameter {
            instance: "nope".into(),
            name: "outlet_pressure".into(),
        },
        Command::RemoveRecycle {
            stream: "nope".into(),
        },
    ] {
        let mut flowsheet = before.clone();
        apply(&mut flowsheet, &palette, &command).expect("the command takes effect");
        assert_eq!(
            flowsheet, before,
            "{command:?} changed a document it named nothing in"
        );
    }

    // The diagnostic is the checker's, and it is already there for the same document.
    let mut flowsheet = before.clone();
    apply(
        &mut flowsheet,
        &palette,
        &Command::Connect {
            from: "nowhere".into(),
            to: "mix1.feed".into(),
        },
    )
    .expect("the connection is written");
    let diagnostics = validate(&flowsheet, &palette);
    assert!(
        diagnostics
            .iter()
            .any(|diagnostic| diagnostic.code() == "unknown_feed"),
        "the checker names what the command could not: {diagnostics:?}"
    );
}

/// **Every edit lands, and what it produced reads back as itself.** This is the editor's round
/// trip executed rather than claimed: the document a command produced is written and reparsed on
/// every command, so a change that lost something cannot reach a caller.
#[test]
fn each_edit_lands_and_the_document_round_trips() {
    let palette = palette();
    for command in every_command() {
        // **`unset_recycle_field` is the one command the shipped document cannot show**: it states
        // none of the tear's seven settings - that is what a silence *is* - so there is nothing for
        // an unset to remove. The test below does it on a document that states one.
        if matches!(command, Command::UnsetRecycleField { .. }) {
            continue;
        }
        let mut flowsheet = demo();
        apply(&mut flowsheet, &palette, &command).expect("the command takes effect");
        assert_ne!(flowsheet, demo(), "{command:?} changed nothing");

        let written = flowsheet.to_toml().expect("it writes");
        let read = Flowsheet::from_toml(&written).expect("it reads back");
        assert_eq!(read, flowsheet, "{command:?} lost something:\n{written}");
        assert_eq!(
            read.to_toml().expect("it writes again"),
            written,
            "{command:?}: writing is not a fixed point"
        );
    }
}

/// **`set` and `unset` are each other's inverse**, and the difference between an unstated setting
/// and a stated one is the reason both exist: `[[recycles]]` writes only what it states, so a form
/// that showed the class's defaults could not tell a declaration that *meant* this number from one
/// that never mentioned it.
#[test]
fn a_recycle_field_can_be_returned_to_silence() {
    let palette = palette();
    let mut flowsheet = demo();
    let field = |name: &str| Command::UnsetRecycleField {
        stream: "recycle_1".into(),
        field: name.into(),
    };

    // Unsetting what is already unstated changes nothing, which is the shell of the same fact.
    let untouched = flowsheet.clone();
    apply(&mut flowsheet, &palette, &field("max_iterations")).expect("it is a no-op");
    assert_eq!(flowsheet, untouched);

    // Stated, written, then unset - and it leaves the document rather than becoming a default.
    flowsheet.recycles[0].max_iterations = Some(25);
    let written = flowsheet.to_toml().expect("it writes");
    assert!(written.contains("max_iterations = 25"), "{written}");
    apply(&mut flowsheet, &palette, &field("max_iterations")).expect("it unsets");
    assert_eq!(flowsheet.recycles[0].max_iterations, None);
    let written = flowsheet.to_toml().expect("it writes");
    assert!(
        !written.contains("max_iterations"),
        "an unset setting is absent from the document, not written as the default: {written}"
    );

    // Every one of the seven can be returned, and the name is checked as `set` checks it.
    for name in azoth_process::middleware::command::RECYCLE_FIELDS {
        apply(&mut flowsheet, &palette, &field(name)).expect("every setting unsets");
    }
    let error = apply(&mut flowsheet, &palette, &field("nope"))
        .expect_err("a field the schema has no word for is refused");
    assert!(error.to_string().contains("nope"), "{error}");
}

#[test]
fn a_removed_instance_takes_its_edges_with_it() {
    let palette = palette();
    let mut flowsheet = demo();
    apply(
        &mut flowsheet,
        &palette,
        &Command::RemoveInstance { id: "hx1".into() },
    )
    .expect("the command takes effect");

    // `hx1` had two connections, and neither is left to name an instance that is gone.
    assert!(
        !flowsheet
            .connections
            .iter()
            .any(|connection| connection.from.starts_with("hx1.")
                || connection.to.starts_with("hx1.")),
        "{:?}",
        flowsheet.connections
    );
    assert_eq!(flowsheet.connections.len(), 3);
    let diagnostics = validate(&flowsheet, &palette);
    assert!(
        !diagnostics
            .iter()
            .any(|diagnostic| diagnostic.code() == "unknown_instance"),
        "removing an instance left dangling endpoints: {diagnostics:?}"
    );
}

/// The kind a parameter's value is written as, which is the one thing the document does not know.
#[test]
fn a_parameter_is_written_as_the_kind_its_declaration_states() {
    let palette = palette();
    let mut flowsheet = demo();

    // A vector parameter takes a list, and a front-end sending a number is refused rather than
    // writing a document the checker would refuse a moment later.
    let vector = Command::SetParameter {
        instance: "sep1".into(),
        name: "gas_in_liquid".into(),
        value: serde_json::json!(0.0),
    };
    apply(&mut flowsheet, &palette, &vector).expect("a quantity takes a number");
    assert_eq!(
        flowsheet.instances[3].parameters.get("gas_in_liquid"),
        Some(&toml::Value::Float(0.0))
    );

    let mut flowsheet = demo();
    flowsheet.instances.push(azoth_process::Instance {
        id: "s1".into(),
        unit: "unit_ops.splitter".into(),
        parameters: BTreeMap::new(),
    });
    let as_list = Command::SetParameter {
        instance: "s1".into(),
        name: "split_factors".into(),
        value: serde_json::json!([0.5, 0.5]),
    };
    apply(&mut flowsheet, &palette, &as_list).expect("a vector takes a list");
    assert_eq!(
        flowsheet.instances[4].parameters.get("split_factors"),
        Some(&toml::Value::Array(vec![
            toml::Value::Float(0.5),
            toml::Value::Float(0.5)
        ]))
    );

    let mut flowsheet = demo();
    flowsheet.instances.push(azoth_process::Instance {
        id: "s1".into(),
        unit: "unit_ops.splitter".into(),
        parameters: BTreeMap::new(),
    });
    let as_number = Command::SetParameter {
        instance: "s1".into(),
        name: "split_factors".into(),
        value: serde_json::json!(0.5),
    };
    let error = apply(&mut flowsheet, &palette, &as_number)
        .expect_err("a vector parameter does not take a number");
    assert!(error.to_string().contains("list"), "{error}");
}

/// TOML has no null, so a parameter cannot be set to one - and the refusal says that rather than
/// writing a zero.
#[test]
fn a_null_parameter_is_refused() {
    let palette = palette();
    let mut flowsheet = demo();
    let error = apply(
        &mut flowsheet,
        &palette,
        &Command::SetParameter {
            instance: "p1".into(),
            name: "outlet_pressure".into(),
            value: serde_json::Value::Null,
        },
    )
    .expect_err("a null is refused");
    assert!(error.to_string().contains("null"), "{error}");
    assert_eq!(flowsheet, demo(), "the document is untouched");
}

#[test]
fn a_position_is_only_placed_on_a_node_the_document_declares() {
    let palette = palette();
    let mut flowsheet = demo();
    apply(
        &mut flowsheet,
        &palette,
        &Command::SetPosition {
            node: "instance:sep1".into(),
            x: 1234.0,
            y: 56.0,
        },
    )
    .expect("a declared node is placed");
    assert_eq!(
        flowsheet
            .layout
            .as_ref()
            .and_then(|layout| layout.position(azoth_process::NodeRole::Instance, "sep1")),
        Some([1234.0, 56.0])
    );

    // A node the document does not declare cannot be placed: a remove of something absent is
    // already what was asked for, and a placement of something absent cannot take effect.
    let error = apply(
        &mut flowsheet,
        &palette,
        &Command::SetPosition {
            node: "instance:nope".into(),
            x: 0.0,
            y: 0.0,
        },
    )
    .expect_err("an undeclared node is refused");
    assert!(error.to_string().contains("nope"), "{error}");

    // And a node id whose role is not one of the three is a malformed command rather than a
    // missing node.
    let error = apply(
        &mut flowsheet,
        &palette,
        &Command::SetPosition {
            node: "widget:sep1".into(),
            x: 0.0,
            y: 0.0,
        },
    )
    .expect_err("an unknown role is refused");
    assert!(error.to_string().contains("node id"), "{error}");
}

#[test]
fn a_recycle_parameter_is_five_numbers_a_count_and_a_name() {
    let palette = palette();
    let mut flowsheet = demo();
    let field = |name: &str, value: serde_json::Value| Command::SetRecycle {
        stream: "recycle_1".into(),
        field: name.into(),
        value,
    };

    apply(
        &mut flowsheet,
        &palette,
        &field("flow_tolerance", serde_json::json!(1e-6)),
    )
    .expect("a tolerance is a number");
    apply(
        &mut flowsheet,
        &palette,
        &field("max_iterations", serde_json::json!(25)),
    )
    .expect("a count is a whole number");
    apply(
        &mut flowsheet,
        &palette,
        &field("acceleration_method", serde_json::json!("wegstein")),
    )
    .expect("the method is a name");

    let recycle = &flowsheet.recycles[0];
    assert_eq!(recycle.flow_tolerance, Some(1e-6));
    assert_eq!(recycle.max_iterations, Some(25));
    assert_eq!(recycle.acceleration_method.as_deref(), Some("wegstein"));

    // Measured: the class's own reason, not a second one.
    apply(
        &mut flowsheet,
        &palette,
        &field("acceleration_method", serde_json::json!("broyden")),
    )
    .expect("the name is written; the checker is what refuses it");
    let diagnostics = validate(&flowsheet, &palette);
    let refusal = diagnostics
        .iter()
        .find(|diagnostic| diagnostic.code() == "acceleration")
        .unwrap_or_else(|| panic!("a broyden tear is not reported: {diagnostics:?}"));
    assert_eq!(refusal.severity(), azoth_process::Severity::Error);
    assert_eq!(refusal.location().section, "recycles");
    assert_eq!(refusal.location().path, "recycle_1");
    assert!(
        refusal.message().contains("broyden") && refusal.message().contains("Measured"),
        "the class's own measured reason, verbatim: {}",
        refusal.message()
    );

    // And a name that is none of the three is the same refusal, not a silent default.
    let mut flowsheet = demo();
    apply(
        &mut flowsheet,
        &palette,
        &field("acceleration_method", serde_json::json!("wegstien")),
    )
    .expect("the name is written");
    let diagnostics = validate(&flowsheet, &palette);
    assert!(
        diagnostics
            .iter()
            .any(|diagnostic| diagnostic.code() == "acceleration"
                && diagnostic.message().contains("wegstien")),
        "{diagnostics:?}"
    );

    // A field the schema has no word for is a malformed command, not a document fact.
    let error = apply(
        &mut flowsheet,
        &palette,
        &field("nope", serde_json::json!(1)),
    )
    .expect_err("an unknown field is refused");
    assert!(error.to_string().contains("nope"), "{error}");
    // And a number where a name belongs.
    let error = apply(
        &mut flowsheet,
        &palette,
        &field("max_iterations", serde_json::json!("many")),
    )
    .expect_err("a name where a count belongs is refused");
    assert!(error.to_string().contains("number"), "{error}");
}

/// **The names a tool offers are the class's own**, which is the half of "not a second surface"
/// that a schema can quietly break: an `enum` in a schema is a list of strings, and a list of
/// strings is a copy the moment nobody holds it to the thing it names.
#[test]
fn the_acceleration_names_a_tool_offers_are_the_recycle_class_s_own() {
    let tools = azoth_process::middleware::tools::tools(&palette());
    let set_recycle = tools
        .iter()
        .find(|tool| tool.name == "set_recycle")
        .expect("set_recycle is a tool");
    assert_eq!(
        set_recycle.input_schema["properties"]["value"]["oneOf"][1]["enum"],
        serde_json::json!(azoth_process::recycle::ACCELERATION_NAMES),
    );
    // And the numeric arm is still there, so a form can tell the two apart rather than treating
    // every setting as a name.
    assert_eq!(
        set_recycle.input_schema["properties"]["value"]["oneOf"][0]["type"],
        "number"
    );
}

/// **The seven settings a tool offers are the command model's own list.** The same gate the
/// acceleration names have, for the same reason: an `enum` in a schema is a copy of a list the
/// moment nothing holds it to the list.
#[test]
fn the_recycle_fields_a_tool_offers_are_the_command_model_s_own() {
    let tools = azoth_process::middleware::tools::tools(&palette());
    for name in ["set_recycle", "unset_recycle_field"] {
        let tool = tools
            .iter()
            .find(|tool| tool.name == name)
            .unwrap_or_else(|| panic!("{name} is a tool"));
        assert_eq!(
            tool.input_schema["properties"]["field"]["enum"],
            serde_json::json!(azoth_process::middleware::command::RECYCLE_FIELDS),
            "{name}'s field enum"
        );
    }
    // And the inverse exists, which is what a control that offers "the default" needs.
    assert!(tools.iter().any(|tool| tool.name == "unset_recycle_field"));
}

/// **The agent's tools are the command model, not a surface beside it.** The names are read off
/// each variant's `Debug` form, so a command cannot arrive without a tool and a tool cannot
/// arrive without a command - and the schema's `command` tag is held to the tool's own name,
/// because a schema that spelled it differently would describe a call the deserialiser refuses
/// with no clue why.
#[test]
fn every_command_has_a_tool_and_no_tool_is_a_second_surface() {
    let tools = azoth_process::middleware::tools::tools(&palette());
    let commands = every_command();
    assert_eq!(
        tools.len(),
        commands.len(),
        "a tool per command, and no more"
    );

    let variants: Vec<String> = commands
        .iter()
        .map(|command| {
            let debug = format!("{command:?}");
            let end = debug.find(['{', '(', ' ']).unwrap_or(debug.len());
            let mut snake = String::new();
            for (index, character) in debug[..end].char_indices() {
                if character.is_uppercase() && index > 0 {
                    snake.push('_');
                }
                snake.extend(character.to_lowercase());
            }
            snake
        })
        .collect();

    assert_eq!(
        tools.iter().map(|tool| tool.name).collect::<Vec<_>>(),
        variants
    );

    for tool in &tools {
        let schema = &tool.input_schema;
        assert_eq!(schema["type"], "object");
        assert_eq!(schema["additionalProperties"], false, "{}", tool.name);
        assert_eq!(
            schema["properties"]["command"]["const"], tool.name,
            "{}'s tag is not its own name",
            tool.name
        );
        assert!(
            schema["required"]
                .as_array()
                .expect("an array")
                .contains(&serde_json::json!("command")),
            "{} does not require its tag",
            tool.name
        );
        assert!(
            !tool.description.is_empty(),
            "{} says nothing about itself",
            tool.name
        );
    }

    // **The schema describes calls the deserialiser accepts.** A tool object with its `command`
    // tag and every required property beside it reads as that command - which is the property
    // that makes the schema usable rather than merely present.
    for (tool, command) in tools.iter().zip(&commands) {
        let mut call = serde_json::Map::new();
        for name in schema_root_properties(tool) {
            call.insert(name.clone(), value_for(command, &name));
        }
        call.insert("command".to_string(), serde_json::json!(tool.name));
        let text = serde_json::Value::Object(call).to_string();
        let read: Command = serde_json::from_str(&text)
            .unwrap_or_else(|error| panic!("{}: {error}\n{text}", tool.name));
        assert_eq!(&read, command, "{}", tool.name);
    }
}

/// The properties a tool's schema requires, `command` included.
fn schema_root_properties(tool: &azoth_process::middleware::tools::Tool) -> Vec<String> {
    tool.input_schema["required"]
        .as_array()
        .expect("an array")
        .iter()
        .filter_map(serde_json::Value::as_str)
        .map(str::to_string)
        .collect()
}

/// A value for one property, read off the command the tool describes.
fn value_for(command: &Command, property: &str) -> serde_json::Value {
    let document = serde_json::to_value(command).expect("a command writes");
    document
        .get(property)
        .cloned()
        .unwrap_or_else(|| serde_json::json!("any"))
}
