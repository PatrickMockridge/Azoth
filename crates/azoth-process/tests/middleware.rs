//! The interoperation surface: what a front-end reads.
//!
//! `docs/src/architecture/middleware.md` states the layer. Its first promise is that a
//! diagnostic reaches a widget as a structured record rather than as a line of text, and the two
//! things that can quietly go wrong with that are a code that stops matching its variant and a
//! target that points at nothing. Both are held here rather than by review.

use azoth_process::middleware::diagnostic::{DiagnosticRecord, records};
use azoth_process::middleware::form::{FormSpec, dimension_id, forms};
use azoth_process::{
    Diagnostic, Flowsheet, NodeRole, Severity, Target, UnitOpSpec, load_palette, split_node_id,
    validate,
};

/// One of every variant.
///
/// **The list is the point.** `Diagnostic::code`'s match is exhaustive, so a new variant is a
/// compile error until it has a code; this list is what catches the other half - a code that was
/// written but does not agree with the variant it was written for. A variant missing from here
/// fails `every_variant_is_in_this_list` rather than passing silently.
fn every_variant() -> Vec<Diagnostic> {
    vec![
        Diagnostic::UnknownDimension {
            unit_op: "unit_ops.pump".into(),
            port: "inlet".into(),
            field: "n".into(),
            dimension: "molar_flow".into(),
        },
        Diagnostic::UnknownParameterUnit {
            unit_op: "unit_ops.pump".into(),
            parameter: "outlet_pressure".into(),
            unit: "psi".into(),
        },
        Diagnostic::DuplicateUnitOpId {
            id: "unit_ops.pump".into(),
        },
        Diagnostic::DuplicatePort {
            unit_op: "unit_ops.pump".into(),
            port: "inlet".into(),
        },
        Diagnostic::DuplicateInstance { id: "p1".into() },
        Diagnostic::DuplicateFeed {
            name: "feed_1".into(),
        },
        Diagnostic::DuplicateProduct {
            name: "vapour_product".into(),
        },
        Diagnostic::NameCollision {
            name: "sep1".into(),
            first: "an instance",
            second: "a feed",
        },
        Diagnostic::UnknownUnitOp {
            instance: "p1".into(),
            unit: "unit_ops.pimp".into(),
        },
        Diagnostic::UnknownParameter {
            instance: "p1".into(),
            parameter: "outlet_presure".into(),
        },
        Diagnostic::MissingParameter {
            instance: "sep1".into(),
            parameter: "gas_in_liquid".into(),
        },
        Diagnostic::UnknownFeed {
            name: "feed_2".into(),
        },
        Diagnostic::UnknownProduct {
            name: "purge".into(),
        },
        Diagnostic::UnknownInstance { name: "p2".into() },
        Diagnostic::UnknownPort {
            instance: "p1".into(),
            port: "discharge".into(),
        },
        Diagnostic::ProducerNotOutlet {
            instance: "p1".into(),
            port: "inlet".into(),
        },
        Diagnostic::ConsumerNotInlet {
            instance: "p1".into(),
            port: "outlet".into(),
        },
        Diagnostic::TypeMismatch {
            from: "p1.outlet".into(),
            to: "hx1.inlet".into(),
            detail: "`h` is a molar_energy and `n` is a molar_flow".into(),
        },
        Diagnostic::OverfedPort {
            instance: "sep1".into(),
            port: "feed".into(),
            count: 2,
        },
        Diagnostic::UnderfedPort {
            instance: "mix1".into(),
            port: "feed".into(),
            count: 0,
        },
        Diagnostic::UnusedFeed {
            name: "feed_1".into(),
        },
        Diagnostic::UnusedProduct {
            name: "vapour_product".into(),
        },
        Diagnostic::InputRecord {
            name: "feed_1".into(),
            detail: "the input `feed_1` names no substance, so it is not a fluid".into(),
        },
        Diagnostic::EndpointIndex {
            endpoint: "split1.products[2]".into(),
            detail: "`products` takes no position at 2".into(),
        },
        Diagnostic::ParameterKind {
            instance: "p1".into(),
            parameter: "outlet_pressure".into(),
            detail: "`outlet_pressure` is a quantity and `\"high\"` is text".into(),
        },
        Diagnostic::UnrecycledLoop {
            detail: "the loop `p1.outlet -> ... -> p1.inlet` declares no recycle".into(),
        },
    ]
}

/// The variant's own name, read off its `Debug` form.
fn variant_name(diagnostic: &Diagnostic) -> String {
    let debug = format!("{diagnostic:?}");
    let end = debug.find(['{', '(', ' ']).unwrap_or(debug.len());
    debug[..end].to_string()
}

/// `UnknownDimension` -> `unknown_dimension`.
fn snake_case(name: &str) -> String {
    let mut out = String::with_capacity(name.len() + 4);
    for (index, character) in name.char_indices() {
        if character.is_uppercase() && index > 0 {
            out.push('_');
        }
        out.extend(character.to_lowercase());
    }
    out
}

#[test]
fn every_diagnostic_code_is_its_variant_name_in_snake_case() {
    for diagnostic in every_variant() {
        assert_eq!(
            diagnostic.code(),
            snake_case(&variant_name(&diagnostic)),
            "the code and the variant disagree about {diagnostic:?}"
        );
    }
}

#[test]
fn every_variant_is_in_this_list() {
    // The count is written down so that a variant added to the enum and forgotten here fails
    // rather than passing a shorter list.
    let diagnostics = every_variant();
    assert_eq!(diagnostics.len(), 26, "the enum and this list have drifted");
    let mut codes: Vec<&str> = diagnostics.iter().map(Diagnostic::code).collect();
    codes.sort_unstable();
    let before = codes.len();
    codes.dedup();
    assert_eq!(before, codes.len(), "two variants share a code");
}

#[test]
fn no_diagnostic_targets_nothing() {
    for diagnostic in every_variant() {
        let subject = match diagnostic.target() {
            Target::Node { id, .. } => id,
            Target::Handle { node, port, .. } => format!("{node}.{port}"),
            Target::Parameter { node, name } => format!("{node}.{name}"),
            Target::Edge { from, to } => format!("{from}->{to}"),
            Target::Endpoint { endpoint } => endpoint,
            Target::Palette { id, .. } => id,
            // The document whole is the one subject that is not a name.
            Target::Document => continue,
        };
        assert!(!subject.is_empty(), "{diagnostic:?} targets an empty name");
    }
}

#[test]
fn a_port_target_and_its_location_name_the_same_port() {
    // `location` is the text and `target` is what a widget draws; where a diagnostic is *located*
    // at a port the two must name the same one, or a mark lands on the wrong handle.
    //
    // **`EndpointIndex` is the case where they differ, and rightly.** It is located in
    // `connections`, because that is where the bad endpoint was written, and points at the port
    // on the instance, because that is what the position belongs to - a different question, and
    // the section is what tells them apart.
    let mut checked = 0;
    for diagnostic in every_variant() {
        let Target::Handle { node, port, .. } = diagnostic.target() else {
            continue;
        };
        let location = diagnostic.location();
        if location.section != "instances" {
            continue;
        }
        assert_eq!(
            format!("{node}.ports.{port}"),
            location.path,
            "{diagnostic:?} points at a different port than it is located at"
        );
        checked += 1;
    }
    assert_eq!(checked, 5, "the port variants are what this test is about");
}

#[test]
fn an_endpoint_index_points_at_the_position_it_names() {
    let diagnostic = Diagnostic::EndpointIndex {
        endpoint: "split1.products[2]".into(),
        detail: String::new(),
    };
    assert_eq!(
        diagnostic.target(),
        Target::Handle {
            role: NodeRole::Instance,
            node: "split1".into(),
            port: "products".into(),
            index: Some(2),
        }
    );
}

#[test]
fn a_node_id_round_trips_through_its_role_prefix() {
    for role in [NodeRole::Instance, NodeRole::Input, NodeRole::Product] {
        let id = role.node_id("sep1");
        assert_eq!(id, format!("{}:sep1", role.name()));
        assert_eq!(split_node_id(&id), Some((role, "sep1")));
    }
    assert_eq!(split_node_id("sep1"), None, "a bare name is not a node id");
    assert_eq!(
        split_node_id("widget:sep1"),
        None,
        "an unknown role is not one"
    );
}

#[test]
fn a_record_is_the_documented_shape() {
    let record = DiagnosticRecord::from(&Diagnostic::OverfedPort {
        instance: "m1".into(),
        port: "feed".into(),
        count: 2,
    });
    let document = serde_json::to_value(&record).expect("a record serialises");
    assert_eq!(
        document,
        serde_json::json!({
            "code": "overfed_port",
            "severity": "error",
            "section": "instances",
            "path": "m1.ports.feed",
            "target": { "kind": "handle", "role": "instance", "node": "m1",
                        "port": "feed", "index": null },
            "message": "`feed` takes one stream and 2 are connected",
            "detail": { "count": 2 },
        })
    );
}

#[test]
fn a_variant_with_no_extra_fields_still_carries_a_detail() {
    let record = DiagnosticRecord::from(&Diagnostic::DuplicateInstance { id: "p1".into() });
    assert_eq!(record.detail, serde_json::json!({}));
    assert_eq!(
        record.target,
        Target::Node {
            role: NodeRole::Instance,
            id: "p1".into(),
        }
    );
}

#[test]
fn a_warning_is_a_warning_on_the_wire() {
    let all = every_variant();
    let document = serde_json::to_value(records(&all)).expect("records serialise");
    let entries = document.as_array().expect("an array");
    let warnings: Vec<&str> = entries
        .iter()
        .filter(|entry| entry["severity"] == "warning")
        .map(|entry| entry["code"].as_str().unwrap_or_default())
        .collect();
    assert_eq!(warnings, ["unused_feed", "unused_product"]);
}

/// The repository root, which is where `specs/` lives.
fn root() -> std::path::PathBuf {
    std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../..")
}

fn palette() -> Vec<UnitOpSpec> {
    load_palette(&root().join("specs/unit_ops")).expect("the shipped palette loads")
}

#[test]
fn every_palette_entry_has_a_form_and_every_form_a_palette_entry() {
    let palette = palette();
    let forms = forms(&palette);
    assert_eq!(forms.len(), palette.len(), "a form per entry, and no more");
    for (spec, form) in palette.iter().zip(&forms) {
        assert_eq!(spec.id, form.id, "the form is the entry's");
        assert_eq!(spec.name, form.name);
        assert_eq!(spec.parameters.len(), form.parameters.len());
        assert_eq!(spec.ports.len(), form.ports.len());
        assert!(form.source.is_some(), "{} has no provenance", form.id);
    }
}

/// **The count of entries a form can render, three ways.**
///
/// Twenty-nine palette entries, twenty-seven with a model, twenty-six with a kernel - and the
/// three that are not runnable are in three different states: one has a model and its ranges, one
/// has parameters and no model, one declares no parameter at all. Collapsing them into one
/// "unsupported" row would lose which of those a caller is looking at.
#[test]
fn the_three_unrunnable_entries_are_three_different_states() {
    let palette = palette();
    let forms = forms(&palette);
    assert_eq!(forms.len(), 29);
    assert_eq!(forms.iter().filter(|f| f.model.is_some()).count(), 27);
    assert_eq!(forms.iter().filter(|f| f.runnable).count(), 26);

    let refused: Vec<&FormSpec> = forms.iter().filter(|f| !f.runnable).collect();
    assert_eq!(refused.len(), 3);
    for form in &refused {
        assert!(
            form.refusal.is_some(),
            "{} is not runnable and says nothing about why",
            form.id
        );
    }

    let by_id = |id: &str| forms.iter().find(|f| f.id == id).expect("declared");
    // A model, its ranges, and a kernel that refuses it: the declaration describes the packing
    // and not the column.
    let packed = by_id("unit_ops.packed_column");
    assert!(packed.model.is_some());
    assert!(!packed.parameters.is_empty());
    // Parameters and no model: nothing declares what kind of thing they are.
    let rate_based = by_id("unit_ops.rate_based_packed_column");
    assert!(rate_based.model.is_none());
    assert!(!rate_based.parameters.is_empty());
    // And one with no parameters at all, only ports.
    let absorber = by_id("unit_ops.simple_absorber");
    assert!(absorber.model.is_none());
    assert!(absorber.parameters.is_empty());
    assert_eq!(absorber.ports.len(), 4);
}

#[test]
fn every_parameter_of_a_modelled_entry_has_a_kind_and_the_other_two_do_not() {
    let mut unknown = Vec::new();
    for form in forms(&palette()) {
        for parameter in &form.parameters {
            if parameter.kind == "unknown" {
                unknown.push(format!("{}.{}", form.id, parameter.name));
            }
        }
    }
    // Measured: the kinds are the five `executor::dispatch` has a reader for, so the kind rule
    // covers every parameter of every entry that has a model - and the two entries without one
    // are exactly the two that say `unknown`.
    assert_eq!(
        unknown,
        [
            "unit_ops.rate_based_packed_column.column_diameter",
            "unit_ops.rate_based_packed_column.column_solver",
            "unit_ops.rate_based_packed_column.convergence_tolerance",
            "unit_ops.rate_based_packed_column.film_model",
            "unit_ops.rate_based_packed_column.heat_transfer_model",
            "unit_ops.rate_based_packed_column.mass_transfer_correlation",
            "unit_ops.rate_based_packed_column.max_iterations",
            "unit_ops.rate_based_packed_column.number_of_segments",
            "unit_ops.rate_based_packed_column.packed_height",
            "unit_ops.rate_based_packed_column.packing_type",
            "unit_ops.rate_based_packed_column.segment_solver",
        ]
    );

    // Every other kind is one the shim reads, which is what makes the rule total.
    let readable = ["quantity", "vector", "boolean", "enum", "string"];
    for form in forms(&palette()) {
        for parameter in &form.parameters {
            if parameter.kind == "unknown" {
                continue;
            }
            assert!(
                readable.contains(&parameter.kind),
                "{}.{} is a {} and no reader takes one",
                form.id,
                parameter.name,
                parameter.kind
            );
        }
    }
}

/// **A model's bounds, held to the palette's names.**
///
/// A model bounds its stream fields as well as its parameters - the compressor's `inlet_t` must
/// be a temperature - and a form has no field for those. They are named rather than dropped,
/// and this is the gate that says the list is exactly the difference and not a sample of it.
#[test]
fn every_model_bound_lands_on_a_parameter_or_in_the_unmodelled_list() {
    let mut checked = 0;
    for form in forms(&palette()) {
        let Some(model) = form.model else {
            assert!(form.unmodelled_ranges.is_empty(), "{}", form.id);
            continue;
        };
        let spec = azoth_process::model_gen::model(model)
            .unwrap_or_else(|| panic!("{} names {model}, which is not a model", form.id));
        let parameters: Vec<&str> = form
            .parameters
            .iter()
            .map(|parameter| parameter.name.as_str())
            .collect();
        // Every bound the model states on an input is either a parameter's range or named as
        // unmodelled - the two together are the model's whole input-bounds list.
        for check in spec.input_checks() {
            let on_a_parameter = parameters.contains(&check.quantity);
            let named = form
                .unmodelled_ranges
                .iter()
                .any(|range| range == check.quantity);
            assert!(
                on_a_parameter || named,
                "{}: the bound on `{}` is on no field and in no list",
                form.id,
                check.quantity
            );
        }
        // And nothing is in the list that is not a bound.
        for ranged in &form.unmodelled_ranges {
            assert!(
                spec.input_checks().any(|check| check.quantity == ranged),
                "{}: `{ranged}` is listed as unmodelled and no bound names it",
                form.id
            );
            assert!(
                !parameters.contains(&ranged.as_str()),
                "{}: `{ranged}` is both a parameter and unmodelled",
                form.id
            );
        }
        checked += spec.input_checks().count();
    }
    assert_eq!(checked, 82, "the bounds this test walked");
}

/// **The generated input table is the spec files' own, read a second time.**
///
/// `gen_model_inputs.py` compiles the input declarations out of the TOML, and nothing but this
/// test reads them again: a generator that skipped a section, dropped an `optional`, or lost an
/// enum's values would otherwise emit a form whose widgets are quietly the wrong ones.
#[test]
fn the_generated_input_table_is_the_models_own() {
    let mut walked = 0;
    for entry in azoth_process::model_inputs_gen::model_inputs() {
        let path = root().join("specs/models/process").join(format!(
            "{}.toml",
            entry.model.trim_start_matches("process.")
        ));
        let text = std::fs::read_to_string(&path).expect("the model spec is there");
        let document: toml::Value = toml::from_str(&text).expect("it parses");
        let declared = document["inputs"].as_table().expect("an `[inputs]` block");
        assert_eq!(
            entry.inputs.len(),
            declared.len(),
            "{}: the table has {} inputs and the spec {}",
            entry.model,
            entry.inputs.len(),
            declared.len()
        );
        for input in entry.inputs {
            let spec = declared
                .get(input.name)
                .unwrap_or_else(|| panic!("{}: `{}` is not an input", entry.model, input.name));
            assert_eq!(
                input.kind,
                spec.get("type")
                    .and_then(toml::Value::as_str)
                    .unwrap_or("quantity"),
                "{}: `{}`'s kind",
                entry.model,
                input.name
            );
            let values: Vec<&str> = spec
                .get("values")
                .and_then(toml::Value::as_array)
                .map(|values| values.iter().filter_map(toml::Value::as_str).collect())
                .unwrap_or_default();
            assert_eq!(
                input.values, values,
                "{}: `{}`'s enum",
                entry.model, input.name
            );
            assert_eq!(
                input.optional,
                spec.get("optional").and_then(toml::Value::as_bool) == Some(true),
                "{}: `{}`'s optional",
                entry.model,
                input.name
            );
            walked += 1;
        }
    }
    assert_eq!(walked, 309, "the inputs this test walked");
}

/// Every bound the models state names an input the models declare - `model_gen`'s `quantity`
/// strings held to `model_inputs_gen`'s names, two generated files from one set of specs.
#[test]
fn every_model_bound_names_an_input_the_model_declares() {
    for entry in azoth_process::model_inputs_gen::model_inputs() {
        let spec = azoth_process::model_gen::model(entry.model).expect("the model is compiled");
        for check in spec.input_checks() {
            assert!(
                entry
                    .inputs
                    .iter()
                    .any(|input| input.name == check.quantity),
                "{} bounds `{}`, which it does not declare",
                entry.model,
                check.quantity
            );
        }
    }
}

#[test]
fn a_dimension_id_is_the_inverse_of_the_exponents_it_names() {
    use azoth_core::unit_vocab_gen::{DIMENSION_EXPONENTS, dimension};

    // Two ids sharing exponents would make the inverse ambiguous, and one would be unreachable.
    for (index, (id, exponents)) in DIMENSION_EXPONENTS.iter().enumerate() {
        for (other, other_exponents) in &DIMENSION_EXPONENTS[index + 1..] {
            assert_ne!(
                exponents, other_exponents,
                "{id} and {other} are one dimension"
            );
        }
    }

    assert_eq!(dimension_id("Pa"), Some("pressure"));
    assert_eq!(dimension_id("K"), Some("thermodynamic_temperature"));
    assert_eq!(dimension_id("mol/s"), Some("molar_flow"));
    assert_eq!(dimension_id("dimensionless"), Some("dimensionless"));
    assert_eq!(
        dimension_id("psi"),
        None,
        "not a unit the vocabulary carries"
    );

    // Every unit resolves to a dimension the vocabulary names.
    for unit in azoth_core::unit_vocab_gen::UNIT_NAMES {
        if dimension(unit).is_some() {
            assert!(dimension_id(unit).is_some(), "`{unit}` has no dimension id");
        }
    }
}

/// **The hole the kind rule closes, driven through the checker rather than through the rule.**
///
/// A parameter that is declared, supplied and of the wrong kind: the name rules are satisfied, so
/// this validated clean and the run refused it - the same shape as the missing-parameter hole
/// before it.
#[test]
fn a_value_of_the_wrong_kind_reaches_the_checker_as_a_diagnostic() {
    let palette = palette();
    let document = |value: &str| {
        format!(
            r#"id = "f"
name = "f"
products = ["out"]

[[inputs]]
name = "in"
components = ["methane"]
n = 1.0
z = [1.0]
P = 5.0e5
T = 300.0

[[instances]]
id = "p1"
unit = "unit_ops.pump"
[instances.parameters]
outlet_pressure = {value}
isentropic_efficiency = 0.75

[[connections]]
from = "in"
to = "p1.inlet"

[[connections]]
from = "p1.outlet"
to = "out"
"#
        )
    };

    let broken = Flowsheet::from_toml(&document("\"high\"")).expect("it parses");
    let diagnostics = validate(&broken, &palette);
    let kind = diagnostics
        .iter()
        .find(|diagnostic| matches!(diagnostic, Diagnostic::ParameterKind { .. }))
        .unwrap_or_else(|| panic!("the kind is not reported: {diagnostics:?}"));
    assert_eq!(kind.code(), "parameter_kind");
    assert_eq!(kind.severity(), Severity::Error);
    assert_eq!(kind.location().section, "instances");
    assert_eq!(kind.location().path, "p1.parameters.outlet_pressure");
    assert_eq!(
        kind.message(),
        "`outlet_pressure` is a string, and its declaration says a number"
    );
    assert_eq!(
        kind.target(),
        Target::Parameter {
            node: "p1".into(),
            name: "outlet_pressure".into()
        }
    );

    // The same document with a number in it is clean, which is the half that says the rule is
    // about the kind and not about the parameter.
    let sound = Flowsheet::from_toml(&document("2.0e6")).expect("it parses");
    let diagnostics = validate(&sound, &palette);
    assert!(diagnostics.is_empty(), "{diagnostics:?}");
}
