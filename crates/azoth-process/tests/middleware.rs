//! The interoperation surface: what a front-end reads.
//!
//! `docs/src/architecture/middleware.md` states the layer. Its first promise is that a
//! diagnostic reaches a widget as a structured record rather than as a line of text, and the two
//! things that can quietly go wrong with that are a code that stops matching its variant and a
//! target that points at nothing. Both are held here rather than by review.

use azoth_process::middleware::diagnostic::{DiagnosticRecord, records};
use azoth_process::{Diagnostic, NodeRole, Target, split_node_id};

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
