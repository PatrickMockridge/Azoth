//! The agent's tool schema, projected from the command model.
//!
//! `middleware.md`: "an agent's tool is `@calc` and its call is `*tool` — the same reflection
//! surface as the GUI, not a third one". What that means mechanically is here: **one tool per
//! [`Command`]**, so the set of things an agent may do and the set of things a canvas may do are
//! the same set, and the test beside this holds the two lists to each other.
//!
//! **Only edits are tools.** The page names `list_unit_ops`, `validate`, `read_stream`, `load` and
//! `save` as the agentic surface, and none of them needs a tool: every call returns the whole
//! [`crate::middleware::envelope::Envelope`], so reading a stream is reading the answer to the
//! call that changed it, and loading, saving and validating are the session's own openings rather
//! than edits to it. A read tool would return the same document as every other tool.
//!
//! **What is derived and what is written.** The names are the enum's, and the test derives them
//! from it; the schemas are written, because a JSON-Schema derive over fourteen commands would be
//! a dependency and an indirection for an object a hand-written one states more plainly. The
//! `unit` an `add_instance` takes is an enum of the palette's own ids, which is the one place the
//! schema reaches past the command to the declaration.

use serde::Serialize;
use serde_json::{Value, json};

use crate::unit_op::UnitOpSpec;

/// One tool an agent may call.
#[derive(Debug, Clone, Serialize)]
pub struct Tool {
    /// The command's own name, which is also the tag its JSON carries.
    pub name: &'static str,
    /// One sentence: what it does, and the rule a caller has to know.
    pub description: &'static str,
    /// The command object's JSON Schema, `command` included.
    pub input_schema: Value,
}

/// The schema of every command, with the palette's ids where a unit op is named.
#[must_use]
pub fn tools(palette: &[UnitOpSpec]) -> Vec<Tool> {
    let units: Vec<&str> = palette.iter().map(|spec| spec.id.as_str()).collect();
    let string = |description: &str| json!({ "type": "string", "description": description });
    let number = |description: &str| json!({ "type": "number", "description": description });

    vec![
        tool(
            "add_instance",
            "Add a unit operation. An id that already exists is replaced, so the call is retryable.",
            vec![
                req("id", string("The instance's name, e.g. `p1`.")),
                req(
                    "unit",
                    json!({ "type": "string", "enum": units,
                            "description": "The palette entry to instantiate." }),
                ),
                opt(
                    "parameters",
                    json!({ "type": "object", "description":
                            "Parameter values, keyed by the entry's own names." }),
                ),
            ],
        ),
        tool(
            "remove_instance",
            "Remove a unit operation, and every connection and recycle that named it.",
            vec![req("id", string("The instance to remove."))],
        ),
        tool(
            "connect",
            "Join a feed or an outlet to an inlet or a product. An outlet is `instance.port`, and \
             a `many` outlet takes a position: `split1.products[0]`.",
            vec![
                req("from", string("Where the stream is produced.")),
                req("to", string("Where it is consumed.")),
            ],
        ),
        tool(
            "disconnect",
            "Remove the first connection or recycle that joins these two endpoints.",
            vec![
                req("from", string("The producing endpoint, as it was written.")),
                req("to", string("The consuming endpoint, as it was written.")),
            ],
        ),
        tool(
            "set_parameter",
            "Set one parameter of one instance. The value's JSON kind follows the declaration: a \
             number for a quantity, a list for a vector, a switch for a boolean, a name for an enum.",
            vec![
                req("instance", string("The instance that owns it.")),
                req("name", string("The parameter's declared name.")),
                req(
                    "value",
                    json!({ "description": "The value, in the declaration's kind." }),
                ),
            ],
        ),
        tool(
            "unset_parameter",
            "Remove a parameter's value, leaving it to the declaration's own default.",
            vec![
                req("instance", string("The instance that owns it.")),
                req("name", string("The parameter's declared name.")),
            ],
        ),
        tool(
            "add_input",
            "Add a boundary feed. `h` is not settable: it is a state function of `T`, `P` and `z` \
             and the kernel derives it.",
            vec![
                req("name", string("The feed's name.")),
                req(
                    "components",
                    json!({ "type": "array", "items": { "type": "string" },
                            "description": "The fluid's substances, by name." }),
                ),
                req("n", number("Molar flow, mol/s.")),
                req(
                    "z",
                    json!({ "type": "array", "items": { "type": "number" },
                            "description": "Mole fractions, one per substance." }),
                ),
                req("P", number("Pressure, Pa.")),
                req("T", number("Temperature, K.")),
            ],
        ),
        tool(
            "remove_input",
            "Remove a boundary feed, and every connection that drew from it.",
            vec![req("name", string("The feed to remove."))],
        ),
        tool(
            "add_product",
            "Declare a boundary product. An output is calculated, so a product is a name only.",
            vec![req("name", string("The product's name."))],
        ),
        tool(
            "remove_product",
            "Remove a boundary product, and every connection that wrote to it.",
            vec![req("name", string("The product to remove."))],
        ),
        tool(
            "add_recycle",
            "Declare a tear: the connection that closes a loop, under a name of its own.",
            vec![
                req("stream", string("The tear's name.")),
                req("from", string("The outlet the loop leaves from.")),
                req("to", string("The inlet the loop returns to.")),
            ],
        ),
        tool(
            "remove_recycle",
            "Remove a tear, leaving its endpoints as an ordinary pair of connections.",
            vec![req("stream", string("The tear to remove."))],
        ),
        tool(
            "set_recycle",
            "Set one convergence parameter of one tear. A tear that does not converge is reported, \
             not raised.",
            vec![
                req("stream", string("The tear's name.")),
                req(
                    "field",
                    json!({ "type": "string",
                            "enum": ["flow_tolerance", "composition_tolerance",
                                     "temperature_tolerance", "pressure_tolerance",
                                     "max_iterations", "minimum_flow", "acceleration_method"],
                            "description": "Which of the tear's seven settings." }),
                ),
                req(
                    "value",
                    json!({ "description":
                            "A number, except for `acceleration_method`, which is a name." }),
                ),
            ],
        ),
        tool(
            "set_position",
            "Place a node on the canvas, by node id: `instance:sep1`, `input:feed_1`, \
             `product:out`.",
            vec![
                req("node", string("The node's id.")),
                req("x", number("The column, in pixels.")),
                req("y", number("The row, in pixels.")),
            ],
        ),
    ]
}

/// One property of a tool's schema, and whether a call must carry it.
struct Property {
    name: &'static str,
    schema: Value,
    required: bool,
}

/// A property a call must carry.
fn req(name: &'static str, schema: Value) -> Property {
    Property {
        name,
        schema,
        required: true,
    }
}

/// A property a call may leave out, which is the same thing the command's own `#[serde(default)]`
/// says about it.
fn opt(name: &'static str, schema: Value) -> Property {
    Property {
        name,
        schema,
        required: false,
    }
}

/// One tool, from its own properties.
///
/// **The `command` tag is added here rather than written fourteen times.** It is the same string
/// as the tool's name, and a schema that spelled it differently would describe a call the
/// deserialiser refuses with no clue why.
fn tool(name: &'static str, description: &'static str, properties: Vec<Property>) -> Tool {
    let mut object = serde_json::Map::new();
    let mut required = vec![json!("command")];
    object.insert(
        "command".to_string(),
        json!({ "const": name, "description": "The command's name." }),
    );
    for Property {
        name: property,
        schema,
        required: needs,
    } in properties
    {
        if needs {
            required.push(json!(property));
        }
        object.insert(property.to_string(), schema);
    }
    Tool {
        name,
        description,
        input_schema: json!({
            "type": "object",
            "properties": object,
            "required": required,
            "additionalProperties": false,
        }),
    }
}
