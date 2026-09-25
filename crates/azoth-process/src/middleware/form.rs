//! A unit operation's declaration, as a form.
//!
//! `docs/src/architecture/middleware.md` asks for "one uniform object per unit op: its declared
//! free variables with units, dimensions and valid ranges", and calls it the adequacy claim made
//! inspectable - a unit-op window *is* its port declaration. This module builds that object.
//!
//! **Three declarations meet here and no new one is invented.** The palette
//! (`specs/unit_ops/`) carries the ports, the parameters, their units and their `required`
//! flags; the model-input table (`model_inputs_gen`) carries the *kind*, which is what chooses a
//! widget; and the model's own bounds (`model_gen`) carry the ranges. A form, the checker and the
//! kernel are therefore reading one declaration between them rather than three that can drift.

use azoth_core::unit_vocab_gen::{DIMENSION_EXPONENTS, dimension};
use serde::Serialize;

use crate::channel::{Direction, FieldType, Multiplicity, Shape};
use crate::executor::{DISPATCH, UNRUNNABLE};
use crate::model_gen;
use crate::model_inputs_gen;
use crate::unit_op::UnitOpSpec;

/// The dimension a unit is at, as a vocabulary id.
///
/// **The inverse of `unit_vocab_gen::dimension_exponents`**, which the generated vocabulary
/// carries as a forward map only: `dimension(unit)` answers exponents, and a form needs the
/// *name* - `pressure` rather than `[-1, 1, -2, 0, 0, 0, 0]`. Equal exponents are the same
/// dimension by definition, so the search is well defined; `middleware.rs` holds that the table
/// has no two ids sharing one.
#[must_use]
pub fn dimension_id(unit: &str) -> Option<&'static str> {
    let exponents = dimension(unit)?;
    DIMENSION_EXPONENTS
        .iter()
        .find(|(_, candidate)| *candidate == exponents)
        .map(|(id, _)| *id)
}

/// One parameter, as a field on a form.
#[derive(Debug, Clone, Serialize)]
pub struct FormParameter {
    pub name: String,
    /// The kind of value: `quantity`, `vector`, `boolean`, `enum`, `components`, `matrix`,
    /// `string`, or `unknown` where no model declares the input.
    ///
    /// **The field neither the palette nor the checker has.** `Param` carries a unit and no
    /// type, so a form given only the palette cannot tell `split_factors` from
    /// `outlet_pressure` - both are a unit and a `required` flag - and cannot choose between a
    /// vector editor and a number box. `unknown` is a real answer rather than a default: a
    /// widget cannot be chosen, and the entry's own refusal says why.
    pub kind: &'static str,
    /// Whether a kernel can run without it; the model's `optional` is its negation, held to
    /// agreement by `tests/palette.rs`.
    pub required: bool,
    /// The canonical unit a form shows and a document states the value in, where the
    /// declaration names one.
    pub unit: Option<String>,
    /// The dimension that unit is at, as a vocabulary id.
    pub dimension: Option<&'static str>,
    /// An `enum`'s allowed values, in the spec's order; empty for every other kind.
    pub values: &'static [&'static str],
    /// The palette's own sentence, which is what a form shows beside the field.
    pub description: String,
    /// The model's bounds on this input, in the spec's order.
    pub range: Vec<FormRange>,
}

/// One bound, as a form draws it.
#[derive(Debug, Clone, Serialize)]
pub struct FormRange {
    pub min: Option<f64>,
    pub min_inclusive: bool,
    pub max: Option<f64>,
    pub max_inclusive: bool,
    /// A single forbidden value, where the bound is an exclusion rather than an interval.
    pub equals: Option<f64>,
    /// `outside` where the interval is allowed and beyond it violates, `inside` where the
    /// interval is the forbidden band.
    pub band: &'static str,
    /// `error` or `warning`, which is what a field marks a violation as.
    pub severity: &'static str,
    /// The stable code the violation carries, e.g. `OUT_OF_VALID_RANGE`.
    pub code: &'static str,
    /// Why the bound exists, which is what the field shows on a violation.
    pub rationale: &'static str,
}

/// One port, as a node's handles come from it.
#[derive(Debug, Clone, Serialize)]
pub struct FormPort {
    pub name: String,
    /// `in` or `out`.
    pub direction: &'static str,
    /// `one` or `many`.
    pub multiplicity: &'static str,
    pub fields: Vec<FormField>,
}

/// One field of a port's record.
#[derive(Debug, Clone, Serialize)]
pub struct FormField {
    pub name: String,
    pub dimension: String,
    /// `scalar` or `vector`.
    pub shape: &'static str,
}

/// One palette entry, as a form.
#[derive(Debug, Clone, Serialize)]
pub struct FormSpec {
    pub id: String,
    pub name: String,
    /// The provenance line the palette carries, e.g. the NeqSim class.
    pub source: Option<String>,
    /// The model id, where one exists.
    pub model: Option<&'static str>,
    /// Whether the executor has a kernel for it.
    pub runnable: bool,
    /// Why not, verbatim from the refusal table, where it does not.
    pub refusal: Option<&'static str>,
    pub ports: Vec<FormPort>,
    pub parameters: Vec<FormParameter>,
    /// The bounds the model puts on inputs no parameter of this entry carries.
    ///
    /// **Named rather than dropped.** A model bounds its stream fields too - the compressor's
    /// `inlet_t` must be a temperature - and a form has no field for those, so a document can be
    /// refused by a rule no field shows. An empty list is the common case.
    pub unmodelled_ranges: Vec<String>,
}

/// The form for one palette entry.
#[must_use]
pub fn form(spec: &UnitOpSpec) -> FormSpec {
    let model_inputs = model_inputs_gen::inputs_for(&spec.id);
    let model = model_inputs.and_then(|entry| model_gen::model(entry.model));
    let checks = || {
        model
            .iter()
            .flat_map(|spec| spec.input_checks())
            .collect::<Vec<_>>()
    };

    let parameters = spec
        .parameters
        .iter()
        .map(|(name, declaration)| {
            let input = model_inputs.and_then(|entry| entry.inputs.iter().find(|i| i.name == name));
            FormParameter {
                name: name.clone(),
                kind: input.map_or("unknown", |input| input.kind),
                required: declaration.required,
                unit: declaration.unit.clone(),
                dimension: declaration.unit.as_deref().and_then(dimension_id),
                values: input.map_or(&[], |input| input.values),
                description: declaration.description.clone(),
                range: checks()
                    .into_iter()
                    .filter(|check| check.quantity == name)
                    .map(range)
                    .collect(),
            }
        })
        .collect();

    let unmodelled_ranges = checks()
        .into_iter()
        .map(|check| check.quantity)
        .filter(|quantity| !spec.parameters.contains_key(*quantity))
        .map(str::to_string)
        .collect();

    FormSpec {
        id: spec.id.clone(),
        name: spec.name.clone(),
        source: spec.source.as_ref().map(|source| source.standard.clone()),
        model: model_inputs.map(|entry| entry.model),
        runnable: DISPATCH.iter().any(|(id, _)| *id == spec.id),
        refusal: UNRUNNABLE
            .iter()
            .find(|(id, _)| *id == spec.id)
            .map(|(_, reason)| *reason),
        ports: spec.ports.iter().map(port).collect(),
        parameters,
        unmodelled_ranges,
    }
}

/// Every entry's form, in the order the palette was given.
#[must_use]
pub fn forms(specs: &[UnitOpSpec]) -> Vec<FormSpec> {
    specs.iter().map(form).collect()
}

fn port(port: &crate::channel::Port) -> FormPort {
    FormPort {
        name: port.name.clone(),
        direction: match port.direction {
            Direction::In => "in",
            Direction::Out => "out",
        },
        multiplicity: match port.multiplicity {
            Multiplicity::One => "one",
            Multiplicity::Many => "many",
        },
        fields: port.fields.iter().map(field).collect(),
    }
}

fn field((name, field): (&String, &FieldType)) -> FormField {
    FormField {
        name: name.clone(),
        dimension: field.dimension.clone(),
        shape: match field.shape {
            Shape::Scalar => "scalar",
            Shape::Vector => "vector",
        },
    }
}

fn range(check: &azoth_core::RangeCheck) -> FormRange {
    FormRange {
        min: check.min,
        min_inclusive: check.min_inclusive,
        max: check.max,
        max_inclusive: check.max_inclusive,
        equals: check.equals,
        band: match check.band {
            azoth_core::Band::Outside => "outside",
            azoth_core::Band::Inside => "inside",
        },
        severity: match check.severity {
            azoth_core::Severity::Warning => "warning",
            azoth_core::Severity::Error => "error",
        },
        code: check.code.as_str(),
        rationale: check.rationale,
    }
}
