//! One calculation, as a sheet's form.
//!
//! The shapes mirror the unit-operation form where a reader already knows it: an input is a
//! `FormParameter` in all but name, so the editor reuses its field. What is added is the
//! calculation's own declarations - its namespace, its outputs, its worked example and the
//! refusal that keeps a table-shaped calculation out of the sheet editor for this tranche.

use serde::Serialize;

/// One worked-example value, as the spec wrote it.
#[derive(Debug, Clone, Copy, Serialize)]
#[serde(untagged)]
pub enum ExampleValue {
    Number(f64),
    Flag(bool),
    Name(&'static str),
}

/// One named worked-example input.
#[derive(Debug, Clone, Copy)]
pub struct CalcExample {
    pub name: &'static str,
    pub value: ExampleValue,
}

/// One bound, as a field draws it.
#[derive(Debug, Clone, Copy, Serialize)]
pub struct CalcRange {
    pub min: Option<f64>,
    pub min_inclusive: bool,
    pub max: Option<f64>,
    pub max_inclusive: bool,
    pub equals: Option<f64>,
    pub band: &'static str,
    pub severity: &'static str,
    pub code: &'static str,
    pub rationale: &'static str,
}

/// One input, with the unit and the kind the field needs.
#[derive(Debug, Clone, Copy, Serialize)]
pub struct CalcInput {
    pub name: &'static str,
    /// `quantity`, `enum`, `boolean` or `string`; the four the sheet editor can draw.
    pub kind: &'static str,
    pub required: bool,
    pub unit: Option<&'static str>,
    pub dimension: Option<&'static str>,
    pub values: &'static [&'static str],
    pub description: &'static str,
    pub range: &'static [CalcRange],
    pub default: Option<ExampleValue>,
}

/// One output, as a result table labels it.
#[derive(Debug, Clone, Copy, Serialize)]
pub struct CalcOutput {
    pub name: &'static str,
    pub kind: &'static str,
    pub unit: Option<&'static str>,
    pub dimension: Option<&'static str>,
    pub values: &'static [&'static str],
    pub description: &'static str,
    pub optional: bool,
}

/// One registered calculation, as the picker and the sheet read it.
#[derive(Debug, Clone, Copy, Serialize)]
pub struct CalcForm {
    pub id: &'static str,
    pub name: &'static str,
    pub namespace: &'static str,
    pub available: bool,
    pub refusal: Option<&'static str>,
    pub equation: Option<&'static str>,
    pub assumptions: &'static [&'static str],
    pub references: &'static [&'static str],
    pub inputs: &'static [CalcInput],
    pub outputs: &'static [CalcOutput],
    /// The worked example, serialised as an object by [`crate::catalogue`].
    #[serde(skip)]
    pub worked_example: &'static [CalcExample],
    /// The example's expected outputs, in the units the spec declares.
    ///
    /// **Not on the wire**: it is the test's pin and not something a sheet shows. The UI shows
    /// what the kernel answered, not what the spec said it should.
    #[serde(skip)]
    pub expected: &'static [CalcExample],
    /// The relative tolerance the worked example is checked at.
    #[serde(skip)]
    pub tolerance: f64,
}
