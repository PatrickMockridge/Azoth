//! The shipped palette, loaded from `specs/unit_ops/`, must be well-formed.
//!
//! This is the CI gate for the palette: every entry parses as a `UnitOpSpec`,
//! names only dimensions the vocabulary knows, carries a provenance line, and has
//! ports with fields. It reads the files from disk rather than an embedded copy,
//! so it validates what a GUI would read, not a snapshot.

use std::collections::HashSet;
use std::path::Path;

use azoth_process::{UnitOpSpec, validate_palette};

fn unit_op_root() -> &'static Path {
    Path::new(concat!(env!("CARGO_MANIFEST_DIR"), "/../../specs/unit_ops"))
}

fn load_palette() -> Vec<UnitOpSpec> {
    let mut specs = Vec::new();

    let mut families: Vec<_> = std::fs::read_dir(unit_op_root())
        .expect("specs/unit_ops must exist")
        .map(|e| e.expect("readable entry").path())
        .filter(|p| p.is_dir())
        .collect();
    families.sort();

    for family in families {
        for entry in std::fs::read_dir(&family).expect("readable family dir") {
            let path = entry.expect("readable entry").path();
            if path.extension().and_then(|e| e.to_str()) != Some("toml") {
                continue;
            }
            let text = std::fs::read_to_string(&path)
                .unwrap_or_else(|e| panic!("{}: {e}", path.display()));
            let spec: UnitOpSpec =
                toml::from_str(&text).unwrap_or_else(|e| panic!("{}: {e}", path.display()));
            specs.push(spec);
        }
    }

    specs
}

#[test]
fn every_palette_entry_is_well_formed() {
    let specs = load_palette();
    assert!(
        !specs.is_empty(),
        "no palette entries found under specs/unit_ops"
    );

    let diags = validate_palette(&specs);
    assert!(diags.is_empty(), "palette defects: {diags:?}");
}

#[test]
fn every_palette_entry_is_named_sourced_and_typed() {
    let mut ids: HashSet<String> = HashSet::new();
    for spec in load_palette() {
        assert!(
            ids.insert(spec.id.clone()),
            "duplicate unit-op id: {}",
            spec.id
        );
        assert!(
            spec.source.is_some(),
            "{} carries no `source`, so its port shape has no provenance",
            spec.id
        );
        assert!(!spec.ports.is_empty(), "{} declares no ports", spec.id);
        for port in &spec.ports {
            assert!(
                !port.fields.is_empty(),
                "{}.{} declares no fields",
                spec.id,
                port.name
            );
        }
    }
}

/// The kernels that exist, by the name a palette entry's id ends with.
fn kernel_names() -> Vec<String> {
    let mut names: Vec<String> = std::fs::read_dir(kernel_root())
        .expect("the kernels directory must exist")
        .map(|e| e.expect("readable entry").path())
        .filter(|path| path.extension().and_then(|e| e.to_str()) == Some("rs"))
        .filter_map(|path| {
            let stem = path.file_stem()?.to_str()?.to_string();
            // `mod.rs` declares the module; it is not a kernel.
            (stem != "mod").then_some(stem)
        })
        .collect();
    names.sort();
    names
}

fn kernel_root() -> &'static Path {
    Path::new(concat!(env!("CARGO_MANIFEST_DIR"), "/src/kernels"))
}

fn model_root() -> &'static Path {
    Path::new(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../../specs/models/process"
    ))
}

/// The inputs a `process.*` model declares, read from its spec.
fn declared_inputs(name: &str) -> Option<HashSet<String>> {
    Some(
        declared_inputs_by_optionality(name)?
            .keys()
            .cloned()
            .collect(),
    )
}

/// The same, with whether each input is optional - which is the model's own word for whether a
/// caller may leave it out.
fn declared_inputs_by_optionality(name: &str) -> Option<toml::map::Map<String, toml::Value>> {
    let path = model_root().join(format!("{name}.toml"));
    let text = std::fs::read_to_string(&path).ok()?;
    let document: toml::Value =
        toml::from_str(&text).unwrap_or_else(|e| panic!("{}: {e}", path.display()));
    let inputs = document.get("inputs")?.as_table()?;
    Some(inputs.clone())
}

/// **The palette's `required` flag, held to the model spec that knows.**
///
/// A model's `[inputs]` block marks an input `optional = true` when a caller may leave it out, and
/// the palette now carries the same fact on `Param::required`. They are two spellings of one
/// thing, and this is the gate that keeps them one thing: without it, an entry could declare a
/// parameter optional that its kernel refuses without, or required one the kernel is happy to
/// default - and the checker's `MissingParameter` would be wrong in the direction that either
/// breaks working documents or misses a broken one.
///
/// **Scoped like its sibling** to the entries that have a kernel, because an entry with no kernel
/// has no model either and its parameters describe a machine that does not exist yet.
#[test]
fn every_required_parameter_is_the_ones_its_model_marks_required() {
    let kernels = kernel_names();
    let mut checked = 0;

    for spec in load_palette() {
        let short = spec.id.trim_start_matches("unit_ops.");
        if !kernels.iter().any(|kernel| kernel == short) {
            continue;
        }
        let inputs = declared_inputs_by_optionality(short).unwrap_or_else(|| {
            panic!(
                "{} has a kernel but no model spec to hold its `required` flags to",
                spec.id
            )
        });

        for (name, declaration) in &spec.parameters {
            let Some(input) = inputs.get(name) else {
                // The sibling gate reports a parameter no model takes.
                continue;
            };
            let optional = input.get("optional").and_then(toml::Value::as_bool) == Some(true);
            assert_eq!(
                declaration.required,
                !optional,
                "{} declares `{name}` with required = {}, and its model at \
                 specs/models/process/{short}.toml says {}",
                spec.id,
                declaration.required,
                if optional { "optional" } else { "required" }
            );
            checked += 1;
        }
    }
    assert!(
        checked > 0,
        "no entry was checked, so this test decided nothing"
    );
}

/// **The process layer's sibling of `test_cross_impl.py::test_signatures_agree_across_languages`.**
///
/// The calc layer has a gate holding a spec's declared parameters against the two
/// implementations' signatures. This tier had none, and the cost of that is on the record:
/// `unit_ops.separator` declared `efficiency`, which `Separator.run` never reads;
/// `unit_ops.throttling_valve` declared `valve_opening`, which only the unported sizing
/// mode reads; `unit_ops.heat_exchanger` declared `duty`, which is that class's *output*.
/// Each was found by hand, one at a time, while writing its kernel.
///
/// The gate is scoped to the entries that have a kernel, because an entry whose kernel is
/// still to come (the 18 outside this workstream) declares parameters nothing reads by
/// construction - that is the tranche's own backlog and not a defect to fail on. A kernel
/// that lands without a model, or a model that lands without the parameters its palette
/// entry declares, fails here.
#[test]
fn every_declared_parameter_reaches_its_kernel() {
    let kernels = kernel_names();
    assert!(!kernels.is_empty(), "no kernels found under src/kernels");

    for spec in load_palette() {
        let short = spec.id.trim_start_matches("unit_ops.");
        if !kernels.iter().any(|kernel| kernel == short) {
            continue;
        }
        let inputs = declared_inputs(short).unwrap_or_else(|| {
            panic!(
                "{} has a kernel at src/kernels/{short}.rs but no spec at \
                 specs/models/process/{short}.toml, so nothing holds its parameters to \
                 what the kernel reads",
                spec.id
            )
        });
        for parameter in spec.parameters.keys() {
            assert!(
                inputs.contains(parameter),
                "{} declares `{parameter}`, and its model at \
                 specs/models/process/{short}.toml does not take it. Either the kernel \
                 reads it and the model must declare it, or nothing reads it and the \
                 palette entry is describing a parameter it cannot honour.",
                spec.id
            );
        }
    }
}
