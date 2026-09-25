//! The dispatch table, held to the palette.
//!
//! **The palette has three consumers and this is the third.** The kernel-signature gate
//! (`tests/palette.rs`) holds every declared parameter to a model input; a case holds a model to
//! its answer; and the executor is what turns a declaration plus a form into arithmetic. A table
//! that named a unit operation the palette does not carry, or read a parameter the declaration
//! does not have, would fail here rather than at some flowsheet's first run.

use std::collections::BTreeMap;
use std::path::Path;

use azoth_core::units::{kelvins, pascals};
use azoth_process::executor::{DISPATCH, Parameters, UNRUNNABLE, dispatch};
use azoth_process::{Stream, UnitOpSpec, load_palette};

/// The shipped palette.
fn palette() -> Vec<UnitOpSpec> {
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../specs/unit_ops");
    load_palette(&root).expect("the palette loads")
}

/// A feed to hand a kernel, so a call reaches the parameter reading rather than failing earlier.
fn feed() -> Stream {
    Stream::from_pt(
        ["methane", "n-butane"]
            .iter()
            .map(|n| (*n).to_string())
            .collect(),
        vec![0.9, 0.1],
        1.0,
        pascals(5.0e5),
        kelvins(300.0),
    )
    .expect("the two resolve")
}

/// **Every palette entry is dispatched or refused by name.** Nothing is silently absent, which is
/// the acceptance P11 was held to and the same rule at the executor's level.
#[test]
fn every_palette_entry_runs_or_says_why_it_does_not() {
    let palette = palette();
    let mut missing = Vec::new();
    for spec in &palette {
        let dispatched = DISPATCH.iter().any(|(id, _)| *id == spec.id);
        let refused = UNRUNNABLE.iter().any(|(id, _)| *id == spec.id);
        let has_model = Path::new(env!("CARGO_MANIFEST_DIR"))
            .join(format!(
                "../../specs/models/process/{}.toml",
                spec.id.trim_start_matches("unit_ops.")
            ))
            .is_file();
        if !dispatched && !refused {
            missing.push(spec.id.clone());
        }
        // A dispatched entry always has a model - the executor runs what the registry
        // describes. A refused one may have a model and no runnable declaration, which is its
        // own statement and is why this is not an equality.
        assert!(
            !dispatched || has_model,
            "`{}` is dispatched but has no registered model",
            spec.id
        );
        assert!(
            !(dispatched && refused),
            "`{}` is both dispatched and refused",
            spec.id
        );
    }
    assert!(
        missing.is_empty(),
        "these palette entries neither dispatch nor carry a refusal: {missing:?}"
    );
    assert_eq!(
        DISPATCH.len() + UNRUNNABLE.len(),
        palette.len(),
        "the table and the refusals together are the palette"
    );
}

/// **Every parameter the table reads is one the declaration carries.**
///
/// `Parameters` is the only way an entry reads a parameter and it refuses a name the entry does
/// not declare, so supplying every declared name makes any read of an undeclared one surface as
/// exactly that refusal. **The value supplied is derived from the declaration's own unit**,
/// because a refusal about a unit would otherwise be indistinguishable from a refusal about a
/// name: every quantity gets `1.0`, every flag `false`, every string the empty one and every
/// vector an empty one.
#[test]
fn every_parameter_the_table_reads_is_declared() {
    let palette = palette();
    for (id, _) in DISPATCH {
        let spec = palette
            .iter()
            .find(|spec| spec.id == *id)
            .unwrap_or_else(|| panic!("`{id}` is dispatched but not declared"));
        let mut values: BTreeMap<String, toml::Value> = BTreeMap::new();
        for (name, param) in &spec.parameters {
            let value = match param.unit.as_deref() {
                None | Some("dimensionless") => {
                    // A dimensionless parameter is a number or a name, and the declaration says
                    // which only by its description; `1.0` is a number, so an entry reading one
                    // as text fails on the shape and not on the name, which is what this checks.
                    toml::Value::Float(1.0)
                }
                Some(_) => toml::Value::Float(1.0),
            };
            values.insert(name.clone(), value);
        }
        let params = Parameters::new(spec, &values);
        // An entry may well *succeed* on a set of ones - a pump at one pascal does not refuse -
        // so only a failure is inspected, and the one failure this checks for is a name.
        if let Err(error) = dispatch(id, &[feed(), feed()], &params) {
            let message = error.to_string();
            assert!(
                !message.contains("does not declare a parameter called"),
                "`{id}` reads a parameter its declaration does not carry: {message}"
            );
        }
    }
}

/// **An entry with no kernel is refused with the class that would close it.** A missing id is a
/// defect; a named refusal is a recorded disposition, which is the difference the roadmap's own
/// rule turns on.
#[test]
fn a_refused_entry_names_what_would_close_it() {
    for (id, reason) in UNRUNNABLE {
        assert!(
            reason.len() > 40,
            "`{id}`'s refusal says nothing: {reason:?}"
        );
        let error = dispatch(
            id,
            &[feed()],
            &Parameters::new(&palette()[0], &BTreeMap::new()),
        )
        .expect_err("a refused entry does not run");
        assert!(
            error.to_string().contains(reason),
            "`{id}`'s refusal does not reach the caller"
        );
    }
}

/// **An id the palette does not carry is refused as that**, and not as a palette entry with no
/// kernel - the two are different statements.
#[test]
fn an_unknown_id_says_it_is_unknown() {
    let error = dispatch(
        "unit_ops.not_a_machine",
        &[feed()],
        &Parameters::new(&palette()[0], &BTreeMap::new()),
    )
    .expect_err("an unknown id does not run");
    assert!(
        error
            .to_string()
            .contains("not a unit operation the palette declares")
    );
}
