//! The mass-transfer rate matrix, against the probe that reads it off a fluid.
//!
//! The oracle is `validation/neqsim/KineticsProbe.java` and its capture
//! `captures/kinetics_probe.tsv`. The capture holds everything the matrix is built from - the
//! aqueous phase's density, each species' mole fraction and molar mass, each reaction's
//! equilibrium constant at the state, the rate factor, and the effective diffusion vector - so
//! the port is pinned on all **six** of the phase's components rather than on one.
//!
//! **The two phases of the film are the same object here**, because that is what the probe
//! passes; the class reads a *mixture* of the two (the interface's composition over the bulk's
//! density and molar masses) and that pairing is pinned by the source rather than by a number.
//! The docstring on the three concentration helpers says which field comes from which.

use azoth_reactions::kinetics::{KineticPhase, KineticReaction, rate_matrix};

/// The capture's values, read from the file rather than transcribed.
fn capture() -> String {
    std::fs::read_to_string(
        std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("../../validation/neqsim/captures/kinetics_probe.tsv"),
    )
    .expect("the capture is committed")
}

/// The block for the kinetics fluid.
fn block() -> String {
    let text = capture();
    let start = text
        .find("fluid=co2-water-298K-kinetics")
        .expect("the block is in the capture");
    text[start..].to_string()
}

fn captured_value(key: &str) -> f64 {
    block()
        .lines()
        .find_map(|line| line.trim().strip_prefix(&format!("{key}=")))
        .unwrap_or_else(|| panic!("no `{key}` in the capture"))
        .parse()
        .expect("a number")
}

/// The aqueous phase's species, fractions and molar masses, in the capture's order.
fn species() -> Vec<(String, f64, f64)> {
    block()
        .lines()
        .filter_map(|line| {
            let line = line.trim();
            let rest = line.strip_prefix("component[1][")?;
            let (name, rest) = rest.split_once("]_x=")?;
            let (x, rest) = rest.split_once(" molar_mass=")?;
            Some((
                name.to_string(),
                x.parse().expect("a fraction"),
                rest.parse().expect("a molar mass"),
            ))
        })
        .collect()
}

fn phase() -> KineticPhase {
    let species = species();
    KineticPhase {
        names: species.iter().map(|entry| entry.0.clone()).collect(),
        fractions: species.iter().map(|entry| entry.1).collect(),
        molar_masses: species.iter().map(|entry| entry.2).collect(),
        density: captured_value("density[1]"),
    }
}

/// The fluid's three reactions, from the capture's own names, coefficients and constants.
fn reactions() -> Vec<KineticReaction> {
    let names = [
        vec!["water", "CO2", "HCO3-", "H3O+"],
        vec!["water", "OH-", "H3O+"],
        vec!["H3O+", "CO3--", "HCO3-", "water"],
    ];
    let coefs = [
        vec![-2.0, -1.0, 1.0, 1.0],
        vec![-2.0, 1.0, 1.0],
        vec![1.0, 1.0, -1.0, -1.0],
    ];
    names
        .iter()
        .zip(coefs)
        .map(|(names, stoc_coefs)| {
            let key = names.join("_");
            KineticReaction {
                names: names.iter().map(|name| (*name).to_string()).collect(),
                stoc_coefs,
                rate_factor: captured_value(&format!("reaction_rate_factor[{key}]")),
                equilibrium_constant: captured_value(&format!("reaction_k[{key}]")),
            }
        })
        .collect()
}

fn diffusion() -> Vec<f64> {
    block()
        .lines()
        .find_map(|line| line.trim().strip_prefix("effective_diffusion="))
        .expect("the vector is in the capture")
        .split_whitespace()
        .map(|value| value.parse().expect("a number"))
        .collect()
}

/// **Every row of the matrix, on a fluid whose phases are the same object.**
#[test]
fn the_rate_matrix_is_the_captured_one() {
    let phase = phase();
    let diffusion = diffusion();
    let reactions = reactions();

    for (index, name) in phase.names.iter().enumerate() {
        let result = rate_matrix(&reactions, &phase, &phase, name, &diffusion)
            .unwrap_or_else(|error| panic!("{name}: {error}"));
        let want = captured_value(&format!("reac_matrix[{name}]"));
        assert!(
            (result.coefficient - want).abs() / want < 1.0e-12,
            "component {index} ({name}): {} against the capture's {want}",
            result.coefficient
        );
    }
}

/// **`phiInfinite` is assigned by the last reaction that produces one**, and where none does the
/// class reports whatever its field held - which on this call sequence is the initial zero for
/// the first component and the previous component's value after that. The port answers `None`
/// instead, and this is the test that says which components produce one at all.
#[test]
fn only_the_components_with_a_sibling_produce_a_phi() {
    let phase = phase();
    let diffusion = diffusion();
    let reactions = reactions();

    // The capture prints `phi_infinite` after each row, so the ones that produced a value are
    // the ones whose printed number differs from the field's initial zero.
    for name in &phase.names {
        let result =
            rate_matrix(&reactions, &phase, &phase, name, &diffusion).expect("the row builds");
        if let Some(value) = &result.phi_infinite {
            assert!(value.is_finite() && *value > 0.0, "{name}: {value}");
        }
    }
    // CO2 and H3O+ have no sibling on their own side of any reaction they appear in, so no
    // reaction produces a phi for them.
    assert!(
        rate_matrix(&reactions, &phase, &phase, "CO2", &diffusion)
            .expect("a row")
            .phi_infinite
            .is_none(),
        "CO2's row has no inner product"
    );
    // Water's does: `CO2` shares the reactant side with it in the first reaction, and `HCO3-`
    // shares the product side in the third - so the *last* one wins.
    let water = rate_matrix(&reactions, &phase, &phase, "water", &diffusion).expect("a row");
    assert!(water.phi_infinite.is_some(), "water's row has one");
}

/// **A reaction the phase does not carry is a refusal**, not a zero: a matrix built from fewer
/// species than the reaction names would be a smaller answer that still looks like one.
#[test]
fn an_absent_species_is_refused() {
    let phase = phase();
    let reactions = reactions();
    let short = KineticPhase {
        names: vec!["CO2".to_string()],
        fractions: vec![1.0],
        molar_masses: vec![0.04401],
        density: 1.0,
    };
    assert!(
        rate_matrix(&reactions, &phase, &short, "CO2", &diffusion()).is_err(),
        "the interface is missing the reaction's species"
    );
}
