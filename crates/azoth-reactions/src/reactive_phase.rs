//! Which of a fluid's phases the reactions are solved in.
//!
//! `ChemicalReactionOperations.getReactivePhaseIndex` scans the system's phases in
//! order and returns the first whose type name is `aqueous`; failing that the first
//! whose type name is `liquid` or `oil`; and failing that `-1`, which makes
//! `solveChemEq` return `false` without attempting anything.
//!
//! **The two-word fallback is a live branch and not a curiosity.** Its own comment says
//! it exists for the state before a flash has assigned phase types, and it does: every
//! one of the five captured fluids takes it before its flash and the `aqueous` branch
//! after, so the same fluid reaches both.
//!
//! **The scan is over labels because azoth has no system phase array.** A calculation
//! here is handed the phase it works on, so the caller runs the scan and
//! [`is_reactive_phase`] answers whether what it found is reactive. The index form is
//! kept beside it because the capture is a list of labels and an index, and the two are
//! the same rule rather than two readings of it.

/// The three phase type names NeqSim solves reactions in.
pub const REACTIVE_PHASE_LABELS: [&str; 3] = ["aqueous", "liquid", "oil"];

/// Whether a phase of this type takes the reaction solve.
///
/// The comparison is case-insensitive, as `equalsIgnoreCase` is. A name that is none of
/// the three - `gas` above all - is the skip.
#[must_use]
pub fn is_reactive_phase(label: &str) -> bool {
    let folded = label.trim();
    REACTIVE_PHASE_LABELS
        .iter()
        .any(|candidate| folded.eq_ignore_ascii_case(candidate))
}

/// The first phase of `labels` that takes the reaction solve.
///
/// **`aqueous` is looked for over every phase before `liquid` or `oil` is looked for
/// over any**, which is two passes and not one: a liquid phase listed ahead of an
/// aqueous one does not win. Within each pass the first match in the caller's order is
/// the answer, so the order is part of the input.
///
/// **`None` is the skip, and the skip is not a failure.** A fluid with no aqueous and no
/// liquid phase has no phase for a water-based equilibrium to be solved in, and NeqSim
/// returns `-1` there for `solveChemEq` to return `false` on. An answer of `None` and a
/// solve that failed are different results and the caller can tell them apart.
#[must_use]
pub fn reactive_phase_index(labels: &[&str]) -> Option<usize> {
    let position =
        |wanted: &dyn Fn(&str) -> bool| labels.iter().position(|label| wanted(label.trim()));
    let aqueous = |label: &str| label.eq_ignore_ascii_case("aqueous");
    let liquid_or_oil =
        |label: &str| label.eq_ignore_ascii_case("liquid") || label.eq_ignore_ascii_case("oil");
    position(&aqueous).or_else(|| position(&liquid_or_oil))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn an_aqueous_phase_wins_over_a_liquid_one_listed_before_it() {
        // Two passes, not one: the fallback does not pre-empt the first pass.
        assert_eq!(reactive_phase_index(&["gas", "liquid", "aqueous"]), Some(2));
    }

    #[test]
    fn the_fallback_takes_liquid_or_oil_in_the_order_they_appear() {
        assert_eq!(reactive_phase_index(&["gas", "oil", "liquid"]), Some(1));
        assert_eq!(reactive_phase_index(&["gas", "liquid", "oil"]), Some(1));
    }

    #[test]
    fn a_fluid_with_neither_is_a_skip_and_not_an_error() {
        assert_eq!(reactive_phase_index(&["gas"]), None);
        assert_eq!(reactive_phase_index(&[]), None);
        assert_eq!(reactive_phase_index(&["gas", "hydrate"]), None);
    }

    #[test]
    fn the_case_of_the_label_does_not_decide_it() {
        assert_eq!(reactive_phase_index(&["GAS", "Aqueous"]), Some(1));
        assert!(is_reactive_phase("AQUEous"));
        assert!(!is_reactive_phase("gas"));
    }
}
