//! The two-film transfer coefficient: the class's scalar form and its Maxwell-Stefan matrix.

use azoth_core::AzothError;
use azoth_reactions::linalg::solve_lu;

use super::fallbacks::{DEFAULT_GAS_DIFFUSIVITY, DEFAULT_LIQUID_DIFFUSIVITY};
use super::phase::PhaseView;

/// `combineFilmFluxes`: the rate-limiting harmonic mean of the two flux densities.
///
/// **Three guards, and the third is the one a simplification destroys.** A non-finite flux is
/// zero; a flux below `1e-30` in magnitude is zero; and **two fluxes of opposite sign are
/// zero** rather than a harmonic mean of them. The last is not defensive: the two films drive
/// the same component in opposite directions exactly where the interface has crossed the bulk,
/// and averaging them there would report a transfer the two films do not agree on.
pub fn combine_film_fluxes(gas_flux: f64, liquid_flux: f64) -> f64 {
    if !gas_flux.is_finite() || !liquid_flux.is_finite() {
        return 0.0;
    }
    if gas_flux.abs() < 1.0e-30 || liquid_flux.abs() < 1.0e-30 {
        return 0.0;
    }
    if gas_flux.signum() != liquid_flux.signum() {
        return 0.0;
    }
    let magnitude = 1.0 / (1.0 / gas_flux.abs() + 1.0 / liquid_flux.abs());
    gas_flux.signum() * magnitude
}

/// `scaleFilmCoefficient`: a film coefficient by the diffusivity's ratio to a reference,
/// clamped to `[0.02, 50]` times the base.
///
/// **Both halves of the ratio fall back on each other**, which is the class's own nesting:
/// `finitePositive(referenceDiffusivity, diffusivity)` and `finitePositive(diffusivity,
/// reference)` - so a state where neither is positive answers a clamped zero rather than a
/// division by nothing.
pub fn scale_film_coefficient(base: f64, diffusivity: f64, reference: f64) -> f64 {
    let reference = finite_positive(reference, diffusivity);
    let scaled = base * finite_positive(diffusivity, reference) / reference;
    clamp(scaled, base * 0.02, base * 50.0)
}

/// `binaryDiffusivity`: the pair coefficient, with the class's two-step fallback.
///
/// Measured on the class's own absorber state, the first step answers: the gas phase's
/// `getDiffusionCoefficient(0, 1)` is `9.457085848059925e-7` m²/s, which is the value
/// `eos.phase_transport`'s own capture holds for methane/CO2 at 313 K and 50 bar. The second
/// step - the *effective* coefficient - is the one that answers `0.0` upstream, because
/// NeqSim's flash path never populates that vector.
pub fn binary_diffusivity(view: &PhaseView, first: usize, second: usize, gas_phase: bool) -> f64 {
    if first == second {
        return 0.0;
    }
    let matrix = view
        .transport
        .d_binary
        .get(first)
        .and_then(|row| row.get(second))
        .copied()
        .unwrap_or(0.0);
    if matrix > 0.0 && matrix.is_finite() {
        return matrix;
    }
    let effective = view
        .transport
        .d_effective
        .get(first)
        .copied()
        .unwrap_or(0.0);
    if effective > 0.0 && effective.is_finite() {
        return effective;
    }
    if gas_phase {
        DEFAULT_GAS_DIFFUSIVITY
    } else {
        DEFAULT_LIQUID_DIFFUSIVITY
    }
}

/// `binaryFilmCoefficient`.
fn binary_film_coefficient(
    view: &PhaseView,
    first: usize,
    second: usize,
    base: f64,
    reference: f64,
    gas_phase: bool,
) -> f64 {
    if first == second {
        return base;
    }
    scale_film_coefficient(
        base,
        binary_diffusivity(view, first, second, gas_phase),
        reference,
    )
}

/// `mixtureDiffusivityForComponent`: the Wilke resistance sum over the other components.
fn mixture_diffusivity(view: &PhaseView, index: usize, reference: f64, gas_phase: bool) -> f64 {
    if view.z.len() <= 1 {
        return reference;
    }
    let mut resistance = 0.0;
    for other in 0..view.z.len() {
        if other == index {
            continue;
        }
        let diffusivity = binary_diffusivity(view, index, other, gas_phase);
        resistance += view.z[other].max(0.0) / diffusivity;
    }
    if resistance > 0.0 && resistance.is_finite() {
        1.0 / resistance
    } else {
        reference
    }
}

/// The diagonal of a matrix's inverse, through the workspace's one LU.
///
/// **`solve_lu` answers `A x = b`, and the class asks for `A⁻¹` and reads its diagonal** -
/// which is `n` unit right-hand sides, the route JAMA's own `inverse` takes
/// (`LUDecomposition.solve(identity)`). Keeping one LU in the workspace and wrapping it is
/// the whole reason this is not a new solver.
///
/// **The refusal is not an error here.** `solve_lu` refuses a zero pivot, and the class
/// *catches* its own `SingularMatrixException` and takes a different fallback from the one it
/// takes for a non-finite answer - so the two must stay apart, and this returns which
/// happened rather than propagating either.
enum Inverted {
    /// The inverse's diagonal.
    Diagonal(Vec<f64>),
    /// A pivot was zero: the class's `SingularMatrixException` branch.
    Singular,
}

fn inverse_diagonal(a: &[Vec<f64>]) -> Inverted {
    let n = a.len();
    let mut diagonal = Vec::with_capacity(n);
    for column in 0..n {
        let mut unit = vec![0.0; n];
        unit[column] = 1.0;
        match solve_lu(a, &unit) {
            Ok(solution) => diagonal.push(solution[column]),
            Err(AzothError::InvalidInput { .. }) => return Inverted::Singular,
            Err(_) => return Inverted::Singular,
        }
    }
    Inverted::Diagonal(diagonal)
}

/// `maxwellStefanFilmCoefficient`.
///
/// The resistance matrix is `(n − 1)²` over the components' binary film coefficients, and
/// only its diagonal is read.
pub fn maxwell_stefan_film_coefficient(
    view: &PhaseView,
    component_index: Option<usize>,
    base: f64,
    reference: f64,
    gas_phase: bool,
) -> f64 {
    if !(base > 0.0 && base.is_finite()) {
        return 0.0;
    }
    let Some(component_index) = component_index else {
        // **The class's `componentIndex(phase, component) < 0` branch**, and it is not a
        // refusal: a component this phase does not carry has no row of the matrix, and
        // `mixtureDiffusivityForComponent` answers the *reference* diffusivity for it - so
        // `scaleFilmCoefficient` returns the base coefficient unscaled. Reaching for another
        // component's row instead is what `unwrap_or(0)` did here, and component zero's
        // diffusivity is not this component's.
        return scale_film_coefficient(base, reference, reference);
    };
    let count = view.z.len();
    let reduced = count.saturating_sub(1);
    // **Two states take the scalar route, and both are the class's own.** A phase of two
    // components has a one-by-one matrix whose inverse is the same scalar, and the *last*
    // component has no row of its own - so the class scales it by its mixture diffusivity
    // instead. Measured, the class's own absorber states are two-component, which is why the
    // matrix below is not on their path at all.
    if component_index >= reduced {
        return scale_film_coefficient(
            base,
            mixture_diffusivity(view, component_index, reference, gas_phase),
            reference,
        );
    }

    let mut matrix = vec![vec![0.0; reduced]; reduced];
    for (row, matrix_row) in matrix.iter_mut().enumerate() {
        let mut row_sum = 0.0;
        let reference_coefficient =
            binary_film_coefficient(view, row, reduced, base, reference, gas_phase);
        for (column, slot) in matrix_row.iter_mut().enumerate() {
            let binary = binary_film_coefficient(view, row, column, base, reference, gas_phase);
            if row != column {
                row_sum += view.z[column].max(0.0) / binary;
            }
            *slot = -view.z[row].max(0.0) * (1.0 / binary - 1.0 / reference_coefficient);
        }
        // **The last component's column has no entry of its own** - the matrix is `(n-1)²` -
        // so it contributes to the row sum and to nothing else. The class's own loop runs to
        // `componentCount` for exactly this reason, and its `column < reduced` guard is the
        // write above.
        let last = binary_film_coefficient(view, row, reduced, base, reference, gas_phase);
        if row != reduced {
            row_sum += view.z[reduced].max(0.0) / last;
        }
        matrix_row[row] += row_sum + view.z[row].max(0.0) / reference_coefficient;
    }

    match inverse_diagonal(&matrix) {
        Inverted::Diagonal(diagonal) => match diagonal.get(component_index).copied() {
            Some(coefficient) if coefficient > 0.0 && coefficient.is_finite() => {
                clamp(coefficient, base * 0.02, base * 50.0)
            }
            // The inverse answered, and the entry is not a positive number: the class falls
            // out of the `if` and returns the base coefficient unchanged.
            _ => base,
        },
        Inverted::Singular => scale_film_coefficient(
            base,
            mixture_diffusivity(view, component_index, reference, gas_phase),
            reference,
        ),
    }
}

/// `calculateFilmCoefficient`: the base coefficient under the scalar model, and the matrix's
/// correction under the Maxwell-Stefan one.
///
/// `component_index` is the component's position **in this phase's own system**, or `None`
/// where that phase does not carry it - the class's own `componentIndex(phase, component)`,
/// which it looks up once per phase rather than once per segment.
pub fn film_coefficient(
    view: &PhaseView,
    component_index: Option<usize>,
    base: f64,
    reference: f64,
    gas_phase: bool,
    matrix_model: bool,
) -> f64 {
    if !matrix_model {
        return base;
    }
    maxwell_stefan_film_coefficient(view, component_index, base, reference, gas_phase)
}

/// `finitePositive`: the value where it is positive and finite, the fallback otherwise.
pub fn finite_positive(value: f64, fallback: f64) -> f64 {
    if value > 0.0 && value.is_finite() {
        value
    } else {
        fallback
    }
}

/// `clamp`.
pub fn clamp(value: f64, min: f64, max: f64) -> f64 {
    value.max(min).min(max)
}

#[cfg(test)]
mod tests {
    use azoth_core::units::{kelvins, pascal_seconds, pascals, watts_per_meter_kelvin};
    use azoth_eos::phase_transport::PhaseKind;
    use azoth_eos::results::PhaseTransportResult;

    use super::*;

    /// A three-component phase whose components 0 and 1 have different binary diffusivities,
    /// so a film coefficient read from component 0's row is visible against one read from
    /// component 1's, and against the reference.
    fn a_three_component_phase() -> PhaseView {
        let d = |first: usize, second: usize| match (first.min(second), first.max(second)) {
            (0, 1) => 1.0e-5,
            (0, 2) => 2.0e-5,
            (1, 2) => 4.0e-5,
            _ => 0.0,
        };
        PhaseView {
            z: vec![0.6, 0.3, 0.1],
            n: 1.0,
            t: kelvins(313.15),
            p: pascals(5.0e6),
            molar_mass: 0.02,
            density: 40.0,
            cp_mass: 2100.0,
            kind: PhaseKind::Gas,
            transport: PhaseTransportResult {
                mu: pascal_seconds(1.0e-5),
                k: watts_per_meter_kelvin(0.02),
                d_binary: (0..3)
                    .map(|row| (0..3).map(|column| d(row, column)).collect())
                    .collect(),
                d_effective: vec![1.0e-5, 2.0e-5, 3.0e-5],
                warnings: Vec::new(),
            },
        }
    }

    /// **A component the phase does not carry is scaled by nothing, not by component zero.**
    ///
    /// `componentIndex(phase, component)` is `-1` for a component the phase's own system does
    /// not have, and the class's `mixtureDiffusivityForComponent` answers the *reference*
    /// diffusivity for it - so the film coefficient comes back as the base coefficient. The
    /// port read component 0's row instead, which is a different number and a plausible one.
    #[test]
    fn a_component_the_phase_does_not_carry_is_the_base_coefficient() {
        let view = a_three_component_phase();
        let base = 100.0;
        let reference = 2.0e-5;

        assert_eq!(
            film_coefficient(&view, None, base, reference, true, true),
            base
        );

        // And the row that `unwrap_or(0)` reached for is not that number, which is what makes
        // the `None` a value rather than a spelling of zero.
        let component_zero = film_coefficient(&view, Some(0), base, reference, true, true);
        assert!(
            (component_zero - base).abs() > 1.0,
            "component zero's row answers {component_zero}, so a fabricated index is visible"
        );
    }
}
