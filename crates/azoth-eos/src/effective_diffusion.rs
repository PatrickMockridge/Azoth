//! `eos.effective_diffusion` - the effective diffusion coefficients of a phase.
//!
//! Spec: `specs/models/eos/effective_diffusion.toml`. A *direct* model: eight lines, no
//! iteration, so no `algorithm` block.
//!
//! ```text
//! D_eff_i = (1 - x_i) / sum_{j != i} x_j / D_ij
//! ```
//!
//! `PhysicalProperties.calcEffectiveDiffusionCoefficients` delegates to the phase's
//! diffusivity model, and every model in the family writes the same loop:
//! `commonphasephysicalproperties/diffusivity/Diffusivity` and
//! `liquidphysicalproperties/diffusivity/Diffusivity` are line for line the same eight.
//!
//! # The input that was recorded as unreachable
//!
//! This port's `reactions.kinetics` takes the effective vector as an input, on the ground
//! that the binary matrix behind it could not be read from outside its class: a probe
//! called `getFickDiffusionCoefficient` - or the *liquid* `getFickBinaryDiffusionCoefficient`
//! - and got a diagonal array.
//!
//! **That was the wrong method, and the matrix was reachable all along.**
//! `DiffusivityInterface.calcDiffusionCoefficients(int, int)` **returns** the matrix, and
//! `PhysicalProperties.diffusivityCalc` is a public field, so no cast is needed either. The
//! diagonal came from the liquid override, which multiplies each entry by a Maxwell-Stefan
//! correction `delta_ij + x_i (d ln phi_j / d n) n_t` - and with the derivative
//! uninitialised that correction is the identity, which leaves the diagonal of the matrix
//! and nothing else. `validation/neqsim/EffectiveDiffusionProbe.java` prints the matrix, the
//! vector the class computes from it, and the vector recomputed from it here, in one run:
//! the last two are equal to the digit.
//!
//! # The matrix is not symmetric
//!
//! `D_ij` is component `i`'s coefficient at infinite dilution in `j`, so the two directions
//! are different numbers - `1.432e-7` against `1.428e-7` on the captured oil. Only row `i`
//! is read for `D_eff_i`, which is what makes an asymmetric matrix the normal case.

use azoth_core::{AzothError, Result, apply_checks};

use crate::results::EffectiveDiffusionResult;

/// The effective diffusion coefficients of a phase.
///
/// `binary_diffusion` is row-major and `binary_diffusion[i][j]` is `D_ij`; `x` is the
/// phase's mole fractions in the same order.
///
/// # Errors
/// * [`AzothError::InvalidInput`] if there are fewer than two components, if the matrix is
///   not square, if it is too small for `x`, or if `x` has a negative entry.
/// * [`AzothError::OutOfRange`] if a pair coefficient a sum divides by is zero - the class
///   does not guard it and the quotient is an infinity, which is not a diffusion
///   coefficient.
///
/// # Example
/// ```
/// use azoth_eos::effective_diffusion::effective_diffusion;
///
/// let r = effective_diffusion(
///     &[vec![1.329_810_109_816_717_7e-7, 1.329_818_260_121_59e-7],
///       vec![1.329_826_410_361_882e-7, 2.056_859_540_349_671_9e-7]],
///     &[1.405_247_510_226_945_8e-5, 0.999_985_947_524_897_8],
/// )?;
/// assert!((r.effective_diffusion[0] - 1.329_818_260_121_59e-7).abs() < 1e-20);
/// # Ok::<(), azoth_core::AzothError>(())
/// ```
pub fn effective_diffusion(
    binary_diffusion: &[Vec<f64>],
    x: &[f64],
) -> Result<EffectiveDiffusionResult> {
    let spec = &crate::model_gen::EFFECTIVE_DIFFUSION_SPEC;
    let mut warnings = Vec::new();

    let count = x.len();
    if count < 2 {
        return Err(AzothError::InvalidInput {
            field: "x".to_string(),
            reason: format!(
                "an effective coefficient needs a second component to diffuse into and \
                 there are {count}"
            ),
        });
    }
    if binary_diffusion.len() != count {
        return Err(AzothError::InvalidInput {
            field: "binary_diffusion".to_string(),
            reason: format!(
                "{} row(s) against {count} component(s)",
                binary_diffusion.len()
            ),
        });
    }
    for (row, values) in binary_diffusion.iter().enumerate() {
        if values.len() != count {
            return Err(AzothError::InvalidInput {
                field: "binary_diffusion".to_string(),
                reason: format!(
                    "row {row} has {} column(s) against {count} component(s)",
                    values.len()
                ),
            });
        }
    }

    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            // The class does not check a negative fraction, and a negative one silently
            // subtracts from a sum, so the spec's bound is what catches it.
            "x" => x.iter().copied().reduce(f64::min),
            _ => None,
        },
        &mut warnings,
    )?;

    let mut effective_diffusion = Vec::with_capacity(count);
    for i in 0..count {
        let mut sum = 0.0_f64;
        for j in 0..count {
            if i == j {
                continue;
            }
            let pair = binary_diffusion[i][j];
            if pair == 0.0 {
                return Err(AzothError::OutOfRange {
                    field: "binary_diffusion".to_string(),
                    value: pair,
                    detail: format!(
                        "D[{i}][{j}] is zero and the sum divides by it; a zero pair \
                         coefficient is not a state this can take a quotient of"
                    ),
                });
            }
            sum += x[j] / pair;
        }
        if sum == 0.0 {
            return Err(AzothError::OutOfRange {
                field: "x".to_string(),
                value: x[i],
                detail: format!(
                    "every component other than {i} carries a zero mole fraction, so there \
                     is nothing for it to diffuse into"
                ),
            });
        }
        effective_diffusion.push((1.0 - x[i]) / sum);
    }

    Ok(EffectiveDiffusionResult {
        effective_diffusion,
        warnings,
    })
}
