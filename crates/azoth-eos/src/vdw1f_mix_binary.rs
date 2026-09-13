//! `eos.vdw1f_mix_binary` - van der Waals one-fluid mixing, for a binary.
//!
//! ```text
//! a_mix = z1**2*a1 + 2*z1*z2*(1 - k12)*sqrt(a1*a2) + z2**2*a2
//! b_mix = z1*b1 + z2*b2
//! ```
//!
//! Spec: `specs/calcs/eos/vdw1f_mix_binary.yaml`
//!
//! # Why two components and not N
//!
//! The quadratic form for `a_mix` is a double sum, so a general version takes a
//! composition vector and a matrix of `kij`. The calc registry is scalar - its
//! inputs are named quantities with units, and there is no vector or matrix type -
//! so two components is the largest it can express. The general case belongs to the
//! model layer, and that layer's spec requires it to reduce to this calc at N = 2.
//!
//! # The two terms a reader should check
//!
//! Written longhand rather than as a loop, so the double sum's three distinct parts
//! are visible: the two pure terms, and the cross term carrying `k12`. An
//! implementation that got the cross term's factor of two wrong, or that transposed
//! `z1` and `z2` in it, would still look like a weighted average at a mid-range
//! composition - which is why the pure-component limits are a spec test.

use azoth_core::{Result, apply_checks};

use crate::results::Vdw1fMixBinaryResult;
use crate::spec_gen;

/// The van der Waals one-fluid mixture parameters for a binary.
///
/// `z1` is the mole fraction of component 1 and component 2 is `1 - z1`; there is no
/// `z2` input, because passing both would allow a pair that does not sum to one.
///
/// `k12` is the binary interaction parameter, fitted per pair and supplied by the
/// caller - this library ships no values for it.
///
/// # Errors
/// * [`azoth_core::AzothError::OutOfRange`] if `z1` is outside `[0, 1]`, or if the
///   resulting `a_mix` is negative - which happens when `k12` is outside `[0, 2]`
///   and the composition is unfavourable, and which `eos.pr_z_factor` would refuse
///   downstream anyway.
///
/// # Example
/// ```
/// use azoth_eos::vdw1f_mix_binary;
///
/// let m = vdw1f_mix_binary(0.6, 0.20206500174625697, 0.08448417260831159,
///                          0.02431127309496514, 0.025932024634629486, 0.05)?;
/// assert!((m.a_mix - 0.14584053499701302).abs() < 1e-15);
/// assert!((m.b_mix - 0.02495957371083088).abs() < 1e-15);
/// # Ok::<(), azoth_core::AzothError>(())
/// ```
pub fn vdw1f_mix_binary(
    z1: f64,
    a1: f64,
    a2: f64,
    b1: f64,
    b2: f64,
    k12: f64,
) -> Result<Vdw1fMixBinaryResult> {
    let spec = &spec_gen::VDW1F_MIX_BINARY_SPEC;
    let mut warnings = Vec::new();

    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "z1" => Some(z1),
            "a1" => Some(a1),
            "a2" => Some(a2),
            "b1" => Some(b1),
            "b2" => Some(b2),
            "k12" => Some(k12),
            _ => None,
        },
        &mut warnings,
    )?;

    let z2 = 1.0 - z1;
    let a_mix = z1 * z1 * a1 + 2.0 * z1 * z2 * (1.0 - k12) * (a1 * a2).sqrt() + z2 * z2 * a2;
    let b_mix = z1 * b1 + z2 * b2;

    apply_checks(
        spec.derived_checks(),
        |quantity| match quantity {
            "a_mix" => Some(a_mix),
            _ => None,
        },
        &mut warnings,
    )?;

    Ok(Vdw1fMixBinaryResult {
        a_mix,
        b_mix,
        warnings,
    })
}
