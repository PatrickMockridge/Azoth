//! `hydraulics.choked_flow_area` - the throat area a choked gas flow needs.
//!
//! ```text
//! A = m_dot / (sqrt(k * rho0 * P0) * (2 / (k + 1))**((k + 1) / (2 * (k - 1))))
//! ```
//!
//! Spec: `specs/calcs/hydraulics/choked_flow_area.toml`
//!
//! # What this is, and what it deliberately is not
//!
//! This is the isentropic critical-flow relation: the mass flux through a throat when
//! the downstream pressure is low enough that the flow reaches sonic velocity there.
//!
//! It is the physical basis of relief valve sizing, and it is **not** relief valve
//! sizing to a standard. API 520 wraps this relation in de-rating coefficients whose
//! values are tabulated in the standard, and this library does not reproduce tables.
//! A caller sizing a relief valve applies those factors themselves, visibly. That is
//! why the calculation is named for what it computes; see the spec.
//!
//! # Why the units work out
//!
//! The unit structure lives in `sqrt(k rho0 P0)`, which is a mass flux:
//! `sqrt(Pa * kg/m**3)` is `kg/(m**2 s)`. Everything else is a function of `k` alone
//! and is dimensionless.

use azoth_core::units::{MassDensity, MassRate, Pressure, square_meters};
use azoth_core::{Result, apply_checks};

use crate::results::ChokedFlowAreaResult;
use crate::spec_gen;

/// Throat area required for a choked mass flow of an ideal gas. `k` is dimensionless.
///
/// # Errors
/// * [`azoth_core::AzothError::OutOfRange`] if `P0` or `rho0` is not positive, if
///   `m_dot` is negative, or if `k` is not greater than 1.
///
/// Outside the range of real substances - `k` above 5/3, the monatomic limit - the
/// area is still returned, carrying an `OutOfValidRange` warning.
///
/// The flow must be choked for this area to be the right one, and this calc cannot
/// check that: the critical pressure ratio depends on `k`, which the spec's typed
/// range checks cannot express. See the spec's assumptions.
///
/// # Example
/// ```
/// use azoth_core::units::{kilograms_per_cubic_meter, kilograms_per_second, pascals};
/// use azoth_hydraulics::choked_flow_area;
///
/// let r = choked_flow_area(
///     kilograms_per_second(1.0),
///     pascals(1.0e6),
///     kilograms_per_cubic_meter(10.0),
///     1.4,
/// )?;
/// assert!((r.a.value - 0.00046182742602466937).abs() < 1e-18);
/// # Ok::<(), azoth_core::AzothError>(())
/// ```
#[allow(non_snake_case)] // `P0` is the symbol in the published relation
pub fn choked_flow_area(
    m_dot: MassRate,
    P0: Pressure,
    rho0: MassDensity,
    k: f64,
) -> Result<ChokedFlowAreaResult> {
    let spec = &spec_gen::CHOKED_FLOW_AREA_SPEC;
    let mut warnings = Vec::new();

    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "m_dot" => Some(m_dot.value),
            "P0" => Some(P0.value),
            "rho0" => Some(rho0.value),
            // Dimensionless, so it arrives as a plain f64 with no unit to extract.
            "k" => Some(k),
            _ => None,
        },
        &mut warnings,
    )?;

    // Guarded by the checks above, so k > 1 and P0, rho0 > 0.
    //
    // The geometric factor is the derived expression rather than a named constant,
    // because it depends on `k` and so cannot be one. At k = 1.4 it is exactly
    // (5/6)**3 = 125/216, which is what makes the spec's worked example checkable by
    // hand rather than only by trusting this code.
    let geometric = (2.0 / (k + 1.0)).powf((k + 1.0) / (2.0 * (k - 1.0)));
    let flux = (k * rho0.value * P0.value).sqrt() * geometric;
    let a = m_dot.value / flux;

    apply_checks(
        spec.derived_checks(),
        |name| (name == "a").then_some(a),
        &mut warnings,
    )?;

    Ok(ChokedFlowAreaResult {
        a: square_meters(a),
        warnings,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use azoth_core::units::{kilograms_per_cubic_meter, kilograms_per_second, pascals};

    #[test]
    fn an_area_of_one_comes_back_when_the_flow_is_the_critical_flux() {
        // At k = 1.4 the exponent (k+1)/(2(k-1)) is exactly 3, so the geometric factor
        // is exactly (5/6)**3 = 125/216, and with rho0 = P0 = 1 the critical flux is
        // exactly sqrt(1.4) * 125/216. Feeding that in as the mass flow must give
        // exactly one square metre.
        //
        // Stated as "the area must be 1" rather than by recomputing the expression:
        // a test that restates the implementation cannot catch the implementation
        // being wrong, and this one moves off 1 for a wrong exponent, a misplaced
        // root, or a unit that did not convert.
        let k = 1.4_f64;
        let critical_flux = k.sqrt() * 125.0 / 216.0;
        let r = choked_flow_area(
            kilograms_per_second(critical_flux),
            pascals(1.0),
            kilograms_per_cubic_meter(1.0),
            k,
        )
        .expect("should compute");
        assert!((r.a.value - 1.0).abs() < 1e-15, "got {} m**2", r.a.value);
    }

    #[test]
    fn agrees_with_the_stagnation_temperature_form() {
        // The relation is usually written with the stagnation temperature and the gas
        // constant: G* = P0 sqrt(k/(R T0)) (2/(k+1))**((k+1)/(2(k-1))). It is the same
        // physics through different variables - `rho0 = P0/(R T0)` - so the two share
        // no arithmetic beyond the geometric factor.
        //
        // Feeding in that textbook flux must give one square metre here too. This is a
        // check against an independent algebraic *form* rather than against the same
        // expression written twice, which is the difference between a test that can
        // fail and one that cannot.
        let (p0, t0, r_gas, k) = (1.0e6_f64, 300.0_f64, 287.0_f64, 1.4_f64);
        let rho0 = p0 / (r_gas * t0);
        let textbook_flux =
            p0 * (k / (r_gas * t0)).sqrt() * (2.0 / (k + 1.0)).powf((k + 1.0) / (2.0 * (k - 1.0)));

        let r = choked_flow_area(
            kilograms_per_second(textbook_flux),
            pascals(p0),
            kilograms_per_cubic_meter(rho0),
            k,
        )
        .expect("should compute");
        assert!((r.a.value - 1.0).abs() < 1e-9, "got {} m**2", r.a.value);
    }
}
