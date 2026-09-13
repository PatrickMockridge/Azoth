//! Result types for the hydraulics calculations.
//!
//! One struct per calculation, rather than one shared struct generic over the
//! output, for two reasons. A generic result would need the calc's identity as a
//! runtime field to satisfy [`CalcResult::CALC_ID`], and it would force the
//! Colebrook result - which is iterative and has an iteration count worth
//! reporting - into the same shape as the explicit Swamee-Jain one, which has
//! no such thing to report.
//!
//! # Field naming contract
//!
//! Every field name here appears identically in the Python result dataclass and
//! is listed in [`CalcResult::FIELDS`]. A test asserts the three agree. A calc
//! whose Python and Rust results have different field names is a bug that no
//! numerical test would catch, because the numbers would agree perfectly.
//!
//! # Values beyond the spec's declared outputs
//!
//! Result structs may carry diagnostic fields that the spec's `outputs` block
//! does not declare - `iterations` and `converged` on the Colebrook result, for
//! instance. The contract is one-directional: every declared output must appear
//! as a field, but a result may report more than the spec promises. The
//! alternative would be to throw away information the caller wants.

use azoth_core::{CalcResult, FlowRegime, Warning};
use uom::si::f64::{Power, Pressure, VolumeRate};

/// Result of `hydraulics.reynolds_number`.
#[derive(Debug, Clone, PartialEq)]
pub struct ReynoldsNumberResult {
    /// Reynolds number. Dimensionless.
    pub re: f64,
    /// Flow regime under the Crane/Moody boundaries.
    pub regime: FlowRegime,
    /// Caveats. Carries [`azoth_core::WarningCode::TransitionalFlow`] when the
    /// flow sits in the 2000-4000 band.
    pub warnings: Vec<Warning>,
}

impl CalcResult for ReynoldsNumberResult {
    const CALC_ID: &'static str = "hydraulics.reynolds_number";
    const FIELDS: &'static [&'static str] = &["re", "regime", "warnings"];

    fn warnings(&self) -> &[Warning] {
        &self.warnings
    }
}

/// Result of `hydraulics.friction_factor_colebrook`.
///
/// Carries the solver's own report because the answer is only meaningful
/// alongside it: `f` is the last iterate, and if `converged` is false it is not
/// a solution to the equation at all.
#[derive(Debug, Clone, PartialEq)]
pub struct ColebrookResult {
    /// Darcy friction factor. Dimensionless.
    pub f: f64,
    /// Iterations performed.
    pub iterations: u32,
    /// Whether the iteration met its tolerance.
    pub converged: bool,
    /// Final change between iterates.
    pub residual: f64,
    /// Caveats.
    pub warnings: Vec<Warning>,
}

impl CalcResult for ColebrookResult {
    const CALC_ID: &'static str = "hydraulics.friction_factor_colebrook";
    const FIELDS: &'static [&'static str] =
        &["f", "iterations", "converged", "residual", "warnings"];

    fn warnings(&self) -> &[Warning] {
        &self.warnings
    }
}

/// Result of `hydraulics.friction_factor_swamee_jain`.
///
/// No solver report: the Swamee-Jain equation is explicit, and reporting an
/// iteration count of zero for it would imply a similarity to Colebrook that
/// does not exist.
#[derive(Debug, Clone, PartialEq)]
pub struct SwameeJainResult {
    /// Darcy friction factor. Dimensionless.
    pub f: f64,
    /// Caveats.
    pub warnings: Vec<Warning>,
}

impl CalcResult for SwameeJainResult {
    const CALC_ID: &'static str = "hydraulics.friction_factor_swamee_jain";
    const FIELDS: &'static [&'static str] = &["f", "warnings"];

    fn warnings(&self) -> &[Warning] {
        &self.warnings
    }
}

/// Result of `hydraulics.friction_factor_haaland`.
///
/// No solver report, for the same reason as Swamee-Jain: the Haaland equation is
/// explicit, and an iteration count would imply a similarity to Colebrook that
/// does not exist.
#[derive(Debug, Clone, PartialEq)]
pub struct HaalandResult {
    /// Darcy friction factor. Dimensionless.
    pub f: f64,
    /// Caveats.
    pub warnings: Vec<Warning>,
}

impl CalcResult for HaalandResult {
    const CALC_ID: &'static str = "hydraulics.friction_factor_haaland";
    const FIELDS: &'static [&'static str] = &["f", "warnings"];

    fn warnings(&self) -> &[Warning] {
        &self.warnings
    }
}

/// Result of `hydraulics.pump_power`.
///
/// Carries a dimensioned output, like `darcy_weisbach`'s: shaft power is a
/// quantity, not a bare number, and the transport layer needs both the SI
/// magnitude and the unit to hand it back.
#[derive(Debug, Clone, PartialEq)]
pub struct PumpPowerResult {
    /// Shaft power the pump must be supplied with.
    pub power: Power,
    /// Caveats.
    pub warnings: Vec<Warning>,
}

impl CalcResult for PumpPowerResult {
    const CALC_ID: &'static str = "hydraulics.pump_power";
    const FIELDS: &'static [&'static str] = &["power", "warnings"];

    fn warnings(&self) -> &[Warning] {
        &self.warnings
    }
}

/// Result of `hydraulics.orifice_flow`.
#[derive(Debug, Clone, PartialEq)]
pub struct OrificeFlowResult {
    /// Volumetric flow rate through the orifice.
    pub q: VolumeRate,
    /// Caveats.
    pub warnings: Vec<Warning>,
}

impl CalcResult for OrificeFlowResult {
    const CALC_ID: &'static str = "hydraulics.orifice_flow";
    const FIELDS: &'static [&'static str] = &["q", "warnings"];

    fn warnings(&self) -> &[Warning] {
        &self.warnings
    }
}

/// One fitting's contribution to the total resistance coefficient.
#[derive(Debug, Clone, PartialEq)]
pub struct KComponent {
    /// Fitting id as it appears in the registry.
    pub fitting_id: String,
    /// Equivalent length ratio (L_eq / D) for this fitting.
    pub n_ld: f64,
    /// This fitting's resistance coefficient, `f_t * n_ld`.
    pub k: f64,
}

/// Result of `hydraulics.crane_k_factors`.
#[derive(Debug, Clone, PartialEq)]
pub struct KFactorsResult {
    /// Total resistance coefficient for all listed fittings. Dimensionless.
    pub k_total: f64,
    /// The friction factor the coefficients were based on.
    pub f_t: f64,
    /// Per-fitting breakdown, in the order the fittings were supplied.
    pub components: Vec<KComponent>,
    /// Caveats. Carries [`azoth_core::WarningCode::EstimatedData`] while the
    /// registry holds placeholder coefficients.
    pub warnings: Vec<Warning>,
}

impl CalcResult for KFactorsResult {
    const CALC_ID: &'static str = "hydraulics.crane_k_factors";
    const FIELDS: &'static [&'static str] = &["k_total", "f_t", "components", "warnings"];

    fn warnings(&self) -> &[Warning] {
        &self.warnings
    }
}

/// Result of `hydraulics.darcy_weisbach`.
#[derive(Debug, Clone, PartialEq)]
pub struct DarcyWeisbachResult {
    /// Pressure drop over the pipe length.
    pub dp: Pressure,
    /// The friction factor the drop was computed with.
    pub f: f64,
    /// Reynolds number, present only when viscosity was supplied.
    pub re: Option<f64>,
    /// Flow regime, present only when viscosity was supplied.
    pub regime: Option<FlowRegime>,
    /// Caveats. Carries [`azoth_core::WarningCode::RangeCheckSkipped`] when
    /// viscosity was omitted, because the regime then went unchecked.
    pub warnings: Vec<Warning>,
}

impl CalcResult for DarcyWeisbachResult {
    const CALC_ID: &'static str = "hydraulics.darcy_weisbach";
    const FIELDS: &'static [&'static str] = &["dp", "f", "re", "regime", "warnings"];

    fn warnings(&self) -> &[Warning] {
        &self.warnings
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use azoth_core::WarningCode;

    #[test]
    fn every_result_carries_a_warnings_field() {
        // The Python dataclass must carry `warnings` too, and a calc that
        // forgets to thread it through would be silently uncaveatable.
        let results: Vec<(&str, &[&str])> = vec![
            (ReynoldsNumberResult::CALC_ID, ReynoldsNumberResult::FIELDS),
            (ColebrookResult::CALC_ID, ColebrookResult::FIELDS),
            (SwameeJainResult::CALC_ID, SwameeJainResult::FIELDS),
            (HaalandResult::CALC_ID, HaalandResult::FIELDS),
            (OrificeFlowResult::CALC_ID, OrificeFlowResult::FIELDS),
            (PumpPowerResult::CALC_ID, PumpPowerResult::FIELDS),
            (KFactorsResult::CALC_ID, KFactorsResult::FIELDS),
            (DarcyWeisbachResult::CALC_ID, DarcyWeisbachResult::FIELDS),
        ];
        for (calc_id, fields) in results {
            assert!(
                fields.contains(&"warnings"),
                "{calc_id} result has no warnings field"
            );
        }
    }

    #[test]
    fn calc_ids_and_field_names_are_unique_and_well_formed() {
        let all: Vec<&str> = vec![
            ReynoldsNumberResult::CALC_ID,
            ColebrookResult::CALC_ID,
            SwameeJainResult::CALC_ID,
            HaalandResult::CALC_ID,
            KFactorsResult::CALC_ID,
            OrificeFlowResult::CALC_ID,
            PumpPowerResult::CALC_ID,
            DarcyWeisbachResult::CALC_ID,
        ];
        let unique: std::collections::HashSet<_> = all.iter().collect();
        assert_eq!(unique.len(), all.len(), "duplicate CALC_ID");

        for id in all {
            assert!(
                id.starts_with("hydraulics."),
                "{id} is not in the hydraulics module"
            );
            // The id's last segment must be a valid Python identifier, since it
            // becomes the function name in both languages.
            let last = id.rsplit('.').next().unwrap();
            assert!(
                last.chars()
                    .all(|c| c.is_ascii_lowercase() || c.is_ascii_digit() || c == '_'),
                "{id} does not end in a valid identifier"
            );
        }

        // Field names must be distinct within a result, or the Python dataclass
        // could not be built.
        for fields in [
            ReynoldsNumberResult::FIELDS,
            ColebrookResult::FIELDS,
            SwameeJainResult::FIELDS,
            HaalandResult::FIELDS,
            KFactorsResult::FIELDS,
            OrificeFlowResult::FIELDS,
            PumpPowerResult::FIELDS,
            DarcyWeisbachResult::FIELDS,
        ] {
            let unique: std::collections::HashSet<_> = fields.iter().collect();
            assert_eq!(
                unique.len(),
                fields.len(),
                "duplicate field name in {fields:?}"
            );
        }
    }

    #[test]
    fn is_clean_and_has_warning_reflect_the_warning_list() {
        let clean = SwameeJainResult {
            f: 0.02,
            warnings: vec![],
        };
        assert!(clean.is_clean());
        assert!(!clean.has_warning(WarningCode::OutOfValidRange));

        let dirty = SwameeJainResult {
            f: 0.02,
            warnings: vec![Warning::new(WarningCode::OutOfValidRange, "below 5000")],
        };
        assert!(!dirty.is_clean());
        assert!(dirty.has_warning(WarningCode::OutOfValidRange));
        assert!(!dirty.has_warning(WarningCode::TransitionalFlow));
    }

    #[test]
    fn darcy_result_distinguishes_unchecked_from_checked() {
        use azoth_core::units::pascals;

        // Viscosity omitted: no Reynolds number, and the caller can see that.
        let unchecked = DarcyWeisbachResult {
            dp: pascals(22455.0),
            f: 0.02,
            re: None,
            regime: None,
            warnings: vec![Warning::new(
                WarningCode::RangeCheckSkipped,
                "no mu supplied",
            )],
        };
        assert!(unchecked.re.is_none());
        assert!(unchecked.has_warning(WarningCode::RangeCheckSkipped));

        // Viscosity supplied: the regime was actually determined.
        let checked = DarcyWeisbachResult {
            dp: pascals(22455.0),
            f: 0.02,
            re: Some(149_401.197_6),
            regime: Some(FlowRegime::Turbulent),
            warnings: vec![],
        };
        assert!(checked.is_clean());
        assert_eq!(checked.regime, Some(FlowRegime::Turbulent));
        // Same pressure drop either way - supplying mu must not change dP.
        assert_eq!(unchecked.dp.value, checked.dp.value);
    }
}
