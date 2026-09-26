//! `unit_ops.rate_based_packed_column` - the segment model's kernel.
//!
//! The arithmetic is [`crate::segment`]'s; what this module adds is the class's own
//! `validateSetup`, the five enum families with the values that are carried and the values that
//! are refused by name, and the shape a flowsheet's four streams take.

use azoth_core::{AzothError, Result};

use crate::segment::{
    Fallbacks, ProfileSettings, SegmentResult, SnapshotSettings, solve_fixed_point_profile,
};
use crate::stream::Stream;

/// `MassTransferCorrelation`: which correlation the wetted area and the film coefficients come
/// from.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MassTransferCorrelation {
    /// `ONDA_1968`, the class's default and the one every state of its own tests runs.
    Onda1968,
    /// `BILLET_SCHULTES_1999`, **not ported**.
    BilletSchultes1999,
}

/// `FilmModel`: the scalar film pair or the Maxwell-Stefan matrix.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FilmModel {
    /// `OVERALL_TWO_RESISTANCE`, the scalar path the matrix one falls back to.
    OverallTwoResistance,
    /// `MAXWELL_STEFAN_MATRIX`, the class's default.
    MaxwellStefanMatrix,
}

/// `HeatTransferModel`: whether the interphase heat step runs at all.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum HeatTransferModel {
    /// `NONE`: the outlets are re-equilibrated after material transfer alone, and both
    /// coefficients are exactly zero.
    None,
    /// `CHILTON_COLBURN_ANALOGY`, the class's default.
    ChiltonColburnAnalogy,
}

/// `SegmentSolver`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SegmentSolver {
    /// `SEQUENTIAL_EXPLICIT`, the class's default and the only one this port carries.
    SequentialExplicit,
    /// `SIMULTANEOUS_RESIDUAL`, **not ported** - and **NeqSim disables its own test** for it.
    SimultaneousResidual,
}

/// `ColumnSolver`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ColumnSolver {
    /// `FIXED_POINT_PROFILE`, the class's default.
    FixedPointProfile,
    /// `EQUATION_ORIENTED`, **not ported**: the column-wide damped Newton with homotopy
    /// continuation.
    EquationOriented,
}

/// A value the class carries but this port does not, named with the class behind it.
fn refused(parameter: &str, value: &str, class: &str, detail: &str) -> AzothError {
    AzothError::invalid_input(
        parameter,
        format!(
            "`{value}` is not ported: `RateBasedPackedColumn.{class}` is the class that would \
             close it. {detail}"
        ),
    )
}

impl MassTransferCorrelation {
    /// Parse the spec's spelling.
    ///
    /// # Errors
    /// [`AzothError::InvalidInput`] for `billet_schultes_1999`, naming the class.
    pub fn parse(value: &str) -> Result<Self> {
        match value {
            "onda_1968" => Ok(Self::Onda1968),
            "billet_schultes_1999" => Err(refused(
                "mass_transfer_correlation",
                value,
                "MassTransferCorrelation.BILLET_SCHULTES_1999",
                "Measured, that value is not a correlation at all: it is a constant multiplier, \
                 `max(0.1, Ch/0.4)` on `kGa` and `max(0.1, Cp)` on `kLa`.",
            )),
            other => Err(AzothError::invalid_input(
                "mass_transfer_correlation",
                format!("`{other}` is not one of onda_1968 or billet_schultes_1999"),
            )),
        }
    }
}

impl FilmModel {
    /// Parse the spec's spelling.
    ///
    /// # Errors
    /// [`AzothError::InvalidInput`] for a value the class does not carry.
    pub fn parse(value: &str) -> Result<Self> {
        match value {
            "overall_two_resistance" => Ok(Self::OverallTwoResistance),
            "maxwell_stefan_matrix" => Ok(Self::MaxwellStefanMatrix),
            other => Err(AzothError::invalid_input(
                "film_model",
                format!("`{other}` is not one of overall_two_resistance or maxwell_stefan_matrix"),
            )),
        }
    }
}

impl HeatTransferModel {
    /// Parse the spec's spelling.
    ///
    /// # Errors
    /// [`AzothError::InvalidInput`] for a value the class does not carry.
    pub fn parse(value: &str) -> Result<Self> {
        match value {
            "none" => Ok(Self::None),
            "chilton_colburn_analogy" => Ok(Self::ChiltonColburnAnalogy),
            other => Err(AzothError::invalid_input(
                "heat_transfer_model",
                format!("`{other}` is not one of none or chilton_colburn_analogy"),
            )),
        }
    }
}

impl SegmentSolver {
    /// Parse the spec's spelling.
    ///
    /// # Errors
    /// [`AzothError::InvalidInput`] for `simultaneous_residual`, naming the class.
    pub fn parse(value: &str) -> Result<Self> {
        match value {
            "sequential_explicit" => Ok(Self::SequentialExplicit),
            "simultaneous_residual" => Err(refused(
                "segment_solver",
                value,
                "SegmentSolver.SIMULTANEOUS_RESIDUAL",
                "It solves the flux residuals and the interfacial heat balance together, and \
                 **NeqSim disables its own test** for the branch: `@Disabled(\"TODO: not working \
                 per 19.06.2060\")`.",
            )),
            other => Err(AzothError::invalid_input(
                "segment_solver",
                format!("`{other}` is not one of sequential_explicit or simultaneous_residual"),
            )),
        }
    }
}

impl ColumnSolver {
    /// Parse the spec's spelling.
    ///
    /// # Errors
    /// [`AzothError::InvalidInput`] for `equation_oriented`, naming the class.
    pub fn parse(value: &str) -> Result<Self> {
        match value {
            "fixed_point_profile" => Ok(Self::FixedPointProfile),
            "equation_oriented" => Err(refused(
                "column_solver",
                value,
                "ColumnSolver.EQUATION_ORIENTED",
                "It is a column-wide damped Newton with homotopy continuation \
                 (`solveEquationOrientedProfile`).",
            )),
            other => Err(AzothError::invalid_input(
                "column_solver",
                format!("`{other}` is not one of fixed_point_profile or equation_oriented"),
            )),
        }
    }
}

/// What the column's solve hands back.
#[derive(Debug, Clone)]
pub struct RateBasedOutcome {
    /// The gas leaving the top segment.
    pub gas_out: Stream,
    /// The liquid leaving the bottom segment.
    pub liquid_out: Stream,
    /// One record per segment, bottom first.
    pub segments: Vec<SegmentResult>,
    /// The passes taken.
    pub iterations: usize,
    /// The outlet residual at the last pass, mol/s.
    pub convergence_residual: f64,
    /// Whether the gate was met.
    pub converged: bool,
    /// The sum of every segment's transfers, magnitudes added, mol/s.
    pub total_absolute_molar_transfer: f64,
    /// Each component's net total, in the order the segments produced them.
    pub component_transfer_totals: Vec<(String, f64)>,
    /// Which of the class's constants stood in for a missing property.
    pub fallbacks: Fallbacks,
}

/// Everything the solve takes.
#[derive(Debug, Clone)]
pub struct RateBasedSetup {
    /// The gas, which enters the bottom segment.
    pub gas: Stream,
    /// The liquid, which enters the top segment.
    pub liquid: Stream,
    pub column_diameter: f64,
    pub packed_height: f64,
    pub number_of_segments: usize,
    pub packing_type: String,
    pub max_iterations: usize,
    pub convergence_tolerance: f64,
    pub mass_transfer_correction: f64,
    pub heat_transfer_correction: f64,
    /// The transfer whitelist; `None` is the union of the two inlets' components.
    pub transfer_components: Option<Vec<String>>,
    pub film_model: FilmModel,
    pub heat_transfer_model: HeatTransferModel,
}

/// Solve a rate-based packed column.
///
/// # Errors
/// [`AzothError::InvalidInput`] for the five conditions `validateSetup` states - a non-positive
/// diameter, a negative height, no segments, a non-positive iteration cap and a non-positive
/// convergence gate - and whatever the profile solve raises.
pub fn rate_based_packed_column(setup: &RateBasedSetup) -> Result<RateBasedOutcome> {
    // `validateSetup`'s own five, in its own order.
    if !(setup.column_diameter > 0.0 && setup.column_diameter.is_finite()) {
        return Err(AzothError::invalid_input(
            "column_diameter",
            format!(
                "a column diameter of {} is not positive, which `validateSetup` refuses with \
                 \"Column diameter must be positive\"",
                setup.column_diameter
            ),
        ));
    }
    if setup.packed_height < 0.0 || !setup.packed_height.is_finite() {
        return Err(AzothError::invalid_input(
            "packed_height",
            format!(
                "a packed height of {} is negative, which `validateSetup` refuses with \"Packed \
                 height can not be negative\". **Zero is accepted** and is a state that \
                 transfers nothing",
                setup.packed_height
            ),
        ));
    }
    if setup.number_of_segments < 1 {
        return Err(AzothError::invalid_input(
            "number_of_segments",
            "a column of no segments has no profile, which `validateSetup` refuses with \"At \
             least one segment is required\"",
        ));
    }
    if setup.max_iterations < 1 {
        return Err(AzothError::invalid_input(
            "max_iterations",
            "the profile loop runs at least once, and a cap of zero has no last iterate to \
             publish",
        ));
    }
    if !(setup.convergence_tolerance > 0.0 && setup.convergence_tolerance.is_finite()) {
        return Err(AzothError::invalid_input(
            "convergence_tolerance",
            format!(
                "a convergence gate of {} is not positive, and zero is a solve that never stops",
                setup.convergence_tolerance
            ),
        ));
    }

    let outcome = solve_fixed_point_profile(
        &setup.gas,
        &setup.liquid,
        &ProfileSettings {
            packed_height: setup.packed_height,
            number_of_segments: setup.number_of_segments,
            tolerance: setup.convergence_tolerance,
            max_iterations: setup.max_iterations,
        },
        &SnapshotSettings {
            column_diameter: setup.column_diameter,
            packed_height: setup.packed_height,
            packing: setup.packing_type.clone(),
            heat_transfer_none: setup.heat_transfer_model == HeatTransferModel::None,
            mass_transfer_correction: setup.mass_transfer_correction,
            heat_transfer_correction: setup.heat_transfer_correction,
        },
        setup.transfer_components.as_deref().unwrap_or(&[]),
        setup.film_model == FilmModel::MaxwellStefanMatrix,
    )?;

    Ok(RateBasedOutcome {
        gas_out: outcome.gas_outlet,
        liquid_out: outcome.liquid_outlet,
        segments: outcome.segments,
        iterations: outcome.iterations,
        convergence_residual: outcome.convergence_residual,
        converged: outcome.converged,
        total_absolute_molar_transfer: outcome.total_absolute_molar_transfer,
        component_transfer_totals: outcome.component_transfer_totals,
        fallbacks: outcome.fallbacks,
    })
}
