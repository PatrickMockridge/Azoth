//! Hydraulics calculations.
//!
//! The slice implemented here runs from the Reynolds number through to the
//! Darcy-Weisbach pressure drop for a straight pipe, with fitting losses
//! available separately:
//!
//! * [`reynolds_number`] - the flow regime, and the input every other calc needs
//! * [`friction_factor_colebrook`] - the implicit, accurate friction factor
//! * [`friction_factor_swamee_jain`] - the explicit approximation to it
//! * [`friction_factor_haaland`] - a second explicit approximation, fitted differently
//! * [`crane_k_factors`] - fitting losses by the equivalent-length method
//! * [`darcy_weisbach`] - pressure drop over a straight pipe
//! * [`pump_power`] - shaft power from flow, head and efficiency
//! * [`orifice_flow`] - flow through an orifice from its pressure difference
//! * [`control_valve_cv`] - liquid flow through a control valve
//!
//! Pipe *with* fittings is a composition of the last two, performed by the
//! `azoth pipe` CLI rather than by a calc of its own, because the two losses
//! are computed by different methods and adding them is a modelling decision the
//! caller should be able to see.
//!
//! # Warning before use
//!
//! `crane_k_factors` reads coefficients from `data/fittings/crane_k_factors.csv`,
//! where every row is currently an **estimated dummy value** - a placeholder for
//! software testing, not engineering data. Results built from it carry an
//! [`azoth_core::WarningCode::EstimatedData`] warning. Nothing in this crate
//! should be used for design work until that file is populated from a primary
//! standard.

pub mod control_valve_cv;
pub mod crane_k_factors;
pub mod darcy_weisbach;
pub mod fittings;
pub mod fluids;
pub mod friction_factor_colebrook;
pub mod friction_factor_haaland;
pub mod friction_factor_swamee_jain;
pub mod orifice_flow;
pub mod provenance;
pub mod pump_power;
pub mod results;
pub mod reynolds_number;
pub mod solver;
pub mod spec_gen;

pub use control_valve_cv::{CV_TO_SI, control_valve_cv};
pub use crane_k_factors::{crane_k_factors, known_fittings};
pub use darcy_weisbach::{add_fitting_loss, darcy_weisbach, propagate_estimated_data};
pub use fluids::{available_fluids, provider_for};
pub use friction_factor_colebrook::{friction_factor_colebrook, fully_rough_limit};
pub use friction_factor_haaland::friction_factor_haaland;
pub use friction_factor_swamee_jain::friction_factor_swamee_jain;
pub use orifice_flow::orifice_flow;
pub use provenance::VerifyStatus;
pub use pump_power::{STANDARD_GRAVITY_M_S2, pump_power};
pub use results::{
    ColebrookResult, ControlValveCvResult, DarcyWeisbachResult, HaalandResult, KComponent,
    KFactorsResult, OrificeFlowResult, PumpPowerResult, ReynoldsNumberResult, SwameeJainResult,
};
pub use reynolds_number::{regime_for, regime_warning, reynolds_number};
