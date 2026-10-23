//! Result types for the characterisation calculations.
//!
//! The same contract as the other namespaces' results: one struct per calculation, field names
//! identical to the Python result dataclass and listed in [`CalcResult::FIELDS`], with a test
//! asserting the three agree. See `crates/azoth-core/src/result.rs` for why the duplication is
//! deliberate.

use azoth_core::{CalcResult, Warning};

use crate::characterise_plus_fraction::PlusModel;
use crate::component_row::CharacterizationComponent;
use uom::si::f64::{MassDensity, MolarMass, Pressure, ThermodynamicTemperature};

/// Result of `characterization.characterise_plus_fraction`.
#[derive(Debug, Clone, PartialEq)]
pub struct CharacterisePlusFractionResult {
    /// Which model actually ran, which is not always the one asked for.
    pub selected_model: PlusModel,
    /// Each lump's share of the plus fraction.
    pub fraction_of_heavy_end: Vec<f64>,
    /// Each lump's mass-weighted mean molar mass.
    pub lump_molar_mass: Vec<MolarMass>,
    /// Each lump's mass-weighted harmonic mean density.
    pub lump_density: Vec<MassDensity>,
    /// Caveats.
    pub warnings: Vec<Warning>,
}

impl CalcResult for CharacterisePlusFractionResult {
    const CALC_ID: &'static str = "characterization.characterise_plus_fraction";
    const FIELDS: &'static [&'static str] = &[
        "selected_model",
        "fraction_of_heavy_end",
        "lump_molar_mass",
        "lump_density",
        "warnings",
    ];

    fn warnings(&self) -> &[Warning] {
        &self.warnings
    }
}

/// Result of `characterization.characterize_to_reference`.
#[derive(Debug, Clone, PartialEq)]
pub struct CharacterizeToReferenceResult {
    /// Which reference cut each non-empty group belongs to.
    pub reference_index: Vec<f64>,
    /// Each group's summed mole amount.
    pub group_moles: Vec<f64>,
    /// Each group's molar mass, its accumulated mass over its accumulated moles.
    pub group_molar_mass: Vec<MolarMass>,
    /// Each group's density, the mass over the volume its rows occupy.
    pub group_density: Vec<MassDensity>,
    /// Caveats.
    pub warnings: Vec<Warning>,
}

impl CalcResult for CharacterizeToReferenceResult {
    const CALC_ID: &'static str = "characterization.characterize_to_reference";
    const FIELDS: &'static [&'static str] = &[
        "reference_index",
        "group_moles",
        "group_molar_mass",
        "group_density",
        "warnings",
    ];

    fn warnings(&self) -> &[Warning] {
        &self.warnings
    }
}

/// Result of `characterization.assay_mass_fractions`.
#[derive(Debug, Clone, PartialEq)]
pub struct AssayMassFractionsResult {
    /// Each cut's resolved mass fraction, which sums to one.
    pub mass_fraction: Vec<f64>,
    /// The declared fractions' own sum, before normalisation.
    pub total_declared_fraction: f64,
    /// The assay's bulk density, the mass-weighted harmonic mean, present only when every cut
    /// carries a density.
    pub bulk_density: Option<MassDensity>,
    /// Caveats.
    pub warnings: Vec<Warning>,
}

impl CalcResult for AssayMassFractionsResult {
    const CALC_ID: &'static str = "characterization.assay_mass_fractions";
    const FIELDS: &'static [&'static str] = &[
        "mass_fraction",
        "total_declared_fraction",
        "bulk_density",
        "warnings",
    ];

    fn warnings(&self) -> &[Warning] {
        &self.warnings
    }
}

/// Result of `characterization.tbp_cut_properties`.
#[derive(Debug, Clone, PartialEq)]
pub struct TbpCutPropertiesResult {
    /// Critical temperature.
    pub tc: ThermodynamicTemperature,
    /// Critical pressure. Every model's correlation is written in bar; this is pascals.
    pub pc: Pressure,
    /// Normal boiling point.
    pub boiling_temperature: ThermodynamicTemperature,
    /// Pitzer's acentric factor, by whichever correlation the model uses.
    pub acentric_factor: f64,
    /// The cubic alpha function's `m`, **absent** for the five models that have no exponent to
    /// offer: `RiaziDaubert` sets `calcm = false` and the other four inherit a method that throws.
    pub attraction_exponent: Option<f64>,
    /// Watson's characterization factor, `(1.8*Tb)^(1/3)/d`.
    pub watson_k: f64,
    /// Caveats.
    pub warnings: Vec<Warning>,
}

impl CalcResult for TbpCutPropertiesResult {
    const CALC_ID: &'static str = "characterization.tbp_cut_properties";
    const FIELDS: &'static [&'static str] = &[
        "tc",
        "pc",
        "boiling_temperature",
        "acentric_factor",
        "attraction_exponent",
        "watson_k",
        "warnings",
    ];

    fn warnings(&self) -> &[Warning] {
        &self.warnings
    }
}

/// Result of `characterization.tbp_closure`.
#[derive(Debug, Clone, PartialEq)]
pub struct TbpClosureResult {
    /// The cut's molar mass, the number a cubic needs.
    pub molar_mass: MolarMass,
    /// Caveats.
    pub warnings: Vec<Warning>,
}

impl CalcResult for TbpClosureResult {
    const CALC_ID: &'static str = "characterization.tbp_closure";
    const FIELDS: &'static [&'static str] = &["molar_mass", "warnings"];

    fn warnings(&self) -> &[Warning] {
        &self.warnings
    }
}

/// Result of `characterization.pedersen_plus_split`.
#[derive(Debug, Clone, PartialEq)]
pub struct PedersenPlusSplitResult {
    /// Each cut's mole fraction, `exp(a + b*CN)`.
    pub cut_z: Vec<f64>,
    /// Each cut's molar mass, its own row of the PVTsim table.
    pub cut_molar_mass: Vec<MolarMass>,
    /// Each cut's normal liquid density, `c + d*ln(CN)`.
    pub cut_density: Vec<MassDensity>,
    /// The solved `a` of `z = exp(a + b*CN)`.
    pub z_intercept: f64,
    /// The solved `b` of `z = exp(a + b*CN)`.
    pub z_slope: f64,
    /// The solved `c` of `rho = c + d*ln(CN)`, in NeqSim's g/cm3 scale.
    pub density_intercept: f64,
    /// The solved `d` of `rho = c + d*ln(CN)`, in NeqSim's g/cm3 scale.
    pub density_slope: f64,
    /// Caveats.
    pub warnings: Vec<Warning>,
}

impl CalcResult for PedersenPlusSplitResult {
    const CALC_ID: &'static str = "characterization.pedersen_plus_split";
    const FIELDS: &'static [&'static str] = &[
        "cut_z",
        "cut_molar_mass",
        "cut_density",
        "z_intercept",
        "z_slope",
        "density_intercept",
        "density_slope",
        "warnings",
    ];

    fn warnings(&self) -> &[Warning] {
        &self.warnings
    }
}

/// Result of `characterization.whitson_gamma_split`.
#[derive(Debug, Clone, PartialEq)]
pub struct WhitsonGammaSplitResult {
    /// Each cut's mole fraction, the gamma density integrated across its window.
    pub cut_z: Vec<f64>,
    /// Each cut's molar mass: its window's first gamma moment, or the window's midpoint.
    pub cut_molar_mass: Vec<MolarMass>,
    /// Each cut's normal liquid density from the selected correlation.
    pub cut_density: Vec<MassDensity>,
    /// The shape parameter actually used, given or estimated.
    pub shape: f64,
    /// The minimum molar mass actually used.
    pub minimum_molar_mass: MolarMass,
    /// The derived gamma scale, `(M_plus - eta)/alpha`.
    pub scale: MolarMass,
    /// Caveats.
    pub warnings: Vec<Warning>,
}

impl CalcResult for WhitsonGammaSplitResult {
    const CALC_ID: &'static str = "characterization.whitson_gamma_split";
    const FIELDS: &'static [&'static str] = &[
        "cut_z",
        "cut_molar_mass",
        "cut_density",
        "shape",
        "minimum_molar_mass",
        "scale",
        "warnings",
    ];

    fn warnings(&self) -> &[Warning] {
        &self.warnings
    }
}

/// Result of `characterization.tbp_grouping`.
#[derive(Debug, Clone, PartialEq)]
pub struct TbpGroupingResult {
    /// Each bin's summed mole fraction, over twenty bins of which the last fourteen can fill.
    pub group_fraction: Vec<f64>,
    /// Caveats.
    pub warnings: Vec<Warning>,
}

impl CalcResult for TbpGroupingResult {
    const CALC_ID: &'static str = "characterization.tbp_grouping";
    const FIELDS: &'static [&'static str] = &["group_fraction", "warnings"];

    fn warnings(&self) -> &[Warning] {
        &self.warnings
    }
}

/// Result of `characterization.lumping`.
#[derive(Debug, Clone, PartialEq)]
pub struct LumpingResult {
    /// Each lump's share of the plus fraction, which sums to one only as far as the split does.
    pub fraction_of_heavy_end: Vec<f64>,
    /// Each lump's accumulated mole fraction, which sums to the cut table's own total.
    pub lump_mole_fraction: Vec<f64>,
    /// Each lump's mass-weighted mean molar mass.
    pub lump_molar_mass: Vec<MolarMass>,
    /// Each lump's mass-weighted harmonic mean density.
    pub lump_density: Vec<MassDensity>,
    /// Caveats.
    pub warnings: Vec<Warning>,
}

impl CalcResult for LumpingResult {
    const CALC_ID: &'static str = "characterization.lumping";
    const FIELDS: &'static [&'static str] = &[
        "fraction_of_heavy_end",
        "lump_mole_fraction",
        "lump_molar_mass",
        "lump_density",
        "warnings",
    ];

    fn warnings(&self) -> &[Warning] {
        &self.warnings
    }
}

/// Result of `characterization.tbp_density`.
#[derive(Debug, Clone, PartialEq)]
pub struct TbpDensityResult {
    /// The cut's normal liquid density at 15 C.
    pub density: MassDensity,
    /// Caveats.
    pub warnings: Vec<Warning>,
}

impl CalcResult for TbpDensityResult {
    const CALC_ID: &'static str = "characterization.tbp_density";
    const FIELDS: &'static [&'static str] = &["density", "warnings"];

    fn warnings(&self) -> &[Warning] {
        &self.warnings
    }
}

/// Result of `characterization.set_heavy_tbp_to_plus`.
#[derive(Debug, Clone, PartialEq)]
pub struct SetHeavyTbpToPlusResult {
    /// Which row the plus end begins at, as a position in the fluid's TBP slate.
    pub plus_index: f64,
    /// That row's name, which is the string its carbon number was read out of.
    pub plus_name: String,
    /// That row's molar mass, which the class publishes as `MPlus`.
    pub plus_molar_mass: MolarMass,
    /// That row's mole fraction, which the class publishes as `zPlus`.
    pub plus_mole_fraction: f64,
    /// That row's normal liquid density, which the class publishes as `densPlus`.
    pub plus_density: MassDensity,
    /// The density of the row **below** that one in index order, which the class publishes as
    /// `densLastTBP`.
    pub last_tbp_density: MassDensity,
    /// `last_tbp_density` plus the class's own step, which it stores at `coefs[2]`.
    pub density_intercept: MassDensity,
    /// Caveats.
    pub warnings: Vec<Warning>,
}

impl CalcResult for SetHeavyTbpToPlusResult {
    const CALC_ID: &'static str = "characterization.set_heavy_tbp_to_plus";
    const FIELDS: &'static [&'static str] = &[
        "plus_index",
        "plus_name",
        "plus_molar_mass",
        "plus_mole_fraction",
        "plus_density",
        "last_tbp_density",
        "density_intercept",
        "warnings",
    ];

    fn warnings(&self) -> &[Warning] {
        &self.warnings
    }
}

/// Result of `characterization.tbp_cut_boundaries`.
#[derive(Debug, Clone, PartialEq)]
pub struct TbpCutBoundariesResult {
    /// Each cut's name, its prefix with its position from one appended.
    pub cut_name: Vec<String>,
    /// Each cut's share of the curve, aligned with `cut_name`.
    pub volume_fraction: Vec<f64>,
    /// Each cut's specific gravity, aligned with `cut_name`.
    pub specific_gravity: Vec<f64>,
    /// Each cut's lower boiling-point boundary, aligned with `cut_name`.
    pub lower_boiling_point: Vec<ThermodynamicTemperature>,
    /// Each cut's upper boiling-point boundary, aligned with `cut_name`.
    pub upper_boiling_point: Vec<ThermodynamicTemperature>,
    /// Each cut's average normal boiling point, the midpoint of its two boundaries.
    pub average_boiling_point: Vec<ThermodynamicTemperature>,
    /// Caveats.
    pub warnings: Vec<Warning>,
}

impl CalcResult for TbpCutBoundariesResult {
    const CALC_ID: &'static str = "characterization.tbp_cut_boundaries";
    const FIELDS: &'static [&'static str] = &[
        "cut_name",
        "volume_fraction",
        "specific_gravity",
        "lower_boiling_point",
        "upper_boiling_point",
        "average_boiling_point",
        "warnings",
    ];

    fn warnings(&self) -> &[Warning] {
        &self.warnings
    }
}

/// Result of `characterization.oil_assay_derivation`.
#[derive(Debug, Clone, PartialEq)]
pub struct OilAssayDerivationResult {
    /// Each derived component's name, in the slate's order.
    pub component_name: Vec<String>,
    /// Each component's resolved mass fraction, aligned with `component_name`.
    pub mass_fraction: Vec<f64>,
    /// Each component's molar mass, aligned with `component_name`.
    pub molar_mass: Vec<MolarMass>,
    /// Each component's mole amount, aligned with `component_name`.
    pub moles: Vec<f64>,
    /// The fluid's total mole amount.
    pub total_moles: f64,
    /// The class's own relative mass balance, which its `1e-10` gate is taken on.
    pub mass_balance_error: f64,
    /// Caveats.
    pub warnings: Vec<Warning>,
}

impl CalcResult for OilAssayDerivationResult {
    const CALC_ID: &'static str = "characterization.oil_assay_derivation";
    const FIELDS: &'static [&'static str] = &[
        "component_name",
        "mass_fraction",
        "molar_mass",
        "moles",
        "total_moles",
        "mass_balance_error",
        "warnings",
    ];

    fn warnings(&self) -> &[Warning] {
        &self.warnings
    }
}

/// One sub-fraction of a delumped lump, as `characterization.delumping` returns it.
///
/// **A `record_list` output's row, so its fields are the spec's own `fields` block.** The
/// spelling here is the transport's - a `uom` alias for a dimensioned field and a bare number
/// for a dimensionless one - and `tools/gen_python_record_outputs.py` emits the pyclass and the
/// Python dataclass from that same declaration rather than from this struct.
#[derive(Debug, Clone, PartialEq)]
pub struct SubFraction {
    /// The sub-fraction's name: its parent's, and the grid position appended to it.
    pub name: String,
    /// Its mole amount.
    pub moles: f64,
    /// Its molar mass.
    pub molar_mass: MolarMass,
    /// Its normal liquid density, which is its parent's unchanged.
    pub density: Option<MassDensity>,
    /// Its normal boiling point, spread with the molar mass - and its parent's where the parent
    /// states none.
    pub boiling_point: Option<f64>,
}

/// Result of `characterization.delumping`.
#[derive(Debug, Clone, PartialEq)]
pub struct DelumpingResult {
    /// The sub-fractions, in ascending sorting key.
    pub fractions: Vec<SubFraction>,
    /// Caveats.
    pub warnings: Vec<Warning>,
}

impl CalcResult for DelumpingResult {
    const CALC_ID: &'static str = "characterization.delumping";
    const FIELDS: &'static [&'static str] = &["fractions", "warnings"];

    fn warnings(&self) -> &[Warning] {
        &self.warnings
    }
}

/// Result of `characterization.combine_reservoir_fluids`.
#[derive(Debug, Clone, PartialEq)]
pub struct CombineReservoirFluidsResult {
    /// The blended fluid's components: its merged base components first, then the combined
    /// pseudo-components in ascending cut order, renamed `PC1_PC`, `PC2_PC` and so on.
    pub components: Vec<CharacterizationComponent>,
    /// The blended fluid's total mole amount, base and pseudo.
    pub total_moles: f64,
    /// Caveats.
    pub warnings: Vec<Warning>,
}

impl CalcResult for CombineReservoirFluidsResult {
    const CALC_ID: &'static str = "characterization.combine_reservoir_fluids";
    const FIELDS: &'static [&'static str] = &["components", "total_moles", "warnings"];

    fn warnings(&self) -> &[Warning] {
        &self.warnings
    }
}
