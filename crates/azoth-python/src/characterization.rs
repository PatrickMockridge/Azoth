//! The characterisation calculations, exposed to Python.
//!
//! Same rule as the other bindings: every parameter below is an SI magnitude, not a `uom`
//! quantity. Unit handling happens once, in Python, before the call crosses this boundary - one
//! conversion site used by both backends cannot disagree with itself about what a number is in.

use azoth_characterization::{
    AssayBasis, AssayCut as KernelAssayCut, PlusModel, TbpClosureKind, TbpModel,
    WhitsonDensityModel, assay_mass_fractions as assay_kernel,
    characterise_plus_fraction as facade_kernel, tbp_closure as closure_kernel,
    tbp_cut_properties as kernel, tbp_density as density_kernel,
    whitson_gamma_split as gamma_kernel,
};
use azoth_core::units::{kelvins, kilograms_per_cubic_meter, kilograms_per_mole};
use pyo3::prelude::*;

use crate::errors::to_pyerr;
use crate::transport_gen::{
    PyAssayMassFractionsResult, PyCharacterisePlusFractionResult, PyTbpClosureResult,
    PyTbpCutPropertiesResult, PyTbpDensityResult, PyWhitsonGammaSplitResult,
};

/// A C7+ end characterised end to end: model, split and lumps.
///
/// All arguments are SI magnitudes. `plus_model` names one of the three the spec declares;
/// `None` takes `pedersen`. **`selected_model` on the result is not always this**: a plus
/// fraction heavier than the requested model's maximum is re-modelled, silently. **This wrapper
/// is hand-written** because an optional enum is a shape `gen_python_wrappers` refuses.
#[pyfunction]
#[pyo3(signature = (molar_mass, density, mole_fraction, first_carbon_number, plus_model=None, number_of_lumps=None))]
#[pyo3(text_signature = "(molar_mass, density, mole_fraction, first_carbon_number, plus_model=None, number_of_lumps=None)")]
pub fn characterise_plus_fraction(
    py: Python<'_>,
    molar_mass: f64,
    density: f64,
    mole_fraction: f64,
    first_carbon_number: usize,
    plus_model: Option<&str>,
    number_of_lumps: Option<usize>,
) -> PyResult<PyCharacterisePlusFractionResult> {
    let parsed: Option<PlusModel> = plus_model.map(|name| name.parse().unwrap_or_default());
    facade_kernel(
        kilograms_per_mole(molar_mass),
        kilograms_per_cubic_meter(density),
        mole_fraction,
        first_carbon_number,
        parsed,
        number_of_lumps,
    )
    .map(|result| PyCharacterisePlusFractionResult::from(&result))
    .map_err(|error| to_pyerr(py, error))
}

/// One cut of an assay, as it crosses the boundary.
///
/// **A record rather than a vector entry**, and the density is why: `Option` per *cut* has no
/// aligned-vector form. Every field is an SI magnitude already - this file's rule - so the struct
/// is a carrier and the kernel builds the quantities.
///
/// It is named by the spec's own `row_type`, which is what `gen_stub` generates the stub's
/// `TypedDict` under and what the agreement test looks the Rust type up by; the name has to be
/// the same on both sides or the three declarations of one row are three shapes.
#[derive(Debug, Clone, Copy, Default, FromPyObject)]
pub struct AssayCut {
    #[pyo3(item)]
    pub declared_fraction: f64,
    #[pyo3(item)]
    pub density: Option<f64>,
}

/// An oil assay's declared yields, resolved to a mass basis.
///
/// All arguments are SI magnitudes. `basis` is one of the two the spec declares, and `cuts` is
/// the assay's rows - each a mapping keyed by the spec's `field` names, with an absent key
/// meaning the cut states none. **This wrapper is hand-written** because a slice of records is a
/// shape `gen_python_wrappers` refuses rather than guesses at.
#[pyfunction]
#[pyo3(signature = (basis, cuts))]
#[pyo3(text_signature = "(basis, cuts)")]
pub fn assay_mass_fractions(
    py: Python<'_>,
    basis: &str,
    cuts: Vec<AssayCut>,
) -> PyResult<PyAssayMassFractionsResult> {
    let parsed: AssayBasis = basis.parse().unwrap_or_default();
    let rows: Vec<KernelAssayCut> = cuts
        .into_iter()
        .map(|cut| KernelAssayCut {
            declared_fraction: cut.declared_fraction,
            density: cut.density.map(kilograms_per_cubic_meter),
        })
        .collect();
    assay_kernel(parsed, &rows)
        .map(|result| PyAssayMassFractionsResult::from(&result))
        .map_err(|error| to_pyerr(py, error))
}

/// A plus fraction split into cuts by Whitson's three-parameter gamma distribution.
///
/// All arguments are SI magnitudes. `eta` is absent by default, taking the class's own 90 g/mol;
/// `density_model` names one of the two gravity correlations, `uop` by default.
#[pyfunction]
#[pyo3(signature = (molar_mass, density, mole_fraction, first_carbon_number, last_carbon_number, alpha=None, eta=None, density_model=None, auto_estimate_shape=None))]
#[pyo3(text_signature = "(molar_mass, density, mole_fraction, first_carbon_number, last_carbon_number, alpha=None, eta=None, density_model=None, auto_estimate_shape=None)")]
#[allow(clippy::too_many_arguments)] // one per declared input, as the spec states them
pub fn whitson_gamma_split(
    py: Python<'_>,
    molar_mass: f64,
    density: f64,
    mole_fraction: f64,
    first_carbon_number: usize,
    last_carbon_number: usize,
    alpha: Option<f64>,
    eta: Option<f64>,
    density_model: Option<&str>,
    auto_estimate_shape: Option<bool>,
) -> PyResult<PyWhitsonGammaSplitResult> {
    let parsed: WhitsonDensityModel = density_model
        .map(|name| name.parse().unwrap_or_default())
        .unwrap_or_default();
    gamma_kernel(
        kilograms_per_mole(molar_mass),
        kilograms_per_cubic_meter(density),
        mole_fraction,
        first_carbon_number,
        last_carbon_number,
        alpha,
        eta.map(kilograms_per_mole),
        parsed,
        auto_estimate_shape.unwrap_or(false),
    )
    .map(|result| PyWhitsonGammaSplitResult::from(&result))
    .map_err(|error| to_pyerr(py, error))
}

/// A cut's molar mass, from its normal boiling point and its specific gravity.
///
/// All arguments are SI magnitudes. `closure` is one of the four members the spec's enum
/// declares; `model` is needed only for `tbp_model`.
#[pyfunction]
#[pyo3(signature = (boiling_point, density, closure, model=None))]
#[pyo3(text_signature = "(boiling_point, density, closure, model=None)")]
pub fn tbp_closure(
    py: Python<'_>,
    boiling_point: f64,
    density: f64,
    closure: &str,
    model: Option<&str>,
) -> PyResult<PyTbpClosureResult> {
    let parsed_model: Option<TbpModel> = model.map(|name| name.parse().unwrap_or_default());
    closure_kernel(
        closure.parse().unwrap_or_default(),
        kelvins(boiling_point),
        kilograms_per_cubic_meter(density),
        parsed_model,
    )
    .map(|result| PyTbpClosureResult::from(&result))
    .map_err(|error| to_pyerr(py, error))
}

/// A cut's normal liquid density, from its normal boiling point and its molar mass.
///
/// Only `riazi_daubert_1980` supports this direction; the other three closures are refused.
#[pyfunction]
#[pyo3(signature = (boiling_point, molar_mass, closure=None))]
#[pyo3(text_signature = "(boiling_point, molar_mass, closure=None)")]
pub fn tbp_density(
    py: Python<'_>,
    boiling_point: f64,
    molar_mass: f64,
    closure: Option<&str>,
) -> PyResult<PyTbpDensityResult> {
    let parsed: Option<TbpClosureKind> = closure.map(|name| name.parse().unwrap_or_default());
    density_kernel(
        parsed.unwrap_or_default(),
        kelvins(boiling_point),
        kilograms_per_mole(molar_mass),
    )
    .map(|result| PyTbpDensityResult::from(&result))
    .map_err(|error| to_pyerr(py, error))
}

/// A TBP cut's critical properties, by any of NeqSim's ten models.
///
/// All arguments are SI magnitudes. `model` is one of the names the spec's enum declares;
/// `None` takes NeqSim's own default. An unrecognised name falls back to that default, which is
/// what the NeqSim class does - the spec's vocabulary check is what refuses.
#[pyfunction]
#[pyo3(signature = (molar_mass, density, model=None, boiling_point=None))]
#[pyo3(text_signature = "(molar_mass, density, model=None, boiling_point=None)")]
pub fn tbp_cut_properties(
    py: Python<'_>,
    molar_mass: f64,
    density: f64,
    model: Option<&str>,
    boiling_point: Option<f64>,
) -> PyResult<PyTbpCutPropertiesResult> {
    let parsed: Option<TbpModel> = model.map(|name| name.parse().unwrap_or_default());
    kernel(
        parsed,
        kilograms_per_mole(molar_mass),
        kilograms_per_cubic_meter(density),
        boiling_point.map(kelvins),
    )
    .map(|result| PyTbpCutPropertiesResult::from(&result))
    .map_err(|error| to_pyerr(py, error))
}
