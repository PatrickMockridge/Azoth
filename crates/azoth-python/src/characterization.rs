//! The characterisation calculations, exposed to Python.
//!
//! Same rule as the other bindings: every parameter below is an SI magnitude, not a `uom`
//! quantity. Unit handling happens once, in Python, before the call crosses this boundary - one
//! conversion site used by both backends cannot disagree with itself about what a number is in.

use azoth_characterization::{
    AssayBasis, PlusModel, TbpClosureKind, TbpModel, WhitsonDensityModel,
    assay_mass_fractions as assay_kernel, characterise_plus_fraction as facade_kernel,
    tbp_closure as closure_kernel,
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

/// An oil assay's declared yields, resolved to a mass basis.
///
/// All arguments are SI magnitudes. `basis` is one of the two the spec declares; `density` is one
/// entry per cut and is needed by a volume basis, and only wanted by a mass basis for the bulk
/// density. **This wrapper is hand-written** because an optional dimensioned slice is a shape
/// `gen_python_wrappers` refuses rather than guesses at.
#[pyfunction]
#[pyo3(signature = (basis, declared_fraction, density=None))]
#[pyo3(text_signature = "(basis, declared_fraction, density=None)")]
pub fn assay_mass_fractions(
    py: Python<'_>,
    basis: &str,
    declared_fraction: Vec<f64>,
    density: Option<Vec<f64>>,
) -> PyResult<PyAssayMassFractionsResult> {
    let parsed: AssayBasis = basis.parse().unwrap_or_default();
    let densities: Option<Vec<_>> = density
        .map(|values| values.into_iter().map(kilograms_per_cubic_meter).collect());
    assay_kernel(parsed, &declared_fraction, densities.as_deref())
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
