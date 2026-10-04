//! The five calculations, exposed to Python.
//!
//! # Why these take plain floats
//!
//! Every parameter below is an SI magnitude, not a `uom` quantity. Unit handling
//! happens once, in Python, in `azoth.hydraulics`, before the call crosses this
//! boundary.
//!
//! That is deliberate rather than a shortcut. If the extension did its own unit
//! conversion, the two backends would each have a unit path, and the failure they
//! could then produce - agreeing on the numbers while disagreeing about what the
//! numbers are in - is precisely the kind of bug this project is built to
//! prevent. One conversion site, used by both backends, cannot disagree with
//! itself.
//!
//! The cost is that `azoth._core` is not units-safe on its own. It is a private
//! module: the public API is `azoth.hydraulics`, which is.

use azoth_core::units::{
    kilograms_per_cubic_meter, kilograms_per_second, meters, millimeters, newtons_per_meter,
    pascal_seconds,
};
use azoth_hydraulics as hyd;
use pyo3::prelude::*;

use crate::errors::to_pyerr;

/// Pressure drop over a length of straight pipe.
///
/// `mu` is optional: it is needed only to check the flow regime. Omit it and the
/// result says the regime went unchecked rather than implying it passed.
/// A packed bed's hydraulics, from the packing's name and the state it runs at.
#[pyfunction]
#[pyo3(signature = (packing, column_diameter, packed_height, vapor_mass_flow, liquid_mass_flow, vapor_density, liquid_density, vapor_viscosity, liquid_viscosity, surface_tension, vapor_diffusivity, liquid_diffusivity, hydraulic_capacity_factor))]
#[pyo3(
    text_signature = "(packing, column_diameter, packed_height, vapor_mass_flow, liquid_mass_flow, vapor_density, liquid_density, vapor_viscosity, liquid_viscosity, surface_tension, vapor_diffusivity, liquid_diffusivity, hydraulic_capacity_factor)"
)]
#[allow(clippy::too_many_arguments)] // One per declared input, and a bed reads many.
pub fn packing_hydraulics(
    py: Python<'_>,
    packing: &str,
    column_diameter: f64,
    packed_height: f64,
    vapor_mass_flow: f64,
    liquid_mass_flow: f64,
    vapor_density: f64,
    liquid_density: f64,
    vapor_viscosity: f64,
    liquid_viscosity: f64,
    surface_tension: f64,
    vapor_diffusivity: f64,
    liquid_diffusivity: f64,
    hydraulic_capacity_factor: f64,
) -> PyResult<crate::transport_gen::PyPackingHydraulicsResult> {
    hyd::packing_hydraulics::packing_hydraulics(
        packing,
        hyd::packing_hydraulics::PackingState {
            column_diameter: meters(column_diameter),
            packed_height,
            vapor_mass_flow: kilograms_per_second(vapor_mass_flow),
            liquid_mass_flow: kilograms_per_second(liquid_mass_flow),
            vapor_density: kilograms_per_cubic_meter(vapor_density),
            liquid_density: kilograms_per_cubic_meter(liquid_density),
            vapor_viscosity: pascal_seconds(vapor_viscosity),
            liquid_viscosity: pascal_seconds(liquid_viscosity),
            surface_tension: newtons_per_meter(surface_tension),
            vapor_diffusivity,
            liquid_diffusivity,
            hydraulic_capacity_factor,
        },
    )
    .map(|r| crate::transport_gen::PyPackingHydraulicsResult::from(&r))
    .map_err(|e| to_pyerr(py, e))
}

/// A tray's flooding, weeping, entrainment, pressure drop and efficiency.
#[pyfunction]
#[pyo3(signature = (tray_type, column_diameter, tray_spacing, weir_height, weir_length, downcommer_area_fraction, hole_diameter, hole_area_fraction, design_flood_fraction, vapor_mass_flow, liquid_mass_flow, vapor_density, liquid_density, liquid_viscosity, surface_tension, relative_volatility))]
#[pyo3(
    text_signature = "(tray_type, column_diameter, tray_spacing, weir_height, weir_length, downcommer_area_fraction, hole_diameter, hole_area_fraction, design_flood_fraction, vapor_mass_flow, liquid_mass_flow, vapor_density, liquid_density, liquid_viscosity, surface_tension, relative_volatility)"
)]
#[allow(clippy::too_many_arguments)] // One per declared input, and a tray reads many.
pub fn tray_hydraulics(
    py: Python<'_>,
    tray_type: &str,
    column_diameter: f64,
    tray_spacing: f64,
    weir_height: f64,
    weir_length: f64,
    downcommer_area_fraction: f64,
    hole_diameter: f64,
    hole_area_fraction: f64,
    design_flood_fraction: f64,
    vapor_mass_flow: f64,
    liquid_mass_flow: f64,
    vapor_density: f64,
    liquid_density: f64,
    liquid_viscosity: f64,
    surface_tension: f64,
    relative_volatility: f64,
) -> PyResult<crate::transport_gen::PyTrayHydraulicsResult> {
    hyd::tray_hydraulics::tray_hydraulics(hyd::tray_hydraulics::TrayHydraulicsState {
        tray_type: tray_type.to_string(),
        column_diameter: meters(column_diameter),
        tray_spacing: meters(tray_spacing),
        weir_height: meters(weir_height),
        weir_length: meters(weir_length),
        downcommer_area_fraction,
        // **Millimetres in the spec, metres on the wire**, which is the same conversion every
        // other length crosses with.
        hole_diameter: millimeters(hole_diameter),
        hole_area_fraction,
        design_flood_fraction,
        vapor_mass_flow: kilograms_per_second(vapor_mass_flow),
        liquid_mass_flow: kilograms_per_second(liquid_mass_flow),
        vapor_density: kilograms_per_cubic_meter(vapor_density),
        liquid_density: kilograms_per_cubic_meter(liquid_density),
        liquid_viscosity: pascal_seconds(liquid_viscosity),
        surface_tension: newtons_per_meter(surface_tension),
        relative_volatility,
    })
    .map(|r| crate::transport_gen::PyTrayHydraulicsResult::from(&r))
    .map_err(|e| to_pyerr(py, e))
}

/// The internal diameter a packed bed needs at a chosen fraction of flood.
#[pyfunction]
#[pyo3(signature = (packing, design_flood_fraction, vapor_mass_flow, liquid_mass_flow, vapor_density, liquid_density, liquid_viscosity, hydraulic_capacity_factor))]
#[pyo3(
    text_signature = "(packing, design_flood_fraction, vapor_mass_flow, liquid_mass_flow, vapor_density, liquid_density, liquid_viscosity, hydraulic_capacity_factor)"
)]
#[allow(clippy::too_many_arguments)] // One per declared input.
pub fn packing_sizing(
    py: Python<'_>,
    packing: &str,
    design_flood_fraction: f64,
    vapor_mass_flow: f64,
    liquid_mass_flow: f64,
    vapor_density: f64,
    liquid_density: f64,
    liquid_viscosity: f64,
    hydraulic_capacity_factor: f64,
) -> PyResult<crate::transport_gen::PyPackingSizingResult> {
    hyd::packing_sizing::packing_sizing(
        packing,
        hyd::packing_sizing::PackingSizingState {
            design_flood_fraction,
            vapor_mass_flow: kilograms_per_second(vapor_mass_flow),
            liquid_mass_flow: kilograms_per_second(liquid_mass_flow),
            vapor_density: kilograms_per_cubic_meter(vapor_density),
            liquid_density: kilograms_per_cubic_meter(liquid_density),
            liquid_viscosity: pascal_seconds(liquid_viscosity),
            hydraulic_capacity_factor,
        },
    )
    .map(|r| crate::transport_gen::PyPackingSizingResult::from(&r))
    .map_err(|e| to_pyerr(py, e))
}
