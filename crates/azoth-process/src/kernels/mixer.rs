//! `unit_ops.mixer` - many streams joined into one, conserving molar flow and
//! enthalpy.

use azoth_core::AzothError;
use azoth_core::units::{Pressure, joules_per_mole};

use crate::stream::Stream;

/// Join several inlets into one outlet, conserving molar flow and enthalpy.
///
/// The outlet composition is the flow-weighted average of the inlets, and its
/// enthalpy is the flow-weighted average too - mixing is isenthalpic. The outlet
/// temperature is the one at which the joined mixture has that enthalpy, found
/// by [`azoth_eos::ph_flash`].
///
/// The outlet pressure is `outlet_pressure` when given, else the lowest inlet
/// pressure: a mixer cannot raise its own pressure.
pub fn mixer(inlets: &[Stream], outlet_pressure: Option<Pressure>) -> azoth_core::Result<Stream> {
    if inlets.is_empty() {
        return Err(AzothError::invalid_input(
            "inlets",
            "a mixer needs at least one feed",
        ));
    }

    let components = inlets[0].components.clone();
    let n_components = components.len();

    let n_total: f64 = inlets.iter().map(|s| s.n).sum();
    let mut z = vec![0.0; n_components];
    let mut h_total = 0.0;
    for inlet in inlets {
        for (i, &zi) in inlet.z.iter().enumerate() {
            z[i] += inlet.n * zi;
        }
        h_total += inlet.n * inlet.h.value;
    }
    for zi in &mut z {
        *zi /= n_total;
    }

    let pressure = outlet_pressure.unwrap_or_else(|| {
        inlets
            .iter()
            .map(|s| s.p)
            .min_by(|a, b| a.value.total_cmp(&b.value))
            .expect("a mixer with inlets has a lowest pressure")
    });

    Stream::from_ph(
        components,
        z,
        n_total,
        pressure,
        joules_per_mole(h_total / n_total),
    )
}
