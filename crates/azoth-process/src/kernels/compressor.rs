//! `unit_ops.compressor` - the isentropic route, through the kernel the pump landed.

use azoth_core::Result;
use azoth_core::units::Pressure;

use crate::kernels::pump::pump;
use crate::stream::Stream;

/// Raise a stream's pressure along an isentrope, dividing the step by the efficiency.
///
/// **This is `unit_ops.pump`'s arithmetic, and NeqSim's own code says so.** `Pump.run` takes
/// the compressor route because `calculateAsCompressor` defaults to `true`, and
/// `Compressor.run`'s no-chart branch is that route verbatim: set the outlet pressure,
/// `PSflash` at the *inlet's own entropy*, and move the enthalpy by
/// `(h_isentropic - h_in) / eta`. Writing the arithmetic out again here would be a second
/// implementation of one step, so this calls it.
///
/// **The entropy is derived, not carried.** The record has five fields and `s` is not one of
/// them: it is a function of `(T, P, z)`, and the flash the kernel already runs at the inlet
/// decides it. The captured isentropic row is the demonstration - at an efficiency of one the
/// outlet's entropy equals the inlet's to `1.3e-14`, and the step's entropy production is
/// `-1.4e-17` kJ/(mol K), which is zero.
///
/// **The shaft power is not an output**, for the reason `s` is not a field: `Compressor.getPower`
/// is `dH = h_out - h_in` over the machine, so the record this returns already carries it -
/// `n * (outlet_h - inlet_h)`, which the captured rows satisfy to their last digit.
///
/// Named out: the compressor chart and catalogue, the polytropic branch (the same rule with
/// the polytropic efficiency on the other side of the division), the gamma-based
/// pressure-ratio mode, the anti-surge machinery, the deposit model, and the
/// `useGERG2008`/`useLeachman`/`useVega` property overrides.
pub fn compressor(feed: &Stream, outlet_pressure: Pressure, efficiency: f64) -> Result<Stream> {
    pump(feed, outlet_pressure, efficiency)
}
