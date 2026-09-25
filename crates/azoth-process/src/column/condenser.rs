//! The top end: `Condenser`, which is a stage with a reflux split.

use azoth_core::units::{Power, Pressure, ThermodynamicTemperature, watts};
use azoth_core::{AzothError, Result};
use azoth_eos::bubble_temperature::bubble_temperature;
use azoth_eos::pv_reflux_flash::{RefluxPhase, pv_reflux_flash};

use crate::kernels::mixer;
use crate::stream::Stream;

use super::phase_fractions;
use super::tray;
use super::tray::SideDraws;

/// How the condenser is specified.
///
/// The three modes are `DistillationColumn.CondenserMode`'s, and the first two differ only in
/// whether a reflux ratio was set - which is what `Condenser.run` tests.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum CondenserMode {
    /// No reflux set: `super.run`, the stage's own flash. An equilibrium **partial**
    /// condenser, and `CondenserMode.PARTIAL` with no ratio.
    Equilibrium,
    /// `PVrefluxflash(ratio, 0)`: a partial condenser at the temperature whose *vapour*
    /// fraction satisfies `ratio = 1/beta_V - 1`. `CondenserMode.PARTIAL` with a ratio.
    RefluxRatio(f64),
    /// A bubble-point **total** condenser whose condensate is split at `ratio/(1 + ratio)`.
    /// `CondenserMode.TOTAL`.
    Total(f64),
    /// A partial condenser with an explicit fixed liquid reflux **flow**, mol/s, taken from
    /// the liquid outlet. `CondenserMode.LIQUID_REFLUX_SPLIT`.
    LiquidRefluxSplit(f64),
}

/// What the top end hands the column.
///
/// **The distillate arrives through NeqSim's *gas* getter in the total mode.** For a total
/// condenser `getGasOutStream()` returns `mixedStreamSplitter.getSplitStream(1)`, which is
/// liquid, and `getLiquidOutStream()` returns `getSplitStream(0)`, which is the reflux. The
/// names here say what the streams *are*; a port that kept the class's names would have a
/// "vapour" that is a liquid.
#[derive(Debug, Clone)]
pub struct CondenserOutcome {
    /// The end's temperature, which is its flash's.
    pub temperature: ThermodynamicTemperature,
    /// The end's pressure.
    pub pressure: Pressure,
    /// The distillate leaving the column, or `None` where the flash found no vapour.
    pub distillate: Option<Stream>,
    /// The reflux returning to the column, or `None` where there is no liquid to return.
    pub reflux: Option<Stream>,
    /// The separate liquid product, where the mode makes one. Only
    /// [`CondenserMode::LiquidRefluxSplit`] does.
    pub liquid_product: Option<Stream>,
    /// The heat the end needs, W. Negative for a condenser: it is the enthalpy the products
    /// carry out less the enthalpy the vapour brought in.
    pub duty: Power,
}

/// The column's top end.
///
/// `Condenser.run`'s four branches, in its own order of testing:
///
/// 1. a reflux ratio **and** `setTotalCondenser(true)` - a bubble-point flash, and then a
///    split of the condensate at `R/(1 + R)` into reflux and distillate;
/// 2. no reflux set - `super.run`, the stage's own flash;
/// 3. a fixed liquid reflux - `super.run`, then a fixed **flow** taken from the liquid;
/// 4. otherwise - `PVrefluxflash(ratio, 0)`, a partial condenser at the ratio.
///
/// **In the total mode both products are liquid.** The flash is a bubble-point one, so there
/// is no vapour at all; the "distillate" is the split fraction the class exposes through its
/// gas getter. This models the four branches as four modes rather than as four flags, because
/// `Condenser.run` is itself a chain of tests on those flags and a mode names which branch was
/// taken.
///
/// # Errors
/// [`azoth_core::AzothError::InvalidInput`] for a negative or non-finite ratio, for a
/// fixed reflux flow that is negative, and whatever the stage's or the flash's refusals are.
pub fn condenser(
    inlets: &[Stream],
    tray_pressure: Option<Pressure>,
    mode: CondenserMode,
) -> Result<CondenserOutcome> {
    let inlets_enthalpy: f64 = inlets.iter().map(|s| s.n * s.h.value).sum();

    match mode {
        CondenserMode::Equilibrium => {
            let out = tray::tray(
                inlets,
                tray_pressure,
                None,
                watts(0.0),
                SideDraws::NONE,
                false,
            )?;
            let duty =
                enthalpy_of(out.gas.as_ref()) + enthalpy_of(out.liquid.as_ref()) - inlets_enthalpy;
            Ok(CondenserOutcome {
                temperature: out.temperature,
                pressure: out.pressure,
                distillate: out.gas,
                reflux: out.liquid,
                liquid_product: None,
                duty: watts(duty),
            })
        }
        CondenserMode::RefluxRatio(ratio) => {
            check_ratio(ratio)?;
            let mixed = mixer(inlets, tray_pressure)?;
            let (mixture, _ideal_gas) = mixed.mixture()?;
            let flash = pv_reflux_flash(
                &mixture,
                mixed.p,
                ratio,
                RefluxPhase::Vapour,
                mixed.t,
                &mixed.z,
            )?;
            let (vapour_fraction, liquid_fraction) = phase_fractions(flash.phase, flash.beta)?;
            let distillate = phase_stream(&mixed, flash.y, vapour_fraction, flash.t)?;
            let reflux = phase_stream(&mixed, flash.x, liquid_fraction, flash.t)?;
            let duty =
                enthalpy_of(distillate.as_ref()) + enthalpy_of(reflux.as_ref()) - inlets_enthalpy;
            Ok(CondenserOutcome {
                temperature: flash.t,
                pressure: mixed.p,
                distillate,
                reflux,
                liquid_product: None,
                duty: watts(duty),
            })
        }
        CondenserMode::Total(ratio) => {
            check_ratio(ratio)?;
            let mixed = mixer(inlets, tray_pressure)?;
            let (mixture, _ideal_gas) = mixed.mixture()?;
            // The class's own `bubblePointTemperatureFlash`, at the mixed composition.
            let bubble = bubble_temperature(&mixture, mixed.p, &mixed.z)?;
            let condensate = Stream::from_pt(
                mixed.components.clone(),
                mixed.z.clone(),
                mixed.n,
                mixed.p,
                bubble.temperature,
            )?;
            // `Splitter.setSplitFactors({R/(1+R), 1 - R/(1+R)})`, in moles: the first
            // fraction is the **reflux**, which `getLiquidOutStream()` returns.
            let reflux_fraction = if ratio <= 0.0 {
                0.0
            } else {
                ratio / (1.0 + ratio)
            };
            let reflux = scale(&condensate, reflux_fraction)?;
            let distillate = scale(&condensate, 1.0 - reflux_fraction)?;
            let duty = condensate.n * condensate.h.value - inlets_enthalpy;
            Ok(CondenserOutcome {
                temperature: bubble.temperature,
                pressure: mixed.p,
                distillate,
                reflux,
                liquid_product: None,
                duty: watts(duty),
            })
        }
        CondenserMode::LiquidRefluxSplit(reflux_flow) => {
            if !reflux_flow.is_finite() || reflux_flow < 0.0 {
                return Err(AzothError::invalid_input(
                    "reflux_flow",
                    format!(
                        "a fixed liquid reflux is a flow in mol/s, and {reflux_flow} is not one"
                    ),
                ));
            }
            let out = tray::tray(
                inlets,
                tray_pressure,
                None,
                watts(0.0),
                SideDraws::NONE,
                false,
            )?;
            let available = out.liquid.as_ref().map_or(0.0, |l| l.n);
            // `Splitter.setFlowRates({reflux, REMAINDER})`: the reflux is what was asked for
            // and the product is the remainder, which is why the class keeps a *shortfall*
            // residual - a reflux larger than the condensate is delivered short rather than
            // refused.
            let delivered = reflux_flow.min(available);
            let reflux = match out.liquid.as_ref() {
                Some(l) => scale(l, delivered / available)?,
                None => None,
            };
            let product_fraction = if available > 0.0 {
                (available - delivered) / available
            } else {
                0.0
            };
            let liquid_product = match out.liquid.as_ref() {
                Some(l) => scale(l, product_fraction)?,
                None => None,
            };
            // **All three products, and the vapour is one of them.** The class's
            // `getMaterialOutletEnthalpy` is the stage's gas and liquid plus the separate
            // liquid product - and in this mode `getLiquidOutStream()` is the *reflux branch*,
            // so the reflux and the product are the liquid's two halves and the vapour is the
            // third stream. A duty that omitted it would be the vapour's enthalpy, which is
            // what this read at `-10.74` W before the measurement caught it.
            let duty = enthalpy_of(out.gas.as_ref())
                + enthalpy_of(reflux.as_ref())
                + enthalpy_of(liquid_product.as_ref())
                - inlets_enthalpy;
            Ok(CondenserOutcome {
                temperature: out.temperature,
                pressure: out.pressure,
                distillate: out.gas,
                reflux,
                liquid_product,
                duty: watts(duty),
            })
        }
    }
}

/// A stream scaled to a fraction of itself, by moles.
///
/// **A zero fraction is `None`, not an empty stream.** NeqSim's `Splitter` publishes a
/// zero-flow branch - the capture's `total_reflux_ratio_zero` row carries a reflux of
/// `1e-50` mol/s - and a stream of no moles is not a state this library can express. The
/// absence of a branch is a fact, so it is spelled as one, the same way the tray spells an
/// absent phase.
fn scale(stream: &Stream, fraction: f64) -> Result<Option<Stream>> {
    if !fraction.is_finite() || fraction < 0.0 {
        return Err(AzothError::invalid_input(
            "split",
            format!("a split fraction is a fraction, and {fraction} is not one"),
        ));
    }
    if fraction == 0.0 {
        return Ok(None);
    }
    Ok(Some(Stream::from_pt(
        stream.components.clone(),
        stream.z.clone(),
        stream.n * fraction,
        stream.p,
        stream.t,
    )?))
}

/// One phase of a flash as a stream, or `None` where it is absent.
fn phase_stream(
    mixed: &Stream,
    z: Vec<f64>,
    fraction: f64,
    temperature: ThermodynamicTemperature,
) -> Result<Option<Stream>> {
    if fraction <= 0.0 {
        return Ok(None);
    }
    Ok(Some(Stream::from_pt(
        mixed.components.clone(),
        z,
        mixed.n * fraction,
        mixed.p,
        temperature,
    )?))
}

/// A ratio has to be a ratio.
fn check_ratio(ratio: f64) -> Result<()> {
    if !ratio.is_finite() || ratio < 0.0 {
        return Err(AzothError::invalid_input(
            "reflux_ratio",
            format!("a reflux ratio is a ratio, and {ratio} is not one"),
        ));
    }
    Ok(())
}

/// The total enthalpy of one stream, W: zero where the stream is absent, which is NeqSim's
/// own `getMaterialStreamEnthalpy` on a zero-flow stream.
fn enthalpy_of(stream: Option<&Stream>) -> f64 {
    stream.map_or(0.0, |s| s.n * s.h.value)
}
