//! `unit_ops.heat_exchanger` - two streams, coupled by an overall conductance.

use azoth_core::units::{
    MolarEnergy, Power, ThermalConductance, ThermodynamicTemperature, joules_per_mole, watts,
    watts_per_kelvin,
};
use azoth_core::{AzothError, Result};

use crate::stream::Stream;

/// The flow arrangements `HeatExchanger` names, in azoth's own spelling.
///
/// NeqSim's are `"concentric tube counterflow"`, `"concentric tube paralellflow"` - its own
/// misspelling - and `"shell and tube"`, and its default is the first. An arrangement it
/// does not recognise takes the counterflow relation through `run`'s `else`; **this refuses
/// one instead**, because a typo that silently becomes counterflow is a plausible answer to
/// a question nobody asked.
pub const FLOW_ARRANGEMENTS: [&str; 3] = ["counterflow", "parallelflow", "shell_and_tube"];

/// Which side of the exchanger a stream is. Only the `outTemperature` mode needs it.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Side {
    Hot,
    Cold,
}

impl Side {
    fn other(self) -> Side {
        match self {
            Side::Hot => Side::Cold,
            Side::Cold => Side::Hot,
        }
    }
}

/// The numbers the effectiveness-NTU rating sized the exchanger by.
///
/// **Present only when a rating was asked for**, which is the default mode and not the only one:
/// a run handed a pinned outlet temperature energy-balances the other side instead, and reaches
/// the duty without any of these. That is why [`HeatExchangerOutcome::rating`] is an `Option` and
/// why the record writes `null` for them - the absence of a rating, not a rating of zero.
pub struct Rating {
    /// `UA / C_min`, dimensionless.
    pub ntu: f64,
    /// The swing the relation scaled, as a fraction of the capacity-limited side's.
    pub effectiveness: f64,
    /// The smaller of the two estimated capacities, W/K.
    pub c_min: ThermalConductance,
    /// The larger of them, W/K.
    pub c_max: ThermalConductance,
    /// `C_min / C_max`, dimensionless.
    pub capacity_ratio: f64,
}

/// What an exchanger did: both outlet records, the duty, and the rating where there was one.
pub struct HeatExchangerOutcome {
    /// The hot outlet record.
    pub hot_out: Stream,
    /// The cold outlet record.
    pub cold_out: Stream,
    /// The heat the hot side released, W - positive when the hot side cools, which is what the
    /// cold side's gain balances. **The same number either mode reaches**, and the one figure a
    /// faceplate shows whatever the exchanger was told to do.
    pub duty: Power,
    /// The rating's own numbers, or `None` where one outlet was pinned instead.
    pub rating: Option<Rating>,
}

impl HeatExchangerOutcome {
    /// The outcome of one solve, with the duty read off the hot side's own two states.
    ///
    /// **The duty is derived rather than carried from the branch that computed it**, because the
    /// two branches hold it in different terms - the rating scales a swing it chose, the pinned
    /// solve releases the swing it pinned - and the hot side's own enthalpy change is the one
    /// expression both reach. It is also what makes the sign one rule instead of two.
    #[must_use]
    fn of(hot: &Stream, hot_out: Stream, cold_out: Stream, rating: Option<Rating>) -> Self {
        let duty = watts(-total_swing(hot, &hot_out));
        Self {
            hot_out,
            cold_out,
            duty,
            rating,
        }
    }
}

/// Exchange heat between two streams.
///
/// Two modes, and which one runs is what the caller states:
///
/// * **`hot_outlet_temperature` or `cold_outlet_temperature`** pins that side's outlet at
///   its inlet pressure, and the other side is energy-balanced against it. Exactly one may
///   be given.
/// * **the effectiveness-NTU rating**, otherwise. `ua` is the overall conductance and
///   `flow_arrangement` picks the relation.
///
/// The rating follows `HeatExchanger.run`: each side's capacity is estimated by flashing it
/// at the **other side's inlet temperature** and dividing the total enthalpy change by the
/// temperature change, `NTU = UA / C_min`, and the effectiveness scales that swing. When
/// the effectiveness comes out at one the side that was not energy-balanced keeps the state
/// it was seeded at, which is the whole of how `run` expresses the infinite-area limit.
///
/// # Errors
/// [`azoth_core::AzothError::InvalidInput`] if neither or both outlet temperatures are
/// given, if the rating is asked for without a `ua`, if the arrangement is not one of
/// [`FLOW_ARRANGEMENTS`], or if the two inlets are at the same temperature so that there is
/// no driving force. Whatever the flashes refuse propagates too.
pub fn heat_exchanger(
    hot: &Stream,
    cold: &Stream,
    ua: Option<ThermalConductance>,
    flow_arrangement: &str,
    hot_outlet_temperature: Option<ThermodynamicTemperature>,
    cold_outlet_temperature: Option<ThermodynamicTemperature>,
) -> Result<HeatExchangerOutcome> {
    if hot_outlet_temperature.is_some() && cold_outlet_temperature.is_some() {
        return Err(AzothError::invalid_input(
            "hot_outlet_temperature",
            "an exchanger pins one outlet at most: with both stated the other side has \
             nothing left to solve for",
        ));
    }
    if let Some(pinned) = hot_outlet_temperature {
        let (hot_out, cold_out) = specified_outlet(hot, cold, Side::Hot, pinned)?;
        return Ok(HeatExchangerOutcome::of(hot, hot_out, cold_out, None));
    }
    if let Some(pinned) = cold_outlet_temperature {
        let (hot_out, cold_out) = specified_outlet(hot, cold, Side::Cold, pinned)?;
        return Ok(HeatExchangerOutcome::of(hot, hot_out, cold_out, None));
    }

    let Some(ua) = ua else {
        return Err(AzothError::invalid_input(
            "ua",
            "the rating needs an overall conductance. NeqSim's `UAvalue` defaults to 500 \
             W/K, and this refuses to guess it: a default that only changes the answer is \
             not one a caller can see",
        ));
    };
    if !FLOW_ARRANGEMENTS.contains(&flow_arrangement) {
        return Err(AzothError::invalid_input(
            "flow_arrangement",
            format!(
                "{flow_arrangement} is not one of {FLOW_ARRANGEMENTS:?}. NeqSim's `run` \
                 falls through to the counterflow relation for an arrangement it does not \
                 know, which would make a misspelling a plausible answer"
            ),
        ));
    }
    if (hot.t.value - cold.t.value).abs()
        <= f64::EPSILON * hot.t.value.abs().max(cold.t.value.abs())
    {
        return Err(AzothError::invalid_input(
            "hot_in",
            format!(
                "both inlets are at {} K, so there is no driving force for the rating to \
                 size against",
                hot.t.value
            ),
        ));
    }

    // The seeded capacities: each side flashed at the *other's* inlet temperature, which
    // is `run`'s own estimate and is what makes the pass count not matter.
    let hot_seeded = Stream::from_pt(hot.components.clone(), hot.z.clone(), hot.n, hot.p, cold.t)?;
    let cold_seeded = Stream::from_pt(
        cold.components.clone(),
        cold.z.clone(),
        cold.n,
        cold.p,
        hot.t,
    )?;

    // **Absolute values, because a swing is signed and a capacity is not.** The hot side's
    // seeded swing is negative - it cools - and taking it at face value would make `C_min`
    // negative and the NTU with it, which is an exponential in the wrong direction rather
    // than an error.
    let span = (hot.t.value - cold.t.value).abs();
    let cold_capacity = total_swing(cold, &cold_seeded).abs() / span;
    let hot_capacity = total_swing(hot, &hot_seeded).abs() / span;
    let (c_min, c_max) = if cold_capacity < hot_capacity {
        (cold_capacity, hot_capacity)
    } else {
        (hot_capacity, cold_capacity)
    };
    let capacity_ratio = c_min / c_max;
    let ntu = ua.value / c_min;
    let effectiveness = effectiveness(ntu, capacity_ratio, flow_arrangement);
    // **The rating's own numbers, kept rather than dropped.** They are what sized the exchanger,
    // and a faceplate that shows a duty and nothing else cannot say whether the answer was
    // capacity-limited or area-limited - which is the difference between a bigger machine helping
    // and it not helping at all.
    let rating = Rating {
        ntu,
        effectiveness,
        c_min: watts_per_kelvin(c_min),
        c_max: watts_per_kelvin(c_max),
        capacity_ratio,
    };

    // `run` energy-balances the side it swapped *away* from, and the side it keeps is the
    // one whose seeded swing is the larger. Both are totals, so the comparison is between
    // the two sides' actual duties and not between their molar capacities.
    let hot_swing = total_swing(hot, &hot_seeded);
    let cold_swing = total_swing(cold, &cold_seeded);
    let set_side = if cold_swing.abs() > hot_swing.abs() {
        Side::Hot
    } else {
        Side::Cold
    };
    let (set, set_seeded) = match set_side {
        Side::Hot => (hot, &hot_seeded),
        Side::Cold => (cold, &cold_seeded),
    };
    let calculated = match set_side.other() {
        Side::Hot => hot,
        Side::Cold => cold,
    };
    let d = effectiveness * total_swing(set, set_seeded);

    let calculated_out = Stream::from_ph(
        calculated.components.clone(),
        calculated.z.clone(),
        calculated.n,
        calculated.p,
        moles_with(calculated, -d),
    )?;
    // **At an effectiveness of one the set side keeps its seeded state.** That is not an
    // omission in `run`: the seeded state is "this side at the other's inlet temperature",
    // which is exactly what an effectiveness of one means, and the branch exists because
    // re-flashing it would find the same state by a longer route.
    let set_out = if (effectiveness - 1.0).abs() > 1e-10 {
        Stream::from_ph(
            set.components.clone(),
            set.z.clone(),
            set.n,
            set.p,
            moles_with(set, d),
        )?
    } else {
        set_seeded.clone()
    };

    let (hot_out, cold_out) = match set_side {
        Side::Hot => (set_out, calculated_out),
        Side::Cold => (calculated_out, set_out),
    };
    Ok(HeatExchangerOutcome::of(
        hot,
        hot_out,
        cold_out,
        Some(rating),
    ))
}

/// The `outTemperature` specification: pin one outlet, energy-balance the other.
fn specified_outlet(
    hot: &Stream,
    cold: &Stream,
    side: Side,
    pinned: ThermodynamicTemperature,
) -> Result<(Stream, Stream)> {
    let (specified, other) = match side {
        Side::Hot => (hot, cold),
        Side::Cold => (cold, hot),
    };
    let specified_out = Stream::from_pt(
        specified.components.clone(),
        specified.z.clone(),
        specified.n,
        specified.p,
        pinned,
    )?;
    // The duty the pinned side released, as a total: `run` compares and applies totals.
    let d = total_swing(specified, &specified_out);
    let other_out = Stream::from_ph(
        other.components.clone(),
        other.z.clone(),
        other.n,
        other.p,
        moles_with(other, -d),
    )?;
    Ok(match side {
        Side::Hot => (specified_out, other_out),
        Side::Cold => (other_out, specified_out),
    })
}

/// A side's total enthalpy change against a second state of the same stream, in W.
fn total_swing(from: &Stream, to: &Stream) -> f64 {
    to.n * (to.h.value - from.h.value)
}

/// A stream's molar enthalpy shifted by `total` watts over its own flow.
fn moles_with(stream: &Stream, total: f64) -> MolarEnergy {
    joules_per_mole(stream.h.value + total / stream.n)
}

/// The effectiveness-NTU relations `HeatExchanger.calcThermalEffectivenes` carries.
///
/// NeqSim spells the counterflow and shell-and-tube cases with the same expression and the
/// parallel-flow one differently; that is reproduced rather than tidied, because it is what
/// the class answers and `shell and tube` is a name a caller can state.
fn effectiveness(ntu: f64, capacity_ratio: f64, arrangement: &str) -> f64 {
    if capacity_ratio == 0.0 {
        return 1.0 - (-ntu).exp();
    }
    if arrangement == "parallelflow" {
        return (1.0 - (-ntu * (1.0 + capacity_ratio)).exp()) / (1.0 + capacity_ratio);
    }
    if arrangement == "counterflow" && capacity_ratio == 1.0 {
        return ntu / (1.0 + ntu);
    }
    let exp = (-ntu * (1.0 - capacity_ratio)).exp();
    (1.0 - exp) / (1.0 - capacity_ratio * exp)
}
