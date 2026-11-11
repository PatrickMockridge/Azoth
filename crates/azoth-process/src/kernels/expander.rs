//! `unit_ops.expander` - the isentropic route with the efficiency on the other side.

use azoth_core::units::{
    Power, Pressure, ThermodynamicTemperature, joules_per_mole, joules_per_mole_kelvin, pascals,
};
use azoth_core::{AzothError, Result};
use azoth_eos::{ph_flash, ps_flash};

use crate::stream::Stream;

/// Drop a stream's pressure along an isentrope, multiplying the step by the efficiency.
///
/// **One operator away from `unit_ops.compressor`, and the operator is the point.**
/// `Expander.run`'s steady-state branch is `Compressor.run`'s with `* isentropicEfficiency`
/// where the compressor has `/ isentropicEfficiency` - and the two are the same physical rule,
/// because an expansion's isentropic enthalpy difference is *negative*: dividing it would
/// make the machine recover more work than the reversible one, with an efficiency below one.
/// At `eta = 1` both reduce to the isentropic step, which is the row each capture pins with an
/// entropy production of zero.
///
/// The route itself is `Pump.run`'s: the inlet's own entropy - derived from `(T, P, z)` and
/// not carried on the record - flashed to the outlet pressure for the reversible enthalpy, and
/// a second flash at the enthalpy the efficiency leaves.
///
/// Named out: the `powerSet` mode that takes a shaft power instead of an efficiency, the
/// `useOutTemperature` mode, and the same chart, anti-surge and property overrides the
/// compressor names. The polytropic branch is [`expander_polytropic`].
pub fn expander(feed: &Stream, outlet_pressure: Pressure, efficiency: f64) -> Result<Stream> {
    if !(efficiency > 0.0 && efficiency <= 1.0) {
        return Err(AzothError::invalid_input(
            "isentropic_efficiency",
            "the expander's isentropic efficiency must be in (0, 1]",
        ));
    }
    expander_replay(feed, outlet_pressure, efficiency)
}

/// The isentropic route without the upper bound, for the turbo-expander's replay.
///
/// **The coupled machine assigns its solved efficiency to the field directly**, bypassing
/// `setIsentropicEfficiency`'s `(0, 1]` clamp because the outlet-temperature inverse needs a
/// required efficiency that can exceed one on the way to a back-calculated design value. The
/// public [`expander`] keeps the clamp; this is the path the machine's replay takes.
///
/// # Errors
/// [`azoth_core::AzothError::InvalidInput`] for a non-positive efficiency, and whatever the
/// flash refuses.
pub fn expander_replay(
    feed: &Stream,
    outlet_pressure: Pressure,
    efficiency: f64,
) -> Result<Stream> {
    if efficiency.is_nan() || efficiency <= 0.0 {
        return Err(AzothError::invalid_input(
            "isentropic_efficiency",
            format!("the replay needs a positive efficiency, and {efficiency} is not"),
        ));
    }

    let (mixture, ideal_gas) = feed.mixture()?;
    let (entropy, _) = ps_flash::entropy_at(&mixture, &ideal_gas, feed.t, feed.p, &feed.z)?;

    let isentropic = ps_flash::ps_flash(
        &mixture,
        &ideal_gas,
        outlet_pressure,
        joules_per_mole_kelvin(entropy),
        &feed.z,
    )?;
    let (h_isentropic, _) = ph_flash::enthalpy_at(
        &mixture,
        &ideal_gas,
        isentropic.temperature,
        outlet_pressure,
        &feed.z,
    )?;

    let dh = (h_isentropic - feed.h.value) * efficiency;
    Stream::from_ph(
        feed.components.clone(),
        feed.z.clone(),
        feed.n,
        outlet_pressure,
        joules_per_mole(feed.h.value + dh),
    )
}

/// Drop a stream's pressure in steps, multiplying each isentropic step by the efficiency.
///
/// **`Expander.run`'s polytropic branch, `Expander.java:620-633`.** The pressure falls in
/// `numbersteps` equal increments; each step takes the *current* entropy and enthalpy, flashes
/// isentropically to the next pressure, and sets the real enthalpy to
/// `h + (h_is - h) eta` before a `PHflash` lands the state. The expander **multiplies** where
/// the compressor divides, because an expansion's isentropic difference is negative.
///
/// **The step count is floored at one**, which is the class's `Math.max(1, ...)`; a count of
/// zero is one step rather than an error.
///
/// # Errors
/// [`azoth_core::AzothError::InvalidInput`] if the efficiency is outside `(0, 1]`, and
/// whatever the flash refuses.
pub fn expander_polytropic(
    feed: &Stream,
    outlet_pressure: Pressure,
    polytropic_efficiency: f64,
    number_of_steps: usize,
) -> Result<Stream> {
    if !(polytropic_efficiency > 0.0 && polytropic_efficiency <= 1.0) {
        return Err(AzothError::invalid_input(
            "polytropic_efficiency",
            "the expander's polytropic efficiency must be in (0, 1]",
        ));
    }

    let steps = number_of_steps.max(1);
    let (mixture, ideal_gas) = feed.mixture()?;
    let increment = (outlet_pressure.value - feed.p.value) / steps as f64;
    let mut state = feed.clone();
    for _ in 0..steps {
        let (entropy, _) = ps_flash::entropy_at(&mixture, &ideal_gas, state.t, state.p, &state.z)?;
        let next_pressure = pascals(state.p.value + increment);
        let isentropic = ps_flash::ps_flash(
            &mixture,
            &ideal_gas,
            next_pressure,
            joules_per_mole_kelvin(entropy),
            &state.z,
        )?;
        let (h_isentropic, _) = ph_flash::enthalpy_at(
            &mixture,
            &ideal_gas,
            isentropic.temperature,
            next_pressure,
            &state.z,
        )?;
        let enthalpy = state.h.value + (h_isentropic - state.h.value) * polytropic_efficiency;
        state = Stream::from_ph(
            state.components.clone(),
            state.z.clone(),
            state.n,
            next_pressure,
            joules_per_mole(enthalpy),
        )?;
    }
    Ok(state)
}

/// `Compressor.MIN_EFFICIENCY`, the setter's lower clamp.
pub const MIN_EFFICIENCY: f64 = 1.0e-6;
/// `Compressor.MAX_EFFICIENCY`, the setter's upper clamp.
pub const MAX_EFFICIENCY: f64 = 1.0;
/// `solveEfficiency`'s outlet-temperature tolerance, K.
pub const EFFICIENCY_SOLVE_TOLERANCE: f64 = 0.01;

/// Drop a stream's pressure to a **stated shaft power** rather than an efficiency.
///
/// **`setPower` states the total enthalpy change, and the class back-computes the efficiency
/// from it**: `eta = dH / (h_isentropic - h_in)`, which can be above one or negative. The
/// reported efficiency is absent where the isentropic step is zero - the class answers `NaN`
/// there, and a value JSON cannot carry is published as the absence it is.
///
/// # Errors
/// [`azoth_core::AzothError::InvalidInput`] for a non-positive flow (the stated power is divided
/// by it), and whatever the flash refuses.
pub fn expander_at_power(
    feed: &Stream,
    outlet_pressure: Pressure,
    power: Power,
) -> Result<(Stream, Option<f64>)> {
    if feed.n <= 0.0 || !feed.n.is_finite() {
        return Err(AzothError::invalid_input(
            "shaft_power",
            format!(
                "the stated power divides by the molar flow, and {} mol/s is not positive",
                feed.n
            ),
        ));
    }
    let (mixture, ideal_gas) = feed.mixture()?;
    let (entropy, _) = ps_flash::entropy_at(&mixture, &ideal_gas, feed.t, feed.p, &feed.z)?;
    let isentropic = ps_flash::ps_flash(
        &mixture,
        &ideal_gas,
        outlet_pressure,
        joules_per_mole_kelvin(entropy),
        &feed.z,
    )?;
    let (h_isentropic, _) = ph_flash::enthalpy_at(
        &mixture,
        &ideal_gas,
        isentropic.temperature,
        outlet_pressure,
        &feed.z,
    )?;

    let dh = power.value / feed.n;
    let isentropic_step = h_isentropic - feed.h.value;
    let efficiency = if isentropic_step == 0.0 {
        None
    } else {
        Some(dh / isentropic_step)
    };
    let outlet = Stream::from_ph(
        feed.components.clone(),
        feed.z.clone(),
        feed.n,
        outlet_pressure,
        joules_per_mole(feed.h.value + dh),
    )?;
    Ok((outlet, efficiency))
}

/// `solveEfficiency`: the efficiency that lands the outlet at a stated temperature.
///
/// **The class's secant, with its first-step and fallback rules.** A run at the starting
/// efficiency gives the first residual; the first iteration always steps `+0.01`, later ones
/// take the secant slope when it is finite and non-singular and the `+-0.01` fallback when it
/// is not, and the efficiency is clamped by the setter's `[1e-6, 1]` window. The loop stops
/// when the temperature is within `0.01` K after at least three iterations, or at fifty.
///
/// An unreachable target therefore lands on the nearest clamp rather than converging - the
/// hot side is the isenthalpic floor and the cold side the reversible outlet.
///
/// # Errors
/// Whatever the flash refuses.
pub fn solve_efficiency(
    feed: &Stream,
    outlet_pressure: Pressure,
    target: ThermodynamicTemperature,
    start_efficiency: f64,
) -> Result<(Stream, f64)> {
    let mut efficiency = if start_efficiency.is_finite() {
        start_efficiency.clamp(MIN_EFFICIENCY, MAX_EFFICIENCY)
    } else {
        MIN_EFFICIENCY
    };
    let mut outlet = expander_replay(feed, outlet_pressure, efficiency)?;
    let mut function_old: f64 = 0.0;
    let mut delta_poly: f64 = 100.0;
    let mut iteration = 0;
    loop {
        iteration += 1;
        let function = outlet.t.value - target.value;
        let slope = if delta_poly.abs() > 1.0e-12
            && delta_poly.is_finite()
            && function.is_finite()
            && function_old.is_finite()
        {
            (function - function_old) / delta_poly
        } else {
            f64::NAN
        };
        let mut candidate = if !slope.is_finite() || slope.abs() < 1.0e-12 {
            efficiency + if function > 0.0 { 0.01 } else { -0.01 }
        } else {
            efficiency - function / slope
        };
        if iteration <= 1 {
            candidate = efficiency + 0.01;
        }
        if !candidate.is_finite() {
            candidate = efficiency + if function > 0.0 { 0.01 } else { -0.01 };
        }
        let old = efficiency;
        efficiency = candidate.clamp(MIN_EFFICIENCY, MAX_EFFICIENCY);
        delta_poly = efficiency - old;
        if delta_poly.abs() < 1.0e-12 {
            delta_poly = if function > 0.0 { 1.0e-12 } else { -1.0e-12 };
        }
        function_old = function;
        outlet = expander_replay(feed, outlet_pressure, efficiency)?;
        let residual = outlet.t.value - target.value;
        if !(residual.abs() > EFFICIENCY_SOLVE_TOLERANCE || iteration < 3) || iteration >= 50 {
            return Ok((outlet, efficiency));
        }
    }
}

/// `CapacityConstraint.MAX_UTILIZATION`, the ceiling on any utilization.
pub const MAX_UTILIZATION: f64 = 9.99;
/// `Expander.initializeCapacityConstraints`' hard band: the maximum is 1.1 x the rating.
pub const RECOVERED_POWER_HARD_MARGIN: f64 = 1.1;

/// `Expander`'s recovered-power constraint and its summaries.
///
/// The parent's speed, surge and stonewall constraints are **disabled** when the machine has
/// no chart, and the class's summaries skip a disabled constraint, so the only enabled entry
/// this surface can produce is `recoveredPower` - and only when a positive rating is set.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct ExpanderCapacity {
    /// The `recoveredPower` constraint's utilization, when a positive rating is set.
    pub recovered_power_utilization: Option<f64>,
    /// `getMaxUtilization`, 0.0 with no enabled constraint.
    pub max_utilization: f64,
    /// `isCapacityExceeded`: utilization above 1.0.
    pub capacity_exceeded: bool,
    /// `isHardLimitExceeded`: the recovered magnitude above 1.1 x the rating.
    pub hard_limit_exceeded: bool,
    /// `getBottleneckConstraint`'s name.
    pub bottleneck_constraint: Option<&'static str>,
}

/// The `recoveredPower` constraint, evaluated on the recovered shaft power.
///
/// The value is `abs(power)`, and a non-finite power answers 0.0, as the class's supplier
/// does. The rating is set to a non-positive value by the class with no constraint created
/// (`ratedRecoveredPower > 0.0`), so a non-positive input is the absent case.
#[must_use]
pub fn recovered_power_capacity(power: Power, rated: Option<Power>) -> ExpanderCapacity {
    let rated = rated.filter(|rated| rated.value > 0.0);
    let Some(rated) = rated else {
        return ExpanderCapacity::default();
    };
    let recovered = if power.value.is_finite() {
        power.value.abs()
    } else {
        0.0
    };
    let utilization = (recovered / rated.value).min(MAX_UTILIZATION);
    ExpanderCapacity {
        recovered_power_utilization: Some(utilization),
        max_utilization: if utilization > 0.0 { utilization } else { 0.0 },
        capacity_exceeded: utilization > 1.0,
        hard_limit_exceeded: recovered > RECOVERED_POWER_HARD_MARGIN * rated.value,
        bottleneck_constraint: (utilization > 0.0).then_some("recoveredPower"),
    }
}

#[cfg(test)]
mod capacity_tests {
    use super::*;
    use azoth_core::units::watts;

    const POWER_KW: f64 = -1.1055373484898778;

    fn power_kw() -> Power {
        watts(POWER_KW * 1000.0)
    }

    #[test]
    fn no_rating_is_the_empty_map() {
        let map = recovered_power_capacity(power_kw(), None);
        assert!(map.recovered_power_utilization.is_none());
        assert_eq!(map.max_utilization, 0.0);
        assert!(!map.capacity_exceeded);
        assert!(!map.hard_limit_exceeded);
        assert!(map.bottleneck_constraint.is_none());
    }

    #[test]
    fn a_rating_above_the_power_is_its_utilization() {
        let map = recovered_power_capacity(power_kw(), Some(watts(2000.0)));
        assert!(
            (map.recovered_power_utilization.expect("present") - 0.5527686742449389).abs() < 1e-12
        );
        assert!(!map.capacity_exceeded);
        assert!(!map.hard_limit_exceeded);
        assert_eq!(map.bottleneck_constraint, Some("recoveredPower"));
    }

    #[test]
    fn a_rating_between_the_power_and_its_hard_band_is_exceeded_but_not_hard() {
        let map = recovered_power_capacity(power_kw(), Some(watts(1200.0)));
        assert!(
            (map.recovered_power_utilization.expect("present") - 0.9212811237415649).abs() < 1e-12
        );
        assert!(!map.capacity_exceeded);
        assert!(!map.hard_limit_exceeded);
    }

    #[test]
    fn a_rating_below_the_power_trips_both_flags() {
        let map = recovered_power_capacity(power_kw(), Some(watts(1000.0)));
        assert!(
            (map.recovered_power_utilization.expect("present") - 1.1055373484898778).abs() < 1e-12
        );
        assert!(map.capacity_exceeded);
        assert!(map.hard_limit_exceeded);
    }

    #[test]
    fn a_non_positive_rating_is_unset() {
        assert_eq!(
            recovered_power_capacity(power_kw(), Some(watts(0.0))),
            ExpanderCapacity::default()
        );
    }

    #[test]
    fn a_non_finite_power_answers_zero() {
        let map = recovered_power_capacity(watts(f64::NAN), Some(watts(1000.0)));
        assert_eq!(map.recovered_power_utilization, Some(0.0));
        assert!(map.bottleneck_constraint.is_none());
    }
}

#[cfg(test)]
mod power_modes_tests {
    use super::*;
    use azoth_core::units::{kelvins, watts};

    fn feed() -> Stream {
        Stream::from_pt(
            vec![
                "methane".to_string(),
                "ethane".to_string(),
                "propane".to_string(),
                "n-butane".to_string(),
            ],
            vec![0.82, 0.10, 0.05, 0.03],
            690.618_418_406_472_6,
            pascals(9_000_000.0),
            kelvins(298.15),
        )
        .expect("the captured state")
    }

    fn close(actual: f64, expected: f64, tolerance: f64, what: &str) {
        let error = (actual - expected).abs() / expected.abs().max(f64::MIN_POSITIVE);
        assert!(
            error <= tolerance,
            "{what}: {actual} against {expected}, relative error {error:.3e} > {tolerance:.3e}"
        );
    }

    #[test]
    fn the_stated_power_is_the_efficiency_s() {
        let (outlet, efficiency) =
            expander_at_power(&feed(), pascals(3_000_000.0), watts(-985_780.729_484_916_9))
                .expect("the captured state");
        close(outlet.t.value, 246.508_123_701_218_98, 5.0e-4, "outlet_t");
        close(
            efficiency.expect("an efficiency"),
            0.800_000_006_064_420_7,
            5.0e-4,
            "efficiency",
        );
    }

    #[test]
    fn zero_power_is_the_isenthalpic_throttle() {
        let (outlet, efficiency) = expander_at_power(&feed(), pascals(3_000_000.0), watts(0.0))
            .expect("the captured state");
        close(outlet.t.value, 263.307_033_214_221_86, 5.0e-4, "outlet_t");
        close(
            efficiency.expect("an efficiency"),
            0.0,
            1.0e-12,
            "efficiency",
        );
    }

    #[test]
    fn a_positive_power_compresses_through_the_expander() {
        let (outlet, efficiency) =
            expander_at_power(&feed(), pascals(3_000_000.0), watts(98_578.068_410_488_08))
                .expect("the captured state");
        close(outlet.t.value, 264.938_253_580_763_96, 5.0e-4, "outlet_t");
        close(
            efficiency.expect("an efficiency"),
            -0.080_000_000_606_442_08,
            5.0e-4,
            "efficiency",
        );
    }

    #[test]
    fn the_outlet_temperature_solve_finds_two_hundred_and_fifty() {
        let (outlet, efficiency) =
            solve_efficiency(&feed(), pascals(3_000_000.0), kelvins(250.0), 0.8)
                .expect("the captured state");
        close(efficiency, 0.635_490_578_922_437_7, 5.0e-4, "efficiency");
        close(outlet.t.value, 250.000_002_665_501_09, 5.0e-4, "outlet_t");
    }

    #[test]
    fn the_unreachable_hot_target_lands_on_the_lower_clamp() {
        let (outlet, efficiency) =
            solve_efficiency(&feed(), pascals(3_000_000.0), kelvins(280.0), 0.8)
                .expect("the captured state");
        close(efficiency, MIN_EFFICIENCY, 5.0e-4, "efficiency");
        close(outlet.t.value, 263.307_012_738_515_37, 5.0e-4, "outlet_t");
    }

    #[test]
    fn the_unreachable_cold_target_lands_on_the_upper_clamp() {
        let (outlet, efficiency) =
            solve_efficiency(&feed(), pascals(3_000_000.0), kelvins(240.0), 0.8)
                .expect("the captured state");
        close(efficiency, 1.0, 5.0e-4, "efficiency");
        close(outlet.t.value, 242.266_840_019_28, 5.0e-4, "outlet_t");
    }
}
