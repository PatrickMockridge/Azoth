//! `unit_ops.heater` - a stream moved to a stated temperature or duty.

use azoth_core::units::{
    Power, Pressure, ThermodynamicTemperature, joules_per_mole, pascals, watts,
};
use azoth_core::{AzothError, Result};

use crate::stream::Stream;

/// What a heater did: the outlet record, and the duty it moved.
pub struct HeaterOutcome {
    /// The outlet record: the inlet's flow and composition at the new state.
    pub outlet: Stream,
    /// The duty it moved, W. Negative is heat removed, which is what a cooler is.
    pub duty: Power,
}

/// Heat or cool a stream to a stated temperature, or by a stated duty.
///
/// **The duty is in W and the enthalpy is per mole, and they meet because a NeqSim process
/// fluid's moles are one second of flow.** `Heater.run` adds `getDuty()` to
/// `system.getEnthalpy()` directly - W to J - which balances only while the fluid's total
/// moles equal its mol/s rate. Here the conversion is explicit, `duty / n` being the shift
/// in J/mol, so the convention is stated rather than inherited.
///
/// **The reported duty is the enthalpy the outlet reached, not the one the caller stated.**
/// `run` overwrites its own `energyInput` with `newH - oldH` after the flash, whatever branch
/// it took, so a temperature-specified heater still reports a duty - the state's, computed
/// from what the flash did.
///
/// **A temperature and a duty together are refused, and that is a measurement rather than a
/// preference.** `setOutletTemperature` sets `setTemperature` and clears `setEnergyInput`;
/// `setEnergyInput` (and so `setDuty`) sets `setEnergyInput` and clears `setTemperature`. The
/// pair is therefore *order-dependent* - the last setter called decides - which is a state
/// the declaration cannot express. Refusing it keeps this model a function of its declared
/// inputs; picking one would be inventing a rule the class does not have.
///
/// **With neither, the drop is isothermal.** `run`'s else branch is `T_in + dT` with `dT`
/// defaulting to zero, so a heater that specifies nothing still moves the pressure.
///
/// # Errors
/// [`azoth_core::AzothError::InvalidInput`] if a temperature and a duty are both given, if the
/// pressure drop would leave a non-positive pressure, or if a duty is asked of a stream
/// carrying no flow, and whatever the flash refuses.
pub fn heater(
    feed: &Stream,
    outlet_temperature: Option<ThermodynamicTemperature>,
    duty: Option<Power>,
    pressure_drop: Option<Pressure>,
) -> Result<HeaterOutcome> {
    if outlet_temperature.is_some() && duty.is_some() {
        return Err(AzothError::invalid_input(
            "duty",
            "a heater holds a temperature specification or a duty specification and not \
             both: `Heater.setOutletTemperature` clears the duty flag and `setDuty` clears \
             the temperature flag, so which one the class runs is the order they were set in \
             and not a function of these inputs",
        ));
    }
    let p_out = match pressure_drop {
        Some(drop) => pascals(feed.p.value - drop.value),
        None => feed.p,
    };
    if p_out.value <= 0.0 {
        return Err(AzothError::invalid_input(
            "pressure_drop",
            format!(
                "a pressure drop of {} Pa takes a {} Pa inlet to {} Pa, which is not a pressure",
                pressure_drop.map_or(0.0, |drop| drop.value),
                feed.p.value,
                p_out.value
            ),
        ));
    }

    let outlet = match (outlet_temperature, duty) {
        (Some(t), _) => Stream::from_pt(feed.components.clone(), feed.z.clone(), feed.n, p_out, t)?,
        (None, Some(power)) => {
            if feed.n == 0.0 {
                return Err(AzothError::invalid_input(
                    "duty",
                    "a duty is moved per mole of flow, and this stream carries none",
                ));
            }
            Stream::from_ph(
                feed.components.clone(),
                feed.z.clone(),
                feed.n,
                p_out,
                joules_per_mole(feed.h.value + power.value / feed.n),
            )?
        }
        (None, None) => Stream::from_pt(
            feed.components.clone(),
            feed.z.clone(),
            feed.n,
            p_out,
            feed.t,
        )?,
    };

    Ok(HeaterOutcome {
        duty: watts(feed.n * (outlet.h.value - feed.h.value)),
        outlet,
    })
}

/// `CapacityConstraint.MAX_UTILIZATION`, the ceiling on any utilization.
pub const MAX_UTILIZATION: f64 = 9.99;

/// The two outlet-temperature constraints the class's limits create, and their summaries.
///
/// Both are HARD. The maximum's utilization is `temperature / limit`; the minimum's is the
/// inverse, `limit / temperature`. The per-constraint fields are `None` when the limit was not
/// set, because `getCapacityConstraints` adds a constraint only for a set limit.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct HeaterCapacity {
    /// `maxOutletTemperature`'s utilization, when the maximum limit is set.
    pub max_outlet_temperature_utilization: Option<f64>,
    /// `minOutletTemperature`'s utilization, when the minimum limit is positive.
    pub min_outlet_temperature_utilization: Option<f64>,
    /// `getMaxUtilization`, 0.0 on an empty map.
    pub max_utilization: f64,
    /// `isCapacityExceeded`: any constraint's utilization above 1.0.
    pub capacity_exceeded: bool,
    /// `isHardLimitExceeded`: only the minimum can trip it - see [`capacity`].
    pub hard_limit_exceeded: bool,
    /// `getBottleneckConstraint`'s class-side name.
    pub bottleneck_constraint: Option<&'static str>,
}

/// `Heater`'s outlet-temperature capacity map, evaluated on the outlet temperature.
///
/// **The maximum's hard-limit check is a no-op in the class.** `isHardLimitExceeded` reads
/// `maxValue`, while `setMaxOutletTemperature` sets `designValue`; `maxValue` stays
/// `Double.MAX_VALUE`, so only the minimum constraint (whose `minValue` is set) can trip it.
///
/// **The limits are the class's sentinels.** A maximum of `f64::MAX` means "not set", and a
/// minimum of zero or less is not a constraint at all: `hasMinOutletTemperatureLimit` is
/// `min > 0.0`. Supplying either setter clears the map; the getters then initialize it with
/// whichever constraints are present, maximum first.
#[must_use]
pub fn capacity(
    outlet_temperature: ThermodynamicTemperature,
    max_outlet_temperature: Option<ThermodynamicTemperature>,
    min_outlet_temperature: Option<ThermodynamicTemperature>,
) -> HeaterCapacity {
    // numerics-ok: the class's unset sentinel, copied from `Double.MAX_VALUE`.
    let max_limit = max_outlet_temperature.filter(|limit| limit.value != f64::MAX);
    let min_limit = min_outlet_temperature.filter(|limit| limit.value > 0.0);

    let max_utilization = max_limit.map(|limit| {
        // `CapacityConstraint.getUtilization` answers zero for a non-positive design value.
        if limit.value <= 0.0 {
            0.0
        } else {
            (outlet_temperature.value / limit.value).min(MAX_UTILIZATION)
        }
    });
    let min_utilization = min_limit.map(|limit| {
        if outlet_temperature.value <= 0.0 {
            MAX_UTILIZATION
        } else {
            (limit.value / outlet_temperature.value).min(MAX_UTILIZATION)
        }
    });

    // The map is insertion-ordered: the maximum first, the minimum second, and the bottleneck
    // is the last one strictly above the running maximum, which starts at zero.
    let mut max_utilization_value = 0.0;
    let mut bottleneck_constraint = None;
    if let Some(value) = max_utilization {
        if value > max_utilization_value {
            max_utilization_value = value;
            bottleneck_constraint = Some("maxOutletTemperature");
        }
    }
    if let Some(value) = min_utilization {
        if value > max_utilization_value {
            max_utilization_value = value;
            bottleneck_constraint = Some("minOutletTemperature");
        }
    }

    HeaterCapacity {
        max_outlet_temperature_utilization: max_utilization,
        min_outlet_temperature_utilization: min_utilization,
        max_utilization: max_utilization_value,
        capacity_exceeded: max_utilization.is_some_and(|value| value > 1.0)
            || min_utilization.is_some_and(|value| value > 1.0),
        hard_limit_exceeded: min_limit.is_some_and(|limit| outlet_temperature.value < limit.value),
        bottleneck_constraint,
    }
}

/// `HeatExchangerInstrumentDesign.HeatExchangerType`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum InstrumentType {
    /// An electric resistance heater; the class detects this for a `Heater`.
    ElectricHeater,
    /// An air cooler's fan instrumentation.
    AirCooler,
    /// A shell-and-tube exchanger's utility-side instrumentation.
    ShellAndTube,
}

impl InstrumentType {
    /// The class's string spelling.
    #[must_use]
    pub fn parse(spelling: &str) -> Option<Self> {
        Some(match spelling {
            "electric_heater" => Self::ElectricHeater,
            "air_cooler" => Self::AirCooler,
            "shell_and_tube" => Self::ShellAndTube,
            _ => return None,
        })
    }
}

/// One `InstrumentSpecification` row, before the list assigns its tag.
#[derive(Debug, Clone, PartialEq)]
pub struct InstrumentRow {
    /// The ISA symbol, e.g. `TT`.
    pub isa_symbol: String,
    /// The service description.
    pub service: String,
    /// The range's lower bound, in `range_unit`.
    pub range_min: f64,
    /// The range's upper bound, in `range_unit`.
    pub range_max: f64,
    /// The range's unit; absent for a digital trip, whose constructor leaves it null.
    pub range_unit: Option<String>,
    /// `AI`, `AO`, `DI` or `DO`.
    pub io_type: String,
    /// The SIL rating; zero for a non-safety instrument.
    pub sil_rating: u32,
    /// Whether the row is safety-related.
    pub safety_related: bool,
}

#[allow(clippy::too_many_arguments)] // one parameter per column of the class's row
fn row(
    isa_symbol: &str,
    service: &str,
    range_min: f64,
    range_max: f64,
    range_unit: Option<&str>,
    io_type: &str,
    sil_rating: u32,
    safety_related: bool,
) -> InstrumentRow {
    InstrumentRow {
        isa_symbol: isa_symbol.to_string(),
        service: service.to_string(),
        range_min,
        range_max,
        range_unit: range_unit.map(str::to_owned),
        io_type: io_type.to_string(),
        sil_rating,
        safety_related,
    }
}

/// `HeatExchangerInstrumentDesign.calcDesign`'s list, in the class's insertion order.
///
/// **The list is static.** The five process instruments and the overtemperature trip are
/// common; the type adds the element monitor, the utility side's three, or the fan's three.
/// The tag is `SYMBOL-(1000 + position)`, which the list assigns on insertion.
#[must_use]
pub fn instruments(
    kind: InstrumentType,
    include_safety: bool,
    sil_rating: u32,
) -> Vec<InstrumentRow> {
    let mut rows = vec![
        row(
            "TT",
            "Process Inlet Temperature",
            -50.0,
            400.0,
            Some("degC"),
            "AI",
            0,
            false,
        ),
        row(
            "TT",
            "Process Outlet Temperature",
            -50.0,
            400.0,
            Some("degC"),
            "AI",
            0,
            false,
        ),
        row(
            "PT",
            "Process Inlet Pressure",
            0.0,
            200.0,
            Some("bara"),
            "AI",
            0,
            false,
        ),
        row(
            "PT",
            "Process Outlet Pressure",
            0.0,
            200.0,
            Some("bara"),
            "AI",
            0,
            false,
        ),
        row(
            "PDT",
            "Process Side dP",
            0.0,
            10.0,
            Some("bar"),
            "AI",
            0,
            false,
        ),
    ];
    if include_safety {
        rows.push(row(
            "TSHH",
            "Overtemperature Trip",
            0.0,
            0.0,
            None,
            "DI",
            sil_rating,
            true,
        ));
    }
    match kind {
        InstrumentType::ElectricHeater => {
            rows.push(row(
                "TT",
                "Heating Element Temperature",
                0.0,
                600.0,
                Some("degC"),
                "AI",
                0,
                false,
            ));
            if include_safety {
                rows.push(row(
                    "TSHH",
                    "Element Overtemperature Trip",
                    0.0,
                    0.0,
                    None,
                    "DI",
                    sil_rating,
                    true,
                ));
            }
        }
        InstrumentType::ShellAndTube => {
            rows.push(row(
                "TT",
                "Utility Inlet Temperature",
                -50.0,
                400.0,
                Some("degC"),
                "AI",
                0,
                false,
            ));
            rows.push(row(
                "TT",
                "Utility Outlet Temperature",
                -50.0,
                400.0,
                Some("degC"),
                "AI",
                0,
                false,
            ));
            rows.push(row(
                "PT",
                "Utility Inlet Pressure",
                0.0,
                50.0,
                Some("bara"),
                "AI",
                0,
                false,
            ));
        }
        InstrumentType::AirCooler => {
            rows.push(row(
                "TT",
                "Ambient Air Temperature",
                -40.0,
                60.0,
                Some("degC"),
                "AI",
                0,
                false,
            ));
            rows.push(row(
                "ST",
                "Fan Speed",
                0.0,
                1800.0,
                Some("rpm"),
                "AI",
                0,
                false,
            ));
            rows.push(row(
                "VT",
                "Fan Vibration",
                0.0,
                100.0,
                Some("mm/s"),
                "AI",
                0,
                false,
            ));
        }
    }
    rows
}

#[cfg(test)]
mod capacity_tests {
    use super::*;
    use azoth_core::units::kelvins;

    const OUTLET_K: f64 = 380.0;

    fn report(max: Option<f64>, min: Option<f64>) -> HeaterCapacity {
        capacity(kelvins(OUTLET_K), max.map(kelvins), min.map(kelvins))
    }

    #[test]
    fn no_limits_is_the_empty_map() {
        let map = report(None, None);
        assert!(map.max_outlet_temperature_utilization.is_none());
        assert!(map.min_outlet_temperature_utilization.is_none());
        assert_eq!(map.max_utilization, 0.0);
        assert!(!map.capacity_exceeded);
        assert!(!map.hard_limit_exceeded);
        assert!(map.bottleneck_constraint.is_none());
    }

    #[test]
    fn the_maximum_is_temperature_over_its_limit() {
        let map = report(Some(400.0), None);
        assert_eq!(map.max_outlet_temperature_utilization, Some(0.95));
        assert_eq!(map.max_utilization, 0.95);
        assert_eq!(map.bottleneck_constraint, Some("maxOutletTemperature"));
        assert!(!map.capacity_exceeded);

        let exceeded = report(Some(350.0), None);
        assert!(exceeded.capacity_exceeded);
        // The class's `isHardLimitExceeded` reads `maxValue`, which this setter never sets.
        assert!(!exceeded.hard_limit_exceeded);
    }

    #[test]
    fn the_minimum_is_its_limit_over_temperature() {
        let map = report(None, Some(300.0));
        assert!(
            (map.min_outlet_temperature_utilization.expect("present") - 300.0 / OUTLET_K).abs()
                < 1e-12
        );
        assert_eq!(map.bottleneck_constraint, Some("minOutletTemperature"));
        assert!(!map.capacity_exceeded);

        let exceeded = report(None, Some(400.0));
        assert!(exceeded.capacity_exceeded);
        assert!(exceeded.hard_limit_exceeded);
    }

    #[test]
    fn a_minimum_of_zero_is_not_a_constraint() {
        let map = report(None, Some(0.0));
        assert!(map.min_outlet_temperature_utilization.is_none());
        assert_eq!(map.max_utilization, 0.0);
        assert!(map.bottleneck_constraint.is_none());
    }

    #[test]
    fn a_maximum_of_the_sentinel_is_unset() {
        // numerics-ok: the sentinel the class reads to mean "no bound".
        let map = report(Some(f64::MAX), None);
        assert!(map.max_outlet_temperature_utilization.is_none());
    }

    #[test]
    fn the_utilization_is_capped_at_9_99() {
        let map = report(None, Some(1.0e6));
        assert_eq!(
            map.min_outlet_temperature_utilization,
            Some(MAX_UTILIZATION)
        );
        assert_eq!(map.max_utilization, MAX_UTILIZATION);
    }

    #[test]
    fn both_constraints_leave_the_larger_utilization_as_the_bottleneck() {
        let map = report(Some(400.0), Some(300.0));
        assert_eq!(map.max_utilization, 0.95);
        assert_eq!(map.bottleneck_constraint, Some("maxOutletTemperature"));
        assert!(!map.capacity_exceeded);
    }
}

#[cfg(test)]
mod instrument_tests {
    use super::*;

    #[test]
    fn the_electric_heater_has_eight_with_safety_and_six_without() {
        let with = instruments(InstrumentType::ElectricHeater, true, 2);
        assert_eq!(with.len(), 8);
        assert_eq!(
            with.iter()
                .map(|row| row.isa_symbol.as_str())
                .collect::<Vec<_>>(),
            vec!["TT", "TT", "PT", "PT", "PDT", "TSHH", "TT", "TSHH"]
        );
        assert_eq!(with[5].sil_rating, 2);
        assert!(with[5].safety_related);
        assert_eq!(with[5].range_unit, None);

        let without = instruments(InstrumentType::ElectricHeater, false, 2);
        assert_eq!(without.len(), 6);
        assert!(without.iter().all(|row| !row.safety_related));
    }

    #[test]
    fn the_other_types_add_three_instruments_each() {
        let shell = instruments(InstrumentType::ShellAndTube, true, 2);
        assert_eq!(shell.len(), 9);
        assert_eq!(shell[6].service, "Utility Inlet Temperature");
        assert_eq!(shell[7].service, "Utility Outlet Temperature");
        assert_eq!(shell[8].service, "Utility Inlet Pressure");

        let air = instruments(InstrumentType::AirCooler, true, 2);
        assert_eq!(air.len(), 9);
        assert_eq!(air[7].isa_symbol, "ST");
        assert_eq!(air[8].isa_symbol, "VT");
    }

    #[test]
    fn the_type_parses_its_three_spellings() {
        assert_eq!(
            InstrumentType::parse("electric_heater"),
            Some(InstrumentType::ElectricHeater)
        );
        assert_eq!(
            InstrumentType::parse("shell_and_tube"),
            Some(InstrumentType::ShellAndTube)
        );
        assert_eq!(InstrumentType::parse("coalescer"), None);
    }
}
