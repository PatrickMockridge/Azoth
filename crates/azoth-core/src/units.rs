//! Units, and the rule that governs how they are used.
//!
//! # The boundary rule
//!
//! `uom` quantities are used at the **public boundary** and nowhere else:
//! functions take dimensioned quantities in, and results carry dimensioned
//! quantities out. Internally, calculations extract the SI base value with
//! `.value` and work in plain `f64`.
//!
//! That is a deliberate choice, not a shortcut, and the reason is that `uom`'s
//! compile-time dimension arithmetic is genuinely valuable at the boundary - it
//! makes `Length + Time` a type error - while being actively unhelpful in the
//! middle of an equation. Dividing two same-dimension quantities in `uom` does
//! not yield a plain scalar, so every Reynolds number, relative roughness and
//! friction factor would need type-level ceremony to express a ratio that is
//! numerically just a division.
//!
//! The cost of the rule is that dimensional correctness is enforced on the way
//! in and out but not step by step through an equation. The benefit is that the
//! bodies of the calculations are readable and identical in structure to the
//! published equations. For short equations transcribed directly from a
//! standard, matching the published form is the more valuable safety property:
//! a reader can check the code against the paper line by line.
//!
//! # Why `SI<f64>` and not a generic `U`
//!
//! Public functions are typed against fixed `SI<f64>` quantities rather than
//! being generic over `U: Units`. Full genericity requires threading
//! `U: Conversion<f64> + Copy` style bounds through every signature, which makes
//! them unreadable and buys nothing a caller wants: callers convert at the
//! boundary either way. Fixed-SI still catches unit errors at compile time,
//! because a caller cannot construct a `Length` from an unconverted number
//! without saying which unit it is in.
//!
//! # Dimensionless quantities
//!
//! Genuinely dimensionless quantities - Reynolds number, relative roughness,
//! friction factor, resistance coefficient - are plain `f64`, here and in the
//! Python API. They carry no unit to be safe about, and a newtype wrapper around
//! a ratio would be an abstraction with exactly one implementation.

pub use uom::si::f64::{
    Area, DynamicViscosity, HeatTransfer, Length, MassDensity, MassRate, MolarEnergy,
    MolarHeatCapacity, MolarMass, MolarVolume, Power, Pressure, SpecificHeatCapacity,
    TemperatureInterval, ThermalConductivity, ThermodynamicTemperature, Velocity, VolumeRate,
};
pub use uom::si::{
    area::square_meter, dynamic_viscosity::pascal_second,
    heat_transfer::watt_per_square_meter_kelvin, length::meter, length::millimeter,
    mass_density::kilogram_per_cubic_meter, mass_rate::kilogram_per_second,
    molar_energy::joule_per_mole, molar_heat_capacity::joule_per_kelvin_mole,
    molar_mass::kilogram_per_mole, molar_volume::cubic_meter_per_mole, power::watt,
    pressure::pascal, specific_heat_capacity::joule_per_kilogram_kelvin,
    temperature_interval::kelvin as kelvin_interval, thermal_conductivity::watt_per_meter_kelvin,
    thermodynamic_temperature::kelvin, velocity::meter_per_second,
    volume_rate::cubic_meter_per_second,
};

/// A length in metres.
#[must_use]
pub fn meters(value: f64) -> Length {
    Length::new::<meter>(value)
}

/// A length in millimetres.
///
/// The one unit in the vocabulary that is not its own SI base unit, which is why
/// it is worth having explicitly: `.value` is still metres, so a caller building a
/// length from millimetres cannot accidentally work in them. Pipe diameters are
/// conventionally quoted in millimetres, so this is the constructor the next
/// hydraulics calcs will reach for.
#[must_use]
pub fn millimeters(value: f64) -> Length {
    Length::new::<millimeter>(value)
}

/// A velocity in metres per second.
#[must_use]
pub fn meters_per_second(value: f64) -> Velocity {
    Velocity::new::<meter_per_second>(value)
}

/// A mass density in kilograms per cubic metre.
#[must_use]
pub fn kilograms_per_cubic_meter(value: f64) -> MassDensity {
    MassDensity::new::<kilogram_per_cubic_meter>(value)
}

/// A dynamic viscosity in pascal seconds.
#[must_use]
pub fn pascal_seconds(value: f64) -> DynamicViscosity {
    DynamicViscosity::new::<pascal_second>(value)
}

/// A pressure in pascals.
#[must_use]
pub fn pascals(value: f64) -> Pressure {
    Pressure::new::<pascal>(value)
}

/// A thermodynamic temperature in kelvin.
///
/// Added when the unit vocabulary was made checkable, not when a calc first
/// needed it: `K` had been a unit the schema permitted and `CANONICAL_UNITS`
/// knew about since the beginning, while this crate had no temperature type at
/// all. Nothing used it, so nothing noticed.
#[must_use]
pub fn kelvins(value: f64) -> ThermodynamicTemperature {
    ThermodynamicTemperature::new::<kelvin>(value)
}

/// A temperature *interval* in kelvin: a difference between two temperatures,
/// which is a different thing from [`kelvins`] and deliberately a different type.
///
/// `uom` separates them because they convert differently, and that difference is
/// the reason to keep the types apart at the boundary. A 30 K interval is a 30 degC
/// interval, but an absolute 30 K is -243.15 degC: treating a difference as an
/// absolute temperature silently adds 273.15, which is a plausible-looking wrong
/// number rather than an error.
///
/// A calc that means a difference therefore takes this type, and cannot be handed
/// the absolute one by accident. `conduction_plane_wall` is the first caller - its
/// `dT` is a difference across a wall, and it has no opinion about either face's
/// absolute temperature.
///
/// The vocabulary's `K` maps to [`kelvins`], the absolute one, because that is what
/// the unit name means on its own. A difference measured in kelvin is the same
/// number either way, so a spec declaring `K` for a difference converts to the same
/// magnitude through either type - the distinction is only enforceable in Rust,
/// where the caller has to choose.
#[must_use]
pub fn kelvin_intervals(value: f64) -> TemperatureInterval {
    TemperatureInterval::new::<kelvin_interval>(value)
}

/// An area in square metres.
#[must_use]
pub fn square_meters(value: f64) -> Area {
    Area::new::<square_meter>(value)
}

/// A volumetric flow rate in cubic metres per second.
#[must_use]
pub fn cubic_meters_per_second(value: f64) -> VolumeRate {
    VolumeRate::new::<cubic_meter_per_second>(value)
}

/// A mass flow rate in kilograms per second.
#[must_use]
pub fn kilograms_per_second(value: f64) -> MassRate {
    MassRate::new::<kilogram_per_second>(value)
}

/// A power in watts.
#[must_use]
pub fn watts(value: f64) -> Power {
    Power::new::<watt>(value)
}

/// A specific heat capacity in joules per kilogram kelvin.
#[must_use]
pub fn joules_per_kilogram_kelvin(value: f64) -> SpecificHeatCapacity {
    SpecificHeatCapacity::new::<joule_per_kilogram_kelvin>(value)
}

/// A thermal conductivity in watts per metre kelvin.
#[must_use]
pub fn watts_per_meter_kelvin(value: f64) -> ThermalConductivity {
    ThermalConductivity::new::<watt_per_meter_kelvin>(value)
}

/// A heat transfer coefficient in watts per square metre kelvin.
#[must_use]
pub fn watts_per_square_meter_kelvin(value: f64) -> HeatTransfer {
    HeatTransfer::new::<watt_per_square_meter_kelvin>(value)
}

/// A molar volume in cubic metres per mole.
///
/// The namespace's first dimensional quantity, and the one that takes an equation
/// of state from a compressibility factor to a volume.
#[must_use]
pub fn cubic_meters_per_mole(value: f64) -> MolarVolume {
    MolarVolume::new::<cubic_meter_per_mole>(value)
}

/// A molar energy in joules per mole.
///
/// Carries an enthalpy, and - as [`molar_heat_capacity`] explains - an entropy too.
/// The model layer is what makes it necessary: every kernel in `eos` returns
/// dimensionless departures, and the multiplication by `R*T` that turns one into
/// joules happens where `R` and `T` are, at the top.
#[must_use]
pub fn joules_per_mole(value: f64) -> MolarEnergy {
    MolarEnergy::new::<joule_per_mole>(value)
}

/// A molar heat capacity in joules per mole kelvin.
///
/// **This is also the carrier for a molar *entropy*.** `uom` has no
/// `MolarEntropy`, and it does not need one: the two are dimensionally identical -
/// `J/(mol*K)` either way - so a second quantity type would be a second name for one
/// dimension. The rereading is worth a comment rather than a silent reuse, because a
/// reader who sees `MolarHeatCapacity` carrying an entropy should be able to find out
/// in one place why that is right.
///
/// Note the conversion path's name: `joule_per_kelvin_mole`, not the
/// `joule_per_mole_kelvin` the unit string suggests. That ordering is `uom`'s, and
/// getting it wrong is a compile error rather than a silent one, which is one of the
/// reasons the boundary uses `uom` at all.
#[must_use]
pub fn joules_per_mole_kelvin(value: f64) -> MolarHeatCapacity {
    MolarHeatCapacity::new::<joule_per_kelvin_mole>(value)
}

/// A molar mass in kilograms per mole.
#[must_use]
pub fn kilograms_per_mole(value: f64) -> MolarMass {
    MolarMass::new::<kilogram_per_mole>(value)
}

/// The canonical unit strings the spec schema permits.
///
/// Generated from `specs/vocabulary/vocabulary.yaml`, which is the one
/// hand-written source of this list - and of the dimension each name carries, the
/// conversion this crate performs for it, and the `pint` name it has on the Python
/// side. The list reaches Python through `azoth._core.unit_names`.
///
/// A name here is a claim that this crate has a *correct* conversion path for it -
/// see [`unit_vocab_gen::CONVERSION_PATHS`], and
/// `every_unit_name_has_a_conversion_path` below, which fails if a name is added
/// without one.
pub use crate::unit_vocab_gen::{CONVERSION_PATHS, SLOTS, UNIT_DIMENSIONS, UNIT_NAMES};

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn constructors_are_in_si_base() {
        // `.value` is always the SI base value, which is what the boundary rule
        // relies on. These assertions are what would fail first if a uom upgrade
        // changed the meaning of the stored value.
        assert_eq!(meters(2.5).value, 2.5);
        assert_eq!(meters_per_second(1.5).value, 1.5);
        assert_eq!(kilograms_per_cubic_meter(998.0).value, 998.0);
        assert_eq!(pascal_seconds(1.002e-3).value, 1.002e-3);
        assert_eq!(pascals(22455.0).value, 22455.0);
    }

    #[test]
    fn unit_conversion_happens_at_construction() {
        // A caller working in US customary converts once, at the boundary, and
        // everything downstream is SI. 1 ft = 0.3048 m exactly.
        let l = Length::new::<uom::si::length::foot>(1.0);
        assert!((l.value - 0.3048).abs() < 1e-15, "got {}", l.value);

        // The same physical state described two ways must be the same quantity.
        let a = meters(0.3048);
        let b = Length::new::<uom::si::length::foot>(1.0);
        assert!((a.value - b.value).abs() < 1e-15);
    }

    #[test]
    fn every_unit_name_has_a_conversion_path() {
        use std::collections::BTreeSet;

        let declared: BTreeSet<&str> = UNIT_NAMES.iter().copied().collect();
        let convertible: BTreeSet<&str> = CONVERSION_PATHS.iter().map(|(n, _)| *n).collect();

        let missing: Vec<_> = declared.difference(&convertible).collect();
        let extra: Vec<_> = convertible.difference(&declared).collect();
        assert!(
            missing.is_empty(),
            "unit name(s) {missing:?} are permitted by the vocabulary but this crate has \
             no conversion for them, so a spec could declare one and a calculation would \
             receive a number in the wrong unit"
        );
        assert!(
            extra.is_empty(),
            "conversion path(s) {extra:?} exist for unit name(s) not in UNIT_NAMES"
        );
    }

    #[test]
    fn every_dimension_has_the_same_width_as_the_slots() {
        // The exponent tuples and the slot names are two halves of one encoding,
        // and a tuple of the wrong length would be read against the wrong slots
        // rather than rejected - `mm` as `[1]` would be a length, and as
        // `[1, 0, 0, 0, 0, 0, 0, 0]` would be nonsense that still compared equal
        // to nothing.
        for (name, exponents) in UNIT_DIMENSIONS {
            assert_eq!(
                exponents.len(),
                SLOTS.len(),
                "{name}: {} exponent(s) for {} slot(s)",
                exponents.len(),
                SLOTS.len()
            );
        }
    }

    #[test]
    fn every_conversion_is_total_and_positive() {
        // A weak statement on purpose: it is the strongest one this side can make.
        //
        // What the conversion should *yield* is a number the units libraries
        // already know, so asserting a magnitude here would put a hand-typed
        // factor back into this repository - which is the defect the generated
        // table exists to remove. The check that a conversion yields the right
        // number is `python/tests/test_units_cross_library.py`, which compares
        // each of these against `pint`'s own answer for the same unit name.
        //
        // What is checkable here is that no entry is a stub: every conversion
        // returns a finite, positive magnitude, so a name added with nothing
        // behind it fails rather than sitting in the vocabulary looking live.
        for (name, convert) in CONVERSION_PATHS {
            let got = convert(1.0);
            assert!(
                got.is_finite() && got > 0.0,
                "{name}: converting 1.0 yielded {got}, which is not a positive finite \
                 magnitude"
            );
        }
    }

    #[test]
    fn unit_names_are_unique() {
        use std::collections::BTreeSet;
        let unique: BTreeSet<&str> = UNIT_NAMES.iter().copied().collect();
        assert_eq!(
            unique.len(),
            UNIT_NAMES.len(),
            "UNIT_NAMES contains a duplicate"
        );
    }
}
