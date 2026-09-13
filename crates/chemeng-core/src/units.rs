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

pub use uom::si::f64::{DynamicViscosity, Length, MassDensity, Pressure, Velocity};
pub use uom::si::{
    dynamic_viscosity::pascal_second, length::meter, mass_density::kilogram_per_cubic_meter,
    pressure::pascal, velocity::meter_per_second,
};

/// A length in metres.
#[must_use]
pub fn meters(value: f64) -> Length {
    Length::new::<meter>(value)
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
}
