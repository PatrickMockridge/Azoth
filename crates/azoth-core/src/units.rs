//! Units, and the rule that governs how they are used.
//!
//! `uom` quantities are used at the **public boundary** and nowhere else: functions take
//! dimensioned quantities in and results carry them out, and internally a calculation
//! extracts the SI base value with `.value` and works in plain `f64`. Dividing two
//! same-dimension quantities in `uom` does not yield a plain scalar, so a Reynolds number
//! or a friction factor expressed as a ratio would need type-level ceremony for what is
//! numerically a division.
//!
//! Public functions are typed against fixed `SI<f64>` quantities rather than being
//! generic over `U: Units`: callers convert at the boundary either way, and fixed-SI
//! still makes `Length + Time` a type error.
//!
//! Dimensionless quantities - Reynolds number, relative roughness, friction factor,
//! resistance coefficient - are plain `f64`, here and in the Python API.

pub use uom::si::f64::{
    AmountOfSubstance, Area, DiffusionCoefficient, DynamicViscosity, ElectricCharge, Energy,
    HeatTransfer, Length, Mass, MassDensity, MassRate, Molality, MolarEnergy, MolarHeatCapacity,
    MolarMass, MolarVolume, Power, Pressure, SpecificHeatCapacity, SurfaceTension,
    TemperatureInterval, ThermalConductance, ThermalConductivity, ThermodynamicTemperature,
    Velocity, VolumeRate,
};
pub use uom::si::{
    acceleration::standard_gravity, amount_of_substance::kilomole, amount_of_substance::mole,
    area::square_meter, diffusion_coefficient::square_meter_per_second,
    dynamic_viscosity::pascal_second, electric_charge::coulomb, energy::btu_it, energy::joule,
    energy::kilojoule, energy::megajoule, heat_transfer::watt_per_square_meter_kelvin,
    length::angstrom, length::centimeter, length::foot, length::inch, length::meter,
    length::millimeter, mass::kilogram, mass::ton, mass_density::kilogram_per_cubic_meter,
    mass_rate::kilogram_per_hour, mass_rate::kilogram_per_second, mass_rate::ton_per_hour,
    molality::mole_per_kilogram, molar_energy::joule_per_mole, molar_energy::kilojoule_per_mole,
    molar_heat_capacity::joule_per_kelvin_mole, molar_mass::kilogram_per_mole,
    molar_volume::cubic_meter_per_mole, power::kilowatt, power::megawatt, power::watt,
    pressure::atmosphere, pressure::bar, pressure::kilopascal, pressure::megapascal,
    pressure::pascal, specific_heat_capacity::joule_per_kilogram_kelvin,
    specific_heat_capacity::kilojoule_per_kilogram_kelvin, surface_tension::newton_per_meter,
    temperature_interval::degree_fahrenheit, temperature_interval::kelvin as kelvin_interval,
    thermal_conductance::watt_per_kelvin, thermal_conductivity::watt_per_meter_kelvin,
    thermodynamic_temperature::kelvin, time::hour, time::minute, velocity::foot_per_second,
    velocity::meter_per_second, volume_rate::cubic_meter_per_hour,
    volume_rate::cubic_meter_per_second, volume_rate::liter_per_minute,
};

/// A length in metres.
#[must_use]
pub fn meters(value: f64) -> Length {
    Length::new::<meter>(value)
}

/// A length in millimetres.
///
/// One of the two units in the vocabulary that are not their own SI base unit, which is
/// why it is worth having explicitly: `.value` is still metres, so a caller building a
/// length from millimetres cannot accidentally work in them. Pipe diameters are
/// conventionally quoted in millimetres, so this is the constructor the next hydraulics
/// calcs will reach for.
#[must_use]
pub fn millimeters(value: f64) -> Length {
    Length::new::<millimeter>(value)
}

/// A length in ångström.
///
/// `.value` is still metres. The electrolyte models carry ion diameters in ångström -
/// NeqSim's `ComponentDesmukhMather` multiplies its own column by `1e-10` at the point of
/// use - and the conversion belongs here rather than at each read.
#[must_use]
pub fn angstroms(value: f64) -> Length {
    Length::new::<angstrom>(value)
}

/// A charge in coulombs.
///
/// SI has no base dimension for charge: the ampere is base and a coulomb is `A s`, which
/// is why this quantity's exponents are the `T` and `I` slots and not a slot of its own.
#[must_use]
pub fn coulombs(value: f64) -> ElectricCharge {
    ElectricCharge::new::<coulomb>(value)
}

/// A velocity in metres per second.
#[must_use]
pub fn meters_per_second(value: f64) -> Velocity {
    Velocity::new::<meter_per_second>(value)
}

/// A diffusion coefficient in square metres per second.
#[must_use]
pub fn square_meters_per_second(value: f64) -> DiffusionCoefficient {
    DiffusionCoefficient::new::<square_meter_per_second>(value)
}

/// A mass density in kilograms per cubic metre.
#[must_use]
pub fn kilograms_per_cubic_meter(value: f64) -> MassDensity {
    MassDensity::new::<kilogram_per_cubic_meter>(value)
}

/// A mass in kilograms.
#[must_use]
pub fn kilograms(value: f64) -> Mass {
    Mass::new::<kilogram>(value)
}

/// An energy in joules, which is a *total* and not a molar quantity.
#[must_use]
pub fn joules(value: f64) -> Energy {
    Energy::new::<joule>(value)
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

/// A thermal conductance in watts per kelvin.
///
/// `UA`, the heat exchanger rating's whole parameter: a heat-transfer coefficient times
/// an area, which is the power the exchanger moves per kelvin of driving temperature.
#[must_use]
pub fn watts_per_kelvin(value: f64) -> ThermalConductance {
    ThermalConductance::new::<watt_per_kelvin>(value)
}

/// A surface tension in newtons per metre.
#[must_use]
pub fn newtons_per_meter(value: f64) -> SurfaceTension {
    SurfaceTension::new::<newton_per_meter>(value)
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

/// A molality in moles per kilogram of solvent.
///
/// `.value` is still mol/kg. This is the scale every activity-coefficient phase works in,
/// and it is per kilogram of *solvent* rather than of solution - a distinction the number
/// cannot carry and the models' own documentation has to.
/// An amount of substance in moles.
///
/// `.value` is the mole number itself. **Added for `eos.hydrate_inhibitor_concentration`**,
/// whose answer is an absolute amount: NeqSim's secant adds moles to a system and reports the
/// inventory it reached, so a model taking a normalised composition would reproduce the
/// equation and not the path.
#[must_use]
pub fn moles(value: f64) -> AmountOfSubstance {
    AmountOfSubstance::new::<mole>(value)
}

/// A molality in moles per kilogram of solvent.
#[must_use]
pub fn moles_per_kilogram(value: f64) -> Molality {
    Molality::new::<mole_per_kilogram>(value)
}

/// A molar mass in kilograms per mole.
#[must_use]
pub fn kilograms_per_mole(value: f64) -> MolarMass {
    MolarMass::new::<kilogram_per_mole>(value)
}

/* --- the engineering units -------------------------------------------------
 *
 * The units a process engineer reads a flowsheet in: bar, psi, kW, Btu, lb/h, a
 * column's reflux in kmol/h. Every one of them is a *scale* of the SI base its
 * dimension already has, so nothing here is a new dimension - and that is the
 * property the generated compile-time assertions check, one unit at a time.
 *
 * **Nothing below is a measured factor.** Where uom carries the unit, the
 * constructor is uom's own definition, named so the vocabulary has something to
 * point at. Where it does not, or where uom's literal is rounded - its imperial
 * units are six or seven significant figures, and `psi` and `hp` are both out by
 * nearly 1e-7 - the constructor is the unit's *definition* assembled from parts
 * uom does define exactly: a pound-force is a pound times standard gravity, and a
 * mechanical horsepower is 550 foot-pounds-force per second. The check that any
 * of this is right is `python/tests/test_units_cross_library.py`, which holds
 * every one of these to `pint`'s own answer at 1e-15 relative - a tolerance the
 * rounded literals cannot pass, which is how they were found.
 */

/// One of a unit, as the SI base magnitude the quantity is stored in.
///
/// The two constructors below need a *part* of a definition rather than a whole
/// quantity - a pound-force is a pound times standard gravity, and a horsepower
/// is 550 of those per second - and `uom` exposes no other way to ask a unit for
/// its coefficient without naming the quantity it measures.
fn si_of<U: uom::Conversion<f64, T = f64>>() -> f64 {
    U::coefficient()
}

/// The international avoirdupois pound in kilograms: `0.453 592 37` exactly, which
/// the 1959 agreement fixed and no measurement decides.
///
/// **This number is here because `uom`'s is wrong at the seventh figure.** uom
/// carries `mass::pound` as `4.535_924_E-1`, 6.6e-8 from the definition, and it
/// carries `volume::cubic_foot` and `volume::gallon` rounded the same way - so
/// `psi`, `hp`, `lb/ft**3`, `gpm`, `Btu/(lb*degF)` and `ft**3/min` are all out by
/// that much if they are built on uom's. The cross-library check
/// (`python/tests/test_units_cross_library.py`) compares every declared unit
/// against `pint` at 1e-15 relative, which is how this was found and is the only
/// reason to trust the replacement: `pint` reads the agreement, uom reads a
/// rounding of it.
///
/// The rest of a pound-derived definition needs no number: a pound-force is this
/// times `standard_gravity`, which uom carries exactly, and a foot, an inch and a
/// Fahrenheit *interval* are exact in uom as well. This constant and the `550`
/// below are the only two quantities in this module that are written out, and both
/// are definitions rather than conversions.
const POUND_KILOGRAMS: f64 = 0.453_592_37;

/// A pressure in bars.
#[must_use]
pub fn bars(value: f64) -> Pressure {
    Pressure::new::<bar>(value)
}

/// A pressure in kilopascals.
#[must_use]
pub fn kilopascals(value: f64) -> Pressure {
    Pressure::new::<kilopascal>(value)
}

/// A pressure in megapascals.
#[must_use]
pub fn megapascals(value: f64) -> Pressure {
    Pressure::new::<megapascal>(value)
}

/// A pressure in standard atmospheres: 101 325 Pa by definition.
#[must_use]
pub fn atmospheres(value: f64) -> Pressure {
    Pressure::new::<atmosphere>(value)
}

/// A pressure in pounds-force per square inch.
///
/// Not uom's `pressure::psi` (`6.894_757_E3`, a rounded literal) and not its
/// `pound_force_per_square_inch` either, which is built on uom's rounded pound. A
/// pound-force is [`POUND_KILOGRAMS`] times `standard_gravity` and an inch is
/// exact, so this is the definition with nothing rounded in it.
#[must_use]
pub fn pounds_per_square_inch(value: f64) -> Pressure {
    let pound_force = POUND_KILOGRAMS * si_of::<standard_gravity>();
    let square_inch = si_of::<inch>() * si_of::<inch>();
    Pressure::new::<pascal>(value * pound_force / square_inch)
}

/// A power in watts, as mechanical horsepower: 550 foot-pounds-force per second.
///
/// Not uom's `power::horsepower` (`7.456_999_E2`) or `foot_pound_per_second`
/// (`1.355_818`), both rounded - and a compressor's shaft power is a number a
/// purchaser reads off a datasheet.
#[must_use]
pub fn mechanical_horsepower(value: f64) -> Power {
    let foot_pound_force = si_of::<foot>() * POUND_KILOGRAMS * si_of::<standard_gravity>();
    Power::new::<watt>(value * 550.0 * foot_pound_force)
}

/// A mass flow rate in kilograms per hour.
#[must_use]
pub fn kilograms_per_hour(value: f64) -> MassRate {
    MassRate::new::<kilogram_per_hour>(value)
}

/// A mass flow rate in tonnes per hour.
#[must_use]
pub fn tonnes_per_hour(value: f64) -> MassRate {
    MassRate::new::<ton_per_hour>(value)
}

/// A mass flow rate in pounds per hour.
///
/// Built on [`POUND_KILOGRAMS`] rather than on uom's `pound_per_hour`, for the
/// reason that constant gives.
#[must_use]
pub fn pounds_per_hour(value: f64) -> MassRate {
    MassRate::new::<kilogram_per_second>(value * POUND_KILOGRAMS / si_of::<hour>())
}

/// A mass in tonnes.
#[must_use]
pub fn tonnes(value: f64) -> Mass {
    Mass::new::<ton>(value)
}

/// A mass in pounds.
#[must_use]
pub fn pounds(value: f64) -> Mass {
    Mass::new::<kilogram>(value * POUND_KILOGRAMS)
}

/// An energy in kilojoules.
#[must_use]
pub fn kilojoules(value: f64) -> Energy {
    Energy::new::<kilojoule>(value)
}

/// An energy in megajoules.
#[must_use]
pub fn megajoules(value: f64) -> Energy {
    Energy::new::<megajoule>(value)
}

/// An energy in British thermal units, International Table.
///
/// `uom` calls this `btu_it`, and its literal `1.055_056_E3` is what `pint`'s
/// `Btu` is to the last digit it carries, so this one is uom's own.
#[must_use]
pub fn british_thermal_units(value: f64) -> Energy {
    Energy::new::<btu_it>(value)
}

/// A power in kilowatts.
#[must_use]
pub fn kilowatts(value: f64) -> Power {
    Power::new::<kilowatt>(value)
}

/// A power in megawatts.
#[must_use]
pub fn megawatts(value: f64) -> Power {
    Power::new::<megawatt>(value)
}

/// A length in feet.
#[must_use]
pub fn feet(value: f64) -> Length {
    Length::new::<foot>(value)
}

/// A length in inches.
#[must_use]
pub fn inches(value: f64) -> Length {
    Length::new::<inch>(value)
}

/// A length in centimetres.
#[must_use]
pub fn centimeters(value: f64) -> Length {
    Length::new::<centimeter>(value)
}

/// A volumetric flow rate in cubic metres per hour.
#[must_use]
pub fn cubic_meters_per_hour(value: f64) -> VolumeRate {
    VolumeRate::new::<cubic_meter_per_hour>(value)
}

/// A volumetric flow rate in litres per minute.
#[must_use]
pub fn liters_per_minute(value: f64) -> VolumeRate {
    VolumeRate::new::<liter_per_minute>(value)
}

/// A volumetric flow rate in US gallons per minute.
///
/// A US liquid gallon is 231 cubic inches **by definition**, and an inch is exact
/// in uom, so the `231` here is the definition and not a factor - uom's own
/// `gallon` is a rounded literal for the same reason its pound is.
#[must_use]
pub fn gallons_per_minute(value: f64) -> VolumeRate {
    let gallon = 231.0 * si_of::<inch>().powi(3);
    VolumeRate::new::<cubic_meter_per_second>(value * gallon / si_of::<minute>())
}

/// A volumetric flow rate in cubic feet per minute.
///
/// uom's `cubic_foot` is `2.831_685_E-2`, 5.9e-8 from the cube of its own exact
/// foot, so this is built from the foot.
#[must_use]
pub fn cubic_feet_per_minute(value: f64) -> VolumeRate {
    let cubic_foot = si_of::<foot>().powi(3);
    VolumeRate::new::<cubic_meter_per_second>(value * cubic_foot / si_of::<minute>())
}

/// A velocity in feet per second.
#[must_use]
pub fn feet_per_second(value: f64) -> Velocity {
    Velocity::new::<foot_per_second>(value)
}

/// A mass density in pounds per cubic foot.
#[must_use]
pub fn pounds_per_cubic_foot(value: f64) -> MassDensity {
    MassDensity::new::<kilogram_per_cubic_meter>(value * POUND_KILOGRAMS / si_of::<foot>().powi(3))
}

/// A specific heat capacity in kilojoules per kilogram kelvin.
#[must_use]
pub fn kilojoules_per_kilogram_kelvin(value: f64) -> SpecificHeatCapacity {
    SpecificHeatCapacity::new::<kilojoule_per_kilogram_kelvin>(value)
}

/// A specific heat capacity in British thermal units per pound degree Fahrenheit.
///
/// The Fahrenheit here is an **interval** and not a temperature: a specific heat
/// capacity is energy per mass per degree of *difference*, and `uom`'s
/// `btu_it_per_pound_degree_fahrenheit` is built that way. An absolute °F would be
/// an offset unit, which is a different thing this crate does not carry - see the
/// note on [`kelvin_intervals`].
#[must_use]
pub fn british_thermal_units_per_pound_degree_fahrenheit(value: f64) -> SpecificHeatCapacity {
    // The Fahrenheit here is uom's *interval*, which is exactly five ninths of a
    // kelvin; only the pound in the denominator is uom's rounded one.
    let fahrenheit_interval = si_of::<degree_fahrenheit>();
    SpecificHeatCapacity::new::<joule_per_kilogram_kelvin>(
        value * si_of::<btu_it>() / (POUND_KILOGRAMS * fahrenheit_interval),
    )
}

/// A molar energy in kilojoules per mole.
#[must_use]
pub fn kilojoules_per_mole(value: f64) -> MolarEnergy {
    MolarEnergy::new::<kilojoule_per_mole>(value)
}

/// A dynamic viscosity in centipoise.
#[must_use]
pub fn centipoise(value: f64) -> DynamicViscosity {
    // Qualified, because the unit and this function are one name: the vocabulary's
    // `rust_ctor` is spelled after the unit, and a `use` of both would be a
    // redefinition rather than an import.
    DynamicViscosity::new::<uom::si::dynamic_viscosity::centipoise>(value)
}

/// An amount of substance in kilomoles.
#[must_use]
pub fn kilomoles(value: f64) -> AmountOfSubstance {
    AmountOfSubstance::new::<kilomole>(value)
}

/// A molar flow rate in kilomoles per hour, **as an SI base magnitude**.
///
/// The one constructor here that returns a plain `f64` rather than a quantity, and
/// the reason is uom's: it carries no molar-flow quantity at all, so there is no
/// `MolarFlow` type to build. The dimension's SI base is `mol/s`, and this is that
/// — uom's kilomole over uom's hour, so the thousand and the three thousand six
/// hundred are both read from the units library rather than written here.
#[must_use]
pub fn kilomoles_per_hour(value: f64) -> f64 {
    value * si_of::<kilomole>() / si_of::<hour>()
}

/// The canonical unit strings the spec schema permits.
///
/// Generated from `specs/vocabulary/vocabulary.toml`, which is the one
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
