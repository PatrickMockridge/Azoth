//! `ColumnInternalsDesigner.calculateTrayed`: the trayed internals tree, as a report.
//!
//! **The class the tray calculator is driven *by*.** It walks the column's trays, picks the
//! **controlling** one by the largest vapour mass flow, sizes a diameter off that tray, then runs
//! one calculator per tray at the diameter it resolved and sums them. Everything the class
//! publishes comes out of that second loop: the total pressure drop, the two percent-flood
//! extremes, the average efficiency and the per-tray results.
//!
//! **Three of the class's own lines are reproduced rather than corrected**, each measured in
//! `validation/neqsim/captures/process_internals_designer.tsv`:
//!
//! **The controlling tray's calculator is discarded.** It exists only to size the diameter from,
//! and it is built **without a relative volatility** while every tray in the summation loop
//! carries one - so the sizing runs at the class's own default `2.0` on trays that read `9.25`
//! and `20.0`.
//!
//! **The extraction's fallbacks are the class's.** `getTrayProperties` reads the tray's own
//! fluid, answers `1.0`/`800.0`/`0.001`/`0.02`/`2.0` when it is absent or single-phase, and on a
//! single-phase tray it takes **phase 0** as the vapour density and substitutes `600.0` for the
//! liquid's when that phase is a gas.
//!
//! **The mass flows are the tray's own two streams.** `getTrayFlows` reads `kg/hr` off
//! `getGasOutStream` and `getLiquidOutStream` and divides by 3600 - so they are not equal, which
//! the capture shows on every tray, and a port that read the tray's *system* would publish one
//! mass flow twice.

use azoth_core::units::{Length, Pressure, meters, pascals};
use azoth_core::{AzothError, Result};
use azoth_hydraulics::tray_hydraulics::{
    TrayHydraulicsState, size_column_diameter, tray_hydraulics,
};

use crate::kernels::distillation_column::{TrayProfile, tray_streams};
use crate::segment::phase::{liquid_view, vapour_view};

/// `ColumnInternalsDesigner`'s own constructor defaults, one per geometry input. **They are the
/// class's initialisers and not correlated values**, and `calcColumnInternals` sets none of
/// them: what a coupling pass designs against is these.
pub const DEFAULT_TRAY_SPACING_M: f64 = 0.6;
pub const DEFAULT_WEIR_HEIGHT_M: f64 = 0.05;
pub const DEFAULT_HOLE_DIAMETER_MM: f64 = 12.7;
pub const DEFAULT_HOLE_AREA_FRACTION: f64 = 0.1;
pub const DEFAULT_DOWNCOMMER_AREA_FRACTION: f64 = 0.1;
pub const DEFAULT_DESIGNER_FLOOD_FRACTION: f64 = 0.8;
/// The internals type the class's own constructor carries.
pub const DEFAULT_INTERNALS_TYPE: &str = "sieve";

/// The relative volatility the class's own default stands at, and the one the *sizing* branch
/// runs at because its calculator is never given the tray's.
pub const DEFAULT_RELATIVE_VOLATILITY: f64 = 2.0;
/// The liquid density a **single gas phase** substitutes for a missing liquid, `getTrayProperties`'s
/// own `600.0`.
pub const SINGLE_PHASE_GAS_LIQUID_DENSITY: f64 = 600.0;
/// The liquid density the class's extraction leaves in place on a single non-gas phase, which is
/// `getTrayProperties`'s own initialiser and not a correlated value.
pub const SINGLE_PHASE_LIQUID_DENSITY_KG_PER_M3: f64 = 800.0;
/// The liquid viscosity the same extraction falls back to where the phase answers nothing.
pub const DEFAULT_LIQUID_VISCOSITY_PA_S: f64 = 0.001;
/// The near-dry surface tension the extraction falls back to when the interface answers nothing.
pub const DESIGNER_SURFACE_TENSION_N_PER_M: f64 = 0.02;
/// The average the class answers when its loop walked no tray at all.
pub const EMPTY_AVERAGE_TRAY_EFFICIENCY: f64 = 0.65;
/// The column diameter that means "size it", `ColumnInternalsDesigner`'s own `-1.0`.
pub const UNSIZED_COLUMN_DIAMETER_M: f64 = -1.0;

/// The designer's own geometry: its setters, and the class's constructor defaults beside them.
#[derive(Debug, Clone)]
pub struct DesignerGeometry {
    /// `sieve`, `valve` or `bubble-cap`. The class's own constructor default is `sieve`.
    pub internals_type: String,
    /// The tray spacing, m. The class's own default is `0.6`.
    pub tray_spacing: Length,
    /// The weir height, m. The class's own default is `0.05`.
    pub weir_height: Length,
    /// The hole diameter, m - **the class's field is millimetres and so is the spec's unit**.
    pub hole_diameter: Length,
    /// The hole area over the active area. The class's own default is `0.1`.
    pub hole_area_fraction: f64,
    /// The downcomer area over the total. The class's own default is `0.1`.
    pub downcommer_area_fraction: f64,
    /// The fraction of flood a tray is sized to. The class's own default is `0.8`.
    pub design_flood_fraction: f64,
    /// A stated diameter, which **replaces the sizing branch entirely** at or below zero.
    pub column_diameter_override: Length,
}

/// One tray's own results, which the designer keeps as a list of calculators.
#[derive(Debug, Clone, PartialEq)]
pub struct DesignerTrayResult {
    /// The tray's load, in per cent of flood.
    pub percent_flood: f64,
    /// The tray's own total pressure drop, Pa.
    pub total_pressure_drop: Pressure,
    /// The tray's efficiency.
    pub tray_efficiency: f64,
    /// Whether this tray is inside its own design window.
    pub design_ok: bool,
}

/// `ColumnInternalsDesigner.calculateTrayed`'s whole answer.
#[derive(Debug, Clone, PartialEq)]
pub struct DesignerReport {
    /// The diameter the sizing resolved, or the stated override where one was.
    pub required_diameter: Length,
    /// **Which tray the diameter was sized from**: the one with the largest vapour mass flow,
    /// the first such tray where two are equal.
    pub controlling_tray_index: usize,
    /// Every tray's own verdict, and-ed together.
    pub design_ok: bool,
    /// The largest per-tray load, in per cent of flood.
    pub max_percent_flood: f64,
    /// **The smallest load above zero**, which is not the smallest - a tray that reads exactly
    /// zero is skipped rather than taken.
    pub min_percent_flood: f64,
    /// The mean of the per-tray efficiencies, or the class's own `0.65` with no tray walked.
    pub average_tray_efficiency: f64,
    /// The trays' pressure drops, summed.
    pub total_pressure_drop: Pressure,
    /// The same sum in millibars, which the class publishes beside it.
    pub total_pressure_drop_mbar: f64,
    /// One entry per tray, in the column's own order.
    pub trays: Vec<DesignerTrayResult>,
}

/// What `getTrayFlows` and `getTrayProperties` read off one tray.
#[derive(Debug, Clone)]
struct TrayReadings {
    vapor_mass_flow: f64,
    liquid_mass_flow: f64,
    vapor_density: f64,
    liquid_density: f64,
    liquid_viscosity: f64,
    surface_tension: f64,
    relative_volatility: f64,
}

/// The class's own two extraction helpers, at one tray.
///
/// **The fallbacks are the class's**: a tray whose two streams are both single-phase takes phase
/// zero's density as the vapour's, and **`600.0`** as the liquid's where that phase is a gas.
fn readings(tray: &TrayProfile, components: &[String]) -> Result<TrayReadings> {
    let (vapour_stream, liquid_stream) = tray_streams(tray, components)?;
    // **The phase each stream carries, and not what a re-flash of it answers.** A pinned end is
    // saturated by construction, so a re-flash of either side lands on the knife edge - see
    // `vapour_view` and `liquid_view`, whose doc says what that cost before they existed.
    let vapour = vapour_view(&vapour_stream)?;
    let liquid = liquid_view(&liquid_stream)?;

    let vapor_mass_flow = tray.gas_n * vapour.molar_mass;
    let liquid_mass_flow = tray.liquid_n * liquid.molar_mass;

    // **The tray's own two traffics decide, and they are the class's `getNumberOfPhases()`**:
    // a tray with no vapour or no liquid has one phase, and one with both has two - which is not
    // the same question as whether either *rebuilt* stream happens to flash into one. Reading it
    // off the two picks is wrong and was measured wrong: the reboiler's vapour stream and its
    // liquid stream both flash to a single phase here, so a kind comparison calls the tray
    // single-phase and substitutes `600.0` for a liquid that has its own density of `435.7`.
    let single_phase = tray.gas_n <= 0.0 || tray.liquid_n <= 0.0;
    let (vapor_density, liquid_density) = if single_phase {
        let density = vapour.density;
        let substituted = if matches!(vapour.kind, azoth_eos::phase_transport::PhaseKind::Gas) {
            SINGLE_PHASE_GAS_LIQUID_DENSITY
        } else {
            // The class leaves its `800.0` initialiser in place when phase zero is not a gas.
            SINGLE_PHASE_LIQUID_DENSITY_KG_PER_M3
        };
        (density, substituted)
    } else {
        (vapour.density, liquid.density)
    };
    let relative_volatility = if single_phase {
        DEFAULT_RELATIVE_VOLATILITY
    } else {
        relative_volatility(&vapour.z, &liquid.z)
    };

    // `getTrayProperties`'s own guard on the viscosity, which is a floor and not a substitution.
    let viscosity = liquid.transport.mu.value;
    let liquid_viscosity = if viscosity <= 0.0 {
        DEFAULT_LIQUID_VISCOSITY_PA_S
    } else {
        viscosity
    };

    Ok(TrayReadings {
        vapor_mass_flow,
        liquid_mass_flow,
        vapor_density,
        liquid_density,
        liquid_viscosity,
        // **The interface route answers nothing on these states** - every tray of the capture
        // prints `0.0` - so the extraction's own fallback is what stands in.
        surface_tension: DESIGNER_SURFACE_TENSION_N_PER_M,
        relative_volatility,
    })
}

/// `getTrayProperties`'s own relative-volatility spread.
///
/// The largest `y_i/x_i` over the smallest, over every component clearing `1e-10` on both sides,
/// **capped at `20.0`**; the class's default `2.0` where no component does.
fn relative_volatility(vapor_z: &[f64], liquid_z: &[f64]) -> f64 {
    let mut max_ratio: f64 = 0.0;
    let mut min_ratio: f64 = 1.0e10;
    for (y, x) in vapor_z.iter().zip(liquid_z) {
        if *x > 1.0e-10 && *y > 1.0e-10 {
            let ratio = y / x;
            if ratio > max_ratio {
                max_ratio = ratio;
            }
            if ratio < min_ratio {
                min_ratio = ratio;
            }
        }
    }
    let relative = if min_ratio > 0.0 {
        max_ratio / min_ratio
    } else {
        DEFAULT_RELATIVE_VOLATILITY
    };
    relative.min(20.0)
}

/// The one calculator state every tray in the summation loop is built from.
fn tray_state(
    geometry: &DesignerGeometry,
    diameter: Length,
    readings: &TrayReadings,
    relative_volatility: f64,
) -> TrayHydraulicsState {
    TrayHydraulicsState {
        tray_type: geometry.internals_type.clone(),
        column_diameter: diameter,
        tray_spacing: geometry.tray_spacing,
        weir_height: geometry.weir_height,
        // **The designer has no weir-length field at all**, so the calculator's own `-1.0` -
        // the rule that derives `0.73 D` - is what runs.
        weir_length: meters(-1.0),
        downcommer_area_fraction: geometry.downcommer_area_fraction,
        hole_diameter: geometry.hole_diameter,
        hole_area_fraction: geometry.hole_area_fraction,
        design_flood_fraction: geometry.design_flood_fraction,
        vapor_mass_flow: azoth_core::units::kilograms_per_second(readings.vapor_mass_flow),
        liquid_mass_flow: azoth_core::units::kilograms_per_second(readings.liquid_mass_flow),
        vapor_density: azoth_core::units::kilograms_per_cubic_meter(readings.vapor_density),
        liquid_density: azoth_core::units::kilograms_per_cubic_meter(readings.liquid_density),
        liquid_viscosity: azoth_core::units::pascal_seconds(readings.liquid_viscosity),
        surface_tension: azoth_core::units::newtons_per_meter(readings.surface_tension),
        relative_volatility,
    }
}

/// `ColumnInternalsDesigner.calculateTrayed`, over a solved column's trays.
///
/// # Errors
/// * [`AzothError::InvalidInput`] where the column carries no tray at all, which the class logs
///   and returns from - here it is a refusal, because a report of nothing is not an answer;
/// * whatever a tray's flash, its two streams or the sizing raise.
pub fn designer_report(
    trays: &[TrayProfile],
    components: &[String],
    geometry: &DesignerGeometry,
) -> Result<DesignerReport> {
    if trays.is_empty() {
        return Err(AzothError::invalid_input(
            "trays",
            "the column has no trays, so its internals have nothing to be designed against",
        ));
    }

    // ---- The controlling tray: the largest vapour mass flow, the first where two are equal.
    let mut controlling = 0;
    let mut largest = f64::NEG_INFINITY;
    for (index, tray) in trays.iter().enumerate() {
        let flow = readings(tray, components)?.vapor_mass_flow;
        if flow > largest {
            largest = flow;
            controlling = index;
        }
    }

    // ---- The diameter, which is the class's sizing branch or the stated override.
    let override_m = geometry.column_diameter_override.value;
    let required_diameter = if override_m > 0.0 {
        meters(override_m)
    } else {
        let controller = readings(&trays[controlling], components)?;
        // **Without the tray's relative volatility**, which is the class's own omission: this
        // calculator is built before the loop and never given one.
        let sized = size_column_diameter(&tray_state(
            geometry,
            meters(1.0),
            &controller,
            DEFAULT_RELATIVE_VOLATILITY,
        ))?;
        meters(sized)
    };

    // ---- The summation loop, at the diameter just resolved.
    let mut results = Vec::with_capacity(trays.len());
    let mut total_pressure_drop = 0.0;
    let mut max_percent_flood: f64 = 0.0;
    let mut min_percent_flood: f64 = 100.0;
    let mut design_ok = true;
    let mut efficiency_sum = 0.0;
    let mut walked = 0;
    for tray in trays {
        let tray_readings = readings(tray, components)?;
        let out = tray_hydraulics(tray_state(
            geometry,
            required_diameter,
            &tray_readings,
            tray_readings.relative_volatility,
        ))?;
        total_pressure_drop += out.total_tray_pressure_drop.value;
        if out.percent_flood > max_percent_flood {
            max_percent_flood = out.percent_flood;
        }
        if out.percent_flood > 0.0 && out.percent_flood < min_percent_flood {
            min_percent_flood = out.percent_flood;
        }
        if !out.design_ok {
            design_ok = false;
        }
        efficiency_sum += out.tray_efficiency;
        walked += 1;
        results.push(DesignerTrayResult {
            percent_flood: out.percent_flood,
            total_pressure_drop: out.total_tray_pressure_drop,
            tray_efficiency: out.tray_efficiency,
            design_ok: out.design_ok,
        });
    }
    let average_tray_efficiency = if walked > 0 {
        efficiency_sum / f64::from(walked)
    } else {
        EMPTY_AVERAGE_TRAY_EFFICIENCY
    };

    Ok(DesignerReport {
        required_diameter,
        controlling_tray_index: controlling,
        design_ok,
        max_percent_flood,
        min_percent_flood,
        average_tray_efficiency,
        total_pressure_drop: pascals(total_pressure_drop),
        total_pressure_drop_mbar: total_pressure_drop / 100.0,
        trays: results,
    })
}

#[cfg(test)]
mod report {
    use azoth_core::units::{kelvins, meters, pascals};

    use crate::models::distillation_column::distillation_column_outcome;

    use super::{DesignerGeometry, designer_report};

    /// The designer's own constructor defaults, which the capture's first row states.
    fn geometry() -> DesignerGeometry {
        DesignerGeometry {
            internals_type: "sieve".to_string(),
            tray_spacing: meters(0.6),
            weir_height: meters(0.05),
            hole_diameter: meters(12.7e-3),
            hole_area_fraction: 0.1,
            downcommer_area_fraction: 0.1,
            design_flood_fraction: 0.8,
            column_diameter_override: meters(-1.0),
        }
    }

    /// **`process_internals_designer.tsv`'s `binary_designer_defaults` row, end to end.**
    ///
    /// The column is solved from the row's own inputs and every quantity the designer publishes
    /// is held to what the class printed: the diameter it sized, the tray it sized from, the two
    /// flood extremes, the average efficiency, the summed pressure drop, and each tray's own
    /// flood, drop, efficiency and verdict.
    ///
    /// **The controlling tray is what makes this a measurement rather than an echo**: the class
    /// picks the largest vapour mass flow, and this row's is tray 2 - the feed stage, at `0.1088`
    /// kg/s against the reboiler's `0.0663`.
    #[test]
    #[allow(clippy::too_many_lines)] // one call, one assertion per captured quantity
    fn the_report_reproduces_the_captures_row() {
        let (outcome, _warnings) = distillation_column_outcome(
            &["methane".to_string(), "n-butane".to_string()],
            7.490_704_036_290_964,
            &[0.5, 0.5],
            pascals(2.0e6),
            kelvins(300.0),
            4,
            2,
            true,
            true,
            pascals(1.9e6),
            pascals(2.0e6),
            Some(kelvins(373.15)),
            Some(kelvins(253.15)),
            1.0e-6,
            200,
            None,
            None,
            Some("direct_substitution"),
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
        )
        .expect("the capture's binary column solves");

        let components = vec!["methane".to_string(), "n-butane".to_string()];
        let report = designer_report(&outcome.trays, &components, &geometry())
            .expect("the trayed report computes");

        // **`1e-5` and not a bit-comparison**, because every tray's state is *rebuilt* from the
        // profile this port carries rather than read off the class's own tray system: the two
        // agree to the last digits a flash can agree on, which is the same agreement the packed
        // column's report has. The discrete answers - the controlling index and the sized
        // diameter - are compared exactly, because the sizing lands on a table entry.
        let close = |got: f64, want: f64| (got - want).abs() / want.abs() < 1e-5;
        assert_eq!(report.controlling_tray_index, 2, "the controlling tray");
        assert!(
            !report.design_ok,
            "the class's own row is outside the window"
        );
        assert!(
            close(report.required_diameter.value, 0.5),
            "required_diameter: {}",
            report.required_diameter.value
        );
        assert!(
            close(report.max_percent_flood, 11.982_384_552_013_03),
            "max_percent_flood: {}",
            report.max_percent_flood
        );
        assert!(
            close(report.min_percent_flood, 4.438_168_765_587_891),
            "min_percent_flood: {}",
            report.min_percent_flood
        );
        assert!(
            close(report.average_tray_efficiency, 0.402_518_231_968_827_9),
            "average_tray_efficiency: {}",
            report.average_tray_efficiency
        );
        assert!(
            close(report.total_pressure_drop.value, 1_780.136_169_410_284_8),
            "total_pressure_drop: {}",
            report.total_pressure_drop.value
        );
        assert!(
            close(report.total_pressure_drop_mbar, 17.801_361_694_102_848),
            "total_pressure_drop_mbar: {}",
            report.total_pressure_drop_mbar
        );

        // **Every tray, in the column's own order** - the reboiler at index 0 to the condenser at
        // index 5, which is the six `getTrays()` walked.
        let floods = [
            11.982_384_552_013_03,
            11.691_032_790_829_645,
            11.526_376_182_416_941,
            6.819_195_914_181_229,
            6.477_112_100_578_219_5,
            4.438_168_765_587_891,
        ];
        let drops = [
            262.370_751_479_392_54,
            302.463_917_751_519_83,
            320.192_556_418_923_6,
            291.208_763_966_257_95,
            293.735_134_375_144_6,
            310.165_045_419_046_46,
        ];
        let efficiencies = [
            0.576_346_633_437_301_9,
            0.417_441_577_700_638_3,
            0.374_165_984_211_667_4,
            0.372_712_899_218_381_8,
            0.365_627_774_305_828,
            0.308_814_522_939_150_46,
        ];
        assert_eq!(report.trays.len(), 6, "the trays the designer walked");
        for (index, tray) in report.trays.iter().enumerate() {
            assert!(
                close(tray.percent_flood, floods[index]),
                "tray {index} percent_flood: {} against {}",
                tray.percent_flood,
                floods[index]
            );
            assert!(
                close(tray.total_pressure_drop.value, drops[index]),
                "tray {index} pressure drop: {} against {}",
                tray.total_pressure_drop.value,
                drops[index]
            );
            assert!(
                close(tray.tray_efficiency, efficiencies[index]),
                "tray {index} efficiency: {} against {}",
                tray.tray_efficiency,
                efficiencies[index]
            );
            assert!(!tray.design_ok, "tray {index} is outside the window");
        }

        // **And the stated override replaces the sizing branch**, which the capture's own
        // `binary_designer_diameter_override_0_5` row shows by reproducing the sized row to the
        // last digit - the sized value *is* `0.5`, so the two rows are one state twice.
        let mut overridden = geometry();
        overridden.column_diameter_override = meters(0.5);
        let other = designer_report(&outcome.trays, &components, &overridden)
            .expect("the overridden report computes");
        assert_eq!(other.required_diameter.value, 0.5);
        assert!(
            close(other.total_pressure_drop.value, 1_780.136_169_410_284_8),
            "an override at the sized value is the sized row: {}",
            other.total_pressure_drop.value
        );
    }
}
