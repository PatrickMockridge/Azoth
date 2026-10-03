//! `DistillationColumnMechanicalDesign.calcDesign`: the vessel around the internals tree.
//!
//! **The class runs two sizings and publishes the second one's diameter beside the first one's
//! numbers**, which is the shape `validation/neqsim/captures/process_column_mechanical_design.tsv`
//! was written to measure. The first pass is Souders-Brown on tray 0's two outlets; at that
//! diameter it also fixes the flooding factor and takes a weir loading from `0.7 D`. Then
//! `calculateContactorCapacity` builds the internals designer, hands it
//! [`resolve_rating_diameter`] as an override, and **replaces the diameter with the designer's
//! answer** - while the three quantities above keep the values they were computed at. So the
//! capture's `2.0` m row reports that diameter beside a flooding factor and a weir loading that
//! are the default row's to the last digit, and a wall thickness recomputed after the swap.
//!
//! **And the designer is never allowed to size.** `resolveRatingDiameter` answers the stated
//! override where there is one and the first pass's own diameter otherwise, and the first pass
//! is never below the class's `0.5` floor - so `getRequiredDiameter()` always returns the
//! override, the sizing branch is dead on this path, and the class's
//! `contactorDesignFloodFraction` (its own `0.70`) is never read. Both are stated here rather
//! than ported as knobs.
//!
//! **One quantity of `calcDesign` is computed and then thrown away**: `totalPressureDrop` is
//! assigned the tray estimate at the class's own line 376 and unconditionally overwritten by the
//! designer's sum at 512, so the port computes only the second.
//!
//! **The two tray-0 density reads are [`equinor/neqsim#4140`]**, and this is where the port
//! states the difference rather than reproducing it: the class takes `getFluid().getDensity()`
//! without ever asking the fluid to initialise, so on an absorber's *liquid* outlet it reads
//! `0.0` - and the weir loading dividing by it is `Infinity` while the tray pressure drop drops
//! its liquid-head term. azoth always has a density for a `(T, P, z)`, so this port publishes
//! the number the capture holds after one `initProperties()`: `645.1756172849634` where the
//! class reads `0.0`, `2.664675513235188` where it reads `Infinity`, `8.164586402782746` mbar
//! where it reads its own `5.0`. **The difference is the defect and not a correction of it** -
//! the capture carries both sides on every row, and the port is held to the repaired one.
//!
//! [`equinor/neqsim#4140`]: https://github.com/equinor/neqsim/issues/4140

use azoth_core::units::{Length, Power, meters};
use azoth_core::{AzothError, Result};
use azoth_hydraulics::tray_hydraulics::round_to_standard_diameter;

use crate::column::designer::{DesignerGeometry, designer_report};
use crate::kernels::distillation_column::{TrayProfile, tray_streams};
use crate::segment::phase::{liquid_view, vapour_view};
use crate::stream::Stream;

/// `trayEfficiency`'s own field initialiser, and the efficiency `calcDesign` divides the tray
/// count by.
pub const DEFAULT_TRAY_EFFICIENCY: f64 = 0.65;
/// `maxFloodingFactor`'s own initialiser - the fraction of flood the design velocity is.
pub const DEFAULT_MAX_FLOODING_FACTOR: f64 = 0.85;
/// `trayType`'s own initialiser. **It answers twice**, which the capture's valve row shows: the
/// Souders-Brown `kFactor` below and the designer's internals type.
pub const DEFAULT_TRAY_TYPE: &str = "sieve";
/// `contactorInternalsType`'s own initialiser.
pub const DEFAULT_CONTACTOR_INTERNALS_TYPE: &str = "auto";
/// `contactorDesignFloodFraction`'s own initialiser. **Carried and not declared**, because the
/// designer it reaches is always given an override and so never runs the branch that reads it.
pub const DEFAULT_CONTACTOR_DESIGN_FLOOD_FRACTION: f64 = 0.7;
/// `materialGrade`'s own initialiser.
pub const DEFAULT_MATERIAL_GRADE: &str = "SA-516-70";
/// `MechanicalDesign.maxOperationPressure`'s own field initialiser, bara.
///
/// **The equipment path overwrites it**: `DistillationColumn` sets it from
/// `getEconomicDesignPressure()` - the larger of the two endpoint pressures, `:6260` - before it
/// calls `calcDesign`, and the capture's rows construct the design directly and so never take
/// that path. This port declares it, so a case can state either.
pub const DEFAULT_MAX_OPERATION_PRESSURE_BARA: f64 = 100.0;
/// `MechanicalDesign.tensileStrength`'s own initialiser, MPa. Not declared: no path this id runs
/// reaches its setter.
pub const TENSILE_STRENGTH_MPA: f64 = 483.0;
/// `MechanicalDesign.jointEfficiency`'s own initialiser, fully radiographed.
pub const JOINT_EFFICIENCY: f64 = 1.0;
/// `MechanicalDesign.corrosionAllowance`'s own initialiser, mm.
pub const CORROSION_ALLOWANCE_MM: f64 = 0.0;
/// `getMaxOperationPressure() * 1.1`, the class's own design margin.
pub const DESIGN_PRESSURE_MARGIN: f64 = 1.1;
/// `getTensileStrength() * 0.4`, the allowable stress in `calculateRequiredWallThickness`.
pub const ALLOWABLE_STRESS_FRACTION: f64 = 0.4;
/// The floor `calculateRequiredWallThickness` takes its answer at, mm.
pub const MINIMUM_WALL_THICKNESS_MM: f64 = 6.0;
/// The Souders-Brown `kFactor` per tray type, m/s, **compared case-sensitively** as the class's
/// own `"valve".equals(trayType)` does.
pub const SIEVE_K_FACTOR: f64 = 0.1;
/// See [`SIEVE_K_FACTOR`].
pub const VALVE_K_FACTOR: f64 = 0.12;
/// See [`SIEVE_K_FACTOR`].
pub const BUBBLE_CAP_K_FACTOR: f64 = 0.08;
/// `dryTrayDp`, the class's own constant for a sieve tray, mbar/tray. It also stands in for the
/// whole tray pressure drop wherever #4140 zeroes the liquid density.
pub const DRY_TRAY_PRESSURE_DROP_MBAR: f64 = 5.0;
/// The weir length as a fraction of the diameter, which is `weirLength = columnDiameter * 0.7`.
pub const WEIR_LENGTH_FRACTION: f64 = 0.7;
/// The vapour disengagement section, m.
pub const TOP_SECTION_M: f64 = 1.0;
/// The liquid holdup section, m.
pub const BOTTOM_SECTION_M: f64 = 2.0;
/// One head, m - taken twice.
pub const HEAD_HEIGHT_M: f64 = 0.5;
/// The diameter the class falls back to where the vapour volume flow is not positive, before
/// rounding it through the table: `roundToStandardDiameter(0.5)`, which is `0.5`.
pub const DEGENERATE_DIAMETER_M: f64 = 0.5;
/// The class's own `500.0` in `calculateRequiredWallThickness`: `D/2` in millimetres.
pub const THICKNESS_DIAMETER_MM_PER_M: f64 = 500.0;

/// `DistillationColumnMechanicalDesign`'s own geometry and design statements.
///
/// **Six of these are the class's fields and the rest are the internals tree's**, which the class
/// passes straight through to the designer it builds - so a caller states each of them once. The
/// designer's own `designFloodFraction` is the exception: the class overwrites it with its own
/// [`DEFAULT_CONTACTOR_DESIGN_FLOOD_FRACTION`], which is why it is not here.
#[derive(Debug, Clone)]
pub struct MechanicalGeometry {
    /// `trayType`: the Souders-Brown `kFactor` and, through `resolveInternalsType`, the
    /// designer's internals type.
    pub tray_type: String,
    /// `contactorInternalsType`. `auto` resolves to [`Self::tray_type`] on this equipment.
    pub contactor_internals_type: String,
    /// `trayEfficiency`, which the tray count is divided by.
    pub tray_efficiency: f64,
    /// `maxFloodingFactor`, the fraction of flood the design velocity is.
    pub max_flooding_factor: f64,
    /// `materialGrade`, a vessel statement the class carries and reports.
    pub material_grade: String,
    /// `getMaxOperationPressure()`, bara - the wall thickness's design basis.
    pub max_operation_pressure_bara: f64,
    /// The tray spacing, m, passed to the designer.
    pub tray_spacing: Length,
    /// The weir height, m: the tray pressure drop's liquid head and the designer's own field.
    pub weir_height: Length,
    /// The sieve hole diameter, m.
    pub hole_diameter: Length,
    /// The hole area over the active area.
    pub hole_area_fraction: f64,
    /// The downcomer area over the total, which the sizing's total area is divided by.
    pub downcommer_area_fraction: f64,
    /// `columnDiameterOverride`, m: a stated diameter, or `-1` for the class's own sizing.
    pub column_diameter_override: Length,
}

/// `calcDesign`'s own report, minus the capacity result, the demister and the costs.
#[derive(Debug, Clone, PartialEq)]
pub struct MechanicalReport {
    /// `getColumnDiameter()`: the diameter the designer answered, which is
    /// [`resolve_rating_diameter`]'s own value - its first pass's rounded diameter where none was
    /// stated.
    pub vessel_diameter: Length,
    /// `getColumnHeight()`: `actualTrays * spacing + 1 + 2 + 2 * 0.5`, m. **The packed branch's
    /// `packedHeight + 4` is not reachable here**, because a packed contactor type is refused.
    pub vessel_height: Length,
    /// `getColumnWallThickness()`, mm - taken at the *final* diameter, after the designer
    /// replaced it.
    pub vessel_wall_thickness_mm: f64,
    /// `getActualTrays()`: `ceil(trays / efficiency)`, the class's own integer cast saturated at
    /// `i32::MAX`.
    pub actual_trays: usize,
    /// `getFloodingFactor()`: the actual velocity over the flooding velocity, **at the first
    /// pass's diameter** - the class's own order of statements, and not at the one it publishes.
    pub flooding_factor: f64,
    /// `getWeirLoading()`: the liquid volume flow over `0.7 D`, `m³/hr per m`, at the same first
    /// pass's diameter.
    pub weir_loading: f64,
    /// `getTrayPressureDrop()`, mbar/tray: the class's own `5.0` constant plus
    /// `weirHeight * liquidDensity * 9.81 / 100`.
    pub tray_pressure_drop_mbar: f64,
    /// `getTotalPressureDrop()`, **bar**: the designer's summed pressure drop over `1e5`, and
    /// not the tray estimate the class computes and discards.
    pub total_pressure_drop_bar: f64,
    /// `getReboilerDuty()`, kW.
    pub reboiler_duty_kw: f64,
    /// `getCondenserDuty()`, kW - **the absolute value**, which is the class's own.
    pub condenser_duty_kw: f64,
    /// `getMaterialGrade()`, carried and reported.
    pub material_grade: String,
}

/// `calcDesign`, over a solved column.
///
/// **The two duties are `Option`s because the class's guards are `getReboiler() != null`**, and a
/// column with no end of that kind reports zero rather than a duty.
///
/// # Errors
/// * [`AzothError::InvalidInput`] where the column carries no tray at all: `calcDesign` returns
///   early and leaves every quantity at its own initialiser, which is a vessel of zero diameter
///   rather than an answer;
/// * [`AzothError::InvalidInput`] where `tray_efficiency` is not positive and finite, because the
///   class's own `(int) Math.ceil(trays / efficiency)` has no guard on it;
/// * [`AzothError::InvalidInput`] for a `contactor_internals_type` this port does not carry;
/// * whatever a tray's flash, its two streams, their densities or the designer raise.
pub fn mechanical_design(
    trays: &[TrayProfile],
    components: &[String],
    distillate: &Stream,
    bottoms: &Stream,
    reboiler_duty: Option<Power>,
    condenser_duty: Option<Power>,
    geometry: &MechanicalGeometry,
) -> Result<MechanicalReport> {
    if trays.is_empty() {
        return Err(AzothError::invalid_input(
            "trays",
            "`DistillationColumnMechanicalDesign.calcDesign` returns early on a column with no \
             tray and leaves every quantity at its own initialiser - a vessel of zero diameter \
             and zero height, which is not a design of anything",
        ));
    }
    if !geometry.tray_efficiency.is_finite() || geometry.tray_efficiency <= 0.0 {
        return Err(AzothError::invalid_input(
            "tray_efficiency",
            format!(
                "a tray efficiency of {} is not positive and finite, and `calcDesign` divides \
                 the tray count by it with no guard - `(int) Math.ceil` answers `{}` there, which \
                 is a saturation of its own cast rather than a tray count",
                geometry.tray_efficiency,
                i32::MAX
            ),
        ));
    }
    let active_type = resolve_internals_type(geometry)?;

    // ---- The tray count, and what the efficiency makes of it.
    let actual_trays = ((trays.len() as f64) / geometry.tray_efficiency)
        .ceil()
        .clamp(0.0, f64::from(i32::MAX)) as usize;

    // ---- Tray 0's two outlets, which is the whole of the class's input surface. **The density
    // is the *phase's*, not the stream's system average**, and that is a measurement rather than
    // a preference: on this row the label rule splits the tray's vapour outlet into two phases
    // while its flash reports a vapour fraction of `1.0`, so the system density is `84.32487`
    // where the class reads `44.17833075904439` - the vapour's own.
    let (tray_zero_gas, tray_zero_liquid) = tray_streams(&trays[0], components)?;
    let (gas_stream, liquid_stream) = (
        fallback(&tray_zero_gas, distillate),
        fallback(&tray_zero_liquid, bottoms),
    );
    let vapor_view = vapour_view(gas_stream)?;
    let liquid_view = liquid_view(liquid_stream)?;
    let vapor_molar_flow = gas_stream.n * 3600.0;
    let liquid_molar_flow = liquid_stream.n * 3600.0;
    // The class's own `getMolarMass() * 1000`, which cancels against the `/ 1000` below and is
    // kept so this reads against the source.
    let vapor_mw = gas_stream.molar_mass()?.value * 1000.0;
    let liquid_mw = liquid_stream.molar_mass()?.value * 1000.0;
    let vapor_density = vapor_view.density;
    let liquid_density = liquid_view.density;
    let vapor_mass_flow = vapor_molar_flow * vapor_mw / 1000.0;
    let liquid_mass_flow = liquid_molar_flow * liquid_mw / 1000.0;
    let vapor_volume_flow = vapor_mass_flow / vapor_density;
    let liquid_volume_flow = liquid_mass_flow / liquid_density;

    // ---- The Souders-Brown pass, whose diameter the second sizing then throws away.
    let k_factor = k_factor(&geometry.tray_type);
    let u_flood_raw = k_factor * ((liquid_density - vapor_density) / vapor_density).sqrt();
    let u_flood = if u_flood_raw.is_nan() || u_flood_raw.is_infinite() || u_flood_raw <= 0.0 {
        k_factor
    } else {
        u_flood_raw
    };
    let u_design = u_flood * geometry.max_flooding_factor;
    let vapor_volume_flow_m3s = vapor_volume_flow / 3600.0;
    let first_pass = if vapor_volume_flow_m3s <= 0.0
        || vapor_volume_flow_m3s.is_nan()
        || vapor_volume_flow_m3s.is_infinite()
        || u_design <= 0.0
    {
        round_to_standard_diameter(DEGENERATE_DIAMETER_M)
    } else {
        let required_vapor_area = vapor_volume_flow_m3s / u_design;
        let total_area = required_vapor_area / (1.0 - geometry.downcommer_area_fraction);
        round_to_standard_diameter((4.0 * total_area / std::f64::consts::PI).sqrt())
    };

    // ---- The three quantities read *at* that diameter and published beside the one below.
    let actual_area = std::f64::consts::PI * (first_pass / 2.0).powi(2);
    let actual_vapor_area = actual_area * (1.0 - geometry.downcommer_area_fraction);
    let flooding_factor = vapor_volume_flow_m3s / actual_vapor_area / u_flood;
    let weir_length = first_pass * WEIR_LENGTH_FRACTION;
    let weir_loading = liquid_volume_flow / weir_length;

    let tray_pressure_drop_mbar =
        DRY_TRAY_PRESSURE_DROP_MBAR + geometry.weir_height.value * liquid_density * 9.81 / 100.0;

    // ---- The internals tree, given the rating diameter as an override so it never sizes.
    let rating_diameter = resolve_rating_diameter(geometry, first_pass);
    let internals = designer_report(
        trays,
        components,
        &DesignerGeometry {
            internals_type: active_type,
            tray_spacing: geometry.tray_spacing,
            weir_height: geometry.weir_height,
            hole_diameter: geometry.hole_diameter,
            hole_area_fraction: geometry.hole_area_fraction,
            downcommer_area_fraction: geometry.downcommer_area_fraction,
            design_flood_fraction: DEFAULT_CONTACTOR_DESIGN_FLOOD_FRACTION,
            column_diameter_override: meters(rating_diameter),
        },
    )?;
    // **The designer's answer replaces the diameter**, which the class guards with `> 0.0`.
    let vessel_diameter = if internals.required_diameter.value > 0.0 {
        internals.required_diameter.value
    } else {
        first_pass
    };

    let actual_trays_f = actual_trays as f64;
    Ok(MechanicalReport {
        vessel_diameter: meters(vessel_diameter),
        vessel_height: meters(
            actual_trays_f * geometry.tray_spacing.value
                + TOP_SECTION_M
                + BOTTOM_SECTION_M
                + 2.0 * HEAD_HEIGHT_M,
        ),
        // Taken at the final diameter, because the class recomputes it there.
        vessel_wall_thickness_mm: required_wall_thickness(
            vessel_diameter,
            geometry.max_operation_pressure_bara,
        ),
        actual_trays,
        flooding_factor,
        weir_loading,
        tray_pressure_drop_mbar,
        total_pressure_drop_bar: internals.total_pressure_drop.value / 1.0e5,
        reboiler_duty_kw: reboiler_duty.map_or(0.0, |duty| duty.value / 1000.0),
        condenser_duty_kw: condenser_duty.map_or(0.0, |duty| duty.value.abs() / 1000.0),
        material_grade: geometry.material_grade.clone(),
    })
}

/// The class's own fallback: where a tray-0 outlet carries no flow, the reading is taken from
/// the column's own product instead. Its guard is the *flow* and not the stream - the product is
/// read whenever the tray's traffic is not positive, even where that is also zero.
fn fallback<'a>(tray_stream: &'a Stream, product: &'a Stream) -> &'a Stream {
    if tray_stream.n <= 0.0 {
        product
    } else {
        tray_stream
    }
}

/// `resolveInternalsType`: the class's own `contactorInternalsType` unless it is `auto`, in
/// which case the equipment's own - which on a `DistillationColumn` is `trayType`.
///
/// # Errors
/// [`AzothError::InvalidInput`] for `packed`. The class carries that branch, and reaching it here
/// would need a packed preset, a bed height and a capacity factor and would put the height on the
/// class's `packedHeight + 4` path instead of the trayed one - a second machine, refused by name
/// rather than half-carried.
fn resolve_internals_type(geometry: &MechanicalGeometry) -> Result<String> {
    if geometry
        .contactor_internals_type
        .eq_ignore_ascii_case("auto")
    {
        return Ok(geometry.tray_type.clone());
    }
    if geometry
        .contactor_internals_type
        .eq_ignore_ascii_case("packed")
    {
        return Err(AzothError::invalid_input(
            "contactor_internals_type",
            "`packed` is refused: this port carries the trayed contactor, and the class's packed \
             branch reads a preset, a bed height and a hydraulic capacity factor this id declares \
             nowhere",
        ));
    }
    Ok(geometry.contactor_internals_type.clone())
}

/// `resolveRatingDiameter`: the stated override where there is one, and otherwise the diameter
/// the Souders-Brown pass resolved.
///
/// **The class's packed arm is not reachable here**, and its `column instanceof PackedColumn`
/// test could not be reached on this equipment even where a packed contactor type were carried -
/// so the two conditions that remain are these.
fn resolve_rating_diameter(geometry: &MechanicalGeometry, first_pass: f64) -> f64 {
    if geometry.column_diameter_override.value > 0.0 {
        geometry.column_diameter_override.value
    } else {
        first_pass
    }
}

/// `calculateRequiredWallThickness`, whose design pressure is `bara`, whose stress is `MPa` and
/// whose answer is millimetres - the class's own unit mixing, reproduced rather than corrected.
fn required_wall_thickness(diameter_m: f64, max_operation_pressure_bara: f64) -> f64 {
    let design_pressure = max_operation_pressure_bara * DESIGN_PRESSURE_MARGIN;
    let allowable_stress = TENSILE_STRENGTH_MPA * ALLOWABLE_STRESS_FRACTION;
    let denominator = allowable_stress * JOINT_EFFICIENCY - 0.6 * design_pressure;
    let calculated = if denominator > 0.0 {
        design_pressure * diameter_m * THICKNESS_DIAMETER_MM_PER_M / denominator
    } else {
        0.0
    };
    (calculated + CORROSION_ALLOWANCE_MM).max(MINIMUM_WALL_THICKNESS_MM)
}

/// The class's three Souders-Brown constants, compared with `String.equals` and so
/// case-sensitively: an unrecognised type is the sieve tray's `0.1`, which is what its own
/// constructor default is.
fn k_factor(tray_type: &str) -> f64 {
    if tray_type == "valve" {
        VALVE_K_FACTOR
    } else if tray_type == "bubble-cap" {
        BUBBLE_CAP_K_FACTOR
    } else {
        SIEVE_K_FACTOR
    }
}

#[cfg(test)]
mod report {
    use azoth_core::units::{kelvins, meters, pascals, watts};

    use crate::models::distillation_column::distillation_column_outcome;

    use super::{MechanicalGeometry, mechanical_design};

    /// The capture's binary column: `col1`, four trays, both ends, at the stated temperatures and
    /// pressures.
    fn binary() -> crate::kernels::distillation_column::ColumnOutcome {
        distillation_column_outcome(
            None,
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
        .expect("the capture's binary column solves")
        .0
    }

    /// The class's own constructor defaults, which the capture's first row states.
    fn geometry() -> MechanicalGeometry {
        MechanicalGeometry {
            tray_type: "sieve".to_string(),
            contactor_internals_type: "auto".to_string(),
            tray_efficiency: 0.65,
            max_flooding_factor: 0.85,
            material_grade: "SA-516-70".to_string(),
            max_operation_pressure_bara: 100.0,
            tray_spacing: meters(0.6),
            weir_height: meters(0.05),
            hole_diameter: meters(12.7e-3),
            hole_area_fraction: 0.1,
            downcommer_area_fraction: 0.1,
            column_diameter_override: meters(-1.0),
        }
    }

    /// **`process_column_mechanical_design.tsv`'s `binary_mechanical_defaults` row.**
    ///
    /// Every quantity the class publishes on this row, held to what it printed. **The two
    /// duties are the column's own over `1000`**, and the condenser's is taken absolute.
    #[test]
    fn the_report_reproduces_the_defaults_row() {
        let outcome = binary();
        let components = vec!["methane".to_string(), "n-butane".to_string()];
        let report = mechanical_design(
            &outcome.trays,
            &components,
            &outcome.distillate,
            &outcome.bottoms,
            Some(watts(outcome.reboiler_duty.value)),
            Some(watts(outcome.condenser_duty.value)),
            &geometry(),
        )
        .expect("the mechanical design computes");

        let close = |got: f64, want: f64, what: &str| {
            assert!(
                (got - want).abs() / want.abs() < 1e-5,
                "{what}: {got} against {want}"
            );
        };
        assert_eq!(report.actual_trays, 10, "the actual trays");
        assert_eq!(report.material_grade, "SA-516-70", "the material grade");
        close(report.vessel_diameter.value, 0.5, "vessel_diameter");
        close(report.vessel_height.value, 10.0, "vessel_height");
        close(
            report.vessel_wall_thickness_mm,
            216.194_968_553_459_12,
            "vessel_wall_thickness_mm",
        );
        close(
            report.flooding_factor,
            0.029_981_966_003_039_984,
            "flooding_factor",
        );
        close(report.weir_loading, 4.992_709_047_038_033, "weir_loading");
        close(
            report.tray_pressure_drop_mbar,
            7.137_331_084_047_494,
            "tray_pressure_drop_mbar",
        );
        close(
            report.total_pressure_drop_bar,
            0.017_801_361_694_102_846,
            "total_pressure_drop_bar",
        );
        // **The two duties are pass-throughs, so their comparison is the *column* port's and not
        // this module's.** The class reads `getReboiler()`/`getCondenser()`; azoth hands over
        // `distillation_column_outcome`'s own duties, which carry the two libraries' enthalpy
        // offset - `tests/column.rs` measures it at `0.39` W of `21323` and `1.06` W of `47786`
        // and holds the same numbers to `5` W. Held to that measured offset rather than to this
        // report's `1e-5`, which would demand `0.48` W and which no duty of this port meets.
        let duty_close = |got: f64, want: f64, what: &str| {
            assert!((got - want).abs() < 5.0e-3, "{what}: {got} against {want}");
        };
        duty_close(
            report.reboiler_duty_kw,
            47.786_583_893_810_15,
            "reboiler_duty_kw",
        );
        duty_close(
            report.condenser_duty_kw,
            21.323_042_789_303_567,
            "condenser_duty_kw",
        );
    }

    /// **The capture's `binary_mechanical_override` row, which is the whole of the finding.**
    ///
    /// A stated `2.0` m diameter replaces the rating diameter the designer is driven at, and it
    /// moves **one** published quantity - the wall thickness, because it is the one recomputed
    /// after the swap. The flooding factor and the weir loading are the defaults row's to the
    /// last digit, because the class computed them one statement earlier at the `0.5` m it then
    /// threw away.
    #[test]
    fn a_stated_diameter_moves_the_thickness_and_not_the_flooding() {
        let outcome = binary();
        let components = vec!["methane".to_string(), "n-butane".to_string()];
        let mut stated = geometry();
        stated.column_diameter_override = meters(2.0);
        let report = mechanical_design(
            &outcome.trays,
            &components,
            &outcome.distillate,
            &outcome.bottoms,
            Some(watts(outcome.reboiler_duty.value)),
            Some(watts(outcome.condenser_duty.value)),
            &stated,
        )
        .expect("the stated-diameter design computes");

        let close = |got: f64, want: f64, what: &str| {
            assert!(
                (got - want).abs() / want.abs() < 1e-5,
                "{what}: {got} against {want}"
            );
        };
        close(report.vessel_diameter.value, 2.0, "vessel_diameter");
        close(
            report.vessel_wall_thickness_mm,
            864.779_874_213_836_5,
            "vessel_wall_thickness_mm",
        );
        close(
            report.total_pressure_drop_bar,
            0.016_680_181_876_092_014,
            "total_pressure_drop_bar",
        );
        // **The two that did not move**, which is the row's reason for existing: they are the
        // defaults row's own numbers, at the diameter the override replaced.
        close(
            report.flooding_factor,
            0.029_981_966_003_039_984,
            "flooding_factor",
        );
        close(report.weir_loading, 4.992_709_047_038_033, "weir_loading");
        close(report.vessel_height.value, 10.0, "vessel_height");
    }

    /// **The capture's `binary_mechanical_valve` row**: one tray type, two answers.
    #[test]
    fn the_valve_row_moves_the_souders_brown_factor_and_the_internals_type() {
        let outcome = binary();
        let components = vec!["methane".to_string(), "n-butane".to_string()];
        let mut valved = geometry();
        valved.tray_type = "valve".to_string();
        let report = mechanical_design(
            &outcome.trays,
            &components,
            &outcome.distillate,
            &outcome.bottoms,
            Some(watts(outcome.reboiler_duty.value)),
            Some(watts(outcome.condenser_duty.value)),
            &valved,
        )
        .expect("the valve design computes");

        let close = |got: f64, want: f64, what: &str| {
            assert!(
                (got - want).abs() / want.abs() < 1e-5,
                "{what}: {got} against {want}"
            );
        };
        close(
            report.flooding_factor,
            0.024_984_971_669_199_99,
            "flooding_factor",
        );
        close(
            report.total_pressure_drop_bar,
            0.017_755_928_607_031_072,
            "total_pressure_drop_bar",
        );
        close(report.vessel_diameter.value, 0.5, "vessel_diameter");
        // The weir loading is untouched by the tray type, as it is untouched by the override: it
        // reads the liquid flow, the diameter and nothing else.
        close(report.weir_loading, 4.992_709_047_038_033, "weir_loading");
    }

    /// **The packed contactor is refused by name rather than half-carried.**
    #[test]
    fn a_packed_contactor_is_refused() {
        let outcome = binary();
        let components = vec!["methane".to_string(), "n-butane".to_string()];
        let mut packed = geometry();
        packed.contactor_internals_type = "packed".to_string();
        let error = mechanical_design(
            &outcome.trays,
            &components,
            &outcome.distillate,
            &outcome.bottoms,
            None,
            None,
            &packed,
        )
        .expect_err("the packed branch is refused");
        assert!(
            error.to_string().contains("contactor_internals_type"),
            "{error}"
        );
    }

    /// **A tray efficiency the class divides by without a guard** is refused here instead, and a
    /// column with no tray is refused rather than answered with zeros.
    #[test]
    fn the_two_partial_inputs_are_refused() {
        let outcome = binary();
        let components = vec!["methane".to_string(), "n-butane".to_string()];
        let mut zero = geometry();
        zero.tray_efficiency = 0.0;
        assert!(
            mechanical_design(
                &outcome.trays,
                &components,
                &outcome.distillate,
                &outcome.bottoms,
                None,
                None,
                &zero,
            )
            .is_err(),
            "a zero efficiency is refused"
        );
        assert!(
            mechanical_design(
                &[],
                &components,
                &outcome.distillate,
                &outcome.bottoms,
                None,
                None,
                &geometry(),
            )
            .is_err(),
            "a column with no tray is refused"
        );
    }
}
