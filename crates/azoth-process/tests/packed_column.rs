//! The packed column, against `validation/neqsim/captures/process_packed_column.tsv`.

use azoth_core::units::{kelvins, pascals};
use azoth_process::kernels::packed_column::stage_count;
use azoth_process::{distillation_column, packed_column};

fn relative(actual: f64, expected: f64, tolerance: f64, what: &str) {
    let scale = 1.0 + actual.abs().max(expected.abs());
    assert!(
        (actual - expected).abs() <= tolerance * scale,
        "{what}: {actual} vs {expected}, {} relative",
        (actual - expected).abs() / scale
    );
}

fn absolute(actual: f64, expected: f64, tolerance: f64, what: &str) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{what}: {actual} vs {expected}, {} absolute",
        (actual - expected).abs()
    );
}

/// **The stage rule, against the capture's eleven stage-count rows.** Each row builds
/// `PackedColumn(name, height, packing, true, true)` and prints the trays it made, with no solve
/// at all - so this is the class's own arithmetic, measured.
///
/// The last two of the three non-positive heights are the interesting ones: `0.0` m and `-1.0` m
/// both give two middle trays rather than an error, because neither `estimateStages` nor
/// `setPackedHeight` has a domain check.
#[test]
fn the_stage_count_is_the_heights_own_arithmetic() {
    let rows = [
        (-1.0, 2usize),
        (0.0, 2),
        (0.2, 2),
        (0.4, 2),
        (0.75, 2),
        (1.0, 2),
        (1.2, 3),
        (2.5, 5),
        (5.0, 10),
        (5.3, 11),
        (10.0, 20),
    ];
    for (height, stages) in rows {
        assert_eq!(
            stage_count(height),
            stages,
            "{height} m is {stages} middle trays in the capture"
        );
    }
}

/// **The id's claim, checked against the base rather than against itself.** `PackedColumn extends
/// DistillationColumn` and its `run` is `super.run(id)`, so a packed column at the height whose
/// constructor makes four middle trays is the base column at four stages - and NeqSim reports the
/// two bit-identical, on all twenty-two captured quantities.
#[test]
fn the_packed_column_is_the_base_column_at_the_heights_stage_count() {
    let packed = packed_column(
        &["methane".to_string(), "n-butane".to_string()],
        7.490704036290964,
        &[0.5, 0.5],
        pascals(2.0e6),
        kelvins(300.0),
        2.0,
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
        None,
        None,
        None,
        None,
        None,
        None,
        None,
    )
    .expect("the packed column converges");

    let base = distillation_column(
        &["methane".to_string(), "n-butane".to_string()],
        7.490704036290964,
        &[0.5, 0.5],
        pascals(2.0e6),
        kelvins(300.0),
        stage_count(2.0),
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
        // No per-stage efficiency overrides: this state states neither field.
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
    .expect("the base column converges");

    assert_eq!(
        packed.tray_temperature, base.tray_temperature,
        "the packing does not move the profile"
    );
    assert_eq!(packed.distillate_n, base.distillate_n);
    assert_eq!(packed.bottoms_n, base.bottoms_n);
    assert_eq!(packed.condenser_duty.value, base.condenser_duty.value);
}

/// **The same row against NeqSim**, which is what the case states: the capture's
/// `packed_distillation_binary_2m` is the base column's own `binary_rigorous` state at a packed
/// height of `2.0` m.
#[test]
fn the_packed_binary_row_reaches_neqsims_profile() {
    let out = packed_column(
        &["methane".to_string(), "n-butane".to_string()],
        7.490704036290964,
        &[0.5, 0.5],
        pascals(2.0e6),
        kelvins(300.0),
        2.0,
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
        None,
        None,
        None,
        None,
        None,
        None,
        None,
    )
    .expect("the packed column converges");

    assert_eq!(
        out.tray_temperature.len(),
        6,
        "four middle trays and two ends"
    );
    // NeqSim takes 15 iterations at this gate and the port 14: where a solve stops is its own.
    assert!(out.iterations <= 25, "in {} iterations", out.iterations);
    let temperatures = [
        373.15,
        336.15382381733946,
        303.0710539116486,
        302.1225394582587,
        296.84208939635556,
        253.15,
    ];
    for (i, (tray, expected)) in out.tray_temperature.iter().zip(temperatures).enumerate() {
        absolute(tray.value, expected, 2.0e-3, &format!("tray {i}"));
    }
    relative(
        out.distillate_n,
        3.7914997294607584,
        1.0e-4,
        "the distillate",
    );
    relative(out.bottoms_n, 3.6992043068302056, 1.0e-4, "the bottoms");
    relative(
        out.distillate_z[0],
        0.9659099052324753,
        1.0e-5,
        "the distillate's methane",
    );
    // **The duties carry the libraries' ideal-gas offset**, as the base column's own case says:
    // measured at `2.2e-5` relative here, which an absolute gate on a wattage would either miss
    // or fail on. The base case's own tolerance is `2e-4` relative.
    relative(
        out.condenser_duty.value,
        -21323.042789303567,
        2.0e-4,
        "the condenser's duty",
    );
    relative(
        out.reboiler_duty.value,
        47786.58389381015,
        2.0e-4,
        "the reboiler's duty",
    );
}

/// **The rule bites on the same feed.** The capture's `packed_distillation_binary_2m3` is
/// `ceil(2.3 / 0.5)` = five middle trays - seven trays with the ends - so a port that solved the
/// base's four stages and ignored `packed_height` would report six.
#[test]
fn the_height_moves_the_stage_count_and_the_profile() {
    let out = packed_column(
        &["methane".to_string(), "n-butane".to_string()],
        7.490704036290964,
        &[0.5, 0.5],
        pascals(2.0e6),
        kelvins(300.0),
        2.3,
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
        None,
        None,
        None,
        None,
        None,
        None,
        None,
    )
    .expect("the packed column converges");

    assert_eq!(
        out.tray_temperature.len(),
        7,
        "five middle trays and two ends"
    );
    absolute(
        out.tray_temperature[4].value,
        301.985281139596,
        2.0e-3,
        "tray 4",
    );
    relative(
        out.distillate_n,
        3.7914996721606786,
        1.0e-4,
        "the distillate",
    );
    relative(
        out.condenser_duty.value,
        -21283.206254634413,
        2.0e-4,
        "the condenser's duty",
    );
}

/// **The one packing parameter the class refuses rather than reads.**
/// `setPackingHydraulicCapacityFactor` throws for a value that is not positive and finite, so the
/// port refuses it too - the other four are carried inertly.
#[test]
fn a_non_positive_capacity_factor_is_refused() {
    for factor in [0.0, -1.0, f64::NAN] {
        let error = packed_column(
            &["methane".to_string(), "n-butane".to_string()],
            7.490704036290964,
            &[0.5, 0.5],
            pascals(2.0e6),
            kelvins(300.0),
            2.0,
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
            None,
            Some(factor),
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
            None,
            None,
            None,
            None,
            None,
        )
        .expect_err("the class's own setter refuses this");
        assert!(
            format!("{error}").contains("packing_hydraulic_capacity_factor"),
            "{error}"
        );
    }
}

/// **The injected packing parameters are declarations, not settings**: a packing type, a
/// flood fraction and a diameter the class would read only in its hydraulics report leave the
/// answer bit-identical.
#[test]
fn the_packing_parameters_do_not_move_the_solve() {
    let plain = packed_column(
        &["methane".to_string(), "n-butane".to_string()],
        7.490704036290964,
        &[0.5, 0.5],
        pascals(2.0e6),
        kelvins(300.0),
        2.0,
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
        None,
        None,
        None,
        None,
        None,
        None,
        None,
    )
    .expect("the packed column converges");

    let dressed = packed_column(
        &["methane".to_string(), "n-butane".to_string()],
        7.490704036290964,
        &[0.5, 0.5],
        pascals(2.0e6),
        kelvins(300.0),
        2.0,
        2,
        true,
        true,
        pascals(1.9e6),
        pascals(2.0e6),
        Some(kelvins(373.15)),
        Some(kelvins(253.15)),
        1.0e-6,
        200,
        Some("Mellapak-250Y"),
        Some(true),
        Some(0.75),
        Some(1.30),
        Some(1.5),
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
        None,
        None,
        None,
        None,
    )
    .expect("the packed column converges");

    assert_eq!(
        plain.tray_temperature, dressed.tray_temperature,
        "the packing group is read by the report after the solve and by nothing in it"
    );
    assert_eq!(plain.distillate_n, dressed.distillate_n);
    assert_eq!(plain.condenser_duty.value, dressed.condenser_duty.value);
}

/// **The reactive section on the packed column, which is the base's own and is the identity on
/// this fluid.**
///
/// `PackedColumn` inherits `setReactive` and does not override it, so the flag reaches the four
/// middle trays the height derives - `reactive_tray_flags=011110`, the two ends never - and
/// NeqSim's own row `packed_distillation_binary_2m_reactive` is **bit-identical to
/// `packed_distillation_binary_2m` across all sixty-nine captured keys**. That is not a
/// coincidence and not a no-op in the port: this fluid has no independent reaction. Methane and
/// n-butane give a carbon-and-hydrogen formula matrix of rank two against two species, so
/// `independent_reactions()` is zero and the reactive route *is* the equilibrium flash - the
/// property the class's own `ReactiveDistillationTest` asserts at `0.01` kg/hr.
///
/// **So the case is the identity**, and what it holds is the *routing*: a
/// `reactive = true` that silently never reached the trays would pass this test and fail the
/// capture's flag line, which is why the row is a case at all.
#[test]
fn the_reactive_section_is_the_delegations_own_on_a_fluid_that_does_not_react() {
    let reactive = |reactive: bool| {
        packed_column(
            &["methane".to_string(), "n-butane".to_string()],
            7.490704036290964,
            &[0.5, 0.5],
            pascals(2.0e6),
            kelvins(300.0),
            2.0,
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
            None,
            None,
            None,
            None,
            Some(reactive),
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
            None,
            None,
        )
        .expect("the packed column converges")
    };
    let plain = reactive(false);
    let dressed = reactive(true);

    // **The same state, and not bit-identical in this port.** NeqSim's two rows *are*
    // bit-identical, because its reactive route delegates to the plain flash when there are no
    // independent reactions; this port runs its own reactive solve, which reaches the same point
    // from a different arithmetic - measured at `3.7e-12` K on the worst tray, a relative
    // `1.1e-14`. Holding the port to NeqSim's bit-identity would be asserting something about
    // the *route* and not about the answer.
    for (plain_t, dressed_t) in plain.tray_temperature.iter().zip(&dressed.tray_temperature) {
        relative(
            dressed_t.value,
            plain_t.value,
            1.0e-11,
            "a tray on a fluid with no independent reaction",
        );
    }
    for (plain_n, dressed_n) in plain.tray_gas_n.iter().zip(&dressed.tray_gas_n) {
        relative(*dressed_n, *plain_n, 1.0e-11, "a tray's vapour");
    }
    relative(
        dressed.distillate_n,
        plain.distillate_n,
        1.0e-12,
        "distillate",
    );
    for (plain_z, dressed_z) in plain.distillate_z.iter().zip(&dressed.distillate_z) {
        relative(
            *dressed_z,
            *plain_z,
            1.0e-12,
            "the distillate's composition",
        );
    }
    relative(
        dressed.condenser_duty.value,
        plain.condenser_duty.value,
        1.0e-10,
        "condenser duty",
    );
    assert_eq!(
        plain.iterations, dressed.iterations,
        "and the two solves take the same number of passes"
    );
}

/// **The three parameters are the class's two forms**, so half a section is refused rather than
/// guessed - the same resolution `process.distillation_column` states, shared by all four ids.
#[test]
fn half_a_reactive_section_is_refused() {
    let error = packed_column(
        &["methane".to_string(), "n-butane".to_string()],
        7.490704036290964,
        &[0.5, 0.5],
        pascals(2.0e6),
        kelvins(300.0),
        2.0,
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
        None,
        None,
        None,
        None,
        Some(true),
        Some(1),
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
        None,
    )
    .expect_err("one bound alone is not a section the class can state");
    assert!(
        error.to_string().contains("both bounds or by neither"),
        "{error}"
    );
}
