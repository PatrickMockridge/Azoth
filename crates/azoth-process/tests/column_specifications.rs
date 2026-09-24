//! The column's product specifications, against `validation/neqsim/captures/process_column.tsv`.
//!
//! `ColumnSpecification` is two degrees of freedom at two ends, and the class splits its five
//! types in two: a **reflux ratio** and a **duty** are written straight onto the end, and the
//! other three are driven by an **outer secant on the end's temperature**. Each row here is one
//! type, on the same binary column the profile tests use, so the specification's effect is
//! measurable against a column whose unsolved answer is already known.
//!
//! **Two of the six captured rows are evidence rather than oracles.** A recovery specification
//! does not converge in NeqSim (its secant drives the column to the whole feed and stops), and a
//! reflux-ratio specification converges *there* on a state this port's `pv_reflux_flash` cannot
//! find. Both are recorded below rather than excused.

use azoth_core::units::{kelvins, pascals};
use azoth_process::Stream;
use azoth_process::kernels::distillation_column::{
    ColumnSetup, SolverType, Specification, SpecificationKind, distillation_column,
};

fn column(top: Option<Specification>, bottom: Option<Specification>, pin: bool) -> ColumnSetup {
    ColumnSetup {
        feed: Stream::from_pt(
            vec!["methane".into(), "n-butane".into()],
            vec![0.5, 0.5],
            7.490704036290964,
            pascals(20.0e5),
            kelvins(300.0),
        )
        .expect("the fluid resolves"),
        feed_stage: 2,
        number_of_stages: 4,
        has_reboiler: true,
        has_condenser: true,
        top_pressure: pascals(19.0e5),
        bottom_pressure: pascals(20.0e5),
        condenser_temperature: if pin { Some(kelvins(253.15)) } else { None },
        reboiler_temperature: Some(kelvins(373.15)),
        temperature_tolerance: 1.0e-6,
        max_iterations: 200,
        top_specification: top,
        bottom_specification: bottom,
        solver_type: SolverType::DirectSubstitution,
    }
}

fn spec(kind: SpecificationKind, target: f64, component: Option<&str>) -> Option<Specification> {
    Some(Specification {
        kind,
        target,
        component: component.map(String::from),
    })
}

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

/// **A top purity specification, driven by the outer secant on the condenser's temperature.**
/// NeqSim searches its way to `241.0` K, where the distillate is `97.9996 %` methane, and the
/// port finds the same temperature to **`1e-9` K** - which is the strongest agreement in this
/// tranche, because the searched quantity is the one thing the specification fixes.
#[test]
fn a_top_purity_specification_searches_the_condenser_temperature() {
    let out = distillation_column(&column(
        spec(SpecificationKind::ProductPurity, 0.98, Some("methane")),
        None,
        true,
    ))
    .expect("the purity specification converges");

    absolute(
        out.trays[5].temperature.value,
        240.99878449426998,
        1.0e-6,
        "the searched condenser temperature",
    );
    relative(out.distillate.n, 3.735722909260181, 1.0e-5, "distillate");
    // The specification is *met*, which is the only thing that makes it a specification: the
    // distillate is 0.98 methane to within the class's own tolerance of 1e-4.
    absolute(
        out.distillate.z[0],
        0.98,
        1.0e-4,
        "the purity it was asked for",
    );
    relative(out.bottoms.n, 3.7549811270307836, 1.0e-5, "bottoms");
    absolute(
        out.condenser_duty.value,
        -23926.07136724649,
        5.0,
        "condenser duty",
    );
}

/// **A top flow-rate specification**, whose target is in **mol/hr** - `defaultTargetUnit`'s own
/// unit for the type. 14000 mol/hr is 3.888888 mol/s, and the port lands on it to nine figures,
/// which is what makes the unit a measurement rather than a convention: a target read in mol/s
/// would be a column asked for a millionth of its flow.
#[test]
fn a_top_flow_rate_specification_is_in_mol_per_hour() {
    let out = distillation_column(&column(
        spec(SpecificationKind::ProductFlowRate, 14000.0, None),
        None,
        true,
    ))
    .expect("the flow-rate specification converges");

    relative(out.distillate.n, 3.888888888888851, 1.0e-6, "distillate");
    absolute(
        out.distillate.n * 3600.0,
        14000.0,
        1.0e-3,
        "the flow it was asked for",
    );
    absolute(
        out.condenser_duty.value,
        -17715.90123296374,
        5.0,
        "condenser duty",
    );
}

/// **A duty specification, applied directly rather than searched for.** Its condenser is not
/// pinned, so its temperature is the flash's own answer and the duty lands on the target to
/// `3e-8` relative.
#[test]
fn a_duty_specification_is_applied_to_the_end() {
    let out = distillation_column(&column(
        spec(SpecificationKind::Duty, -20000.0, None),
        None,
        false,
    ))
    .expect("the duty specification converges");

    relative(
        out.condenser_duty.value,
        -20000.000040496427,
        1.0e-6,
        "the duty it was asked for",
    );
    relative(out.distillate.n, 3.8249569643784094, 1.0e-5, "distillate");
}

/// **A duty specification under a temperature pin is inert, and the capture measures it.**
/// The condenser's own `outTemperature` is stated, so its flash is a `TPflash` and
/// `setHeatInput` is never read: the row reports `-21323.04` W for a specification of `-20000`.
///
/// **Nothing in the class reports this.** A `DUTY` is not *adjusted*, so
/// `needsAdjustment` is false, no error is evaluated and `getLastTopSpecificationResidual`
/// returns `0.0` - a zero that means "not checked" rather than "satisfied". The port reproduces
/// the silence: its duty is the state's, and the specification's target is unreachable.
#[test]
fn a_duty_specification_under_a_pin_moves_nothing() {
    let out = distillation_column(&column(
        spec(SpecificationKind::Duty, -20000.0, None),
        None,
        true,
    ))
    .expect("the column solves");

    absolute(
        out.condenser_duty.value,
        -21323.042789303567,
        5.0,
        "the pinned duty",
    );
    relative(
        out.distillate.n,
        3.7914997294607584,
        1.0e-5,
        "the unsold column's distillate",
    );
    // **The pin wins, and that is NeqSim's own order of tests.** `SimpleTray.run` checks its
    // stated outlet temperature before anything else and takes a `TPflash` there, so the heat
    // input a `DUTY` specification wrote is never read. Reproducing the silence is the port
    // rule; honouring the duty here would be an improvement NeqSim does not make.
    // And the specification's own target is *not* met, which is the finding: a value nine per
    // cent away reads as satisfied because nothing evaluates it.
    assert!(
        (out.condenser_duty.value - (-20000.0)).abs() > 1000.0,
        "the duty is {} W against a target of -20000",
        out.condenser_duty.value
    );
}

/// **A reflux-ratio specification is not reproduced, and the reason is measured.**
///
/// The class applies it directly to the condenser, whose own mode is `PVrefluxflash(ratio, 0)`,
/// a temperature search for the state whose vapour fraction is `1/(1 + ratio)`. On this row
/// NeqSim reports `402.24` K, which for a 19 bara mixture that is 82 % methane is past its dew
/// point: a state there is single-phase, and this port's `pv_reflux_flash` therefore refuses to
/// converge rather than landing on it.
///
/// So the row is *not* an oracle and is declared uncased. What is asserted here is the refusal,
/// because a spec type that quietly returned the wrong state would be worse than one that says
/// it cannot.
#[test]
fn a_reflux_ratio_specification_is_refused_rather_than_guessed() {
    let out = distillation_column(&column(
        spec(SpecificationKind::RefluxRatio, 1.5, None),
        None,
        false,
    ));

    assert!(
        out.is_err(),
        "the reflux-ratio row must refuse rather than land on NeqSim's state: {:?}",
        out.map(|o| (o.distillate.n, o.trays[5].temperature.value))
    );
}

/// **The same purity specification at the other location**, which is `setBottomSpecification`.
///
/// It drives the **reboiler's** temperature where the top's drives the condenser's - `374.80` K
/// against the `373.15` it started from - and the product it constrains is the bottoms. So this
/// is the row that says a location is a degree of freedom rather than a mirrored spelling of the
/// top's, and the port lands on the searched temperature to `6e-10` K.
#[test]
fn a_bottom_purity_specification_searches_the_reboiler_temperature() {
    let out = distillation_column(&column(
        None,
        spec(SpecificationKind::ProductPurity, 0.98, Some("n-butane")),
        true,
    ))
    .expect("the bottom purity specification converges");

    absolute(
        out.trays[0].temperature.value,
        374.80351966552143,
        1.0e-6,
        "the searched reboiler temperature",
    );
    absolute(
        out.bottoms.z[1],
        0.98,
        1.0e-4,
        "the purity it was asked for",
    );
    relative(out.bottoms.n, 3.6898141185673485, 1.0e-6, "bottoms");
    relative(out.distillate.n, 3.800889917723617, 1.0e-6, "distillate");
    absolute(
        out.reboiler_duty.value,
        49560.82240422304,
        5.0,
        "reboiler duty",
    );
}

/// A specification that constrains a component must name it.
#[test]
fn a_purity_specification_needs_its_component() {
    assert!(
        distillation_column(&column(
            spec(SpecificationKind::ProductPurity, 0.98, None),
            None,
            true
        ))
        .is_err()
    );
}

/// **A recovery specification is the row the capture declares uncased**, and the port's answer
/// to it is a refusal rather than a state.
///
/// NeqSim's secant walks the condenser towards the whole feed - `D = 7.25` mol/s against a
/// `7.49` feed - and stops at a residual of `2.8e-2` against its own `1e-4` tolerance. The same
/// walk here reaches a tray composition that sums to `1.026`, which the library refuses rather
/// than renormalises, so the row ends in an error.
///
/// **The Python twin refuses it through a different door**: its walk reaches
/// `_flash_newton.split_at`'s `1 - beta` division on the vapour side, where the Rust's
/// feasibility gate keeps the split out of that state, and raises `ZeroDivisionError` instead.
/// Both sides decline the row; the divergence is `azoth-eos`'s and is recorded here rather than
/// patched from this crate.
#[test]
fn a_recovery_specification_is_the_row_neqsim_cannot_drive() {
    let out = distillation_column(&column(
        spec(SpecificationKind::ComponentRecovery, 0.97, Some("methane")),
        None,
        false,
    ));
    assert!(
        out.is_err(),
        "the recovery row must refuse rather than return a state: {:?}",
        out.map(|o| (o.distillate.n, o.bottoms.n))
    );
}
