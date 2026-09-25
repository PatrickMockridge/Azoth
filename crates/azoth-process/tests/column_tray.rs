//! The column's stage, against `validation/neqsim/captures/process_column_tray.tsv`.
//!
//! Every row below drives one `SimpleTray` in the probe and one [`tray`] here, on the same
//! four-component feed. The stage is pinned before a column exists to contain it, which is
//! the point of D2: a column's conversation is between stages, and a stage that is wrong
//! makes that conversation meaningless.

use azoth_core::units::{kelvins, pascals, watts};
use azoth_process::Stream;
use azoth_process::column::{SideDraws, TrayOutcome, tray};

/// The capture's fluid: a four-component mixture whose split exercises both ends of the
/// relative-volatility range.
fn feed(t: f64, p_bar: f64, n: f64) -> Stream {
    Stream::from_pt(
        vec![
            "methane".into(),
            "ethane".into(),
            "propane".into(),
            "n-butane".into(),
        ],
        vec![0.1, 0.3, 0.4, 0.2],
        n,
        pascals(p_bar * 1.0e5),
        kelvins(t),
    )
    .expect("the fluid resolves")
}

/// The molar enthalpy the stage's flash was asked for: the inlets' own records plus the
/// duty, over the flow.
fn flash_target(inlets: &[Stream], duty: f64) -> f64 {
    let n: f64 = inlets.iter().map(|s| s.n).sum();
    (inlets.iter().map(|s| s.n * s.h.value).sum::<f64>() + duty) / n
}

/// **The stage is checked against its own enthalpy balance, because NeqSim's `PHflash` does
/// not keep one.**
///
/// `ProcessProbe tray`'s last rows ask NeqSim's `PHflash(h, 0)` for a molar enthalpy on a
/// **fresh fluid**, with no tray involved, and read the state's own enthalpy back: the
/// answers are `+2.5e-9`, `+179.48`, `-27.29` and `+0.00002` J/mol out for requests of `0`,
/// `+2500`, `+5000` and `-2000` J/mol. A 179 J/mol error on a 2500 J/mol step is a solver's
/// stopping rule and not its equations of state.
///
/// So the capture's temperature is not an oracle - it is the state NeqSim's solve landed on,
/// not the state its enthalpy implies - and the port is self-consistent instead. Everything
/// the two libraries *do* agree on is oracled: the outlet flows, compositions and enthalpies,
/// to `0.06` J/mol. This checks the balance and holds the port within a measured band of the
/// capture, so NeqSim repairing its flash fails the test rather than passing it quietly.
fn assert_consistent_with_its_enthalpy(
    out: &TrayOutcome,
    inlets: &[Stream],
    duty: f64,
    neqsim_temperature: f64,
    band: f64,
) {
    let target = flash_target(inlets, duty);
    let (mut moles, mut enthalpy) = (0.0, 0.0);
    for phase in [out.gas.as_ref(), out.liquid.as_ref()]
        .into_iter()
        .flatten()
    {
        moles += phase.n;
        enthalpy += phase.n * phase.h.value;
    }
    let weighted = enthalpy / moles;
    // `1e-4` J/mol rather than the flash's own `5.6e-7`: each outlet is a phase extract
    // rebuilt at (T, P), so the round trip through `Stream::from_pt` costs about `7e-6`
    // J/mol, measured on these rows.
    assert!(
        (weighted - target).abs() < 1.0e-4,
        "the stage's outlets weigh to {weighted} J/mol against the {target} J/mol its flash \
         was given; azoth's own PH flash is consistent to 5.6e-7 J/mol, measured"
    );
    // The band is that row's own enthalpy drift converted at its own `dT/dh`, so it is a
    // statement about NeqSim's solver rather than a tolerance chosen to fit: it fails the
    // moment NeqSim's flash improves, and fails the other way if the port's balance moves.
    assert!(
        (out.temperature.value - neqsim_temperature).abs() < band,
        "the port's temperature {} against NeqSim's {neqsim_temperature}, band {band} K: the \
         capture's number is the state NeqSim's flash landed on and not the one its enthalpy \
         implies",
        out.temperature.value
    );
}

fn relative(actual: f64, expected: f64, tolerance: f64, what: &str) {
    let scale = 1.0 + actual.abs().max(expected.abs());
    assert!(
        (actual - expected).abs() <= tolerance * scale,
        "{what}: {actual} vs {expected}, {} relative",
        (actual - expected).abs() / scale
    );
}

/// An absolute comparison, for a quantity that is a difference and can sit near zero.
fn absolute(actual: f64, expected: f64, tolerance: f64, what: &str) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{what}: {actual} vs {expected}, {} absolute",
        (actual - expected).abs()
    );
}

/// **One two-phase feed, which is the equilibrium split at the tray's own pressure.** With
/// no duty and no stated temperature the flash is at the mixed enthalpy, and the feed
/// arrives already at equilibrium - so the tray reproduces its own inlet split, measured at
/// `14` digits of the temperature.
#[test]
fn a_tray_splits_one_two_phase_feed() {
    let out = tray(
        &[feed(300.0, 20.0, 1.0)],
        None,
        None,
        watts(0.0),
        SideDraws::NONE,
        false,
    )
    .expect("the tray solves");

    absolute(
        out.temperature.value,
        300.0000000000223,
        1e-9,
        "tray temperature",
    );
    assert_eq!(out.pressure.value, 20.0e5, "the inlet's pressure");

    let gas = out.gas.expect("a two-phase feed has vapour");
    let liquid = out.liquid.expect("a two-phase feed has liquid");
    relative(gas.n, 0.3106080504317247, 1e-9, "gas n");
    relative(liquid.n, 0.6893919495682753, 1e-9, "liquid n");

    let gas_z = [
        0.24442604577015606,
        0.4215571797096646,
        0.27634084105196494,
        0.05767593346821439,
    ];
    let liquid_z = [
        0.034928319752273916,
        0.24523197053756648,
        0.4557150838538851,
        0.26412462585627466,
    ];
    for i in 0..4 {
        relative(gas.z[i], gas_z[i], 1e-6, "gas z");
        relative(liquid.z[i], liquid_z[i], 1e-6, "liquid z");
    }
    // The phase enthalpies are the state's at the tray's temperature and pressure. The
    // libraries' ideal-gas offset is the whole of the gap, and it is absolute.
    absolute(gas.h.value, 34.64421491209969, 0.1, "gas h");
    absolute(liquid.h.value, -13696.958208647447, 0.1, "liquid h");

    relative(gas.n + liquid.n, 1.0, 1e-12, "the mole balance");
}

/// **Two inlets, which is what a column stage actually receives**: vapour rising from the
/// stage below at 320 K and liquid falling from the stage above at 290 K. The tray mixes
/// them and flashes the mixture, so the answer is neither inlet's.
#[test]
fn a_tray_mixes_two_inlets_before_it_flashes() {
    let out = tray(
        &[feed(320.0, 20.0, 1.0), feed(290.0, 20.0, 1.0)],
        None,
        None,
        watts(0.0),
        SideDraws::NONE,
        false,
    )
    .expect("the tray solves");

    let inlets = [feed(320.0, 20.0, 1.0), feed(290.0, 20.0, 1.0)];
    let gas = out.gas.as_ref().expect("vapour");
    let liquid = out.liquid.as_ref().expect("liquid");
    // **A band and not agreement, because the split follows the temperature.** NeqSim's
    // flash landed 5.24 J/mol from the enthalpy it was asked for on this row, so its vapour
    // fraction is the one at *its* temperature; measured `3.3e-4` relative. The TP row
    // below is where the split itself is pinned, at a temperature both libraries agree on.
    relative(gas.n, 0.8723103388402136, 1e-3, "gas n");
    relative(liquid.n, 1.1276896611597869, 1e-3, "liquid n");
    relative(gas.n + liquid.n, 2.0, 1e-12, "the mole balance");
    // Neither inlet's temperature, and neither inlet's composition.
    assert!(out.temperature.value > 290.0 && out.temperature.value < 320.0);
    // This row's flash drifted `+5.24` J/mol, which is `0.016` K here.
    assert_consistent_with_its_enthalpy(&out, &inlets, 0.0, 307.3701867657665, 0.05);
}

/// A stated outlet temperature takes the class's other branch: a `TPflash` at that
/// temperature rather than a flash at the mixed enthalpy.
#[test]
fn a_stated_outlet_temperature_replaces_the_enthalpy_flash() {
    let out = tray(
        &[feed(300.0, 20.0, 1.0)],
        None,
        Some(kelvins(320.0)),
        watts(0.0),
        SideDraws::NONE,
        false,
    )
    .expect("the tray solves");

    assert_eq!(out.temperature.value, 320.0, "the stated temperature, held");
    relative(
        out.gas.expect("vapour").n,
        0.7134488551785124,
        1e-9,
        "gas n",
    );
    relative(
        out.liquid.expect("liquid").n,
        0.2865511448214876,
        1e-9,
        "liquid n",
    );
}

/// **The duty raises the flash's enthalpy, and it is divided by the flow.** The same 5000 W
/// is run against one mol/s and against two, and the temperature rises 15.3 K and 8.9 K
/// respectively - a *total* enthalpy, which is `PHflash(h, 0)`'s own reading. This is the
/// pair of rows that makes the tier's duty convention a measurement rather than a habit:
/// every other duty in this tier has been cased at one mol/s, where a molar reading and a
/// total one are the same number.
#[test]
fn a_trays_duty_is_a_total_enthalpy_and_the_flow_divides_it() {
    let one_in = [feed(300.0, 20.0, 1.0)];
    let two_in = [feed(300.0, 20.0, 2.0)];
    let one =
        tray(&one_in, None, None, watts(5000.0), SideDraws::NONE, false).expect("the tray solves");
    let two =
        tray(&two_in, None, None, watts(5000.0), SideDraws::NONE, false).expect("the tray solves");

    // The two-mol row's vapour fraction is the oracle for the division: a molar reading of
    // the duty would have left it at the one-mol row's 0.6018.
    // **The two vapour fractions are not oracled, and why is a measurement.** NeqSim's
    // flashes on these rows landed `-27.3` and `+179.5` J/mol from the enthalpy they were
    // asked for, which moves the vapour fraction by `1.0e-3` and `1.1e-2` relative - a band
    // too wide to be an oracle. What the pair proves instead is the division itself: a
    // *molar* reading of the duty would leave the two-mol row at the one-mol row's fraction,
    // and it is far below it.
    // The *fraction*, not the flow: the two-mol row vaporises more moles because it has
    // more, and a fraction of `0.454` against `0.603` is the division showing itself.
    let one_fraction = one.gas.as_ref().expect("vapour").n;
    let two_fraction = two.gas.as_ref().expect("vapour").n / 2.0;
    assert!(
        two_fraction < one_fraction - 0.1,
        "twice the flow on the same duty vaporises far less of itself: {two_fraction} against \
         {one_fraction}"
    );
    // A smaller rise, and **more** than half of it: halving the molar duty does not halve
    // the temperature rise, because the enthalpy is not linear in temperature - which is
    // exactly why the vapour fractions above are not linear in the duty either.
    let (one_rise, two_rise) = (one.temperature.value - 300.0, two.temperature.value - 300.0);
    assert!(
        two_rise < one_rise && two_rise > 0.5 * one_rise,
        "twice the flow on the same duty: {two_rise} K against {one_rise} K"
    );
    // `-27.3` J/mol drifted, which is `0.073` K here.
    assert_consistent_with_its_enthalpy(&one, &one_in, 5000.0, 315.31522969151587, 0.1);
    // `+179.5` J/mol drifted, which is `0.54` K here - the widest in the capture.
    assert_consistent_with_its_enthalpy(&two, &two_in, 5000.0, 308.84975823529163, 0.6);
}

/// The tray's own pressure overrides the inlet's, and the flash follows it down.
#[test]
fn a_tray_pressure_overrides_the_inlets() {
    let out = tray(
        &[feed(300.0, 20.0, 1.0)],
        Some(pascals(15.0e5)),
        None,
        watts(0.0),
        SideDraws::NONE,
        false,
    )
    .expect("the tray solves");

    let inlets = [feed(300.0, 20.0, 1.0)];
    assert_eq!(out.pressure.value, 15.0e5);
    // Measured `1.4e-5` relative, from the `-0.38` J/mol this row's flash drifted.
    relative(
        out.gas.as_ref().expect("vapour").n,
        0.3788143213813363,
        1e-4,
        "gas n",
    );
    // This row's flash drifted `-0.38` J/mol, which is `0.0012` K here.
    assert_consistent_with_its_enthalpy(&out, &inlets, 0.0, 291.58970263779287, 0.01);
}

/// **The absent phase, which is `None` here and a zero-flow stream in NeqSim.**
///
/// The class's own outlet getters look for a phase *type*: a subcooled feed has no gas
/// phase, so its vapour outlet is a zero-flow stream carrying the tray's composition whose
/// enthalpy is `-Infinity` - the capture's own number, because NeqSim divides an empty
/// system's enthalpy by its zero moles. A `Stream` cannot hold that, and fabricating a
/// finite one would be a state that looks real. So the absence is `None`.
///
/// **`beta` is not the split.** On the subcooled row the capture's `mixed_beta` is `1.0`
/// while the vapour outlet's flow is `0.0`: NeqSim's `getBeta()` answers one for any single
/// phase, so reading the split off it would give a subcooled feed all of its flow as vapour.
#[test]
fn an_absent_phase_is_none_and_beta_does_not_say_so() {
    let subcooled = tray(
        &[feed(220.0, 20.0, 1.0)],
        None,
        None,
        watts(0.0),
        SideDraws::NONE,
        false,
    )
    .expect("the tray");
    assert!(subcooled.gas.is_none(), "a subcooled feed has no vapour");
    let liquid = subcooled.liquid.expect("a subcooled feed is all liquid");
    relative(liquid.n, 1.0, 1e-12, "the liquid carries the whole flow");
    assert_eq!(liquid.z, feed(220.0, 20.0, 1.0).z, "the tray's composition");
    absolute(liquid.h.value, -20429.951874378614, 0.1, "liquid h");

    let superheated = tray(
        &[feed(400.0, 20.0, 1.0)],
        None,
        None,
        watts(0.0),
        SideDraws::NONE,
        false,
    )
    .expect("the tray");
    assert!(
        superheated.liquid.is_none(),
        "a superheated feed has no liquid"
    );
    relative(
        superheated.gas.expect("all vapour").n,
        1.0,
        1e-12,
        "the vapour carries the whole flow",
    );
}

/// No inlets is not a tray.
#[test]
fn a_tray_with_no_inlets_is_refused() {
    assert!(tray(&[], None, None, watts(0.0), SideDraws::NONE, false).is_err());
}

/// **A side draw is a split of the tray's own outlet phase, and the split is an identity.**
///
/// `SimpleTraySideDrawTest.gasSideDrawSplitsTrayOutletFlow` runs two trays on one feed - a
/// reference and one drawing a quarter of its vapour - and asserts that the reference's gas flow
/// equals the drawing tray's own gas plus its draw, and that the draw is the fraction of the
/// reference. Both are identities rather than numbers, so they hold whatever the flash answers,
/// and asserting them here is what makes the *routing* right rather than the arithmetic lucky.
#[test]
fn a_gas_side_draw_is_a_split_of_the_trays_own_vapour() {
    let inlets = [feed(300.0, 20.0, 1.0)];
    let plain = tray(&inlets, None, None, watts(0.0), SideDraws::NONE, false).expect("it solves");
    let drawn = tray(
        &inlets,
        None,
        None,
        watts(0.0),
        SideDraws {
            gas: 0.25,
            ..SideDraws::NONE
        },
        false,
    )
    .expect("it solves");

    let reference = plain.gas.as_ref().expect("the plain tray holds vapour").n;
    let kept = drawn.gas.as_ref().expect("the drawing tray keeps vapour").n;
    let side = drawn
        .gas_side_draw
        .as_ref()
        .expect("a quarter is withdrawn")
        .n;

    relative(
        kept + side,
        reference,
        1.0e-12,
        "the draw and the outlet against the whole",
    );
    relative(
        side,
        0.25 * reference,
        1.0e-12,
        "the draw against a quarter of it",
    );
    relative(
        plain.liquid.as_ref().expect("liquid").n,
        drawn.liquid.as_ref().expect("liquid").n,
        1.0e-12,
        "the liquid, which the gas draw does not touch",
    );
    // **A draw carries the tray's own state**: the same composition, temperature and pressure.
    let phase = drawn.gas.as_ref().expect("vapour");
    assert_eq!(drawn.gas_side_draw.as_ref().expect("the draw").z, phase.z);
    assert_eq!(
        drawn.gas_side_draw.as_ref().expect("the draw").t.value,
        drawn.temperature.value
    );
}

/// **The two liquid fractions are bounded together**, which is the class's own
/// `validateLiquidSplitFractions`: a side draw of 0.6 and a pumparound of 0.5 withdraw more than
/// the tray's liquid, and `SimpleTraySideDrawTest` asserts the refusal at the setter.
#[test]
fn a_liquid_draw_and_a_pumparound_may_not_withdraw_more_than_the_liquid() {
    let inlets = [feed(300.0, 20.0, 1.0)];
    let refused = tray(
        &inlets,
        None,
        None,
        watts(0.0),
        SideDraws {
            liquid: 0.6,
            pumparound: 0.5,
            ..SideDraws::NONE
        },
        false,
    );
    assert!(refused.is_err(), "the pair sums above one");

    // And a fraction outside [0, 1] is refused too, as `validateSideDrawFraction` refuses it.
    let out_of_range = tray(
        &inlets,
        None,
        None,
        watts(0.0),
        SideDraws {
            gas: 1.1,
            ..SideDraws::NONE
        },
        false,
    );
    assert!(out_of_range.is_err(), "a fraction above one");
}
