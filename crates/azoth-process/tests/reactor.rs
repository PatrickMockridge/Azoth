//! The reactor tier's internals: the axial stepper, on its own.
//!
//! **A stepper's oracle is an integral with a closed form.** Every other kernel in this crate
//! is pinned against a NeqSim capture, but a capture of a *stepper* would only be another
//! stepper's answer - so this holds `reactor::stepper::march` to arithmetic instead: the same
//! equations integrated exactly, at two step sizes, and the two schemes held to their own
//! orders. It is the one test in the tier that no NeqSim row can replace.

use azoth_process::Stream;
use azoth_process::reactor::stepper::{Scheme, march};

/// **The plug-flow reactor against its own capture, station by station.**
///
/// The capture is `validation/neqsim/captures/process_plug_flow_reactor.tsv`, row `rk4_default`:
/// methane `0.05`, oxygen `0.10`, nitrogen `0.85` at `600 K`, `5` bara and `1` mol/s; a power-law
/// combustion with `A = 1e4`, `Ea = 80000` J/mol and `ΔH = -802000` J/mol; five metres of
/// `0.10` m tube in a hundred RK4 steps with the property state refreshed every ten.
///
/// The answer is the whole profile, so this checks the class's own numbers at the inlet, the
/// first station and the outlet - a `dT/dz` that is right only at the ends would not be the same
/// curve.
#[test]
fn the_plug_flow_reactor_reproduces_the_captured_profile() {
    use azoth_core::units::{kelvins, pascals};
    use azoth_process::kernels::plug_flow_reactor::{ReactorSetup, plug_flow_reactor};
    use azoth_process::reactor::kinetic_reaction::KineticReaction;

    let feed = Stream::from_pt(
        ["methane", "oxygen", "nitrogen"]
            .iter()
            .map(|n| (*n).to_string())
            .collect(),
        vec![0.05, 0.10, 0.85],
        1.0,
        pascals(5.0e5),
        kelvins(600.0),
    )
    .expect("the three resolve");

    let mut reaction = KineticReaction::new("methanecombustion");
    reaction.add_reactant("methane", 1.0, 1.0);
    reaction.add_reactant("oxygen", 2.0, 1.0);
    reaction.add_product("CO2", 1.0);
    reaction.add_product("water", 2.0);
    reaction.pre_exponential_factor = 1.0e4;
    reaction.activation_energy = 80_000.0;
    reaction.heat_of_reaction = -802_000.0;

    let setup = ReactorSetup {
        reactions: vec![reaction],
        ..ReactorSetup::default()
    };
    let (_outlet, numbers, profile) =
        plug_flow_reactor(&feed, &setup).expect("the reactor marches");

    println!(
        "components={:?}\nconversion={} (neqsim 0.10666898902165212)\noutlet_t={} (neqsim 733.4597018618274)\n\
         drop={} (neqsim 1.1979021268260226e-5)\nresidence={} (neqsim 3.3649329465814315)\n\
         t1={} (neqsim 600.5421)\nt100={} (neqsim 733.4597)\nrate1={} (neqsim 5.428858e-2)",
        profile.components,
        numbers.conversion,
        numbers.outlet_temperature,
        numbers.pressure_drop_bar,
        numbers.residence_time,
        profile.temperatures[1],
        profile.temperatures[100],
        profile.rates[1],
    );

    assert_eq!(
        profile.components,
        vec!["methane", "oxygen", "nitrogen", "CO2", "water"],
        "the products are appended in the reaction's stoichiometry order"
    );
    assert_eq!(
        profile.positions.len(),
        101,
        "one station per step, plus the inlet"
    );

    // **The rate at a station is the tightest check, and it is the one that says the inputs are
    // right.** It reads the composition, the corrected molar density, the temperature and the
    // Arrhenius constant at once, with no integration to accumulate anything: `5e-6` relative.
    let rate = (profile.rates[1] - 5.428_858e-2).abs() / 5.428_858e-2;
    assert!(
        rate < 1e-4,
        "the first station's rate is {} against the capture's 5.428858e-2, {rate} relative",
        profile.rates[1]
    );

    // **The temperatures drift, and the drift is the heat capacity.** `dT/dz` divides by
    // `Cp * ΣF`, and `eos`'s Cp is `6.85e-4` below the class's (the divergence
    // `the_reactor_heat_capacity_is_the_eos_one_and_diverges_by_seven_per_ten_thousand` pins):
    // the first station is right to `7e-7` and the outlet is `4.2e-4` out, which is that
    // divergence integrated over a hundred steps and not an error in the march.
    let first = (profile.temperatures[1] - 600.5421).abs() / 600.5421;
    assert!(
        first < 1e-5,
        "the first station is {} K against the capture's 600.5421, {first} relative",
        profile.temperatures[1]
    );
    let outlet = (numbers.outlet_temperature - 733.459_701_861_827_4).abs() / 733.459_701_861_827_4;
    assert!(
        outlet < 1e-3,
        "outlet {} K against the capture's 733.4597018618274, {outlet} relative",
        numbers.outlet_temperature
    );

    // The conversion inherits the temperature drift through the rate, so it is looser than the
    // rate it comes from and tighter than the temperature it follows.
    let conversion =
        (numbers.conversion - 0.106_668_989_021_652_12).abs() / 0.106_668_989_021_652_12;
    assert!(
        conversion < 5e-3,
        "conversion {} against the capture's 0.10666898902165212, {conversion} relative",
        numbers.conversion
    );

    // **The residence time is the quirk, and it is the reason this test exists.**
    // `calculateResidenceTime` runs before the class's final state update, so it reads the
    // last *refreshed* state - a ninety-step temperature on a hundred-step march - and not the
    // outlet's. Using the outlet gives `3.2132` s; the capture says `3.3649`.
    let residence =
        (numbers.residence_time - 3.364_932_946_581_431_5).abs() / 3.364_932_946_581_431_5;
    assert!(
        residence < 1e-3,
        "residence {} s against the capture's 3.3649329465814315, {residence} relative - if this \
         is near 4.5 per cent, the port is reading the outlet state where the class reads the \
         last refresh",
        numbers.residence_time
    );

    // And the pressure row is the empty tube's, the class's default with no bed set: `4.8e-4`.
    let drop =
        (numbers.pressure_drop_bar - 1.197_902_126_826_022_6e-5).abs() / 1.197_902_126_826_022_6e-5;
    assert!(
        drop < 5e-3,
        "the drop is {} bar against the capture's 1.1979021268260226e-5, {drop} relative",
        numbers.pressure_drop_bar
    );
}

/// **The reactor's heat capacity exists, and it is `6.85e-4` below the class's.**
///
/// `PlugFlowReactor.calculateDerivatives` divides its heat generation by
/// `getCp("J/molK") * ΣF`, and the capture's `cp_probe` measures that value at the probe's own
/// feed state: `31.537822487323226`, which the probe also shows is NeqSim's own `(∂H/∂T)_P` to
/// eight significant figures.
///
/// **The surface is not missing.** `eos.molar_enthalpy_entropy` returns `cp` beside `h` and `s`,
/// as an ideal-gas part plus an analytic departure, and this pins what it gives at that state:
/// `31.516222490658407`, of which `0.053517` is the departure. Nothing has to be built.
///
/// **What the two do not do is agree**, and the gap is `6.85e-4` relative - larger than the
/// `0.00275` J/(mol·K) ideal-gas difference the process tier's enthalpy offset records, so it is
/// not that difference alone. It is the same family of library divergence as the offset itself:
/// invisible where a case pins a difference, and setting the reactor's temperature tolerance
/// because the whole of `dT/dz` is proportional to it. Asserting the divergence rather than a
/// band around the class's number is how it stays visible - the shape
/// `azoth-process/tests/stream.rs` uses for the aqueous viscosity that is still `38` per cent out.
#[test]
fn the_reactor_heat_capacity_is_the_eos_one_and_diverges_by_seven_per_ten_thousand() {
    use azoth_core::units::{kelvins, pascals};
    use azoth_eos::{Cubic, RootSide, databank, molar_enthalpy_entropy};

    let names = ["methane", "oxygen", "nitrogen"];
    let (mixture, ideal_gas) =
        databank::mixture_of(&names, Cubic::Pr, None).expect("the databank carries the three");
    let z = vec![0.05, 0.10, 0.85];
    let (t, p) = (kelvins(600.0), pascals(5.0e5));
    let reduced = mixture.reduced_parameters(t, p).expect("a state");
    let root = mixture
        .phase_state(&reduced, &z, RootSide::Vapour)
        .expect("a single vapour root")
        .z;
    let state = molar_enthalpy_entropy(&mixture, &ideal_gas, t, p, &z, root).expect("a state");

    let neqsim = 31.537_822_487_323_226;
    let relative = (state.cp.value - neqsim).abs() / neqsim;
    println!(
        "azoth cp={} (ideal {} + departure {}) neqsim={neqsim} relative={relative}",
        state.cp.value, state.cp_ideal.value, state.cp_departure.value
    );

    // It is the real-mixture capacity and not the ideal-gas one, so the departure is carried.
    assert!(
        state.cp_departure.value.abs() > 1e-3,
        "the departure is {} and would be zero if this were an ideal-gas Cp",
        state.cp_departure.value
    );
    // And it is the measured distance from the class, not a band around it: landing inside
    // `1e-4` would mean the divergence moved and this test must not silently accept that.
    assert!(
        (relative - 6.85e-4).abs() < 2e-5,
        "the libraries' heat capacities were measured {0:e} apart; this run is at {relative:e}. If an eos \
         change moved it, say so and re-measure the reactor's tolerance rather than widening this bound.",
        6.85e-4
    );
}

/// The exact solution of `dy/dz = -k y`, `y(0) = y0`.
fn decaying(z: f64, k: f64, y0: f64) -> f64 {
    y0 * (-k * z).exp()
}

/// March `dy/dz = -k y` over `length` in `steps` steps and return the final `y`.
fn integrate(k: f64, y0: f64, length: f64, steps: usize, scheme: Scheme) -> f64 {
    let mut state = vec![y0];
    let dz = length / steps as f64;
    march(
        &mut state,
        steps,
        dz,
        scheme,
        |state| vec![-k * state[0]],
        |_| {},
    );
    state[0]
}

/// **RK4 is fourth order and Euler is first, and the port's weights are the classical ones.**
///
/// If the `1:2:2:1` weights or the `1/6` were wrong, RK4 would fall to first or second order
/// and the ratio below would not be `16`.
#[test]
fn rk4_is_fourth_order_and_euler_is_first() {
    let (k, y0, length) = (0.7, 2.0, 1.0);
    let exact = decaying(length, k, y0);

    let coarse = (integrate(k, y0, length, 8, Scheme::Rk4) - exact).abs();
    let fine = (integrate(k, y0, length, 16, Scheme::Rk4) - exact).abs();
    let ratio = coarse / fine;
    assert!(
        (ratio - 16.0).abs() < 1.0,
        "halving the step should cut RK4's error sixteenfold, not {ratio} (coarse {coarse:e}, fine {fine:e})"
    );

    let euler_coarse = (integrate(k, y0, length, 8, Scheme::Euler) - exact).abs();
    let euler_fine = (integrate(k, y0, length, 16, Scheme::Euler) - exact).abs();
    let euler_ratio = euler_coarse / euler_fine;
    assert!(
        (euler_ratio - 2.0).abs() < 0.2,
        "halving the step should halve Euler's error, not scale it by {euler_ratio}"
    );

    // The two schemes are not the same arithmetic, which is the point of carrying both.
    assert!(
        euler_coarse > 100.0 * coarse,
        "Euler ({euler_coarse:e}) should be far behind RK4 ({coarse:e}) on the same step"
    );
}

/// **A march terminates on its own step count, and the length is the steps times the step.**
///
/// The class computes `dz = length / numberOfSteps` and loops exactly `numberOfSteps` times, so
/// the profile ends at `length` and nowhere else.
#[test]
fn the_profile_ends_at_the_length_it_was_given() {
    let mut state = vec![0.0];
    let mut stations = vec![0.0];
    let (length, steps) = (2.5, 5);
    let dz = length / steps as f64;

    march(
        &mut state,
        steps,
        dz,
        Scheme::Rk4,
        |_| vec![1.0],
        |state| stations.push(state[0]),
    );

    assert_eq!(
        stations.len(),
        steps + 1,
        "one station per step, plus the start"
    );
    // `dy/dz = 1` from `y(0) = 0` is `y = z`, so the state is the position.
    assert!(
        (state[0] - length).abs() < 1e-12,
        "the march should end at {length}, not {}",
        state[0]
    );
    assert!(
        (stations[steps] - length).abs() < 1e-12,
        "and the last station should be the last state"
    );
}

/// **`settle` runs after every advance, and what it writes is what the next step reads.**
///
/// The reactor clamps its molar flows at zero and its pressure at `0.1` there, and overrides an
/// isothermal reactor's temperature - so a `settle` that did not take effect before the next
/// derivative would silently change the answer.
#[test]
fn settle_rewrites_the_state_before_the_next_step() {
    let mut state = vec![0.0];
    let mut saw = Vec::new();

    // `dy/dz = 1` makes the state the step index; `settle` then doubles it, so each step starts
    // from twice the last. Holding 1, 3, 7, 15 is that rule and not the bare march's 1, 2, 3.
    march(
        &mut state,
        4,
        1.0,
        Scheme::Euler,
        |state| {
            saw.push(state[0]);
            vec![1.0]
        },
        |state| state[0] *= 2.0,
    );

    assert_eq!(
        saw,
        vec![0.0, 2.0, 6.0, 14.0],
        "the next derivative must see the settled state"
    );
    assert_eq!(state[0], 30.0, "and the last settle is the answer");
}

/// **The scheme name follows the class's setter, fallback and all.**
///
/// `PlugFlowReactor.setIntegrationMethod` matches `"EULER"` and answers RK4 otherwise - an
/// unknown name is the default rather than a refusal, which is a behaviour and not an
/// oversight.
#[test]
fn an_unknown_scheme_name_is_the_default_and_not_a_refusal() {
    assert_eq!(Scheme::named("EULER"), Scheme::Euler);
    assert_eq!(Scheme::named("RK4"), Scheme::Rk4);
    assert_eq!(Scheme::named("runge"), Scheme::Rk4);
    assert_eq!(
        Scheme::named("euler"),
        Scheme::Rk4,
        "the class's match is case-sensitive"
    );
}
