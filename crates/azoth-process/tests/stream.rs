//! The stream's own properties, held to the process layer's oracle.
//!
//! `validation/neqsim/captures/process_stream_properties.tsv`, from
//! `ProcessProbe stream`. The numbers below are the capture's, at the precision it prints
//! them - `tools/oracle_sweep.py` reads Rust literals at this precision for exactly that
//! reason, so a pin move that moved one of these is reported rather than silently kept.
//!
//! **These are the quantities no port carries.** `h` is on the record and is measured by
//! every other capture; `s` is derived from `(T, P, z)` the same way, and density comes from
//! the state rather than from the record.
//!
//! # The pair is not a pair: the viscosity is not ported, and this is why
//!
//! **There is no `Stream::viscosity()`.** `eos.viscosity` is the PFCT *heavy-oil*
//! correlation, which its spec names as "NeqSim's liquid default" - and that is true of a
//! hydrocarbon liquid and false in general. On the capture's water row azoth gives
//! `5.302926524345695e-4` Pa·s where NeqSim's `getViscosity("kg/msec")` gives
//! `8.547585317553409e-4`, and water's viscosity at 300 K is `8.53e-4`: a crude-oil
//! correlation is 38% low on water, and NeqSim does not use it there.
//!
//! NeqSim's `SystemInterface.getViscosity` **dispatches on the phase's type**, and
//! `azoth-eos` carries one branch of that dispatch - the hydrocarbon liquid's. A `Stream`
//! has no phase type, so nothing here can choose, and returning the one correlation for
//! every stream would be a method named `viscosity` that is 38% wrong for water. The
//! dispatch is `azoth-eos`'s work; the four kernels that need a viscosity need it there,
//! not a wrapper here that hides which correlation answered.

use azoth_core::units::{kelvins, pascals};
use azoth_process::Stream;

fn stream(names: &[&str], z: &[f64], t: f64, p: f64, n: f64) -> Stream {
    Stream::from_pt(
        names.iter().map(|name| (*name).to_string()).collect(),
        z.to_vec(),
        n,
        pascals(p),
        kelvins(t),
    )
    .expect("the fluid resolves and flashes")
}

fn close(actual: f64, expected: f64, what: &str) {
    let scale = expected.abs().max(1.0);
    let relative = (actual - expected).abs() / scale;
    assert!(
        relative < 1e-4,
        "{what}: {actual} against the capture's {expected} is {relative:e} relative"
    );
}

#[test]
fn butane_liquid() {
    let s = stream(&["n-butane"], &[1.0], 300.0, 10.0e5, 1.0);
    close(
        s.molar_mass().unwrap().value,
        0.058122999999999994,
        "molar_mass",
    );
    close(s.mass_flow().unwrap(), 0.058122999999999994, "mass_flow");
    close(s.entropy().unwrap(), -70.09304618085763, "entropy");
    // The capture's `cubic_density`, not its `corrected_density`: NeqSim's
    // `getPhase(0).getDensity("kg/m3")` adds a volume correction on top of the cubic, and
    // `Stream` does not carry the component volume shifts that correction is built from.
    close(s.density().unwrap(), 601.2649457697453, "density");
    // No viscosity: see the module doc. The capture records NeqSim's for the day the
    // dispatch lands.
}

/// **The one row where the flash root disagrees, and it is not this accessor's.**
///
/// Three of the five rows match the capture's `cubic_density` in every digit. This one does
/// not: azoth's cubic root for 0.9/0.1 methane/n-butane at 320 K and 30 bar is `Z = 0.87296`
/// and NeqSim's is `Z = 0.91976`, so the same `M / (Z R T / P)` gives `26.157` and `24.826`.
///
/// **And NeqSim's own stream fluid disagrees with itself there.** Its
/// `getPhase(0).getZ()` is `0.91976` and its `getPhase(0).getDensity()` is `24.826126577519663`,
/// which is that root - while the *same object's* `getEntropy()/getTotalNumberOfMoles()` is
/// `-20.429350395514717`, which is the value a `Z` of `0.87296` gives. Azoth is
/// self-consistent at `0.87296` on all three: density `26.157222505467086`, entropy
/// `-20.429765016975505`, which is the capture's own entropy to `2e-5`. A fresh NeqSim fluid
/// at the same state has `Z = 0.91976` throughout and an entropy of `-19.843548034475464`.
///
/// The divergence is a root-selection question in `eos.pt_flash` for a state no `eos` case
/// pins, and locating it is that tier's business rather than this one's - `W2`'s layered
/// harness exists to say which layer moved. The assertion below is **that the divergence is
/// still there**, so that fixing `eos` fails this test and points back at this comment
/// instead of leaving a silently stale expectation.
#[test]
fn methane_butane_gas() {
    let s = stream(&["methane", "n-butane"], &[0.9, 0.1], 320.0, 30.0e5, 1.0);
    close(
        s.molar_mass().unwrap().value,
        0.020250999999999998,
        "molar_mass",
    );
    close(s.mass_flow().unwrap(), 0.020250999999999998, "mass_flow");
    close(s.entropy().unwrap(), -20.429765016975505, "entropy");

    let density = s.density().unwrap();
    assert!(
        (density - 24.826126577519663).abs() / 24.826126577519663 > 0.05,
        "the cubic root now agrees with the capture's {density}; if `eos.pt_flash` was fixed, \
         update this test and the comment above rather than the divergence"
    );
}

#[test]
fn water_liquid() {
    let s = stream(&["water"], &[1.0], 300.0, 1.0e5, 1.0);
    close(s.molar_mass().unwrap().value, 0.018015, "molar_mass");
    close(s.mass_flow().unwrap(), 0.018015, "mass_flow");
    close(s.entropy().unwrap(), -119.71049634390958, "entropy");
    // `848.2314181081963` is the cubic's; the capture's `corrected_density` is 983.93, and
    // the 14% between them is the volume correction `Stream` does not carry. Liquid water
    // is where that correction is largest, which is why this row is here.
    close(s.density().unwrap(), 848.2314181081963, "density");
}

#[test]
fn methane_co2_gas() {
    let s = stream(&["methane", "CO2"], &[0.7, 0.3], 300.0, 50.0e5, 1.0);
    close(s.molar_mass().unwrap().value, 0.0244331, "molar_mass");
    close(s.mass_flow().unwrap(), 0.0244331, "mass_flow");
    close(s.entropy().unwrap(), -26.80314870961037, "entropy");
    close(s.density().unwrap(), 56.64846407158036, "density");
}

/// **A two-phase stream has no one density, and this refuses rather than averaging.**
///
/// The capture's fifth row is the same fluid as `process.separator`'s feed, and NeqSim
/// reports `30.080813997105757` kg/m3 for it without saying that the number is a
/// volume-weighted average of a vapour and a liquid. A caller with a two-phase line needs
/// to say which density they mean; a plausible number that is neither is the failure this
/// refuses.
#[test]
fn a_two_phase_stream_refuses_a_density() {
    let s = stream(&["methane", "n-butane"], &[0.7, 0.3], 300.0, 20.0e5, 1.0);
    assert!(
        s.density().is_err(),
        "a two-phase stream has no one density"
    );
    // The three quantities that *are* defined for any stream still are.
    close(
        s.molar_mass().unwrap().value,
        0.028666999999999998,
        "molar_mass",
    );
    close(s.entropy().unwrap(), -26.962728758469062, "entropy");
}

/// **The measurement the module doc rests on, asserted so it cannot go stale.**
///
/// `eos.viscosity` is the PFCT heavy-oil correlation and returns `5.302926524345695e-4` Pa·s
/// for water at 300 K and 1 bar, where the capture has NeqSim's `8.547585317553409e-4` and
/// the true value is `8.53e-4`. When the per-phase dispatch lands and this stops being true,
/// this test fails and says so - which is the point, because the alternative is a comment
/// that quietly stops describing the code.
#[test]
fn the_heavy_oil_viscosity_is_still_wrong_for_water() {
    let (mixture, _) = azoth_eos::databank::mixture_of(&["water"], azoth_eos::Cubic::Pr, None)
        .expect("water resolves");
    let mu = azoth_eos::viscosity(&mixture, kelvins(300.0), pascals(1.0e5), &[1.0])
        .expect("the correlation answers")
        .mu
        .value;
    let relative = (mu - 8.547585317553409e-4).abs() / 8.547585317553409e-4;
    assert!(
        relative > 0.3,
        "the heavy-oil correlation now gives {mu} for water, {relative:e} from the capture's \
         8.547585317553409e-4. If the per-phase dispatch landed, read the module doc above \
         and add `Stream::viscosity()` rather than updating this expectation"
    );
}
