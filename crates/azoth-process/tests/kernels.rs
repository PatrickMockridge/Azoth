//! The kernels' balance invariants: moles and enthalpy are conserved.

use azoth_core::units::{kelvins, pascals};
use azoth_process::{Stream, mixer, separator, splitter};

fn close(a: f64, b: f64) {
    let scale = 1.0 + a.abs() + b.abs();
    assert!((a - b).abs() < 1e-9 * scale, "{a} vs {b}");
}

fn binary(z_methane: f64, n: f64, p: f64, t: f64) -> Stream {
    Stream::from_pt(
        vec!["methane".into(), "n-butane".into()],
        vec![z_methane, 1.0 - z_methane],
        n,
        pascals(p),
        kelvins(t),
    )
    .expect("methane / n-butane resolves and flashes")
}

#[test]
fn a_splitter_conserves_moles_and_state() {
    let feed = binary(0.5, 100.0, 1e5, 300.0);
    let outs = splitter(&feed, &[0.3, 0.7]).expect("splitter");

    assert_eq!(outs.len(), 2);
    close(outs[0].n + outs[1].n, feed.n);
    for out in &outs {
        assert_eq!(out.z, feed.z);
        assert_eq!(out.t, feed.t);
        assert_eq!(out.p, feed.p);
        assert_eq!(out.h, feed.h);
    }
}

#[test]
fn a_mixer_conserves_moles_and_enthalpy() {
    let a = binary(1.0, 50.0, 1e5, 300.0);
    let b = binary(0.0, 50.0, 1e5, 300.0);
    let m = mixer(&[a.clone(), b.clone()], None).expect("mixer");

    close(m.n, 100.0);
    close(m.z[0], 0.5);
    close(m.z[1], 0.5);
    close(m.h.value * m.n, a.h.value * a.n + b.h.value * b.n);
}

#[test]
fn a_separator_conserves_moles() {
    let feed = binary(0.5, 100.0, 5e5, 270.0);
    let (vapour, liquid) = separator(&feed, kelvins(270.0)).expect("separator");

    close(vapour.n + liquid.n, feed.n);
    if vapour.n > 0.0 && liquid.n > 0.0 {
        for i in 0..2 {
            close(
                vapour.z[i] * vapour.n + liquid.z[i] * liquid.n,
                feed.z[i] * feed.n,
            );
        }
    }
}
