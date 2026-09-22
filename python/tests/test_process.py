"""The process layer's Python binding: streams, kernels and the checker.

`azoth.process.kernels` is a binding to `crates/azoth-process`, not a second
implementation, so these tests require the compiled extension and exercise the
binding's balance invariants rather than a reference.
"""

from __future__ import annotations

from pathlib import Path

import pytest

import azoth
from azoth import process
from azoth.process import kernels

REPO_ROOT = Path(__file__).resolve().parents[2]

pytestmark = pytest.mark.requires_rust


def _binary(z_methane: float, n: float, p: float, t: float) -> process.Stream:
    q = azoth.ureg.Quantity
    return process.Stream.from_pt(
        ["methane", "n-butane"],
        [z_methane, 1.0 - z_methane],
        n,
        q(p, "Pa"),
        q(t, "K"),
    )


def _close(actual: float, expected: float, tol: float = 1e-6) -> None:
    assert abs(actual - expected) < tol * (1.0 + abs(expected)), f"{actual} vs {expected}"


def test_a_splitter_conserves_moles_and_state() -> None:
    feed = _binary(0.5, 100.0, 1e5, 300.0)
    outs = kernels.splitter(feed, [0.3, 0.7])

    assert len(outs) == 2
    _close(outs[0].n + outs[1].n, feed.n)
    for out in outs:
        assert out.z == feed.z
        assert out.p == feed.p
        assert out.t == feed.t
        assert out.h == feed.h


def test_a_mixer_conserves_moles_and_enthalpy() -> None:
    a = _binary(1.0, 50.0, 1e5, 300.0)
    b = _binary(0.0, 50.0, 1e5, 300.0)
    m = kernels.mixer([a, b])

    _close(m.n, 100.0)
    _close(m.z[0], 0.5)
    _close(m.h.magnitude * m.n, a.h.magnitude * a.n + b.h.magnitude * b.n)


def test_a_separator_conserves_moles() -> None:
    feed = _binary(0.5, 100.0, 5e5, 270.0)
    vapour, liquid = kernels.separator(feed, azoth.ureg.Quantity(270.0, "K"))

    _close(vapour.n + liquid.n, feed.n)


def test_validate_holds_the_demo_and_rejects_a_double_fed_inlet() -> None:
    palette = str(REPO_ROOT / "specs" / "unit_ops")
    demo = str(REPO_ROOT / "specs" / "flowsheets" / "demo.toml")
    assert process.load_flowsheet(demo, palette) == []

    broken = (
        'id = "broken"\n'
        'name = "broken"\n'
        'feeds = ["f1", "f2"]\n'
        "products = []\n"
        "[[instances]]\n"
        'id = "p1"\n'
        'unit = "unit_ops.pump"\n'
        "[[connections]]\n"
        'from = "f1"\n'
        'to = "p1.inlet"\n'
        "[[connections]]\n"
        'from = "f2"\n'
        'to = "p1.inlet"\n'
    )
    diagnostics = process.validate(broken, palette)
    assert any("OverfedPort" in line for line in diagnostics)
