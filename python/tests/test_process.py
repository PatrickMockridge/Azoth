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
    q = azoth.ureg.Quantity
    feed = _binary(0.5, 100.0, 5e5, 270.0)
    vapour, liquid = kernels.separator(feed, q(0.0, "Pa"))

    _close(vapour.n + liquid.n, feed.n)
    # Energy too, and only because nothing was asked of the vessel: the flash is at the
    # feed's own temperature and pressure. A pressure drop moves it, which is a case's
    # business rather than an invariant's.
    _close(
        vapour.h.magnitude * vapour.n + liquid.h.magnitude * liquid.n,
        feed.h.magnitude * feed.n,
    )


def test_a_separator_carries_vapour_into_the_liquid() -> None:
    q = azoth.ureg.Quantity
    feed = _binary(0.5, 100.0, 5e5, 270.0)
    plain, _ = kernels.separator(feed, q(0.0, "Pa"))
    carried, liquid = kernels.separator(feed, q(0.0, "Pa"), gas_in_liquid=0.1)

    _close(carried.n, plain.n * 0.9)
    _close(carried.n + liquid.n, feed.n)


def test_validate_holds_the_demo_and_rejects_a_double_fed_inlet() -> None:
    palette = str(REPO_ROOT / "specs" / "unit_ops")
    demo = str(REPO_ROOT / "specs" / "flowsheets" / "demo.toml")
    assert process.load_flowsheet(demo, palette) == []

    broken = (
        'id = "broken"\n'
        'name = "broken"\n'
        "products = []\n"
        "[[inputs]]\n"
        'name = "f1"\n'
        'components = ["methane"]\n'
        "n = 1.0\n"
        "z = [1.0]\n"
        "P = 1.0e5\n"
        "T = 300.0\n"
        "[[inputs]]\n"
        'name = "f2"\n'
        'components = ["methane"]\n'
        "n = 1.0\n"
        "z = [1.0]\n"
        "P = 1.0e5\n"
        "T = 300.0\n"
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


def test_run_flowsheet_runs_the_shipped_flowsheet() -> None:
    """The Python surface is the executor's, and the shipped document runs through it.

    **What this test is allowed to pin, and what it is not.** The executor's arithmetic is
    Rust's and is checked where it lives - `tests/executor_calls_kernels.rs` holds each
    instance's outlet to that unit's own registered model, measured at `5.7e-12` relative on
    the pump. So a divergence from the capture here would not be the bridge's, and the one
    number pinned below is pinned as *the capture's*, with the band its measurement gives.
    """
    palette = str(REPO_ROOT / "specs" / "unit_ops")
    demo = (REPO_ROOT / "specs" / "flowsheets" / "demo.toml").read_text()

    # **No feeds: the document declares its own input**, which is what a self-contained
    # flowsheet means. `test_a_supplied_feed_replaces_the_document_s` is the override.
    result = process.run_flowsheet(demo, palette_dir=palette)

    assert result.flowsheet == "flowsheets.demo"
    assert result.converged
    assert result.iterations == 2, "a recycle takes a second observation before it converges"
    assert set(result.streams) == {
        "feed_1",
        "mix1.product",
        "p1.outlet",
        "hx1.outlet",
        "sep1.vapour",
        "sep1.liquid",
        "vapour_product",
    }

    # **The shipped document's recycle is switched off rather than closed**, which is measured
    # rather than incidental: at 20 bar and 320 K that feed is all vapour, so the tear carries
    # nothing, `Recycle.run` takes its low-flow branch, and the four residuals are *declared* zero
    # rather than measured. `active=False` beside `solved=True` is that state, and a reader who saw
    # only `solved` would take a loop carrying nothing for one that closed.
    #
    # `iterations` is the tear's own count and not the loop's: the class skips an inactive unit
    # rather than running it, so the second pass `report` still takes does not advance it. The
    # capture's `demo` row reads `recycle_iterations=1`, and `result.iterations` is 2.
    (tear,) = result.tears
    assert tear.stream == "recycle_1"
    assert tear.iterations == 1, "the tear's own count, not the loop's"
    assert result.iterations == 2, "the outer loop ran two passes"
    assert tear.solved
    assert not tear.active, "solved by being switched off"
    assert tear.residuals is not None
    assert (tear.residuals.flow, tear.residuals.composition) == (0.0, 0.0)
    assert (tear.residuals.temperature, tear.residuals.pressure) == (0.0, 0.0)
    assert result.streams["sep1.liquid"].n.magnitude == 0.0

    # `p1.outlet.T`, against `captures/process_flowsheet.tsv`'s `416.011567832174`. The measured
    # difference is `0.012 K`, and it is upstream of the execution layer rather than in it.
    pump_out = result.streams["p1.outlet"]
    _close(pump_out.n.magnitude, 1.0)
    _close(pump_out.p.to("bar").magnitude, 20.0)
    assert abs(pump_out.t.magnitude - 416.011567832174) < 0.05


def test_the_document_is_the_one_the_fields_were_read_from() -> None:
    """**One codec, and no second writer on the Python side.**

    `FlowsheetResult.document` is byte-identical to what the extension returned, so a
    front-end handed the document and a caller reading `.streams` cannot disagree about a
    quantity - there is one text and the dataclasses are a reading of it. A Python
    re-encode of the parsed fields would be a second writer with its own key order and its
    own number formatting, which is the drift `executor::json` exists to make impossible.
    """
    import json

    from azoth import _core

    palette = str(REPO_ROOT / "specs" / "unit_ops")
    demo = (REPO_ROOT / "specs" / "flowsheets" / "demo.toml").read_text()

    result = process.run_flowsheet(demo, palette_dir=palette)
    assert result.document == _core.run_flowsheet(demo, None, palette)

    raw = json.loads(result.document)
    assert raw["flowsheet"] == result.flowsheet
    assert raw["converged"] is result.converged
    assert raw["iterations"] == result.iterations
    assert set(raw["streams"]) == set(result.streams)
    for endpoint, record in raw["streams"].items():
        stream = result.streams[endpoint]
        assert record["n"]["magnitude_si"] == stream.n.magnitude

    # And the document is read-only by the name it is held through: a frozen result whose
    # mapping could be mutated in place would be frozen in name only.
    with pytest.raises(TypeError):
        result.streams["feed_1"] = result.streams["p1.outlet"]  # type: ignore[index]


def test_a_supplied_feed_replaces_the_document_s() -> None:
    """**The override, and the refusal beside it.**

    A document is self-contained, so a supplied feed replaces the value of a boundary it
    already declares rather than supplying a missing one. The measurement is the molar flow:
    the heater pins the separator's inlet temperature, so a hotter feed reaches the same
    state, while doubling the flow doubles what the separator hands the tear. A name the
    document does not declare is **refused** rather than added - a boundary nothing consumes
    is a wiring error, not an extra feed.
    """
    palette = str(REPO_ROOT / "specs" / "unit_ops")
    demo = (REPO_ROOT / "specs" / "flowsheets" / "demo.toml").read_text()

    as_written = process.run_flowsheet(demo, palette_dir=palette)
    declared = as_written.streams["sep1.liquid"].n.magnitude

    doubled = process.run_flowsheet(demo, {"feed_1": _binary(0.9, 2.0, 5e5, 300.0)}, palette)
    assert doubled.streams["sep1.liquid"].n.magnitude == pytest.approx(2.0 * declared)

    with pytest.raises(azoth.InvalidInputError, match="feed_9"):
        process.run_flowsheet(demo, {"feed_9": _binary(0.9, 1.0, 5e5, 300.0)}, palette)


def test_run_flowsheet_refuses_a_flowsheet_it_cannot_read() -> None:
    """The two refusals the parser and the order check make, each named rather than guessed."""
    palette = str(REPO_ROOT / "specs" / "unit_ops")
    demo = (REPO_ROOT / "specs" / "flowsheets" / "demo.toml").read_text()

    # `useGraphBasedExecution` is a flag, so there are two orders and not a menu.
    with pytest.raises(ValueError, match="insertion"):
        process.run_flowsheet(demo, palette_dir=palette, execution_order="kahn")


def _demo() -> str:
    return (REPO_ROOT / "specs" / "flowsheets" / "demo.toml").read_text()


def test_the_palette_comes_back_as_a_form_per_entry() -> None:
    """The declaration a unit-op window *is*, and the embedded bundle is the shipped palette."""
    entries = process.forms()
    assert len(entries) == 29
    assert sum(1 for entry in entries if entry["model"]) == 27
    assert sum(1 for entry in entries if entry["runnable"]) == 26

    pump = next(entry for entry in entries if entry["id"] == "unit_ops.pump")
    pressure = next(p for p in pump["parameters"] if p["name"] == "outlet_pressure")
    # The kind is what picks a widget, and it is the model's own declaration - not in the palette.
    assert pressure["kind"] == "quantity"
    assert pressure["unit"] == "Pa"
    assert pressure["dimension"] == "pressure"
    assert pressure["required"] is True

    efficiency = next(p for p in pump["parameters"] if p["name"] == "isentropic_efficiency")
    assert efficiency["kind"] == "quantity"
    assert efficiency["range"][0]["min"] == 0.0
    assert efficiency["range"][0]["min_inclusive"] is False
    assert efficiency["range"][0]["max"] == 1.0
    assert efficiency["range"][0]["rationale"]


def test_the_agent_tools_are_the_command_model() -> None:
    names = [tool["name"] for tool in process.tools()]
    assert len(names) == 14
    assert names[:4] == ["add_instance", "remove_instance", "connect", "disconnect"]

    add = process.tools()[0]["input_schema"]
    assert add["additionalProperties"] is False
    assert add["properties"]["command"]["const"] == "add_instance"
    assert len(add["properties"]["unit"]["enum"]) == 29


def test_a_session_holds_a_document_and_marks_its_values_stale() -> None:
    """The state `run_flowsheet` cannot give: an edit re-checks, and the values go older."""
    session = process.Session(_demo())
    assert session.ok is True
    assert session.dirty is True
    assert session.paths == []

    ran = session.run()
    assert ran["session"]["converged"] is True
    # The flags are read off the envelope, which is what a front-end is handed.
    assert ran["dirty"] is False, "the values are the document's now"
    assert "p1.outlet.P" in session.paths
    _close(session.value("p1.outlet.P"), 2.0e6, 1e-9)

    edited = session.apply(
        {
            "command": "set_parameter",
            "instance": "p1",
            "name": "outlet_pressure",
            "value": 4.0e6,
        }
    )
    assert edited["ok"] is True
    assert edited["dirty"] is True, "the values are older than the document"
    session.run()
    _close(session.value("p1.outlet.P"), 4.0e6, 1e-9)


def test_a_broken_edit_is_reported_by_the_check_that_follows_it() -> None:
    session = process.Session(_demo())
    envelope = session.apply({"command": "remove_instance", "id": "hx1"})

    assert envelope["ok"] is False
    codes = [diagnostic["code"] for diagnostic in envelope["diagnostics"]]
    assert codes == ["underfed_port", "underfed_port"]
    # Each diagnostic says what to put a mark on, and where in the text to find it.
    target = envelope["diagnostics"][0]["target"]
    assert target == {
        "kind": "handle",
        "role": "instance",
        "node": "p1",
        "port": "outlet",
        "index": None,
    }
    assert session.ok is False


def test_the_graph_is_the_editor_s_own_node_and_edge_document() -> None:
    session = process.Session(_demo())
    graph = session.graph
    assert [node["id"] for node in graph["nodes"]] == [
        "instance:mix1",
        "instance:p1",
        "instance:hx1",
        "instance:sep1",
        "input:feed_1",
        "product:vapour_product",
    ]
    assert len(graph["edges"]) == 6
    # A handle's id is the stream's own path, so an edge and the session agree without a mapping.
    loop_back = graph["edges"][-1]
    assert loop_back["data"]["kind"] == "recycle"
    assert loop_back["data"]["path"] == "recycle_1"
    assert loop_back["sourceHandle"] == "sep1.liquid"


def test_a_command_that_is_not_one_raises_and_the_document_that_is_not_one_raises() -> None:
    session = process.Session(_demo())
    with pytest.raises(ValueError, match="delete_everything"):
        session.apply({"command": "delete_everything"})
    with pytest.raises(ValueError):
        process.Session("this is not TOML")


def test_the_document_a_session_holds_writes_and_reads_back() -> None:
    session = process.Session(_demo())
    session.apply({"command": "add_product", "name": "purge"})
    document = process.Session(session.document)
    assert "purge" in document.document
    assert document.ok is True
