"""The process tier's layers, compared against the captures they came from.

**A case compares a total, and a total that is 4% out says nothing about where.** The
`eos` tier has had the sibling of this since the beginning -
`tools/neqsim_layer_diff.py` puts a model's layers against its capture key by key and names
the first that moved - and P11's W2 is where the process tier gets the same instrument.

What this checks is narrow and worth stating: that every layer of a unit operation whose
key the capture *also* prints agrees with it, and that at least one such key exists. A
layer the capture has no counterpart for is not compared - `Pump` exposes no getter for its
isentropic outlet, so that interior is reached through the shaft power and the entropy the
step produced - and a dump that shared no key with its capture would be a dumper reading
the wrong thing rather than a passing test.

**The three rules a layer may be held to are `layers.compare`'s**, not this file's: the
case's relative tolerance, a declared absolute bound where the oracle is zero by
construction, and a declared [`layers.Divergence`] where the port is deliberately not
NeqSim's answer. `tools/neqsim_layer_diff.py` applies the same three, so the gate and the
tool a developer reaches for when a case fails cannot come to different conclusions.
"""

from __future__ import annotations

from typing import Any

import pytest

from azoth import _models_gen
from azoth.process import layers


def _case_inputs(model_id: str, case_id: str) -> tuple[dict[str, Any], float]:
    """A case's inputs and tolerance, from the generated registry."""
    try:
        return layers.case_inputs(model_id, case_id)
    except KeyError as missing:
        raise AssertionError(f"{model_id} has no case {case_id!r}") from missing


def test_every_process_model_says_which_kind_it_is() -> None:
    """Each is either dumped or declared to have no interior, and never both.

    **The harness's own version of "nothing is silently absent".** A model added to
    `specs/models/process/` with no entry here would be one the layer diff never sees, and
    the silence would look exactly like agreement - which is the failure mode the whole
    tier's record is written against. So the partition has to be total, and a model that
    has no layer between its input and its output has to say so in `NO_INTERIOR` rather
    than be left out.
    """
    declared = {m["id"] for m in _models_gen.MODELS if str(m["id"]).startswith("process.")}
    dumped = set(layers.DUMPERS)
    interiorless = set(layers.NO_INTERIOR)
    assert not dumped & interiorless, (
        f"{sorted(dumped & interiorless)} is both dumped and declared to have no interior"
    )
    assert declared == dumped | interiorless, (
        f"neither dumped nor declared interiorless: {sorted(declared - dumped - interiorless)}; "
        f"declared for an id that is not a process model: "
        f"{sorted((dumped | interiorless) - declared)}"
    )


def _models() -> list[str]:
    seen: list[str] = []
    for layer_case in layers.LAYER_CASES:
        if layer_case.model not in seen:
            seen.append(layer_case.model)
    return seen


@pytest.mark.parametrize("model", _models(), ids=lambda name: name)
def test_a_capture_has_one_block_per_case(model: str) -> None:
    """The pairing is positional, so the two sequences have to be the same length.

    A probe row added without a case, or a case whose row was reordered, would otherwise
    pair a dump against another state's numbers - which compares and fails, or worse,
    compares and passes by coincidence.

    A block a case deliberately does not cover - a row that is evidence for a divergence
    rather than an oracle for a state the port reproduces - is declared in
    `LAYER_CASES`'s sibling [`layers.UNCASED_ROWS`] rather than excused by widening the
    comparison: the check stays exact, and the block count it holds to is the capture's
    own. The blocks the cases *do* name are checked one by one below, so the two together
    say which rows the cases cover and not merely how many.
    """
    declared = next(m for m in _models_gen.MODELS if m["id"] == model)["cases"]
    captures = {c.capture for c in layers.LAYER_CASES if c.model == model}
    assert len(captures) == 1, f"{model} reads more than one capture: {sorted(captures)}"
    capture = captures.pop()
    blocks = layers.capture_blocks(capture)
    uncased = layers.UNCASED_ROWS.get(capture, 0)
    assert len(blocks) == len(declared) + uncased, (
        f"{capture} has {len(blocks)} block(s), {model} has {len(declared)} case(s) and "
        f"{uncased} block(s) are declared uncased; the pairing is positional, so the "
        f"three move together"
    )
    named = [c.block for c in layers.LAYER_CASES if c.model == model]
    assert len(set(named)) == len(named), f"{model} names a block twice: {named}"
    assert len(set(range(len(blocks))) - set(named)) == uncased, (
        f"{capture} has {sorted(set(range(len(blocks))) - set(named))} block(s) no case "
        f"names and {uncased} are declared uncased"
    )


@pytest.mark.parametrize(
    ("layer_case"),
    layers.LAYER_CASES,
    ids=lambda case: f"{case.model}::{case.case}",
)
def test_a_models_layers_agree_with_its_capture(layer_case: layers.LayerCase) -> None:
    inputs, tolerance = _case_inputs(layer_case.model, layer_case.case)
    block = layers.capture_blocks(layer_case.capture)[layer_case.block]
    id_key, id_value = layer_case.identified_by
    assert block.get(id_key) == id_value, (
        f"{layer_case.capture} block {layer_case.block} is not the one "
        f"{layer_case.case} names: its {id_key} is {block.get(id_key)!r}, not {id_value!r}"
    )

    diffs = layers.compare(layer_case, inputs, tolerance=tolerance)
    assert diffs, (
        f"{layer_case.model}'s dump shares no key with {layer_case.capture} block "
        f"{layer_case.block}; layers carries {sorted(layers.rows(layer_case.model, inputs))} "
        f"and the capture {sorted(k for k in block if not k.endswith('_z'))}"
    )
    for diff in diffs:
        assert diff.ok, f"{layer_case.model}::{layer_case.case}{diff.describe()}"
