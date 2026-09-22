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
"""

from __future__ import annotations

from typing import Any

import pytest

from azoth import _models_gen
from azoth.process import layers


def _case_inputs(model_id: str, case_id: str) -> tuple[dict[str, Any], float]:
    """A case's inputs and tolerance, from the generated registry."""
    for model in _models_gen.MODELS:
        if model["id"] != model_id:
            continue
        for case in model["cases"]:
            if case["id"] == case_id:
                return dict(case["inputs"]), float(case["tolerance"])
    raise AssertionError(f"{model_id} has no case {case_id!r}")


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
    """
    declared = next(m for m in _models_gen.MODELS if m["id"] == model)["cases"]
    captures = {c.capture for c in layers.LAYER_CASES if c.model == model}
    assert len(captures) == 1, f"{model} reads more than one capture: {sorted(captures)}"
    capture = captures.pop()
    blocks = layers.capture_blocks(capture)
    assert len(blocks) == len(declared), (
        f"{capture} has {len(blocks)} block(s) and {model} has {len(declared)} case(s); "
        f"the pairing is positional, so the two move together"
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

    dumped = layers.rows(layer_case.model, inputs)
    absolute = dict(layer_case.diagnostic)
    compared = 0
    for key, value in dumped.items():
        if key not in block:
            continue
        expected = float(block[key])
        if key in absolute:
            assert abs(value - expected) <= absolute[key], (
                f"{layer_case.model}::{layer_case.case}.{key}: {value} against the "
                f"capture's {expected} is {abs(value - expected):.3e} absolute, over the "
                f"{absolute[key]:.1e} this layer is compared on"
            )
        else:
            scale = abs(expected)
            relative = abs(value - expected) / scale if scale else abs(value - expected)
            assert relative <= tolerance, (
                f"{layer_case.model}::{layer_case.case}.{key}: {value} against the "
                f"capture's {expected} is {relative:.3e} relative, over the case's "
                f"{tolerance:.1e}"
            )
        compared += 1
    assert compared > 0, (
        f"{layer_case.model}'s dump shares no key with {layer_case.capture} block "
        f"{layer_case.block}; layers carries {sorted(dumped)} and the capture "
        f"{sorted(k for k in block if not k.endswith('_z'))}"
    )
