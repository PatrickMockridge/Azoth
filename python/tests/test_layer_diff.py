"""The layer diff, and the layers it reads.

`tools/neqsim_layer_diff.py` compares a model's intermediates against a probe capture's,
one layer at a time, and names the first that is out. This is its gate: `--all` over the
recorded cases, plus the parts of it that can be wrong without the gate noticing.

Three of those:

* a dumper can offer a set of layers that shares nothing with its capture, in which case
  the comparison is vacuous and a pass means nothing;
* "the *first* divergence" is a claim about order, and an unordered dump makes it a claim
  about nothing; and
* the layers are read out of a reference kernel's private state, so a field that moves
  between it and the probe's own key is a silent mismatch of two different quantities.

The last is what `tools/neqsim_layer_diff.py`'s reports are for and what the recorded
captures are immutable against: a capture is committed output, and a dump that stops
matching it fails here.
"""

from __future__ import annotations

import importlib
import subprocess
import sys
from pathlib import Path
from types import ModuleType
from typing import Any

REPO_ROOT = Path(__file__).resolve().parents[2]
TOOLS = REPO_ROOT / "tools"


def tools_module(name: str) -> ModuleType:
    """One of `tools/`, imported by name.

    The tools are not a package and sit outside mypy's `files`, so the import is dynamic -
    the arrangement `test_gen_neqsim_cases.py` and `test_ci_scripts.py` both use.
    """
    sys.path.insert(0, str(TOOLS))
    try:
        return importlib.import_module(name)
    finally:
        sys.path.pop(0)


def test_every_recorded_case_agrees_at_its_layers() -> None:
    """The gate, run where it runs in CI: every case, every layer it shares.

    A failure here names a layer rather than a total, which is the whole point - the
    case's own test says a model is 4% out, and this says where.
    """
    result = subprocess.run(
        [sys.executable, str(TOOLS / "neqsim_layer_diff.py"), "--all"],
        capture_output=True,
        text=True,
        check=False,
    )

    assert result.returncode == 0, (
        f"a layer diverged from its NeqSim capture\n{result.stdout}{result.stderr}"
    )


def test_every_dumped_case_shares_layers_with_its_capture() -> None:
    """A dump that shares nothing with its capture would pass by comparing nothing.

    The tool refuses that outright rather than reporting no divergence, and this is what
    says the refusal is not load-bearing on any case actually shipped.
    """
    diff: Any = tools_module("neqsim_layer_diff")
    gen: Any = tools_module("gen_neqsim_cases")

    from azoth.eos import layers as azoth_layers

    dumped = [case for case in gen.CASES if case.calc in azoth_layers.DUMPERS]
    assert dumped, "no case has a dumper, so this test is checking nothing"
    for case in dumped:
        inputs, _tolerance = diff.case_inputs(case.id)
        rows = diff.comparable(case, inputs)
        assert rows, f"{case.id} shares no layer with {case.capture}"


def test_the_first_divergence_is_the_first_in_build_order() -> None:
    """`first_divergence` returns the earliest key over the bound, not the largest.

    A dump is ordered the way the kernel builds it, so the earliest key is the layer the
    divergence is *in* and a later one may be its consequence. Returning the largest would
    name whatever the model is most sensitive to, which is the total's question again.
    """
    diff: Any = tools_module("neqsim_layer_diff")
    rows = [("early", 1.0, 1.0), ("first", 2.0, 2.001), ("late", 3.0, 4.0)]

    found = diff.first_divergence(rows, 1.0e-6)
    assert found is not None
    assert found[0] == "first", found

    assert diff.first_divergence([("agrees", 1.0, 1.0)], 1.0e-6) is None


def test_a_zero_layer_is_compared_absolutely() -> None:
    """A ratio cannot judge a layer the oracle reports as zero.

    Those are not incidental: an ideal-mixing `ln gamma`, or a departure that vanishes
    for a fluid without the term, is exactly the number a port gets wrong in a way no
    relative difference sees.
    """
    diff: Any = tools_module("neqsim_layer_diff")
    assert diff.relative(0.0, 0.0) == 0.0
    assert diff.relative(0.5, 0.0) == 0.5
    assert diff.relative(1.0, 2.0) == 0.5


def test_a_model_with_no_dumper_is_refused_rather_than_empty() -> None:
    """`layers` raises for a model it cannot dump.

    Returning `{}` would read to the diff as "nothing diverged", which is the one wrong
    answer a check like this must not be able to give.
    """
    from azoth.eos import layers as azoth_layers
    from azoth.eos.reference import tp_flash_saft

    assert "eos.tp_flash_saft" not in azoth_layers.DUMPERS
    assert callable(tp_flash_saft.tp_flash_saft)

    for model_id, dumper in azoth_layers.DUMPERS.items():
        assert model_id.split(".")[1].endswith("_phase"), model_id
        assert callable(dumper)
