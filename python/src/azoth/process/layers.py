"""The intermediates a unit operation's kernel already computes, as named rows.

    from azoth.process import layers
    layers.rows("process.pump", {"components": ["n-butane"], "inlet_n": 1.0, ...})

**A case compares a total, and a total that is 4% out says nothing about where.** A process
capture carries the machine's own intermediates beside the record - the shaft power, the
entropy the step produced - and those already exist inside the kernel on the way to the
answer. This exposes them by name so that the two sides can be put side by side and the
first layer that moved can be named.

# The sibling of `azoth.eos.layers`, and different in one way

The `eos` module reads a *reference* kernel's own solved state, because the reference is
where that kernel's arithmetic is written. A process model has two implementations too, and
its Python reference is the twin of the Rust kernel rather than the original - so this reads
the reference's intermediates through a factored entry point (`_states`, `_route`) rather
than re-deriving them. **A layer computed by a route the kernel does not take would be a
second implementation**, and a divergence in it would be a disagreement with ourselves.

# What counts as a layer

Only what the oracle also reports. `ProcessProbe` prints the record at every port and the
machine's own numbers beside it, and a layer whose key the capture has no counterpart for is
simply not compared - `tools/neqsim_layer_diff.py` says so rather than passing it.

That constraint is the point rather than a limitation: `Pump` exposes no getter for the
isentropic outlet, so the pump's interior is oracled through the two quantities it *does*
expose - the shaft power and the entropy the step produced - and the second of those is the
second-law statement of the same head.
"""

from __future__ import annotations

from collections.abc import Callable, Mapping
from dataclasses import dataclass
from pathlib import Path
from typing import Any

#: A model's declared inputs, keyed as a case states them, to its layers. The inputs are
#: bare numbers in the spec's units - the same ones `specs/cases/process/*.toml` carries -
#: so a caller states a state the way every other test in the suite does.
Dumper = Callable[[Mapping[str, Any]], dict[str, float]]

ROOT = Path(__file__).resolve().parents[4]
CAPTURES = ROOT / "validation" / "neqsim" / "captures"


def _pump(inputs: Mapping[str, Any]) -> dict[str, float]:
    """`process.pump`'s layers, from the reference's own factored arithmetic."""
    from azoth.process.reference.pump import _states

    states = _states(
        [str(name) for name in inputs["components"]],
        float(inputs["inlet_t"]),
        float(inputs["inlet_p"]),
        [float(v) for v in inputs["inlet_z"]],
        float(inputs["outlet_pressure"]),
        float(inputs["isentropic_efficiency"]),
    )
    n = float(inputs["inlet_n"])
    return {
        "inlet_h": states.h_in,
        "inlet_s": states.s_in,
        "outlet_h": states.h_out,
        "outlet_T": states.outlet_t.to("K").magnitude,
        # The two quantities `Pump` exposes. `kJ/molK` and `kW` are the capture's units:
        # the probe reads `getEntropyProduction` and `getPower`, so the dump has to be in
        # the same ones or the comparison is a unit conversion nobody wrote down.
        "power_kW": n * states.dh_actual / 1000.0,
        "entropy_production_kJ_per_molK": (states.s_out - states.s_in) / 1000.0,
    }


#: One model's layers, keyed by the model id a case names.
DUMPERS: dict[str, Dumper] = {
    "process.pump": _pump,
}


@dataclass(frozen=True)
class LayerCase:
    """One capture block, and the model case it is the oracle for."""

    model: str
    case: str
    capture: str
    #: Which block of the capture, counted from zero in the order the probe printed them.
    #: Positional because that is the convention the `eos` cases already use, and because a
    #: probe that named its own case ids would couple the oracle to the spec tree.
    block: int
    #: The one row the probe prints in that block that identifies it, as `(key, value)`.
    #: Checked rather than documented: the pairing is positional, so a probe row inserted
    #: or reordered would otherwise pair a dump against another state's numbers.
    identified_by: tuple[str, str]
    #: Layers compared on an **absolute** bound, with the bound, because the oracle is zero
    #: there by construction and a ratio divides by the vanishing thing.
    #:
    #: The same rule `python/tests/_helpers.py::_assert_diagnostic` applies to a solver
    #: residual, and for the same reason: a reversible pump's entropy production is zero by
    #: definition, so the capture has `-6.5780714209040525e-15` kJ/(mol*K) for the row at an
    #: efficiency of one and azoth has `9.18e-12`, which is 1400 times it and also zero. No
    #: relative tolerance can be met there, at any value.
    diagnostic: tuple[tuple[str, float], ...] = ()


#: The process cases a layer diff can read. **One entry per case**, so a probe row added
#: without a case - or a case whose probe row was reordered - is a mismatch
#: `python/tests/test_process_layer_diff.py` fails on rather than a silent mis-pairing.
LAYER_CASES: tuple[LayerCase, ...] = (
    LayerCase(
        model="process.pump",
        case="butane_5_to_20_bara_at_250_k",
        capture="process_pump.tsv",
        block=0,
        identified_by=("isentropic_efficiency", "0.75"),
    ),
    LayerCase(
        model="process.pump",
        case="the_work_divides_by_the_efficiency",
        capture="process_pump.tsv",
        block=1,
        identified_by=("isentropic_efficiency", "1.0"),
        # 1e-9 kJ/(mol*K) is 1e-6 J/(mol*K): four orders below any physical irreversibility
        # and four above the double-precision rounding of an enthalpy near 2.5e4 J/mol.
        diagnostic=(("entropy_production_kJ_per_molK", 1e-9),),
    ),
)


def capture_blocks(capture: str) -> list[dict[str, str]]:
    """A capture's blocks, each as its `key=value` rows.

    The process probes print one block per case: a label line, the ports' records, the
    machine's own numbers, and a blank line. The label is the block's first line and is not
    a `key=value` row, so it is dropped - what identifies a block for a reader is a line
    *inside* it, which is what `LayerCase.identified_by` is for.
    """
    text = (CAPTURES / capture).read_text(encoding="utf-8")
    blocks: list[dict[str, str]] = []
    for chunk in text.strip().split("\n\n"):
        rows: dict[str, str] = {}
        for line in chunk.strip().splitlines():
            if "=" in line:
                key, _, value = line.partition("=")
                rows[key.strip()] = value.strip()
        if rows:
            blocks.append(rows)
    return blocks


def rows(model: str, inputs: Mapping[str, Any]) -> dict[str, float]:
    """A model's layers, under the capture's own key names."""
    dumper = DUMPERS.get(model)
    if dumper is None:
        raise KeyError(f"no layer dump for {model!r}; DUMPERS carries {sorted(DUMPERS)}")
    return dumper(inputs)
