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

#: The key a block's own label is kept under. `#` is not a character a probe key starts
#: with, so it cannot collide with a row.
LABEL_KEY = "#label"


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


def _heat_exchanger(inputs: Mapping[str, Any]) -> dict[str, float]:
    """`process.heat_exchanger`'s layers, from the reference's own factored arithmetic.

    **The rating's interior is here and the capture does not have it.** `NTU` is a
    package-private field of `HeatExchanger` with no getter, so the probe cannot print it -
    but `getDuty` and `getThermalEffectiveness` are two numbers it can, and those two pin
    every number below: the kept side is the one whose seeded swing is larger, so
    `duty = effectiveness * C_max * span`, and the effectiveness relation then leaves
    `C_min` as its only unknown. A wrong `c_min` cannot reach the same duty and the same
    effectiveness at the same `UA`, which is why the un-oracled four are worth dumping
    beside the two that are compared.
    """
    from azoth.process.reference.heat_exchanger import _states

    def optional(key: str) -> float | None:
        value = inputs.get(key)
        return None if value is None else float(value)

    states = _states(
        [str(name) for name in inputs["hot_components"]],
        float(inputs["hot_in_n"]),
        [float(v) for v in inputs["hot_in_z"]],
        float(inputs["hot_in_p"]),
        float(inputs["hot_in_t"]),
        [str(name) for name in inputs["cold_components"]],
        float(inputs["cold_in_n"]),
        [float(v) for v in inputs["cold_in_z"]],
        float(inputs["cold_in_p"]),
        float(inputs["cold_in_t"]),
        optional("ua"),
        str(inputs.get("flow_arrangement", "counterflow")),
        optional("hot_outlet_temperature"),
        optional("cold_outlet_temperature"),
    )
    # In the order the kernel computes them - the inlets, then the rating, then the
    # outlets - because the harness reports the *first* layer that moved, and the first one
    # to move should be the cause rather than the answer it produced.
    dump = {
        "hot_in_h": states.hot_in_h,
        "cold_in_h": states.cold_in_h,
    }
    if states.rating is not None:
        rating = states.rating
        dump |= {
            "hot_capacity": rating.hot_capacity,
            "cold_capacity": rating.cold_capacity,
            "c_min": rating.c_min,
            "c_max": rating.c_max,
            "capacity_ratio": rating.capacity_ratio,
            "ntu": rating.ntu,
            "effectiveness": rating.effectiveness,
            "duty_W": states.duty,
        }
    dump |= {
        "hot_out_h": states.hot_out_h,
        "cold_out_h": states.cold_out_h,
        "hot_out_T": states.hot_out_t.to("K").magnitude,
        "cold_out_T": states.cold_out_t.to("K").magnitude,
    }
    return dump


#: One model's layers, keyed by the model id a case names.
DUMPERS: dict[str, Dumper] = {
    "process.pump": _pump,
    "process.heat_exchanger": _heat_exchanger,
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
    #: or reordered would otherwise pair a dump against another state's numbers. Name the
    #: block's [`LABEL_KEY`] where there is one - a label is an identity, a state is a
    #: thing the comparison is already testing.
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
    # The exchanger's five cases sit on blocks 0, 1, 2, 5 and 6. Blocks 3 and 4 are
    # NeqSim's own `runSpecifiedStream` rows, which the port deliberately does not
    # reproduce - see `UNCASED_ROWS`.
    LayerCase(
        model="process.heat_exchanger",
        case="the_effectiveness_ntu_rating",
        capture="process_heat_exchanger.tsv",
        block=0,
        identified_by=(LABEL_KEY, "ua_rating_counterflow"),
    ),
    LayerCase(
        model="process.heat_exchanger",
        case="the_rating_answers_the_flow",
        capture="process_heat_exchanger.tsv",
        block=1,
        identified_by=(LABEL_KEY, "ua_rating_hot_flow_2"),
    ),
    LayerCase(
        model="process.heat_exchanger",
        case="the_arrangement_moves_less_heat",
        capture="process_heat_exchanger.tsv",
        block=2,
        identified_by=(LABEL_KEY, "ua_rating_parallelflow"),
    ),
    LayerCase(
        model="process.heat_exchanger",
        case="a_pinned_hot_outlet",
        capture="process_heat_exchanger.tsv",
        block=5,
        identified_by=(LABEL_KEY, "reference_pin_hot_350"),
    ),
    LayerCase(
        model="process.heat_exchanger",
        case="a_pinned_cold_outlet",
        capture="process_heat_exchanger.tsv",
        block=6,
        identified_by=(LABEL_KEY, "reference_pin_cold_320"),
    ),
)

#: Probe rows that are **deliberately not cases**, per capture, by the count they take.
#:
#: A probe row and a case are not the same thing: a row can be evidence for a divergence
#: rather than an oracle for a state the port reproduces. The exchanger's two
#: `out_temperature_pins_*` rows are that - NeqSim's own `runSpecifiedStream` lands on an
#: enthalpy its own `PHflash` places at another temperature, so the two pinned cases are
#: held to the `reference_pin_*` rows instead. Counting them here rather than leaving the
#: block count open keeps the check exact: a probe row added without a case still fails.
UNCASED_ROWS: dict[str, int] = {
    "process_heat_exchanger.tsv": 2,
}


def capture_blocks(capture: str) -> list[dict[str, str]]:
    """A capture's blocks, each as its `key=value` rows plus its label.

    The process probes print one block per case: a label line, the ports' records, the
    machine's own numbers, and a blank line. The label is the block's first line and is not
    a `key=value` row, so it is kept under [`LABEL_KEY`] rather than dropped. **It is what
    `LayerCase.identified_by` should name**, because it is the block's identity rather than
    a state: keying the pairing on a number would make an identification failure reachable
    by the very arithmetic the comparison exists to test, and two rows that differ only in
    an argument the ports do not carry - the exchanger's counterflow and parallel-flow rows
    - have no other distinguishing line between them.
    """
    text = (CAPTURES / capture).read_text(encoding="utf-8")
    blocks: list[dict[str, str]] = []
    for chunk in text.strip().split("\n\n"):
        rows: dict[str, str] = {}
        for line in chunk.strip().splitlines():
            line = line.strip()
            if not line:
                continue
            if "=" in line:
                key, _, value = line.partition("=")
                rows[key.strip()] = value.strip()
            else:
                rows[LABEL_KEY] = line
        if rows:
            blocks.append(rows)
    return blocks


def rows(model: str, inputs: Mapping[str, Any]) -> dict[str, float]:
    """A model's layers, under the capture's own key names."""
    dumper = DUMPERS.get(model)
    if dumper is None:
        raise KeyError(f"no layer dump for {model!r}; DUMPERS carries {sorted(DUMPERS)}")
    return dumper(inputs)
