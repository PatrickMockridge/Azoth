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


def _splitter(inputs: Mapping[str, Any]) -> dict[str, float]:
    """`process.splitter`'s layers, from the reference's own factored arithmetic.

    **The branches' enthalpies are the layers where NeqSim and this library part**, and
    they are dumped so that the parting is measured rather than described. Two of the three
    are declared divergences; see `LAYER_CASES`.

    Nothing here reads a state the kernel did not compute. A branch's pressure, temperature
    and composition are the feed's *by construction* - the kernel copies them - so dumping
    them would be dumping the input, and the entropy is a function of `(T, P, z)` the
    splitter never forms, which would make this a second implementation of
    `Stream::entropy` rather than a layer of the split.
    """
    from azoth.process.reference.splitter import _route

    states = _route(
        [str(name) for name in inputs["components"]],
        float(inputs["feed_t"]),
        float(inputs["feed_p"]),
        [float(v) for v in inputs["feed_z"]],
        [float(f) for f in inputs["split_factors"]],
    )
    n = float(inputs["feed_n"])
    dump = {"feed_h": states.h_in}
    for index, fraction in enumerate(states.fractions):
        dump[f"products{index}_n"] = n * fraction
        dump[f"products{index}_h"] = states.h_in
    # The balance a split owes, under the same name the probe prints it by. It is the
    # diverging one: `run`'s branches carry the defect, so their weighted enthalpy is not
    # the feed's.
    dump["molar_enthalpy_out_weighted"] = sum(
        fraction * states.h_in for fraction in states.fractions
    )
    return dump


def _mixer(inputs: Mapping[str, Any]) -> dict[str, float]:
    """`process.mixer`'s layers, from the reference's own factored arithmetic.

    **`mixed_enthalpy_W` is the layer worth having.** An outlet's molar enthalpy is the
    weighted total over the total flow, and a mean of the right magnitude hides a wrong
    weighting - the two inlets here are at `-20689.86` and `-22027.34` J/mol, so weighting
    them by moles or by mass gives answers `0.2%` apart while the totals are `65 kJ` apart
    and cannot be confused. `Mixer.calcMixStreamEnthalpy` is the class's own name for it.

    The joined composition and the outlet pressure are on the record and are not dumped: a
    layer the port already carries would be the validation case over again.
    """
    from azoth.process.reference.mixer import _route

    states = _route(
        [str(name) for name in inputs["components"]],
        [float(v) for v in inputs["feed_n"]],
        [[float(v) for v in row] for row in inputs["feed_z"]],
        [float(v) for v in inputs["feed_p"]],
        [float(v) for v in inputs["feed_t"]],
        None if inputs.get("outlet_pressure") is None else float(inputs["outlet_pressure"]),
    )
    return {"mixed_enthalpy_W": states.h_total}


#: One model's layers, keyed by the model id a case names.
DUMPERS: dict[str, Dumper] = {
    "process.pump": _pump,
    "process.splitter": _splitter,
    "process.mixer": _mixer,
    "process.heat_exchanger": _heat_exchanger,
}

#: The process models whose kernel forms **nothing the port records do not carry**, and why.
#:
#: The harness exists to localise a wrong answer, and localising needs a layer between the
#: input and the output. Some kernels do not have one - their whole arithmetic *is* the
#: outlet record, so every number they form is on a port and the validation case already
#: holds it to the same capture. Declaring them here rather than leaving them absent is the
#: rule the palette follows too: nothing is silently missing, and a new process model has to
#: say which side of this it is on.
#:
#: **This is not "the model is too simple to check".** Each entry names the arithmetic and
#: what would have to exist for there to be a layer.
NO_INTERIOR: dict[str, str] = {
    "process.throttling_valve": (
        "the kernel is one ``Stream::from_ph`` - the outlet pressure the valve applied, and "
        "the flash it runs at the inlet's enthalpy. The applied pressure *is* the outlet "
        "record, and `getDeltaPressure` and `getEntropyProduction` are both differences of "
        "the two records' own fields, so nothing is left over. A layer would need the valve "
        "to size itself, which is the `MechanicalDesign.calcValveSize()` path the palette "
        "entry's `valve_opening` was withdrawn for"
    ),
    "process.cooler": (
        "the same kernel as ``process.heater`` and therefore the same answer: `Cooler` "
        "overrides `runTransient` and some getters and not `run`, and the probe's two "
        "captures are byte-identical. A dumper here would be a second dump of the heater's "
        "arithmetic under another name, which is the thing this dict exists to refuse"
    ),
    "process.heater": (
        "the kernel is one ``Stream::from_pt`` at the stated temperature or one "
        "``Stream::from_ph`` at the enthalpy a stated duty implies, at the pressure the drop "
        "leaves - and the outlet record is that flash's answer in full. The capture's extra "
        "line is ``duty_W``, which is ``Heater.run``'s own ``newH - oldH`` and therefore a "
        "difference of the two records' enthalpies over the inlet's flow; the model reports "
        "it as ``outlet_duty`` for exactly that reason. `getEnergyInput` is no way out of "
        "that: `run` overwrites the field with the same difference, so it reads back the "
        "number the record already implies rather than a layer behind it"
    ),
    "process.filter": (
        "the kernel is one ``enthalpy_at`` at the pressure the drop leaves, and the temperature "
        "is the inlet's - so the outlet record is that state in full. The capture's two extra "
        "lines are the applied drop, which the model reports as ``applied_drop``, and "
        "``Cv = sqrt(dP) / massFlow``, which the class computes from the two states and this "
        "port leaves out on purpose: it is a function of the applied drop and the inlet's mass "
        "flow, both of which the result already carries, and its unit is one this library has "
        "no dimension for. **The capture's third row is uncased rather than dumped**: the "
        "clamp lands the outlet at a microbar, where this library's cubic refuses to converge "
        "and NeqSim's extrapolates - so there is no layer to compare, and that is the finding"
    ),
    "process.separator": (
        "the kernel forms the flash's ``beta``, the two phase compositions and the moles the "
        "entrainment moves, and every one of the four is on the two outlet records - the "
        "vapour fraction is their flow ratio and the compositions are their ``z``. "
        "`Separator` exposes no getter for the state *before* the entrainment either: "
        "`getThermoSystem()` is the working system with the fraction already applied. "
        "**That state is not missing from the oracle, it is in the block beside it.** "
        "Entrainment is a move applied after the flash, and the capture's first row runs "
        "the same feed at the same conditions with no entrainment - so the fourth row's "
        "flash is the first row's answer exactly, and it is: the vapour leaves at "
        "`0.7773101311100368`, which is the first row's `0.8182211906421439` less the "
        "stated `0.05`, at that row's composition and enthalpy to `1.3e-15`. The flash's "
        "arithmetic is therefore already held by the three rows that state no entrainment, "
        "and reaching for it again here would be a second implementation of the ported "
        "arithmetic rather than a layer of it"
    ),
}


@dataclass(frozen=True)
class Divergence:
    """A layer the port **deliberately** does not reproduce, and how far from it it is.

    **Sometimes the layer that moved is supposed to.** `Splitter.run` does not conserve
    enthalpy: it clones the inlet fluid, zeroes it with `init(0)` and adds each component
    back phase by phase, so its outlets weight to `2.3` times what came in. A split changes
    no state, so the port answers with the feed's enthalpy - and a comparison against the
    capture's number can only ever read as a failure, at any tolerance.

    Reading it as a failure would be the wrong answer, and reading it as nothing would let
    the port drift. So the divergence is **declared and bounded**: the capture's value is
    asserted to be at least `at_least` times the port's. NeqSim repairing `Splitter.run`, or
    the port starting to agree, moves the ratio out of that band and fails - which is what
    makes this a measurement rather than a comment that quietly stops being true. The same
    shape as `crates/azoth-process/tests/stream.rs`'s
    `the_heavy_oil_viscosity_is_still_wrong_for_water`.
    """

    key: str
    #: The least the capture's value, divided by the port's, is allowed to be. A bound
    #: rather than an equality: the ratio is a ratio of two flashes' arithmetic, and the
    #: point is which of them it is near, not the digits.
    at_least: float
    #: Why the two are apart, in the shape the spec's `source` states it.
    reason: str


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
    #: Layers held to a declared [`Divergence`] instead of to agreement.
    divergence: tuple[Divergence, ...] = ()


#: The splitter's three divergent layers, declared once because both its rows carry the
#: same split and therefore the same bands.
#:
#: The capture has `-93092.87`, `-28142.87` and a weighted `-47627.87` against a feed of
#: `-20689.86`, which is `4.50`, `1.36` and `2.30` times the port's answer. The bounds sit
#: a little under each, so that the bands are the *ratios* and not their last digits: what
#: is asserted is that the capture is still several times away, not how many.
_SPLITTER_DIVERGENCE: tuple[Divergence, ...] = (
    Divergence(
        key="products0_h",
        at_least=3.0,
        reason="`Splitter.run` reaches its branch state through `init(0)` and "
        "`addComponent(int, double)` over every phase, so the first branch's enthalpy is "
        "not the feed's and the port's is: a split changes no state",
    ),
    Divergence(
        key="products1_h",
        at_least=1.2,
        reason="the same defect on the branch that carries most of the flow, where it is "
        "the smallest - the stale phase's share of the total is what moves",
    ),
    Divergence(
        key="molar_enthalpy_out_weighted",
        at_least=1.8,
        reason="the two branches weighted by their fractions, which is the balance a "
        "splitter owes and the row where NeqSim's 2.3 times the feed is visible",
    ),
)

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
    LayerCase(
        model="process.mixer",
        case="two_feeds_at_different_pressures",
        capture="process_mixer.tsv",
        block=0,
        identified_by=("specified_outlet_pressure_bara", "null"),
    ),
    LayerCase(
        model="process.mixer",
        case="a_stated_outlet_pressure_overrides_the_lowest_feed",
        capture="process_mixer.tsv",
        block=1,
        identified_by=("specified_outlet_pressure_bara", "8.0"),
    ),
    # The splitter's two rows carry the same numbers, so the pairing is the label's and
    # not a state's.
    LayerCase(
        model="process.splitter",
        case="a_two_way_split_of_a_liquid",
        capture="process_splitter.tsv",
        block=0,
        identified_by=("split_factors", "0.3 0.7"),
        divergence=_SPLITTER_DIVERGENCE,
    ),
    LayerCase(
        model="process.splitter",
        case="the_factors_are_normalised",
        capture="process_splitter.tsv",
        block=1,
        identified_by=("split_factors", "3.0 7.0"),
        divergence=_SPLITTER_DIVERGENCE,
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


def case_inputs(model: str, case: str) -> tuple[dict[str, Any], float]:
    """A case's declared inputs and tolerance, from the generated registry.

    The registry rather than the case file: `specs/cases/process/*.toml` is what a case
    *is*, and `azoth._models_gen` is its generated form - the one the reference kernels and
    the generated stub already read.
    """
    from azoth import _models_gen

    for declared in _models_gen.MODELS:
        if declared["id"] != model:
            continue
        for entry in declared["cases"]:
            if entry["id"] == case:
                return dict(entry["inputs"]), float(entry["tolerance"])
    raise KeyError(f"{model} has no case {case!r}")


@dataclass(frozen=True)
class LayerDiff:
    """One layer both sides have, and whether it is where it is declared to be.

    The three rules live here rather than in the gate that applies them, so that the gate
    and `tools/neqsim_layer_diff.py` cannot come to different answers about one layer.
    """

    key: str
    capture: float
    azoth: float
    #: Which rule holds it: `"case"` for the case's relative tolerance, `"absolute"` for a
    #: declared [`LayerCase.diagnostic`], `"divergence"` for a declared [`Divergence`].
    by: str
    #: The bound the rule compares against.
    bound: float
    #: What the rule measured - a relative gap, an absolute one, or a ratio.
    gap: float
    ok: bool

    def describe(self) -> str:
        """The layer, both numbers, and the rule it missed."""
        if self.by == "divergence":
            return (
                f".{self.key}: the capture's {self.capture} is {self.gap:.3f} times azoth's "
                f"{self.azoth}, under the {self.bound} this divergence is declared to be at "
                f"least"
            )
        if self.by == "absolute":
            return (
                f".{self.key}: {self.azoth} against the capture's {self.capture} is "
                f"{self.gap:.3e} absolute, over the {self.bound:.1e} this layer is compared on"
            )
        return (
            f".{self.key}: {self.azoth} against the capture's {self.capture} is "
            f"{self.gap:.3e} relative, over the case's {self.bound:.1e}"
        )


def compare(
    layer_case: LayerCase,
    inputs: Mapping[str, Any],
    *,
    tolerance: float | None = None,
) -> list[LayerDiff]:
    """Every layer the dump and the capture share, in the dump's own order, with its verdict.

    **In the dumper's order and not the capture's**, because "the first layer that moved" is
    a statement about the build order of the model - the reason the dumpers put the cause
    before the answer it produced, and the reason this returns the list rather than the
    first failure.
    """
    block = capture_blocks(layer_case.capture)[layer_case.block]
    dumped = rows(layer_case.model, inputs)
    absolute = dict(layer_case.diagnostic)
    divergence = {declared.key: declared for declared in layer_case.divergence}
    out: list[LayerDiff] = []
    for key, mine in dumped.items():
        if key not in block:
            continue
        theirs = float(block[key])
        if key in divergence:
            declared = divergence[key]
            ratio = abs(theirs) / abs(mine) if mine else float("inf")
            out.append(
                LayerDiff(
                    key=key,
                    capture=theirs,
                    azoth=mine,
                    by="divergence",
                    bound=declared.at_least,
                    gap=ratio,
                    ok=ratio >= declared.at_least,
                )
            )
        elif key in absolute:
            bound = absolute[key]
            out.append(
                LayerDiff(
                    key=key,
                    capture=theirs,
                    azoth=mine,
                    by="absolute",
                    bound=bound,
                    gap=abs(mine - theirs),
                    ok=abs(mine - theirs) <= bound,
                )
            )
        else:
            bound = tolerance if tolerance is not None else case_tolerance(layer_case)
            scale = abs(theirs)
            gap = abs(mine - theirs) / scale if scale else abs(mine - theirs)
            out.append(
                LayerDiff(
                    key=key,
                    capture=theirs,
                    azoth=mine,
                    by="case",
                    bound=bound,
                    gap=gap,
                    ok=gap <= bound,
                )
            )
    return out


def case_tolerance(layer_case: LayerCase) -> float:
    """The tolerance the case itself compares with."""
    return case_inputs(layer_case.model, layer_case.case)[1]


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
