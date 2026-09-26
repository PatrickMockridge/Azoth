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


def _compressor(inputs: Mapping[str, Any]) -> dict[str, float]:
    """`process.compressor`'s layers, from the reference's own factored arithmetic.

    **The isentropic outlet is in the dump and not in the capture**, and that is worth
    stating: `Compressor` keeps its reversible enthalpy in a local, so the probe cannot print
    it - but the capture's efficiency-of-one row *is* that state, which is why that row is in
    the case set at all. A wrong isentropic step cannot reach both rows' enthalpies and both
    rows' entropy productions at once.
    """
    from azoth.process.reference.compressor import _states

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
        "isentropic_h": states.h_isentropic,
        "outlet_h": states.h_out,
        "outlet_T": states.outlet_t.to("K").magnitude,
        # `getPower` is `dH = h_out - h_in` over the machine, which is what the class's own
        # `power_kW` line prints - in kW against the record's J/mol.
        "power_kW": n * states.dh_actual / 1000.0,
        "entropy_production_kJ_per_molK": (states.s_out - states.s_in) / 1000.0,
    }


def _expander(inputs: Mapping[str, Any]) -> dict[str, float]:
    """`process.expander`'s layers, from the reference's own factored arithmetic.

    The same six as the compressor's plus the isentropic enthalpy, and the same reasoning:
    `Expander.run`'s second-law statement is the only other number it exposes, and the
    efficiency-of-one row is where the reversible state is the answer.
    """
    from azoth.process.reference.expander import _states

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
        "isentropic_h": states.h_isentropic,
        "outlet_h": states.h_out,
        "outlet_T": states.outlet_t.to("K").magnitude,
        # Negative, and it should be: `getPower` is the work the machine does, and an
        # expander produces it.
        "power_kW": n * states.dh_actual / 1000.0,
        "entropy_production_kJ_per_molK": (states.s_out - states.s_in) / 1000.0,
    }


def _pipe(inputs: Mapping[str, Any]) -> dict[str, float]:
    """`process.pipe`'s layers, from the reference's own factored arithmetic.

    **Every number the capture prints about the line's interior is here**, including the
    pressure drop in the capture's own unit: the probe reports `bara` and the arithmetic is
    in pascals, so a dump that crossed the SI figure would be a unit conversion nobody wrote
    down. `passes` is the one layer with no counterpart - the capture cannot print how many
    passes the class's loop took - and it is dumped anyway so that the *path* is visible
    beside the answer it produced.
    """
    from azoth.process.reference.pipe import _route

    states = _route(
        [str(name) for name in inputs["components"]],
        float(inputs["inlet_n"]),
        [float(v) for v in inputs["inlet_z"]],
        float(inputs["inlet_p"]),
        float(inputs["inlet_t"]),
        float(inputs["length"]),
        float(inputs["diameter"]),
        float(inputs["roughness"]),
    )
    return {
        "velocity_m_per_s": states.velocity,
        "reynolds_number": states.reynolds,
        "friction_factor": states.friction_factor,
        "pressure_drop_bara": (float(inputs["inlet_p"]) - states.p_out) / 1.0e5,
        "passes": float(states.passes),
    }


def _manifold(inputs: Mapping[str, Any]) -> dict[str, float]:
    """`process.manifold`'s layers, from the two routes it composes.

    **The mixture's own rows and each branch's flow**, which is what the manifold adds to
    its two children: `Manifold.run` is `localmixer.run()` then `localsplitter.run()` over
    the mixture, and the capture prints the mixture beside the branches. The branch
    *enthalpies* are deliberately not dumped: `process.splitter`'s own case declares them
    as a divergence, and the manifold's branches carry the same one - stated once, in the
    model that owns it, rather than twice.
    """
    from azoth.process.reference.manifold import DEFAULT_MINIMUM_FLOW_KG_PER_HOUR, _mass_flow
    from azoth.process.reference.mixer import _route as mix_route
    from azoth.process.reference.splitter import _route as split_route

    components = [str(name) for name in inputs["components"]]
    feed_n = [float(value) for value in inputs["feed_n"]]
    feed_z = [[float(x) for x in row] for row in inputs["feed_z"]]
    feed_p = [float(value) for value in inputs["feed_p"]]
    feed_t = [float(value) for value in inputs["feed_t"]]
    factors = [float(value) for value in inputs["split_factors"]]

    keep = [
        index
        for index in range(len(feed_n))
        if _mass_flow(components, feed_z[index], feed_n[index]) * 3600.0
        > DEFAULT_MINIMUM_FLOW_KG_PER_HOUR
    ]
    mixture = mix_route(
        components,
        [feed_n[index] for index in keep],
        [feed_z[index] for index in keep],
        [feed_p[index] for index in keep],
        [feed_t[index] for index in keep],
        None,
    )
    branches = split_route(
        components,
        mixture.product_t.to("K").magnitude,
        mixture.pressure,
        list(mixture.z),
        factors,
    )

    dump = {
        "product_n": mixture.n_total,
        "product_h": mixture.h_out,
        "product_T": mixture.product_t.to("K").magnitude,
    }
    for index, fraction in enumerate(branches.fractions):
        dump[f"products{index}_n"] = mixture.n_total * fraction
    return dump


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
def _shortcut_distillation_column(inputs: Mapping[str, Any]) -> dict[str, float]:
    """`process.shortcut_distillation_column`'s layers, from the reference's own arithmetic.

    **The layers are the K-values and the relative volatilities, because everything
    downstream of `alpha` is rearrangement.** Fenske, Underwood, Molokanov and Kirkbride are
    closed form in `alpha_i = K_i / K_HK`, so a divergence in the stage count or the feed
    tray is located either in the flash that produced the K-values or in the correlations -
    and this dump separates the two. It is also the row where the port and NeqSim are closest,
    measured at `3e-12` on the captured propane/n-butane row, so a gap here is a plumbing bug
    rather than a physics disagreement.

    `alpha_lk_hk_from_flash` is dumped beside the per-component vectors rather than instead of
    them because the probe prints it twice - once from the flash it re-runs, once from the
    class's own `getRelativeVolatility` - and the two agreeing is what makes the capture's
    K-values the class's own rather than a second flash's.
    """
    from azoth.process.reference.shortcut_distillation_column import _states

    components = [str(name) for name in inputs["components"]]
    states = _states(
        components,
        float(inputs["feed_n"]),
        float(inputs["feed_t"]),
        float(inputs["feed_p"]),
        [float(v) for v in inputs["feed_z"]],
        str(inputs["light_key"]),
        str(inputs["heavy_key"]),
        float(inputs["light_key_recovery_distillate"]),
        float(inputs["heavy_key_recovery_bottoms"]),
        float(inputs["reflux_ratio_multiplier"]),
        None if inputs.get("condenser_pressure") is None else float(inputs["condenser_pressure"]),
        None if inputs.get("reboiler_pressure") is None else float(inputs["reboiler_pressure"]),
    )
    layers: dict[str, float] = {"alpha_lk_hk_from_flash": states.relative_volatility}
    for index, name in enumerate(components):
        layers[f"k_{name}"] = states.k_values[index]
        layers[f"alpha_{name}"] = states.alpha[index]
    return layers


def _component_splitter(inputs: Mapping[str, Any]) -> dict[str, float]:
    """`process.component_splitter`'s layers, from the reference's own arithmetic.

    **The flows, temperatures and compositions are compared; the two enthalpies are declared
    divergences**, for the reason the sets beside `LAYER_CASES` give. The capture's own
    units, as bare magnitudes: the harness compares numbers.
    """
    from azoth.core.units import ureg
    from azoth.process.reference.component_splitter import component_splitter

    result = component_splitter(
        [str(name) for name in inputs["components"]],
        ureg.Quantity(float(inputs["feed_n"]), "mol/s"),
        [float(v) for v in inputs["feed_z"]],
        ureg.Quantity(float(inputs["feed_p"]), "Pa"),
        ureg.Quantity(float(inputs["feed_t"]), "K"),
        [float(v) for v in inputs["split_factors"]],
    )
    return {
        "overhead_n": result.overhead_n.to("mol/s").magnitude,
        "overhead_T": result.overhead_t.to("K").magnitude,
        "overhead_h": result.overhead_h.to("J/mol").magnitude,
        "bottoms_n": result.bottoms_n.to("mol/s").magnitude,
        "bottoms_T": result.bottoms_t.to("K").magnitude,
        "bottoms_h": result.bottoms_h.to("J/mol").magnitude,
    }


def _absorption_column(inputs: Mapping[str, Any]) -> dict[str, float]:
    """`process.absorption_column`'s layers: the profile, under the capture's own key names.

    The same layer as the column's and for the same reason - a column's response is a profile
    rather than a scalar - and the capture prints the same four keys per tray, so the first tray
    that moved is the layer the diff names.
    """
    from azoth.process.reference.distillation_column import _feed, _states

    components = [str(name) for name in inputs["gas_components"]]
    states = _states(
        components,
        float(inputs["gas_n"]),
        [float(v) for v in inputs["gas_z"]],
        float(inputs["gas_t"]),
        float(inputs["gas_p"]),
        int(inputs["number_of_stages"]),
        0,
        False,
        False,
        float(inputs["top_pressure"]),
        float(inputs["bottom_pressure"]),
        None,
        None,
        float(inputs["temperature_tolerance"]),
        int(inputs["max_iterations"]),
        None,
        None,
        str(inputs["solver_type"]) if "solver_type" in inputs else None,
        top_feed=_feed(
            [str(name) for name in inputs["solvent_components"]],
            float(inputs["solvent_n"]),
            [float(v) for v in inputs["solvent_z"]],
            float(inputs["solvent_t"]),
            float(inputs["solvent_p"]),
        ),
        tray_temperatures=(
            None
            if "tray_temperatures" not in inputs
            else tuple(float(v) for v in inputs["tray_temperatures"])
        ),
    )
    layers: dict[str, float] = {}
    for i in range(len(states.tray_temperature)):
        layers[f"tray{i}_temperature_K"] = states.tray_temperature[i]
        layers[f"tray{i}_pressure_bara"] = states.tray_pressure[i] / 1.0e5
        layers[f"tray{i}_gas_n"] = states.tray_gas_n[i]
        layers[f"tray{i}_liquid_n"] = states.tray_liquid_n[i]
    return layers


def _stripping_column(inputs: Mapping[str, Any]) -> dict[str, float]:
    """`process.stripping_column`'s layers: the profile, under the capture's own key names.

    The same layer as the absorber's, because the class is the same machine: `StrippingColumn`
    renames the two inlets and the two products and adds no equations.
    """
    from azoth.process.reference.distillation_column import _feed, _states

    components = [str(name) for name in inputs["stripping_gas_components"]]
    states = _states(
        components,
        float(inputs["stripping_gas_n"]),
        [float(v) for v in inputs["stripping_gas_z"]],
        float(inputs["stripping_gas_t"]),
        float(inputs["stripping_gas_p"]),
        int(inputs["number_of_stages"]),
        0,
        False,
        False,
        float(inputs["top_pressure"]),
        float(inputs["bottom_pressure"]),
        None,
        None,
        float(inputs["temperature_tolerance"]),
        int(inputs["max_iterations"]),
        None,
        None,
        str(inputs["solver_type"]) if "solver_type" in inputs else None,
        top_feed=_feed(
            [str(name) for name in inputs["rich_liquid_components"]],
            float(inputs["rich_liquid_n"]),
            [float(v) for v in inputs["rich_liquid_z"]],
            float(inputs["rich_liquid_t"]),
            float(inputs["rich_liquid_p"]),
        ),
        tray_temperatures=(
            None
            if "tray_temperatures" not in inputs
            else tuple(float(v) for v in inputs["tray_temperatures"])
        ),
    )
    layers: dict[str, float] = {}
    for i in range(len(states.tray_temperature)):
        layers[f"tray{i}_temperature_K"] = states.tray_temperature[i]
        layers[f"tray{i}_pressure_bara"] = states.tray_pressure[i] / 1.0e5
        layers[f"tray{i}_gas_n"] = states.tray_gas_n[i]
        layers[f"tray{i}_liquid_n"] = states.tray_liquid_n[i]
    return layers


def _spec_of(inputs: Mapping[str, Any], which: str) -> Any:
    """One end's specification from a case's own inputs, or ``None`` where none is stated."""
    from azoth.process.reference.distillation_column import Specification

    kind = inputs.get(f"{which}_specification_type")
    if kind is None:
        return None
    return Specification(
        kind=str(kind),
        target=float(inputs[f"{which}_specification_target"]),
        component=(
            None
            if inputs.get(f"{which}_specification_component") is None
            else str(inputs[f"{which}_specification_component"])
        ),
    )


def _distillation_column(inputs: Mapping[str, Any]) -> dict[str, float]:
    """`process.distillation_column`'s layers, from the reference's own factored arithmetic.

    **The profile is the layer, and it is the whole answer.** A column's response is not a
    scalar that a divergence can be attributed to; it is a temperature and two traffic rates on
    every tray, and the capture prints exactly those. So the dump is the profile under the
    capture's own key names - `tray3_temperature_K` and not `tray_temperature[3]` - and the
    first tray that moved is the layer the diff names.

    That is also why this model needs no divergence declaration: the two solvers reach the same
    fixed point one iteration apart, so the profile is compared directly rather than with a
    declared band.
    """
    from azoth.process.reference.distillation_column import _states

    # **Every declared input the case states is passed on, the section and the draws included.**
    # A dump that dropped one would compare a different column's profile against the capture and
    # call the difference a divergence: measured, the side-draw row's dump without its draw is
    # `4.9e-2` from the capture on tray 0.
    section = None
    if inputs.get("reactive"):
        start, end = inputs.get("reactive_start_tray"), inputs.get("reactive_end_tray")
        section = (-1, -1) if start is None or end is None else (int(start), int(end))
    draws = None
    if any(
        name in inputs
        for name in (
            "gas_side_draw_fractions",
            "liquid_side_draw_fractions",
            "pumparound_fractions",
        )
    ):

        def vector(name: str) -> tuple[float, ...] | None:
            # Absent and present-but-null are the same here: a case that states one draw kind
            # leaves the other two out entirely.
            values = inputs.get(name)
            return tuple(float(v) for v in values) if values is not None else None

        draws = (
            vector("gas_side_draw_fractions"),
            vector("liquid_side_draw_fractions"),
            vector("pumparound_fractions"),
        )

    states = _states(
        [str(name) for name in inputs["components"]],
        float(inputs["feed_n"]),
        [float(v) for v in inputs["feed_z"]],
        float(inputs["feed_t"]),
        float(inputs["feed_p"]),
        int(inputs["number_of_stages"]),
        int(inputs["feed_stage"]),
        bool(inputs["has_reboiler"]),
        bool(inputs["has_condenser"]),
        float(inputs["top_pressure"]),
        float(inputs["bottom_pressure"]),
        float(inputs["reboiler_temperature"]) if "reboiler_temperature" in inputs else None,
        float(inputs["condenser_temperature"]) if "condenser_temperature" in inputs else None,
        float(inputs["temperature_tolerance"]),
        int(inputs["max_iterations"]),
        _spec_of(inputs, "top"),
        _spec_of(inputs, "bottom"),
        str(inputs["solver_type"]) if "solver_type" in inputs else None,
        reactive=section,
        draws=draws,
    )
    layers: dict[str, float] = {}
    for i in range(len(states.tray_temperature)):
        layers[f"tray{i}_temperature_K"] = states.tray_temperature[i]
        layers[f"tray{i}_pressure_bara"] = states.tray_pressure[i] / 1.0e5
        layers[f"tray{i}_gas_n"] = states.tray_gas_n[i]
        layers[f"tray{i}_liquid_n"] = states.tray_liquid_n[i]
    # The three draws the capture states, summed: **its side-draw rows draw on one tray**, so the
    # sum is that tray's draw and the diff holds the mechanism's own number.
    layers["gas_side_draw_n"] = sum(states.gas_side_draw_n)
    layers["liquid_side_draw_n"] = sum(states.liquid_side_draw_n)
    layers["pumparound_n"] = sum(states.pumparound_n)
    return layers


def _packed_column(inputs: Mapping[str, Any]) -> dict[str, float]:
    """`process.packed_column`'s layers: the base column's profile at the derived stage count.

    The same layer as the distillation column's, because the class's arithmetic is: the packing
    reaches the solve only through the stage count a constructor derives from the packed height,
    and the rest of the packing is a hydraulics report on the far side of the converged column.
    """
    from azoth.process.reference.distillation_column import _states
    from azoth.process.reference.packed_column import stage_count

    packed_height = float(inputs["packed_height"])
    states = _states(
        [str(name) for name in inputs["components"]],
        float(inputs["feed_n"]),
        [float(v) for v in inputs["feed_z"]],
        float(inputs["feed_t"]),
        float(inputs["feed_p"]),
        stage_count(packed_height),
        int(inputs["feed_stage"]),
        bool(inputs["has_reboiler"]),
        bool(inputs["has_condenser"]),
        float(inputs["top_pressure"]),
        float(inputs["bottom_pressure"]),
        float(inputs["reboiler_temperature"]) if "reboiler_temperature" in inputs else None,
        float(inputs["condenser_temperature"]) if "condenser_temperature" in inputs else None,
        float(inputs["temperature_tolerance"]),
        int(inputs["max_iterations"]),
        _spec_of(inputs, "top"),
        _spec_of(inputs, "bottom"),
        str(inputs["solver_type"]) if "solver_type" in inputs else None,
    )
    layers: dict[str, float] = {}
    for i in range(len(states.tray_temperature)):
        layers[f"tray{i}_temperature_K"] = states.tray_temperature[i]
        layers[f"tray{i}_pressure_bara"] = states.tray_pressure[i] / 1.0e5
        layers[f"tray{i}_gas_n"] = states.tray_gas_n[i]
        layers[f"tray{i}_liquid_n"] = states.tray_liquid_n[i]
    return layers


def _rate_based_packed_column(inputs: Mapping[str, Any]) -> dict[str, float]:
    """`process.rate_based_packed_column`'s layers: the segment profile itself.

    **A rate-based column's four ports being out says nothing about which segment moved**, and
    the segment profile is the layer - the same reason `_distillation_column` dumps the tray
    profile. The keys are the capture's own, so a divergence names the segment it is in.

    Three mixtures are formed per segment - the two systems and an interface mix - so this is
    the closest layer the model has to its own arithmetic.

    **The two densities are deliberately not dumped.** They are the snapshot's *inputs*, read
    from NeqSim's `getDensity("kg/m3")` - the cubic's volume with the Peneloux shift - and the
    two implementations' shift paths agree to `6e-8` relative rather than to the dump's own
    band. A key that never agrees is noise in the one instrument whose job is to say which
    segment moved; the wetted area and the two film coefficients already carry what the density
    feeds.
    """
    from azoth.process.reference.rate_based_packed_column import DEFAULTS, _states

    def value(name: str, fallback: float) -> float:
        return float(inputs[name]) if name in inputs else fallback

    states = _states(
        [str(name) for name in inputs["gas_components"]],
        float(inputs["gas_n"]),
        [float(v) for v in inputs["gas_z"]],
        float(inputs["gas_p"]),
        float(inputs["gas_t"]),
        [str(name) for name in inputs["liquid_components"]],
        float(inputs["liquid_n"]),
        [float(v) for v in inputs["liquid_z"]],
        float(inputs["liquid_p"]),
        float(inputs["liquid_t"]),
        [str(name) for name in inputs["transfer_components"]]
        if "transfer_components" in inputs
        else None,
        value("column_diameter", DEFAULTS["column_diameter"]),
        value("packed_height", DEFAULTS["packed_height"]),
        int(value("number_of_segments", DEFAULTS["number_of_segments"])),
        str(inputs["packing_type"]) if "packing_type" in inputs else DEFAULTS["packing_type"],
        int(value("max_iterations", DEFAULTS["max_iterations"])),
        value("convergence_tolerance", DEFAULTS["convergence_tolerance"]),
        value("mass_transfer_correction", DEFAULTS["mass_transfer_correction"]),
        str(inputs["mass_transfer_correlation"])
        if "mass_transfer_correlation" in inputs
        else DEFAULTS["mass_transfer_correlation"],
        str(inputs["film_model"]) if "film_model" in inputs else DEFAULTS["film_model"],
        str(inputs["heat_transfer_model"])
        if "heat_transfer_model" in inputs
        else DEFAULTS["heat_transfer_model"],
        str(inputs["segment_solver"]) if "segment_solver" in inputs else DEFAULTS["segment_solver"],
        str(inputs["column_solver"]) if "column_solver" in inputs else DEFAULTS["column_solver"],
    )
    layers: dict[str, float] = {}
    for segment in states["segments"]:
        index = segment["number"]
        layers[f"segment{index}_gas_temperature_K"] = segment["gas_temperature"]
        layers[f"segment{index}_liquid_temperature_K"] = segment["liquid_temperature"]
        layers[f"segment{index}_wetted_area_m2_per_m3"] = segment["wetted_area"]
        layers[f"segment{index}_kga"] = segment["k_ga"]
        layers[f"segment{index}_kla"] = segment["k_la"]
        layers[f"segment{index}_interface_temperature_K"] = segment["interface_temperature"]
        layers[f"segment{index}_heat_transfer_rate_W"] = segment["heat_transfer_rate"]
        layers[f"segment{index}_percent_flood"] = segment["percent_flood"]
        layers[f"segment{index}_net_molar_transfer_mol_per_s"] = segment["net_molar_transfer"]
    return layers


DUMPERS: dict[str, Dumper] = {
    "process.rate_based_packed_column": _rate_based_packed_column,
    "process.absorption_column": _absorption_column,
    "process.packed_column": _packed_column,
    "process.stripping_column": _stripping_column,
    "process.component_splitter": _component_splitter,
    "process.distillation_column": _distillation_column,
    "process.expander": _expander,
    "process.manifold": _manifold,
    "process.pipe": _pipe,
    "process.pump": _pump,
    "process.splitter": _splitter,
    "process.mixer": _mixer,
    "process.compressor": _compressor,
    "process.heat_exchanger": _heat_exchanger,
    "process.shortcut_distillation_column": _shortcut_distillation_column,
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
    "process.gas_scrubber": (
        "the same kernel as ``process.separator`` and therefore the same answer: `GasScrubber` "
        "does not override `run`, so the two entries' captures carry the same rows. A dumper "
        "here would be a second dump of the separator's arithmetic under another name, which "
        "is the thing this dict exists to refuse - and the class's own arithmetic, the "
        "Souders-Brown capacity metric, is mechanical design that no port field reaches"
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
    "process.stirred_tank_reactor": (
        "the kernel forms the reacted composition, the outlet pressure and the flashed state, "
        "and **every one of them is on the record** - the composition is the product's ``z``, the "
        "pressure its ``P``, the state its ``T`` and ``h``. The extent is the one intermediate "
        "the record does not carry, and it is not hidden: the composition *is* the extent, "
        "applied. The duty the isothermal branch reports is a difference of the two records' "
        "enthalpies over the flow, which is an output rather than a layer"
    ),
    "process.flare": (
        "the kernel is a *pass-through*: `run` clones the inlet into the outlet, so the whole "
        "record is the answer and there is no intermediate between the two ports. Its two "
        "reported numbers - the duty and the CO2 emission - are outputs rather than layers, "
        "and the capture prints them beside the record exactly as the model does"
    ),
    "process.ejector": (
        "the kernel forms the mixing pressure, four velocities, two efficiencies' enthalpy "
        "drops and the joined fluid's static enthalpy - and the outlet record carries the "
        "answer to all of them, because the last flash *is* the state the diffuser lands on. "
        "**The one intermediate not on the record is the class's own instrumentation**: the "
        "entrainment, compression and expansion ratios and the area ratio are ratios of the "
        "record's fields (`n`, `P` and `P`), and the Mach numbers need a speed of sound the "
        "port would have to model separately rather than derive. A dumper would be dumping "
        "arithmetic the capture's own ratios already state"
    ),
    "process.three_phase_separator": (
        "the kernel forms the multiphase flash's three fractions, the three compositions and "
        "the moles the entrainment moves, and every one of them is on the three outlet "
        "records - a phase's share is its outlet's flow over the feed's, and its composition "
        "is that outlet's ``z``. The **labels** are the one intermediate the records do not "
        "carry, and a swap is not a hidden layer: the phases are named by "
        "`PhaseEos.init`'s rule and routed to the outlet that names them, so an oil sent to "
        "the water outlet is a wrong ``z`` in the capture rather than a state beside it. And "
        "the flash itself is held by the first row, which states the same feed with no "
        "entrainment - the other five are that answer with material moved, and where it moved "
        "to is the difference between their flows and its"
    ),
    "process.tank": (
        "the same kernel as ``process.separator`` at zero drop, zero entrainment and no "
        "heat input, and the capture's two-phase row is the separator's first row to the "
        "last digit. A dumper here would be a second dump of the separator's arithmetic "
        "under another name, which is the thing this dict exists to refuse. **The capture's "
        "other rows are the reason this is not a case rather than a dump**: two of them are "
        "single-phase, where `run`'s absent-oil branch writes `gasOutStream` instead of the "
        "liquid one and the two outlets come back swapped against `Separator`'s on the same "
        "feed - a construction-time state a kernel that is a pure function of its inlets "
        "has no way to reproduce"
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


#: The component splitter's diverging layers, **per row**, because how far apart the class's
#: enthalpies are from the state's is a property of the state.
#:
#: **`ComponentSplitter.run` reports outlet enthalpies that are not the states it reports.**
#: The first row's bottoms reads `-14431.31737097011` where a *fresh* NeqSim fluid at that
#: outlet's own composition, temperature and pressure gives `-20309.24832914077`, and a fresh
#: fluid is what this port computes. The bounds sit a little under each measured ratio, so
#: what is asserted is that the divergence is still there and at least as wide - not its last
#: digits. The first row's overhead is the narrowest at `1.06`, which is worth seeing: the
#: defect is large where the outlet crosses a phase boundary and small where it does not.
_COMPONENT_SPLITTER_DIVERGENCE: tuple[Divergence, ...] = (
    Divergence(
        key="overhead_h",
        at_least=1.05,
        reason="`run` reports an enthalpy the outlet's own composition does not have; a "
        "fresh NeqSim fluid at that state gives the port's number, and on this row the two "
        "are 6% apart",
    ),
    Divergence(
        key="bottoms_h",
        at_least=1.3,
        reason="the same, and wider here: `-14431.32` against the fresh fluid's `-20309.25`",
    ),
)

#: The even-routing row, where both outlets carry the *feed's* composition and the class's
#: enthalpy is more than twice the state's.
_COMPONENT_SPLITTER_EVEN_DIVERGENCE: tuple[Divergence, ...] = (
    Divergence(
        key="overhead_h",
        at_least=2.0,
        reason="both outlets are the feed, so the state's enthalpy is the feed's `-8984.69` "
        "and the class reports `-18823.74`",
    ),
    Divergence(
        key="bottoms_h",
        at_least=2.0,
        reason="the same number on the other outlet",
    ),
)

#: The all-of-one-component row, where the overhead loses a component entirely and the two
#: enthalpies part by more than forty times.
_COMPONENT_SPLITTER_SEPARATION_DIVERGENCE: tuple[Divergence, ...] = (
    Divergence(
        key="overhead_h",
        at_least=40.0,
        reason="the overhead takes all of the methane and none of the pentane, and the class "
        "reports `-45658.96` where the state is `-1032.85`",
    ),
    Divergence(
        key="bottoms_h",
        at_least=1.3,
        reason="`-16148.11` against the state's `-21393.95`",
    ),
)


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
        model="process.rate_based_packed_column",
        case="co2_water_absorber",
        capture="process_rate_based_packed_column.tsv",
        block=0,
        identified_by=("#label", "co2_water_absorber"),
        # **The heat rate is the model's loosest layer, and it is bounded absolutely.**
        # Measured, the four segments differ from the capture by `1.955e-4` relative - the
        # largest of any dumped key, because the heat step's capacity rates go through the
        # *phase mass* (`n M`) and the class's own `getMass()` carries the flash's own rounding.
        # A band of 1 W on a rate of order 2731 W is `3.7e-4` relative and names the key rather
        # than loosening the case.
        diagnostic=(
            ("segment1_heat_transfer_rate_W", 1.0),
            ("segment2_heat_transfer_rate_W", 1.0),
            ("segment3_heat_transfer_rate_W", 1.0),
            ("segment4_heat_transfer_rate_W", 1.0),
        ),
    ),
    LayerCase(
        model="process.rate_based_packed_column",
        case="heat_transfer_disabled",
        capture="process_rate_based_packed_column.tsv",
        block=6,
        identified_by=("#label", "heat_transfer_disabled"),
    ),
    LayerCase(
        model="process.rate_based_packed_column",
        case="zero_packed_height",
        capture="process_rate_based_packed_column.tsv",
        block=7,
        identified_by=("#label", "zero_packed_height"),
    ),
    # The specification rows. **Four are cases and three are declared uncased**, and the split
    # is the measurement: a purity, a flow rate, a duty and a purity at the *bottom* location
    # converge, while a recovery specification does not converge in NeqSim at all, a reflux
    # ratio lands on a state this port's own flash cannot find, and a duty under a temperature
    # pin is inert.
    LayerCase(
        model="process.distillation_column",
        case="spec_top_purity",
        capture="process_column.tsv",
        block=4,
        identified_by=("#label", "spec_top_purity_0_98_methane"),
    ),
    LayerCase(
        model="process.distillation_column",
        case="spec_top_flow_rate",
        capture="process_column.tsv",
        block=6,
        identified_by=("#label", "spec_top_flow_rate_14000_mol_per_hour"),
    ),
    LayerCase(
        model="process.distillation_column",
        case="spec_top_duty",
        capture="process_column.tsv",
        block=9,
        identified_by=("#label", "spec_top_duty_minus_20000"),
    ),
    # **The one row at the other location**, whose secant drives the reboiler's temperature
    # rather than the condenser's - so it is the row that fails if a location is ever treated
    # as a field rather than the slot a specification sits in.
    LayerCase(
        model="process.distillation_column",
        case="spec_bottom_purity",
        capture="process_column.tsv",
        block=10,
        identified_by=("#label", "spec_bottom_purity_0_98_n_butane"),
    ),
    # **The second solve's row**, on the column capture's own last block: the same binary
    # column under `solver_type = "naphtali_sandholm"`, whose layer is the same profile read
    # off a mesh solve's own state. The ladder capture holds the identical row among the ten
    # strategies - a measurement, not an oracle, so it is not a case.
    LayerCase(
        model="process.distillation_column",
        case="binary_mesh_solve",
        capture="process_column.tsv",
        block=11,
        identified_by=("#label", "binary_methane_butane_mesh_solve"),
    ),
    # The absorber's two rows: the oracle - unpinned at `1e-4` K - and the class's own
    # isothermal case, which pins every tray and therefore stops after one sweep.
    LayerCase(
        model="process.absorption_column",
        case="lean_oil_absorber",
        capture="process_absorber.tsv",
        block=2,
        identified_by=("#label", "lean_oil_absorber"),
    ),
    LayerCase(
        model="process.stripping_column",
        case="hydrocarbon_stripper",
        capture="process_absorber.tsv",
        block=5,
        identified_by=("#label", "hydrocarbon_stripper"),
    ),
    # The column's converged rows. **Both are cases, and the second is there because a looser
    # gate is a different measurement** rather than a sloppier version of the first: where a
    # solve stops is what its answer is.
    LayerCase(
        model="process.distillation_column",
        case="binary_rigorous",
        capture="process_column.tsv",
        block=2,
        identified_by=("#label", "binary_methane_butane_4_stages"),
    ),
    LayerCase(
        model="process.distillation_column",
        case="binary_looser_gate",
        capture="process_column.tsv",
        block=3,
        identified_by=("#label", "binary_methane_butane_loose"),
    ),
    # The component splitter's three rows. The two enthalpies are declared divergences - see
    # `_COMPONENT_SPLITTER_DIVERGENCE` - and the flows, compositions and temperatures are
    # compared, which is what the class gets right.
    LayerCase(
        model="process.component_splitter",
        case="near_total_separation",
        capture="process_component_splitter.tsv",
        block=0,
        identified_by=("split_factors", "0.98 0.05 0.02"),
        divergence=_COMPONENT_SPLITTER_DIVERGENCE,
    ),
    LayerCase(
        model="process.component_splitter",
        case="an_even_routing_is_the_feed_twice",
        capture="process_component_splitter.tsv",
        block=1,
        identified_by=("split_factors", "0.5 0.5 0.5"),
        divergence=_COMPONENT_SPLITTER_EVEN_DIVERGENCE,
    ),
    LayerCase(
        model="process.component_splitter",
        case="all_of_one_component",
        capture="process_component_splitter.tsv",
        block=2,
        identified_by=("split_factors", "1.0 0.5 0.0"),
        divergence=_COMPONENT_SPLITTER_SEPARATION_DIVERGENCE,
    ),
    # The shortcut column's four rows. **Their block order is the probe's**, and the two
    # degenerate rows are the gap between the fourth and the fifth: `UNCASED_ROWS` declares
    # them, and this tuple names the other four by position.
    LayerCase(
        model="process.shortcut_distillation_column",
        case="propane_nbutane_split",
        capture="process_shortcut_distillation_column.tsv",
        block=0,
        identified_by=("#label", "propane_nbutane_300K_20bara"),
    ),
    LayerCase(
        model="process.shortcut_distillation_column",
        case="pressures_stated_move_the_product_phase",
        capture="process_shortcut_distillation_column.tsv",
        block=1,
        identified_by=("#label", "pressures_stated_18_21_bara"),
    ),
    LayerCase(
        model="process.shortcut_distillation_column",
        case="wilson_fallback_superheated_gas",
        capture="process_shortcut_distillation_column.tsv",
        block=3,
        identified_by=("#label", "wilson_fallback_superheated_gas"),
    ),
    LayerCase(
        model="process.shortcut_distillation_column",
        case="binary_methane_nbutane",
        capture="process_shortcut_distillation_column.tsv",
        block=4,
        identified_by=("#label", "binary_methane_nbutane"),
    ),
    # The isentropic pair. **Both rows of each are cases, and the efficiency-of-one row is
    # the one that matters**: its reversible step *is* its answer, so the capture states the
    # isentropic enthalpy the class keeps in a local - and its entropy production is zero by
    # construction, which is why it is compared on an absolute bound.
    LayerCase(
        model="process.compressor",
        case="the_reversible_limit_at_efficiency_one",
        capture="process_compressor.tsv",
        block=1,
        identified_by=("isentropic_efficiency", "1.0"),
        # The pump's bound and the pump's reason: 1e-9 kJ/(mol*K) is 1e-6 J/(mol*K), four
        # orders below any physical irreversibility and four above the double-precision
        # rounding of an enthalpy near 2.5e4. **Measured on these rows**: the compressor's two
        # sides are `-1.4e-17` and `1.4e-15`, and the expander's `-4.8e-12` and `1.4e-12` -
        # both zero, and neither of them equal.
        diagnostic=(("entropy_production_kJ_per_molK", 1e-9),),
    ),
    LayerCase(
        model="process.compressor",
        case="the_efficiency_divides_the_step",
        capture="process_compressor.tsv",
        block=0,
        identified_by=("isentropic_efficiency", "0.75"),
    ),
    LayerCase(
        model="process.expander",
        case="the_reversible_limit_at_efficiency_one",
        capture="process_expander.tsv",
        block=1,
        identified_by=("isentropic_efficiency", "1.0"),
        # The pump's bound and the pump's reason: 1e-9 kJ/(mol*K) is 1e-6 J/(mol*K), four
        # orders below any physical irreversibility and four above the double-precision
        # rounding of an enthalpy near 2.5e4. **Measured on these rows**: the compressor's two
        # sides are `-1.4e-17` and `1.4e-15`, and the expander's `-4.8e-12` and `1.4e-12` -
        # both zero, and neither of them equal.
        diagnostic=(("entropy_production_kJ_per_molK", 1e-9),),
    ),
    LayerCase(
        model="process.expander",
        case="the_efficiency_multiplies_the_step_and_not_divides_it",
        capture="process_expander.tsv",
        block=0,
        identified_by=("isentropic_efficiency", "0.75"),
    ),
    # The pipe's three rows, one per branch and one that diverges. The water row is a case
    # *because the case records the divergence*, not because the model reproduces it: the
    # layers that moved there are the Reynolds number and the friction factor, and this is
    # where a reader sees by how much.
    # The manifold's three rows: two outlets, three, and one feed at zero flow. It is the
    # first id whose capture prints *both* of its composed models' rows - the mixture and
    # the branches - so the dumper reports both and the comparison is the composition.
    LayerCase(
        model="process.manifold",
        case="two_feeds_two_outlets",
        capture="process_manifold.tsv",
        block=0,
        identified_by=("split_factors", "0.25 0.75"),
    ),
    LayerCase(
        model="process.manifold",
        case="two_feeds_three_outlets",
        capture="process_manifold.tsv",
        block=1,
        identified_by=("split_factors", "0.2 0.3 0.5"),
    ),
    LayerCase(
        model="process.manifold",
        case="a_zero_flow_feed_is_dropped",
        capture="process_manifold.tsv",
        block=2,
        identified_by=("split_factors", "0.5 0.5"),
    ),
    LayerCase(
        model="process.pipe",
        case="gas_methane_co2_1000m",
        capture="process_pipe.tsv",
        block=0,
        identified_by=(LABEL_KEY, "gas_methane_co2_1000m"),
    ),
    LayerCase(
        model="process.pipe",
        case="liquid_n_butane_1000m",
        capture="process_pipe.tsv",
        block=1,
        identified_by=(LABEL_KEY, "liquid_n_butane_1000m"),
    ),
    # **The water row, which was three declared divergences and is now a case like the
    # others.** NeqSim gives an aqueous phase `WaterPhysicalProperties` and the liquid
    # `Viscosity` correlation, which `eos.aqueous_viscosity` ports; before that id existed
    # this port used the PFCT form its gas and oil branches take and was `1.61` out on the
    # Reynolds number and the drop. The row is kept because it is the measurement of that.
    LayerCase(
        model="process.pipe",
        case="water_1000m",
        capture="process_pipe.tsv",
        block=2,
        identified_by=(LABEL_KEY, "liquid_water_1000m"),
    ),
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
    # **The reactive section's block**, on the binary state the model's own case is held to:
    # NeqSim's reactive column there is bit-identical to its standard twin, which is the class's
    # own `NR = 0` claim. It sits in this capture because the layer diff pairs one capture with
    # one model, and this model's case is two blocks above.
    LayerCase(
        model="process.distillation_column",
        case="reactive_section_is_the_delegations_own",
        capture="process_column.tsv",
        block=13,
        identified_by=("#label", "binary_reactive_column_pr"),
    ),
    # **The side-draw block**, four below the reactive one and the last state in this capture: the
    # binary column with tray 3 drawing a quarter of its vapour. The three blocks above it are the
    # class's own one-tray columns, which this port refuses, and the block below it is the same
    # column drawing liquid and a pumparound, which NeqSim reconciles - see `UNCASED_ROWS`.
    LayerCase(
        model="process.distillation_column",
        case="a_tray_withdraws_a_quarter_of_its_vapour",
        capture="process_column.tsv",
        block=17,
        identified_by=("#label", "side_draw_column_binary_gas_quarter"),
    ),
    # The packed column's two cases sit on blocks 12 and 13. **Its first block is the base
    # column's own `binary_rigorous` state**, because the packing does not change the separation:
    # NeqSim's `PackedColumn` at `2.0` m and its `DistillationColumn` at four stages report
    # bit-identical answers on all twenty-two captured quantities. The second is the same feed at
    # `2.3` m, where `ceil(4.6)` makes it a five-stage column instead - so it is the row that
    # fails if the height stops moving the stage count. See `UNCASED_ROWS` for the other thirteen.
    LayerCase(
        model="process.packed_column",
        case="packed_binary_2m",
        capture="process_packed_column.tsv",
        block=12,
        identified_by=(LABEL_KEY, "packed_distillation_binary_2m"),
    ),
    LayerCase(
        model="process.packed_column",
        case="packed_binary_2m3_is_one_stage_more",
        capture="process_packed_column.tsv",
        block=13,
        identified_by=(LABEL_KEY, "packed_distillation_binary_2m3"),
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
    # The rate-based column's seven. `_srk` is the class's own state on the cubic its test file
    # uses, kept beside the PR row that the port is held to so the pair measures what the cubic
    # moved; the two TEG rows are `SystemSrkCPAstatoil` with mixing rule 10, which PR cannot
    # represent - so they are evidence for a refusal rather than an oracle for a state - and the
    # three remaining CO2/water rows are the stripper and the two heights, which the port's own
    # test asserts as *directions* rather than as digits.
    "process_rate_based_packed_column.tsv": 7,
    "process_heat_exchanger.tsv": 2,
    # The shortcut column's two degenerate rows. **Both are evidence rather than oracles**:
    # a reflux multiplier of exactly one leaves Gilliland's `X` at zero and the class returns
    # `actual_stages = Infinity` with a feed tray of zero - `(int) Math.round(Infinity) + 1`
    # wrapping through `Integer.MIN_VALUE` - and swapped keys give the class's own
    # `solved = false` with every answer at its field initialiser. The port **refuses** both,
    # so neither has a state to be its case; they stay in the capture as the measurement the
    # refusal rests on.
    "process_shortcut_distillation_column.tsv": 2,
    # The column's two deethanizer rows. **Both are the class's own failure to converge** - the
    # same column at a hard-capped 80 iterations and at 200, ending `7.3e-3` and `2.9e-2` K from
    # a `1e-5` gate, the second worse than the first. The port converges on the same column, so
    # its state is not NeqSim's and the rows are evidence rather than oracles: what the port's
    # own test asserts about them is self-consistency, and this is why they are declared here
    # rather than paired with a case.
    # The column's own, one more: the reactive block's *standard* twin, which is the same state
    # at the same gate as the case two blocks above and is kept as the measurement the reactive
    # row is identical to.
    # And the side-draw capture's four: the class's own three one-tray columns, and the binary
    # column drawing liquid and a pumparound. **The three are the class's own states for this
    # mechanism** - `columnReportsSideDrawAsOutletStream` and `columnEnergyBalanceIncludes-
    # SideDrawStreams` run a one-tray column with no ends and a pure feed, and the class is happy
    # with the zero-flow product that leaves: the methane rows report a bottoms of `n = 0` with an
    # enthalpy of `-Infinity` (`getMaterialOutletEnthalpy` divides an empty system's enthalpy by
    # its zero moles) and a mass balance of exactly `0.0`. **The port refuses all three by name**,
    # because its products must exist: the two methane rows have no liquid for the bottom to
    # publish, and the pentane row's one tray is all liquid, so the column's two products would be
    # the same stream and the closure counts it twice - refused as `SolverNotConverged` with a
    # residual of `0.75`.
    # The fourth is a divergence rather than a refusal. **NeqSim does not close its own balance
    # there**: it answers `RECONCILED_PRODUCTS` with a mass balance of `-2.449` kg/hr against a
    # 1000 kg/hr feed and an energy error of `6.4e-2`, and its two products come out `0.0200`
    # mol/s larger than the port's. The port closes to `9.7e-9` and treats the drawn liquid as a
    # net outlet, which is what `SimpleTray`'s own split does - a pumparound fraction *without* a
    # return. What the class's reconciliation restores is the pumparound's return, which is
    # `ColumnPumparound`'s subject and is named out of this tranche.
    "process_column.tsv": 10,
    # One capture for two ids, because the two machines it drives are one class with two names,
    # and five of its six rows are uncased for each of them. **The pinned pair is the classes'
    # own isothermal case**: `setOutletTemperature` on every stage makes the base's gate exactly
    # zero, so the solve stops after its first sweep - and the port **refuses** the state it
    # lands on, because that state's energy closure is `0.28` against the class's own default
    # gate of `1.6e-2`. The classes' *tests* accept it only because they loosened that gate to
    # `5e-2`, a settable this port has not carried. The two loose-gate rows are the same
    # measurement one step out, and the stripper's pinned row states `MESH_RESIDUAL` besides.
    "process_absorber.tsv": 5,
    # The packed column's thirteen: eleven stage-count rows and two solved rows. **The eleven are
    # the constructor's own arithmetic** - `PackedColumn(name, height, packing, true, true)` with
    # no feed and no solve, printing the trays it made - so there is no model state for them to be
    # a case of; the port's own test asserts the rule against all eleven. The two solved rows are
    # evidence rather than oracles: **the class's own test case converges to a state that is not a
    # column**, reporting `RIGOROUS_CONVERGED` with residuals of exactly zero while the trays
    # below the feed carry no traffic at all, and the contactor row reports
    # `FALLBACK_PRODUCTS`, which the class's own warning says means the products are a single
    # flash of the mixed feeds rather than the tray solution.
    "process_packed_column.tsv": 13,
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
            # **The distance between the two, in whichever direction it runs.** A divergence
            # whose sign flips from row to row cannot be stated as `capture / azoth >= n`,
            # because half its rows sit below one - `process.component_splitter`'s outlet
            # enthalpies are `0.71` and `0.75` on two rows and `2.09` on the third, the class
            # reporting an enthalpy below the state's on one state and above on another.
            # `max(r, 1/r)` says what the older rule said wherever the capture was the larger
            # - `process.splitter`'s three ratios are all above one, so its bounds read
            # identically - and says something on the rows where it was not.
            ratio = abs(theirs) / abs(mine) if mine else float("inf")
            ratio = max(ratio, 1.0 / ratio) if ratio > 0.0 else float("inf")
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
