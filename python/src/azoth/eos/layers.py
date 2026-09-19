"""The intermediates a model's reference kernel already computes, as named rows.

    from azoth.eos import layers
    layers.rows("eos.saft_vr_mie_phase", {"components": ["methane"], "T": 300.0, ...})

**A case compares a total, and a total that is 4% out says nothing about where.** A probe
capture carries the oracle's intermediates - the segment diameters, the packing fraction,
the chain contact value, the three dispersion terms, the mixture parameters - and those
already exist here, inside the reference kernel, on the way to the answer. This exposes
them by name so that `tools/neqsim_layer_diff.py` can put the two side by side and say
which layer moved first.

# The key names are the capture's

Every layer is spelled the way the probe prints it, so a capture and a dump meet without a
mapping table. That is deliberate: a mapping table is where a wrong-layer comparison hides,
and the probes' own headers already declare the correspondence key by key for the same
reason.

The consequence is that a layer set is not a model's whole state - a probe prints what its
author was looking for. A key the dump does not carry is simply not compared, and a model
whose dump shares nothing with a capture is an error rather than a pass.

# What counts as a layer

Each entry is one of three things, and nothing else:

* a field of the reference kernel's own solved state - `a_hs`, `g_hs`, `i1`, `alpha_mix`;
* the value of one of its own functions - `f()`, `pressure_over_rt()`;
* its own decomposition of the Helmholtz energy, applied to those fields, where the kernel
  computes the pieces inside a function rather than exposing them. PC-SAFT's `f()` spells
  out its hard-chain and two dispersion brackets and keeps none of them, and SAFT-VR-Mie's
  are `MieState::f_hc` and `f_disp` in the Rust kernel.

Nothing here re-derives a model. A layer computed by a route the kernel does not take would
be a second implementation, and a divergence in it would be a disagreement with ourselves.

# No model gains an output

This is an introspection surface and not API: nothing is added to a spec's `[outputs]`, and
a caller who wants an answer calls the model.
"""

from __future__ import annotations

import math
from collections.abc import Callable, Mapping
from typing import Any

#: A model's inputs, keyed as a case states them, to its layers. The inputs are bare
#: numbers in the spec's units - the same ones `validation/eos/*.json` carries - so a
#: caller states a state the way every other test in the suite does.
Dumper = Callable[[Mapping[str, Any]], dict[str, float]]


#: NeqSim's internal scale for `a` and `b`, `Component.java:526-531`.
_NEQSIM_INTERNAL = 1.0e5


def _side(inputs: Mapping[str, Any]) -> bool:
    """Whether a case's `compressed_phase` names the liquid root."""
    return str(inputs["compressed_phase"]) == "liquid"


def _dimensional(a: float, b: float, t: float, p: float) -> tuple[float, float]:
    """`(a, b)` in SI, from the reduced pair a cubic is solved in.

    **`ReducedParameters.a` and `.b` are the cubic's dimensionless `A_i` and `B_i`**,
    `a_i P/(R T)^2` and `b_i P/(R T)` - the docstring says so, they are functions of the
    temperature and pressure alone, and the reduced pair is what makes one set of them
    serve both roots. A probe prints SI, because that is what a reader can check against a
    published parameter, so the dump divides the reduction back out.
    """
    rt = 8.3144621 * t
    return a * rt * rt / p, b * rt / p


def _internal(a: float, b: float, t: float, p: float) -> tuple[float, float]:
    """The same pair in NeqSim's own scale, which is SI times `1e5`.

    `UmrCpaProbe` prints `phase.getA()` and `phase.getB()` raw where `CpaSweep` divides
    both by `1e5`, so the same model's mixture parameters are two different numbers in two
    captures. A layer carries the scale its own capture does, because the alternative is a
    conversion table between captures - which is the thing this module exists not to have.
    """
    si_a, si_b = _dimensional(a, b, t, p)
    return si_a * _NEQSIM_INTERNAL, si_b * _NEQSIM_INTERNAL


def _srk_cpa(inputs: Mapping[str, Any]) -> dict[str, float]:
    """The CPA family's layers, under the names `validation/neqsim/CpaSweep.java` prints.

    **Only the cubic and the mixture parameters, not the association's own.** The site
    fractions and the contact value live in `_association`, a kernel of its own and a
    different reach; `CpaSweep` prints them and nothing here compares them.
    """
    from azoth.eos import components as databank
    from azoth.eos.reference._mixture_state import phase_state, reduced_parameters

    mixture = databank.from_names(list(inputs["components"]), eos="srk", associating=True)
    t_si = float(inputs["T"])
    p_si = float(inputs["P"])
    reduced = reduced_parameters(mixture, t_si, p_si)
    state = phase_state(reduced, mixture.kij, list(inputs["z"]), liquid=_side(inputs))

    t_si = float(inputs["T"])
    p_si = float(inputs["P"])
    out: dict[str, float] = {}
    for i, value in enumerate(reduced.a):
        out[f"aT[{i}]"] = _dimensional(value, 0.0, t_si, p_si)[0]
    for i, value in enumerate(reduced.b):
        out[f"b[{i}]"] = _dimensional(0.0, value, t_si, p_si)[1]
    out["A_mix"] = _dimensional(state.a_mix, 0.0, t_si, p_si)[0]
    out["B_mix"] = _dimensional(0.0, state.b_mix, t_si, p_si)[1]
    out["Z"] = state.z
    for i, value in enumerate(state.ln_phi):
        out[f"lnPhi[{i}]"] = value
    return out


def _umr_cpa(inputs: Mapping[str, Any]) -> dict[str, float]:
    """The UMR-CPA layers, under the names `validation/neqsim/UmrCpaProbe.java` prints.

    **`alpha_mix` is the one this exists for.** The universal rule's mixing is where this
    model's new arithmetic is, and a divergence in it otherwise shows up in `Z` as a number
    rather than as a rule.
    """
    from azoth.eos.components import umr_cpa_mixture_of
    from azoth.eos.reference._mixture_state import (
        _umr_alpha_mix,
        phase_state,
        reduced_parameters,
    )
    from azoth.eos.reference.umr_cpa_phase import R

    t_si = float(inputs["T"])
    fluid, _ = umr_cpa_mixture_of(list(inputs["components"]))
    p_si = float(inputs["P"])
    reduced = reduced_parameters(fluid, t_si, p_si)
    x = [float(value) for value in inputs["z"]]
    state = phase_state(reduced, fluid.kij, x, liquid=_side(inputs))

    out: dict[str, float] = {}
    alpha_mix = _umr_alpha_mix(reduced, x)
    if alpha_mix is not None:
        out["alpha_mix"] = alpha_mix
    out["A_attraction"] = _internal(state.a_mix, 0.0, t_si, p_si)[0]
    out["B_covolume"] = _internal(0.0, state.b_mix, t_si, p_si)[1]
    out["Z"] = state.z
    for i, value in enumerate(state.ln_phi):
        out[f"lnPhi[{i}]"] = value
    out["HresTP_J_per_mol"] = state.h_dep_rt * R * t_si
    out["SresTP_J_per_molK"] = state.s_dep_r * R
    return out


def _pcsaft_kij(names: list[str]) -> list[float]:
    """The interaction matrix, as `pcsaft_rahmat_phase`'s own entry function builds it.

    Repeated here rather than lifted out of it because the matrix is an *adapter* - the
    databank gives the upper triangle and the dispersion sums run over ordered pairs - and
    the entry function returns the answer rather than the state it built on the way.
    """
    from azoth.eos.components import pcsaft_kij_for

    pairs = pcsaft_kij_for(tuple(names))
    n = len(names)
    kij = [0.0] * (n * n)
    for (i, j), value in pairs.items():
        kij[i * n + j] = value
        kij[j * n + i] = value
    return kij


def _pcsaft(inputs: Mapping[str, Any]) -> dict[str, float]:
    """The PC-SAFT layers, under the names `validation/neqsim/PcsaftProbe.java` prints.

    **`v` is not offered**, and that is the probe's defect rather than a choice: its `v`
    line is `PhasePCSAFT.getMolarVolume()`, which returns `Infinity` at every state tried,
    so a dump carrying a `v` layer would differ from the capture by an infinity everywhere.
    `volumeSAFT` is the molar volume the class's own solve converged to, and is compared.
    """
    from azoth.eos.reference import pcsaft_rahmat_phase as ref

    names = [str(name) for name in inputs["components"]]
    components = [ref._Component.of(name) for name in names]
    x = [float(value) for value in inputs["z"]]
    t_si = float(inputs["T"])
    p_si = float(inputs["P"])
    side = "liquid" if _side(inputs) else "vapour"
    kij = _pcsaft_kij(names)
    v, _z, _iterations = ref._molar_volume(components, kij, x, t_si, p_si, side)
    state = ref._State(components, kij, x, t_si, v)

    out: dict[str, float] = {
        "m": state.m_bar,
        "mmin1": state.m_minus_1,
        "md": state.md3,
        "volumeSAFT": state.v,
        "nSAFT": state.eta,
        "aHS": state.a_hs,
        "gHS": state.g_hs,
        "f1sum": state.s1,
        "f2sum": state.s2,
        "I1": state.i1,
        "I2": state.i2,
        "C1": state.c1,
        # `f()`'s own three terms, which it computes and keeps none of.
        "F_hc": state.m_bar * state.a_hs - state.m_minus_1 * math.log(state.g_hs),
        "F_disp1": -2.0 * math.pi * state.rho * state.s1 * state.i1,
        "F_disp2": -math.pi * state.m_bar * state.rho * state.s2 * state.i2 * state.c1,
        "F": state.f(),
        "Z": state.pressure_over_rt() * state.v,
    }
    for i, value in enumerate(state.ln_fugacity_coefficients()):
        out[f"lnPhi[{i}]"] = value
    return out


def _saft_vr_mie(inputs: Mapping[str, Any]) -> dict[str, float]:
    """The SAFT-VR-Mie layers, under the names `SaftVrMieProbe.java` prints.

    **`v` is not offered**, for the same reason as PC-SAFT's: that probe prints
    `getVolume()/n` in the class's internal scale, so a dump carrying it would be comparing
    a molar volume against a number 1e5 times it.
    """
    from azoth.eos.reference import saft_vr_mie_phase as ref

    components = [ref._Component.of(name) for name in inputs["components"]]
    x = [float(value) for value in inputs["z"]]
    t_si = float(inputs["T"])
    p_si = float(inputs["P"])
    side = "liquid" if _side(inputs) else "vapour"
    v, _z, _iterations = ref._molar_volume(components, x, t_si, p_si, side)
    state = ref._state(components, x, t_si, v)

    out: dict[str, float] = {
        "volumeSAFT": v,
        "nSAFT": state.eta,
        "aHS": state.a_hs,
        "gHS": state.g_hs,
        "A1Disp": state.a1,
        "A2Disp": state.a2,
        "A3Disp": state.a3,
        # `MieState::f_hc` and `f_disp`, which the reference computes inside `_energy`.
        "F_hc": state.m_bar * state.a_hs - state.m_minus_1 * math.log(state.g_hs),
        "F_disp": state.m_bar * (state.a1 + state.a2 + state.a3),
        "F": ref._energy(state),
        "Z": ref._pressure_over_rt(components, x, t_si, v) * v,
    }
    for i, value in enumerate(ref._ln_phi(components, x, t_si, v)):
        out[f"lnPhi[{i}]"] = value
    return out


#: The models this can dump, which are the ones with a capture whose keys they can meet.
#:
#: `eos.pr_cpa_phase` is absent because NeqSim has no PR-CPA answer to compare it with -
#: its `SystemPrCPA` carries the SRK `Omega` constants and the `aCPA_SRK` fit rather than a
#: Peng-Robinson's - and `eos.tp_flash_saft` because its probe prints prose rather than
#: rows a capture reader can key.
DUMPERS: dict[str, Dumper] = {
    "eos.srk_cpa_phase": _srk_cpa,
    "eos.pcsaft_rahmat_phase": _pcsaft,
    "eos.saft_vr_mie_phase": _saft_vr_mie,
    "eos.umr_cpa_phase": _umr_cpa,
}


def layers(model_id: str, inputs: Mapping[str, Any]) -> dict[str, float]:
    """One model's layers at one state, keyed as its probe prints them.

    Args:
        model_id: a model with a dumper - see [`DUMPERS`].
        inputs: the state, keyed as a case states it.

    Raises:
        KeyError: if the model has no dumper. **A refusal rather than an empty mapping**,
            because a caller comparing an empty mapping against a capture would be told
            that nothing diverged.
    """
    return DUMPERS[model_id](inputs)


def rows(model_id: str, inputs: Mapping[str, Any]) -> list[tuple[str, float]]:
    """The same, as rows in the order the kernel builds them.

    The order is what makes "the first layer that diverges" mean something, and it is the
    order the probes print their own rows in.
    """
    return list(layers(model_id, inputs).items())


__all__ = ["DUMPERS", "layers", "rows"]
