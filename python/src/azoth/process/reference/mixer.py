"""``process.mixer`` - several feeds blended into one.

Spec: ``specs/models/process/mixer.yaml``, which carries this model's provenance and the
two places its arithmetic diverges from NeqSim's.

The outlet is at the **lowest** inlet pressure - a mixer is a vessel, and nothing in it
can be above the pressure any feed arrives at - and its temperature is an isenthalpic
flash of the flow-weighted blend. There is no ``stream`` type: the six things a stream is
are the arguments a unit operation takes, and a type would be a second way to state them.
"""

from __future__ import annotations

from azoth import _models_gen
from azoth.core.errors import InvalidInputError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import MixerResult
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning
from azoth.eos.mixture import Mixture
from azoth.eos.reference.molar_enthalpy_entropy import IdealGasModel
from azoth.eos.reference.ph_flash import enthalpy_at, ph_flash

MODEL_ID = "process.mixer"


def mixer(
    mixture: Mixture,
    ideal_gas: IdealGasModel,
    T: list[Q],
    P: list[Q],
    n: list[Q],
    z: list[list[float]],
) -> MixerResult:
    """Several feeds blended into one, at the lowest inlet pressure.

    Args:
        mixture: the components and their interaction parameters, **one mixture for every
            inlet**. All the feeds share a component set here, which is what makes the
            accumulation a vector sum.
        ideal_gas: the datum the enthalpy balance is measured from.
        T: each inlet's absolute temperature.
        P: each inlet's absolute pressure. The outlet is the minimum.
        n: each inlet's molar flow rate.
        z: each inlet's mole fractions, one row per inlet.

    Returns:
        The blend's state, its flow and its composition.

    Raises:
        InvalidInputError: if the inlet vectors disagree about how many streams there
            are, if ``z`` does not carry one composition per inlet, or if any inlet's
            composition is not a composition.
        OutOfRangeError: if an inlet state is non-positive, or an inlet carries no flow.
        SolverNotConvergedError: if the blend's enthalpy is outside the flash's bracket.
    """
    spec = _models_gen.model(MODEL_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    t_si = [input_to_si(spec, "T", value) for value in T]
    p_si = [input_to_si(spec, "P", value) for value in P]
    n_si = [input_to_si(spec, "n", value) for value in n]

    streams = len(t_si)
    components = len(mixture.components)
    if len(p_si) != streams or len(n_si) != streams:
        raise InvalidInputError(
            "P",
            f"a mixer's inlets must be the same number of streams; got {streams} "
            f"temperature(s), {len(p_si)} pressure(s) and {len(n_si)} flow(s)",
        )
    rows = [list(row) for row in z]
    if len(rows) != streams or any(len(row) != components for row in rows):
        raise InvalidInputError(
            "z",
            f"`z` must carry one composition per inlet - {streams} rows of {components} - "
            f"and carries {len(rows)} row(s) of {[len(r) for r in rows]}",
        )

    # A vector input has no single value for a range check, so the binding one is
    # reported: the smallest inlet flow, the lowest inlet pressure and the lowest inlet
    # temperature are the ones that would produce a meaningless blend, which is what the
    # checks are for. All three of the spec's bounds are resolved here, and a fourth
    # declared bound resolving to nothing would report itself as unevaluated rather than
    # pass.
    apply_checks(
        checks.on_input,
        {"n": min(n_si), "P": min(p_si), "T": min(t_si)}.get,
        warnings,
    )

    total_flow = sum(n_si)
    outlet_pressure = min(p_si)

    z_out = [0.0] * components
    total_enthalpy = 0.0
    for stream in range(streams):
        enthalpy, _ = enthalpy_at(mixture, ideal_gas, t_si[stream], p_si[stream], rows[stream])
        total_enthalpy += n_si[stream] * enthalpy
        for index in range(components):
            z_out[index] += n_si[stream] * rows[stream][index]

    blended: list[float] = [value / total_flow for value in z_out]
    molar_enthalpy = total_enthalpy / total_flow

    flash = ph_flash(
        mixture,
        ideal_gas,
        from_si(outlet_pressure, "Pa"),
        from_si(molar_enthalpy, "J/mol"),
        blended,
    )
    warnings.extend(flash.warnings)

    return MixerResult(
        T=flash.T,
        P=from_si(outlet_pressure, "Pa"),
        flow=total_flow,
        z_out=tuple(blended),
        beta=flash.beta,
        phase=flash.phase,
        iterations=flash.iterations,
        warnings=tuple(warnings),
    )
