"""``process.compressor`` - a pressure rise at a stated isentropic efficiency.

Spec: ``specs/models/process/compressor.yaml``, which carries the provenance and what is
not ported. The procedure is :mod:`azoth.process.reference._isentropic`'s - the three
machines share it and differ only in which way the efficiency scales the ideal enthalpy
change. This module is the argument handling and the result, and the same is true of the
pump and the expander.
"""

from __future__ import annotations

from azoth import _models_gen
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import CompressorResult
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning
from azoth.eos.mixture import Mixture
from azoth.eos.reference.molar_enthalpy_entropy import IdealGasModel
from azoth.process.reference._isentropic import Direction, isentropic_run

MODEL_ID = "process.compressor"


def compressor(
    mixture: Mixture,
    ideal_gas: IdealGasModel,
    T: Q,
    P: Q,
    n: Q,
    z: list[float],
    outlet_pressure: Q,
    efficiency: float,
) -> CompressorResult:
    """A stream compressed to a stated outlet pressure at a stated isentropic efficiency.

    Args:
        mixture: the components and their interaction parameters.
        ideal_gas: the datum the enthalpy and the entropy are measured from.
        T: the inlet absolute temperature.
        P: the inlet absolute pressure. The outlet must be **above** it.
        n: the inlet molar flow rate, which turns the molar enthalpy change into a power.
        z: the inlet mole fractions.
        outlet_pressure: the pressure the compressor delivers.
        efficiency: the isentropic efficiency, in ``(0, 1]``. **Required rather than
            defaulted**, and range-checked: an efficiency of one is a machine with no
            losses, and a compressor assumed ideal is one that understates every duty it
            is asked for.

    Returns:
        The outlet temperature and phase, the shaft power, and the temperature a perfect
        machine would have reached.

    Raises:
        InvalidInputError: if ``outlet_pressure`` is not above the inlet pressure.
        OutOfRangeError: if an input is outside the spec's declared range, including an
            efficiency above one.
        SolverNotConvergedError: if either flash cannot reach its target.
    """
    spec = _models_gen.model(MODEL_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    t_si = input_to_si(spec, "T", T)
    p_si = input_to_si(spec, "P", P)
    n_si = input_to_si(spec, "n", n)
    p_out_si = input_to_si(spec, "outlet_pressure", outlet_pressure)
    apply_checks(
        checks.on_input,
        {
            "T": t_si,
            "P": p_si,
            "n": n_si,
            "outlet_pressure": p_out_si,
            "efficiency": efficiency,
        }.get,
        warnings,
    )

    solved = isentropic_run(
        mixture, ideal_gas, t_si, p_si, p_out_si, n_si, list(z), efficiency, Direction.CONSUMING
    )
    warnings.extend(solved.warnings)

    return CompressorResult(
        T=solved.T,
        P=from_si(p_out_si, "Pa"),
        power=solved.power,
        beta=solved.beta,
        phase=solved.phase,
        isentropic_temperature=solved.isentropic_temperature,
        iterations=solved.iterations,
        warnings=tuple(warnings),
    )
