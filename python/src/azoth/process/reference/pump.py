"""``process.pump`` - a pressure rise in a liquid, at a stated isentropic efficiency.

Spec: ``specs/models/process/pump.yaml``, which carries the provenance, which of the four
paths through NeqSim's ``Pump.run`` this is, and why a liquid pump runs an isentropic
flash at all. The procedure is :mod:`azoth.process.reference._isentropic`'s.
"""

from __future__ import annotations

from azoth import _models_gen
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import PumpResult
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning
from azoth.eos.mixture import Mixture
from azoth.eos.reference.molar_enthalpy_entropy import IdealGasModel
from azoth.process.reference._isentropic import Direction, isentropic_run

MODEL_ID = "process.pump"


def pump(
    mixture: Mixture,
    ideal_gas: IdealGasModel,
    T: Q,
    P: Q,
    n: Q,
    z: list[float],
    outlet_pressure: Q,
    efficiency: float,
) -> PumpResult:
    """A stream pumped to a stated outlet pressure at a stated isentropic efficiency.

    Args:
        mixture: the components and their interaction parameters.
        ideal_gas: the datum the enthalpy and the entropy are measured from.
        T: the inlet absolute temperature.
        P: the inlet absolute pressure. The outlet must be **above** it.
        n: the inlet molar flow rate.
        z: the inlet mole fractions.
        outlet_pressure: the pressure the pump delivers.
        efficiency: the isentropic efficiency, in ``(0, 1]``. Required rather than
            defaulted, for the reason :mod:`azoth.process.reference.compressor` gives: a
            silently ideal pump is a pump that understates its power.

    Returns:
        The outlet temperature and phase, the shaft power, and the temperature a perfect
        machine would have reached.

    Raises:
        InvalidInputError: if ``outlet_pressure`` is not above the inlet pressure.
        OutOfRangeError: if an input is outside the spec's declared range.
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

    return PumpResult(
        T=solved.T,
        P=from_si(p_out_si, "Pa"),
        power=solved.power,
        beta=solved.beta,
        phase=solved.phase,
        isentropic_temperature=solved.isentropic_temperature,
        iterations=solved.iterations,
        warnings=tuple(warnings),
    )
