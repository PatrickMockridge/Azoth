"""``process.expander`` - a pressure drop that produces work.

Spec: ``specs/models/process/expander.yaml``

The port source is ``neqsim.process.equipment.expander.Expander``, ``run(UUID)`` at lines
608-668 of a 669-line file. ``Expander`` **extends ``Compressor``** and overrides ``run``,
and the difference from the compressor is one expression:

    Compressor.java:1738   dH = (H(P_out, s_in) - H_in) / isentropicEfficiency
    Expander.java:653      dH = (H(P_out, s_in) - H_in) * isentropicEfficiency

The enthalpy change of an expansion is negative, so the multiplication makes it less
negative and the real outlet **warmer** than the ideal one. That is the physical statement:
an expander cannot deliver more work than the isentropic drop contains, so its outlet
cannot be colder than the isentropic temperature.

The procedure itself is in :mod:`azoth.process.reference._isentropic`, and the only thing
this module passes that the compressor does not is ``Direction.PRODUCING``.
"""

from __future__ import annotations

from azoth import _models_gen
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import ExpanderResult
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning
from azoth.eos.mixture import Mixture
from azoth.eos.reference.molar_enthalpy_entropy import IdealGasModel
from azoth.process.reference._isentropic import Direction, isentropic_run

MODEL_ID = "process.expander"


def expander(
    mixture: Mixture,
    ideal_gas: IdealGasModel,
    T: Q,
    P: Q,
    n: Q,
    z: list[float],
    outlet_pressure: Q,
    efficiency: float,
) -> ExpanderResult:
    """A stream expanded to a stated outlet pressure, producing work.

    Args:
        mixture: the components and their interaction parameters.
        ideal_gas: the datum the enthalpy and the entropy are measured from.
        T: the inlet absolute temperature.
        P: the inlet absolute pressure. The outlet must be **below** it.
        n: the inlet molar flow rate.
        z: the inlet mole fractions.
        outlet_pressure: the pressure the expander discharges to.
        efficiency: the isentropic efficiency, in ``(0, 1]``. A *smaller* number gives a
            warmer outlet and less power, which is the opposite of the compressor's sense.

    Returns:
        The outlet temperature and phase, the shaft power - **negative**, because the
        fluid is doing the work - and the coldest temperature this pressure drop could
        reach.

    Raises:
        InvalidInputError: if ``outlet_pressure`` is not below the inlet pressure.
        OutOfRangeError: if an input is outside the spec's declared range.
        SolverNotConvergedError: if either flash cannot reach its target.

    Note:
        ``power`` is negative here and positive for a compressor, on the same convention:
        positive means energy into the fluid. NeqSim reports the opposite sign on its
        energy port (``Expander.java:661``, ``setDuty(-dH)``), which is worth knowing
        before comparing the two.
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
        mixture, ideal_gas, t_si, p_si, p_out_si, n_si, list(z), efficiency, Direction.PRODUCING
    )
    warnings.extend(solved.warnings)

    return ExpanderResult(
        T=solved.T,
        P=from_si(p_out_si, "Pa"),
        power=solved.power,
        beta=solved.beta,
        phase=solved.phase,
        isentropic_temperature=solved.isentropic_temperature,
        iterations=solved.iterations,
        warnings=tuple(warnings),
    )
