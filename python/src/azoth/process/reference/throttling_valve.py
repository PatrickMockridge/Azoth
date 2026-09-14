"""``process.throttling_valve`` - pressure dropped at constant enthalpy.

Spec: ``specs/models/process/throttling_valve.yaml``

The port source is ``neqsim.process.equipment.valve.ThrottlingValve``, ``run(UUID)`` at
lines 363-453 of the 3.20.0 tree, with the isenthalpic flash in
``runPHflashWithNaNRetry`` at ``:467-495``. Two statements of physics:

    P_out = P_in - deltaP        (:407-412)
    T_out = PHflash(P_out, H_in) (:470)

**The temperature really does change.** A throttling valve is the classic Joule-Thomson
device, and a model that returned the inlet temperature would be wrong for every real gas.
The ``PHflash`` is what makes it right, and it is the reason this unit operation could not
be written before ``eos.ph_flash`` existed.

**There is no flow rate in this signature**, and that is a statement about the physics
rather than a saving: an isenthalpic flash is a molar property, so the outlet state does
not depend on the flow at all. NeqSim computes a molar flow here too and applies it only
on the transient path.
"""

from __future__ import annotations

from azoth import _models_gen
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import ThrottlingValveResult
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning
from azoth.eos.mixture import Mixture
from azoth.eos.reference.molar_enthalpy_entropy import IdealGasModel
from azoth.eos.reference.ph_flash import enthalpy_at, ph_flash

MODEL_ID = "process.throttling_valve"


def throttling_valve(
    mixture: Mixture,
    ideal_gas: IdealGasModel,
    T: Q,
    P: Q,
    z: list[float],
    pressure_drop: Q,
) -> ThrottlingValveResult:
    """A stream's pressure dropped at constant enthalpy.

    Args:
        mixture: the components and their interaction parameters.
        ideal_gas: the datum the enthalpy is measured from.
        T: the feed's absolute temperature.
        P: the feed's absolute pressure, before the drop.
        z: the feed's mole fractions.
        pressure_drop: the pressure the valve drops, subtracted from ``P``.

    Returns:
        The outlet temperature and pressure, and the phase there.

    Raises:
        OutOfRangeError: if an input is outside the spec's range, or the drop takes the
            outlet to or below zero - the flash catches that last one.
        InvalidInputError: if ``z`` is not a composition.
        SolverNotConvergedError: if the outlet state is outside the flash's bracket.
    """
    spec = _models_gen.model(MODEL_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    t_si = input_to_si(spec, "T", T)
    p_si = input_to_si(spec, "P", P)
    drop_si = input_to_si(spec, "pressure_drop", pressure_drop)
    apply_checks(
        checks.on_input,
        {"T": t_si, "P": p_si, "pressure_drop": drop_si}.get,
        warnings,
    )

    p_out_si = p_si - drop_si

    # The enthalpy the fluid arrived with. This is the whole of the model: what leaves has
    # the same enthalpy and less pressure, and everything else follows from that.
    h_in, _ = enthalpy_at(mixture, ideal_gas, t_si, p_si, list(z))
    flash = ph_flash(mixture, ideal_gas, from_si(p_out_si, "Pa"), from_si(h_in, "J/mol"), list(z))
    warnings.extend(flash.warnings)

    return ThrottlingValveResult(
        T=flash.T,
        P=from_si(p_out_si, "Pa"),
        phase=flash.phase,
        beta=flash.beta,
        iterations=flash.iterations,
        warnings=tuple(warnings),
    )
