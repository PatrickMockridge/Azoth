"""``process.heater`` - a duty applied at a fixed pressure.

Spec: ``specs/models/process/heater.yaml``, which carries the provenance and the reason
one of the four specifications is ported and the others are not.

A cooler is this model with a negative ``heat_duty``; the sign is the whole difference.
"""

from __future__ import annotations

from azoth import _models_gen
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import HeaterResult
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning
from azoth.eos.mixture import Mixture
from azoth.eos.reference.molar_enthalpy_entropy import IdealGasModel
from azoth.eos.reference.ph_flash import enthalpy_at, ph_flash

MODEL_ID = "process.heater"


def heater(
    mixture: Mixture,
    ideal_gas: IdealGasModel,
    T: Q,
    P: Q,
    n: Q,
    z: list[float],
    pressure_drop: Q,
    heat_duty: Q,
) -> HeaterResult:
    """A duty applied to a stream at a fixed pressure.

    Args:
        mixture: the components and their interaction parameters.
        ideal_gas: the datum the enthalpy is measured from.
        T: the feed's absolute temperature.
        P: the feed's absolute pressure, before the drop.
        n: the feed's molar flow rate. Needed rather than decorative - the duty is an
            *extensive* quantity and the flash inverts a *molar* enthalpy, so the
            conversion ``Q / n`` is part of the model. The spec's range check puts it
            strictly above zero, which is what makes the division safe.
        z: the feed's mole fractions.
        pressure_drop: the pressure the exchanger drops, subtracted from ``P`` first.
        heat_duty: the duty in watts, **signed**. Positive raises the temperature;
            negative is a cooler.

    Returns:
        The outlet temperature - the model's answer - and the pressure and phase there.

    Raises:
        OutOfRangeError: if an input is outside the spec's range, including a zero flow,
            or if the flash cannot find the requested enthalpy on its bracket.
        InvalidInputError: if ``z`` is not a composition.
    """
    spec = _models_gen.model(MODEL_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    t_si = input_to_si(spec, "T", T)
    p_si = input_to_si(spec, "P", P)
    n_si = input_to_si(spec, "n", n)
    drop_si = input_to_si(spec, "pressure_drop", pressure_drop)
    duty_si = input_to_si(spec, "heat_duty", heat_duty)
    apply_checks(
        checks.on_input,
        {"T": t_si, "P": p_si, "n": n_si, "pressure_drop": drop_si}.get,
        warnings,
    )

    p_out_si = p_si - drop_si

    h_in, _ = enthalpy_at(mixture, ideal_gas, t_si, p_si, list(z))
    flash = ph_flash(
        mixture,
        ideal_gas,
        from_si(p_out_si, "Pa"),
        from_si(h_in + duty_si / n_si, "J/mol"),
        list(z),
    )
    warnings.extend(flash.warnings)

    return HeaterResult(
        T=flash.T,
        P=from_si(p_out_si, "Pa"),
        phase=flash.phase,
        beta=flash.beta,
        iterations=flash.iterations,
        warnings=tuple(warnings),
    )
