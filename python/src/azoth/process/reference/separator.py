"""``process.separator`` - one feed split into a gas and a liquid at a single state.

Spec: ``specs/models/process/separator.yaml``

The first unit operation ported from NeqSim, and the pure-Python half of it. Its
physics is one flash and the arithmetic that follows from it: the vapour fraction
splits the molar flow, the equilibrium compositions become the two outlets', and both
outlets leave at the flash's temperature and pressure because that is what a vessel
held at one state does.

The port source is ``neqsim.process.equipment.separator.Separator``, ``run(UUID)`` at
lines 674-782 of the 3.20.0 tree - 109 lines of a 4,511-line file. What is taken from
them is the pressure drop, the choice between an isothermal and an isenthalpic flash,
and building an outlet from a flashed phase. What is not is the internal mixer, the
memoization guard, the low-flow bypass and the entrainment model; the spec's notes
record each with its reason.

**The split is decided on ``phase`` and never on ``beta``.** A single-phase feed has no
vapour fraction, and the number the flash would report for one is the *extrapolated*
split - values outside ``[0, 1]`` are ordinary. Multiplying the feed by 1.888 would give
a wrong answer shaped exactly like a right one.
"""

from __future__ import annotations

from azoth import _models_gen
from azoth.core.errors import InvalidInputError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import Phase, SeparatorResult
from azoth.core.units import Q, from_si, input_to_si, to_si
from azoth.core.warnings import Warning
from azoth.eos.mixture import Mixture
from azoth.eos.reference.molar_enthalpy_entropy import IdealGasModel
from azoth.eos.reference.ph_flash import enthalpy_at, ph_flash
from azoth.eos.reference.pt_flash import pt_flash

MODEL_ID = "process.separator"


def separator(
    mixture: Mixture,
    ideal_gas: IdealGasModel,
    T: Q,
    P: Q,
    n: Q,
    z: list[float],
    pressure_drop: Q,
    heat_duty: Q,
) -> SeparatorResult:
    """One feed split into a gas and a liquid at a single temperature and pressure.

    Args:
        mixture: the components and their interaction parameters.
        ideal_gas: the ``Cp/R`` polynomial and reference state the duty branch measures
            its enthalpy from. Required even when ``heat_duty`` is zero, because the
            model's signature is fixed and the duty branch is the same model.
        T: the feed's absolute temperature. It is the outlet temperature too unless a
            duty is set.
        P: the feed's absolute pressure, before the pressure drop.
        n: the feed's molar flow rate. Every answer is linear in it.
        z: the feed's mole fractions. Checked rather than renormalised.
        pressure_drop: the pressure the vessel drops, subtracted from ``P`` first.
        heat_duty: the duty applied to the feed, in watts. Zero means the separator is
            isothermal at its feed temperature and the procedure is a ``TPflash``.

    Returns:
        Both outlets at one temperature and pressure, with the flow and composition of
        each.

    Raises:
        InvalidInputError: if ``z`` is not a composition, or if the flash reports a
            two-phase split with no vapour fraction - a contradiction rather than a
            state.
        OutOfRangeError: if an input is outside the spec's declared range - a zero
            flow, a negative pressure drop, a non-positive absolute state.
        SolverNotConvergedError: if the flash hits its cap, or a duty puts the answer
            outside the bracket ``eos.ph_flash`` searches.
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

    if duty_si == 0.0:
        tp = pt_flash(mixture, from_si(t_si, "K"), from_si(p_out_si, "Pa"), list(z))
        temperature_si = t_si
        beta, phase, iterations = tp.beta, tp.phase, tp.iterations
        x, y = list(tp.x), list(tp.y)
        warnings.extend(tp.warnings)
    else:
        # The feed's enthalpy at its *own* state, evaluated explicitly. NeqSim adds the
        # duty to whatever the last flash left on the object; this depends on the
        # inputs and on nothing else. See the spec's notes.
        h_in, _ = enthalpy_at(mixture, ideal_gas, t_si, p_si, list(z))
        ph = ph_flash(
            mixture,
            ideal_gas,
            from_si(p_out_si, "Pa"),
            from_si(h_in + duty_si / n_si, "J/mol"),
            list(z),
        )
        temperature_si = to_si(ph.T, "K", "T")
        beta, phase, iterations = ph.beta, ph.phase, ph.iterations
        x, y = list(ph.x), list(ph.y)
        warnings.extend(ph.warnings)

    if phase == Phase.TWO_PHASE:
        if beta is None:
            raise InvalidInputError(
                "beta",
                f"{MODEL_ID} reports a two-phase split with no vapour fraction, so "
                "there is nothing to split the feed by",
            )
        gas_flow = beta * n_si
        liquid_flow = (1.0 - beta) * n_si
        gas_z, liquid_z = y, x
    elif phase == Phase.ALL_LIQUID:
        gas_flow, liquid_flow, gas_z, liquid_z = 0.0, n_si, list(z), x
    else:
        # ALL_VAPOUR, and TRIVIAL - where the flash converged to `x = y = z` and
        # established that the feed is single phase but not which one. The whole feed
        # goes to the gas outlet, which is a convention rather than a finding; `phase`
        # is what tells a caller it happened.
        gas_flow, liquid_flow, gas_z, liquid_z = n_si, 0.0, y, list(z)

    return SeparatorResult(
        T=from_si(temperature_si, "K"),
        P=from_si(p_out_si, "Pa"),
        beta=beta,
        gas_flow=gas_flow,
        gas_z=tuple(gas_z),
        liquid_flow=liquid_flow,
        liquid_z=tuple(liquid_z),
        phase=phase,
        iterations=iterations,
        warnings=tuple(warnings),
    )
