"""``process.pump`` - the pump's kernel.

Spec: ``specs/models/process/pump.toml``

The Python twin of ``crates/azoth-process/src/models/pump.rs`` over
``crates/azoth-process/src/kernels/pump.rs``, written to mirror them.

# The work is an isentropic head, not ``v (P_out - P_in)``

NeqSim's ``Pump.run`` clones the fluid to the outlet pressure and runs a ``PSflash`` at the
inlet entropy, because ``calculateAsCompressor`` defaults to ``true`` and selects that
branch. So the work is ``(H(P_out, s_in) - H(P_in)) / eta`` - the enthalpy a real fluid
needs along an isentrope - and the incompressible form is its linearisation. Measured on the
captured fluid the two differ by ``1.2e-3``, which is why this is the isentropic route.

# The stream record

An inlet carries ``n, z, P, T`` and an outlet carries ``n, z, P, T, h``. ``h`` is a state
function of ``(T, P, z)`` - the same flash decides it - so accepting one as well would let a
case hand over a state that does not exist.
"""

from __future__ import annotations

from dataclasses import dataclass

from azoth.core.errors import InvalidInputError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import PumpResult
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning
from azoth.eos import components as _components

# **The modules' functions, imported by name.** `azoth.eos.reference` re-exports `ph_flash`
# and `ps_flash` as *functions*, so reaching the modules through the package gives the
# function where a module is meant. The `from ... import` form resolves the module path
# itself and is unambiguous.
from azoth.eos.reference.ph_flash import enthalpy_at
from azoth.eos.reference.ph_flash import ph_flash as ph_flash_solve
from azoth.eos.reference.ps_flash import entropy_at
from azoth.eos.reference.ps_flash import ps_flash as ps_flash_solve


def pump(
    components: list[str],
    inlet_n: Q,
    inlet_z: list[float],
    inlet_p: Q,
    inlet_t: Q,
    outlet_pressure: Q,
    isentropic_efficiency: float,
) -> PumpResult:
    """Raise a stream's pressure, adding the pump's work as enthalpy.

    Args:
        components: the fluid's substances, by name.
        inlet_n: the inlet's molar flow.
        inlet_z: the inlet composition.
        inlet_p: the inlet pressure.
        inlet_t: the inlet temperature.
        outlet_pressure: the pressure the pump raises the stream to.
        isentropic_efficiency: the pump's isentropic efficiency, in ``(0, 1]``.

    Returns:
        The outlet's record: its flow and composition (the inlet's), its pressure, and the
        temperature and molar enthalpy the isentropic head implies.

    Raises:
        InvalidInputError: where the shapes disagree or the efficiency is outside ``(0, 1]``.
        OutOfRangeError: for a temperature the spec bounds.

    Example:
        >>> import azoth
        >>> q = azoth.ureg.Quantity
        >>> r = pump(
        ...     ["n-butane"],
        ...     1.0,
        ...     q(1.0, "mol/s"),
        ...     q(5.0, "bar"),
        ...     q(250.0, "K"),
        ...     q(20.0, "bar"),
        ...     0.75,
        ... )
        >>> round(r.outlet_t.to("K").magnitude, 4)
        250.7621
    """
    spec = _spec()
    checks = checks_for(spec)
    warnings: list[Warning] = []

    n = input_to_si(spec, "inlet_n", inlet_n)
    t = input_to_si(spec, "inlet_t", inlet_t)
    p = input_to_si(spec, "inlet_p", inlet_p)
    p_out = input_to_si(spec, "outlet_pressure", outlet_pressure)

    apply_checks(
        checks.on_input,
        {"inlet_t": t, "isentropic_efficiency": isentropic_efficiency}.get,
        warnings,
    )

    if not 0.0 < isentropic_efficiency <= 1.0:
        raise InvalidInputError(
            "isentropic_efficiency",
            f"the pump's isentropic efficiency must be in (0, 1], and {isentropic_efficiency} "
            f"is not",
        )

    states = _states(components, t, p, inlet_z, p_out, isentropic_efficiency)

    return PumpResult(
        outlet_n=from_si(n, "mol/s"),
        outlet_z=tuple(inlet_z),
        outlet_p=outlet_pressure,
        outlet_t=states.outlet_t,
        outlet_h=from_si(states.h_out, "J/mol"),
        warnings=tuple(warnings),
    )


@dataclass(frozen=True, slots=True)
class PumpStates:
    """The pump's intermediates, in the order it computes them, in SI.

    A record rather than a dict so that a layer the harness reads and a field the result
    carries are the same typed thing, and neither can acquire a name the other lacks.
    """

    #: The inlet's molar enthalpy and entropy, J/mol and J/(mol*K).
    h_in: float
    s_in: float
    #: The isentropic outlet: the temperature at the raised pressure on the inlet's
    #: entropy, and its molar enthalpy.
    t_isentropic: Q
    h_isentropic: float
    #: The isentropic head, and the head the efficiency actually asks for.
    dh_isentropic: float
    dh_actual: float
    #: The outlet: its molar enthalpy, its temperature, and its molar entropy - which is
    #: the second-law statement of the step, since it rises iff the step was irreversible.
    h_out: float
    outlet_t: Q
    s_out: float


def _states(
    components: list[str],
    t: float,
    p: float,
    inlet_z: list[float],
    p_out: float,
    isentropic_efficiency: float,
) -> PumpStates:
    """The pump's intermediates, in the order it computes them, in SI.

    **One arithmetic, two consumers.** :func:`pump` builds the result from this and
    :mod:`azoth.process.layers` reports the same numbers against a NeqSim capture, so a
    layer the harness compares is a layer the model actually took - which is the rule
    `azoth.eos.layers` states and the reason this is factored rather than restated.
    """
    mixture, ideal_gas = _components.mixture_of(components, eos="pr")

    # The inlet's own state, which is where the entropy and the enthalpy come from.
    h_in, _ = enthalpy_at(mixture, ideal_gas, t, p, inlet_z)
    entropy, _ = entropy_at(mixture, ideal_gas, t, p, inlet_z)

    isentropic = ps_flash_solve(
        mixture, ideal_gas, from_si(p_out, "Pa"), from_si(entropy, "J/(mol*K)"), inlet_z
    )
    h_isentropic, _ = enthalpy_at(
        mixture, ideal_gas, isentropic.T.to("K").magnitude, p_out, inlet_z
    )

    dh_isentropic = h_isentropic - h_in
    dh_actual = dh_isentropic / isentropic_efficiency
    h_out = h_in + dh_actual
    outlet_t = ph_flash_solve(
        mixture, ideal_gas, from_si(p_out, "Pa"), from_si(h_out, "J/mol"), inlet_z
    )

    # The outlet's own entropy, which is the second-law statement of the head: it rises
    # across an irreversible step and is unchanged across a reversible one. `Pump` exposes
    # it as `getEntropyProduction` and nothing else about the machine's interior, which is
    # why the harness localises a wrong head through it.
    s_out, _ = entropy_at(mixture, ideal_gas, outlet_t.T.to("K").magnitude, p_out, inlet_z)

    return PumpStates(
        h_in=h_in,
        s_in=entropy,
        t_isentropic=isentropic.T,
        h_isentropic=h_isentropic,
        dh_isentropic=dh_isentropic,
        dh_actual=dh_actual,
        h_out=h_out,
        outlet_t=outlet_t.T,
        s_out=s_out,
    )


def _spec() -> dict[str, object]:
    """The generated spec, imported here so the module has no import-time cycle."""
    from azoth._models_gen import model

    return model("process.pump")
