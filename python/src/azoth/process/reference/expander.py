"""``process.expander`` - the isentropic route, multiplying by the efficiency.

Spec: ``specs/models/process/expander.toml``

The Python twin of ``crates/azoth-process/src/models/expander.rs`` over
``crates/azoth-process/src/kernels/expander.rs``, written to mirror them.

# The step is isentropic, and the inlet's entropy is derived

NeqSim reads ``getEntropy()`` off the inlet system, sets the outlet pressure, runs a
``PSflash`` on that entropy and takes the enthalpy it lands on as the reversible outlet. The
record has no entropy field - ``s = f(T, P, z)``, the same flash that gives ``h`` - so a port
that agrees with the captured efficiency-of-one rows has computed it rather than carried it.

# The efficiency multiplies the isentropic step

``h_out = h_in + (h_isentropic - h_in) * eta``, so an efficiency below one recovers
less work than the reversible step. That is ``Compressor.run``'s rule with the efficiency on
the other side of the division, and it is the same physical rule: an expansion's difference is
negative, so dividing would make the machine beat the reversible one.
"""

from __future__ import annotations

from dataclasses import dataclass

from azoth.core.errors import InvalidInputError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import ExpanderResult
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning
from azoth.eos import components as _components
from azoth.eos.reference.ph_flash import enthalpy_at
from azoth.eos.reference.ph_flash import ph_flash as ph_flash_solve
from azoth.eos.reference.ps_flash import entropy_at
from azoth.eos.reference.ps_flash import ps_flash as ps_flash_solve


def expander(
    components: list[str],
    inlet_n: Q,
    inlet_z: list[float],
    inlet_p: Q,
    inlet_t: Q,
    outlet_pressure: Q,
    isentropic_efficiency: float,
) -> ExpanderResult:
    """Drop a stream's pressure along an isentrope, multiplying the step by the efficiency.

    Args:
        components: the fluid's substances, by name.
        inlet_n: the inlet's molar flow.
        inlet_z: the inlet composition.
        inlet_p: the inlet pressure.
        inlet_t: the inlet temperature.
        outlet_pressure: the pressure the machine takes the stream to.
        isentropic_efficiency: the machine's isentropic efficiency, in ``(0, 1]``.

    Returns:
        The outlet's record.

    Example:
        >>> import azoth
        >>> q = azoth.ureg.Quantity
        >>> r = expander(
        ...     ["methane", "n-butane"],
        ...     1.0,
        ...     q([0.9, 0.1]),
        ...     q(60.0, "bar"),
        ...     q(320.0, "K"),
        ...     q(30.0, "bar"),
        ...     0.75,
        ... )
        >>> round(r.outlet_h.to("J/mol").magnitude, 3)
        -618.933
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
            f"the machine's isentropic efficiency must be in (0, 1], and "
            f"{isentropic_efficiency} is not",
        )

    states = _states(components, t, p, inlet_z, p_out, isentropic_efficiency)

    return ExpanderResult(
        outlet_n=from_si(n, "mol/s"),
        outlet_z=tuple(inlet_z),
        outlet_p=outlet_pressure,
        outlet_t=states.outlet_t,
        outlet_h=from_si(states.h_out, "J/mol"),
        warnings=tuple(warnings),
    )


@dataclass(frozen=True, slots=True)
class ExpanderStates:
    """The machine's intermediates, in the order it computes them, in SI."""

    #: The inlet's molar enthalpy and entropy, J/mol and J/(mol*K).
    h_in: float
    s_in: float
    #: The reversible outlet: the temperature at the outlet pressure on the inlet's entropy,
    #: and its molar enthalpy.
    t_isentropic: Q
    h_isentropic: float
    #: The isentropic step, and the step the efficiency actually asks for.
    dh_isentropic: float
    dh_actual: float
    #: The outlet: its molar enthalpy, its temperature, and its molar entropy - which rises
    #: across an irreversible step and is unchanged across a reversible one.
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
) -> ExpanderStates:
    """The machine's intermediates, in SI, in the order it computes them."""
    mixture, ideal_gas = _components.mixture_of(components, eos="pr")

    h_in, _ = enthalpy_at(mixture, ideal_gas, t, p, inlet_z)
    entropy, _ = entropy_at(mixture, ideal_gas, t, p, inlet_z)

    isentropic = ps_flash_solve(
        mixture, ideal_gas, from_si(p_out, "Pa"), from_si(entropy, "J/(mol*K)"), inlet_z
    )
    h_isentropic, _ = enthalpy_at(
        mixture, ideal_gas, isentropic.T.to("K").magnitude, p_out, inlet_z
    )

    dh_isentropic = h_isentropic - h_in
    dh_actual = dh_isentropic * isentropic_efficiency
    h_out = h_in + dh_actual
    outlet_t = ph_flash_solve(
        mixture, ideal_gas, from_si(p_out, "Pa"), from_si(h_out, "J/mol"), inlet_z
    )

    s_out, _ = entropy_at(mixture, ideal_gas, outlet_t.T.to("K").magnitude, p_out, inlet_z)

    return ExpanderStates(
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

    return model("process.expander")
