"""The procedure a compressor, a pump and an expander have in common.

NeqSim has three classes for the three machines and the same eleven lines in each. Read
out of ``Compressor.java:1721-1776``, ``Pump.java:558-572`` and ``Expander.java:646-660``,
all three do:

    s_in  = S(T_in, P_in)                        the entropy the fluid arrives with
    T_is  = PSflash(P_out, s_in)                 where it would get to ideally
    h_is  = H(T_is, P_out)                       the ideal outlet enthalpy
    h_out = h_in + (h_is - h_in) / eta           compressor and pump  (Compressor.java:1738)
    h_out = h_in + (h_is - h_in) * eta           expander             (Expander.java:653)
    T_out = PHflash(P_out, h_out)                the real outlet state

The only difference between the three is **which way the efficiency scales the ideal
enthalpy change**, and that is what :class:`Direction` names. Dividing gives an outlet
further from the inlet than ideal, which is a machine consuming work; multiplying gives
one closer to the inlet, which is a machine producing it.

The mirror of ``crates/azoth-process/src/isentropic.rs``, and it exists for the same
reason that one does: three copies of eleven lines is three places for the efficiency to
end up on the wrong side of the division.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import StrEnum

from azoth.core.errors import InvalidInputError
from azoth.core.result import Phase
from azoth.core.units import Q, from_si, to_si
from azoth.core.warnings import Warning
from azoth.eos.mixture import Mixture
from azoth.eos.reference.molar_enthalpy_entropy import IdealGasModel
from azoth.eos.reference.ph_flash import enthalpy_at, ph_flash
from azoth.eos.reference.ps_flash import entropy_at, ps_flash


class Direction(StrEnum):
    """Which way a machine moves the pressure, and so which way the efficiency scales."""

    #: Work in. The efficiency **divides**, so the real outlet is further from the inlet
    #: than the ideal one.
    CONSUMING = "consuming"

    #: Work out. The efficiency **multiplies**, so the real outlet is closer to the inlet
    #: than the ideal one.
    PRODUCING = "producing"


@dataclass(frozen=True, slots=True)
class Solved:
    """Everything the three machines compute, before it is put into their own results.

    A dataclass rather than a dict, because the three callers each read six of these
    fields and a dict would make every read an unchecked lookup - which is exactly the
    kind of silent wrong-field mistake the three result classes exist to avoid.
    """

    #: The real outlet temperature.
    T: Q
    #: The shaft power, signed: positive into the fluid.
    power: float
    #: The vapour fraction at the outlet.
    beta: float | None
    #: Which phase the stream is in at the outlet.
    phase: Phase
    #: The temperature at unit efficiency.
    isentropic_temperature: Q
    #: Iterations, summed over the flashes this ran.
    iterations: int
    #: Caveats from the flashes.
    warnings: tuple[Warning, ...]


def isentropic_run(
    mixture: Mixture,
    ideal_gas: IdealGasModel,
    t_si: float,
    p_si: float,
    p_out_si: float,
    n_si: float,
    z: list[float],
    efficiency: float,
    direction: Direction,
) -> Solved:
    """Run one machine and return the six things its three result types are built from.

    ``efficiency`` is the **isentropic** efficiency, in ``(0, 1]``. NeqSim allows exactly
    that range too, clamped rather than refused (``Compressor.java:2118``); here the
    spec's range check refuses it, because a clamp turns a caller's mistake into a
    slightly different answer that looks deliberate.

    Raises:
        InvalidInputError: if the pressure moves the wrong way for the direction claimed.
            The schema's range checks are per quantity and cannot express a relation
            *between* two inputs, so this is the only place the check can live.
        OutOfRangeError: from the flashes, for a state with no root.
        SolverNotConvergedError: if a flash cannot reach its target on its bracket.
    """
    consistent = p_out_si > p_si if direction is Direction.CONSUMING else p_out_si < p_si
    if not consistent:
        machine = "a compressor or a pump" if direction is Direction.CONSUMING else "an expander"
        sense = "above" if direction is Direction.CONSUMING else "below"
        raise InvalidInputError(
            "outlet_pressure",
            f"{machine} must leave at a pressure {sense} its inlet; this one goes from "
            f"{p_si} Pa to {p_out_si} Pa. A machine that moves the pressure the other way "
            f"is a different unit operation, and admitting it here would return an answer "
            f"for a process nobody can build",
        )

    # The state the fluid arrives in. Both properties come from the same model at the same
    # state, so they describe one inlet and not two.
    h_in, _ = enthalpy_at(mixture, ideal_gas, t_si, p_si, z)
    s_in, _ = entropy_at(mixture, ideal_gas, t_si, p_si, z)

    # Where it would get to if the machine were perfect.
    ideal = ps_flash(mixture, ideal_gas, from_si(p_out_si, "Pa"), from_si(s_in, "J/(mol*K)"), z)
    t_is = to_si(ideal.T, "K", "T")
    h_ideal, _ = enthalpy_at(mixture, ideal_gas, t_is, p_out_si, z)

    # The whole difference between the three machines, at NeqSim's own line numbers.
    if direction is Direction.CONSUMING:
        d_h = (h_ideal - h_in) / efficiency
    else:
        d_h = (h_ideal - h_in) * efficiency

    flash = ph_flash(mixture, ideal_gas, from_si(p_out_si, "Pa"), from_si(h_in + d_h, "J/mol"), z)

    warnings: list[Warning] = list(ideal.warnings)
    warnings.extend(flash.warnings)

    return Solved(
        T=flash.T,
        power=n_si * d_h,
        beta=flash.beta,
        phase=flash.phase,
        isentropic_temperature=ideal.T,
        iterations=ideal.iterations + flash.iterations,
        warnings=tuple(warnings),
    )
