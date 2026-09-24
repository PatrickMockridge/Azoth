"""``process.gas_scrubber`` - the separator's kernel under the other entry.

Spec: ``specs/models/process/gas_scrubber.toml``

The Python twin of ``crates/azoth-process/src/models/separator.rs`` over
``crates/azoth-process/src/kernels/separator.rs``, written to mirror them.

# The vessel holds the temperature, and that is not a throttling

``Separator.run`` sets the pressure to ``P - pressureDrop`` and runs a ``TPflash``, so the
outlets come out at the **feed's** temperature - not at a temperature a caller states, and
not at the one an isenthalpic drop would give. Measured on the captured fluid, a 2 bar drop
takes -2705.35 J/mol in to a weighted -2364.95 out. That is the class's own model of the
vessel, and reproducing it means the balance is open whenever a pressure drop is asked for.
A ``heat_input`` is the one thing that moves the flash, and it moves it by ``Q / n`` J/mol -
``Heater.run``'s one-second convention, where the duty is W and the flow is mol/s.

# Entrainment moves material and not energy

``gas_in_liquid`` carries that fraction of the vapour's moles into the liquid outlet,
component by component, which is NeqSim's ``addPhaseFractionToPhase`` on its
``feed``/``mole`` basis. Both outlets are then at the same temperature, so the transfer
does not close the balance either - and it acts only where both phases exist, because a
from-phase or a to-phase that is not there has nothing to move.
"""

from __future__ import annotations

from azoth.core.errors import InvalidInputError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import GasScrubberResult, Phase
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning
from azoth.eos import components as _components
from azoth.eos.reference.ph_flash import enthalpy_at
from azoth.eos.reference.ph_flash import ph_flash as ph_flash_solve
from azoth.eos.reference.pt_flash import pt_flash
from azoth.process.reference._phase_outlet import phase_enthalpy


def gas_scrubber(
    components: list[str],
    feed_n: Q,
    feed_z: list[float],
    feed_p: Q,
    feed_t: Q,
    pressure_drop: Q,
    gas_in_liquid: float,
    heat_input: Q | None = None,
) -> GasScrubberResult:
    """Flash a stream into vapour and liquid outlets.

    Args:
        components: the fluid's substances, by name.
        feed_n: the inlet's molar flow.
        feed_z: the inlet composition.
        feed_p: the inlet pressure.
        feed_t: the inlet temperature, which the vessel holds.
        pressure_drop: the pressure drop across the vessel.
        heat_input: heat added before the flash, or ``None`` for a flash at the feed's
            temperature.
        gas_in_liquid: the fraction of the vapour's moles carried into the liquid outlet.

    Returns:
        Both outlets' records.

    Raises:
        InvalidInputError: where the pressure drop takes the outlet below zero or the
            entrainment fraction is outside ``[0, 1]``.
        OutOfRangeError: for a pressure drop or a fraction the spec bounds.

    Example:
        >>> import azoth
        >>> q = azoth.ureg.Quantity
        >>> r = separator(
        ...     ["methane", "n-butane"],
        ...     q(1.0, "mol/s"),
        ...     [0.7, 0.3],
        ...     q(20.0, "bar"),
        ...     q(300.0, "K"),
        ...     q(0.0, "Pa"),
        ... )
        >>> round(r.vapour_n.to("mol/s").magnitude, 6)
        0.818221
    """
    spec = _spec()
    checks = checks_for(spec)
    warnings: list[Warning] = []

    n = input_to_si(spec, "feed_n", feed_n)
    t = input_to_si(spec, "feed_t", feed_t)
    p = input_to_si(spec, "feed_p", feed_p)
    drop = input_to_si(spec, "pressure_drop", pressure_drop)

    apply_checks(
        checks.on_input,
        {"pressure_drop": drop, "gas_in_liquid": gas_in_liquid}.get,
        warnings,
    )

    if not 0.0 <= gas_in_liquid <= 1.0:
        raise InvalidInputError(
            "gas_in_liquid",
            f"an entrainment fraction is a fraction, and {gas_in_liquid} is not in [0, 1]",
        )
    p_out = p - drop
    if p_out <= 0.0:
        raise InvalidInputError(
            "pressure_drop",
            f"a pressure drop of {drop} Pa takes a {p} Pa inlet to {p_out} Pa, which is "
            f"not a pressure",
        )

    mixture, ideal_gas = _components.mixture_of(components, eos="pr")

    if heat_input is None:
        # The vessel holds the feed's temperature, so the flash is a TP one and the
        # outlets share the temperature they came in at.
        flash = pt_flash(mixture, feed_t, from_si(p_out, "Pa"), feed_z)
        temperature = feed_t
        x, y = list(flash.x), list(flash.y)
        beta = _vapour_fraction(flash.phase, flash.beta)
    else:
        # W over (mol/s) is J/mol, for the reason `Heater.run` gives.
        h_in, _ = enthalpy_at(mixture, ideal_gas, t, p, feed_z)
        duty = input_to_si(spec, "heat_input", heat_input)
        moved = ph_flash_solve(
            mixture,
            ideal_gas,
            from_si(p_out, "Pa"),
            from_si(h_in + duty / n, "J/mol"),
            feed_z,
        )
        temperature = moved.T
        x, y = list(moved.x), list(moved.y)
        beta = _vapour_fraction(moved.phase, moved.beta)

    n_vapour = n * beta
    n_liquid = n * (1.0 - beta)
    # Only where both phases exist: a from-phase or a to-phase that is not there has
    # nothing to move, which is NeqSim's own guard.
    carried = n_vapour * gas_in_liquid if n_vapour > 0.0 and n_liquid > 0.0 else 0.0

    liquid_n = n_liquid + carried
    liquid_z = (
        [(n_liquid * x[i] + carried * y[i]) / liquid_n for i in range(len(x))]
        if carried > 0.0
        else x
    )

    # **The outlet is a phase, except where the class re-runs it as a stream.** `GasScrubber`
    # inherits `Separator.run`, so the same branch decides the liquid's enthalpy: a phase root
    # while nothing was carried into it, a re-flash of the carried composition once something
    # was. See `_phase_outlet`, which is where the difference is measured.
    t_out = temperature.to("K").magnitude
    vapour_h = phase_enthalpy(mixture, ideal_gas, t_out, p_out, y, liquid=False)
    if gas_in_liquid != 0.0:
        liquid_h, _ = enthalpy_at(mixture, ideal_gas, t_out, p_out, list(liquid_z))
    else:
        liquid_h = phase_enthalpy(mixture, ideal_gas, t_out, p_out, list(liquid_z), liquid=True)

    return GasScrubberResult(
        vapour_n=from_si(n_vapour - carried, "mol/s"),
        vapour_z=tuple(y),
        vapour_p=from_si(p_out, "Pa"),
        vapour_t=temperature,
        vapour_h=from_si(vapour_h, "J/mol"),
        liquid_n=from_si(liquid_n, "mol/s"),
        liquid_z=tuple(liquid_z),
        liquid_p=from_si(p_out, "Pa"),
        liquid_t=temperature,
        liquid_h=from_si(liquid_h, "J/mol"),
        warnings=tuple(warnings),
    )


def _vapour_fraction(phase: Phase, beta: float | None) -> float:
    """The vapour fraction a flash's phase and ``beta`` imply.

    ``beta`` is absent for a single-phase answer - the flash reports ``None`` rather than a
    number outside ``[0, 1]``, because there is no split to report - so the phase is what
    says which side of the interval the outlet is on.
    """
    if phase in (Phase.ALL_VAPOUR, Phase.TRIVIAL):
        return 1.0
    return 0.0 if beta is None else float(beta)


def _spec() -> dict[str, object]:
    """The generated spec, imported here so the module has no import-time cycle."""
    from azoth._models_gen import model

    return model("process.gas_scrubber")
