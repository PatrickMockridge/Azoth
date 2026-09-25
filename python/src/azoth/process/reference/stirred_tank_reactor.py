"""``process.stirred_tank_reactor`` - a stoichiometric conversion, then a flash.

Spec: ``specs/models/process/stirred_tank_reactor.toml``

The Python twin of ``crates/azoth-process/src/models/stirred_tank_reactor.rs`` over
``crates/azoth-process/src/kernels/stirred_tank_reactor.rs``, written to mirror them.

# The reaction is a movement of moles

``StoichiometricReaction.react`` takes ``moles(limiting) * conversion`` and adjusts every
component by ``that * coeff / |coeff(limiting)|`` - so the extent is the *limiting reactant's*
and every other coefficient is scaled by it. The stoichiometry comes from
``data/reactions/stoichiometry.csv``, resolved by the reaction's id.

**A row naming a substance the feed does not carry is skipped**, which is the class's own
behaviour: ``addComponent`` throws for a component the fluid has no slot for and the catch
drops it. A reaction whose products are not in the feed therefore moves less material than its
stoichiometry says, and this does not repair it.

# Then the flash, and what the duty is

``run`` sets the outlet pressure - a stated reactor pressure *replaces* the feed's rather than
dropping from it - holds the temperature when it is isothermal, and flashes: a ``TPflash`` at
the held temperature, or a ``PHflash`` at the **feed's own enthalpy** when it is adiabatic.
There is no energy balance over the reaction itself; the heat released is the temperature the
adiabatic flash lands on, and the duty the isothermal branch reports is
``outletEnthalpy - inletEnthalpy``.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

from azoth.core.errors import InvalidInputError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import StirredTankReactorResult
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning
from azoth.eos import components as _components
from azoth.eos.reference.ph_flash import enthalpy_at
from azoth.eos.reference.ph_flash import ph_flash as ph_flash_solve
from azoth.eos.reference.pt_flash import pt_flash
from azoth.reactions.reference import _tables


def stirred_tank_reactor(
    components: list[str],
    feed_n: Q,
    feed_z: list[float],
    feed_p: Q,
    feed_t: Q,
    reaction: str,
    limiting_reactant: str,
    conversion: float,
    isothermal: bool,
    reactor_temperature: Q | None = None,
    reactor_pressure: Q | None = None,
    pressure_drop: Q | None = None,
) -> StirredTankReactorResult:
    """React a feed stoichiometrically and flash the product.

    Args:
        components: the fluid's substances, by name.
        feed_n: the feed's molar flow.
        feed_z: the feed's composition.
        feed_p: the feed's pressure.
        feed_t: the feed's temperature.
        reaction: the reaction's id in the reaction data.
        limiting_reactant: the component whose moles times ``conversion`` set the extent.
        conversion: the fractional conversion of the limiting reactant, in ``[0, 1]``.
        isothermal: whether the vessel holds its temperature.
        reactor_temperature: the temperature an isothermal vessel holds.
        reactor_pressure: the pressure the vessel holds, or ``None`` for the feed's less the drop.
        pressure_drop: the pressure drop the outlet takes where no reactor pressure is stated.

    Returns:
        The product's record and the duty the vessel supplied.

    Raises:
        InvalidInputError: for a reaction the data does not carry, a limiting reactant that is
            not a feed component, a conversion outside ``[0, 1]``, or a conversion that would
            take a substance below zero.

    Example:
        >>> import azoth
        >>> q = azoth.ureg.Quantity
        >>> r = stirred_tank_reactor(
        ...     ["methane", "oxygen", "CO2", "water", "nitrogen"],
        ...     q(1.0, "mol/s"),
        ...     [0.05, 0.10, 0.02, 0.03, 0.80],
        ...     q(5.0, "bar"),
        ...     q(500.0, "K"),
        ...     "methanecombustion",
        ...     "oxygen",
        ...     0.5,
        ...     False,
        ... )
        >>> round(r.product_t.to("K").magnitude, 3)
        498.537
    """
    spec = _spec()
    checks = checks_for(spec)
    warnings: list[Warning] = []

    n = input_to_si(spec, "feed_n", feed_n)
    p = input_to_si(spec, "feed_p", feed_p)
    drop = 0.0 if pressure_drop is None else input_to_si(spec, "pressure_drop", pressure_drop)
    held = (
        None
        if reactor_temperature is None
        else input_to_si(spec, "reactor_temperature", reactor_temperature)
    )
    held_p = (
        None
        if reactor_pressure is None
        else input_to_si(spec, "reactor_pressure", reactor_pressure)
    )

    apply_checks(
        checks.on_input,
        {"conversion": conversion, "pressure_drop": drop}.get,
        warnings,
    )

    states = _route(
        components,
        n,
        feed_z,
        p,
        feed_t,
        reaction,
        limiting_reactant,
        conversion,
        isothermal,
        held,
        held_p,
        drop,
    )

    return StirredTankReactorResult(
        product_n=from_si(states.product_n, "mol/s"),
        product_z=states.product_z,
        product_p=from_si(states.product_p, "Pa"),
        product_t=states.product_t,
        product_h=from_si(states.product_h, "J/mol"),
        heat_duty=from_si(states.heat_duty, "W"),
        warnings=tuple(warnings),
    )


@dataclass(frozen=True, slots=True)
class ReactorStates:
    """The product's record and the duty, in SI."""

    product_n: float
    product_z: tuple[float, ...]
    product_p: float
    product_t: Q
    product_h: float
    heat_duty: float


def _route(
    components: list[str],
    feed_n: float,
    feed_z: list[float],
    feed_p: float,
    feed_t: Q,
    reaction: str,
    limiting_reactant: str,
    conversion: float,
    isothermal: bool,
    reactor_temperature: float | None,
    reactor_pressure: float | None,
    pressure_drop: float,
) -> ReactorStates:
    """The reactor's arithmetic, in SI, in the order it computes it.

    **One arithmetic, two consumers**: :func:`stirred_tank_reactor` builds the result from this
    and it is the twin of ``crates/azoth-process/src/kernels/stirred_tank_reactor.rs``. The
    refusals live here for the reason Rust puts them in the kernel.
    """
    if not 0.0 <= conversion <= 1.0:
        raise InvalidInputError(
            "conversion",
            f"a conversion is a fraction of the limiting reactant, and {conversion} is "
            "not in [0, 1]",
        )
    rows = _tables.stoichiometry(reaction)
    if not rows:
        raise InvalidInputError(
            "reaction",
            f"`{reaction}` is not a reaction the data carries, so there is no "
            "stoichiometry to apply",
        )
    limiting = next(
        ((name, coeff) for name, coeff in rows if name.lower() == limiting_reactant.lower()), None
    )
    if limiting is None:
        raise InvalidInputError(
            "limiting_reactant",
            f"`{limiting_reactant}` is not a substance of reaction `{reaction}`",
        )
    indices = {name.lower(): index for index, name in enumerate(components)}
    if limiting[0].lower() not in indices:
        raise InvalidInputError(
            "components",
            f"the limiting reactant `{limiting[0]}` is not one of the feed's substances",
        )

    # `moles(limiting) * conversion`, then every coefficient scaled by the limiting one's.
    reacted = feed_n * feed_z[indices[limiting[0].lower()]] * conversion
    scale = reacted / abs(limiting[1])
    amounts = [fraction * feed_n for fraction in feed_z]
    for name, coefficient in rows:
        index = indices.get(name.lower())
        if index is not None:
            amounts[index] += scale * coefficient
    if any(amount < 0.0 for amount in amounts):
        bad = next(name for name, amount in zip(components, amounts, strict=True) if amount < 0.0)
        raise InvalidInputError(
            "conversion",
            f"this conversion takes {bad} below zero, so the reaction would consume more of it "
            f"than the feed carries - the extent is the limiting reactant's, and the "
            f"stoichiometry is what overdraws it",
        )
    total = sum(amounts)
    if total <= 0.0:
        raise InvalidInputError(
            "conversion",
            "the reaction leaves no moles in the reactor, so there is no state to flash",
        )
    composition = [amount / total for amount in amounts]

    p_out = feed_p - pressure_drop if reactor_pressure is None else reactor_pressure
    if p_out <= 0.0:
        raise InvalidInputError(
            "pressure_drop",
            f"a pressure drop of {pressure_drop} Pa takes a {feed_p} Pa feed to {p_out} Pa, "
            f"which is not a pressure",
        )

    mixture, ideal_gas = _components.mixture_of(components, eos="pr")
    if isothermal:
        if reactor_temperature is None:
            raise InvalidInputError(
                "reactor_temperature",
                "an isothermal vessel holds a temperature, and none was stated",
            )
        pt_flash(mixture, from_si(reactor_temperature, "K"), from_si(p_out, "Pa"), composition)
        temperature = from_si(reactor_temperature, "K")
    else:
        h_in, _ = enthalpy_at(mixture, ideal_gas, feed_t.to("K").magnitude, feed_p, feed_z)
        moved = ph_flash_solve(
            mixture, ideal_gas, from_si(p_out, "Pa"), from_si(h_in, "J/mol"), composition
        )
        temperature = moved.T

    product_h = float(
        enthalpy_at(mixture, ideal_gas, temperature.to("K").magnitude, p_out, composition)[0]
    )
    duty = (
        (product_h * total - feed_n * _enthalpy_in(feed_t, feed_p, feed_z, mixture, ideal_gas))
        if isothermal
        else 0.0
    )

    return ReactorStates(
        product_n=total,
        product_z=tuple(composition),
        product_p=p_out,
        product_t=temperature,
        product_h=product_h,
        heat_duty=duty,
    )


def _enthalpy_in(
    feed_t: Q, feed_p: float, feed_z: list[float], mixture: Any, ideal_gas: Any
) -> float:
    """The feed's own molar enthalpy, on the same basis as the product's."""
    return float(enthalpy_at(mixture, ideal_gas, feed_t.to("K").magnitude, feed_p, feed_z)[0])


def _spec() -> dict[str, object]:
    """The generated spec, imported here so the module has no import-time cycle."""
    from azoth._models_gen import model

    return model("process.stirred_tank_reactor")
