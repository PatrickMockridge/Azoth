"""``process.gas_scrubber`` - the separator's kernel under the other entry, and the one
number the class adds to it.

Spec: ``specs/models/process/gas_scrubber.toml``

The Python twin of ``crates/azoth-process/src/models/gas_scrubber.rs`` over
``crates/azoth-process/src/kernels/gas_scrubber.rs``, written to mirror them: the steady state
is the separator's flash, and ``getCapacityUtilization`` is the whole of the difference.

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

# The capacity metric reads the *entrained* system

``getCapacityUtilization`` is ``Q / (K sqrt((rho_l - rho_g)/rho_g) A)`` with ``A = pi d**2/4``.
``Separator.run`` applies the entrainment to the system and only then does the getter read it, so
``Q`` is the vapour that is *left* and ``rho_l`` is the oil that has absorbed the carried moles -
measured on the capture's pair, ``0.0033609430039118087`` against the dry ``0.003341057093792472``.
Two rules of the class's own are worth naming because the columns' capacity families do neither:
``rho_l`` is the ``oil`` phase's density, ``aqueous``'s where there is no oil and ``1000.0`` where
there is neither (no ``10.0`` floor), and a stream with no gas phase answers ``0.0``. NeqSim
answers ``NaN`` when the liquid is no denser than the gas; this refuses that state rather than
carrying it.
"""

from __future__ import annotations

import math

from azoth.core.errors import InvalidInputError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import GasScrubberResult, Phase
from azoth.core.units import Q, from_si, input_to_si, quantity
from azoth.core.warnings import Warning
from azoth.eos import components as _components
from azoth.eos.reference.hydrate_inhibitor_wt import GAS
from azoth.eos.reference.ph_flash import enthalpy_at
from azoth.eos.reference.ph_flash import ph_flash as ph_flash_solve
from azoth.eos.reference.pt_flash import pt_flash
from azoth.process.reference._phase_outlet import phase_enthalpy
from azoth.process.reference.capacity import _labelled_phases, _mass_density, _molar_mass

#: ``Separator.getMaxAllowableGasVelocity``'s own substitution where the system has no liquid
#: phase at all - a bare ``1000.0`` beside the ``50.0`` it puts on an absent gas phase, which a
#: stream with no gas never reaches because ``getCapacityUtilization`` answers zero first.
DEFAULT_LIQUID_DENSITY = 1000.0


def gas_scrubber(
    components: list[str],
    feed_n: Q,
    feed_z: list[float],
    feed_p: Q,
    feed_t: Q,
    pressure_drop: Q,
    gas_in_liquid: float,
    internal_diameter: Q | None = None,
    design_gas_load_factor: Q | None = None,
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
        gas_in_liquid: the fraction of the vapour's moles carried into the liquid outlet.
        internal_diameter: the vessel's internal diameter, or ``None`` for no capacity metric.
        design_gas_load_factor: the Souders-Brown ``K`` the vessel is designed to, or ``None``.
        heat_input: heat added before the flash, or ``None`` for a flash at the feed's
            temperature.

    Returns:
        Both outlets' records.

    Raises:
        InvalidInputError: where the pressure drop takes the outlet below zero, the
            entrainment fraction is outside ``[0, 1]``, only one of the two mechanical
            inputs is stated, or the state makes the capacity metric undefined.
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
        {
            "pressure_drop": drop,
            "gas_in_liquid": gas_in_liquid,
            "internal_diameter": (
                None
                if internal_diameter is None
                else input_to_si(spec, "internal_diameter", internal_diameter)
            ),
            "design_gas_load_factor": (
                None
                if design_gas_load_factor is None
                else input_to_si(spec, "design_gas_load_factor", design_gas_load_factor)
            ),
        }.get,
        warnings,
    )

    if not 0.0 <= gas_in_liquid <= 1.0:
        raise InvalidInputError(
            "gas_in_liquid",
            f"an entrainment fraction is a fraction, and {gas_in_liquid} is not in [0, 1]",
        )
    # **Both or neither.** `getCapacityUtilization` is a function of the vessel, and a caller
    # who states half of the vessel has stated nothing: the metric is not computed at all, and
    # the half that was given is refused rather than quietly ignored.
    if (internal_diameter is None) != (design_gas_load_factor is None):
        named, missing = (
            ("internal_diameter", "design_gas_load_factor")
            if internal_diameter is not None
            else ("design_gas_load_factor", "internal_diameter")
        )
        raise InvalidInputError(
            missing,
            f"`{named}` was stated and `{missing}` was not, and the Souders-Brown capacity "
            f"metric needs both - a vessel with a size but no design basis, or the reverse, "
            f"has nothing to check",
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
        beta = _vapour_fraction(flash.phase, flash.vapour_fraction)
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
        beta = _vapour_fraction(moved.phase, moved.vapour_fraction)

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

    # The class reads the system `run` left, so the metric is at the *outlet* state - the feed's
    # temperature unless a `heat_input` moved the flash.
    capacity = (
        None
        if internal_diameter is None or design_gas_load_factor is None
        else _capacity_utilization(
            components,
            t_out,
            p_out,
            list(feed_z),
            n_vapour - carried,
            list(liquid_z),
            input_to_si(spec, "internal_diameter", internal_diameter),
            input_to_si(spec, "design_gas_load_factor", design_gas_load_factor),
        )
    )

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
        capacity_utilization=capacity,
        warnings=tuple(warnings),
    )


def _capacity_utilization(
    components: list[str],
    t: float,
    p: float,
    z: list[float],
    gas_n: float,
    liquid_z: list[float],
    diameter_m: float,
    design_gas_load_factor: float,
) -> float:
    """``Separator.getCapacityUtilization``, on the system the run left behind.

    ``z`` is the feed the flash started from, because the phase *types* the class asks for are
    the ones that flash assigned - ``hasPhaseType`` reads a label that ``addPhaseFractionToPhase``
    does not revisit - while ``liquid_z`` is what the oil phase *holds* now, because
    ``initPhysicalProperties`` does recompute its density. The two are the same on a run with no
    entrainment, which is what the capture's dry row fixes.
    """
    mixture, labelled, _shares = _labelled_phases(components, t, p, z)
    gas = next((entry for entry in labelled if entry[0] == GAS), None)
    if gas is None:
        # The class's own first answer: a stream with no gas phase uses none of its capacity.
        return 0.0
    _gas_label, gas_z, gas_root = gas
    gas_density = _mass_density(mixture, t, p, gas_z, gas_root)
    gas_molar_mass = _molar_mass(mixture, gas_z)

    if any(entry[0] != GAS for entry in labelled):
        # The oil's density is the *entrained* composition's, at the root that composition's own
        # flash puts it on - `Separator.run` adds the carried moles before the getter reads.
        entrained = pt_flash(mixture, quantity(t, "K"), quantity(p, "Pa"), list(liquid_z))
        liquid_density = _mass_density(mixture, t, p, list(liquid_z), entrained.z_liquid)
    else:
        liquid_density = DEFAULT_LIQUID_DENSITY

    if math.isnan(liquid_density) or math.isnan(gas_density) or liquid_density <= gas_density:
        raise InvalidInputError(
            "capacity_utilization",
            f"the liquid phase is {liquid_density} kg/m3 against the gas's {gas_density}, so "
            f"`rho_l - rho_g` is not positive and the Souders-Brown velocity is the square root "
            f"of a negative number. NeqSim's `getCapacityUtilization` answers NaN there; this "
            f"refuses the state rather than answering it",
        )

    area = math.pi * diameter_m * diameter_m / 4.0
    velocity = design_gas_load_factor * math.sqrt((liquid_density - gas_density) / gas_density)
    volumetric_flow = gas_n * gas_molar_mass / gas_density
    return volumetric_flow / (velocity * area)


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
