"""``process.plug_flow_reactor`` - a fixed bed marched along its own length.

Spec: ``specs/models/process/plug_flow_reactor.toml``

The Python twin of ``crates/azoth-process/src/models/plug_flow_reactor.rs`` over
``crates/azoth-process/src/kernels/plug_flow_reactor.rs``, written to mirror them.

# The march

``PlugFlowReactor.run`` advances ``[F1..Fn, T, P]`` in ``numberOfSteps`` fixed steps of
``dz = length / steps`` - the classical fourth-order Runge-Kutta by default, forward Euler on
request - and **there is no step-size control of any kind**. The state's pressure is in
**bara**, because ``calculateDerivatives`` writes its own pressure row as ``-dPdz / 1e5``;
pascals appear only at the boundary.

After every step the class clamps the molar flows at zero and the pressure at ``0.1`` bar,
overrides the temperature back to the inlet when it is isothermal, and records the station.

# Two volumes, and they are different volumes

The rate law's concentration is ``x_i * getPhase(0).getDensity("mol/m3")``, which is the
**Peneloux-corrected** molar density. The superficial velocity is
``getVolume("m3") / getNumberOfMoles() * getTotalNumberOfMoles() / A``, and ``getVolume`` is
the **untranslated cubic's**. Measured on the capture's own feed state they differ by
``4.2e-4`` relative - small, and not zero, so both are carried.

# Properties are frozen by default

``thermodynamicCoupling`` is ``FROZEN_PROPERTIES``: the rate law reads the last re-flash and
never one of the RK4 sub-states, and the re-flash happens at the *top* of step ``s`` when
``s % propertyUpdateFrequency == 0 && s > 0``. Two consequences the capture measures: the
class's RK4 and Euler agree to thirteen digits on its own default, because a frozen march is
near enough to a straight line; and ``calculateResidenceTime``, called **before** the class's
final state update, reads a ninety-step temperature on a hundred-step march.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any

from azoth.core.errors import InvalidInputError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import Phase, PlugFlowReactorResult
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning
from azoth.eos import components as _databank
from azoth.eos.reference.molar_enthalpy_entropy import molar_enthalpy_entropy
from azoth.eos.reference.pr_molar_volume import pr_molar_volume
from azoth.eos.reference.pt_flash import pt_flash as pt_flash_solve
from azoth.eos.reference.viscosity import viscosity as viscosity_solve
from azoth.reactions.reference import _tables

#: The class's own gas constant, `KineticReaction.R_GAS`, J/(mol*K).
R_GAS = 8.31446

#: `updateSystemState`'s floor on a component's moles, mol/s.
MINIMUM_MOLES = 1.0e-30

#: The pressure floor the loop holds, bara.
MINIMUM_PRESSURE = 0.1

#: The temperature floor `calculateDerivatives` reads with, K.
MINIMUM_TEMPERATURE = 1.0

#: The empty tube's friction factor, where no bed is configured.
EMPTY_TUBE_FRICTION_FACTOR = 0.01

#: `CatalystBed`'s defaults, which a supplied bulk density does not change.
PARTICLE_DIAMETER = 0.003
VOID_FRACTION = 0.40
PARTICLE_POROSITY = 0.50
TORTUOSITY = 3.0

#: The class's default molecular diffusivity, m**2/s.
MOLECULAR_DIFFUSIVITY = 1.0e-5


@dataclass(frozen=True, slots=True)
class State:
    """One flashed property state, and the properties the derivatives read off it."""

    temperature: float
    pressure: float
    composition: tuple[float, ...]
    molar_mass: float
    corrected_density: float
    cubic_density: float
    cp: float
    viscosity: float

    def concentration(self, index: int) -> float:
        """A component's molar concentration, mol/m**3, on the **corrected** density."""
        return self.composition[index] * self.corrected_density / self.molar_mass

    def velocity(self, total_flow: float, total_area: float) -> float:
        """The superficial velocity, m/s, on the **untranslated cubic**."""
        molar_volume = self.molar_mass / self.cubic_density
        return total_flow * molar_volume / total_area


@dataclass(frozen=True, slots=True)
class Bed:
    """A `CatalystBed`, with the class's defaults for what a caller does not set."""

    bulk_density: float
    activity_factor: float
    particle_diameter: float
    void_fraction: float

    def pressure_drop(self, velocity: float, density: float, viscosity: float) -> float:
        """The Ergun drop per unit length, Pa/m, positive for a loss."""
        eps = self.void_fraction
        dp = self.particle_diameter
        u = abs(velocity)
        one_minus = 1.0 - eps
        eps_cubed = eps * eps * eps
        viscous = 150.0 * viscosity * one_minus * one_minus * u / (dp * dp * eps_cubed)
        inertial = 1.75 * density * one_minus * u * u / (dp * eps_cubed)
        return viscous + inertial

    def effective_diffusivity(self, molecular: float) -> float:
        """`D * eps_particle / tau`."""
        return molecular * PARTICLE_POROSITY / TORTUOSITY

    def thiele_modulus(self, rate_constant: float, diffusivity: float) -> float:
        """`(R/3) * sqrt(|k| / D_eff)`, with the class's `1e6` sentinel."""
        if diffusivity <= 0.0:
            return 1.0e6
        radius = self.particle_diameter / 2.0
        return (radius / 3.0) * math.sqrt(abs(rate_constant) / diffusivity)

    def effectiveness_factor(self, thiele: float) -> float:
        """`(1/phi) * [1/tanh(3 phi) - 1/(3 phi)]`, with the class's two cut-offs."""
        if thiele < 0.01:
            return 1.0
        if thiele > 500.0:
            return 1.0 / (3.0 * thiele)
        three = 3.0 * thiele
        return (1.0 / thiele) * (1.0 / math.tanh(three) - 1.0 / three)


@dataclass(frozen=True, slots=True)
class March:
    """The march's own record, which `layers` dumps and the result is built from."""

    components: tuple[str, ...]
    positions: tuple[float, ...]
    temperatures: tuple[float, ...]
    pressures: tuple[float, ...]
    conversions: tuple[float, ...]
    rates: tuple[float, ...]
    outlet_z: tuple[float, ...]
    outlet_n: float
    outlet_t: float
    outlet_p: float
    outlet_h: float
    conversion: float
    pressure_drop_bar: float
    residence_time: float
    heat_duty: float


@dataclass(frozen=True, slots=True)
class _Reaction:
    """The rate law's parameters, which the declaration carries."""

    stoich: tuple[tuple[str, float], ...]
    orders: dict[str, float]
    rate_type: str
    pre_exponential_factor: float
    activation_energy: float
    temperature_exponent: float
    heat_of_reaction: float
    bed: Bed | None
    effectiveness_enabled: bool
    molecular_diffusivity: float

    def rate_constant(self, temperature: float) -> float:
        """`A * T**n * exp(-Ea/(R*T))`, with the class's zero-exponent guard."""
        t_power = (
            temperature**self.temperature_exponent
            if abs(self.temperature_exponent) > 1.0e-10
            else 1.0
        )
        return (
            self.pre_exponential_factor
            * t_power
            * math.exp(-self.activation_energy / (R_GAS * temperature))
        )

    def volumetric_rate(self, state: State, species: list[str]) -> float:
        """`KineticReaction.calculateRate` in mol/(m**3 reactor * s)."""
        k = self.rate_constant(state.temperature)
        forward = k
        for name, order in self.orders.items():
            c = state.concentration(species.index(name)) if name in species else 0.0
            if c <= 0.0 and order > 0.0:
                return 0.0
            forward *= max(c, 0.0) ** order
        if self.rate_type == "lhhw":
            # **A declaration carries no adsorption terms**, so the denominator is
            # `(1 + 0) ** m`, which is one - and the class behaves the same way when
            # `addAdsorptionTerm` is never called.
            return forward
        if self.bed is not None:
            forward *= self.bed.activity_factor
            if self.effectiveness_enabled:
                effective = self.bed.effective_diffusivity(self.molecular_diffusivity)
                forward *= self.bed.effectiveness_factor(self.bed.thiele_modulus(k, effective))
        return forward


def _state_at(
    components: list[str], flows: list[float], temperature: float, pressure: float
) -> State:
    """The property state a flash at `(temperature, pressure_bar)` gives."""
    floored = [max(f, MINIMUM_MOLES) for f in flows]
    total = sum(floored)
    composition = [f / total for f in floored]
    mixture, ideal_gas = _databank.mixture_of(components, eos="pr")
    t_si, p_si = from_si(temperature, "K"), from_si(pressure * 1.0e5, "Pa")
    # **The single-phase root, which refuses a two-phase state.** `Stream::single_phase_root`
    # is the Rust side of the same rule: a stream that is two phases has no one density, so
    # nothing in the derivative below is defined and the kernel says so rather than picking one.
    flash = pt_flash_solve(mixture, t_si, p_si, composition)
    if flash.phase == Phase.TWO_PHASE:
        raise InvalidInputError(
            "stream",
            f"this state is two phases at {temperature} K and {pressure} bar, so it has no one "
            "density or viscosity; a caller has to say which phase's they mean",
        )
    z_factor = flash.z_vapour if flash.phase == Phase.ALL_VAPOUR else flash.z_liquid
    masses = [
        zi * component.molar_mass.to("kg/mol").magnitude  # type: ignore[union-attr]
        for zi, component in zip(composition, mixture.components, strict=True)
    ]
    molar_mass = sum(masses)
    cubic_volume = pr_molar_volume(z_factor, t_si, p_si).v.to("m**3/mol").magnitude
    shift = mixture.volume_shift(composition)
    props = molar_enthalpy_entropy(mixture, ideal_gas, t_si, p_si, composition, z_factor)
    mu = viscosity_solve(mixture, t_si, p_si, composition).mu.to("Pa*s").magnitude
    return State(
        temperature=temperature,
        pressure=pressure,
        composition=tuple(composition),
        molar_mass=molar_mass,
        corrected_density=molar_mass / (cubic_volume - shift),
        cubic_density=molar_mass / cubic_volume,
        cp=props.cp.to("J/(mol*K)").magnitude,
        viscosity=mu,
    )


def plug_flow_reactor(
    components: list[str],
    feed_n: Q,
    feed_z: list[float],
    feed_p: Q,
    feed_t: Q,
    length: Q,
    diameter: Q,
    number_of_tubes: float,
    energy_mode: str,
    coolant_temperature: Q,
    overall_heat_transfer_coefficient: Q,
    number_of_steps: float,
    integration_method: str,
    property_update_frequency: float,
    thermodynamic_coupling: str,
    reaction: str,
    reaction_orders: list[float],
    rate_type: str,
    pre_exponential_factor: float,
    activation_energy: Q,
    temperature_exponent: float,
    heat_of_reaction: Q,
    catalyst_bulk_density: Q | None = None,
    catalyst_activity_factor: float | None = None,
    catalyst_particle_diameter: Q | None = None,
    catalyst_void_fraction: float | None = None,
    catalyst_molecular_diffusivity: Q | None = None,
    catalyst_effectiveness_enabled: bool | None = None,
    key_component: str | None = None,
) -> PlugFlowReactorResult:
    """March a feed along a tube, reacting as it goes.

    See :func:`azoth.process.plug_flow_reactor` for what the machine is; this is the
    arithmetic.
    """
    warnings: list[Warning] = []
    spec = _spec()
    checks = checks_for(spec)
    # **A unit-carrying input arrives as a quantity and a dimensionless one as a bare float**,
    # which is the split `input_to_si` refuses to guess across - so the two are read apart and
    # everything below this line is SI magnitudes.
    length_si = input_to_si(spec, "length", length)
    diameter_si = input_to_si(spec, "diameter", diameter)
    activation_si = input_to_si(spec, "activation_energy", activation_energy)
    heat_of_reaction_si = input_to_si(spec, "heat_of_reaction", heat_of_reaction)
    particle_si = (
        None
        if catalyst_particle_diameter is None
        else input_to_si(spec, "catalyst_particle_diameter", catalyst_particle_diameter)
    )
    feed_n_si = input_to_si(spec, "feed_n", feed_n)
    feed_p_si = input_to_si(spec, "feed_p", feed_p)
    feed_t_si = input_to_si(spec, "feed_t", feed_t)
    apply_checks(
        checks.on_input,
        {
            "length": length_si,
            "diameter": diameter_si,
            "number_of_steps": number_of_steps,
            "property_update_frequency": property_update_frequency,
            "catalyst_activity_factor": catalyst_activity_factor,
        }.get,
        warnings,
    )

    if energy_mode not in ("adiabatic", "isothermal", "coolant"):
        raise InvalidInputError(
            "energy_mode",
            f"`{energy_mode}` is not one of adiabatic, isothermal or coolant",
        )
    if thermodynamic_coupling not in ("frozen_properties", "fully_coupled"):
        raise InvalidInputError(
            "thermodynamic_coupling",
            f"`{thermodynamic_coupling}` is not one of frozen_properties or fully_coupled",
        )
    if rate_type not in ("power_law", "lhhw", "equilibrium"):
        raise InvalidInputError(
            "rate_type", f"`{rate_type}` is not one of power_law, lhhw or equilibrium"
        )

    rows = _tables.stoichiometry(reaction)
    if not rows:
        raise InvalidInputError("reaction", f"`{reaction}` is not a reaction the data carries")
    reactants = sum(1 for _, coefficient in rows if coefficient < 0.0)
    if len(reaction_orders) != reactants:
        raise InvalidInputError(
            "reaction_orders",
            f"`{reaction}` names {reactants} reactants, so it needs {reactants} orders, "
            f"but {len(reaction_orders)} were given",
        )

    orders: dict[str, float] = {}
    next_order = iter(reaction_orders)
    for name, coefficient in rows:
        if coefficient < 0.0:
            orders[name] = next(next_order, 0.0)

    bed = (
        None
        if catalyst_bulk_density is None
        else Bed(
            bulk_density=input_to_si(spec, "catalyst_bulk_density", catalyst_bulk_density),
            activity_factor=catalyst_activity_factor
            if catalyst_activity_factor is not None
            else 1.0,
            particle_diameter=particle_si if particle_si is not None else PARTICLE_DIAMETER,
            void_fraction=catalyst_void_fraction
            if catalyst_void_fraction is not None
            else VOID_FRACTION,
        )
    )
    kinetics = _Reaction(
        stoich=tuple(rows),
        orders=orders,
        rate_type=rate_type,
        pre_exponential_factor=pre_exponential_factor,
        activation_energy=activation_si,
        temperature_exponent=temperature_exponent,
        heat_of_reaction=heat_of_reaction_si,
        bed=bed,
        effectiveness_enabled=bool(catalyst_effectiveness_enabled),
        molecular_diffusivity=(
            MOLECULAR_DIFFUSIVITY
            if catalyst_molecular_diffusivity is None
            else input_to_si(spec, "catalyst_molecular_diffusivity", catalyst_molecular_diffusivity)
        ),
    )

    march = _march(
        components,
        feed_n_si,
        feed_z,
        # **The march carries bara**, as `calculateDerivatives` writes its own pressure row,
        # so the pascals stop here and only the two boundary fields convert back.
        feed_p_si / 1.0e5,
        feed_t_si,
        length=length_si,
        diameter=diameter_si,
        number_of_tubes=max(1, int(number_of_tubes)),
        energy_mode=energy_mode,
        coolant_temperature=input_to_si(spec, "coolant_temperature", coolant_temperature),
        overall_heat_transfer_coefficient=input_to_si(
            spec, "overall_heat_transfer_coefficient", overall_heat_transfer_coefficient
        ),
        number_of_steps=max(1, int(number_of_steps)),
        integration_method=integration_method,
        update_every=max(1, int(property_update_frequency)),
        fully_coupled=thermodynamic_coupling == "fully_coupled",
        kinetics=kinetics,
        key_component=key_component,
    )

    return PlugFlowReactorResult(
        product_n=from_si(march.outlet_n, "mol/s"),
        product_z=march.outlet_z,
        product_p=from_si(march.outlet_p, "Pa"),
        product_t=from_si(march.outlet_t, "K"),
        product_h=from_si(march.outlet_h, "J/mol"),
        conversion=march.conversion,
        pressure_drop=from_si(march.pressure_drop_bar * 1.0e5, "Pa"),
        outlet_temperature=from_si(march.outlet_t, "K"),
        heat_duty=from_si(march.heat_duty, "W"),
        positions=march.positions,
        temperature_profile=march.temperatures,
        pressure_profile=tuple(p * 1.0e5 for p in march.pressures),
        conversion_profile=march.conversions,
        warnings=tuple(warnings),
    )


def _march(
    feed_components: list[str],
    feed_n: float,
    feed_z: list[float],
    feed_p_bar: float,
    feed_t: float,
    *,
    length: float,
    diameter: float,
    number_of_tubes: int,
    energy_mode: str,
    coolant_temperature: float,
    overall_heat_transfer_coefficient: float,
    number_of_steps: int,
    integration_method: str,
    update_every: int,
    fully_coupled: bool,
    kinetics: _Reaction,
    key_component: str | None,
) -> March:
    """The whole march: the kernel's ``plug_flow_reactor``."""
    # `ensureProductComponentsExist`: species a reaction names and the feed does not carry are
    # appended at 1e-20 mol, in the reaction's own stoichiometry order.
    species = list(feed_components)
    flows = [z * feed_n for z in feed_z]
    for name, _ in kinetics.stoich:
        if name not in species:
            species.append(name)
            flows.append(1.0e-20)
    n = len(species)

    if key_component is not None:
        if key_component not in species:
            raise InvalidInputError(
                "key_component", f"`{key_component}` is not in the reactor's species set"
            )
        key_index = species.index(key_component)
    else:
        reactant = next((name for name, c in kinetics.stoich if c < 0.0), None)
        if reactant is None:
            raise InvalidInputError(
                "key_component",
                "no reaction names a reactant, so no conversion can be reported",
            )
        key_index = species.index(reactant)
    initial_key = flows[key_index]

    total_area = math.pi * diameter * diameter / 4.0 * number_of_tubes
    perimeter = math.pi * diameter
    dz = length / number_of_steps

    temperature, pressure = feed_t, feed_p_bar
    inlet_temperature = temperature
    state = [*flows, temperature, pressure]
    properties = _state_at(species, flows, temperature, pressure)

    positions, temperatures, pressures = [0.0], [temperature], [pressure]
    conversions, rates = [0.0], [0.0]

    def derivatives(trial: list[float], current: State) -> list[float]:
        return _derivatives(
            trial,
            species,
            current,
            kinetics,
            total_area=total_area,
            perimeter=perimeter,
            number_of_tubes=number_of_tubes,
            diameter=diameter,
            energy_mode=energy_mode,
            coolant_temperature=coolant_temperature,
            overall_heat_transfer_coefficient=overall_heat_transfer_coefficient,
            n=n,
        )

    for step in range(1, number_of_steps + 1):
        if integration_method == "EULER":
            k1 = derivatives(state, properties)
            for i in range(len(state)):
                state[i] += k1[i] * dz
        else:
            k1 = derivatives(state, properties)
            s2 = [state[i] + 0.5 * dz * k1[i] for i in range(len(state))]
            k2 = derivatives(s2, properties)
            s3 = [state[i] + 0.5 * dz * k2[i] for i in range(len(state))]
            k3 = derivatives(s3, properties)
            s4 = [state[i] + dz * k3[i] for i in range(len(state))]
            k4 = derivatives(s4, properties)
            for i in range(len(state)):
                state[i] += dz / 6.0 * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i])

        for i in range(n):
            state[i] = max(state[i], 0.0)
        state[n + 1] = max(state[n + 1], MINIMUM_PRESSURE)
        if energy_mode == "isothermal":
            state[n] = inlet_temperature
        temperature, pressure = state[n], state[n + 1]

        # **The class's refresh schedule, including where it does not fire.** At the default
        # frequency of ten a hundred-step march last re-flashes after step ninety, so the
        # `state` at the end of the loop is *not* what the properties hold - which is what
        # makes `calculateResidenceTime`, called before the class's final update, read a
        # ninety-step temperature.
        if not fully_coupled and step % update_every == 0 and step < number_of_steps:
            properties = _state_at(species, state[:n], temperature, pressure)

        conversions.append(
            min(1.0, max(0.0, 1.0 - state[key_index] / initial_key))
            if initial_key > 1.0e-30
            else 0.0
        )
        rates.append(abs(kinetics.volumetric_rate(properties, species)))
        positions.append(step * dz)
        temperatures.append(temperature)
        pressures.append(pressure)

    final = state[:n]
    total = sum(final)
    volumetric = properties.velocity(total, total_area) * total_area
    conversion = (
        min(1.0, max(0.0, 1.0 - state[key_index] / initial_key)) if initial_key > 1.0e-30 else 0.0
    )
    duty = 0.0
    if energy_mode == "isothermal":
        for name, coefficient in kinetics.stoich:
            if coefficient < 0.0 and name in feed_components:
                index = species.index(name)
                initial = feed_z[feed_components.index(name)] * feed_n
                duty += -kinetics.heat_of_reaction * (initial - final[index]) / abs(coefficient)
                break
    outlet = _state_at(species, final, temperature, pressure)
    mixture, ideal_gas = _databank.mixture_of(species, eos="pr")
    t_si, p_si = from_si(temperature, "K"), from_si(pressure * 1.0e5, "Pa")
    composition = list(outlet.composition)
    # **The single-phase root, which refuses a two-phase state.** `Stream::single_phase_root`
    # is the Rust side of the same rule: a stream that is two phases has no one density, so
    # nothing in the derivative below is defined and the kernel says so rather than picking one.
    flash = pt_flash_solve(mixture, t_si, p_si, composition)
    if flash.phase == Phase.TWO_PHASE:
        raise InvalidInputError(
            "stream",
            f"this state is two phases at {temperature} K and {pressure} bar, so it has no one "
            "density or viscosity; a caller has to say which phase's they mean",
        )
    z_factor = flash.z_vapour if flash.phase == Phase.ALL_VAPOUR else flash.z_liquid
    outlet_h = (
        molar_enthalpy_entropy(mixture, ideal_gas, t_si, p_si, composition, z_factor)
        .h.to("J/mol")
        .magnitude
    )

    return March(
        components=tuple(species),
        positions=tuple(positions),
        temperatures=tuple(temperatures),
        pressures=tuple(pressures),
        conversions=tuple(conversions),
        rates=tuple(rates),
        outlet_z=tuple(f / total for f in final),
        outlet_n=total,
        outlet_t=temperature,
        outlet_p=pressure * 1.0e5,
        outlet_h=outlet_h,
        conversion=conversion,
        pressure_drop_bar=feed_p_bar - pressure,
        residence_time=total_area * length / volumetric if volumetric > 0.0 else 0.0,
        heat_duty=duty,
    )


def states(
    feed_components: list[str],
    feed_n: float,
    feed_z: list[float],
    feed_p_bar: float,
    feed_t: float,
    **kwargs: Any,
) -> March:
    """The march's intermediates, for the layer diff.

    Factored out so :func:`plug_flow_reactor` and the dumper read one arithmetic - the recipe
    `layers.py` states, and the reason the profile is checkable station by station.
    """
    return _march(feed_components, feed_n, feed_z, feed_p_bar, feed_t, **kwargs)


def _derivatives(
    trial: list[float],
    species: list[str],
    properties: State,
    kinetics: _Reaction,
    *,
    total_area: float,
    perimeter: float,
    number_of_tubes: int,
    diameter: float,
    energy_mode: str,
    coolant_temperature: float,
    overall_heat_transfer_coefficient: float,
    n: int,
) -> list[float]:
    """`calculateDerivatives`: `[dF1/dz..dFn/dz, dT/dz, dP/dz]`."""
    derivs = [0.0] * (n + 2)
    temperature = max(trial[n], MINIMUM_TEMPERATURE)
    total_mol_flow = sum(max(f, 0.0) for f in trial[:n])
    if total_mol_flow < MINIMUM_MOLES:
        return derivs

    volumetric = kinetics.volumetric_rate(properties, species)
    for name, coefficient in kinetics.stoich:
        derivs[species.index(name)] += total_area * coefficient * volumetric
    heat_generation = volumetric * kinetics.heat_of_reaction

    capacity_flow = properties.cp * total_mol_flow
    if energy_mode == "adiabatic":
        if abs(capacity_flow) > 1.0e-20:
            derivs[n] = -heat_generation * total_area / capacity_flow
    elif energy_mode == "coolant":
        transfer = (
            overall_heat_transfer_coefficient
            * perimeter
            * number_of_tubes
            * (coolant_temperature - temperature)
        )
        if abs(capacity_flow) > 1.0e-20:
            derivs[n] = (-heat_generation * total_area + transfer) / capacity_flow

    velocity = properties.velocity(total_mol_flow, total_area)
    if kinetics.bed is not None:
        drop = kinetics.bed.pressure_drop(velocity, properties.cubic_density, properties.viscosity)
    else:
        drop = (
            EMPTY_TUBE_FRICTION_FACTOR
            * properties.cubic_density
            * velocity
            * velocity
            / (2.0 * diameter)
        )
    derivs[n + 1] = -drop / 1.0e5
    return derivs


def _spec() -> dict[str, Any]:
    """The model spec's checks, which both implementations apply."""
    from azoth._models_gen import model

    return model("process.plug_flow_reactor")
