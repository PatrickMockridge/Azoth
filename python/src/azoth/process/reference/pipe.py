"""``process.pipe`` - a line's pressure drop, with its outlet pressure solved.

Spec: ``specs/models/process/pipe.toml``

The Python twin of ``crates/azoth-process/src/models/pipe.rs`` over
``crates/azoth-process/src/kernels/pipe.rs``, written to mirror them.

# Two equations, one per phase

``AdiabaticPipe.calcPressureOut`` branches on the phase's type. A **gas** takes the
compressible ``P1^2 - P2^2`` form - ``Z``, molar mass and temperature, and no density at
all - and anything else takes Darcy-Weisbach. The branch is the phase label, which is
NeqSim's ``PhaseEos.init`` rule and which this library ports for the hydrate dose pair.

# Two densities, which is the quirk

The liquid branch's *velocity* takes ``getPhysicalProperties().getDensity()`` - the volume
translation applied - and its *Reynolds number* takes ``getKinematicViscosity()``, which
divides the same viscosity by the **untranslated** cubic ``M / (Z R T / P)``. Measured on the
captured butane row the two are ``567.33`` and ``603.86`` kg/m3, 6% apart, and a port that
used one density for both is 6% out on ``Re``.

# The viscosity follows the phase type

``getPhysicalProperties()`` dispatches on it: a gas or a hydrocarbon liquid takes
``PFCTViscosityMethodHeavyOil`` - ``eos.viscosity`` - and an **aqueous** phase takes
``WaterPhysicalProperties``, whose correlation is the liquid ``Viscosity`` class
(``eos.aqueous_viscosity``). Water at 300 K is ``8.5510e-4`` Pa s through the second and
``5.3097e-4`` through the first, so taking one for both puts a water line's Reynolds number
``1.61`` out - which is what this row did until that id existed.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from azoth.core.errors import InvalidInputError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import Phase, PipeResult
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning
from azoth.eos import components as _databank
from azoth.eos.mixture import Mixture
from azoth.eos.reference._mixture_state import reduced_parameters
from azoth.eos.reference.aqueous_viscosity import aqueous_viscosity as aqueous_viscosity_solve
from azoth.eos.reference.hydrate_inhibitor_wt import AQUEOUS, GAS, phase_label
from azoth.eos.reference.ph_flash import enthalpy_at
from azoth.eos.reference.pr_molar_volume import pr_molar_volume
from azoth.eos.reference.pt_flash import pt_flash as pt_flash_solve
from azoth.eos.reference.viscosity import viscosity as viscosity_solve

#: The gas constant, as `eos` uses it everywhere.
R = 8.31446261815324

#: `AdiabaticPipe.run`'s own convergence: `1e-2` bar, and `iter < 25`.
PRESSURE_TOLERANCE = 1.0e3
MAX_PASSES = 25


def pipe(
    components: list[str],
    inlet_n: Q,
    inlet_z: list[float],
    inlet_p: Q,
    inlet_t: Q,
    length: Q,
    diameter: Q,
    roughness: Q,
) -> PipeResult:
    """Drop a stream's pressure along a line.

    Args:
        components: the fluid's substances, by name.
        inlet_n: the inlet's molar flow.
        inlet_z: the inlet composition.
        inlet_p: the inlet pressure.
        inlet_t: the inlet temperature, which is also the outlet's.
        length: the pipe's length.
        diameter: the internal diameter.
        roughness: the sand-grain roughness height.

    Returns:
        The outlet's record and the pressure drop the line implies.

    Raises:
        InvalidInputError: where the geometry is not positive.
        OutOfRangeError: for a temperature the spec bounds.

    Example:
        >>> import azoth
        >>> q = azoth.ureg.Quantity
        >>> r = pipe(
        ...     ["methane", "CO2"],
        ...     q(1.0, "mol/s"),
        ...     [0.7, 0.3],
        ...     q(50.0, "bar"),
        ...     q(300.0, "K"),
        ...     length=q(1000.0, "m"),
        ...     diameter=q(0.1, "m"),
        ...     roughness=q(1e-5, "m"),
        ... )
        >>> round(r.pressure_drop.to("Pa").magnitude, 4)
        21.7722
    """
    spec = _spec()
    checks = checks_for(spec)
    warnings: list[Warning] = []

    n = input_to_si(spec, "inlet_n", inlet_n)
    t = input_to_si(spec, "inlet_t", inlet_t)
    p = input_to_si(spec, "inlet_p", inlet_p)
    length_si = input_to_si(spec, "length", length)
    diameter_si = input_to_si(spec, "diameter", diameter)
    roughness_si = input_to_si(spec, "roughness", roughness)

    apply_checks(
        checks.on_input,
        {
            "inlet_t": t,
            "length": length_si,
            "diameter": diameter_si,
            "roughness": roughness_si,
        }.get,
        warnings,
    )

    states = _route(components, n, inlet_z, p, t, length_si, diameter_si, roughness_si)

    return PipeResult(
        outlet_n=from_si(n, "mol/s"),
        outlet_z=tuple(inlet_z),
        outlet_p=from_si(states.p_out, "Pa"),
        outlet_t=from_si(t, "K"),
        outlet_h=from_si(states.h_out, "J/mol"),
        pressure_drop=from_si(p - states.p_out, "Pa"),
        warnings=tuple(warnings),
    )


@dataclass(frozen=True, slots=True)
class PipeStates:
    """The line's intermediates, in the order it computes them, in SI."""

    #: The outlet pressure the loop settled on.
    p_out: float
    #: The outlet molar enthalpy.
    h_out: float
    #: The velocity the drop was computed at, m/s.
    velocity: float
    #: The Reynolds number it was computed at.
    reynolds: float
    #: The Darcy friction factor.
    friction_factor: float
    #: `laminar`, `transition` or `turbulent`.
    regime: str
    #: The passes the outlet pressure took to settle.
    passes: int


@dataclass(frozen=True, slots=True)
class _Pass:
    """One pass's state, before the loop's pressure has been settled."""

    is_gas: bool
    z_factor: float
    molar_mass: float
    temperature: float
    moles: float
    density: float
    velocity: float
    reynolds: float
    friction_factor: float
    regime: str


def _route(
    components: list[str],
    inlet_n: float,
    inlet_z: list[float],
    inlet_p: float,
    inlet_t: float,
    length: float,
    diameter: float,
    roughness: float,
) -> PipeStates:
    """The line's intermediates, in SI, in the order it computes them."""
    for name, value in (
        ("length", length),
        ("diameter", diameter),
        ("roughness", roughness),
    ):
        if not value > 0.0:
            raise InvalidInputError(
                name, f"a pipe's {name} is a positive length, and {value} is not"
            )

    mixture, ideal_gas = _databank.mixture_of(components, eos="pr")
    area = math.pi / 4.0 * diameter * diameter
    relative_roughness = roughness / diameter

    pressure = inlet_p
    passes = 0
    state = _pass(mixture, inlet_n, inlet_z, inlet_t, pressure, area, diameter, relative_roughness)
    while True:
        passes += 1
        nxt = _outlet_pressure(state, inlet_p, length, diameter)
        settled = abs(nxt - pressure) <= PRESSURE_TOLERANCE
        pressure = nxt
        if settled or passes >= MAX_PASSES:
            break
        state = _pass(
            mixture, inlet_n, inlet_z, inlet_t, pressure, area, diameter, relative_roughness
        )

    # The outlet is a flash at the pressure the loop settled on, at the inlet's temperature:
    # `run` ends with a `TPflash` over the system it has been moving.
    h_out, _ = enthalpy_at(mixture, ideal_gas, inlet_t, pressure, inlet_z)

    return PipeStates(
        p_out=pressure,
        h_out=h_out,
        velocity=state.velocity,
        reynolds=state.reynolds,
        friction_factor=state.friction_factor,
        regime=state.regime,
        passes=passes,
    )


def _pass(
    mixture: Mixture,
    inlet_n: float,
    inlet_z: list[float],
    inlet_t: float,
    pressure: float,
    area: float,
    diameter: float,
    relative_roughness: float,
) -> _Pass:
    """Evaluate the line at one pressure."""
    flash = pt_flash_solve(mixture, from_si(inlet_t, "K"), from_si(pressure, "Pa"), inlet_z)
    # Phase 0 of the flashed system: the liquid when there are two, the single phase when
    # there is one. `pt_flash` reports `x` for the liquid root and `y` for the vapour one.
    if flash.phase == Phase.TWO_PHASE:
        z_factor, composition = flash.z_liquid, list(flash.x)
    elif flash.phase == Phase.ALL_VAPOUR:
        z_factor, composition = flash.z_vapour, list(flash.y)
    else:
        z_factor, composition = flash.z_liquid, list(flash.x)

    reduced = reduced_parameters(mixture, inlet_t, pressure)
    label = phase_label(mixture, reduced, composition, z_factor)
    is_gas = label == GAS

    # In kg/mol, as a bare float: the rest of this arithmetic is SI magnitudes, and a
    # quantity that leaked into it would surface as a comparison against a float.
    masses = [
        zi * component.molar_mass.to("kg/mol").magnitude  # type: ignore[union-attr]
        for zi, component in zip(inlet_z, mixture.components, strict=True)
    ]
    molar_mass = sum(masses)
    # **Two densities.** The translated one is what the liquid branch's velocity takes; the
    # cubic one is what `getKinematicViscosity` divides by. See the module doc.
    volume = pr_molar_volume(z_factor, from_si(inlet_t, "K"), from_si(pressure, "Pa")).v
    density = molar_mass / (volume.to("m**3/mol").magnitude - mixture.volume_shift(inlet_z))
    cubic_density = molar_mass / (z_factor * R * inlet_t / pressure)

    # **The viscosity is the phase's, and NeqSim dispatches on the phase type.** A gas or a
    # hydrocarbon liquid takes the PFCT form; an **aqueous** phase takes
    # `WaterPhysicalProperties` and the liquid `Viscosity` correlation. Water at 300 K is
    # `8.5510e-4` Pa s through the second and `5.3097e-4` through the first.
    if label == AQUEOUS:
        mu = (
            aqueous_viscosity_solve(
                mixture, from_si(inlet_t, "K"), from_si(pressure, "Pa"), composition
            )
            .viscosity.to("Pa*s")
            .magnitude
        )
    else:
        mu = (
            viscosity_solve(mixture, from_si(inlet_t, "K"), from_si(pressure, "Pa"), composition)
            .mu.to("Pa*s")
            .magnitude
        )
    kinematic = mu / cubic_density

    if is_gas:
        volume_rate = inlet_n * z_factor * R * inlet_t / pressure
        velocity = volume_rate / area
    else:
        velocity = inlet_n * molar_mass / density / area

    reynolds = velocity * diameter / kinematic
    friction, regime = _friction_factor(reynolds, relative_roughness)

    return _Pass(
        is_gas=is_gas,
        z_factor=z_factor,
        molar_mass=molar_mass,
        temperature=inlet_t,
        moles=inlet_n,
        density=density,
        velocity=velocity,
        reynolds=reynolds,
        friction_factor=friction,
        regime=regime,
    )


def _outlet_pressure(state: _Pass, inlet_pressure: float, length: float, diameter: float) -> float:
    """The outlet pressure a pass implies, in Pa."""
    if state.is_gas:
        # The compressible general flow equation, `P1^2 - P2^2`, exactly as `run` writes it -
        # including that its leading term is the *mass* flow, `4 n M / pi`.
        leading = (4.0 * state.moles * state.molar_mass / math.pi) ** 2
        dp = (
            leading
            * state.friction_factor
            * length
            * state.z_factor
            * R
            / state.molar_mass
            * state.temperature
            / diameter**5
        )
        return math.sqrt(max(0.0, inlet_pressure**2 - dp))
    # No gravity term: the palette declares no elevations, and `run`'s `dp_gravity` is
    # `rho g (z_in - z_out)` - zero for a level line.
    dp = (
        state.friction_factor
        * length
        / diameter
        * state.density
        * state.velocity
        * state.velocity
        / 2.0
    )
    return inlet_pressure - dp


def _friction_factor(reynolds: float, relative_roughness: float) -> tuple[float, str]:
    """The Darcy friction factor and its regime, as `calcWallFrictionFactor` states them."""
    re = abs(reynolds)
    if re < 1e-10:
        return 0.0, "laminar"
    if re < 2300.0:
        return 64.0 / re, "laminar"
    if re < 4000.0:
        laminar = 64.0 / 2300.0
        turbulent = _haaland(4000.0, relative_roughness)
        return laminar + (turbulent - laminar) * (re - 2300.0) / 1700.0, "transition"
    return _haaland(re, relative_roughness), "turbulent"


def _haaland(reynolds: float, relative_roughness: float) -> float:
    """Haaland's explicit approximation to Colebrook, as the class writes it."""
    inner = 6.9 / reynolds + (relative_roughness / 3.7) ** 1.11
    return (1.0 / (-1.8 * math.log10(inner))) ** 2


def _spec() -> dict[str, object]:
    """The generated spec, imported here so the module has no import-time cycle."""
    from azoth._models_gen import model

    return model("process.pipe")
