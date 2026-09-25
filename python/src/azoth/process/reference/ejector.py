"""``process.ejector`` - the ejector's kernel.

Spec: ``specs/models/process/ejector.toml``

The Python twin of ``crates/azoth-process/src/models/ejector.rs`` over
``crates/azoth-process/src/kernels/ejector.rs``, written to mirror them.

# The machine in one paragraph

``run`` is a quasi one-dimensional route. Each stream is expanded **isentropically** to the
mixing pressure, then dropped to the enthalpy its own nozzle efficiency leaves, and given the
velocity ``sqrt(2 dh)`` that drop implies; the two momenta mix with the mixing chamber's
efficiency on the ideal velocity; the static enthalpy that leaves is flashed at the mixing
pressure; the diffuser recovers ``eta * v^2/2`` up to the discharge pressure; and the
diffuser's own design velocity takes one more ``v^2/2`` off before the final flash.

# What the class estimates rather than being told

**The mixing pressure.** ``setMixingPressure`` exists and the palette does not declare it, so
``estimateDefaultMixingPressure`` decides - and it is written **in bara**, down to a ``0.01``
floor that is a bar and not a pascal. It is ported in that unit and converted at the boundary.

**The two design velocities.** Each is a blend of ``sqrt(2 target/rho)`` with a logarithmic
flow term, clamped - ``[25, 120]`` m/s at the suction, ``[10, 60]`` at the diffuser outlet -
and the outlet's enthalpy is what the diffuser velocity leaves. They are the class's own
correlations, not a standard's, which is why they are ported rather than replaced.

# The mass basis

``run`` works in ``kg/sec`` and ``J/kg`` throughout; this library's record is molar. The two
cross through the mixture's molar mass, once per flash, and the density the estimates are
built from is the **corrected** one - the cubic's volume with the Peneloux translation, which
is what NeqSim's ``getDensity("kg/m3")`` reports.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

from azoth.core.errors import InvalidInputError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import EjectorResult
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning
from azoth.eos import components as _components
from azoth.eos.reference.ph_flash import enthalpy_at
from azoth.eos.reference.ph_flash import ph_flash as ph_flash_solve
from azoth.eos.reference.pr_mass_density import pr_mass_density
from azoth.eos.reference.pr_molar_volume import pr_molar_volume
from azoth.eos.reference.ps_flash import entropy_at, ps_flash

#: The bar the class's own estimates are written in - not this library's unit.
BAR = 1.0e5


def ejector(
    motive_components: list[str],
    motive_n: Q,
    motive_z: list[float],
    motive_p: Q,
    motive_t: Q,
    suction_components: list[str],
    suction_n: Q,
    suction_z: list[float],
    suction_p: Q,
    suction_t: Q,
    discharge_pressure: Q,
    motive_nozzle_efficiency: float,
    suction_nozzle_efficiency: float,
    mixing_efficiency: float,
    diffuser_efficiency: float,
) -> EjectorResult:
    """Expand a motive stream, entrain a suction stream with it, and diffuse the mixture.

    Args:
        motive_components: the motive stream's substances, by name.
        motive_n: the motive inlet's molar flow.
        motive_z: the motive inlet's composition.
        motive_p: the motive inlet's pressure, which drives the machine.
        motive_t: the motive inlet's temperature.
        suction_components: the suction stream's substances, by name.
        suction_n: the suction inlet's molar flow.
        suction_z: the suction inlet's composition.
        suction_p: the suction inlet's pressure.
        suction_t: the suction inlet's temperature.
        discharge_pressure: the pressure the diffuser discharges at.
        motive_nozzle_efficiency: the motive nozzle's isentropic efficiency.
        suction_nozzle_efficiency: the suction nozzle's isentropic efficiency.
        mixing_efficiency: the mixing chamber's momentum-transfer efficiency.
        diffuser_efficiency: the diffuser's pressure-recovery efficiency.

    Returns:
        The outlet's record.

    Raises:
        InvalidInputError: where an efficiency is outside ``(0, 1]`` or the discharge
            pressure is not positive.
        OutOfRangeError: for a state the spec bounds.

    Example:
        >>> import azoth
        >>> q = azoth.ureg.Quantity
        >>> r = ejector(
        ...     ["methane", "n-butane"],
        ...     q(1.0, "mol/s"),
        ...     [0.9, 0.1],
        ...     q(30.0, "bar"),
        ...     q(400.0, "K"),
        ...     ["methane", "n-butane"],
        ...     q(0.5, "mol/s"),
        ...     [0.9, 0.1],
        ...     q(5.0, "bar"),
        ...     q(300.0, "K"),
        ...     q(10.0, "bar"),
        ...     0.75,
        ...     0.90,
        ...     0.85,
        ...     0.80,
        ... )
        >>> round(r.outlet_n.to("mol/s").magnitude, 6)
        1.5
    """
    spec = _spec()
    checks = checks_for(spec)
    warnings: list[Warning] = []

    efficiencies = {
        "motive_nozzle_efficiency": motive_nozzle_efficiency,
        "suction_nozzle_efficiency": suction_nozzle_efficiency,
        "mixing_efficiency": mixing_efficiency,
        "diffuser_efficiency": diffuser_efficiency,
    }
    apply_checks(
        checks.on_input,
        {
            "discharge_pressure": input_to_si(spec, "discharge_pressure", discharge_pressure),
            **efficiencies,
        }.get,
        warnings,
    )

    states = _route(
        motive_components,
        input_to_si(spec, "motive_n", motive_n),
        motive_z,
        input_to_si(spec, "motive_p", motive_p),
        motive_t,
        suction_components,
        input_to_si(spec, "suction_n", suction_n),
        suction_z,
        input_to_si(spec, "suction_p", suction_p),
        suction_t,
        input_to_si(spec, "discharge_pressure", discharge_pressure),
        efficiencies,
    )

    return EjectorResult(
        outlet_n=from_si(states.n, "mol/s"),
        outlet_z=states.z,
        outlet_p=from_si(states.p, "Pa"),
        outlet_t=states.temperature,
        outlet_h=from_si(states.h, "J/mol"),
        warnings=tuple(warnings),
    )


@dataclass(frozen=True, slots=True)
class EjectorStates:
    """The outlet's state, in SI."""

    n: float
    z: tuple[float, ...]
    p: float
    temperature: Q
    h: float


def _route(
    motive_components: list[str],
    motive_n: float,
    motive_z: list[float],
    motive_p: float,
    motive_t: Q,
    suction_components: list[str],
    suction_n: float,
    suction_z: list[float],
    suction_p: float,
    suction_t: Q,
    discharge_pressure: float,
    efficiencies: dict[str, float],
) -> EjectorStates:
    """The ejector's arithmetic, in SI, in the order it computes it.

    **One arithmetic, two consumers**: :func:`ejector` builds the result from this, and it is
    the twin of ``kernels/ejector.rs``. The kernel-level refusals live here for the reason
    Rust puts them in the kernel.
    """
    for name, value in efficiencies.items():
        if not 0.0 < value <= 1.0:
            raise InvalidInputError(
                name, f"an efficiency is a fraction in (0, 1], and {value} is not"
            )
    if discharge_pressure <= 0.0:
        raise InvalidInputError(
            "discharge_pressure",
            f"a discharge pressure is positive, and {discharge_pressure} Pa is not",
        )

    motive, motive_gas = _fluid(motive_components)
    suction, suction_gas = _fluid(suction_components)
    motive_mass = motive_n * _molar_mass(motive, motive_z)
    suction_mass = suction_n * _molar_mass(suction, suction_z)
    total_mass = motive_mass + suction_mass
    suction_bar = suction_p / BAR
    discharge_bar = discharge_pressure / BAR

    # `run` returns the motive stream unchanged when nothing flows.
    if total_mass <= 0.0:
        motive_h, _ = enthalpy_at(
            motive, motive_gas, motive_t.to("K").magnitude, motive_p, motive_z
        )
        return EjectorStates(
            n=motive_n,
            z=tuple(motive_z),
            p=motive_p,
            temperature=motive_t,
            h=float(motive_h),
        )

    mixing_bar = _estimate_mixing_pressure(suction_bar, discharge_bar, motive_mass, suction_mass)
    if mixing_bar > suction_bar:
        mixing_bar = suction_bar
    mixing_p = mixing_bar * BAR

    # **Every velocity in this block is a mass-basis quantity.** `run` works in `J/kg`, and
    # `sqrt(2 dh)` is a velocity only when `dh` is per kilogram - so the molar enthalpies this
    # library's flashes return are divided by the mixture's molar mass here and multiplied back
    # wherever a flash needs one. The Python twin of `kernels/ejector.rs`, variable for
    # variable.
    motive_molar_mass = _molar_mass(motive, motive_z)
    suction_molar_mass = _molar_mass(suction, suction_z)

    # The motive nozzle: isentropic to the mixing pressure, then the efficiency's drop.
    motive_h = _molar_enthalpy(motive, motive_gas, motive_t, motive_p, motive_z)
    motive_mass_h = motive_h / motive_molar_mass
    motive_s, _ = entropy_at(motive, motive_gas, motive_t.to("K").magnitude, motive_p, motive_z)
    motive_isentropic_h = _isentropic_enthalpy(motive, motive_gas, mixing_p, motive_s, motive_z)
    motive_actual_mass_h = motive_mass_h - efficiencies["motive_nozzle_efficiency"] * (
        motive_mass_h - motive_isentropic_h / motive_molar_mass
    )
    nozzle_velocity = (2.0 * max(motive_mass_h - motive_actual_mass_h, 0.0)) ** 0.5

    # The suction nozzle, the same route with its own efficiency.
    suction_h = _molar_enthalpy(suction, suction_gas, suction_t, suction_p, suction_z)
    suction_mass_h = suction_h / suction_molar_mass
    suction_s, _ = entropy_at(
        suction, suction_gas, suction_t.to("K").magnitude, suction_p, suction_z
    )
    suction_isentropic_h = _isentropic_enthalpy(
        suction, suction_gas, mixing_p, suction_s, suction_z
    )
    suction_actual_mass_h = suction_mass_h - efficiencies["suction_nozzle_efficiency"] * (
        suction_mass_h - suction_isentropic_h / suction_molar_mass
    )
    drawn_t = _flash_at_enthalpy(
        suction, suction_gas, mixing_p, suction_actual_mass_h * suction_molar_mass, suction_z
    )
    suction_from_enthalpy = (2.0 * max(suction_mass_h - suction_actual_mass_h, 0.0)) ** 0.5
    drawn_density = max(_corrected_density(suction, drawn_t, mixing_p, suction_z), 1.0e-9)
    design_suction_velocity = max(
        suction_from_enthalpy,
        _estimate_suction_velocity(suction_bar, discharge_bar, drawn_density, suction_mass),
    )

    # The momentum balance, in mass enthalpies, and the static enthalpy it leaves.
    total_motive_h = motive_actual_mass_h + 0.5 * nozzle_velocity**2
    total_suction_h = suction_actual_mass_h + 0.5 * design_suction_velocity**2
    ideal_mixing_velocity = (
        motive_mass * nozzle_velocity + suction_mass * design_suction_velocity
    ) / total_mass
    mixing_velocity = efficiencies["mixing_efficiency"] * ideal_mixing_velocity
    mixed_total_h = (motive_mass * total_motive_h + suction_mass * total_suction_h) / total_mass
    mixed_static_h = mixed_total_h - 0.5 * mixing_velocity**2

    # The joined fluid at the mixing pressure, then across the diffuser's pressure rise.
    joined_components, joined_z, joined_n = _join(
        motive_components, motive_n, motive_z, suction_components, suction_n, suction_z
    )
    joined, joined_gas = _fluid(joined_components)
    joined_molar_mass = _molar_mass(joined, joined_z)
    _flash_at_enthalpy(joined, joined_gas, mixing_p, mixed_static_h * joined_molar_mass, joined_z)

    recovered = efficiencies["diffuser_efficiency"] * 0.5 * mixing_velocity**2
    before_diffuser = mixed_static_h + recovered
    at_discharge_t = _flash_at_enthalpy(
        joined, joined_gas, discharge_pressure, before_diffuser * joined_molar_mass, joined_z
    )
    diffuser_density = max(
        _corrected_density(joined, at_discharge_t, discharge_pressure, joined_z), 1.0e-9
    )
    design_diffuser_velocity = _estimate_diffuser_velocity(
        mixing_bar, discharge_bar, diffuser_density, total_mass
    )
    final_mass_h = before_diffuser - 0.5 * design_diffuser_velocity**2

    temperature = _flash_at_enthalpy(
        joined, joined_gas, discharge_pressure, final_mass_h * joined_molar_mass, joined_z
    )
    final_state_h, _ = enthalpy_at(
        joined, joined_gas, temperature.to("K").magnitude, discharge_pressure, joined_z
    )
    return EjectorStates(
        n=joined_n,
        z=tuple(joined_z),
        p=discharge_pressure,
        temperature=temperature,
        h=float(final_state_h),
    )


def _fluid(components: list[str]) -> tuple[Any, Any]:
    """The mixture and the ideal-gas model for a set of component names.

    Both, because every state in this module needs one of each and resolving the pair twice
    would be two chances to disagree about which databank rows answered.
    """
    return _components.mixture_of(components, eos="pr")


def _molar_mass(mixture: Any, z: list[float]) -> float:
    """The mixture's molar mass, **in kg/mol as a bare number**.

    Bare, because every use of it here is inside an estimate written in kilograms; a quantity
    that leaked into one of those would make the comparison a unit error rather than an
    arithmetic one.
    """
    total = 0.0
    for component, zi in zip(mixture.components, z, strict=True):
        if component.molar_mass is None:
            raise InvalidInputError(
                "components",
                "a component carries no molar mass, and the machine's velocities are per kilogram",
            )
        total += zi * component.molar_mass.to("kg/mol").magnitude
    return float(total)


def _molar_enthalpy(
    mixture: Any, ideal_gas: Any, temperature: Q, p_si: float, z: list[float]
) -> float:
    """The state's molar enthalpy, J/mol."""
    h, _ = enthalpy_at(mixture, ideal_gas, temperature.to("K").magnitude, p_si, z)
    return float(h)


def _isentropic_enthalpy(
    mixture: Any, ideal_gas: Any, p_out: float, s_in: float, z: list[float]
) -> float:
    """The molar enthalpy after an isentropic expansion to ``p_out``.

    The entropy is the caller's, taken at the stream's *own* state - this only moves the
    pressure.
    """
    moved = ps_flash(mixture, ideal_gas, from_si(p_out, "Pa"), from_si(s_in, "J/(mol*K)"), z)
    h, _ = enthalpy_at(mixture, ideal_gas, moved.T.to("K").magnitude, p_out, z)
    return float(h)


def _flash_at_enthalpy(
    mixture: Any, ideal_gas: Any, p_out: float, h_molar: float, z: list[float]
) -> Q:
    """The temperature a state takes at a molar enthalpy and pressure."""
    moved = ph_flash_solve(mixture, ideal_gas, from_si(p_out, "Pa"), from_si(h_molar, "J/mol"), z)
    return moved.T


def _corrected_density(mixture: Any, temperature: Q, p_si: float, z: list[float]) -> float:
    """The state's mass density with the Peneloux translation, kg/m³.

    The same arithmetic ``Stream::corrected_density`` uses, which is what NeqSim's
    ``getDensity("kg/m3")`` reports - the class's own estimates are built from it.
    """
    from azoth.core.result import Phase as _Phase
    from azoth.eos.reference.pt_flash import pt_flash

    t_si = temperature.to("K").magnitude
    flash = pt_flash(mixture, temperature, from_si(p_si, "Pa"), z)
    # **A two-phase state has no one density**, and `Stream::corrected_density` refuses it for
    # that reason; the same refusal is here rather than a choice between the roots.
    if flash.phase is _Phase.TWO_PHASE:
        raise InvalidInputError(
            "state", f"a two-phase state at {t_si} K has no one density for the estimates"
        )
    root = flash.z_vapour if flash.phase is _Phase.ALL_VAPOUR else flash.z_liquid
    volume = pr_molar_volume(root, from_si(t_si, "K"), from_si(p_si, "Pa")).v
    molar_mass = _molar_mass(mixture, z)
    corrected = volume.to("m**3/mol").magnitude - mixture.volume_shift(z)
    density = pr_mass_density(from_si(molar_mass, "kg/mol"), from_si(corrected, "m**3/mol")).rho
    return float(density.to("kg/m**3").magnitude)


def _join(
    motive_components: list[str],
    motive_n: float,
    motive_z: list[float],
    suction_components: list[str],
    suction_n: float,
    suction_z: list[float],
) -> tuple[list[str], list[float], float]:
    """The two streams as one, by moles - ``SystemThermo.addFluid``'s arithmetic.

    The two inlets need not name the same substances, so the union is what the joined fluid
    is; a component only one side carries arrives at zero from the other.
    """
    names = list(dict.fromkeys([*motive_components, *suction_components]))
    amounts = [0.0] * len(names)
    for components, z, n in (
        (motive_components, motive_z, motive_n),
        (suction_components, suction_z, suction_n),
    ):
        for name, zi in zip(components, z, strict=True):
            amounts[names.index(name)] += n * zi
    total = sum(amounts)
    return names, [amount / total for amount in amounts], total


def _estimate_mixing_pressure(
    suction_pressure: float, discharge_pressure: float, motive_mass: float, suction_mass: float
) -> float:
    """The mixing pressure ``run`` estimates when ``setMixingPressure`` was never called.

    In bara - the floor is ``0.01`` bar and the margin is 3 % of the suction pressure.
    """
    if suction_pressure <= 0.0:
        return max(discharge_pressure, 0.0)
    entrainment_ratio = max(suction_mass, 0.0) / motive_mass if motive_mass > 1.0e-9 else 1.0
    clamped_ratio = _clamp(entrainment_ratio, 0.0, 5.0)
    pressure_lift = max(discharge_pressure - suction_pressure, 0.0)
    pressure_drop = pressure_lift * (0.1 + 0.03 * clamped_ratio)
    suction_margin = suction_pressure * 0.03
    estimated = suction_pressure - max(pressure_drop, suction_margin)
    return _clamp(estimated, max(0.01, suction_pressure * 0.4), suction_pressure)


def _estimate_suction_velocity(
    suction_pressure: float, discharge_pressure: float, suction_density: float, suction_mass: float
) -> float:
    """The suction velocity ``run`` designs to, in m/s."""
    density = max(suction_density, 1.0e-6)
    available_lift = max((discharge_pressure - suction_pressure) * BAR, 0.0)
    target_dynamic = max(available_lift * 0.02, 500.0)
    baseline = (2.0 * target_dynamic / density) ** 0.5
    volumetric_flow = suction_mass / density if suction_mass > 0.0 else 0.0
    flow_scaling = (
        50.0 + 30.0 * _log10(1.0 + volumetric_flow * 5.0) if volumetric_flow > 0.0 else 50.0
    )
    return _clamp(0.6 * baseline + 0.4 * flow_scaling, 25.0, 120.0)


def _estimate_diffuser_velocity(
    mixing_pressure: float, discharge_pressure: float, diffuser_density: float, total_mass: float
) -> float:
    """The diffuser outlet velocity ``run`` designs to, in m/s."""
    density = max(diffuser_density, 1.0e-6)
    pressure_recovery = max((discharge_pressure - mixing_pressure) * BAR, 0.0)
    target_dynamic = max(pressure_recovery * 0.01, 250.0)
    baseline = (2.0 * target_dynamic / density) ** 0.5
    volumetric_flow = total_mass / density if total_mass > 0.0 else 0.0
    flow_scaling = (
        20.0 + 15.0 * _log10(1.0 + volumetric_flow * 4.0) if volumetric_flow > 0.0 else 25.0
    )
    return _clamp(0.5 * baseline + 0.5 * flow_scaling, 10.0, 60.0)


def _log10(value: float) -> float:
    """The base-ten logarithm, spelled out so the twin says which one it means."""
    from math import log10

    return log10(value)


def _clamp(value: float, low: float, max_: float) -> float:
    """``Math.max(min, Math.min(value, max))``, in the argument order the class uses."""
    return max(low, min(value, max_))


def _spec() -> dict[str, object]:
    """The generated spec, imported here so the module has no import-time cycle."""
    from azoth._models_gen import model

    return model("process.ejector")
