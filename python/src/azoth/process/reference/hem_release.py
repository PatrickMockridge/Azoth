"""``process.hem_release`` - NeqSim's homogeneous-equilibrium orifice release.

Spec: ``specs/models/process/hem_release.toml``

The Python twin of ``crates/azoth-process/src/models/hem_release.rs``, written to mirror it.
The model samples an isentropic expansion from the upstream stagnation state down to the
receiving pressure, maximises ``rho sqrt(2 (h0 - h))`` over pressure, and reports the throat
as the flow-limiting station. The throat's equilibrium sound speed is
:func:`azoth.eos.reference.equilibrium_sound_speed.equilibrium_sound_speed`.

The fixed-entropy sample is the shared :mod:`azoth.eos.reference._entropy_state` arithmetic
with this class's root policy, and the closure thresholds are the class's own.
"""

from __future__ import annotations

import math
from typing import Any, NamedTuple

from azoth.core.errors import InvalidInputError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import HemReleaseDiagnostic, HemReleaseResult, HemReleaseStation
from azoth.core.units import Q, from_si, input_to_si, quantity
from azoth.core.warnings import Warning
from azoth.eos import components as _components
from azoth.eos.reference._entropy_state import (
    HEM,
    EntropyState,
    EntropyStateFailure,
    Evaluations,
    flashed_state,
    isentropic_state,
    molar_mass_of,
)
from azoth.eos.reference.equilibrium_sound_speed import equilibrium_sound_speed
from azoth.eos.reference.pt_flash import check_composition

MODEL_ID = "process.hem_release"

#: The class's pressure samples between the receiving pressure and the upstream state.
SAMPLES = 48
#: ``PRESSURE_TOLERANCE``: the golden section's bracket as a fraction of the upstream pressure.
PRESSURE_TOLERANCE = 1.0e-7
#: The golden section's iteration cap.
REFINE_ITERATIONS = 80
#: ``sample()``'s accepted specific-entropy residual, J/(kg K).
ENTROPY_CLOSURE = 1.0e-5
#: ``checkEquilibrium``'s component-inventory residual, over the total moles.
BALANCE_CLOSURE = 1.0e-8
#: ``checkEquilibrium``'s interphase log-fugacity residual.
FUGACITY_CLOSURE = 1.0e-5
#: The smallest enthalpy drop that is still an expansion, J/kg.
ENTHALPY_FLOOR = -1.0e-5
#: The choked-flow test's margin over the receiving pressure.
CHOKED_MARGIN = 1.0e-5
#: The class's Mach band for the throat sound speed.
MACH_BAND = 0.05
#: The CO2 solid-capable guard's temperature, K.
CO2_SOLID_K = 216.592
#: The relative pressure step the throat sound speed is evaluated at.
SOUND_SPEED_RELATIVE_STEP = 1.0e-3

STATIONS = ("upstream_stagnation", "throat_critical", "orifice_exit", "ambient_expanded")
PHASES = ("all_vapour", "all_liquid", "two_phase")


class _Failure(Exception):
    """A failure the class reports inside its result rather than as an exception."""

    def __init__(self, message: str, *, unsupported: bool = False) -> None:
        super().__init__(message)
        self.unsupported = unsupported


class _Sample(NamedTuple):
    """One sample of the isentropic expansion."""

    state: EntropyState
    velocity: float
    flux: float


class _Station(NamedTuple):
    """One station row, in SI magnitudes."""

    station: str
    pressure: float
    temperature: float
    density: float
    molar_enthalpy: float
    molar_entropy: float
    velocity: float
    mass_flux: float
    phase: str
    gas_mass_fraction: float


class _States(NamedTuple):
    """Everything the class forms, in SI magnitudes.

    ``upstream``, ``samples`` and ``best`` are carried for ``azoth.process.layers``, which
    reads the reference's own factored entry point rather than re-deriving the stations.
    """

    mass_flow: float | None
    choked: bool
    sound_speed: float | None
    status: str
    component_mole_fractions: tuple[float, ...]
    component_mass_fractions: tuple[float, ...]
    stations: tuple[_Station, ...]
    diagnostics: tuple[tuple[str, str], ...]
    upstream: EntropyState | None
    samples: tuple[_Sample, ...]
    best: _Sample | None


def hem_release(
    components: list[str],
    upstream_z: list[float],
    upstream_t: Q,
    upstream_p: Q,
    diameter: Q,
    discharge_coefficient: float,
    back_pressure: Q,
) -> HemReleaseResult:
    """The release rate of a homogeneous-equilibrium short-opening orifice.

    Args:
        components: the fluid's substances, by name.
        upstream_z: the upstream stagnation composition, in ``components``' order.
        upstream_t: the upstream stagnation temperature.
        upstream_p: the upstream stagnation pressure.
        diameter: the circular opening's diameter.
        discharge_coefficient: the effective-area factor in ``(0, 1]``.
        back_pressure: the receiving absolute pressure.

    Raises:
        InvalidInputError: for a state or opening the class refuses at construction, a
            composition that is not one entry per component and sums to one, or whatever
            the databank refuses.
    """
    spec = _spec()
    checks = checks_for(spec)
    warnings: list[Warning] = []

    t = input_to_si(spec, "upstream_t", upstream_t)
    p = input_to_si(spec, "upstream_p", upstream_p)
    diameter_m = input_to_si(spec, "diameter", diameter)
    back = input_to_si(spec, "back_pressure", back_pressure)
    apply_checks(
        checks.on_input,
        {
            "upstream_t": t,
            "upstream_p": p,
            "diameter": diameter_m,
            "discharge_coefficient": discharge_coefficient,
            "back_pressure": back,
        }.get,
        warnings,
    )
    if not 0.0 < discharge_coefficient <= 1.0:
        raise InvalidInputError(
            "discharge_coefficient",
            f"the effective-area factor must be in (0, 1], and {discharge_coefficient} is not",
        )

    try:
        states = _states(
            components,
            list(upstream_z),
            t,
            p,
            diameter_m,
            float(discharge_coefficient),
            back,
        )
    except _Failure as failure:
        status = "unsupported" if failure.unsupported else "invalid"
        code = "UNSUPPORTED_PHYSICS" if failure.unsupported else "CALCULATION_FAILED"
        return HemReleaseResult(
            status=status,
            mass_flow_rate=None,
            choked=False,
            throat_sound_speed=None,
            component_mole_fractions=(),
            component_mass_fractions=(),
            stations=(),
            diagnostics=(HemReleaseDiagnostic(code=code, message=str(failure)),),
            diagnostic_codes=(code,),
            warnings=tuple(warnings),
        )

    return HemReleaseResult(
        status=states.status,
        mass_flow_rate=(None if states.mass_flow is None else from_si(states.mass_flow, "kg/s")),
        choked=states.choked,
        throat_sound_speed=(
            None if states.sound_speed is None else from_si(states.sound_speed, "m/s")
        ),
        component_mole_fractions=states.component_mole_fractions,
        component_mass_fractions=states.component_mass_fractions,
        stations=tuple(
            HemReleaseStation(
                station=row.station,
                pressure=from_si(row.pressure, "Pa"),
                temperature=from_si(row.temperature, "K"),
                density=from_si(row.density, "kg/m**3"),
                molar_enthalpy=from_si(row.molar_enthalpy, "J/mol"),
                molar_entropy=from_si(row.molar_entropy, "J/(mol*K)"),
                velocity=from_si(row.velocity, "m/s"),
                mass_flux=from_si(row.mass_flux, "kg/(m**2*s)"),
                phase=row.phase,
                gas_mass_fraction=row.gas_mass_fraction,
            )
            for row in states.stations
        ),
        diagnostics=tuple(
            HemReleaseDiagnostic(code=code, message=message) for code, message in states.diagnostics
        ),
        diagnostic_codes=tuple(code for code, _ in states.diagnostics),
        warnings=tuple(warnings),
    )


def _states(
    components: list[str],
    z: list[float],
    t: float,
    p: float,
    diameter: float,
    discharge_coefficient: float,
    back: float,
) -> _States:
    """The whole calculation in SI magnitudes, for the public function and the layer diff."""
    mixture, ideal_gas = _components.mixture_of(components, eos="pr")
    if len(z) != len(mixture.components):
        raise InvalidInputError(
            "upstream_z",
            f"a mixture of {len(mixture.components)} components needs that many mole "
            f"fractions, and {z!r} has {len(z)}",
        )
    check_composition(list(z), len(mixture.components), "upstream_z")
    mass = molar_mass_of(mixture, list(z))
    mass_fractions = _component_mass_fractions(mixture, z)
    names = mixture.names if mixture.names is not None else tuple(components)
    _support_inputs(t, names, z)

    evaluations = Evaluations()
    diagnostics = [
        (
            "MODEL_ASSUMPTIONS",
            "Equilibrium, zero slip, isentropic acceleration; no pipe friction, heat "
            "transfer or solids",
        ),
        ("EXIT_ALIASES_THROAT", "Zero-length opening: exit equals accepted throat"),
        (
            "ENTROPY_SOLVER",
            "Pure-fluid PS flash retains saturation quality"
            if len(mixture.components) == 1
            else "Cold TP flashes in bracketed mixture entropy root",
        ),
    ]

    try:
        upstream = flashed_state(mixture, ideal_gas, z, mass, p, t, evaluations)
        _closures(upstream)
        _support(upstream, names, z)
        target = upstream.molar_entropy / mass
        context = _Context(
            mixture=mixture,
            ideal_gas=ideal_gas,
            z=z,
            names=names,
            molar_mass=mass,
            target=target,
            upstream_h=upstream.molar_enthalpy,
            upstream_pressure=p,
            start_temperature=t,
        )
        upstream_sample = _Sample(upstream, 0.0, 0.0)

        if p <= back:
            diagnostics.append(
                ("NO_FORWARD_FLOW", "Upstream pressure is at or below receiving pressure")
            )
            rows = tuple(
                _station(upstream_sample, name)
                for name in (
                    "upstream_stagnation",
                    "throat_critical",
                    "orifice_exit",
                )
            )
            return _States(
                mass_flow=0.0,
                choked=False,
                sound_speed=None,
                status="valid",
                component_mole_fractions=tuple(z),
                component_mass_fractions=mass_fractions,
                stations=rows,
                diagnostics=tuple(diagnostics),
                upstream=upstream,
                samples=(),
                best=None,
            )

        samples: list[_Sample] = []
        for index in range(SAMPLES + 1):
            if index == SAMPLES:
                samples.append(upstream_sample)
                continue
            pressure = back * math.exp(math.log(p / back) * index / SAMPLES)
            samples.append(context.sample(pressure, evaluations))

        best = samples[0]
        for index in range(1, SAMPLES):
            if (
                samples[index].flux >= samples[index - 1].flux
                and samples[index].flux >= samples[index + 1].flux
            ):
                candidate = context.refine(
                    samples[index - 1].state.pressure,
                    samples[index + 1].state.pressure,
                    evaluations,
                )
                if candidate.flux > best.flux:
                    best = candidate
        for sample in samples:
            if sample.flux > best.flux:
                best = sample

        choked = best.state.pressure > back * (1.0 + CHOKED_MARGIN)
        sound_speed: float | None = None
        warning = False
        if choked:
            try:
                acoustic = equilibrium_sound_speed(
                    components,
                    quantity(best.state.temperature, "K"),
                    quantity(best.state.pressure, "Pa"),
                    z,
                    relative_pressure_step=SOUND_SPEED_RELATIVE_STEP,
                    specific_entropy=quantity(target, "J/(kg*K)"),
                )
            except Exception as failure:
                raise _Failure(str(failure)) from failure
            if acoustic.status == "converged" and acoustic.sound_speed is not None:
                sound_speed = acoustic.sound_speed.to("m/s").magnitude
                mach = best.velocity / sound_speed
                if abs(mach - 1.0) > MACH_BAND:
                    warning = True
                    diagnostics.append(
                        (
                            "THROAT_MACH_MISMATCH",
                            f"Equilibrium Mach={mach:.6}; phase-boundary maximum may be nonsmooth",
                        )
                    )
            else:
                warning = True
                diagnostics.append(
                    (
                        "ACOUSTIC_UNAVAILABLE",
                        f"{acoustic.status}: {acoustic.message}",
                    )
                )

        rows = (
            _station(upstream_sample, "upstream_stagnation"),
            _station(best, "throat_critical"),
            _station(best, "orifice_exit"),
            _station(samples[0], "ambient_expanded"),
        )
        rate = discharge_coefficient * math.pi * diameter * diameter / 4.0 * best.flux
        return _States(
            mass_flow=rate,
            choked=choked,
            sound_speed=sound_speed,
            status="valid_with_warnings" if warning else "valid",
            component_mole_fractions=tuple(z),
            component_mass_fractions=mass_fractions,
            stations=rows,
            diagnostics=tuple(diagnostics),
            upstream=upstream,
            samples=tuple(samples),
            best=best,
        )
    except EntropyStateFailure as failure:
        raise _Failure(str(failure)) from failure


class _Context(NamedTuple):
    """The state one release sample reads."""

    mixture: Any
    ideal_gas: Any
    z: list[float]
    names: tuple[str, ...]
    molar_mass: float
    target: float
    upstream_h: float
    upstream_pressure: float
    start_temperature: float

    def sample(self, pressure: float, evaluations: Evaluations) -> _Sample:
        """One isentropic sample at a pressure."""
        try:
            state = isentropic_state(
                self.mixture,
                self.ideal_gas,
                self.z,
                self.molar_mass,
                self.target,
                pressure,
                self.start_temperature,
                evaluations,
                HEM,
            )
        except EntropyStateFailure as failure:
            raise _Failure(str(failure)) from failure
        _closures(state)
        _support(state, self.names, self.z)
        drop = (self.upstream_h - state.molar_enthalpy) / self.molar_mass
        if not math.isfinite(drop) or drop < ENTHALPY_FLOOR:
            raise _Failure("Expansion increased enthalpy or returned a nonfinite value")
        velocity = math.sqrt(2.0 * max(0.0, drop))
        flux = state.density * velocity
        if not math.isfinite(flux) or flux < 0.0:
            raise _Failure(f"the mass flux at {pressure} Pa is {flux} kg/(m2 s)")
        return _Sample(state, velocity, flux)

    def refine(self, low: float, high: float, evaluations: Evaluations) -> _Sample:
        """The golden section over a sampled local maximum, in linear pressure."""
        ratio = (math.sqrt(5.0) - 1.0) / 2.0
        left = self.sample(high - ratio * (high - low), evaluations)
        right = self.sample(low + ratio * (high - low), evaluations)
        for _ in range(REFINE_ITERATIONS):
            if high - low <= PRESSURE_TOLERANCE * self.upstream_pressure:
                return left if left.flux > right.flux else right
            if left.flux < right.flux:
                low = left.state.pressure
                left = right
                right = self.sample(low + ratio * (high - low), evaluations)
            else:
                high = right.state.pressure
                right = left
                left = self.sample(high - ratio * (high - low), evaluations)
        raise _Failure("Critical pressure search did not converge")


def _closures(state: EntropyState) -> None:
    """The class's per-sample closure checks."""
    if abs(state.entropy_residual) > ENTROPY_CLOSURE:
        raise _Failure(f"Specific entropy closure failed: {state.entropy_residual} J/(kg K)")
    if state.component_residual > BALANCE_CLOSURE:
        raise _Failure("Component inventory closure failed")
    if state.phase_composition_residual > BALANCE_CLOSURE:
        raise _Failure("Phase composition does not close")
    if state.fugacity_residual > FUGACITY_CLOSURE:
        raise _Failure("Interphase fugacity closure failed")
    if (
        not state.density > 0.0
        or not math.isfinite(state.density)
        or not state.temperature > 0.0
        or not math.isfinite(state.temperature)
        or not state.pressure > 0.0
        or not math.isfinite(state.pressure)
    ):
        raise _Failure(f"the sample at {state.pressure:.6f} Pa is not a finite positive state")


def _support_inputs(temperature: float, names: tuple[str, ...], z: list[float]) -> None:
    """The class's upstream unsupported-regime check, before the first flash."""
    for name, fraction in zip(names, z, strict=True):
        if name.strip().lower() == "co2" and fraction > 0.0 and temperature < CO2_SOLID_K:
            raise _Failure(
                "CO2 below 216.592 K requires a separately assessed solid-capable model",
                unsupported=True,
            )


def _support(state: EntropyState, names: tuple[str, ...], z: list[float]) -> None:
    """The class's unsupported-regime checks this port can reach."""
    if state.phase == "trivial":
        raise _Failure(
            f"the flash found a trivial solution at {state.pressure:.6f} Pa; the phase "
            "is not established"
        )
    for name, fraction in zip(names, z, strict=True):
        if name.strip().lower() == "co2" and fraction > 0.0 and state.temperature < CO2_SOLID_K:
            raise _Failure(
                "CO2 below 216.592 K requires a separately assessed solid-capable model",
                unsupported=True,
            )


def _station(sample: _Sample, station: str) -> _Station:
    """One station row from a sample."""
    state = sample.state
    return _Station(
        station=station,
        pressure=state.pressure,
        temperature=state.temperature,
        density=state.density,
        molar_enthalpy=state.molar_enthalpy,
        molar_entropy=state.molar_entropy,
        velocity=sample.velocity,
        mass_flux=state.density * sample.velocity,
        phase=state.phase,
        gas_mass_fraction=state.gas_mass_fraction,
    )


def _component_mass_fractions(mixture: Any, z: list[float]) -> tuple[float, ...]:
    """The upstream overall mass fractions."""
    masses: list[float] = []
    total = 0.0
    for index, (component, fraction) in enumerate(zip(mixture.components, z, strict=True)):
        if component.molar_mass is None:
            raise InvalidInputError("molar_mass", f"component {index} carries no molar mass")
        mass = fraction * component.molar_mass.to("kg/mol").magnitude
        total += mass
        masses.append(mass)
    if not total > 0.0 or not math.isfinite(total):
        raise InvalidInputError("upstream_z", f"the mixture's total mass is {total} kg/mol")
    return tuple(mass / total for mass in masses)


def _spec() -> dict[str, object]:
    """The generated spec, imported here so the module has no import-time cycle."""
    from azoth._models_gen import model

    return model(MODEL_ID)
