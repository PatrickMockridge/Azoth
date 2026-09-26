"""``eos.phase_transport`` - a phase's transport properties, by NeqSim's phase-type dispatch.

Spec: ``specs/models/eos/phase_transport.toml``

The Python twin of ``crates/azoth-eos/src/phase_transport.rs``, written to mirror it line for
line: it composes the same reference kernels, in the same order, and resolves the same
component data.

# The dispatch is the answer

``PhysicalProperties`` is built by phase type and each subclass assigns its own models, so one
fluid's gas and liquid answer their properties from *different correlations*:

===========  ==========================  ==============================  ==================
phase        viscosity                   conductivity                    diffusivity
===========  ==========================  ==============================  ==================
gas          PFCT (heavy oil)            PFCT                            Chapman-Enskog
oil          PFCT                        PFCT                            Siddiqi-Lucas
aqueous      the polynom ``Viscosity``   the polynom ``Conductivity``     Siddiqi-Lucas
===========  ==========================  ==============================  ==================

# The two ladders, and the fallback that hides one of them

The liquid diffusivity divides by a **pure-component** viscosity, which comes from whichever
LIQVISC ladder the phase's viscosity class carries. An oil's PFCT class inherits the common-phase
one, whose model 2 branch is **empty** - and ``SiddiqiLucasMethod`` then falls back to the
*phase's own* viscosity, so the hole is the normal path rather than an error.
"""

from __future__ import annotations

from collections.abc import Sequence

from azoth.core.errors import InvalidInputError, OutOfRangeError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import PhaseTransportResult
from azoth.core.units import Q, from_si, input_to_si, ureg
from azoth.core.warnings import Warning
from azoth.eos.components import (
    entry,
    lennard_jones_pair,
    mixture_of,
    normal_boiling_molar_volume,
)
from azoth.eos.reference.aqueous_viscosity import aqueous_viscosity
from azoth.eos.reference.chapman_enskog_diffusivity import chapman_enskog_diffusivity
from azoth.eos.reference.effective_diffusion import effective_diffusion
from azoth.eos.reference.liquid_conductivity_polynom import liquid_conductivity_polynom
from azoth.eos.reference.liquid_viscosity_pure import liquid_viscosity_pure
from azoth.eos.reference.siddiqi_lucas_diffusivity import siddiqi_lucas_diffusivity
from azoth.eos.reference.thermal_conductivity import thermal_conductivity
from azoth.eos.reference.viscosity import viscosity

CALC_ID = "eos.phase_transport"

#: NeqSim's `PhaseType` values, as the spec spells them.
KINDS = ("gas", "oil", "aqueous")

#: `SiddiqiLucasMethod`'s own floor on the pure-component viscosity, in cP.
MIN_ETA_CP = 0.01

#: The ladder each liquid kind's viscosity class carries.
LADDERS = {"oil": "common_phase", "aqueous": "liquid"}


def phase_transport(
    components: Sequence[str],
    phase: str,
    T: Q,
    P: Q,
    z: Sequence[float],
) -> PhaseTransportResult:
    """A phase's viscosity, conductivity, binary diffusivity matrix and effective diffusivities.

    Args:
        components: the phase's substances, by name.
        phase: ``gas``, ``oil`` or ``aqueous`` - **which decides the correlations**.
        T: absolute temperature.
        P: absolute pressure.
        z: the phase's mole fractions.

    Returns:
        The four properties the segment model's transport snapshot reads.

    Raises:
        InvalidInputError: for a phase kind that is not one of the three, a composition of the
            wrong length, or a phase with fewer than two components - the effective diffusivity
            divides by an empty sum there and NeqSim answers a ``NaN``.
        OutOfRangeError: from whichever correlation refuses the state.

    Example:
        >>> import azoth
        >>> q = azoth.ureg.Quantity
        >>> r = phase_transport(
        ...     ["CO2", "water"], "aqueous", q(313.15, "K"), q(50.0e5, "Pa"),
        ...     [0.0006691762234084198, 0.9993308237765915])
        >>> round(r.k.magnitude, 13)
        0.6344057039895
    """
    from azoth._models_gen import model

    spec = model(CALC_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    if phase not in KINDS:
        raise InvalidInputError(
            "phase",
            f"`{phase}` is not one of NeqSim's three phase kinds: {', '.join(KINDS)}",
        )
    count = len(components)
    if len(z) != count:
        raise InvalidInputError(
            "z", f"a mixture of {count} component(s) needs {count} mole fractions, got {len(z)}"
        )
    if count < 2:
        raise InvalidInputError(
            "components",
            "one component: the effective diffusivity divides by the sum over the *other* "
            "components, so NeqSim answers a NaN here and this refuses instead of reporting one",
        )

    temperature = input_to_si(spec, "T", T)
    pressure = input_to_si(spec, "P", P)
    apply_checks(
        checks.on_input,
        {
            "T": temperature,
            "P": pressure,
            "z": min(float(value) for value in z),
        }.get,
        warnings,
    )

    mixture, ideal_gas = mixture_of(list(components))
    entries = [entry(name) for name in components]
    molar_mass: list[Q] = []
    for record in entries:
        if record.molar_mass is None:
            raise InvalidInputError(
                "components",
                "a component carries no molar mass, and every correlation here is a "
                "function of one",
            )
        molar_mass.append(record.molar_mass)

    # ---- The viscosity and the conductivity, by phase kind.
    if phase == "aqueous":
        mu = aqueous_viscosity(mixture, T, P, list(z)).viscosity
        k = liquid_conductivity_polynom(
            [record.liquid_conductivity for record in entries],
            molar_mass,
            list(z),
            T,
        ).k
    else:
        mu = viscosity(mixture, T, P, list(z)).mu
        k = thermal_conductivity(mixture, ideal_gas, T, P, list(z)).k

    # ---- The binary diffusivity matrix, by phase kind.
    matrix: list[list[Q]] = [[from_si(0.0, "m**2/s") for _ in range(count)] for _ in range(count)]
    for i in range(count):
        for j in range(count):
            if i == j:
                continue
            if phase == "gas":
                sigma, eps, _ = lennard_jones_pair(
                    entries[i].lennard_jones_diameter,
                    entries[i].lennard_jones_energy,
                    entries[j].lennard_jones_diameter,
                    entries[j].lennard_jones_energy,
                    molar_mass[i].to("kg/mol").magnitude,
                    molar_mass[j].to("kg/mol").magnitude,
                )
                matrix[i][j] = chapman_enskog_diffusivity(
                    molar_mass[i],
                    molar_mass[j],
                    ureg.Quantity(sigma * 1.0e-10, "m"),
                    ureg.Quantity(eps, "K"),
                    T,
                    P,
                ).d
            else:
                # **The aqueous form for every liquid**: the method's `autoSelectCorrelation` is
                # false, so an oil takes the aqueous correlation too.
                va = normal_boiling_molar_volume(
                    _liquid_density(entries[i]),
                    molar_mass[i].to("kg/mol").magnitude,
                    _critical_volume(entries[i]),
                )
                vb = normal_boiling_molar_volume(
                    _liquid_density(entries[j]),
                    molar_mass[j].to("kg/mol").magnitude,
                    _critical_volume(entries[j]),
                )
                # The pure-component viscosity the correlation divides by, from **this phase
                # kind's ladder** - and the class's own fallback where it is not positive, which
                # on an oil is the normal path: the common-phase ladder's model 2 branch is
                # empty, so the phase's own viscosity is what the correlation takes.
                pure = (
                    liquid_viscosity_pure(
                        LADDERS[phase],
                        entries[j].liqvisc_model,
                        *entries[j].liqvisc,
                        entries[j].Tc,
                        entries[j].Pc,
                        entries[j].omega,
                        T,
                        P,
                    )
                    .mu.to("Pa*s")
                    .magnitude
                    * 1000.0
                )
                eta_cp = pure if pure > 0.0 else mu.to("Pa*s").magnitude * 1000.0
                matrix[i][j] = siddiqi_lucas_diffusivity(
                    "aqueous",
                    ureg.Quantity(va * 1.0e-6, "m**3/mol"),
                    ureg.Quantity(vb * 1.0e-6, "m**3/mol"),
                    T,
                    ureg.Quantity(max(eta_cp, MIN_ETA_CP) * 1.0e-3, "Pa*s"),
                ).d

    # The assembly takes quantities, not magnitudes: it is the same kernel the id exposes.
    #
    # **A phase that is one substance has no effective diffusivity, and that is not this id's
    # refusal.** The assembly divides by the *other* components' fractions, so on a pure phase
    # there is nothing for the one component to diffuse into and ``eos.effective_diffusion``
    # refuses - rightly, for the question *it* asks. ``RateBasedPackedColumnTest``'s own lean
    # solvent is pure water, and NeqSim answers a **zero** vector there rather than refusing;
    # the Rust twin carries the same narrowing for the same measurement.
    try:
        assembled = effective_diffusion(matrix, list(z))
        d_effective = assembled.effective_diffusion
        warnings.extend(assembled.warnings)
    except OutOfRangeError as refusal:
        if refusal.input_field != "x":
            raise
        d_effective = tuple(from_si(0.0, "m**2/s") for _ in z)

    apply_checks(
        checks.derived,
        lambda name: (
            mu.to("Pa*s").magnitude
            if name == "mu"
            else (k.to("W/(m*K)").magnitude if name == "k" else None)
        ),
        warnings,
    )

    return PhaseTransportResult(
        mu=mu,
        k=k,
        d_binary=tuple(tuple(value for value in row) for row in matrix),
        d_effective=d_effective,
        warnings=tuple(warnings),
    )


def _liquid_density(record: object) -> float:
    """A component's normal liquid density in kg/m**3, zero where the row carries none.

    **Zero is the table's own spelling of absence** - `getNormalLiquidDensity` estimates from
    the critical volume there - so this passes it through rather than substituting one.
    """
    value = getattr(record, "liquid_density", None)
    return 0.0 if value is None else float(value.to("kg/m**3").magnitude)


def _critical_volume(record: object) -> float | None:
    """A component's critical volume in m**3/mol, `None` where the row carries none."""
    value = getattr(record, "critical_volume", None)
    return None if value is None else float(value.to("m**3/mol").magnitude)
