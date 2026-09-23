"""``reactions.kinetics`` - the Krishna-Standart mass-transfer rate matrix.

Spec: ``specs/models/reactions/kinetics.toml``

The Python twin of ``crates/azoth-reactions/src/kinetics.rs``, written to mirror it line
for line.

# One component's row

```text
coefficient = sum_r  k_r * prod_j c_j^(-nu_j) * prod_{k: nu_k nu_j > 0} c_k^(nu_k/nu_j)
```

where ``c_j = x_j * rho / M_j``, ``k`` is the reaction's rate factor, and the second
product runs only over the siblings on the **same side** of the reaction as the component
being asked about, excluding water and excluding the component itself.

# Three things the class does that a formula does not say

**The class builds three concentrations from two phases at once, and pairs them three
different ways.** :func:`_reaction_concentration` takes the interface's mole fraction over
the *bulk* phase's density and molar mass; :func:`_reactive_concentration` takes the
interface's fraction and density over the bulk's molar mass; and the sibling concentration
is the bulk's throughout. Kept as three functions rather than one so the pairing is
visible.

**The irreversibility test is ``1/K`` scaled by the reactant concentrations**, and it is
computed on every reaction whether or not the component appears in it. It sets a flag
rather than changing the coefficient.

**``phiInfinite`` is assigned, not accumulated**, so where a reaction has more than one
same-side sibling the *last* one in the reaction's own name order wins - which is why the
reactions cross as concatenated name lists and not as a matrix over the phase's species.
Where none produces one the class reads back the field it was constructed with, zero; this
model reports that same zero per call.
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from typing import Any

from azoth.core.errors import InvalidInputError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import KineticsResult
from azoth.core.units import Q, input_to_si
from azoth.core.warnings import Warning

#: The floor on ``|1/K|`` scaled by the reactant concentrations below which the class calls
#: a reaction irreversible, from ``if (Math.abs(irr) < 1e-3)``.
IRREVERSIBLE_THRESHOLD = 1e-3


def kinetics(
    components: Sequence[str],
    reaction_components: Sequence[str],
    reaction_lengths: Sequence[float],
    reaction_coefficients: Sequence[float],
    rate_factors: Sequence[float],
    equilibrium_constants: Sequence[float],
    fractions: Sequence[float],
    molar_masses: Sequence[float],
    density: Q,
    inter_fractions: Sequence[float],
    inter_density: Q,
    diffusion: Sequence[float],
) -> KineticsResult:
    """The whole phase's mass-transfer rate matrix.

    Args:
        components: the bulk phase's substances, by name, in the order ``fractions``,
            ``molar_masses`` and ``diffusion`` are given in.
        reaction_components: **the reactions' own species, concatenated in reaction order
            and in each reaction's own order**, ``reaction_lengths[i]`` of them for
            reaction ``i``. A name list rather than a stoichiometric matrix because a
            reaction's own order decides which sibling's phi it keeps.
        reaction_lengths: how many species each reaction names.
        reaction_coefficients: the signed coefficients, aligned one for one with
            ``reaction_components``.
        rate_factors: each reaction's rate factor at the interface's temperature.
        equilibrium_constants: each reaction's ``K`` at the bulk phase's state.
        fractions: the bulk phase's mole fractions.
        molar_masses: the bulk phase's molar masses, in kg/mol.
        density: the bulk phase's mass density.
        inter_fractions: the interface phase's mole fractions.
        inter_density: the interface phase's mass density.
        diffusion: the effective diffusion coefficients, in m**2/s, one per component.

    Returns:
        Each component's coefficient, its ``phiInfinite`` - zero where no reaction produced
        one - and a mask saying whether the irreversibility test fired while its row was
        built.

    Raises:
        InvalidInputError: where the lengths disagree, or where a reaction names a species
            the phase does not carry.

    Example:
        >>> import azoth
        >>> q = azoth.ureg.Quantity
        >>> r = kinetics(
        ...     ["water", "CO2"],
        ...     ["water", "CO2"],
        ...     [2.0],
        ...     [-2.0, -1.0],
        ...     [0.004327147905576517],
        ...     [4.412340363006617e-07],
        ...     [0.99, 0.01],
        ...     [0.018015, 0.04401],
        ...     q(1000.0, "kg/m**3"),
        ...     [0.99, 0.01],
        ...     q(1000.0, "kg/m**3"),
        ...     [1.0e-9, 1.0e-9],
        ... )
        >>> r.irreversible
        (0.0, 0.0)
    """
    spec = _spec()
    checks = checks_for(spec)
    warnings: list[Warning] = []
    bulk_density = input_to_si(spec, "density", density)
    interface_density = input_to_si(spec, "inter_density", inter_density)

    apply_checks(
        checks.on_input,
        {"density": bulk_density, "inter_density": interface_density}.get,
        warnings,
    )

    count = len(components)
    for field, values in [
        ("fractions", fractions),
        ("molar_masses", molar_masses),
        ("diffusion", diffusion),
        ("inter_fractions", inter_fractions),
    ]:
        if len(values) != count:
            raise InvalidInputError(field, f"{len(values)} value(s) against {count} component(s)")

    if len(reaction_components) != len(reaction_coefficients):
        raise InvalidInputError(
            "reaction_coefficients",
            f"{len(reaction_components)} name(s) and {len(reaction_coefficients)} coefficient(s)",
        )
    if len(rate_factors) != len(reaction_lengths) or len(equilibrium_constants) != len(
        reaction_lengths
    ):
        raise InvalidInputError(
            "rate_factors",
            f"{len(reaction_lengths)} reaction(s) and {len(rate_factors)} rate factor(s) "
            f"against {len(equilibrium_constants)} equilibrium constant(s)",
        )

    # Cut the concatenated lists into the reactions the class iterates, each with its own
    # names in its own order.
    reactions: list[tuple[list[str], list[float], float, float]] = []
    cursor = 0
    for index, raw in enumerate(reaction_lengths):
        if not math.isfinite(raw) or raw < 1.0 or raw != int(raw):
            raise InvalidInputError(
                "reaction_lengths", f"reaction {index} names {raw} species, which is not a count"
            )
        length = int(raw)
        if cursor + length > len(reaction_components):
            raise InvalidInputError(
                "reaction_lengths",
                f"reaction {index} runs past the {len(reaction_components)}-species reaction "
                f"list at {cursor}",
            )
        reactions.append(
            (
                list(reaction_components[cursor : cursor + length]),
                list(reaction_coefficients[cursor : cursor + length]),
                rate_factors[index],
                equilibrium_constants[index],
            )
        )
        cursor += length
    if cursor != len(reaction_components):
        raise InvalidInputError(
            "reaction_lengths",
            f"the lengths sum to {cursor} and the reaction list has "
            f"{len(reaction_components)} species",
        )

    # **The boundary is mixed by design**: a vector whose spec declares a unit arrives as a
    # pint quantity and has to be converted, and a dimensionless one arrives as the bare
    # number it is. `molar_masses` and `diffusion` are the two here.
    masses = [_si(spec, "molar_masses", value) for value in molar_masses]
    diffusivities = [_si(spec, "diffusion", value) for value in diffusion]

    names = list(components)
    bulk = (names, list(fractions), masses, bulk_density)
    interface = (names, list(inter_fractions), masses, interface_density)

    coefficient: list[float] = []
    phi_infinite: list[float] = []
    irreversible: list[float] = []
    for name in components:
        row = _rate_matrix(reactions, bulk, interface, name, diffusivities)
        coefficient.append(row[0])
        phi_infinite.append(0.0 if row[1] is None else row[1])
        irreversible.append(1.0 if row[2] else 0.0)

    return KineticsResult(
        coefficient=tuple(coefficient),
        phi_infinite=tuple(phi_infinite),
        irreversible=tuple(irreversible),
        warnings=tuple(warnings),
    )


def _index(phase: tuple[list[str], list[float], list[float], float], name: str) -> int:
    try:
        return phase[0].index(name)
    except ValueError:
        raise InvalidInputError("names", f"`{name}` is not in the phase") from None


def _reaction_concentration(
    inter_phase: tuple[list[str], list[float], list[float], float],
    phase: tuple[list[str], list[float], list[float], float],
    name: str,
) -> float:
    """The interface's mole fraction over the **bulk** phase's density and molar mass."""
    here = _index(inter_phase, name)
    mass = _index(phase, name)
    return inter_phase[1][here] * phase[3] / phase[2][mass]


def _reactive_concentration(
    inter_phase: tuple[list[str], list[float], list[float], float],
    phase: tuple[list[str], list[float], list[float], float],
    name: str,
) -> float:
    """The interface's fraction and density over the **bulk** phase's molar mass."""
    here = _index(inter_phase, name)
    mass = _index(phase, name)
    return inter_phase[1][here] * inter_phase[3] / phase[2][mass]


def _bulk_concentration(
    phase: tuple[list[str], list[float], list[float], float], name: str
) -> float:
    """The bulk phase's own ``x rho / M``."""
    index = _index(phase, name)
    return phase[1][index] * phase[3] / phase[2][index]


def _rate_matrix(
    reactions: list[tuple[list[str], list[float], float, float]],
    phase: tuple[list[str], list[float], list[float], float],
    inter_phase: tuple[list[str], list[float], list[float], float],
    component: str,
    diffusion: list[float],
) -> tuple[float, float | None, bool]:
    """``calcReacMatrix``: one component's row, as ``(coefficient, phi, irreversible)``."""
    if len(diffusion) != len(phase[0]):
        raise InvalidInputError(
            "diffusion", f"{len(diffusion)} coefficient(s) against {len(phase[0])} species"
        )
    component_index = _index(phase, component)

    coefficient = 0.0
    phi_infinite: float | None = None
    irreversible = False

    for names, stoc_coefs, rate_factor, equilibrium_constant in reactions:
        if len(names) != len(stoc_coefs):
            raise InvalidInputError(
                "reactions",
                f"`{' '.join(names)}` has {len(names)} name(s) and "
                f"{len(stoc_coefs)} coefficient(s)",
            )
        ktemp = rate_factor

        irr = 1.0 / equilibrium_constant
        for j, name in enumerate(names):
            irr *= _reaction_concentration(inter_phase, phase, name) ** (-stoc_coefs[j])
        if abs(irr) < IRREVERSIBLE_THRESHOLD:
            irreversible = True

        for j, name in enumerate(names):
            if name != component:
                continue
            for k, sibling in enumerate(names):
                # The same side of the reaction, not the component itself, and never water.
                if stoc_coefs[k] * stoc_coefs[j] <= 0.0 or k == j or sibling == "water":
                    continue
                exponent = stoc_coefs[k] / stoc_coefs[j]
                reactive = _reactive_concentration(inter_phase, phase, component)
                sibling_concentration = _bulk_concentration(phase, sibling)
                ktemp *= sibling_concentration**exponent

                sibling_index = _index(phase, sibling)
                forward = math.sqrt(diffusion[component_index] / diffusion[sibling_index])
                backward = math.sqrt(diffusion[sibling_index] / diffusion[component_index])
                phi_infinite = forward + backward * sibling_concentration / (exponent * reactive)
        coefficient += ktemp

    return coefficient, phi_infinite, irreversible


def _si(spec: dict[str, Any], name: str, value: float | Q) -> float:
    """A declared input as its SI magnitude, quantity or not."""
    # A quantity or a bare number. `Q` is a type alias rather than a class, so the check
    # is made the other way round: the numeric branch is the one `isinstance` can name.
    if isinstance(value, int | float):
        return float(value)
    return input_to_si(spec, name, value)


def _spec() -> dict[str, Any]:
    from azoth._models_gen import model

    return model("reactions.kinetics")
