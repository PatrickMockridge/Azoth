"""``reactions.equilibrium_constant`` - one reaction's ``K`` and its derivative.

```text
ln K        = K1 + K2/T + K3*ln T + K4*T
d(ln K)/dT  = -K1/T**2 + K2/T + K3
dH          = d(ln K)/dT * R * T**2
```

Spec: ``specs/calcs/reactions/equilibrium_constant.toml``

NeqSim: ``ChemicalReaction.getK(PhaseInterface)`` and
``getReactionHeat(PhaseInterface)``. The Python twin of
``crates/azoth-reactions/src/equilibrium_constant.rs``, written to mirror it line for
line.

# What the correlation is

A van't Hoff term, a heat-capacity term and a linear term, fitted per reaction and per
source. ``K3`` is not decoration: for ``CO2water`` it is ``-39.440767``, so the ``ln T``
term contributes about ``-225`` to a ``ln K`` of ``-14.6``. A port that dropped it would
still be close at 298.15 K, where the fit's constants nearly absorb it, and wrong
everywhere else.

# What is not read

The row's ``TREF``, which is the fit's anchor and not a point where ``K`` takes any
particular value, and the row's ``ACTENERGY`` and rate factor, which parameterise the
rate law rather than the constant.
"""

from __future__ import annotations

import math

from azoth._registry_gen import spec as _spec_for
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import EquilibriumConstantResult
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning
from azoth.reactions.reference import _tables
from azoth.reactions.reference._tables import source_path

CALC_ID = "reactions.equilibrium_constant"

#: `ThermodynamicConstantsInterface.R`, in J/(mol*K). The table's `ACTENERGY` is not in
#: J/mol, and nothing in the correlation below reads it, so the two never meet here.
GAS_CONSTANT = 8.3144621


def equilibrium_constant(source: str, reaction: str, T: Q) -> EquilibriumConstantResult:
    """One reaction's equilibrium constant, derivative and heat of reaction.

    Args:
        source: which of the three reaction tables to read - ``"standard"``,
            ``"pitzer"`` or ``"kent-eisenberg"``. Required rather than defaulting,
            because the three are different standard states whose answers part company
            away from 298 K.
        reaction: the reaction's name as its table spells it, e.g. ``"CO2water"``.
        T: absolute temperature at which the constant is evaluated.

    Raises:
        InvalidInputError: if ``source`` is not one of the three.
        PropertyUnavailableError: if the source carries no row by that name. This is a
            different failure from a bad argument: a caller has to be able to tell "you
            gave me nonsense" from "I do not have that reaction".
        OutOfRangeError: if ``T`` is not positive. The correlation takes ``ln T``.

    Example:
        >>> import azoth
        >>> q = azoth.ureg.Quantity
        >>> r = equilibrium_constant("standard", "CO2water", q(298.15, "K"))
        >>> round(r.ln_k, 6)
        -14.63369
    """
    spec = _spec_for(CALC_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    values = {"T": input_to_si(spec, "T", T)}
    apply_checks(checks.on_input, values.get, warnings)

    t = values["T"]

    # Refused before the row is looked up, so a bad source is reported as a bad source
    # rather than as a missing reaction.
    source_path(source)

    row = _tables.reaction(source, reaction)

    k1, k2, k3, k4 = row.coefficients

    ln_k = k1 + k2 / t + k3 * math.log(t) + k4 * t
    k = math.exp(ln_k)
    ln_k_derivative = -k2 / (t * t) + k3 / t + k4
    reaction_heat = ln_k_derivative * t * t * GAS_CONSTANT

    apply_checks(
        checks.derived,
        lambda name: {"ln_k": ln_k, "k": k}.get(name),
        warnings,
    )

    return EquilibriumConstantResult(
        ln_k=ln_k,
        k=k,
        ln_k_derivative=from_si(ln_k_derivative, "1/K"),
        reaction_heat=from_si(reaction_heat, "J/mol"),
        reference=row.reference,
        warnings=tuple(warnings),
    )
