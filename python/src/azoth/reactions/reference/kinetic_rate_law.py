"""``reactions.kinetic_rate_law`` - the rate factor a reaction's own law answers.

Spec: ``specs/models/reactions/kinetic_rate_law.toml``

The Python twin of ``crates/azoth-reactions/src/kinetic_rate_law.rs``, written to mirror
it line for line.

# One method, two laws, and a selector decides which

```text
legacy      2.576e9 * exp(-6024 / T) / 1000
arrhenius   rate * exp(-Ea / R * (1/T - 1/T_ref))
```

**The absent selector answers the legacy law.** ``getKineticRateLaw`` returns
``LEGACY_TEMPERATURE_CORRELATION`` when its field is unset, and every reaction
``chemicalReactionInit`` builds is in that state - so the legacy law is what a fluid's own
reactions run.

**The legacy law reads none of the three parameters.** ``2.576e9`` and ``6024`` are
literals in the method, so every reaction of a fluid gets the same rate factor at a given
temperature and the stored ``rateFactor`` and ``ACTENERGY`` columns are never evaluated
there. That is why the three Arrhenius parameters are required inputs the legacy branch
*ignores* rather than optional ones.

The gas constant is :data:`.equilibrium_constant.GAS_CONSTANT`, which is
``ThermodynamicConstantsInterface.R`` - the same ``8.3144621`` the equilibrium constant is
built with, and not the ``8.314462`` the RAND solver carries.
"""

from __future__ import annotations

import math
from typing import Any

from azoth.core.errors import InvalidInputError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import KineticRateLawResult
from azoth.core.units import Q, input_to_si
from azoth.core.warnings import Warning
from azoth.reactions.reference.equilibrium_constant import GAS_CONSTANT

#: The legacy correlation's prefactor, from the method's own literal.
LEGACY_PREFACTOR = 2.576e9

#: The temperature it is written against, in K - ``6024`` as the method spells it.
LEGACY_ACTIVATION_TEMPERATURE = 6024.0

#: The divisor the result is scaled by, from the method's ``/ 1000.0``.
LEGACY_DIVISOR = 1000.0

#: The two laws, and the spellings the selector accepts for them.
_LAWS = {
    "legacy": "legacy",
    "legacy_temperature_correlation": "legacy",
    "arrhenius": "arrhenius",
    "reference_arrhenius": "arrhenius",
}


def kinetic_rate_law(
    law: str,
    T: Q,
    reference_rate: float,
    activation_energy: Q,
    reference_temperature: Q,
) -> KineticRateLawResult:
    """The reaction's rate factor at ``T``, by the selected law.

    Args:
        law: the selector's name, ``legacy`` or ``arrhenius``. The class answers the legacy
            law when its field is absent, which is the state a fluid's own reactions are in.
        T: absolute temperature.
        reference_rate: the rate at ``reference_temperature``. **Not read by the legacy
            branch**, which is why the cases set it and still get the method's literal.
        activation_energy: the activation energy, in J/mol.
        reference_temperature: the temperature ``reference_rate`` is stated at.

    Returns:
        The rate factor, and any caveats.

    Raises:
        InvalidInputError: for a law this library does not carry, a temperature that is not
            finite and positive, and on the reference branch an Arrhenius parameter the law
            cannot use - a negative or non-finite rate, a non-finite energy, or a
            non-positive reference temperature.

    Example:
        >>> import azoth
        >>> q = azoth.ureg.Quantity
        >>> r = kinetic_rate_law(
        ...     "legacy", q(298.15, "K"), 1e-3, q(50000.0, "J/mol"), q(298.15, "K")
        ... )
        >>> round(r.rate_factor, 12)
        0.004327147906
    """
    spec = _spec()
    checks = checks_for(spec)
    warnings: list[Warning] = []
    temperature = input_to_si(spec, "T", T)
    rate = float(reference_rate)
    energy = input_to_si(spec, "activation_energy", activation_energy)
    reference = input_to_si(spec, "reference_temperature", reference_temperature)

    apply_checks(
        checks.on_input,
        {
            "T": temperature,
            "reference_rate": rate,
            "reference_temperature": reference,
        }.get,
        warnings,
    )

    # The selector is parsed before the temperature is judged, which is the order the Rust
    # twin takes: a law that does not exist is refused whether or not the state is a state.
    folded = _LAWS.get(law.strip().lower())
    if folded is None:
        raise InvalidInputError("rate_law", f"`{law}` is neither `legacy` nor `arrhenius`")

    if not math.isfinite(temperature) or temperature <= 0.0:
        raise InvalidInputError(
            "temperature", f"{temperature} is not a finite, positive temperature"
        )

    if folded == "legacy":
        return KineticRateLawResult(
            rate_factor=LEGACY_PREFACTOR
            * math.exp(-LEGACY_ACTIVATION_TEMPERATURE / temperature)
            / LEGACY_DIVISOR,
            warnings=tuple(warnings),
        )

    if (
        not math.isfinite(rate)
        or rate < 0.0
        or not math.isfinite(energy)
        or not math.isfinite(reference)
        or reference <= 0.0
    ):
        raise InvalidInputError(
            "kinetics",
            "the reference law needs a finite nonnegative rate, a finite energy and a "
            f"positive temperature, and it has {rate}, {energy} and {reference}",
        )
    if rate == 0.0:
        return KineticRateLawResult(rate_factor=0.0, warnings=tuple(warnings))

    return KineticRateLawResult(
        rate_factor=rate * math.exp(-energy / GAS_CONSTANT * (1.0 / temperature - 1.0 / reference)),
        warnings=tuple(warnings),
    )


def _spec() -> dict[str, Any]:
    from azoth._models_gen import model

    return model("reactions.kinetic_rate_law")
