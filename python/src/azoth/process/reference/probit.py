"""``process.probit`` - the dose-response curve, as NeqSim's ``ProbitModel`` computes it.

Spec: ``specs/models/process/probit.toml``

The Python twin of ``crates/azoth-process/src/models/probit.rs``: `Y = a + b n ln(D)`, the
standard normal CDF over Abramowitz & Stegun 7.1.26, and the five factory triples.
"""

from __future__ import annotations

import math

from azoth.core.errors import InvalidInputError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import ProbitResult
from azoth.core.warnings import Warning

MODEL_ID = "process.probit"


def _erf(x: float) -> float:
    """Abramowitz & Stegun 7.1.26, the class's approximation."""
    sign = -1.0 if x < 0.0 else 1.0
    x = abs(x)
    t = 1.0 / (1.0 + 0.3275911 * x)
    y = 1.0 - (
        (((1.061405429 * t - 1.453152027) * t + 1.421413741) * t - 0.284496736) * t + 0.254829592
    ) * t * math.exp(-x * x)
    return sign * y


def _probability(probit: float) -> float:
    """The standard normal CDF over :func:`_erf`."""
    return 0.5 * (1.0 + _erf((probit - 5.0) / math.sqrt(2.0)))


def _spec() -> dict[str, object]:
    """The generated spec, imported here so the module has no import-time cycle."""
    from azoth._models_gen import model

    return model(MODEL_ID)


def probit(
    a: float,
    b: float,
    n: float,
    dose: float | None = None,
    probit_value: float | None = None,
) -> ProbitResult:
    """The probit of a hazard dose, and the probability a probit value implies.

    Args:
        a: the probit intercept.
        b: the probit slope.
        n: the dose exponent, ``1.0`` where the published ``(a, b)`` absorb it.
        dose: the hazard's dose in its own units; ``probit`` is published from it.
        probit_value: a probit ``Y`` stated directly, the class's ``probability(Y)``
            argument. When stated it is the value ``probability`` is evaluated at;
            otherwise the dose's probit is.

    Returns:
        The probit, absent without a dose, and the probability.

    Raises:
        InvalidInputError: when neither a dose nor a probit value is stated, or when the
            stated dose is not positive.
    """
    spec = _spec()
    checks = checks_for(spec)
    warnings: list[Warning] = []

    duty = None if dose is None else float(dose)
    stated = None if probit_value is None else float(probit_value)
    apply_checks(checks.on_input, {"dose": duty}.get, warnings)
    if duty is None and stated is None:
        raise InvalidInputError(
            "dose", "a probit needs a dose or a probit value, and neither was stated"
        )
    if math.isnan(a) or math.isnan(b) or math.isnan(n):
        raise InvalidInputError("a", "the probit constants must be numbers, and a NaN is not")
    if stated is not None and math.isnan(stated):
        raise InvalidInputError(
            "probit_value", "a stated probit value must be a number, and a NaN is not"
        )

    computed = None if duty is None else a + b * n * math.log(duty)
    evaluated = stated if stated is not None else computed
    return ProbitResult(
        probit=computed,
        probability=_probability(evaluated if evaluated is not None else math.nan),
        warnings=tuple(warnings),
    )
