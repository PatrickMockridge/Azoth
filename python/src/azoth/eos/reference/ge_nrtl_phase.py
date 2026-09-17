"""``eos.ge_nrtl_phase`` - the fugacity coefficients of an NRTL activity-coefficient
liquid.

Spec: ``specs/models/eos/ge_nrtl_phase.toml``. A *direct* model: no iteration, so no
algorithm block.

This is the pure-Python reference: a second, independent expression of the same
physics as the Rust kernel, held to it by ``test_cross_impl.py``.
"""

from __future__ import annotations

import math
from collections.abc import Sequence

from azoth import _models_gen
from azoth.core.errors import InvalidInputError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import GeNrtlPhaseResult
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning
from azoth.eos.components import GeNrtlPhaseParameters
from azoth.eos.reference.antoine_vapor_pressure import antoine_vapor_pressure
from azoth.eos.reference.nrtl_activity_coefficients import nrtl_activity_coefficients

MODEL_ID = "eos.ge_nrtl_phase"

#: How many Antoine coefficients each component carries. Named because both the
#: resolver and this reader slice the flat vector by it.
ANTOINE_COEFFICIENTS = 5


def ge_nrtl_phase(
    params: GeNrtlPhaseParameters,
    T: Q,
    P: Q,
    x: Sequence[float],
) -> GeNrtlPhaseResult:
    """The fugacity coefficients of a liquid whose non-ideality is NRTL's.

    ``phi_i = gamma_i P0_i / P``, which is NeqSim's ``ComponentGE.fugcoef`` for this
    phase. ``gamma_i`` comes from NRTL on ``params`` at ``T`` and ``x``; ``P0_i`` is
    Antoine's correlation on the component's own coefficients at ``T``. ``params`` is
    resolved by name through
    :func:`azoth.eos.components.ge_nrtl_phase_parameters`; ``x`` is checked rather than
    renormalised.

    Args:
        params: the resolved phase parameters, by component.
        T: absolute temperature.
        P: absolute pressure.
        x: the liquid's mole fractions; non-negative and summing to one.

    Returns:
        The activity coefficients, their logarithms, the fugacity coefficients and the
        pure-component saturation pressures.

    Raises:
        InvalidInputError: if ``params`` and ``x`` disagree in length, or ``x`` is not a
            composition.
        OutOfRangeError: if ``T`` or ``P`` is not positive.

    Example:
        >>> import azoth
        >>> from azoth.eos.components import ge_nrtl_phase_parameters
        >>> q = azoth.ureg.Quantity
        >>> r = ge_nrtl_phase(
        ...     ge_nrtl_phase_parameters(["methanol", "water"]), q(298.15, "K"),
        ...     q(100000.0, "Pa"), [0.5, 0.5],
        ... )
        >>> round(r.gamma[0], 12)
        1.233178856168
    """
    spec = _models_gen.model(MODEL_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    t_si = input_to_si(spec, "T", T)
    p_si = input_to_si(spec, "P", P)
    apply_checks(checks.on_input, {"T": t_si, "P": p_si}.get, warnings)

    x = list(x)
    n = len(x)
    if n == 0:
        raise InvalidInputError("x", "a mixture of zero components has no fugacity coefficient")
    if len(params.antoine_type) != n:
        raise InvalidInputError(
            "components",
            f"the phase has {len(params.antoine_type)} component(s) but `x` has {n} entries",
        )
    if any(value < 0.0 for value in x):
        bad = next(i for i, value in enumerate(x) if value < 0.0)
        raise InvalidInputError("x", f"x[{bad}] is {x[bad]} but a mole fraction cannot be negative")
    if abs(sum(x) - 1.0) > 1.0e-09:
        raise InvalidInputError(
            "x",
            f"the mole fractions sum to {sum(x)}, not to one. Renormalising them here "
            "would make a composition error invisible in every number downstream, so "
            "it is refused instead",
        )

    # The activity coefficients, from the same matrices `eos.nrtl_activity_coefficients`
    # uses - resolved through the same function, so the two cannot disagree.
    from azoth.eos.components import NrtlParameters

    gamma = list(
        nrtl_activity_coefficients(NrtlParameters(alpha=params.alpha, dij=params.dij), T, x).gamma
    )

    # The pure-component saturation pressures, one Antoine evaluation each.
    p_sat: list[float] = []
    for i in range(n):
        start = i * ANTOINE_COEFFICIENTS
        [a, b, c, d, e] = params.antoine_coefficients[start : start + ANTOINE_COEFFICIENTS]
        saturated = antoine_vapor_pressure(
            a,
            b,
            c,
            d,
            e,
            params.antoine_type[i],
            from_si(params.antoine_tc[i], "K"),
            from_si(params.antoine_pc[i], "Pa"),
            T,
        )
        warnings.extend(saturated.warnings)
        p_sat.append(saturated.p_sat.to_base_units().magnitude)

    ln_gamma = [math.log(value) for value in gamma]
    ln_phi = [lg + math.log(p0 / p_si) for lg, p0 in zip(ln_gamma, p_sat, strict=True)]

    return GeNrtlPhaseResult(
        gamma=tuple(gamma),
        ln_gamma=tuple(ln_gamma),
        ln_phi=tuple(ln_phi),
        p_sat=tuple(from_si(value, "Pa") for value in p_sat),
        warnings=tuple(warnings),
    )
