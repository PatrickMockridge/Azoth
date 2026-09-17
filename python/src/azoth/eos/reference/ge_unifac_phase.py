"""``eos.ge_unifac_phase`` - the fugacity coefficients of a UNIFAC liquid.

Spec: ``specs/models/eos/ge_unifac_phase.toml``. A *direct* model: no iteration, so no
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
from azoth.core.result import GeUnifacPhaseResult
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning
from azoth.eos.components import GeUnifacPhaseParameters
from azoth.eos.reference._ge_phase import ge_fugacities
from azoth.eos.reference.unifac_activity_coefficients import unifac_activity_coefficients

MODEL_ID = "eos.ge_unifac_phase"


def ge_unifac_phase(
    params: GeUnifacPhaseParameters,
    T: Q,
    P: Q,
    x: Sequence[float],
) -> GeUnifacPhaseResult:
    """The fugacity coefficients of a liquid whose non-ideality is UNIFAC's.

    ``phi_i = gamma_i P0_i / P``, which is NeqSim's ``ComponentGE.fugcoef`` for this
    phase. ``gamma_i`` comes from UNIFAC on ``params``' group tables at ``T`` and ``x``;
    ``P0_i`` is Antoine's correlation on the component's own coefficients at ``T``.
    ``params`` is resolved by name through
    :func:`azoth.eos.components.ge_unifac_phase_parameters`; ``x`` is checked rather than
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
        >>> from azoth.eos.components import ge_unifac_phase_parameters
        >>> q = azoth.ureg.Quantity
        >>> r = ge_unifac_phase(
        ...     ge_unifac_phase_parameters(["methanol", "water"]), q(298.15, "K"),
        ...     q(100000.0, "Pa"), [0.5, 0.5],
        ... )
        >>> round(r.gamma[0], 12)
        1.115681506247
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

    # The activity coefficients, from the same group tables
    # `eos.unifac_activity_coefficients` uses - resolved through the same function, so the
    # two cannot disagree.
    from azoth.eos.components import UnifacParameters

    activity = unifac_activity_coefficients(
        UnifacParameters(
            groups=params.groups,
            group_r=params.group_r,
            group_q=params.group_q,
            aij=params.aij,
        ),
        T,
        x,
    )
    warnings.extend(activity.warnings)
    gamma = list(activity.gamma)

    saturated = ge_fugacities(gamma, params, t_si, p_si)
    warnings.extend(saturated.warnings)

    ln_gamma = [math.log(value) for value in gamma]

    return GeUnifacPhaseResult(
        gamma=tuple(gamma),
        ln_gamma=tuple(ln_gamma),
        ln_phi=tuple(saturated.ln_phi),
        p_sat=tuple(from_si(value, "Pa") for value in saturated.p_sat),
        warnings=tuple(warnings),
    )
