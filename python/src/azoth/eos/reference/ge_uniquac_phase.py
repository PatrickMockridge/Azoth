"""``eos.ge_uniquac_phase`` - the fugacity coefficients of a UNIQUAC liquid.

Spec: ``specs/models/eos/ge_uniquac_phase.toml``. A *direct* model: no iteration, so no
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
from azoth.core.result import GeUniquacPhaseResult
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning
from azoth.eos.components import GeUniquacPhaseParameters
from azoth.eos.reference._ge_phase import ge_fugacities
from azoth.eos.reference.uniquac_activity_coefficients import uniquac_activity_coefficients

MODEL_ID = "eos.ge_uniquac_phase"


def ge_uniquac_phase(
    params: GeUniquacPhaseParameters,
    T: Q,
    P: Q,
    x: Sequence[float],
    aij: Sequence[Sequence[Q]],
) -> GeUniquacPhaseResult:
    """The fugacity coefficients of a liquid whose non-ideality is UNIQUAC's.

    ``phi_i = gamma_i P0_i / P``, which is NeqSim's ``ComponentGE.fugcoef`` for this
    phase. ``gamma_i`` comes from UNIQUAC on ``params``' ``r`` and ``q`` with ``aij``;
    ``P0_i`` is Antoine's correlation on the component's own coefficients at ``T``.

    ``aij`` is the caller's - no upstream table carries a UNIQUAC interaction matrix. The
    public wrapper resolves it from a keycard when it is omitted; this reference is handed
    it.

    Args:
        params: the resolved phase parameters, by component.
        T: absolute temperature.
        P: absolute pressure.
        x: the liquid's mole fractions; non-negative and summing to one.
        aij: the interaction matrix, in kelvin, directional with a zero diagonal.

    Returns:
        The activity coefficients, their logarithms, the fugacity coefficients and the
        pure-component saturation pressures.

    Raises:
        InvalidInputError: if ``params`` and ``x`` disagree in length, if ``x`` is not a
            composition, or if ``aij`` is not an ``N x N`` matrix.
        OutOfRangeError: if ``T`` or ``P`` is not positive.

    Example:
        >>> import azoth
        >>> from azoth.eos.components import ge_uniquac_phase_parameters
        >>> q = azoth.ureg.Quantity
        >>> r = ge_uniquac_phase(
        ...     ge_uniquac_phase_parameters(["methanol", "water"]), q(298.15, "K"),
        ...     q(100000.0, "Pa"), [0.5, 0.5],
        ...     aij=[[q(0.0, "K"), q(-71.0, "K")], [q(209.0, "K"), q(0.0, "K")]],
        ... )
        >>> round(r.gamma[0], 12)
        1.218544184873
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

    if len(aij) != n or any(len(row) != n for row in aij):
        raise InvalidInputError(
            "aij",
            f"the interaction matrix is {len(aij)} x "
            f"{len(aij[0]) if aij else 0} but the mixture has {n} components",
        )
    # The activity coefficients, from the same r and q and the same interaction matrix
    # `eos.uniquac_activity_coefficients` uses - resolved through the same function, so the
    # two cannot disagree.
    from azoth.eos.components import UniquacParameters

    activity = uniquac_activity_coefficients(
        UniquacParameters(r=params.r, q=params.q), T, x, aij=aij
    )
    warnings.extend(activity.warnings)
    gamma = list(activity.gamma)

    saturated = ge_fugacities(gamma, params, t_si, p_si)
    warnings.extend(saturated.warnings)

    ln_gamma = [math.log(value) for value in gamma]

    return GeUniquacPhaseResult(
        gamma=tuple(gamma),
        ln_gamma=tuple(ln_gamma),
        ln_phi=tuple(saturated.ln_phi),
        p_sat=tuple(from_si(value, "Pa") for value in saturated.p_sat),
        warnings=tuple(warnings),
    )
