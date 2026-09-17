"""``eos.ge_van_laar_acid_phase`` - the fugacity coefficients of the
water-nitric-sulfuric acid liquid.

Spec: ``specs/models/eos/ge_van_laar_acid_phase.toml``. A *direct* model: no iteration, so
no algorithm block.

This is the pure-Python reference: a second, independent expression of the same physics as
the Rust kernel, held to it by ``test_cross_impl.py``.
"""

from __future__ import annotations

import math
from collections.abc import Sequence

from azoth import _models_gen
from azoth.core.errors import InvalidInputError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import GeVanLaarAcidPhaseResult
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning
from azoth.eos.components import GeVanLaarAcidPhaseParameters, VanLaarAcidParameters
from azoth.eos.reference._ge_phase import combine, saturation
from azoth.eos.reference.nitric_sulfuric_acid_vapor_pressure import (
    nitric_sulfuric_acid_vapor_pressure,
)
from azoth.eos.reference.van_laar_acid_activity_coefficients import (
    van_laar_acid_activity_coefficients,
)

MODEL_ID = "eos.ge_van_laar_acid_phase"

#: The fugacity coefficient NeqSim gives a species the model does not cover. A *penalty*
#: rather than physics, and finite rather than infinite so a flash can still iterate while
#: it drives the component out.
NON_MODELED_COMPONENT_PENALTY = 1.0e12


def ge_van_laar_acid_phase(
    params: GeVanLaarAcidPhaseParameters,
    T: Q,
    P: Q,
    x: Sequence[float],
) -> GeVanLaarAcidPhaseResult:
    """The fugacity coefficients of the water-nitric-sulfuric acid liquid.

    ``phi_i = gamma_i P0_i / P``, where ``gamma_i`` is the Taleb Van Laar activity
    coefficient and ``P0_i`` is the acid correlation for the three modelled species and the
    database Antoine for anything else. This is the identity
    ``ComponentGEVanLaarAcid.fugcoef`` enforces as ``f_i = gamma_i x_i P0_i``, and it
    **ignores** ``referenceStateType`` - both acids are tagged ``solute``, so the inherited
    method would give them a Henry's-law coefficient instead.

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
        OutOfRangeError: if ``T`` or ``P`` is not positive, or outside the acid
            correlation's own bounds.

    Example:
        >>> import azoth
        >>> from azoth.eos.components import ge_van_laar_acid_phase_parameters
        >>> q = azoth.ureg.Quantity
        >>> r = ge_van_laar_acid_phase(
        ...     ge_van_laar_acid_phase_parameters(["water", "nitric acid"]), q(250.0, "K"),
        ...     q(100000.0, "Pa"), [0.6, 0.4],
        ... )
        >>> round(r.gamma[0], 12)
        0.157903921115
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

    # The activity coefficients, from the same acid identities and the same renormalised
    # acid basis `eos.van_laar_acid_activity_coefficients` uses.
    activity = van_laar_acid_activity_coefficients(
        VanLaarAcidParameters(acid_index=params.acid_index), T, x
    )
    warnings.extend(activity.warnings)
    gamma = list(activity.gamma)

    # `P0` per component: the acid correlation for the three modelled species, the database
    # Antoine for anything else. `saturation` evaluates all of them and the acid entries'
    # results are then discarded - their stored rows are the all-zero ones.
    antoine_p_sat, antoine_warnings = saturation(params, t_si)
    warnings.extend(antoine_warnings)
    acids = nitric_sulfuric_acid_vapor_pressure(from_si(t_si, "K"))
    warnings.extend(acids.warnings)

    acid_pa = {
        1: acids.p_water.to_base_units().magnitude,
        2: acids.p_nitric_acid.to_base_units().magnitude,
        3: acids.p_sulfuric_acid.to_base_units().magnitude,
    }
    p_sat = [acid_pa.get(index, antoine_p_sat[i]) for i, index in enumerate(params.acid_index)]

    # The penalty replaces the whole coefficient, not just `P0`: NeqSim returns it before
    # computing anything else.
    composed = combine(gamma, p_sat, p_si)
    ln_phi = [
        math.log(NON_MODELED_COMPONENT_PENALTY) if index == 0 else composed[i]
        for i, index in enumerate(params.acid_index)
    ]
    ln_gamma = [math.log(value) for value in gamma]

    return GeVanLaarAcidPhaseResult(
        gamma=tuple(gamma),
        ln_gamma=tuple(ln_gamma),
        ln_phi=tuple(ln_phi),
        p_sat=tuple(from_si(value, "Pa") for value in p_sat),
        warnings=tuple(warnings),
    )
