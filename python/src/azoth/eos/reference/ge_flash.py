"""``eos.ge_nrtl_flash`` - the isothermal flash of a cubic vapour over an NRTL liquid.

Spec: ``specs/models/eos/ge_nrtl_flash.toml``

NeqSim's ``SystemNRTL`` pairs ``PhaseSrkEos`` with ``PhaseGENRTL``, and the K-value its
gamma-phi iteration updates is ``K_i = phi_i^L / phi_i^V`` - see ``TPflash``'s
``sucsSubsGammaPhi``. For an activity-coefficient liquid ``phi_i^L`` is
``gamma_i P0_i / P``, which is :mod:`azoth.eos.reference.ge_nrtl_phase` exactly, so the
liquid half of this model is that model and the vapour half is the cubic.

The scaffold - Wilson K-values, the Rachford-Rice bracket and bisection, the
convergence measure, what a trivial solution is - is
:mod:`azoth.eos.reference._mixture_state`'s, shared with
:mod:`azoth.eos.reference.pt_flash`. What is here is the one line that differs.

This is the pure-Python reference: a second, independent expression of the same
physics as the Rust kernel, held to it by ``test_cross_impl.py``.
"""

from __future__ import annotations

import math
from collections.abc import Callable, Sequence

from azoth import _models_gen
from azoth.core.errors import InvalidInputError, SolverNotConvergedError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import GeFlashResult, Phase
from azoth.core.units import Q, input_to_si
from azoth.core.warnings import Warning
from azoth.eos.components import (
    ge_nrtl_phase_parameters,
    ge_unifac_phase_parameters,
    ge_van_laar_acid_phase_parameters,
    ge_wilson_phase_parameters,
    mixture_of,
)
from azoth.eos.reference._mixture_state import (
    compositions,
    is_trivial,
    negative_flash_warning,
    phase_state,
    rachford_rice,
    rachford_rice_bounds,
    reduced_parameters,
    rms_delta,
    single_phase_warning,
    trivial_warning,
    wilson_k,
)
from azoth.eos.reference.ge_nrtl_phase import ge_nrtl_phase
from azoth.eos.reference.ge_unifac_phase import ge_unifac_phase
from azoth.eos.reference.ge_van_laar_acid_phase import ge_van_laar_acid_phase
from azoth.eos.reference.ge_wilson_phase import ge_wilson_phase

MODEL_ID = "eos.ge_flash"

#: The ``ln K`` below which the iteration has found the trivial solution.
TRIVIAL_TOLERANCE = 1.0e-08


def ge_flash(
    components: Sequence[str],
    cubic: str,
    liquid_model: str,
    T: Q,
    P: Q,
    z: list[float],
) -> GeFlashResult:
    """The isothermal flash of a mixture whose liquid is an NRTL phase.

    ``K_i = phi_i^L / phi_i^V``, where ``phi_i^L = gamma_i P0_i / P`` comes from
    :func:`azoth.eos.reference.ge_nrtl_phase` and ``phi_i^V`` from the cubic the
    mixture carries.

    Args:
        params: the resolved NRTL phase parameters, by component.
        cubic: the cubic the vapour is, by name. A gamma-phi flash is defined by the
            pair, so this is part of the model rather than the caller's convenience.
        T: absolute temperature.
        P: absolute pressure.
        z: the feed's mole fractions. Checked rather than renormalised.

    Returns:
        The vapour fraction, both compositions, the K-values, both sets of fugacity
        coefficients, the vapour root and the convergence report.

    Raises:
        InvalidInputError: if ``params``, ``mixture`` and ``z`` disagree in length, or
            ``z`` is not a composition.
        OutOfRangeError: if ``T`` or ``P`` is not positive.
        SolverNotConvergedError: if the iteration hits its cap or a K-value becomes
            non-finite.

    Example:
        >>> import azoth
        >>> from azoth.eos.components import mixture_of
        >>> q = azoth.ureg.Quantity
        >>> names = ["CO2", "water", "nitric acid", "sulfuric acid"]
        >>> r = ge_flash(
        ...     names,
        ...     "srk",
        ...     "van_laar_acid",
        ...     q(273.15, "K"),
        ...     q(1.0e5, "Pa"),
        ...     [0.9090909090909091, 0.06363636363636363, 0.013636363636363636,
        ...      0.013636363636363636],
        ... )
        >>> r.phase.value
        'two_phase'
    """
    spec = _models_gen.model(MODEL_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    t_si = input_to_si(spec, "T", T)
    p_si = input_to_si(spec, "P", P)
    apply_checks(checks.on_input, {"T": t_si, "P": p_si}.get, warnings)

    names = list(components)
    mixture, _ = mixture_of(names, eos=cubic)
    n = len(mixture.components)
    # The liquid's own parameters, resolved by whichever model is named. A model this
    # library does not carry is refused rather than defaulted.
    liquid_ln_phi: Callable[[list[float]], Sequence[float]]
    if liquid_model == "nrtl":
        nrtl = ge_nrtl_phase_parameters(names)
        liquid_ln_phi = lambda x: ge_nrtl_phase(nrtl, T, P, x).ln_phi  # noqa: E731
    elif liquid_model == "unifac":
        unifac = ge_unifac_phase_parameters(names)
        liquid_ln_phi = lambda x: ge_unifac_phase(unifac, T, P, x).ln_phi  # noqa: E731
    elif liquid_model == "wilson":
        wilson = ge_wilson_phase_parameters(names)
        liquid_ln_phi = lambda x: ge_wilson_phase(wilson, mixture, T, P, x).ln_phi  # noqa: E731
    elif liquid_model == "van_laar_acid":
        acid = ge_van_laar_acid_phase_parameters(names)
        liquid_ln_phi = lambda x: ge_van_laar_acid_phase(acid, T, P, x).ln_phi  # noqa: E731
    else:
        raise InvalidInputError(
            "liquid_model",
            f"`{liquid_model}` is not an activity-coefficient phase this library carries; "
            "the four are nrtl, unifac, wilson and van_laar_acid",
        )
    z = list(z)
    if len(z) != n:
        raise InvalidInputError("z", f"a feed for {n} components has {len(z)} entries")
    if any(value < 0.0 for value in z):
        bad = next(i for i, value in enumerate(z) if value < 0.0)
        raise InvalidInputError("z", f"z[{bad}] is {z[bad]} but a mole fraction cannot be negative")
    if abs(sum(z) - 1.0) > 1.0e-09:
        raise InvalidInputError(
            "z",
            f"the feed's mole fractions sum to {sum(z)}, not to one. Renormalising them "
            "here would make a composition error invisible in every number downstream, "
            "so it is refused instead",
        )

    min_t_over_tc = min(
        t_si / component.Tc.to_base_units().magnitude for component in mixture.components
    )
    apply_checks(checks.derived, {"min_t_over_tc": min_t_over_tc}.get, warnings)

    algorithm = spec["algorithm"]
    inner = algorithm["inner"]
    tolerance = float(algorithm["tolerance"])

    reduced = reduced_parameters(mixture, t_si, p_si)
    warnings.extend(reduced.warnings)
    kij = mixture.kij

    k = wilson_k(mixture, t_si, p_si)
    iterations = 0
    residual = math.nan
    settled: str | None = None

    for step in range(1, int(algorithm["max_iterations"]) + 1):
        iterations = step

        # Both guards run *before* the Rachford-Rice solve, because both are cases
        # where the bracket is a division by zero: `K_i = 1` puts a pole at infinity,
        # and K-values all on one side of one put a bracket end there.
        if is_trivial(k, TRIVIAL_TOLERANCE):
            settled = "trivial"
            break
        # A K-value that is not finite and positive is the iteration diverging rather
        # than a state to diagnose: every path below treats the K-values as a
        # composition ratio, and `nan > 1.0` being false would otherwise report a
        # diverged iteration as a single-phase feed.
        if any(not math.isfinite(value) or value <= 0.0 for value in k):
            raise SolverNotConvergedError(iterations, residual, tolerance)
        bounds = rachford_rice_bounds(k)
        if bounds is None:
            settled = "all_vapour" if all(value > 1.0 for value in k) else "all_liquid"
            break

        beta = rachford_rice(z, k, bounds, inner["tolerance"], inner["max_iterations"])
        x, y = compositions(z, k, beta)

        # The one line that makes this a gamma-phi flash rather than `eos.pt_flash`:
        # the liquid's coefficients are the activity model's and the vapour's are the
        # cubic's.
        liquid_ln_phi_values = liquid_ln_phi(x)
        vapour = phase_state(reduced, kij, y, liquid=False)

        ln_k_new = [
            left - right for left, right in zip(liquid_ln_phi_values, vapour.ln_phi, strict=True)
        ]
        residual = rms_delta(ln_k_new, k)
        k = [math.exp(value) for value in ln_k_new]

        if residual <= tolerance:
            break

    beta_out: float | None
    if settled == "trivial" or (settled is None and is_trivial(k, TRIVIAL_TOLERANCE)):
        # Stated *exactly* - x = y = z and K = 1 - rather than as the last iterate
        # that approached it, which is a function of where the bisection stopped on an
        # identically-zero function and is not reproducible.
        k = [1.0] * n
        x, y = list(z), list(z)
        phase = Phase.TRIVIAL
        beta_out = None
        warnings.append(trivial_warning())
    elif settled is not None:
        x, y = list(z), list(z)
        phase = Phase.ALL_VAPOUR if settled == "all_vapour" else Phase.ALL_LIQUID
        beta_out = None
        warnings.append(single_phase_warning(phase))
    elif residual > tolerance:
        raise SolverNotConvergedError(iterations, residual, tolerance)
    else:
        bounds = rachford_rice_bounds(k)
        if bounds is None:  # pragma: no cover - guarded above
            raise SolverNotConvergedError(iterations, residual, tolerance)
        beta_out = rachford_rice(z, k, bounds, inner["tolerance"], inner["max_iterations"])
        x, y = compositions(z, k, beta_out)
        if beta_out < 0.0:
            phase = Phase.ALL_LIQUID
        elif beta_out > 1.0:
            phase = Phase.ALL_VAPOUR
        else:
            phase = Phase.TWO_PHASE
        if phase is not Phase.TWO_PHASE:
            warnings.append(negative_flash_warning(phase, beta_out))

    # The one state the result reports, at the compositions actually returned.
    liquid_ln_phi_values = liquid_ln_phi(x)
    vapour = phase_state(reduced, kij, y, liquid=False)

    return GeFlashResult(
        beta=beta_out,
        x=tuple(x),
        y=tuple(y),
        k=tuple(k),
        ln_phi_liquid=tuple(liquid_ln_phi_values),
        ln_phi_vapour=tuple(vapour.ln_phi),
        z_vapour=vapour.z,
        min_t_over_tc=min_t_over_tc,
        phase=phase,
        iterations=iterations,
        residual=residual,
        warnings=tuple(warnings),
    )
