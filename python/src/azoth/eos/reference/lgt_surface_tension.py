"""``eos.lgt_surface_tension`` - the linear-gradient-theory surface tension.

Spec: ``specs/models/eos/lgt_surface_tension.toml``, which records NeqSim's
``LGTSurfaceTension.calcSurfaceTension``.
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from typing import Any

from azoth.core.errors import InvalidInputError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import LgtSurfaceTensionResult
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning
from azoth.eos.reference._mixture_state import (
    helmholtz_gradient,
    reduced_parameters,
)
from azoth.eos.reference.pt_flash import pt_flash

CALC_ID = "eos.lgt_surface_tension"

#: NeqSim's gas constant.
R = 8.3144621
#: The class's fixed number of steps along the reference component's density.
STEPS = 200


def lgt_surface_tension(
    components: Sequence[str],
    eos: str,
    influence_parameters: Sequence[float],
    reference_component: float,
    pair_interactions: Sequence[Sequence[float]],
    T: Q,
    P: Q,
    z: Sequence[float],
) -> LgtSurfaceTensionResult:
    """The linear-gradient-theory interface tension of a flashed fluid.

    The density path between the two phases is a straight line in molar density over
    ``STEPS`` steps; each step evaluates the T,V chemical potential
    ``mu_i/(RT) = ln(rho_i R T) + dA_R/dn_i``, and the tension is the sum of
    ``sqrt(2 kappa mu_times_rho)`` along the path. The path crosses the cubic's own
    loop, where the pressure is negative; the T,V form is finite there.

    Args:
        components: the fluid's substances, by name.
        eos: the cubic the fluid runs - ``pr``, ``srk`` or ``rk``.
        influence_parameters: each component's gradient-theory influence parameter.
        reference_component: the 0-based index of the highest-boiling component.
        pair_interactions: the class's water/MEG ``interact`` matrix.
        T: absolute temperature.
        P: absolute pressure.
        z: the overall mole fractions.

    Returns:
        The interface tension, in ``N/m``.

    Raises:
        InvalidInputError: for a single-phase feed, an absent component, or a step whose
            integrand is not positive.
        OutOfRangeError: if an input is outside its declared range.

    Example:
        >>> import azoth
        >>> q = azoth.ureg.Quantity
        >>> r = lgt_surface_tension(
        ...     ["methane", "n-butane"], "srk",
        ...     [9.634812734820117e-21, 2.7571945977326497e-19], 1.0,
        ...     [[0.0, 0.0], [0.0, 0.0]], q(300.0, "K"), q(20.0, "bar"),
        ...     [0.5, 0.5])
        >>> abs(r.sigma.magnitude - 0.010011889708138906) < 1e-6
        True
    """
    from azoth._models_gen import model
    from azoth.eos import components as databank

    spec: dict[str, Any] = model(CALC_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    count = len(z)
    if (
        len(components) != count
        or len(influence_parameters) != count
        or len(pair_interactions) != count
        or any(len(row) != count for row in pair_interactions)
    ):
        raise InvalidInputError(
            "influence_parameters",
            f"{count} component(s) against {len(influence_parameters)} influence "
            f"parameter(s), {len(pair_interactions)} interaction row(s) and "
            f"{len(components)} name(s): every one is one entry per component",
        )
    if (
        not math.isfinite(reference_component)
        or reference_component < 0.0
        or reference_component != int(reference_component)
        or int(reference_component) >= count
    ):
        raise InvalidInputError(
            "reference_component",
            f"{reference_component} is not a whole index in 0..{count}: it is the "
            f"component the profile integrates along",
        )
    reference = int(reference_component)

    temperature = input_to_si(spec, "T", T)
    pressure = input_to_si(spec, "P", P)
    influences = [float(value) for value in influence_parameters]
    apply_checks(
        checks.on_input,
        {
            "T": temperature,
            "P": pressure,
            "influence_parameters": min(influences) if influences else None,
            "z": min(float(value) for value in z) if z else None,
        }.get,
        warnings,
    )

    mixture, _ = databank.mixture_of(list(components), eos=eos)
    flash = pt_flash(
        mixture, from_si(temperature, "K"), from_si(pressure, "Pa"), [float(v) for v in z]
    )
    if flash.phase != "two_phase":
        raise InvalidInputError(
            "z",
            "the feed is single phase at this state, so there is no interface to have a "
            "tension; a gradient-theory profile runs between two phases",
        )
    if flash.z_vapour <= 0.0 or flash.z_liquid <= 0.0:
        raise InvalidInputError(
            "z", "a phase's compressibility root is not positive, so it has no density"
        )
    for fraction in list(flash.y) + list(flash.x):
        if fraction <= 0.0:
            raise InvalidInputError(
                "z",
                "a component is absent from one of the phases, so its chemical potential "
                "would be `ln(0)`; this port refuses the partial integral the class would "
                "return at the resulting NaN",
            )

    rt = R * temperature
    rho_vapour = [flash.y[i] * pressure / (flash.z_vapour * rt) for i in range(count)]
    rho_liquid = [flash.x[i] * pressure / (flash.z_liquid * rt) for i in range(count)]
    reduced = reduced_parameters(mixture, temperature, pressure)
    kij = mixture.kij
    gradient_vapour = helmholtz_gradient(reduced, kij, list(flash.y), flash.z_vapour)
    mu_equi = [rt * (math.log(rho_vapour[i] * rt) + gradient_vapour[i]) for i in range(count)]

    a_abs = [reduced.a[i] * rt * rt / pressure for i in range(count)]
    b_abs = [reduced.b[i] * rt / pressure for i in range(count)]
    delta1 = mixture.cubic.delta1
    delta2 = mixture.cubic.delta2

    delta_rho = [(rho_liquid[i] - rho_vapour[i]) / STEPS for i in range(count)]
    if delta_rho[reference] == 0.0:
        raise InvalidInputError(
            "reference_component",
            "the reference component's density is the same in both phases, so the profile "
            "has no integration variable",
        )

    sigma = 0.0
    for step in range(1, STEPS):
        fraction = step / STEPS
        rho = [rho_vapour[i] + (rho_liquid[i] - rho_vapour[i]) * fraction for i in range(count)]
        total = sum(rho)
        x = [value / total for value in rho]
        volume = 1.0 / total
        b_mix = sum(x[i] * b_abs[i] for i in range(count))
        a_mix = sum(
            x[i] * x[j] * (1.0 - kij[i][j]) * math.sqrt(a_abs[i] * a_abs[j])
            for i in range(count)
            for j in range(count)
        )
        p_step = rt / (volume - b_mix) - a_mix / (
            (volume + delta1 * b_mix) * (volume + delta2 * b_mix)
        )
        compressibility = p_step * volume / rt
        ratio = p_step / pressure
        reduced_step = reduced._replace(
            a=[value * ratio for value in reduced.a],
            b=[value * ratio for value in reduced.b],
        )
        gradient = helmholtz_gradient(reduced_step, kij, x, compressibility)
        mu = [rt * (math.log(rho[i] * rt) + gradient[i]) for i in range(count)]

        kappa = 0.0
        for i in range(count):
            ki = 1.0 if i == reference else delta_rho[i] / delta_rho[reference]
            for k in range(count):
                kk = 1.0 if k == reference else delta_rho[k] / delta_rho[reference]
                kappa += (
                    math.sqrt(influences[i] * influences[k])
                    * ki
                    * kk
                    * (1.0 - float(pair_interactions[i][k]))
                )
        mu_times_rho = sum(rho[i] * (mu[i] - mu_equi[i]) for i in range(count)) - (
            p_step - pressure
        )
        if kappa <= 0.0 or mu_times_rho <= 0.0:
            raise InvalidInputError(
                "z",
                f"the gradient-theory integrand is not positive at step {step} "
                f"(kappa = {kappa}, mu_times_rho = {mu_times_rho}); the class takes its "
                f"square root and stops at the NaN, and this port refuses a partial "
                f"integral",
            )
        sigma += math.sqrt(2.0 * kappa * mu_times_rho) * delta_rho[reference]

    return LgtSurfaceTensionResult(sigma=from_si(abs(sigma), "N/m"), warnings=tuple(warnings))
