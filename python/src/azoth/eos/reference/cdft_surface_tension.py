"""``eos.cdft_surface_tension`` - NeqSim's ``CDFTSurfaceTension``.

Spec: ``specs/models/eos/cdft_surface_tension.toml``. A tanh density profile on a
2048-point grid, convolved with a step kernel whose width is one molecular diameter; the
interface width is chosen by a fixed scan and golden section, and the tension is the
integral of the excess grand potential.
"""

from __future__ import annotations

import math
from collections.abc import Callable, Sequence
from typing import Any

from azoth.core.errors import InvalidInputError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import CdftSurfaceTensionResult
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning
from azoth.eos.reference._mixture_state import reduced_parameters

CALC_ID = "eos.cdft_surface_tension"

R = 8.3144621
NA = 6.023e23
NGRID = 2048
DOMAIN_HALF_WIDTH_IN_D = 15.0
ATTRACTIVE_RANGE = 0.5
SCAN_PURE = 80
GOLDEN_PURE = 60
SCAN_MIXTURE = 60
GOLDEN_MIXTURE = 50


def _mu_repulsive_excess(rho: float, t: float, b: float) -> float:
    b_rho = min(b * rho, 0.999)
    return R * t * (-math.log(1.0 - b_rho) + b_rho / (1.0 - b_rho))


def _mu_attractive(rho: float, a: float, b: float, delta1: float, delta2: float) -> float:
    u = b * rho
    denom1 = 1.0 + delta1 * u
    denom2 = 1.0 + delta2 * u
    return a / ((delta1 - delta2) * b) * math.log(denom2 / denom1) - a * rho / (denom1 * denom2)


def _f_attractive_local(rho: float, a: float, b: float, delta1: float, delta2: float) -> float:
    u = b * rho
    return rho * a / ((delta1 - delta2) * b) * math.log((1.0 + delta2 * u) / (1.0 + delta1 * u))


def _pressure(rho: float, t: float, a: float, b: float, delta1: float, delta2: float) -> float:
    b_rho = b * rho
    p_repulsive = R * t * rho / (1.0 - b_rho)
    denominator = (1.0 + delta1 * b_rho) * (1.0 + delta2 * b_rho)
    return p_repulsive - a * rho * rho / denominator


def _step_convolve(rho: Sequence[float], n_half: int, w: float, dz: float) -> list[float]:
    n = len(rho)
    result = []
    for i in range(n):
        total = 0.0
        for k in range(-n_half, n_half + 1):
            index = min(max(i + k, 0), n - 1)
            total += rho[index]
        result.append(w * total * dz)
    return result


def _pure_sigma_for_width(
    rho_l: float,
    rho_v: float,
    t: float,
    a: float,
    b: float,
    delta1: float,
    delta2: float,
    d_mol: float,
    delta: float,
    mu_eq: float,
    p_bulk: float,
) -> float:
    half = DOMAIN_HALF_WIDTH_IN_D * d_mol
    dz = 2.0 * half / (NGRID - 1)
    n_half = math.ceil(ATTRACTIVE_RANGE * d_mol / dz)
    n_points = 2 * n_half + 1
    w = -2.0 * a / (n_points * dz)

    average = 0.5 * (rho_l + rho_v)
    difference = rho_l - rho_v
    rho = [average + 0.5 * difference * math.tanh(-(-half + i * dz) / delta) for i in range(NGRID)]
    convolution = _step_convolve(rho, n_half, w, dz)

    sigma = 0.0
    for i in range(NGRID):
        ri = max(rho[i], 1.0e-10)
        if b * ri >= 0.999:
            ri = 0.999 / b
        f_repulsive = ri * R * t * (math.log(ri) - 1.0 - math.log(1.0 - b * ri))
        f_attractive_nonlocal = 0.5 * ri * convolution[i]
        f_correction = _f_attractive_local(ri, a, b, delta1, delta2) + a * ri * ri
        sigma += (f_repulsive + f_attractive_nonlocal + f_correction - mu_eq * ri + p_bulk) * dz
    return sigma


def _mixture_sigma_for_width(
    rho_liq: Sequence[float],
    rho_vap: Sequence[float],
    t: float,
    aij: Sequence[Sequence[float]],
    b: Sequence[float],
    delta1: float,
    delta2: float,
    d_mol: float,
    delta: float,
    mu_eq: Sequence[float],
    p_bulk: float,
) -> float:
    count = len(rho_liq)
    half = DOMAIN_HALF_WIDTH_IN_D * d_mol
    dz = 2.0 * half / (NGRID - 1)
    n_half = math.ceil(ATTRACTIVE_RANGE * d_mol / dz)
    n_points = 2 * n_half + 1

    profile = [[0.0] * NGRID for _ in range(count)]
    total = [0.0] * NGRID
    for ic in range(count):
        average = 0.5 * (rho_liq[ic] + rho_vap[ic])
        difference = rho_liq[ic] - rho_vap[ic]
        for i in range(NGRID):
            z = -half + i * dz
            value = average + 0.5 * difference * math.tanh(-z / delta)
            if value < 1.0e-30:
                value = 1.0e-30
            profile[ic][i] = value
            total[i] += value
    averaged = [[0.0] * NGRID for _ in range(count)]
    for jc in range(count):
        for i in range(NGRID):
            total_row = 0.0
            for k in range(-n_half, n_half + 1):
                index = min(max(i + k, 0), NGRID - 1)
                total_row += profile[jc][index]
            averaged[jc][i] = total_row / n_points

    sigma = 0.0
    for i in range(NGRID):
        rho_total = max(total[i], 1.0e-10)
        b_mix = 0.0
        a_mix = 0.0
        for ic in range(count):
            xi = profile[ic][i] / rho_total
            b_mix += xi * b[ic]
            for jc in range(count):
                xj = profile[jc][i] / rho_total
                a_mix += xi * xj * aij[ic][jc]
        if b_mix * rho_total >= 0.999:
            rho_total = 0.999 / b_mix
        f_ideal = 0.0
        for ic in range(count):
            ri = profile[ic][i]
            if ri > 1.0e-30:
                f_ideal += ri * R * t * (math.log(ri) - 1.0)
        f_repulsive = -rho_total * R * t * math.log(1.0 - b_mix * rho_total)
        f_attractive_nonlocal = 0.0
        for ic in range(count):
            for jc in range(count):
                f_attractive_nonlocal -= aij[ic][jc] * profile[ic][i] * averaged[jc][i]
        f_correction = (
            _f_attractive_local(rho_total, a_mix, b_mix, delta1, delta2)
            + a_mix * rho_total * rho_total
        )
        mu_rho = 0.0
        for ic in range(count):
            mu_rho += mu_eq[ic] * profile[ic][i]
        sigma += (
            f_ideal + f_repulsive + f_attractive_nonlocal + f_correction - mu_rho + p_bulk
        ) * dz
    return sigma


def _variational_minimum(
    width: Callable[[float], float], d_mol: float, scan_points: int, golden_steps: int
) -> float:
    lo, hi = 0.3 * d_mol, 30.0 * d_mol
    best_delta = d_mol
    best = math.inf
    for k in range(scan_points):
        delta = lo + (hi - lo) * k / (scan_points - 1)
        sigma = width(delta)
        if sigma < best:
            best = sigma
            best_delta = delta
    ref_lo = max(lo, best_delta - 3.0 * d_mol)
    ref_hi = min(hi, best_delta + 3.0 * d_mol)
    ratio = (math.sqrt(5.0) + 1.0) / 2.0
    for _ in range(golden_steps):
        c = ref_hi - (ref_hi - ref_lo) / ratio
        d = ref_lo + (ref_hi - ref_lo) / ratio
        if width(c) < width(d):
            ref_hi = d
        else:
            ref_lo = c
    return width(0.5 * (ref_lo + ref_hi))


def _mixture_free_energy_density(
    rho: Sequence[float],
    aij: Sequence[Sequence[float]],
    b: Sequence[float],
    t: float,
    delta1: float,
    delta2: float,
) -> float:
    count = len(rho)
    rho_total = sum(rho)
    if rho_total < 1.0e-20:
        return 0.0
    b_mix = 0.0
    a_mix = 0.0
    for i in range(count):
        xi = rho[i] / rho_total
        b_mix += xi * b[i]
        for j in range(count):
            xj = rho[j] / rho_total
            a_mix += xi * xj * aij[i][j]
    if b_mix * rho_total >= 0.999:
        rho_total = 0.999 / b_mix
    f_ideal = 0.0
    for i in range(count):
        if rho[i] > 1.0e-30:
            f_ideal += rho[i] * R * t * (math.log(rho[i]) - 1.0)
    f_repulsive = -rho_total * R * t * math.log(1.0 - b_mix * rho_total)
    f_attractive = _f_attractive_local(rho_total, a_mix, b_mix, delta1, delta2)
    return f_ideal + f_repulsive + f_attractive


def _mixture_pressure(
    rho: Sequence[float],
    aij: Sequence[Sequence[float]],
    b: Sequence[float],
    t: float,
    delta1: float,
    delta2: float,
) -> float:
    count = len(rho)
    rho_total = sum(rho)
    b_mix = 0.0
    a_mix = 0.0
    for i in range(count):
        xi = rho[i] / rho_total
        b_mix += xi * b[i]
        for j in range(count):
            xj = rho[j] / rho_total
            a_mix += xi * xj * aij[i][j]
    return _pressure(rho_total, t, a_mix, b_mix, delta1, delta2)


def cdft_surface_tension(
    components: Sequence[str],
    eos: str,
    T: Q,
    x_liquid: Sequence[float],
    x_vapour: Sequence[float],
    molar_volume_liquid: Q,
    molar_volume_vapour: Q,
) -> CdftSurfaceTensionResult:
    """The classical density functional theory's interface tension between two states.

    Args:
        components: the fluid's substances, by name.
        eos: the cubic whose attraction and covolume each component carries.
        T: absolute temperature.
        x_liquid: the liquid phase's mole fractions.
        x_vapour: the vapour phase's mole fractions.
        molar_volume_liquid: the liquid phase's molar volume.
        molar_volume_vapour: the vapour phase's molar volume.

    Returns:
        The interface tension, in ``N/m``.

    Raises:
        InvalidInputError: for a composition of the wrong length or one that does not sum
            to one, or a variational minimum that is not positive.
        OutOfRangeError: if an input is outside its declared range.

    Example:
        >>> import azoth
        >>> q = azoth.ureg.Quantity
        >>> r = cdft_surface_tension(
        ...     ["methane"], "pr", q(150.0, "K"), [1.0], [1.0000000011004573],
        ...     q(4.128785631897625e-5, "m**3/mol"), q(9.709562218655069e-4, "m**3/mol"))
        >>> abs(r.sigma.magnitude - 0.002875574857984931) < 1e-12
        True
    """
    from azoth._models_gen import model
    from azoth.eos import components as databank

    spec: dict[str, Any] = model(CALC_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    count = len(components)
    if len(x_liquid) != count or len(x_vapour) != count:
        raise InvalidInputError(
            "x_liquid",
            f"{count} component(s) against {len(x_liquid)} liquid and "
            f"{len(x_vapour)} vapour mole fraction(s): each is one entry per component",
        )
    for name, composition in (("x_liquid", x_liquid), ("x_vapour", x_vapour)):
        total = sum(composition)
        if abs(total - 1.0) > 1.0e-6:
            raise InvalidInputError(name, f"the mole fractions sum to {total}, not one")

    temperature = input_to_si(spec, "T", T)
    volume_liquid = input_to_si(spec, "molar_volume_liquid", molar_volume_liquid)
    volume_vapour = input_to_si(spec, "molar_volume_vapour", molar_volume_vapour)
    liquid = [float(value) for value in x_liquid]
    vapour = [float(value) for value in x_vapour]
    apply_checks(
        checks.on_input,
        {
            "T": temperature,
            "x_liquid": min(liquid) if liquid else None,
            "x_vapour": min(vapour) if vapour else None,
            "molar_volume_liquid": volume_liquid,
            "molar_volume_vapour": volume_vapour,
        }.get,
        warnings,
    )

    mixture, _ = databank.mixture_of(list(components), eos=eos)
    # The cubic's a_i(T) and b_i are pressure-independent; the reduced parameters are read
    # at one bar and converted back, so the model needs no pressure.
    reference = 1.0e5
    reduced = reduced_parameters(mixture, temperature, reference)
    rt = R * temperature
    a = [reduced.a[i] * rt * rt / reference for i in range(count)]
    b = [reduced.b[i] * rt / reference for i in range(count)]
    if mixture.cubic.name == "srk":
        delta1, delta2 = 1.0, 0.0
    else:
        delta1, delta2 = 1.0 + math.sqrt(2.0), 1.0 - math.sqrt(2.0)

    rho_liquid = [liquid[i] / volume_liquid for i in range(count)]
    rho_vapour = [vapour[i] / volume_vapour for i in range(count)]
    if sum(rho_liquid) < sum(rho_vapour):
        rho_liquid, rho_vapour = rho_vapour, rho_liquid

    if count == 1:
        rho_l, rho_v = rho_liquid[0], rho_vapour[0]
        d_mol = (6.0 * b[0] / (math.pi * NA)) ** (1.0 / 3.0)
        mu_pure = (
            R * temperature * math.log(rho_l)
            + _mu_repulsive_excess(rho_l, temperature, b[0])
            + _mu_attractive(rho_l, a[0], b[0], delta1, delta2)
        )
        p_bulk = _pressure(rho_l, temperature, a[0], b[0], delta1, delta2)
        sigma = _variational_minimum(
            lambda delta: _pure_sigma_for_width(
                rho_l,
                rho_v,
                temperature,
                a[0],
                b[0],
                delta1,
                delta2,
                d_mol,
                delta,
                mu_pure,
                p_bulk,
            ),
            d_mol,
            SCAN_PURE,
            GOLDEN_PURE,
        )
    else:
        aij = [[math.sqrt(a[i] * a[j]) for j in range(count)] for i in range(count)]
        total = sum(rho_liquid)
        b_average = sum(rho_liquid[i] / total * b[i] for i in range(count))
        d_mol = (6.0 * b_average / (math.pi * NA)) ** (1.0 / 3.0)
        eps = total * 1.0e-7
        mu_eq: list[float] = []
        for i in range(count):
            plus = list(rho_liquid)
            minus = list(rho_liquid)
            plus[i] += eps
            minus[i] -= eps
            f_plus = _mixture_free_energy_density(plus, aij, b, temperature, delta1, delta2)
            f_minus = _mixture_free_energy_density(minus, aij, b, temperature, delta1, delta2)
            mu_eq.append((f_plus - f_minus) / (2.0 * eps))
        p_bulk = _mixture_pressure(rho_liquid, aij, b, temperature, delta1, delta2)
        sigma = _variational_minimum(
            lambda delta: _mixture_sigma_for_width(
                rho_liquid,
                rho_vapour,
                temperature,
                aij,
                b,
                delta1,
                delta2,
                d_mol,
                delta,
                mu_eq,
                p_bulk,
            ),
            d_mol,
            SCAN_MIXTURE,
            GOLDEN_MIXTURE,
        )

    if sigma <= 0.0:
        raise InvalidInputError(
            "molar_volume_liquid",
            f"the variational search's minimum is {sigma:.6e} N/m, so the excess grand "
            f"potential never rises above the bulk; the class clamps the answer to 0.0, and "
            f"this port refuses a state it cannot give a tension to",
        )
    return CdftSurfaceTensionResult(sigma=from_si(sigma, "N/m"), warnings=tuple(warnings))
