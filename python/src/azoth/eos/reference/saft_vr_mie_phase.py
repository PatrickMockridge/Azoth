"""``eos.saft_vr_mie_phase`` - the SAFT-VR-Mie phase state, the pure-Python reference.

The second, independent expression of the physics the Rust ``azoth_eos::saft_vr_mie_phase``
computes: Lafitte 2013's statistical associating fluid theory for chain molecules formed
from Mie segments, without its association term.

Spec: ``specs/models/eos/saft_vr_mie_phase.toml``, which carries why the exponents are not
the absence marker, why the chain term's contact value is a weighted geometric mean, and
why the last digits of anything here are the model's own differences.

**The numbers this is checked against are NeqSim's.**
``validation/neqsim/SaftVrMieProbe.java`` prints ``PhaseSAFTVRMie``'s intermediates one per
line - the parameter set, the effective diameter, the chain RDF's own layers, the three
dispersion terms and the state - and
``crates/azoth-eos/tests/saft_vr_mie.rs`` reproduces each of them, so a disagreement
between the kernels here is a disagreement about a layer rather than about a total.
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from typing import NamedTuple

from azoth.core.errors import InvalidInputError, OutOfRangeError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import SaftVrMiePhaseResult
from azoth.core.units import Q, input_to_si, ureg
from azoth.core.warnings import Warning
from azoth.eos.components import entry

MODEL_ID = "eos.saft_vr_mie_phase"

#: Avogadro's constant, **NeqSim's value**, in 1/mol - the packing fraction is built on it.
AVOGADRO = 6.023e23

#: NeqSim's ``ThermodynamicConstantsInterface.R``.
R = 8.3144621

#: The 10-point Gauss-Legendre nodes on ``[0, 1]``, **NeqSim's own rounded literals**.
_GAUSS_NODES = (
    0.01304673574,
    0.06746831665,
    0.16029521585,
    0.28330230294,
    0.42556283050,
    0.57443716950,
    0.71669769706,
    0.83970478415,
    0.93253168335,
    0.98695326426,
)

#: The matching weights. **Ten points carry this integrand poorly and that is the model's
#: choice**: measured against a converged integral the weak branch is ``1e-4`` out at
#: ``theta = 0.9`` and ``4e-3`` at ``0.05``. Reproducing NeqSim means reproducing the rule.
_GAUSS_WEIGHTS = (
    0.03333567215,
    0.07472567458,
    0.10954318126,
    0.13463335965,
    0.14776211236,
    0.14776211236,
    0.13463335965,
    0.10954318126,
    0.07472567458,
    0.03333567215,
)

#: Lafitte 2013's effective-packing-fraction coefficients, ``[coefficient][power of 1/lambda]``.
_ETA_EFF_COEFFS = (
    (0.81096, 1.7888, -37.578, 92.284),
    (1.0205, -19.341, 151.26, -463.50),
    (-1.9057, 22.845, -228.14, 973.92),
    (1.08850, -6.1962, 106.98, -677.64),
)

#: Lafitte 2013's Pade coefficients, ``[coefficient][function]``.
_PHI_PADE = (
    (7.5365557, -359.440, 1550.9, -1.199320, -1911.2800, 9236.9),
    (-37.604630, 1825.60, -5070.1, 9.063632, 21390.175, -129430.0),
    (71.745953, -3168.00, 6534.6, -17.94820, -51320.700, 357230.0),
    (-46.835520, 1884.20, -3288.7, 11.34027, 37064.540, -315530.0),
    (-2.4679820, -0.82376, -2.7171, 20.52142, 1103.7420, 1390.2),
    (-0.5027200, -3.19350, 2.0883, -56.63770, -3264.6100, -4518.2),
    (8.0956883, 3.70900, 0.0000, 40.53683, 2556.1810, 4241.6),
)


class _Component(NamedTuple):
    """One component's SAFT-VR-Mie set, resolved from the databank."""

    m: float
    lambda_r: float
    lambda_a: float
    sigma: float
    epsik: float

    @classmethod
    def of(cls, name: str) -> _Component:
        record = entry(name)
        return cls(
            record.m_mie,
            record.lambda_r_mie,
            record.lambda_a_mie,
            record.sigma_mie,
            record.epsik_mie,
        )

    def has_parameters(self) -> bool:
        """Whether the table gave this component a set.

        **``m`` is the marker and the exponents are not.** The table carries the standard
        ``12``/``6`` on all 286 rows, so ``lambda_r > 0`` says nothing; ``m``, ``sigma`` and
        ``epsik`` are zero together on the 274 without one.
        """
        return self.m > 0.0 and self.sigma > 0.0 and self.epsik > 0.0


def _exp(x: float) -> float:
    """``exp``, saturating to infinity as Java and Rust do.

    **Not a nicety: it is what makes this the same function the other kernel computes.**
    ``Math.exp`` of a large number is ``Infinity``, and Rust's ``f64::exp`` the same, while
    Python raises ``OverflowError`` - so a trial volume the solve walks through on its way
    to a root, where the packing fraction is past anything the model admits, would abort
    here and be carried past there. The comparison would be between a model and an
    exception.
    """
    try:
        return math.exp(x)
    except OverflowError:
        return math.inf


def _log(x: float) -> float:
    """``log``, giving NaN outside its domain as Java and Rust do. See [`_exp`]."""
    try:
        return math.log(x)
    except ValueError:
        return math.nan


def mie_prefactor(lambda_r: float, lambda_a: float) -> float:
    """``C = lambda_r/(lambda_r - lambda_a) (lambda_r/lambda_a)^(lambda_a/(lambda_r-lambda_a))``."""
    return (
        lambda_r
        / (lambda_r - lambda_a)
        * math.pow(lambda_r / lambda_a, lambda_a / (lambda_r - lambda_a))
    )


def barker_henderson(theta: float, lambda_r: float, lambda_a: float) -> float:
    """The Barker-Henderson integral, by NeqSim's two branches.

    For ``theta <= 1`` the ten points carry the whole interval; above it the well is deep
    enough that ``exp(-u*)`` underflows near the origin, so NeqSim treats a clamped
    ``[0, x_cut]`` as exactly ``x_cut`` and puts the same ten points on the rest. Both the
    split and the rule's coarseness are the model's.
    """

    def reduced(x: float) -> float:
        return theta * (math.pow(x, -lambda_r) - math.pow(x, -lambda_a))

    if theta <= 1.0:
        total = 0.0
        for node, weight in zip(_GAUSS_NODES, _GAUSS_WEIGHTS, strict=True):
            if node < 1.0e-20:
                total += weight
                continue
            total += weight * (1.0 - math.exp(-reduced(node)))
        return total

    x_cut = min(max(math.exp(-math.log(theta + 20.0) / lambda_r), 0.5), 0.999)
    low, high = x_cut, 1.0
    half, middle = (high - low) / 2.0, (high + low) / 2.0
    total = x_cut
    for node, weight in zip(_GAUSS_NODES, _GAUSS_WEIGHTS, strict=True):
        x = middle + half * (2.0 * node - 1.0)
        if x < 1.0e-20:
            total += weight * (high - low)
            continue
        total += weight * (1.0 - math.exp(-reduced(x))) * (high - low)
    return total


def barker_henderson_d_theta(theta: float, lambda_r: float, lambda_a: float) -> float:
    """``d/dtheta`` of :func:`barker_henderson`, on the same branch.

    **The derivative of the approximation rather than of the integral**, which is what makes
    it the analytic companion of NeqSim's own central difference: the cut-off ``x_cut``
    moves with ``theta`` and that motion is part of both.
    """

    def powers(x: float) -> tuple[float, float]:
        inv = 1.0 / x
        return math.pow(inv, lambda_r), math.pow(inv, lambda_a)

    if theta <= 1.0:
        total = 0.0
        for node, weight in zip(_GAUSS_NODES, _GAUSS_WEIGHTS, strict=True):
            if node < 1.0e-20:
                continue
            repulsive, attractive = powers(node)
            u = theta * (repulsive - attractive)
            total += weight * _exp(-u) * (repulsive - attractive)
        return total

    raw_cut = _exp(-math.log(theta + 20.0) / lambda_r)
    x_cut = min(max(raw_cut, 0.5), 0.999)
    # Where the clamp holds the cut is constant, and the derivative of a constant is zero.
    d_cut = -x_cut / (lambda_r * (theta + 20.0)) if 0.5 < raw_cut < 0.999 else 0.0
    width = 1.0 - x_cut
    d_width = -d_cut

    covered = 0.0
    moving = 0.0
    for node, weight in zip(_GAUSS_NODES, _GAUSS_WEIGHTS, strict=True):
        x = 1.0 - width * (1.0 - node)
        if x < 1.0e-20:
            continue
        repulsive, attractive = powers(x)
        u = theta * (repulsive - attractive)
        dx = -(1.0 - node) * d_width
        du = (repulsive - attractive) + theta * (
            -lambda_r * repulsive / x + lambda_a * attractive / x
        ) * dx
        covered += weight * (1.0 - _exp(-u))
        moving += weight * _exp(-u) * du
    return d_width * (covered - 1.0) + width * moving


def effective_diameter(component: _Component, t: float) -> float:
    """``d = sigma * integral_0^1 (1 - exp(-u*(x))) dx``.

    Raises:
        OutOfRangeError: if ``t`` is not positive, or the exponents do not make a potential.
    """
    return effective_diameter_and_log_derivative(component, t)[0]


def effective_diameter_and_log_derivative(component: _Component, t: float) -> tuple[float, float]:
    """The diameter, and ``T d ln d/dT`` beside it.

    ``d = sigma BH(theta)`` with ``theta = C eps/(k T)``, so ``T d ln d/dT`` is
    ``-theta d ln BH/dtheta`` - one number for every layer that needs the diameter's
    temperature dependence.

    Raises:
        OutOfRangeError: as :func:`effective_diameter`, or if the integral is not positive.
    """
    if not t > 0.0:
        raise OutOfRangeError("T", t, "the effective diameter is a function of temperature")
    if component.lambda_r <= component.lambda_a:
        raise InvalidInputError(
            "lambda_r",
            f"the repulsive exponent {component.lambda_r} is not above the attractive one "
            f"{component.lambda_a}, so there is no potential well to average",
        )
    theta = mie_prefactor(component.lambda_r, component.lambda_a) * component.epsik / t
    integral = barker_henderson(theta, component.lambda_r, component.lambda_a)
    if integral <= 0.0:
        raise OutOfRangeError(
            "T", t, f"the Barker-Henderson integral is {integral} at this temperature"
        )
    slope = barker_henderson_d_theta(theta, component.lambda_r, component.lambda_a)
    return component.sigma * integral, -theta * slope / integral


def eta_effective(eta: float, lambda_: float) -> float:
    """``eta_eff = c1 eta + c2 eta^2 + c3 eta^3 + c4 eta^4``."""
    inv = 1.0 / lambda_
    total = 0.0
    power = eta
    for coefficient in _ETA_EFF_COEFFS:
        c = (
            coefficient[0]
            + coefficient[1] * inv
            + coefficient[2] * inv * inv
            + coefficient[3] * inv * inv * inv
        )
        total += c * power
        power *= eta
    return total


def a_s1_bare(eta: float, lambda_: float) -> float:
    """The bare ``aS1`` a Mie exponent contributes."""
    effective = eta_effective(eta, lambda_)
    one = 1.0 - effective
    return -(1.0 - effective / 2.0) / one**3 / (lambda_ - 3.0)


def b_bare(eta: float, lambda_: float, x0: float) -> float:
    """The bare ``B`` a Mie exponent contributes, ``x0 = sigma/d``."""
    x0_3l = math.pow(x0, 3.0 - lambda_)
    cap_i = (1.0 - x0_3l) / (lambda_ - 3.0)
    cap_j = (1.0 - (lambda_ - 3.0) * math.pow(x0, 4.0 - lambda_) + (lambda_ - 4.0) * x0_3l) / (
        (lambda_ - 3.0) * (lambda_ - 4.0)
    )
    one = 1.0 - eta
    return cap_i * (1.0 - eta / 2.0) / one**3 - 9.0 * cap_j * eta * (eta + 1.0) / (2.0 * one**3)


def k_hs(eta: float) -> float:
    """The hard-sphere isothermal compressibility."""
    one = 1.0 - eta
    return one**4 / (1.0 + 4.0 * eta + 4.0 * eta * eta - 4.0 * eta**3 + eta**4)


def mie_alpha(lambda_r: float, lambda_a: float) -> float:
    """Lafitte's ``alpha``, the well's softness."""
    return mie_prefactor(lambda_r, lambda_a) * (1.0 / (lambda_a - 3.0) - 1.0 / (lambda_r - 3.0))


def contact_value_0(eta: float, x0: float) -> float:
    """The Mie contact value's quartic in ``x0``."""
    one = 1.0 - eta
    eta2, eta3, eta4 = eta * eta, eta**3, eta**4
    om3 = one**3
    k0 = -_log(one) + (42.0 * eta - 39.0 * eta2 + 9.0 * eta3 - 2.0 * eta4) / (6.0 * om3)
    k1 = (-12.0 * eta + 6.0 * eta2 + eta4) / (2.0 * om3)
    k2 = -3.0 * eta2 / (8.0 * one * one)
    k3 = (3.0 * eta + 3.0 * eta2 - eta4) / (6.0 * om3)
    return _exp(k0 + k1 * x0 + k2 * x0 * x0 + k3 * x0 * x0 * x0)


#: The relative step the central differences here are taken at, ``1e-5``, which is
#: NeqSim's own. Small enough that the truncation is invisible and the cancellation is not.
_RELATIVE_STEP = 1.0e-5


def _eta_step(eta: float) -> tuple[float, float, float]:
    """NeqSim's ``eta`` difference: its step, its floor, and the halved span."""
    step = max(abs(eta) * 1.0e-5, 1.0e-12)
    high = eta + step
    low = max(eta - step, 1.0e-15)
    return high, low, (high - low) / 2.0


def chain_g1(eta: float, lambda_r: float, lambda_a: float, c_mie: float, x0: float) -> float:
    """``g1``, the first-order chain perturbation."""

    def bare(e: float, lambda_: float) -> float:
        return a_s1_bare(e, lambda_) + b_bare(e, lambda_, x0)

    high, low, step = _eta_step(eta)

    def d_full(lambda_: float) -> float:
        return (high * bare(high, lambda_) - low * bare(low, lambda_)) / (2.0 * step)

    da1_drho = c_mie * (
        math.pow(x0, lambda_a) * d_full(lambda_a) - math.pow(x0, lambda_r) * d_full(lambda_r)
    )
    return 3.0 * da1_drho - c_mie * (
        lambda_a * math.pow(x0, lambda_a) * bare(eta, lambda_a)
        - lambda_r * math.pow(x0, lambda_r) * bare(eta, lambda_r)
    )


def chain_g2(
    eta: float,
    zeta_st: float,
    lambda_r: float,
    lambda_a: float,
    eps_over_kt: float,
    c_mie: float,
    x0: float,
) -> float:
    """``g2``, the second-order chain perturbation with Lafitte's ``gamma_c``."""

    def bare(e: float, lambda_: float) -> float:
        return a_s1_bare(e, lambda_) + b_bare(e, lambda_, x0)

    def inner(e: float) -> float:
        return (
            math.pow(x0, 2.0 * lambda_a) * bare(e, 2.0 * lambda_a)
            - 2.0 * math.pow(x0, lambda_a + lambda_r) * bare(e, lambda_a + lambda_r)
            + math.pow(x0, 2.0 * lambda_r) * bare(e, 2.0 * lambda_r)
        )

    high, low, step = _eta_step(eta)
    da2_product = (high * k_hs(high) * inner(high) - low * k_hs(low) * inner(low)) / (2.0 * step)
    da2_drho = 0.5 * c_mie * c_mie * da2_product
    g_mca2 = 3.0 * da2_drho - k_hs(eta) * c_mie * c_mie * (
        lambda_r * math.pow(x0, 2.0 * lambda_r) * bare(eta, 2.0 * lambda_r)
        - (lambda_a + lambda_r) * math.pow(x0, lambda_a + lambda_r) * bare(eta, lambda_a + lambda_r)
        + lambda_a * math.pow(x0, 2.0 * lambda_a) * bare(eta, 2.0 * lambda_a)
    )
    alpha = mie_alpha(lambda_r, lambda_a)
    theta = _exp(eps_over_kt) - 1.0
    gamma_c = (
        10.0
        * (-math.tanh(10.0 * (0.57 - alpha)) + 1.0)
        * zeta_st
        * theta
        * math.exp(-6.7 * zeta_st - 8.0 * zeta_st * zeta_st)
    )
    return (1.0 + gamma_c) * g_mca2


def chain_contact_value(
    components: Sequence[_Component],
    x: Sequence[float],
    t: float,
    eta: float,
    diameters: Sequence[float],
) -> float:
    """The mixture's chain contact value.

    **A fluid with no chain keeps the Carnahan-Starling value**: ``w_i = x_i (m_i - 1)`` is
    zero for a one-segment molecule, so there is no weight to take a mean over and NeqSim
    skips its block outright.
    """
    # The same two guards as the Rust kernel: outside ``(1e-10, 0.7)`` the contact value is
    # the Carnahan-Starling one, because the quartic the mean would need carries
    # ``log(1 - eta)``, which past one is not a number. A flash's trial compositions walk
    # through that region on the way to a liquid root.
    if not math.isfinite(eta) or eta < 1.0e-10 or eta > 0.7:
        return (1.0 - eta / 2.0) / (1.0 - eta) ** 3
    weight = sum(xi * (c.m - 1.0) for xi, c in zip(x, components, strict=True))
    if weight <= 1.0e-10:
        return (1.0 - eta / 2.0) / (1.0 - eta) ** 3
    weighted = 0.0
    for xi, component, d in zip(x, components, diameters, strict=True):
        w = xi * (component.m - 1.0)
        if w < 1.0e-30:
            continue
        x0 = component.sigma / d if d > 0.0 else 1.0
        eps_over_kt = component.epsik / t
        c_mie = mie_prefactor(component.lambda_r, component.lambda_a)
        zeta_st = eta * x0 * x0 * x0
        g0 = contact_value_0(eta, x0)
        g1 = chain_g1(eta, component.lambda_r, component.lambda_a, c_mie, x0)
        g2 = chain_g2(eta, zeta_st, component.lambda_r, component.lambda_a, eps_over_kt, c_mie, x0)
        weighted += w * _log(g0 * _exp(eps_over_kt * (g1 + eps_over_kt * g2) / g0))
    return math.exp(weighted / weight)


def pade_f(index: int, alpha: float) -> float:
    """Lafitte 2013's Pade approximant of the ``index``th function."""
    a2 = alpha * alpha
    numerator = (
        _PHI_PADE[0][index]
        + _PHI_PADE[1][index] * alpha
        + _PHI_PADE[2][index] * a2
        + _PHI_PADE[3][index] * a2 * alpha
    )
    denominator = (
        1.0
        + _PHI_PADE[4][index] * alpha
        + _PHI_PADE[5][index] * a2
        + _PHI_PADE[6][index] * a2 * alpha
    )
    return numerator / denominator


def a1_sutherland(eta: float, lambda_: float, eps_over_kt: float) -> float:
    """Sutherland's first-order attractive term."""
    effective = eta_effective(eta, lambda_)
    one = 1.0 - effective
    return -12.0 * eps_over_kt * eta / (lambda_ - 3.0) * (1.0 - effective / 2.0) / one**3


def b_correction(eta: float, lambda_: float, eps_over_kt: float, x0: float) -> float:
    """The correction the potential's softness adds to Sutherland's term."""
    x0_3l = math.pow(x0, 3.0 - lambda_)
    x0_4l = math.pow(x0, 4.0 - lambda_)
    cap_i = (1.0 - x0_3l) / (lambda_ - 3.0)
    cap_j = (1.0 - (lambda_ - 3.0) * x0_4l + (lambda_ - 4.0) * x0_3l) / (
        (lambda_ - 3.0) * (lambda_ - 4.0)
    )
    one = 1.0 - eta
    om3 = one**3
    b_bar = (1.0 - eta / 2.0) / om3 * cap_i - 9.0 * eta * (1.0 + eta) / (2.0 * om3) * cap_j
    return 12.0 * eta * eps_over_kt * b_bar


def a1_mie(
    eta: float, lambda_r: float, lambda_a: float, eps_over_kt: float, c_mie: float, x0: float
) -> float:
    """``a_1``, the mean-attractive term."""

    def bare(lambda_: float) -> float:
        return a1_sutherland(eta, lambda_, eps_over_kt) + b_correction(
            eta, lambda_, eps_over_kt, x0
        )

    return c_mie * (
        math.pow(x0, lambda_a) * bare(lambda_a) - math.pow(x0, lambda_r) * bare(lambda_r)
    )


def mie_chi(zeta_st: float, lambda_r: float, lambda_a: float) -> float:
    """``chi``, the second-order term's density correction."""
    alpha = mie_alpha(lambda_r, lambda_a)
    return (
        pade_f(0, alpha) * zeta_st + pade_f(1, alpha) * zeta_st**5 + pade_f(2, alpha) * zeta_st**8
    )


def a2_mie(
    eta: float,
    zeta_st: float,
    lambda_r: float,
    lambda_a: float,
    eps_over_kt: float,
    c_mie: float,
    x0: float,
) -> float:
    """``a_2``, the second-order perturbation term."""

    def bare(lambda_: float) -> float:
        return a1_sutherland(eta, lambda_, eps_over_kt) + b_correction(
            eta, lambda_, eps_over_kt, x0
        )

    inner = (
        math.pow(x0, 2.0 * lambda_a) * bare(2.0 * lambda_a)
        - 2.0 * math.pow(x0, lambda_a + lambda_r) * bare(lambda_a + lambda_r)
        + math.pow(x0, 2.0 * lambda_r) * bare(2.0 * lambda_r)
    )
    return (
        0.5
        * k_hs(eta)
        * (1.0 + mie_chi(zeta_st, lambda_r, lambda_a))
        * eps_over_kt
        * c_mie
        * c_mie
        * inner
    )


def a3_mie(zeta_st: float, lambda_r: float, lambda_a: float, eps_over_kt: float) -> float:
    """``a_3``, the third-order term."""
    alpha = mie_alpha(lambda_r, lambda_a)
    return (
        -(eps_over_kt**3)
        * pade_f(3, alpha)
        * zeta_st
        * _exp(pade_f(4, alpha) * zeta_st + pade_f(5, alpha) * zeta_st * zeta_st)
    )


def _b_bar_d_x0(eta: float, lambda_: float, x0: float) -> float:
    """``d b_bar/d x0``, the derivative of :func:`b_correction`'s bracket.

    ``d cap_i/d x0 = x0^(2-lambda)`` and ``d cap_j/d x0 = x0^(3-lambda) - x0^(2-lambda)``:
    the two exponents bring opposite signs down with them, so the second is a difference.
    """
    one = 1.0 - eta
    om3 = one**3
    contact = (1.0 - eta / 2.0) / om3
    cap_i = math.pow(x0, 2.0 - lambda_)
    cap_j = math.pow(x0, 3.0 - lambda_) - cap_i
    return contact * cap_i - 9.0 * eta * (1.0 + eta) / (2.0 * om3) * cap_j


def _fluid_dispersion(
    t: float, eta: float, component: _Component, d: float, t_d_ln_d: float
) -> tuple[float, float, float, float]:
    """One Mie fluid's three dispersion terms and ``T d/dT`` of their sum at a fixed ``eta``.

    The temperature reaches the dispersion through ``x0 = sigma/d(T)`` and
    ``beta = eps/(k T)`` and not through the packing fraction:

    ``T d a_k/dT = -lam (d a_k/d ln x0) - beta (d a_k/d beta)``

    with ``lam = T d ln d/dT``, the minus sign because ``x0`` and ``d`` move opposite ways.
    The second term is free - ``a_1``, ``a_2`` and ``a_3`` are homogeneous in ``beta`` of
    degree one, two and three - and the first is written out.
    """
    lambda_r, lambda_a = component.lambda_r, component.lambda_a
    beta = component.epsik / t
    x0 = component.sigma / d if d > 0.0 else 1.0
    c_mie = mie_prefactor(lambda_r, lambda_a)
    alpha = mie_alpha(lambda_r, lambda_a)
    zeta = eta * x0 * x0 * x0

    def bare(lambda_: float) -> float:
        return a1_sutherland(eta, lambda_, beta) + b_correction(eta, lambda_, beta, x0)

    def bare_d_x0(lambda_: float) -> float:
        return 12.0 * eta * beta * _b_bar_d_x0(eta, lambda_, x0)

    a1 = a1_mie(eta, lambda_r, lambda_a, beta, c_mie, x0)
    a2 = a2_mie(eta, zeta, lambda_r, lambda_a, beta, c_mie, x0)
    a3 = a3_mie(zeta, lambda_r, lambda_a, beta)

    d_a1_d_x0 = c_mie * (
        lambda_a * math.pow(x0, lambda_a - 1.0) * bare(lambda_a)
        + math.pow(x0, lambda_a) * bare_d_x0(lambda_a)
        - lambda_r * math.pow(x0, lambda_r - 1.0) * bare(lambda_r)
        - math.pow(x0, lambda_r) * bare_d_x0(lambda_r)
    )

    def inner(lambda_: float) -> float:
        return math.pow(x0, lambda_) * bare(lambda_)

    def d_inner_d_x0(lambda_: float) -> float:
        return lambda_ * math.pow(x0, lambda_ - 1.0) * bare(lambda_) + math.pow(
            x0, lambda_
        ) * bare_d_x0(lambda_)

    d_a2_d_x0 = (
        0.5
        * k_hs(eta)
        * beta
        * c_mie
        * c_mie
        * (
            (1.0 + mie_chi(zeta, lambda_r, lambda_a))
            * (
                d_inner_d_x0(2.0 * lambda_a)
                - 2.0 * d_inner_d_x0(lambda_a + lambda_r)
                + d_inner_d_x0(2.0 * lambda_r)
            )
            + (inner(2.0 * lambda_a) - 2.0 * inner(lambda_a + lambda_r) + inner(2.0 * lambda_r))
            * (
                pade_f(0, alpha)
                + 5.0 * pade_f(1, alpha) * zeta**4
                + 8.0 * pade_f(2, alpha) * zeta**7
            )
            * 3.0
            * eta
            * x0
            * x0
        )
    )

    exponent = pade_f(4, alpha) * zeta + pade_f(5, alpha) * zeta * zeta
    d_a3_d_zeta = (
        -(beta**3)
        * pade_f(3, alpha)
        * _exp(exponent)
        * (1.0 + zeta * (pade_f(4, alpha) + 2.0 * pade_f(5, alpha) * zeta))
    )

    # ``T d x0/dT = -x0 lam`` and ``T d zeta/dT = -3 zeta lam``.
    t_d_x0 = -x0 * t_d_ln_d
    return (
        a1,
        a2,
        a3,
        d_a1_d_x0 * t_d_x0
        - a1
        + d_a2_d_x0 * t_d_x0
        - 2.0 * a2
        - d_a3_d_zeta * 3.0 * zeta * t_d_ln_d
        - 3.0 * a3,
    )


def _segment_fractions(components: Sequence[_Component], x: Sequence[float]) -> list[float]:
    """``xi_i = x_i m_i / m_bar``, the segment fractions the dispersion sums over."""
    m_bar = sum(xi * c.m for xi, c in zip(x, components, strict=True))
    if m_bar <= 0.0:
        raise InvalidInputError(
            "components",
            "the mixture's segment number is zero, so no segment fraction is defined",
        )
    return [xi * c.m / m_bar for xi, c in zip(x, components, strict=True)]


def _pair_terms(
    first: _Component, second: _Component, t: float, eta: float
) -> tuple[float, float, float, float]:
    """One pair's three dispersion terms and their ``T d/dT`` at a packing fraction."""
    sigma_ij = 0.5 * (first.sigma + second.sigma)
    sigma3 = first.sigma**3 * second.sigma**3
    eps_ij = math.sqrt(first.epsik * second.epsik) * math.sqrt(sigma3) / sigma_ij**3
    lambda_r = 3.0 + math.sqrt((first.lambda_r - 3.0) * (second.lambda_r - 3.0))
    lambda_a = 3.0 + math.sqrt((first.lambda_a - 3.0) * (second.lambda_a - 3.0))
    cross = _Component(1.0, lambda_r, lambda_a, sigma_ij, eps_ij)
    d_ij, t_d_ln_d = effective_diameter_and_log_derivative(cross, t)
    return _fluid_dispersion(t, eta, cross, d_ij, t_d_ln_d)


def dispersion_at(
    components: Sequence[_Component],
    x: Sequence[float],
    t: float,
    eta: float,
    diameters: Sequence[float],
) -> tuple[float, float, float]:
    """The three dispersion terms, by NeqSim's branch.

    **One component evaluates directly and a mixture sums over pairs**, and the direct
    expression is not the pair sum at ``n = 1`` - it is a different expression whose cross
    parameters happen to be the pure ones.
    """
    if len(components) != 1:
        segment = _segment_fractions(components, x)
        a1 = a2 = a3 = 0.0
        for i, first in enumerate(components):
            for j, second in enumerate(components):
                weight = segment[i] * segment[j]
                if weight < 1.0e-30:
                    continue
                one, two, three, _ = _pair_terms(first, second, t, eta)
                a1 += weight * one
                a2 += weight * two
                a3 += weight * three
        return a1, a2, a3

    component = components[0]
    d = diameters[0] if diameters else 0.0
    _, t_d_ln_d = effective_diameter_and_log_derivative(component, t)
    one, two, three, _ = _fluid_dispersion(t, eta, component, d, t_d_ln_d)
    return one, two, three


def dispersion_t_d_dt(
    components: Sequence[_Component],
    x: Sequence[float],
    t: float,
    eta: float,
    diameters: Sequence[float],
) -> float:
    """``T d(a_1 + a_2 + a_3)/dT`` at a fixed packing fraction, on NeqSim's branch."""
    if len(components) != 1:
        segment = _segment_fractions(components, x)
        total = 0.0
        for i, first in enumerate(components):
            for j, second in enumerate(components):
                weight = segment[i] * segment[j]
                if weight < 1.0e-30:
                    continue
                total += weight * _pair_terms(first, second, t, eta)[3]
        return total

    component = components[0]
    d = diameters[0] if diameters else 0.0
    _, t_d_ln_d = effective_diameter_and_log_derivative(component, t)
    return _fluid_dispersion(t, eta, component, d, t_d_ln_d)[3]


def dispersion_row_sums(
    components: Sequence[_Component], x: Sequence[float], t: float, eta: float
) -> list[float]:
    """``sum_l xi_l a^{il}`` for each component, which the composition derivative needs."""
    segment = _segment_fractions(components, x)
    out: list[float] = []
    for first in components:
        total = 0.0
        for weight, second in zip(segment, components, strict=True):
            if weight < 1.0e-30:
                continue
            one, two, three, _ = _pair_terms(first, second, t, eta)
            total += weight * (one + two + three)
        out.append(total)
    return out


class _State(NamedTuple):
    """Every layer the Helmholtz energy is built from, at one state."""

    #: Each component's effective segment diameter, in metres.
    d: tuple[float, ...]
    #: ``sum_i x_i m_i``, the mixture's segment number.
    m_bar: float
    #: ``sum_i x_i (m_i - 1)``, the chain term's weight.
    m_minus_1: float
    #: ``sum_i x_i m_i d_i^3``, in m^3/mol.
    md3: float
    #: The packing fraction.
    eta: float
    #: The hard-sphere compressibility.
    a_hs: float
    #: The chain term's contact value.
    g_hs: float
    #: The three dispersion terms, per mole and unscaled.
    a1: float
    a2: float
    a3: float


def _state(components: Sequence[_Component], x: Sequence[float], t: float, v: float) -> _State:
    """Every layer the Helmholtz energy is built from, at one state."""
    for i, component in enumerate(components):
        if not component.has_parameters():
            raise InvalidInputError(
                "components",
                f"component {i} carries no SAFT-VR-Mie set (m = {component.m}, sigma = "
                f"{component.sigma}, epsilon/k = {component.epsik}). The table spells an "
                f"absent set as zeros rather than a blank, so this is a fluid with no "
                f"segments",
            )
    if not t > 0.0:
        raise OutOfRangeError("T", t, "the effective diameter is a function of temperature")
    if not v > 0.0:
        raise OutOfRangeError("v", v, "the packing fraction is a volume fraction")

    diameters = [effective_diameter(c, t) for c in components]
    m_bar = sum(xi * c.m for xi, c in zip(x, components, strict=True))
    m_minus_1 = sum(xi * (c.m - 1.0) for xi, c in zip(x, components, strict=True))
    md3 = sum(xi * c.m * d**3 for xi, c, d in zip(x, components, diameters, strict=True))
    eta = math.pi / 6.0 * AVOGADRO * md3 / v
    if eta >= 1.0:
        raise OutOfRangeError(
            "v",
            v,
            f"the packing fraction at this volume is {eta}, and the hard-sphere terms "
            f"diverge at one",
        )

    one = 1.0 - eta
    a_hs = (4.0 * eta - 3.0 * eta * eta) / (one * one)
    g_hs = chain_contact_value(components, x, t, eta, diameters)
    a1, a2, a3 = dispersion_at(components, x, t, eta, diameters)
    return _State(tuple(diameters), m_bar, m_minus_1, md3, eta, a_hs, g_hs, a1, a2, a3)


def chain_contact_eta(
    components: Sequence[_Component],
    x: Sequence[float],
    t: float,
    eta: float,
    diameters: Sequence[float],
) -> float:
    """``d g_hs/d eta`` at a fixed temperature, by NeqSim's central difference."""
    high, low, step = _eta_step(eta)
    return (
        chain_contact_value(components, x, t, high, diameters)
        - chain_contact_value(components, x, t, low, diameters)
    ) / (2.0 * step)


def dispersion_eta(
    components: Sequence[_Component],
    x: Sequence[float],
    t: float,
    eta: float,
    diameters: Sequence[float],
) -> float:
    """``d(a_1 + a_2 + a_3)/d eta`` at a fixed temperature, by the same difference."""

    def total(e: float) -> float:
        a1, a2, a3 = dispersion_at(components, x, t, e, diameters)
        return a1 + a2 + a3

    high, low, step = _eta_step(eta)
    return (total(high) - total(low)) / (2.0 * step)


def t_d_helmholtz_rt_dt(
    components: Sequence[_Component], x: Sequence[float], t: float, state: _State
) -> float:
    """``T d(A^R/(R T))/dT`` at constant volume, per mole.

    **NeqSim publishes this and its value is wrong for any fluid with a chain.**
    ``PhaseSAFTVRMie.dF_HC_SAFTdT`` differentiates ``g_hs`` as though it moved with the
    packing fraction alone; it is the Mie-weighted contact value, so it carries
    ``beta = eps/(k T)`` as well, and the class is missing that route and the diameter's.
    For methane (``m = 1``) its chain block is skipped and the two agree; for n-butane they
    disagree in sign. Its ``dF_DISP_SAFTdT`` is right, and its volume derivatives are right.

    The check is a difference of NeqSim's own ``F`` at a pinned molar volume, so the ``F``
    below is with the diameters held and the third route is the diameters' own.
    """
    t_d_ln_d = []
    t_d_md3 = 0.0
    for i, component in enumerate(components):
        _, lam = effective_diameter_and_log_derivative(component, t)
        t_d_ln_d.append(lam)
        t_d_md3 += x[i] * component.m * 3.0 * state.d[i] ** 3 * lam
    t_d_eta = state.eta * t_d_md3 / state.md3

    eta = state.eta
    one = 1.0 - eta
    a_hs_d_eta = ((4.0 - 6.0 * eta) * one + 2.0 * (4.0 * eta - 3.0 * eta * eta)) / one**3

    g_eta = chain_contact_eta(components, x, t, eta, state.d)
    t_step = t * _RELATIVE_STEP
    g_t = (
        chain_contact_value(components, x, t + t_step, eta, state.d)
        - chain_contact_value(components, x, t - t_step, eta, state.d)
    ) / (2.0 * t_step)
    # Each diameter scaled by its own ``lam_i``, so the difference is ``sum_i d_i lam_i
    # d g/d d_i`` - which is exactly what the diameter's route contributes.
    high = [d * (1.0 + lam * _RELATIVE_STEP) for d, lam in zip(state.d, t_d_ln_d, strict=True)]
    low = [d * (1.0 - lam * _RELATIVE_STEP) for d, lam in zip(state.d, t_d_ln_d, strict=True)]
    g_d = (
        chain_contact_value(components, x, t, eta, high)
        - chain_contact_value(components, x, t, eta, low)
    ) / (2.0 * _RELATIVE_STEP)
    t_d_g_hs = g_eta * t_d_eta + t * g_t + g_d

    t_d_f_hc = state.m_bar * a_hs_d_eta * t_d_eta - state.m_minus_1 * t_d_g_hs / state.g_hs

    dispersion_eta_value = dispersion_eta(components, x, t, eta, state.d)
    t_d_dispersion = dispersion_t_d_dt(components, x, t, eta, state.d)
    return t_d_f_hc + state.m_bar * (t_d_dispersion + dispersion_eta_value * t_d_eta)


def _departure(
    components: Sequence[_Component], x: Sequence[float], t: float, v: float, z: float
) -> tuple[float, float]:
    """``(H^R/(RT), S^R/R)``, from ``Z - 1 - T dF/dT`` and the ``P``-to-``V`` conversion."""
    state = _state(components, x, t, v)
    t_d_f = t_d_helmholtz_rt_dt(components, x, t, state)
    return z - 1.0 - t_d_f, _log(z) - t_d_f - _energy(state)


def _energy(state: _State) -> float:
    """``A^R/(RT)`` per mole.

    **The dispersion carries the segment number as well as the mole number**, so
    ``a_1 + a_2 + a_3`` is exact only for a fluid whose ``m_bar`` is one.
    """
    dispersion = state.a1 + state.a2 + state.a3
    return (
        state.m_bar * state.a_hs - state.m_minus_1 * math.log(state.g_hs) + state.m_bar * dispersion
    )


def _pressure_over_rt(
    components: Sequence[_Component], x: Sequence[float], t: float, v: float
) -> float:
    """``P/(RT) = 1/v + eta f_eta/v``, with NeqSim's ``eta`` differences."""
    state = _state(components, x, t, v)
    eta = state.eta
    high, low, step = _eta_step(eta)

    g_eta = (
        chain_contact_value(components, x, t, high, state.d)
        - chain_contact_value(components, x, t, low, state.d)
    ) / (2.0 * step)

    def dispersion(e: float) -> float:
        a1, a2, a3 = dispersion_at(components, x, t, e, state.d)
        return a1 + a2 + a3

    dispersion_eta = (dispersion(high) - dispersion(low)) / (2.0 * step)
    a_hs_eta = (4.0 - 2.0 * eta) / (1.0 - eta) ** 3
    f_eta = (
        state.m_bar * a_hs_eta - state.m_minus_1 * g_eta / state.g_hs + state.m_bar * dispersion_eta
    )
    return 1.0 / v + eta * f_eta / v


#: The step the pressure's slope is differenced at, as a fraction of the volume. Ten times
#: NeqSim's own ``eta`` step, because this model's pressure is already a central difference
#: and differencing that again at ``1e-5`` puts the slope in cancellation.
_SLOPE_STEP = 1.0e-4


def _molar_volume(
    components: Sequence[_Component],
    x: Sequence[float],
    t: float,
    p: float,
    side: str,
) -> tuple[float, float, int]:
    """Solve ``P_calc(v) = p`` for ``(v, Z, iterations)``, as ``volume_solve.rs`` does."""
    if not p > 0.0:
        raise OutOfRangeError("P", p, "the volume solve divides by the pressure")

    rt = R * t
    ideal = rt / p
    md3 = sum(xi * c.m * effective_diameter(c, t) ** 3 for xi, c in zip(x, components, strict=True))
    floor = math.pi / 6.0 * AVOGADRO * md3

    def residual(v: float) -> tuple[float, float]:
        pressure = _pressure_over_rt(components, x, t, v) * rt
        step = _SLOPE_STEP * v
        high, low = v + step, v - step
        if low <= floor:
            # One-sided near the floor: a difference across a divergence is not a slope.
            upper = _pressure_over_rt(components, x, t, high) * rt
            slope = (upper - pressure) / step
        else:
            upper = _pressure_over_rt(components, x, t, high) * rt
            lower = _pressure_over_rt(components, x, t, low) * rt
            slope = (upper - lower) / (2.0 * step)
        return pressure - p, slope

    def value(v: float) -> float:
        return residual(v)[0]

    def bisect(lo: float, hi: float) -> float:
        sign = value(lo) < 0.0
        for _ in range(200):
            mid = 0.5 * (lo + hi)
            if (value(mid) < 0.0) == sign:
                lo = mid
            else:
                hi = mid
        return 0.5 * (lo + hi)

    def brackets() -> list[tuple[float, float]]:
        start = floor * (1.0 + 1.0e-9)
        top = 1.0e4 * ideal
        out: list[tuple[float, float]] = []
        previous = start
        previous_value = value(previous)
        for step in range(1, 129):
            v = start * math.pow(top / start, step / 128.0)
            current = value(v)
            if (previous_value < 0.0) != (current < 0.0):
                out.append((previous, v))
            previous, previous_value = v, current
        return out

    if side == "vapour":
        v = ideal
        converged = None
        for step in range(100):
            delta = -0.9 * value(v) / residual(v)[1]
            if not math.isfinite(delta) or v + delta <= floor:
                break
            relative = abs(delta) / v
            v += delta
            if relative < 1.0e-10:
                converged = step + 1
                break
        if converged is not None:
            return v, p * v / rt, converged
        found = brackets()
        if not found:
            raise OutOfRangeError("P", p, f"the vapour branch has no zero at {t} K and {p} Pa")
        lo, hi = found[-1]
    else:
        found = brackets()
        if not found:
            raise OutOfRangeError("P", p, f"the liquid branch has no zero at {t} K and {p} Pa")
        lo, hi = found[0]
    v = bisect(lo, hi)
    return v, p * v / rt, 0


def _ln_phi(
    components: Sequence[_Component], x: Sequence[float], t: float, v: float
) -> list[float]:
    """``ln phi_i = d(nF)/dn_i - ln Z``, on NeqSim's branch."""
    state = _state(components, x, t, v)
    z = _pressure_over_rt(components, x, t, v) * v
    if not z > 0.0:
        raise OutOfRangeError("Z", z, "the fugacity coefficients carry `- ln Z`")

    if len(components) == 1:
        return [_energy(state) + z - 1.0 - math.log(z)]

    eta = state.eta
    high, low, step = _eta_step(eta)
    a_hs_eta = (4.0 - 2.0 * eta) / (1.0 - eta) ** 3
    g_eta = (
        chain_contact_value(components, x, t, high, state.d)
        - chain_contact_value(components, x, t, low, state.d)
    ) / (2.0 * step)

    def dispersion(e: float) -> float:
        a1, a2, a3 = dispersion_at(components, x, t, e, state.d)
        return a1 + a2 + a3

    dispersion_eta = (dispersion(high) - dispersion(low)) / (2.0 * step)
    rows = dispersion_row_sums(components, x, t, eta)
    a_disp = state.a1 + state.a2 + state.a3

    out: list[float] = []
    for i, component in enumerate(components):
        eta_i = math.pi / 6.0 * AVOGADRO * component.m * state.d[i] ** 3 / v
        hard_chain = (
            component.m * state.a_hs
            + state.m_bar * a_hs_eta * eta_i
            - (component.m - 1.0) * math.log(state.g_hs)
            - state.m_minus_1 * g_eta / state.g_hs * eta_i
        )
        dispersion_i = component.m * (2.0 * rows[i] - a_disp) + state.m_bar * dispersion_eta * eta_i
        out.append(hard_chain + dispersion_i - math.log(z))
    return out


def saft_vr_mie_phase(
    components: list[str],
    T: Q,
    P: Q,
    z: list[float],
    compressed_phase: str,
) -> SaftVrMiePhaseResult:
    """One SAFT-VR-Mie phase's state at a temperature, pressure and composition.

    Args:
        components: the substance names, resolved against the databank with their
            ``m``/``sigma``/``epsilon/k`` and their two exponents.
        T: absolute temperature.
        P: absolute pressure.
        z: the mole fractions, checked rather than renormalised.
        compressed_phase: ``"liquid"`` or ``"vapour"``.

    Raises:
        InvalidInputError: if a component has no SAFT-VR-Mie set, if ``z`` is not a
            composition of the right length, or if ``compressed_phase`` is neither name.
        OutOfRangeError: if ``T`` or ``P`` is not positive, or the wanted branch has no root.

    See :func:`azoth.eos.saft_vr_mie_phase`.
    """
    from azoth import _models_gen

    spec = _models_gen.model(MODEL_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    t_si = input_to_si(spec, "T", T)
    p_si = input_to_si(spec, "P", P)
    apply_checks(checks.on_input, {"T": t_si, "P": p_si}.get, warnings)

    if compressed_phase.strip().lower() not in ("liquid", "vapour", "vapor"):
        raise InvalidInputError(
            "compressed_phase",
            f'{compressed_phase!r} is not a side. The spec\'s values are "liquid" and "vapour"',
        )
    side = "liquid" if compressed_phase.strip().lower() == "liquid" else "vapour"

    if len(components) != len(z):
        raise InvalidInputError(
            "z", f"{len(components)} components but {len(z)} fractions; the two must match"
        )
    total = sum(z)
    if abs(total - 1.0) > 1.0e-9:
        raise InvalidInputError(
            "z",
            f"the mole fractions sum to {total}, not to one. Renormalising them here would "
            f"make a composition error invisible in every number downstream, so it is "
            f"refused instead",
        )

    resolved = [_Component.of(name) for name in components]
    v, z_factor, _ = _molar_volume(resolved, z, t_si, p_si, side)
    h_over_rt, s_over_r = _departure(resolved, z, t_si, v, z_factor)
    return SaftVrMiePhaseResult(
        z_factor=z_factor,
        ln_phi=tuple(_ln_phi(resolved, z, t_si, v)),
        v=ureg.Quantity(v, "m**3/mol"),
        h_res=ureg.Quantity(h_over_rt * R * t_si, "J/mol"),
        s_res=ureg.Quantity(s_over_r * R, "J/(mol*K)"),
        warnings=tuple(warnings),
    )


__all__ = ["saft_vr_mie_phase"]
