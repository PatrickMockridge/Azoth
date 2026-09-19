"""``eos.pcsaft_rahmat_phase`` - the PC-SAFT phase state, the pure-Python reference.

The second, independent expression of the physics the Rust ``azoth_eos::pcsaft_rahmat_phase``
computes: Gross and Sadowski's perturbed-chain statistical associating fluid theory
without its association term, at a temperature, a pressure and a composition.

Spec: ``specs/models/eos/pcsaft_rahmat_phase.toml``, which carries why the two constants are
NeqSim's rounded literals rather than CODATA's, why a component with no parameter set is
refused, and why the branch is a caller's statement.

**The numbers this is checked against are NeqSim's.** ``SystemPCSAFT`` builds
``PhasePCSAFTRahmat``, and ``validation/neqsim/PcsaftProbe.java`` prints that class's
intermediates one per line: the segment diameters, ``md3``, ``m_bar``, the packing
fraction, both dispersion sums, ``I1``, ``I2``, ``C1``, the three Helmholtz terms, the
compressibility and every fugacity coefficient. Both kernels reproduce them, so a
disagreement between them is a disagreement about a layer rather than about a total.
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from typing import NamedTuple

from azoth.core.errors import InvalidInputError, OutOfRangeError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import PcsaftRahmatPhaseResult
from azoth.core.units import Q, input_to_si, ureg
from azoth.core.warnings import Warning
from azoth.eos.components import entry, pcsaft_kij_for

MODEL_ID = "eos.pcsaft_rahmat_phase"

#: Avogadro's constant, **NeqSim's value**, in 1/mol. ``6.023e23`` is the pre-2019
#: definition, a relative ``1.4e-4`` from CODATA's; the packing fraction and both
#: dispersion terms are built on it, so the modern one would diverge from every number
#: the oracle prints. ``6.02214076e23`` is the correct constant and the wrong port.
AVOGADRO = 6.023e23

#: NeqSim's ``ThermodynamicConstantsInterface.R``, for the same reason: the volume the
#: oracle converges to is set by ``P/(RT)``.
R = 8.3144621

#: Gross-Sadowski's ``I1`` constants, ``[order in m][power of eta]``.
A_CONST = (
    (
        0.9105631445,
        0.6361281449,
        2.6861347891,
        -26.547362491,
        97.759208784,
        -159.59154087,
        91.297774084,
    ),
    (
        -0.3084016918,
        0.1860531159,
        -2.5030047259,
        21.419793629,
        -65.255885330,
        83.318680481,
        -33.746922930,
    ),
    (
        -0.0906148351,
        0.4527842806,
        0.5962700728,
        -1.7241829131,
        -4.1302112531,
        13.776631870,
        -8.6728470368,
    ),
)

#: Gross-Sadowski's ``I2`` constants, same layout.
B_CONST = (
    (
        0.7240946941,
        2.2382791861,
        -4.0025849485,
        -21.003576815,
        26.855641363,
        206.55133841,
        -355.60235612,
    ),
    (
        -0.5755498075,
        0.6995095521,
        3.8925673390,
        -17.215471648,
        192.67226447,
        -161.826465,
        -165.20769346,
    ),
    (
        0.0976883116,
        -0.2557574982,
        -9.1558561530,
        20.642075974,
        -38.804430052,
        93.626774077,
        -29.666905585,
    ),
)


class _Component(NamedTuple):
    """One component's PC-SAFT set, resolved from the databank."""

    #: The number of segments.
    m: float
    #: The temperature-independent segment diameter, in metres.
    sigma: float
    #: The segment energy over Boltzmann's constant, in K.
    epsik: float

    @classmethod
    def of(cls, name: str) -> _Component:
        record = entry(name)
        return cls(record.m_saft, record.sigma_saft, record.epsik_saft)

    def has_parameters(self) -> bool:
        """Whether the table gave this component a set: zero in all three means absent."""
        return self.m > 0.0 and self.sigma > 0.0 and self.epsik > 0.0


def _segment_diameter(sigma: float, epsik: float, t: float) -> float:
    """``d_i = sigma_i (1 - 0.12 exp(-3 (epsilon_i/k)/T))``, in metres."""
    return sigma * (1.0 - 0.12 * math.exp(-3.0 * epsik / t))


def _segment_diameter_d_t(sigma: float, epsik: float, t: float) -> float:
    """``d d_i/dT``, the derivative of :func:`_segment_diameter`."""
    return -0.36 * sigma * epsik * math.exp(-3.0 * epsik / t) / (t * t)


def _series(
    order: Sequence[Sequence[float]], m_bar: float, eta: float
) -> tuple[float, float, float, float]:
    """One series and its derivatives: value, in ``eta`` twice, and in ``m_bar``.

    ``a_i`` is linear in ``(m_bar-1)/m_bar`` and ``(m_bar-1)(m_bar-2)/m_bar^2``, so the
    ``m_bar`` derivative is the two ratios' own derivative against the coefficients they
    carry - and ``a_0`` is a constant, whose derivative is zero.
    """
    one = (m_bar - 1.0) / m_bar
    two = one * (m_bar - 2.0) / m_bar
    value = d_eta = d2_eta = d_m_bar = 0.0
    for i, (a0, a1, a2) in enumerate(zip(*order, strict=True)):
        a_i = a0 + one * a1 + two * a2
        power = eta**i
        value += a_i * power
        if i:
            d_eta += i * a_i * eta ** (i - 1)
        if i > 1:
            d2_eta += i * (i - 1) * a_i * eta ** (i - 2)
        d_m_bar += (a1 / m_bar**2 + a2 * (3.0 * m_bar - 4.0) / m_bar**3) * power
    return value, d_eta, d2_eta, d_m_bar


def _c1_terms(eta: float) -> tuple[float, float, float, float, float, float]:
    """``A``, ``B`` and their first two derivatives in ``eta``.

    ``C1 = 1/(1 + m_bar A + (1 - m_bar) B)``. Each has the shape ``(num)/(1-eta)^k``, or
    ``num/P^2`` with ``P = (1-eta)(2-eta)``, so each derivative is one quotient rule.
    """
    one = 1.0 - eta
    u_a = 8.0 + 20.0 * eta - 4.0 * eta * eta
    a = (8.0 * eta - 2.0 * eta * eta) / one**4
    a_d_eta = u_a / one**5
    a_d2_eta = ((20.0 - 8.0 * eta) * one + 5.0 * u_a) / one**6

    p = (1.0 - eta) * (2.0 - eta)
    p_d_eta = 2.0 * eta - 3.0
    num = 20.0 * eta - 27.0 * eta * eta + 12.0 * eta**3 - 2.0 * eta**4
    num_d_eta = 20.0 - 54.0 * eta + 36.0 * eta * eta - 8.0 * eta**3
    num_d2_eta = -54.0 + 72.0 * eta - 24.0 * eta * eta
    f = num_d_eta * p - 2.0 * num * p_d_eta
    f_d_eta = num_d2_eta * p - num_d_eta * p_d_eta - 4.0 * num

    return a, a_d_eta, a_d2_eta, num / p**2, f / p**3, (f_d_eta * p - 3.0 * f * p_d_eta) / p**4


class _State:
    """Every layer PC-SAFT's Helmholtz energy is built from, at one state."""

    __slots__ = (
        "a_hs",
        "c1",
        "d",
        "epsik",
        "eta",
        "g_hs",
        "i1",
        "i2",
        "kij",
        "m",
        "m_bar",
        "m_minus_1",
        "md3",
        "rho",
        "s1",
        "s2",
        "sigma",
        "t",
        "v",
        "x",
    )

    def __init__(
        self,
        components: Sequence[_Component],
        kij: Sequence[float],
        x: Sequence[float],
        t: float,
        v: float,
    ) -> None:
        n = len(components)
        if not t > 0.0:
            raise OutOfRangeError("T", t, "the segment diameter is a function of temperature")
        if not v > 0.0:
            raise OutOfRangeError("v", v, "the packing fraction is a volume fraction")
        for i, c in enumerate(components):
            if not c.has_parameters():
                raise InvalidInputError(
                    "components",
                    f"component {i} carries no PC-SAFT set (m = {c.m}, sigma = {c.sigma}, "
                    f"epsilon/k = {c.epsik}). The table spells an absent set as zeros rather "
                    f"than a blank, so this is a fluid with no segments",
                )

        self.m = [c.m for c in components]
        self.sigma = [c.sigma for c in components]
        self.epsik = [c.epsik for c in components]
        self.d = [_segment_diameter(c.sigma, c.epsik, t) for c in components]
        self.kij, self.x, self.t, self.v = kij, x, t, v
        self.m_bar = sum(xi * c.m for xi, c in zip(x, components, strict=True))
        self.m_minus_1 = sum(xi * (c.m - 1.0) for xi, c in zip(x, components, strict=True))
        self.md3 = sum(xi * c.m * di**3 for xi, c, di in zip(x, components, self.d, strict=True))
        self.eta = math.pi / 6.0 * AVOGADRO * self.md3 / v
        if self.eta >= 1.0:
            raise OutOfRangeError(
                "v",
                v,
                f"the packing fraction at this volume is {self.eta}, and the hard-sphere "
                f"terms diverge at one",
            )

        eta = self.eta
        self.a_hs = (4.0 * eta - 3.0 * eta * eta) / (1.0 - eta) ** 2
        self.g_hs = (1.0 - eta / 2.0) / (1.0 - eta) ** 3

        s1 = s2 = 0.0
        for i in range(n):
            for j in range(n):
                sigma_ij = 0.5 * (components[i].sigma + components[j].sigma)
                e_ij = math.sqrt(components[i].epsik / t * (components[j].epsik / t))
                weight = x[i] * x[j] * components[i].m * components[j].m * sigma_ij**3
                one_minus_k = 1.0 - kij[i * n + j]
                s1 += weight * e_ij * one_minus_k
                s2 += weight * e_ij * e_ij * one_minus_k * one_minus_k
        self.s1, self.s2 = s1, s2

        a, _, _, b, _, _ = _c1_terms(eta)
        self.c1 = 1.0 / (1.0 + self.m_bar * a + (1.0 - self.m_bar) * b)
        self.i1 = _series(A_CONST, self.m_bar, eta)[0]
        self.i2 = _series(B_CONST, self.m_bar, eta)[0]
        self.rho = AVOGADRO / v

    def t_d_helmholtz_rt_dt(self) -> float:
        """``T d(A^R/(RT))/dT`` at constant volume, per mole.

        **NeqSim publishes this and its PC-SAFT value cannot be used.**
        ``PhasePCSAFT.getdDSAFTdT`` multiplies the chain rule for the segment diameter by
        an extra ``3 d_i**2``, so every ``eta``-dependent part of its ``dFdT`` is about
        ``1e-18`` too small. The write-up is at
        ``~/Desktop/neqsim-pcsaft-hard-chain-temperature-derivative.md``; the check this
        is held to is a finite difference of NeqSim's own ``F`` at fixed volume.

        ``eta`` is the only place the volume and the temperature meet, so at constant
        ``v`` everything moves through ``T d eta/dT``; the two dispersion sums each carry
        a ``1/T`` of their own, which is the ``-S`` beside each.
        """
        t_d_md3 = sum(
            self.x[i]
            * self.m[i]
            * 3.0
            * self.d[i] ** 2
            * (self.t * _segment_diameter_d_t(self.sigma[i], self.epsik[i], self.t))
            for i in range(len(self.x))
        )
        t_d_eta = self.eta * t_d_md3 / self.md3

        a_hs_d_eta, _, g_hs_d_eta, _, c1_d_eta, _ = self._eta_layers()
        t_d_f_hc = (
            self.m_bar * a_hs_d_eta * t_d_eta - self.m_minus_1 * g_hs_d_eta / self.g_hs * t_d_eta
        )
        t_d_c1 = c1_d_eta * t_d_eta
        i1_d_eta = _series(A_CONST, self.m_bar, self.eta)[1]
        i2_d_eta = _series(B_CONST, self.m_bar, self.eta)[1]
        t_d_f_disp1 = -2.0 * math.pi * self.rho * self.s1 * (i1_d_eta * t_d_eta - self.i1)
        t_d_f_disp2 = (
            -math.pi
            * self.m_bar
            * self.rho
            * self.s2
            * (self.i2 * (-self.c1 + t_d_c1) + self.c1 * (i2_d_eta * t_d_eta - self.i2))
        )
        return t_d_f_hc + t_d_f_disp1 + t_d_f_disp2

    def departure(self) -> tuple[float, float]:
        """``(h_res/(R T), s_res/R)`` at this state.

        ``Hres/(R T) = Z - 1 - T dF/dT``: the two ``F``-terms in NeqSim's
        ``AresTV + T SresTV + P V - n R T`` cancel, and the entropy then carries the
        ``P``-to-``V`` conversion, ``SresTP = SresTV + n R ln Z``.
        """
        z = self.pressure_over_rt() * self.v
        t_d_f = self.t_d_helmholtz_rt_dt()
        return z - 1.0 - t_d_f, math.log(z) - t_d_f - self.f()

    def f(self) -> float:
        """``A^R/(RT)`` per mole: the hard-sphere, chain and two dispersion terms."""
        return (
            self.m_bar * self.a_hs
            - self.m_minus_1 * math.log(self.g_hs)
            - 2.0 * math.pi * self.rho * self.s1 * self.i1
            - math.pi * self.m_bar * self.rho * self.s2 * self.i2 * self.c1
        )

    def pressure_over_rt(self) -> float:
        """``1/v - eta F_eta/v``, with ``F_eta`` taken along the volume path."""
        f_eta, _ = self._f_eta_path()
        return 1.0 / self.v + self.eta * f_eta / self.v

    def d_pressure_over_rt_dv(self) -> float:
        """``-(1 + 2 eta F_eta + eta^2 F_etaeta)/v^2``, the volume solve's denominator."""
        f_eta, f_eta2 = self._f_eta_path()
        eta = self.eta
        return -(1.0 + 2.0 * eta * f_eta + eta * eta * f_eta2) / (self.v * self.v)

    def _eta_layers(self) -> tuple[float, float, float, float, float, float]:
        """``a_hs'``, ``a_hs''``, ``g_hs'``, ``g_hs''``, ``C1'``, ``C1''`` at this state.

        The volume enters the energy only through the packing fraction, and the number
        density is proportional to it, so the whole volume path is travelled in ``eta``
        with ``rho/eta`` constant - which is why one derivative serves both callers.
        """
        eta, one = self.eta, 1.0 - self.eta
        a_hs_d_eta = (4.0 - 2.0 * eta) / one**3
        a_hs_d2_eta = (10.0 - 4.0 * eta) / one**4
        g_hs_d_eta = (2.5 - eta) / one**4
        g_hs_d2_eta = (9.0 - 3.0 * eta) / one**5

        _, a_d_eta, a_d2_eta, _, b_d_eta, b_d2_eta = _c1_terms(eta)
        w_d_eta = self.m_bar * a_d_eta + (1.0 - self.m_bar) * b_d_eta
        w_d2_eta = self.m_bar * a_d2_eta + (1.0 - self.m_bar) * b_d2_eta
        c1_d_eta = -self.c1 * self.c1 * w_d_eta
        c1_d2_eta = 2.0 * self.c1**3 * w_d_eta * w_d_eta - self.c1 * self.c1 * w_d2_eta
        return a_hs_d_eta, a_hs_d2_eta, g_hs_d_eta, g_hs_d2_eta, c1_d_eta, c1_d2_eta

    def _f_eta_path(self) -> tuple[float, float]:
        """``F_eta`` and ``F_etaeta`` along the volume path, where ``rho`` moves with eta."""
        eta = self.eta
        a_hs_d_eta, a_hs_d2_eta, g_hs_d_eta, g_hs_d2_eta, c1_d_eta, c1_d2_eta = self._eta_layers()
        _, i1_d_eta, i1_d2_eta, _ = _series(A_CONST, self.m_bar, eta)
        _, i2_d_eta, i2_d2_eta, _ = _series(B_CONST, self.m_bar, eta)

        # Both dispersion brackets are `(X + eta X')` with `X = I1` and `X = I2 C1`, and
        # each differentiates to `(2 X' + eta X'')`.
        d_d_eta = i2_d_eta * self.c1 + self.i2 * c1_d_eta
        d_d2_eta = i2_d2_eta * self.c1 + 2.0 * i2_d_eta * c1_d_eta + self.i2 * c1_d2_eta
        rho_over_eta = self.rho / eta

        f_eta = (
            self.m_bar * a_hs_d_eta
            - self.m_minus_1 * g_hs_d_eta / self.g_hs
            - 2.0 * math.pi * rho_over_eta * self.s1 * (self.i1 + eta * i1_d_eta)
            - math.pi * self.m_bar * rho_over_eta * self.s2 * (self.i2 * self.c1 + eta * d_d_eta)
        )
        f_d2_eta = (
            self.m_bar * a_hs_d2_eta
            - self.m_minus_1 * (g_hs_d2_eta / self.g_hs - (g_hs_d_eta / self.g_hs) ** 2)
            - 2.0 * math.pi * rho_over_eta * self.s1 * (2.0 * i1_d_eta + eta * i1_d2_eta)
            - math.pi * self.m_bar * rho_over_eta * self.s2 * (2.0 * d_d_eta + eta * d_d2_eta)
        )
        return f_eta, f_d2_eta

    def ln_fugacity_coefficients(self) -> list[float]:
        """``ln phi_i = d(nF)/dn_i - ln Z``, at fixed temperature and volume.

        A different derivative from the volume one. ``rho`` and ``eta`` both scale with
        ``1/v``, so the volume path can be travelled in ``eta`` alone - but at fixed
        volume the composition moves them apart: ``rho`` through the mole numbers and
        ``eta`` through ``sum n_i m_i d_i^3``. The packing fraction's own derivative is
        ``(pi/6) N_A m_k d_k^3/v``, **not** ``eta/v``.
        """
        eta = self.eta
        a_hs_d_eta = (4.0 - 2.0 * eta) / (1.0 - eta) ** 3
        g_hs_d_eta = (2.5 - eta) / (1.0 - eta) ** 4
        _, a_d_eta, _, _, b_d_eta, _ = _c1_terms(eta)
        c1_d_eta = -self.c1 * self.c1 * (self.m_bar * a_d_eta + (1.0 - self.m_bar) * b_d_eta)
        a_series = _series(A_CONST, self.m_bar, eta)
        b_series = _series(B_CONST, self.m_bar, eta)
        a, _, _, b, _, _ = _c1_terms(eta)
        c1_d_m_bar = -self.c1 * self.c1 * (a - b)

        # `f`'s partials in each thing the composition moves.
        f_rho = (
            -2.0 * math.pi * self.s1 * self.i1 - math.pi * self.m_bar * self.s2 * self.i2 * self.c1
        )
        f_eta = (
            self.m_bar * a_hs_d_eta
            - self.m_minus_1 * g_hs_d_eta / self.g_hs
            - 2.0 * math.pi * self.rho * self.s1 * a_series[1]
            - math.pi
            * self.m_bar
            * self.rho
            * self.s2
            * (b_series[1] * self.c1 + self.i2 * c1_d_eta)
        )
        f_m_bar = (
            self.a_hs
            - 2.0 * math.pi * self.rho * self.s1 * a_series[3]
            - math.pi
            * self.rho
            * self.s2
            * (self.i2 * self.c1 + self.m_bar * (b_series[3] * self.c1 + self.i2 * c1_d_m_bar))
        )
        f_m_minus_1 = -math.log(self.g_hs)
        f_s1 = -2.0 * math.pi * self.rho * self.i1
        f_s2 = -math.pi * self.m_bar * self.rho * self.i2 * self.c1

        z = self.pressure_over_rt() * self.v
        if not z > 0.0:
            raise OutOfRangeError("Z", z, "the fugacity coefficients are `d(nF)/dn_i - ln Z`")

        n = len(self.x)
        out: list[float] = []
        for k in range(n):
            # `S1` and `S2` are quadratic in the mole fractions, so their derivatives are
            # row sums: `-2 S1 + 2 sum_j x_j A_kj`.
            s1_k = -2.0 * self.s1
            s2_k = -2.0 * self.s2
            for j in range(n):
                sigma_ij = 0.5 * (self.sigma[k] + self.sigma[j])
                e_ij = math.sqrt(self.epsik[k] / self.t * (self.epsik[j] / self.t))
                weight = self.x[j] * self.m[k] * self.m[j] * sigma_ij**3
                one_minus_k = 1.0 - self.kij[k * n + j]
                s1_k += 2.0 * weight * e_ij * one_minus_k
                s2_k += 2.0 * weight * e_ij * e_ij * one_minus_k * one_minus_k

            eta_k = math.pi / 6.0 * AVOGADRO * self.m[k] * self.d[k] ** 3 / self.v
            d = (
                self.f()
                + self.rho * f_rho
                + eta_k * f_eta
                + f_m_bar * (self.m[k] - self.m_bar)
                + f_m_minus_1 * ((self.m[k] - 1.0) - self.m_minus_1)
                + f_s1 * s1_k
                + f_s2 * s2_k
            )
            out.append(d - math.log(z))
        return out


def _molar_volume(
    components: Sequence[_Component],
    kij: Sequence[float],
    x: Sequence[float],
    t: float,
    p: float,
    side: str,
) -> tuple[float, float, int]:
    """Solve ``P_calc(v) = p`` for ``(v, Z, iterations)``.

    NeqSim's ``PhasePCSAFTRahmat.calcVolume`` is a damped Newton on the volume:
    ``v += 0.9 (p - P_calc)/(dP_calc/dv)``, stopping at a relative step of ``1e-10`` or a
    hundred steps. The packing fraction reaches one at ``(pi/6) N_A md3``, where every
    hard-sphere term diverges, so that is the floor.

    The vapour branch is seeded from the ideal gas. **The liquid branch cannot be** - at a
    pressure where the isotherm has one root, that root is the vapour one, so a dilute seed
    converges there and stops - so it walks up from the floor and takes the lowest zero.
    """
    if not p > 0.0:
        raise OutOfRangeError("P", p, "the volume solve divides by the pressure")

    rt = R * t
    ideal = rt / p
    md3 = _State(components, kij, x, t, ideal).md3
    floor = math.pi / 6.0 * AVOGADRO * md3

    def h(v: float) -> float:
        return _State(components, kij, x, t, v).pressure_over_rt() * rt - p

    def brackets() -> list[tuple[float, float]]:
        start = floor * (1.0 + 1.0e-9)
        top = 1.0e4 * ideal
        steps = 128
        out: list[tuple[float, float]] = []
        previous = start
        previous_value = h(previous)
        for step in range(1, steps + 1):
            v = start * (top / start) ** (step / steps)
            value = h(v)
            if (previous_value < 0.0) != (value < 0.0):
                out.append((previous, v))
            previous, previous_value = v, value
        return out

    def bisect(lo: float, hi: float) -> float:
        sign = h(lo) < 0.0
        for _ in range(200):
            mid = 0.5 * (lo + hi)
            if (h(mid) < 0.0) == sign:
                lo = mid
            else:
                hi = mid
        return 0.5 * (lo + hi)

    if side == "vapour":
        v = ideal
        converged = None
        for step in range(100):
            state = _State(components, kij, x, t, v)
            # ``v += 0.9 (P - P_calc)/(dP_calc/dv)``, NeqSim's step: ``h`` is
            # ``P_calc - P``, so the sign is Newton's with the residual negated.
            denominator = state.d_pressure_over_rt_dv() * rt
            if denominator == 0.0:
                break
            delta = -0.9 * h(v) / denominator
            if not math.isfinite(delta) or v + delta <= floor:
                break
            relative = abs(delta) / v
            v += delta
            if relative < 1.0e-10:
                converged = step + 1
                break
        if converged is None:
            found = brackets()
            if not found:
                raise OutOfRangeError("P", p, f"the vapour branch has no zero at {t} K and {p} Pa")
            lo, hi = found[-1]
            return bisect(lo, hi), p * bisect(lo, hi) / rt, 0
        return v, p * v / rt, converged

    found = brackets()
    if not found:
        raise OutOfRangeError("P", p, f"the liquid branch has no zero at {t} K and {p} Pa")
    lo, hi = found[0]
    v = bisect(lo, hi)
    return v, p * v / rt, 0


def pcsaft_rahmat_phase(
    components: list[str],
    T: Q,
    P: Q,
    z: list[float],
    compressed_phase: str,
) -> PcsaftRahmatPhaseResult:
    """One PC-SAFT phase's state at a temperature, pressure and composition.

    Args:
        components: the substance names, resolved against the databank with their
            ``mSAFT``/``sigmaSAFT``/``epsikSAFT`` set and the ``KIJPCSAFT`` column.
        T: absolute temperature.
        P: absolute pressure.
        z: the mole fractions, checked rather than renormalised.
        compressed_phase: ``"liquid"`` or ``"vapour"``.

    Raises:
        InvalidInputError: if a component has no PC-SAFT set, if ``z`` is not a
            composition of the right length, or if ``compressed_phase`` is neither name.
        OutOfRangeError: if ``T`` or ``P`` is not positive, or the wanted branch has no
            root at this state.
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
    # `pcsaft_kij_for` gives the upper triangle; the matrix is symmetric, and both
    # directions are filled because the dispersion sums run over every ordered pair.
    pairs = pcsaft_kij_for(tuple(components))
    n = len(components)
    kij = [0.0] * (n * n)
    for (i, j), value in pairs.items():
        kij[i * n + j] = value
        kij[j * n + i] = value

    v, z_factor, _ = _molar_volume(resolved, kij, z, t_si, p_si, side)
    state = _State(resolved, kij, z, t_si, v)
    h_over_rt, s_over_r = state.departure()
    return PcsaftRahmatPhaseResult(
        z_factor=z_factor,
        ln_phi=tuple(state.ln_fugacity_coefficients()),
        v=ureg.Quantity(v, "m**3/mol"),
        h_res=ureg.Quantity(h_over_rt * R * t_si, "J/mol"),
        s_res=ureg.Quantity(s_over_r * R, "J/(mol*K)"),
        warnings=tuple(warnings),
    )


__all__ = ["pcsaft_rahmat_phase"]
