"""The model layer's own arithmetic, checked by reduction.

The mixture fugacity coefficient and the mixture departure functions have **no
registered spec**: the registry's inputs are scalars, and there is nowhere in it to put
a composition vector. That makes them the one piece of arithmetic in this package that
no kernel checks directly, and the discipline that replaces a spec is *reduction* -
each must reproduce a registered calculation where the registered calculation applies.

This file is where those reductions live. Asserting the mixture form against values it
produced itself would check nothing.
"""

from __future__ import annotations

import math

import pytest

import _helpers as h
from azoth import ureg
from azoth.eos import Component, from_names, mixture, pr_departure, pr_kappa
from azoth.eos.components import component
from azoth.eos.mixture import Mixture
from azoth.eos.reference._mixture_state import (
    PhaseState,
    ReducedParameters,
    criticality_matrix,
    helmholtz_energy,
    helmholtz_hessian,
    mixture_parameters,
    phase_derivatives,
    phase_state,
    phase_state_at,
    reduced_parameters,
)
from azoth.eos.reference.pr_z_factor import pr_z_factor

Q = ureg.Quantity

# Resolved through the databank rather than typed here, like every other fixture in
# this test tree. The three written out longhand had drifted from NeqSim's COMP.csv:
# methane at 0.01142 and 4 599 200 Pa, where the table says 0.0115 and 4 599 000.
PROPANE = component("propane")
BUTANE = component("n-butane")
METHANE = component("methane")


def methane_butane() -> Mixture:
    """The pair, with the interaction parameter NeqSim fits for it."""
    return from_names(["methane", "n-butane"])


def ternary() -> Mixture:
    """Three components, their interaction parameters from the same table."""
    return from_names(["methane", "propane", "n-butane"])


@pytest.mark.parametrize(
    ("substance", "t_c", "p_pa"),
    [
        (PROPANE, 300.0, 1.0e6),
        (PROPANE, 350.0, 2.0e6),
        (BUTANE, 350.0, 1.0e6),
        (METHANE, 200.0, 3.0e6),
    ],
    ids=["propane-300", "propane-350", "butane-350", "methane-200"],
)
def test_the_departures_reduce_to_pr_departure_at_one_component(
    substance: Component, t_c: float, p_pa: float
) -> None:
    """At one component the mixture's departure functions *are* ``eos.pr_departure``.

    The strongest check available, and an exact one rather than a tolerance: the
    mixture form is the pure form with ``psi`` replaced by a weighted average, and at
    ``N = 1`` that average is the single component's own ``psi``. ``h_dep_rt`` is
    asserted **bit-identical**; ``s_dep_r`` gets a tolerance of one ulp, because it is
    formed by subtracting the Gibbs term and the pure calc forms it by adding two
    others - a different summation order and nothing more.
    """
    fluid = mixture([substance])
    reduced = reduced_parameters(fluid, t_c, p_pa)
    state = phase_state(reduced, fluid.kij, [1.0], liquid=False)

    kappa = pr_kappa(substance.omega).kappa
    pure = pr_departure(
        reduced.a[0], reduced.b[0], state.z, kappa, t_c / substance.Tc.to_base_units().magnitude
    )

    assert state.psi_bar == reduced.psi[0], f"T={t_c}: psi_bar should be the component's own psi"
    assert state.h_dep_rt == pure.h_dep_rt, (
        f"T={t_c}: the mixture departure enthalpy should be bit-identical"
    )
    assert abs(state.s_dep_r - pure.s_dep_r) < 1e-15, (
        f"T={t_c}: the entropy departures differ by {state.s_dep_r - pure.s_dep_r:e}, "
        f"more than the one ulp a different summation order explains"
    )


def test_the_gibbs_identity_holds_for_a_mixture() -> None:
    """``h_dep_rt - s_dep_r = sum_i z_i ln phi_i``, exactly, and it defines ``s_dep_r``.

    ``G = H - TS`` makes the departure Gibbs energy the composition-weighted log of the
    fugacity coefficients, so a disagreement is a defect rather than a residual to be
    tolerated.
    """
    for fluid, t_c, p_pa, z in (
        (methane_butane(), 330.0, 2.5e6, [0.6, 0.4]),
        (methane_butane(), 300.0, 3.0e6, [0.1, 0.9]),
        (ternary(), 320.0, 2.0e6, [0.5, 0.3, 0.2]),
    ):
        reduced = reduced_parameters(fluid, t_c, p_pa)
        for liquid in (True, False):
            state = phase_state(reduced, fluid.kij, z, liquid=liquid)
            g_dep_rt = sum(zi * lp for zi, lp in zip(z, state.ln_phi, strict=True))
            h.assert_close(
                state.h_dep_rt - state.s_dep_r,
                g_dep_rt,
                1e-14,
                f"T={t_c}, z={z}, liquid={liquid}",
            )


def test_psi_bar_lies_between_the_components_psi() -> None:
    """It is a weighted average, so it cannot leave the range of what it averages.

    A property rather than a value, and it is what would fail if the weights were
    dropped, transposed, or summed over ``i`` alone.
    """
    fluid = ternary()
    for t_c, p_pa in ((320.0, 2.0e6), (350.0, 5.0e6)):
        reduced = reduced_parameters(fluid, t_c, p_pa)
        state = phase_state(reduced, fluid.kij, [0.5, 0.3, 0.2], liquid=False)
        assert min(reduced.psi) - 1e-15 <= state.psi_bar <= max(reduced.psi) + 1e-15, (
            f"T={t_c}: psi_bar is {state.psi_bar}, outside [{min(reduced.psi)}, {max(reduced.psi)}]"
        )


def test_the_departures_are_proportional_to_pressure_at_low_pressure() -> None:
    """At low pressure the departure goes to zero, and the *shape* is the check.

    Measured on methane/n-butane at 330 K, ``h_dep_rt`` is ``-1.961381e-04`` at 1 kPa
    and ``-1.962095e-03`` at 10 kPa, a ratio of **10.0036**; at 100 kPa and 1 MPa the
    ratios are 10.037 and 10.397, drifting up as the ideal-gas limit is left behind.
    Asserting the shape rather than an arbitrary smallness is what makes this a check
    rather than a tolerance.
    """
    fluid = methane_butane()

    def at(p_pa: float) -> PhaseState:
        reduced = reduced_parameters(fluid, 330.0, p_pa)
        return phase_state(reduced, fluid.kij, [0.6, 0.4], liquid=False)

    low, mid = at(1.0e3), at(1.0e4)
    for name in ("h_dep_rt", "s_dep_r"):
        ratio = getattr(mid, name) / getattr(low, name)
        h.assert_close(ratio, 10.0, 1e-2, f"{name} ratio over a tenfold pressure rise")


def compressibility_of(fluid: Mixture, t_c: float, p_pa: float, n: list[float]) -> float:
    """The mixture's own root at a state and composition.

    Computed rather than written down. Both tests below took a literal ``Z``, and it
    was a literal for a fluid that no longer exists: 0.8274482588400789 belonged to
    the hand-typed methane/n-butane with an illustrative ``kij`` of 0.05, and the
    fixture now resolves that pair through the databank. The energy identity below is
    a statement about the *state*, so a ``Z`` that is not that state's makes it fail by
    percent for a reason that looks like a defect.
    """
    reduced = reduced_parameters(fluid, t_c, p_pa)
    a_mix, b_mix = mixture_parameters(reduced.a, reduced.b, fluid.kij, n)
    return pr_z_factor(a_mix, b_mix).z_max


def test_the_energy_differentiates_to_the_fugacity_coefficient() -> None:
    """``d(A^R/RT)/dn_i`` must be ``ln phi_i + ln Z``.

    The right-hand side comes from :func:`phase_state_at`, which reaches it by
    differentiating a departure function; the left is a central difference of the
    energy directly. Two routes to one quantity, so the agreement is evidence rather
    than a restatement - and it is the identity that ties the critical point's
    machinery to the flash's.
    """
    fluid = methane_butane()
    reduced = reduced_parameters(fluid, 330.0, 2_500_000.0)
    n = [0.6, 0.4]
    z = compressibility_of(fluid, 330.0, 2_500_000.0, n)
    state = phase_state_at(reduced, fluid.kij, n, z)

    step = 1e-6
    for i in range(2):
        up, down = list(n), list(n)
        up[i] += step
        down[i] -= step
        gradient = (
            helmholtz_energy(reduced, fluid.kij, up, z)
            - helmholtz_energy(reduced, fluid.kij, down, z)
        ) / (2.0 * step)
        h.assert_close(
            gradient,
            state.ln_phi[i] + math.log(z),
            1e-9,
            f"component {i}: the energy's gradient against `ln phi + ln Z`",
        )


def test_the_hessian_is_the_second_derivative_of_the_energy() -> None:
    """Checked by central finite difference, at a step the two error terms agree on.

    The step is ``1e-3`` and the tolerance loose because a central second difference
    in ``f64`` cannot do better: truncation falls as ``h**2`` and round-off rises as
    ``eps/h**2``, and the measured minimum of the sum is about ``3e-8`` there. A
    tighter tolerance would be a claim about the difference quotient rather than
    about the Hessian - and it would still catch an error in any term, which moves an
    entry by order ``0.1``.
    """
    fluid = methane_butane()
    reduced = reduced_parameters(fluid, 330.0, 2_500_000.0)
    n, z = [0.6, 0.4], 0.8274482588400789
    hessian = helmholtz_hessian(reduced, fluid.kij, n, z)

    step = 1e-3

    def corner(di: float, dj: float, i: int, j: int) -> float:
        shifted = list(n)
        shifted[i] += di * step
        shifted[j] += dj * step
        return helmholtz_energy(reduced, fluid.kij, shifted, z)

    for i in range(2):
        for j in range(2):
            difference = (
                corner(1.0, 1.0, i, j)
                - corner(1.0, -1.0, i, j)
                - corner(-1.0, 1.0, i, j)
                + corner(-1.0, -1.0, i, j)
            ) / (4.0 * step * step)
            h.assert_close(hessian[i][j], difference, 1e-6, f"H[{i}][{j}] against the difference")


def test_the_hessian_and_the_criticality_matrix_are_exactly_symmetric() -> None:
    """Exactly, not closely.

    The expression is symmetric term by term, so the two must agree bit for bit. It
    matters because the eigenvalues are only real if the matrix is symmetric, and an
    asymmetry that appeared here would stay invisible until it produced a complex
    eigenvector somewhere downstream.
    """
    fluid = methane_butane()
    reduced = reduced_parameters(fluid, 330.0, 2_500_000.0)
    n, z = [0.6, 0.4], 0.8274482588400789

    hessian = helmholtz_hessian(reduced, fluid.kij, n, z)
    assert hessian[0][1] == hessian[1][0], "the Hessian is not symmetric"
    matrix = criticality_matrix(reduced, fluid.kij, n, z)
    assert matrix[0][1] == matrix[1][0], "Q is not symmetric"


#: The triple-root solution of the Peng-Robinson cubic, which this library does *not*
#: ship - see `eos.pr_alpha_ab`'s assumptions and the constants' own comments.
#:
#: Named here because these two tests are statements about the cubic rather than about
#: the constants: `Q` vanishes, and is minimised, at the critical point of a cubic whose
#: `Omega` pair satisfies the triple-root condition. NeqSim's pair does not, so the
#: shipped cubic's critical point is 4.8e-5 away from `Tr = Pr = 1` and its `Q` there is
#: 1.37e-4 rather than zero.
TRIPLE_ROOT_A = 0.4572355289213822
TRIPLE_ROOT_B = 0.07779607390388846


def triple_root_state() -> ReducedParameters:
    """A one-component reduced state at the cubic's own critical point.

    Built directly rather than through a substance: the claim is about ``(A, B, Z)``,
    which is the point the other tests in this file make about the construction
    carrying no dimensioned quantity.
    """
    return ReducedParameters(
        a=[TRIPLE_ROOT_A],
        b=[TRIPLE_ROOT_B],
        psi=[0.0],
        psi_t=[0.0],
        warnings=[],
    )


def test_the_criticality_matrix_vanishes_at_the_cubics_critical_point() -> None:
    """The check a port cannot inherit, against an answer known in closed form.

    Heidemann & Khalil's first condition is that the smallest eigenvalue of ``Q``
    reaches zero. For one component ``Q`` is ``1 x 1`` and equals ``H_11 + 1``, and
    the critical point of a Peng-Robinson cubic is analytic: ``Tr = Pr = 1``,
    ``A = Omega_a``, ``B = Omega_b``, ``Z_c = (1 - Omega_b)/3``. So the whole
    construction - the constant-volume Hessian, the ideal part, the scaling - is
    checked against a closed form rather than against another implementation.

    **At the triple-root pair, which is what makes the closed form hold.** This test
    used to reach that state through a real substance at its own tabulated ``Tc`` and
    ``Pc``, which was the same state only for as long as the library shipped the
    triple-root ``Omega`` pair. It does not - it ships NeqSim's, deliberately, and
    `eos.pr_alpha_ab`'s assumptions record why - so the substance route now lands
    4.8e-5 off the cubic's critical point and ``Q`` there is 1.37e-4. Asserting the
    closed form requires the pair the closed form is about, and building the reduced
    state directly is the honest way to ask for it.
    """
    value = criticality_matrix(triple_root_state(), ((0.0,),), [1.0], (1.0 - TRIPLE_ROOT_B) / 3.0)[
        0
    ][0]
    assert abs(value) < 1e-14, f"Q is {value!r} at the critical point, and it must be zero"


def test_the_shipped_omegas_leave_a_critical_point_residue() -> None:
    """How far the shipped pair puts the cubic's critical point from ``Tr = Pr = 1``.

    Recorded rather than asserted loosely, because it is the one number that says what
    carrying NeqSim's literals costs: at the shipped pair the same construction gives
    ``Q = 1.37e-4`` at ``(A, B) = (Omega_a, Omega_b)``, and the cubic is at a triple
    root at ``Tr = 1 + 4.8e-5`` instead. Bounded on both sides so a change to the
    constants is visible here as well as in `eos.pr_alpha_ab`.
    """
    from azoth.eos.reference.pr_alpha_ab import OMEGA_A, OMEGA_B

    state = ReducedParameters(a=[OMEGA_A], b=[OMEGA_B], psi=[0.0], psi_t=[0.0], warnings=[])
    value = criticality_matrix(state, ((0.0,),), [1.0], (1.0 - OMEGA_B) / 3.0)[0][0]
    assert 1e-5 < value < 1e-3, f"Q at the shipped pair's Tr = Pr = 1 is {value!r}"


def test_the_criticality_matrix_is_minimised_at_the_critical_temperature() -> None:
    """The residue is a minimum, and clearly positive on both sides of it.

    Along ``P = Pc`` the critical point is the temperature at which ``Q`` reaches its
    least value, and it is orders of magnitude larger either side. Without this, a
    ``Q`` that was uniformly near zero everywhere would pass the test above - which is
    what a missing diagonal term would give.

    **The minimum is not zero, because the shipped ``Omega`` pair is not the
    triple-root one** - measured at 2.5e-3 for propane at ``T/Tc = 1``. What survives
    the port is where the minimum sits (``T/Tc = 1``, to five figures) and how sharply
    it rises away from it.
    """
    from azoth.eos.reference.pr_z_factor import pr_z_factor

    fluid = mixture([PROPANE])
    t_c = PROPANE.Tc.to_base_units().magnitude
    p_c = PROPANE.Pc.to_base_units().magnitude

    def q_at(offset: float) -> float:
        reduced = reduced_parameters(fluid, t_c * offset, p_c)
        roots = pr_z_factor(reduced.a[0], reduced.b[0])
        root = roots.z_min if offset < 1.0 else roots.z_max
        return criticality_matrix(reduced, fluid.kij, [1.0], root)[0][0]

    at_critical = q_at(1.0)
    assert 1e-3 < at_critical < 1e-2, f"Q at T/Tc = 1 is {at_critical!r}"

    for offset in (0.98, 0.99, 0.999, 1.001, 1.01, 1.02):
        value = q_at(offset)
        assert value > 1e-2, (
            f"T/Tc = {offset}: Q is {value}, but away from the critical point it must "
            f"be clearly positive"
        )
        assert value > at_critical, (
            f"T/Tc = {offset}: Q is {value}, below the value at the critical point - "
            f"the critical temperature is not the minimum"
        )


def test_the_helmholtz_layer_matches_its_recorded_values() -> None:
    """A change detector: the Python results against recorded values, at 1e-12.

    The tolerance is not bit-equality because ``ln`` is not correctly rounded; the
    arithmetic other than that is ``+ - * /`` and ``sqrt``, which would have permitted
    one.

    This is not a cross-language test. It calls one implementation.
    """
    fluid = methane_butane()
    reduced = reduced_parameters(fluid, 330.0, 2_500_000.0)
    n = [0.6, 0.4]
    z = compressibility_of(fluid, 330.0, 2_500_000.0, n)

    h.assert_close(
        helmholtz_energy(reduced, fluid.kij, n, z),
        -0.18925799467861787,
        1e-12,
        "the residual Helmholtz energy",
    )
    hessian = helmholtz_hessian(reduced, fluid.kij, n, z)
    for (i, j), wanted in {
        (0, 0): -0.07053560045586016,
        (0, 1): -0.28355937321520475,
        (1, 0): -0.28355937321520475,
        (1, 1): -1.0641939064423525,
    }.items():
        h.assert_close(hessian[i][j], wanted, 1e-12, f"the Hessian at {i},{j}")

    matrix = criticality_matrix(reduced, fluid.kij, n, z)
    for (i, j), wanted in {
        (0, 0): 0.9576786397264839,
        (0, 1): -0.1389151552321342,
        (1, 0): -0.1389151552321342,
        (1, 1): 0.574322437423059,
    }.items():
        h.assert_close(matrix[i][j], wanted, 1e-12, f"Q at {i},{j}")


def _ln_phi_of(
    mixture: Mixture, reduced: ReducedParameters, n: list[float], *, liquid: bool
) -> list[float]:
    """``ln phi_i`` as a function of mole numbers, which is what the surface differentiates.

    The analytic derivative is the one of the *intensive* function ``ln phi_i(T, P, x)``
    with ``x = n / sum(n)``, so the difference quotient has to be taken of the same
    function. Calling :func:`phase_state` with a raw perturbed vector would differentiate
    a different function - one whose ``A`` is ``sum sum n_i n_j A_ij`` rather than the
    same over ``N**2`` - and the two agree only on the sum-to-one surface.
    """
    total = sum(n)
    x = [value / total for value in n]
    return phase_state(reduced, mixture.kij, x, liquid=liquid).ln_phi


def test_the_composition_derivative_is_the_fugacity_coefficients_own() -> None:
    """A central difference of ``ln phi_i(n)`` at the ``N = 1`` the surface is stated at."""
    mix = methane_butane()
    t, p, n = 330.0, 2_500_000.0, [0.6, 0.4]
    reduced = reduced_parameters(mix, t, p)
    state = phase_state(reduced, mix.kij, n, liquid=True)
    d = phase_derivatives(reduced, mix.kij, n, state.z, temperature=t, pressure=p)

    step = 1e-6
    for j in range(2):
        up, down = list(n), list(n)
        up[j] += step
        down[j] -= step
        above = _ln_phi_of(mix, reduced, up, liquid=True)
        below = _ln_phi_of(mix, reduced, down, liquid=True)
        for i in range(2):
            difference = (above[i] - below[i]) / (2.0 * step)
            assert abs(d.d_ln_phi_dn[i][j] - difference) < 1e-8, (
                f"d ln phi_{i} / d n_{j}: analytic {d.d_ln_phi_dn[i][j]}, difference {difference}"
            )


def test_the_state_derivatives_are_the_fugacity_coefficients_own() -> None:
    """Both at constant composition, which is what the surface states them at."""
    mix = methane_butane()
    t, p, n = 330.0, 2_500_000.0, [0.6, 0.4]

    def at(temperature: float, pressure: float) -> list[float]:
        return _ln_phi_of(mix, reduced_parameters(mix, temperature, pressure), n, liquid=False)

    reduced = reduced_parameters(mix, t, p)
    roots = pr_z_factor(*mixture_parameters(reduced.a, reduced.b, mix.kij, n))
    d = phase_derivatives(reduced, mix.kij, n, roots.z_max, temperature=t, pressure=p)

    (step_t, step_p) = (1e-3, 1e-2)
    above, below = at(t + step_t, p), at(t - step_t, p)
    for i in range(2):
        difference = (above[i] - below[i]) / (2.0 * step_t)
        assert abs(d.d_ln_phi_dt[i] - difference) < 1e-7, (
            f"d ln phi_{i} / dT: analytic {d.d_ln_phi_dt[i]}, difference {difference}"
        )
    above, below = at(t, p + step_p), at(t, p - step_p)
    for i in range(2):
        difference = (above[i] - below[i]) / (2.0 * step_p)
        assert abs(d.d_ln_phi_dp[i] - difference) < 1e-9, (
            f"d ln phi_{i} / dP: analytic {d.d_ln_phi_dp[i]}, difference {difference}"
        )


@pytest.mark.parametrize(
    ("t", "p", "n", "liquid"),
    [
        (330.0, 2_500_000.0, [0.6, 0.4], True),
        (330.0, 2_500_000.0, [0.6, 0.4], False),
        (300.0, 3_000_000.0, [0.1, 0.9], True),
        (400.0, 1_000_000.0, [0.5, 0.5], False),
    ],
    ids=["330-liquid", "330-vapour", "300-liquid", "400-vapour"],
)
def test_the_composition_derivative_obeys_gibbs_duhem_and_is_symmetric(
    t: float, p: float, n: list[float], liquid: bool
) -> None:
    """Two structural properties that no difference quotient can show.

    At constant temperature and pressure ``sum_i n_i d ln phi_i = 0``, so every column
    sums to zero against the composition; and the matrix is ``(1/RT)`` times a Gibbs
    energy's second derivative, so it is symmetric. A matrix can satisfy either and fail
    the other, and together they pin the whole block rather than a sum of it.
    """
    mix = methane_butane()
    reduced = reduced_parameters(mix, t, p)
    state = phase_state(reduced, mix.kij, n, liquid=liquid)
    d = phase_derivatives(reduced, mix.kij, n, state.z, temperature=t, pressure=p)

    # Absolute, not relative: the targets are zero and a symmetry, so a relative test
    # against them would be measuring the round-off of a subtraction of near-equals.
    for j in range(len(n)):
        total = sum(n[i] * d.d_ln_phi_dn[i][j] for i in range(len(n)))
        assert abs(total) < 1e-9, f"Gibbs-Duhem column {j} sums to {total}, not to zero"
    for i in range(len(n)):
        for j in range(len(n)):
            gap = abs(d.d_ln_phi_dn[i][j] - d.d_ln_phi_dn[j][i])
            assert gap < 1e-12, f"symmetry {i},{j} is off by {gap}"
