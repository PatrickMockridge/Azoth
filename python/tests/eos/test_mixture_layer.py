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

import pytest

import _helpers as h
from azoth import ureg
from azoth.eos import Component, mixture, pr_departure, pr_kappa
from azoth.eos.mixture import Mixture
from azoth.eos.reference._mixture_state import PhaseState, phase_state, reduced_parameters

Q = ureg.Quantity

PROPANE = Component(Q(369.83, "K"), Q(4_248_000.0, "Pa"), 0.1523)
BUTANE = Component(Q(425.12, "K"), Q(3_796_000.0, "Pa"), 0.2002)
METHANE = Component(Q(190.56, "K"), Q(4_599_200.0, "Pa"), 0.01142)


def methane_butane() -> Mixture:
    return mixture([METHANE, BUTANE], kij={(0, 1): 0.05})


def ternary() -> Mixture:
    return mixture([METHANE, PROPANE, BUTANE])


@pytest.mark.parametrize(
    ("component", "t_c", "p_pa"),
    [
        (PROPANE, 300.0, 1.0e6),
        (PROPANE, 350.0, 2.0e6),
        (BUTANE, 350.0, 1.0e6),
        (METHANE, 200.0, 3.0e6),
    ],
    ids=["propane-300", "propane-350", "butane-350", "methane-200"],
)
def test_the_departures_reduce_to_pr_departure_at_one_component(
    component: Component, t_c: float, p_pa: float
) -> None:
    """At one component the mixture's departure functions *are* ``eos.pr_departure``.

    The strongest check available, and an exact one rather than a tolerance: the
    mixture form is the pure form with ``psi`` replaced by a weighted average, and at
    ``N = 1`` that average is the single component's own ``psi``. ``h_dep_rt`` is
    asserted **bit-identical**; ``s_dep_r`` gets a tolerance of one ulp, because it is
    formed by subtracting the Gibbs term and the pure calc forms it by adding two
    others - a different summation order and nothing more.
    """
    fluid = mixture([component])
    reduced = reduced_parameters(fluid, t_c, p_pa)
    state = phase_state(reduced, fluid.kij, [1.0], liquid=False)

    kappa = pr_kappa(component.omega).kappa
    pure = pr_departure(
        reduced.a[0], reduced.b[0], state.z, kappa, t_c / component.Tc.to_base_units().magnitude
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

    Measured on methane/n-butane at 330 K, ``h_dep_rt`` is ``-1.926095e-04`` at 1 kPa
    and ``-1.926771e-03`` at 10 kPa, a ratio of **10.0035**; at 100 kPa and 1 MPa the
    ratios are 10.035 and 10.382, drifting up as the ideal-gas limit is left behind.
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
