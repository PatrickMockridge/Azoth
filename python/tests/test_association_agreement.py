"""The Wertheim association kernel, against NeqSim's own numbers.

Two kernels are meant to be compared case by case, and these two cannot be called from
one process: `azoth_eos::association` is not exposed across the boundary, and it is not
exposed *because* nothing yet reaches it from Python - the model that will is
`eos.srk_cpa_phase`, and `test_cross_impl.py` compares that. Until it lands, the two
kernels are held together the way every other port here is: both are pinned to the same
external oracle's numbers, printed by `validation/neqsim/CpaSweep.java`.

    java -cp .:neqsim-f0c7436.jar CpaSweep 300 100 0.6

That is a weaker arrangement than comparing the two implementations, and deliberately so
rather than accidentally: what it cannot catch is a mistake the Rust kernel and NeqSim
make together, which is exactly what `test_cross_impl` exists for. The docstring says so
rather than leaving the reader to assume the stronger check.
"""

from __future__ import annotations

import math

import pytest

from azoth.eos import components
from azoth.eos.mixture import Mixture
from azoth.eos.reference._association import (
    ELLIOTT,
    Association,
    AssociationComponent,
    R,
    Rdf,
    SiteDerivatives,
    SiteScheme,
    SiteState,
    delta_nog,
)

#: Water and methanol as `SystemSrkCPA` builds them: the fitted SRK set, already in SI.
WATER = AssociationComponent(SiteScheme.FOUR_C, 16655.0, 0.0692)
METHANOL = AssociationComponent(SiteScheme.TWO_B, 24591.0, 0.0161)
COVOLUMES = (1.4515e-5, 3.0978e-5)
MOLES = (0.6, 0.4)
#: The state `CpaProbe` was first run at, and the volume it reported.
VOLUME = 2.620_326_748_288_29e-5
TEMPERATURE = 300.0


def _water_methanol() -> Association:
    return Association([WATER, METHANOL])


def test_the_kernel_reproduces_neqsims_cpa_probe() -> None:
    """The site fractions, the fugacity term and the Helmholtz energy, all at one state.

    The same four numbers the Rust kernel asserts, from the same probe: this is the
    agreement between the two kernels stated as agreement with the thing they were both
    ported from.
    """
    state = _water_methanol().solve(COVOLUMES, MOLES, VOLUME, TEMPERATURE)

    for site, want in ((0, 0.101_330_316_299_160), (4, 0.031_557_041_194_167_5)):
        assert state.fractions[site] == pytest.approx(want, rel=1.0e-9), f"xsite[{site}]"
    for i, want in ((0, -9.807_081_619_647_33), (1, -8.298_303_821_393_14)):
        assert state.ln_phi[i] == pytest.approx(want, rel=1.0e-9), f"dFCPAdN[{i}]"
    assert state.helmholtz_rt == pytest.approx(-6.793_473_163_493_31, rel=1.0e-9), "FCPA"


def test_the_site_fractions_satisfy_the_substitution() -> None:
    """`X_i = 1/(1 + S_i)` at the answer, checked against the definition rather than a probe.

    An oracle pins a state; this pins the *solve*, and it would hold at a state no probe
    was ever run at.
    """
    association = _water_methanol()
    state = association.solve(COVOLUMES, MOLES, VOLUME, TEMPERATURE)
    sites = association.site_count
    rdf = Rdf(COVOLUMES, MOLES, VOLUME)
    delta = association.delta_matrix(COVOLUMES, TEMPERATURE, rdf)
    for i in range(sites):
        s_i = (
            sum(
                MOLES[association.component_of_site(k)] * delta[i * sites + k] * state.fractions[k]
                for k in range(sites)
            )
            / VOLUME
        )
        assert state.fractions[i] == pytest.approx(1.0 / (1.0 + s_i), rel=1.0e-9)


def test_the_temperature_derivative_matches_a_finite_difference() -> None:
    """`d(A/RT)/dT`, against a central difference of the energy the solve returns.

    The value surface reads this and nothing else of the derivative family, so it is the
    one that has to be checked before a phase model can report an enthalpy departure.
    """
    association = _water_methanol()
    state = association.solve(COVOLUMES, MOLES, VOLUME, TEMPERATURE)
    analytic = association.temperature_derivative(COVOLUMES, MOLES, VOLUME, TEMPERATURE, state)

    step = 1.0e-4
    up = association.solve(COVOLUMES, MOLES, VOLUME, TEMPERATURE + step).helmholtz_rt
    down = association.solve(COVOLUMES, MOLES, VOLUME, TEMPERATURE - step).helmholtz_rt
    numerical = (up - down) / (2.0 * step)

    assert analytic == pytest.approx(numerical, rel=1.0e-6)


def test_a_two_a_component_associates_not_at_all() -> None:
    """A scheme whose sites cannot bond leaves every fraction at one and the energy at zero.

    The finding the port is built around: the table gives CO2 and H2S a fitted association
    energy under schemes whose charge vectors carry one sign, so NeqSim's product test is
    positive for every pair and no calculation reads the parameter.
    """
    association = Association([AssociationComponent(SiteScheme.TWO_A, 5000.0, 0.001_160_489)])
    state = association.solve((2.9e-5,), (1.0,), 4.0e-5, 300.0)

    assert state.fractions == (1.0, 1.0)
    assert state.helmholtz_rt == 0.0
    assert state.unbonded_sites == 0.0
    assert state.ln_phi == (0.0,)
    assert state.d_helmholtz_dv == 0.0


def test_a_mixture_without_sites_is_a_no_op() -> None:
    """A CPA mixture of light hydrocarbons on a zero-site basis, which is not an error."""
    association = Association([AssociationComponent(SiteScheme.NON_ASSOCIATING, 0.0, 0.0)] * 2)
    state = association.solve((3.0e-5, 3.0e-5), (0.5, 0.5), 6.0e-5, 300.0)

    assert state.fractions == ()
    assert state.helmholtz_rt == 0.0
    assert state.ln_phi == (0.0, 0.0)
    assert not association.has_bonds()


def test_a_negative_association_parameter_is_refused() -> None:
    """Both enter an exponential and a square root, so a negative one is not weaker."""
    from azoth.core.errors import OutOfRangeError

    with pytest.raises(OutOfRangeError):
        Association([AssociationComponent(SiteScheme.TWO_B, -1.0, 0.01)])


def test_the_elliott_rule_is_the_geometric_mean_of_the_pure_strengths() -> None:
    """`Delta_ij = sqrt(Delta_ii Delta_jj)` under the rule NeqSim defaults every pair to."""
    diagonal = delta_nog(WATER, COVOLUMES[0], WATER, COVOLUMES[0], ELLIOTT, TEMPERATURE)
    other = delta_nog(METHANOL, COVOLUMES[1], METHANOL, COVOLUMES[1], ELLIOTT, TEMPERATURE)
    cross = delta_nog(WATER, COVOLUMES[0], METHANOL, COVOLUMES[1], ELLIOTT, TEMPERATURE)

    assert cross == pytest.approx(math.sqrt(diagonal * other), rel=1.0e-12)


# ---------------------------------------------------------------------------
# The mixture layer: the same kernel reached through `_mixture_state`
# ---------------------------------------------------------------------------
#
# `_association` is the kernel; `_mixture_state` is where a *phase model* meets it -
# the fitted substitution, the associating volume root, and the association's share of
# `ln phi` and the departure enthalpy. The two are separate because a caller can reach
# the kernel directly and does not have to go through a cubic.


from azoth.eos.reference._mixture_state import (  # noqa: E402
    phase_derivatives,
    phase_state,
    reduced_parameters,
)


def _water_methanol_mixture() -> Mixture:
    """Water/methanol as `SystemSrkCPA` resolves it: SRK, the CPA column, association on."""
    return components.from_names(["water", "methanol"], eos="srk", associating=True)


def test_the_associating_liquid_root_matches_neqsim() -> None:
    """NeqSim's own `Z` at 300 K and 100 bar, which is the state the root was solved at.

    The association carries a pressure, so this is not the cubic's root: NeqSim's own `A`
    and `B` give a cubic root of 0.15229, and the association moves it to 0.10505 - 31%.
    """
    from azoth.eos.reference._mixture_state import phase_state, reduced_parameters

    mixture = _water_methanol_mixture()
    reduced = reduced_parameters(mixture, 300.0, 1.0e7)
    state = phase_state(reduced, mixture.kij, [0.6, 0.4], liquid=True)

    assert state.z == pytest.approx(0.105_050_962_879_418, rel=1.0e-9)


def test_the_substitution_replaces_the_cubics_attraction_and_covolume() -> None:
    """The fitted `a` and `b`, not the ones `Tc` and `Pc` derive.

    Water's fitted covolume is 1.4515e-5 m3/mol against `0.08664 R Tc/Pc`'s 2.11e-5, so a
    mixture that quietly used the cubic's would be a different fluid - and its alpha
    coefficient is fitted too, so the whole `a_i(T)` differs rather than one factor in it.
    """
    from azoth.eos.reference._mixture_state import reduced_parameters

    mixture = _water_methanol_mixture()
    associating = reduced_parameters(mixture, 356.0, 1.0e5)
    classical = components.from_names(["water", "methanol"], eos="srk")
    plain = reduced_parameters(classical, 356.0, 1.0e5)

    assert associating.b[0] != pytest.approx(plain.b[0])
    assert associating.a[0] != pytest.approx(plain.a[0])
    # `b_i` is `b * P/(R T)`, so this recovers the fitted covolume and checks it against
    # the published SI value rather than against the other call.
    recovered = associating.b[0] * R * 356.0 / 1.0e5
    assert recovered == pytest.approx(1.4515e-5, rel=1.0e-12)


def test_the_low_pressure_liquid_root_is_found_and_is_the_lower_one() -> None:
    """The regime the root finder used to get wrong, at 356 K and 1 bar.

    At low pressure the substituted cubic has **one** real root, so `z_min` and `z_max` are
    the same number and Newton from it lands on the vapour root. The liquid root exists only
    because the association's pressure creates it, a fraction of a per cent above the
    covolume, so it has to be found rather than seeded - and before the geometric walk it
    was not, which made every low-pressure CPA flash report a single phase, silently.
    """
    from azoth.eos.reference._mixture_state import phase_state, reduced_parameters

    mixture = _water_methanol_mixture()
    reduced = reduced_parameters(mixture, 356.0, 1.0e5)
    vapour = phase_state(reduced, mixture.kij, [0.6, 0.4], liquid=False)
    liquid = phase_state(reduced, mixture.kij, [0.6, 0.4], liquid=True)

    assert vapour.z == pytest.approx(0.949_619_265_899_294_9, rel=1.0e-12)
    assert liquid.z == pytest.approx(0.000_935_582_331_200_141_4, rel=1.0e-9)
    assert liquid.z < vapour.z / 100.0


def test_the_association_is_added_to_the_fugacity_and_to_the_enthalpy_together() -> None:
    """`ln phi` and `h_dep_rt` move together, because `s_dep_r` is the Gibbs identity.

    `s_dep_r` is `h_dep_rt - sum_i x_i ln phi_i`, not a second formula, so an association
    added to one and not the other would be folded silently into the entropy. This checks
    the identity rather than the two values.
    """
    from azoth.eos.reference._mixture_state import phase_state, reduced_parameters

    mixture = _water_methanol_mixture()
    reduced = reduced_parameters(mixture, 300.0, 1.0e7)
    state = phase_state(reduced, mixture.kij, [0.6, 0.4], liquid=True)

    assert state.s_dep_r == pytest.approx(
        state.h_dep_rt - sum(x * lp for x, lp in zip([0.6, 0.4], state.ln_phi, strict=True)),
        rel=1.0e-12,
    )
    # And the association is actually in there: a mixture that did not run it disagrees.
    classical = components.from_names(["water", "methanol"], eos="srk")
    plain_reduced = reduced_parameters(classical, 300.0, 1.0e7)
    plain = phase_state(plain_reduced, classical.kij, [0.6, 0.4], liquid=True)
    assert state.ln_phi[0] != pytest.approx(plain.ln_phi[0])


# ---------------------------------------------------------------------------
# The implicit derivative surface
# ---------------------------------------------------------------------------
#
# Every one of these is a derivative of a *solution*, so each is checked against a finite
# difference of the solve rather than against the derivation. That distinction is the
# whole reason this tranche's kernel is trustworthy and its earlier algebra was not: a
# consistently partial derivation satisfies a consistently partial identity, so the only
# check that finds a missing term is one that never sees the identity at all.

#: A volume and temperature at which the liquid branch is well conditioned.
_STATE = (2.620_326_748_288_29e-5, 300.0)


def _state_and_derivatives() -> tuple[Association, SiteState, SiteDerivatives]:
    """Water/methanol at the probe's state, solved, with its derivative surface."""
    association = _water_methanol()
    v, t = _STATE
    state = association.solve(COVOLUMES, MOLES, v, t)
    return association, state, association.derivatives(COVOLUMES, MOLES, v, t, state)


def test_the_site_fraction_volume_derivative_matches_the_solve() -> None:
    """`dX/dV`, against a central difference of the site fractions the solve returns."""
    association, _, derivatives = _state_and_derivatives()
    v, t = _STATE
    h = 1.0e-12

    for site in range(association.site_count):
        up = association.solve(COVOLUMES, MOLES, v + h, t).fractions[site]
        down = association.solve(COVOLUMES, MOLES, v - h, t).fractions[site]
        assert derivatives.d_fractions_dv[site] == pytest.approx(
            (up - down) / (2.0 * h), rel=1.0e-5
        ), f"dX_{site}/dV"


def test_the_second_volume_derivative_matches_a_finite_difference() -> None:
    """`d^2(A/RT)/dV^2`, against a difference of the *first* derivative.

    `d_helmholtz_dv` is carried on the solve, so differencing it is an independent voice
    where differencing the energy twice is not.
    """
    association, _, derivatives = _state_and_derivatives()
    v, t = _STATE
    h = 1.0e-9
    up = association.solve(COVOLUMES, MOLES, v + h, t).d_helmholtz_dv
    down = association.solve(COVOLUMES, MOLES, v - h, t).d_helmholtz_dv

    assert derivatives.d2_helmholtz_dv2 == pytest.approx((up - down) / (2.0 * h), rel=1.0e-6)


def test_the_mixed_composition_derivative_matches_a_finite_difference() -> None:
    """`d^2(A/RT)/dV dn_j`, the derivative that cost the most to get right.

    Its missing term was that `F_{X_c}` holds `X_a m_c Delta_ac / V`, and `m_c` *is* `n_b`
    when `c` is one of b's sites - so perturbing a mole number moves that factor inside the
    product as well as beside it. Three analytic routes disagreed on this number before a
    finite difference of the solve settled it.
    """
    association, _, derivatives = _state_and_derivatives()
    v, t = _STATE
    step = 1.0e-9

    for j in range(2):

        def energy(sign: float, j: int = j) -> float:
            moles = list(MOLES)
            moles[j] += sign * step
            return association.solve(COVOLUMES, moles, v, t).d_helmholtz_dv

        numerical = (energy(1.0) - energy(-1.0)) / (2.0 * step)
        assert derivatives.d2_helmholtz_dv_dn[j] == pytest.approx(numerical, rel=1.0e-5)


def test_the_mixed_temperature_derivative_matches_a_finite_difference() -> None:
    """`d^2(A/RT)/dV dT`, checked at two step sizes so convergence is visible.

    `T` enters through `Delta` alone, but `Delta'` carries the distribution function, so
    `F_{VT}` is not `-F_T/V` - the term a first reading of this derivative misses.
    """
    association, _, derivatives = _state_and_derivatives()
    v, t = _STATE

    for step in (1.0e-4, 1.0e-3):
        up = association.solve(COVOLUMES, MOLES, v, t + step).d_helmholtz_dv
        down = association.solve(COVOLUMES, MOLES, v, t - step).d_helmholtz_dv
        assert derivatives.d2_helmholtz_dv_dt == pytest.approx(
            (up - down) / (2.0 * step), rel=1.0e-6
        ), f"at step {step}"


def test_the_fugacity_composition_derivative_matches_its_definition() -> None:
    """`d ln phi_i/dn_j`, against the derivative of the fugacity term it comes from.

    The identity is `d(ln phi_i)/dn_j = sum_{A in i} X_A^(n_j)/X_A - (1/2) d(h calc_lngi_i)/dn_j`,
    and it is checked here against the site fractions the solve returns rather than against
    the assembled derivative, which is the only way a missing product-rule term shows.
    """
    association, _, derivatives = _state_and_derivatives()
    v, t = _STATE
    step = 1.0e-9

    for j in range(2):

        def fugacity(sign: float, j: int = j) -> tuple[float, ...]:
            moles = list(MOLES)
            moles[j] += sign * step
            return association.solve(COVOLUMES, moles, v, t).ln_phi

        up, down = fugacity(1.0), fugacity(-1.0)
        for i in range(2):
            numerical = (up[i] - down[i]) / (2.0 * step)
            assert derivatives.d_ln_phi_dn[i][j] == pytest.approx(numerical, rel=1.0e-5), (
                f"d ln phi_{i}/dn_{j}"
            )


def test_the_associating_derivative_surface_matches_the_phase_state() -> None:
    """`phase_derivatives` for an associating mixture, against finite differences of the state.

    The same oracle the Rust twin uses, and the one that matters: `phase_state` re-solves the
    associating root at every perturbed state, so it sees a term missing from the root's
    sensitivities - which is where the difference lives and where a check that shared the
    derivation could not look.

    The composition family is the one that caught a real defect here: an edit that added the
    association's contribution to `d_ln_phi_dn` silently failed to apply, and the surface
    then returned the *cubic's* derivative - agreeing with a cubic-only finite difference
    exactly, which is why nothing but a comparison against the full state could see it.
    """
    from azoth.eos import components as databank

    mixture = databank.from_names(["water", "methanol"], eos="srk", associating=True)
    t, p = (356.0, 1.0e5)
    x = [0.6, 0.4]
    reduced = reduced_parameters(mixture, t, p)
    state = phase_state(reduced, mixture.kij, x, liquid=False)
    surface = phase_derivatives(reduced, mixture.kij, x, state.z, temperature=t, pressure=p)

    for i in range(2):
        assert surface.ln_phi[i] == pytest.approx(state.ln_phi[i], rel=1.0e-12)

    def ln_phi_at(temperature: float, pressure: float, composition: list[float]) -> list[float]:
        return phase_state(
            reduced_parameters(mixture, temperature, pressure),
            mixture.kij,
            composition,
            liquid=False,
        ).ln_phi

    # Composition, at constant T and P, with the total held at one - the normalised partial.
    h = 1.0e-7
    for j in range(2):
        up, down = list(x), list(x)
        up[j] += h
        down[j] -= h
        up = [v / (1.0 + h) for v in up]
        down = [v / (1.0 - h) for v in down]
        up_ln_phi, down_ln_phi = ln_phi_at(t, p, up), ln_phi_at(t, p, down)
        for i in range(2):
            numerical = (up_ln_phi[i] - down_ln_phi[i]) / (2.0 * h)
            assert surface.d_ln_phi_dn[i][j] == pytest.approx(numerical, rel=1.0e-5), (
                f"d ln phi_{i}/dn_{j}"
            )

    step = 1.0e-3
    up, down = ln_phi_at(t + step, p, x), ln_phi_at(t - step, p, x)
    for i in range(2):
        assert surface.d_ln_phi_dt[i] == pytest.approx(
            (up[i] - down[i]) / (2.0 * step), rel=1.0e-6
        ), f"d ln phi_{i}/dT"

    step = 100.0
    up, down = ln_phi_at(t, p + step, x), ln_phi_at(t, p - step, x)
    for i in range(2):
        assert surface.d_ln_phi_dp[i] == pytest.approx(
            (up[i] - down[i]) / (2.0 * step), rel=1.0e-6
        ), f"d ln phi_{i}/dP"
