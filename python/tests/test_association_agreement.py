"""The Wertheim association kernel, against NeqSim's own numbers.

Two kernels are meant to be compared case by case, and these two cannot be called from
one process: `azoth_eos::association` is not exposed across the boundary, and it is not
exposed *because* nothing yet reaches it from Python - the model that will is
`eos.srk_cpa_phase`, and `test_cross_impl.py` compares that. Until it lands, the two
kernels are held together the way every other port here is: both are pinned to the same
external oracle's numbers, printed by `validation/neqsim/CpaSweep.java`.

    java -cp .:neqsim-3.20.0.jar CpaSweep 300 100 0.6

That is a weaker arrangement than comparing the two implementations, and deliberately so
rather than accidentally: what it cannot catch is a mistake the Rust kernel and NeqSim
make together, which is exactly what `test_cross_impl` exists for. The docstring says so
rather than leaving the reader to assume the stronger check.
"""

from __future__ import annotations

import math

import pytest

from azoth.eos.reference._association import (
    ELLIOTT,
    Association,
    AssociationComponent,
    Rdf,
    SiteScheme,
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
