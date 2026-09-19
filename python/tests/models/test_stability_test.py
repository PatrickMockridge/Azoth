"""Spec-driven tests for the ``eos.stability_test`` model.

The spec's cases pin the two implementations to each other. These tests are for the
things a case *cannot* say: the cross-model check against ``eos.pt_flash`` - which is
the scientific content of the pair - the identities that hold at any answer, and the
two decisions that were made rather than transcribed, each of which is a wrong answer
that looks like a right one.

``w`` is asserted here and not in the Rust test, and the reason is a generator gap
rather than a choice: a case's ``expected`` can carry scalars and vectors, and the
generated ``TestCase`` has no field for a *matrix* one, so the ``w`` block never
reaches Rust. This file reads the raw spec, so it is where the compositions are pinned
exactly - and the cross-backend comparison below carries them to the Rust side, since
both backends' ``w`` must agree on every case.
"""

from __future__ import annotations

import math
from functools import partial
from typing import Any

import pytest

import _helpers as h
from azoth import _models_gen, ureg, use_backend
from azoth.core.errors import InvalidInputError, OutOfRangeError, SolverNotConvergedError
from azoth.core.result import StabilityTestResult, StabilityVerdict
from azoth.eos import (
    Mixture,
    component,
    from_names,
    mixture,
    pr_z_factor,
    pt_flash,
    stability_test,
)
from azoth.eos.reference._mixture_state import (
    ReducedParameters,
    helmholtz_energy,
    mixture_parameters,
    phase_state,
    phase_state_at,
    reduced_parameters,
)

MODEL_ID = "eos.stability_test"
Q = ureg.Quantity

SPEC = _models_gen.model(MODEL_ID)
CASES = SPEC["cases"]

# The substances the cases and the identities below use, resolved through the
# databank rather than typed here. A `Component` written out longhand is a second
# copy of NeqSim's table, and this one had drifted: its methane was 0.01142 and
# 4 599 200 Pa, where COMP.csv says 0.0115 and 4 599 000.
METHANE = component("methane")
BUTANE = component("n-butane")
PROPANE = component("propane")


def methane_butane() -> Mixture:
    """The methane/n-butane pair, with the interaction parameter NeqSim fits for it.

    Resolved by name rather than assembled here, so the pair a sweep runs and the pair
    a case runs are the same fluid. The `kij` used to be an "illustrative" 0.05 - a
    number in a test file that described no fluid, and that disagreed with `INTER.csv`.
    """
    return from_names(["methane", "n-butane"])


def call(case: dict[str, Any]) -> StabilityTestResult:
    """Run one case declared in the model spec.

    Through :func:`_helpers.model_kwargs`, which is the one place a case's
    declared inputs become arguments: it resolves `components` against the
    databank and hands over the mixture and the ideal-gas model the function
    takes. A hand-built mixture here would be a second fluid, described by the
    case file rather than by NeqSim's tables.
    """
    return stability_test(**h.model_kwargs(SPEC, case["inputs"]))


def gibbs_at_root(
    reduced: ReducedParameters,
    kij: tuple[tuple[float, ...], ...],
    z: list[float],
    compressibility: float,
) -> float:
    """``A^R/RT - ln Z + Z`` at one root, which is ``G/RT`` up to a constant.

    The comparison the model places the feed by, written out here so a test can ask
    which root each *form* of it picks. The constant it drops is the same at both
    roots - same ``T``, ``P`` and composition - which is the whole reason this form
    is the right one.
    """
    return (
        helmholtz_energy(reduced, kij, z, compressibility)
        - math.log(compressibility)
        + compressibility
    )


def test_the_energy_differentiates_to_the_fugacity_it_is_placed_by() -> None:
    """`d(A^R/RT)/dn_i = ln phi_i + ln Z`, finite-differenced, for an associating fluid.

    This is the identity the whole surface rests on: the energy and the fugacity
    coefficient are one function twice, and the root is chosen by the energy while the
    trials are placed by the fugacity. An association in one and not the other would put
    them on different fluids.

    **It is the sharp test of the association's branch here, and the plain root-choice
    test is not**: at 356 K and 1 bar the cubic-only energy still orders the two roots the
    same way, so dropping the branch changes every value and no decision. The identity
    notices: without the branch the two sides disagree by a factor of three.
    """
    from azoth.eos import components as databank

    fluid = databank.from_names(["water", "methanol"], eos="srk", associating=True)
    reduced = reduced_parameters(fluid, 356.0, 1.0e5)
    z = [0.6, 0.4]
    step = 1.0e-6
    for liquid in (True, False):
        state = phase_state(reduced, fluid.kij, z, liquid=liquid)
        for i in range(len(z)):
            up, down = list(z), list(z)
            up[i] += step
            down[i] -= step
            # `n` at a fixed compressibility is a perturbation at fixed volume: this
            # function's volume is `Z R T/P` whatever the moles.
            difference = (
                helmholtz_energy(reduced, fluid.kij, up, state.z)
                - helmholtz_energy(reduced, fluid.kij, down, state.z)
            ) / (2.0 * step)
            expected = state.ln_phi[i] + math.log(state.z)
            assert abs(difference / expected - 1.0) < 1.0e-7, (
                f"{'liquid' if liquid else 'vapour'} root, component {i}: the energy "
                f"differentiates to {difference} where the fugacity is {expected}"
            )


def test_the_root_the_energy_picks_is_the_root_physics_picks() -> None:
    """`gibbs_at_root` must place the feed on the same root its own Gibbs energy does.

    The model chooses the feed's root by comparing `A^R/RT - ln Z + Z` at the two cubic
    roots, and everything downstream - the tangent-plane distances, the trial phases,
    the verdict - is measured from that root. So a root placed wrong is a verdict
    measured from the wrong fluid.

    **It is a guard rather than a pin on the association's branch.** Measured: dropping
    the branch from `helmholtz_energy` leaves this passing at both states here, because
    the cubic-only energy still orders the two roots the same way - it changes every
    value and no decision. What it does catch is a change to `gibbs_at_root` or to the
    energy that ever *does* reverse them, which is the mistake the association's absence
    was one state away from.

    The control is the plain cubic, where the two have always agreed.
    """
    from azoth.eos import components as databank

    for label, fluid, temperature, pressure in (
        (
            "water/methanol, associating",
            databank.from_names(["water", "methanol"], eos="srk", associating=True),
            356.0,
            1.0e5,
        ),
        (
            "methane/n-butane, a plain cubic",
            databank.from_names(["methane", "n-butane"]),
            300.0,
            2.0e6,
        ),
    ):
        reduced = reduced_parameters(fluid, temperature, pressure)
        z = [0.6, 0.4]
        by_energy = {}
        by_physics = {}
        for root, liquid in (("liquid", True), ("vapour", False)):
            state = phase_state(reduced, fluid.kij, z, liquid=liquid)
            by_energy[root] = gibbs_at_root(reduced, fluid.kij, z, state.z)
            # The physical Gibbs energy, up to the standard state, which is the same at
            # both roots - the same reason `gibbs_at_root` may drop its constant.
            by_physics[root] = sum(
                z_i * (math.log(z_i) + ln_phi_i)
                for z_i, ln_phi_i in zip(z, state.ln_phi, strict=True)
            )
        chosen = min(by_energy, key=lambda root: by_energy[root])
        truth = min(by_physics, key=lambda root: by_physics[root])
        assert chosen == truth, (
            f"{label}: the energy places the feed on the {chosen} root and the physical "
            f"Gibbs energy on the {truth}. {by_energy} against {by_physics}"
        )


def assert_tm(actual: float, expected: float, tolerance: float, context: str) -> None:
    """Compare a trial distance, where an expectation of zero has no scale.

    ``h.assert_close`` is *relative*, deliberately - the numbers in this registry
    span nine orders of magnitude - and a relative comparison against an expected
    zero divides by the smallest positive double. A trivial trial's distance is
    exactly that: ``tm = 1 - sum(W)`` where ``sum(W)`` is one to rounding, so the
    distance is ``0.0`` or ``-4.4e-16`` depending on which way the last ``exp``
    rounded - and the spec's second case records one on each side. Comparing two
    such numbers relatively is comparing round-off.

    So an expectation at or below ``1e-12`` - which no non-trivial distance is - is
    compared absolutely, at the case's own tolerance.
    """
    if abs(expected) <= 1e-12:
        assert abs(actual) <= tolerance, (
            f"{context}: got {actual}, expected a distance of zero to within {tolerance:e}"
        )
    else:
        h.assert_close(actual, expected, tolerance, context)


@pytest.mark.parametrize("case", CASES, ids=lambda c: c["id"])
def test_spec_case(case: dict[str, Any]) -> None:
    """Run one case declared in the model spec."""
    result = call(case)
    tolerance = case["tolerance"]
    expected = case["expected"]

    assert len(result.tm) == len(expected["tm"]), f"{case['id']}: tm length"
    for i, (got, want) in enumerate(zip(result.tm, expected["tm"], strict=True)):
        assert_tm(got, want, tolerance, f"{case['id']} (tm[{i}])")

    assert result.iterations == tuple(expected["iterations"]), f"{case['id']}: iterations"

    # The compositions: two rows, always, each a composition.
    assert len(result.w) == 2, f"{case['id']}: two trials, so two rows"
    assert len(expected["w"]) == 2, f"{case['id']}: the case should declare two rows"
    for i, (row, want) in enumerate(zip(result.w, expected["w"], strict=True)):
        assert len(row) == len(want), f"{case['id']}: w[{i}] length"
        for j, (got, want_value) in enumerate(zip(row, want, strict=True)):
            h.assert_close(got, want_value, tolerance, f"{case['id']} (w[{i}][{j}])")
        assert abs(sum(row) - 1.0) <= 1e-12, f"{case['id']}: w[{i}] is not a composition"

    # The verdict is an enum, so it cannot be written into `expected` - it is
    # asserted in the cross-model test below, which is where it means anything.
    h.assert_consistent(result, case["id"])


def test_every_case_ran() -> None:
    """Guard against a spec edit that silently removes every case."""
    assert len(CASES) >= 2, f"expected several cases, found {len(CASES)}"


def test_the_verdict_agrees_with_the_flash_phase_on_every_state() -> None:
    """The scientific content of the pair: this model's verdict against the flash's.

    A feed a flash splits into two phases is unstable as a single phase - that is
    what the split *means* - and a feed it reports as one phase is not. So over a
    sweep the verdict must be ``UNSTABLE`` exactly where the phase is ``two_phase``,
    and ``STABLE`` everywhere else: ``all_liquid``, ``all_vapour`` and ``trivial``
    are all single-phase readings.

    ``trivial`` is the one worth stating. The flash reached ``x = y = z`` and cannot
    say which single phase the feed is; this model *can*, and says it is one phase.
    The two are consistent rather than redundant, which is why the pair exists.

    Nothing here is a value comparison: the two answers are enums, computed by
    separate iterations over separate equations, and the claim is that they never
    contradict each other.
    """
    fluids = (
        methane_butane(),
        mixture([PROPANE, BUTANE], kij={(0, 1): 0.02}),
    )

    splits = 0
    checked = 0
    for fluid in fluids:
        for t_c in (250.0, 280.0, 300.0, 330.0, 350.0, 400.0):
            for p_pa in (1e5, 1e6, 3e6, 1e7, 2e7, 5e7):
                for z in ([0.1, 0.9], [0.5, 0.5], [0.9, 0.1]):
                    try:
                        flash = pt_flash(fluid, T=Q(t_c, "K"), P=Q(p_pa, "Pa"), z=z)
                        verdict = stability_test(fluid, T=Q(t_c, "K"), P=Q(p_pa, "Pa"), z=z).verdict
                    except SolverNotConvergedError:
                        continue
                    checked += 1

                    expected = StabilityVerdict.UNSTABLE
                    if flash.phase.value == "two_phase":
                        splits += 1
                    else:
                        expected = StabilityVerdict.STABLE
                    assert verdict is expected, (
                        f"T={t_c}, P={p_pa}, z={z}: the flash says `{flash.phase.value}` and "
                        f"this model says `{verdict.value}`. A state that splits is unstable "
                        f"as a single phase, and one that does not is not - a disagreement "
                        f"here is a contradiction, not a tolerance"
                    )

    assert splits > 20, f"the sweep should find plenty of two-phase states, found {splits}"
    assert checked > 150, f"the sweep should cover plenty of states, covered {checked}"


def test_a_trivial_flash_is_not_a_stable_feed() -> None:
    """The contrast that is the model's reason to exist, on the spec's second case.

    At 430 K and 60 bar the flash converges to ``x = y = z`` and reports ``trivial``,
    which says the feed is single phase and *not which one* - the K-values straddled
    one throughout, so nothing in that model ever proved what the feed was. This model
    answers the question it could not ask, and answers ``stable``.

    Both trials converge to the feed here, so the verdict is the weakest kind of
    ``stable``: two trials that found nothing, not two trials that separated. The test
    asserts the shape of that - both distances at zero to rounding, and therefore the
    threshold that decides ``stable`` has to tolerate a *negative* zero, which is why
    it is ``-1e-8`` and not zero.
    """
    fluid = methane_butane()
    z = [0.6, 0.4]
    flash = pt_flash(fluid, T=Q(430.0, "K"), P=Q(60.0, "bar"), z=z)
    assert flash.phase.value == "trivial", "the contrast needs a flash that cannot answer"

    result = stability_test(fluid, T=Q(430.0, "K"), P=Q(60.0, "bar"), z=z)
    assert result.verdict is StabilityVerdict.STABLE
    for i, distance in enumerate(result.tm):
        assert abs(distance) <= 1e-8, (
            f"trial {i} reached a stationary point at tm = {distance}, so this is not "
            f"the two-trivial-trials state the case describes"
        )
    # The *sign* of that round-off is deliberately not asserted. It says which way the
    # last `exp` rounded, and this test used to pin it as negative - a measurement of
    # the old illustrative `kij` rather than a property of the model, and one the
    # databank's fitted parameter rounds the other way. What the threshold has to
    # tolerate is a distance below zero, and that a run ever produces one is a fact
    # about single-precision cancellation rather than about 430 K.


def test_the_feed_is_placed_on_its_lower_gibbs_root() -> None:
    """The decision this model is most likely to get quietly wrong.

    The comparison is ``A^R/RT - ln Z + Z`` and **not** ``A^R/RT + Z``: the ideal
    part of ``G/RT`` carries ``-ln Z``, and ``V = Z R T / P`` differs between the two
    roots, so only ``-ln Z`` separates them at one ``(T, P, n)``.

    Pure methane at 150 K and 1 bar is superheated vapour - its saturation pressure
    there is about 10 bar - so the vapour-like root is the state the feed is in. The
    two comparisons disagree about which root that is, and the numbers below are the
    spec's measurement of it. The second half is the consequence rather than the
    rule: with the right comparison the feed is a stable single phase, and with the
    wrong one both trials measure their distance from a state the feed is not in and
    the model calls a plain vapour unstable.
    """
    fluid = mixture([METHANE], kij={})
    reduced = reduced_parameters(fluid, 150.0, 100_000.0)
    z = [1.0]
    a_mix, b_mix = mixture_parameters(reduced.a, reduced.b, fluid.kij, z)
    roots = pr_z_factor(a_mix, b_mix)

    def forms(compressibility: float) -> tuple[float, float, float]:
        residual = helmholtz_energy(reduced, fluid.kij, z, compressibility)
        return (
            residual,
            residual + compressibility,
            residual - math.log(compressibility) + compressibility,
        )

    liquid_residual, liquid_naive, liquid_right = forms(roots.z_min)
    vapour_residual, vapour_naive, vapour_right = forms(roots.z_max)

    h.assert_close(liquid_residual, -2.5578138051194212, 1e-9, "liquid A^R/RT")
    h.assert_close(vapour_residual, -0.015548927641895367, 1e-9, "vapour A^R/RT")
    h.assert_close(liquid_right, 3.1458414529963945, 1e-9, "liquid A^R/RT - ln Z + Z")
    h.assert_close(vapour_right, 0.9845725795485407, 1e-9, "vapour A^R/RT - ln Z + Z")

    assert liquid_naive < vapour_naive, (
        "the wrong comparison is only worth recording if it picks the other root"
    )
    assert vapour_right < liquid_right, (
        "A^R/RT - ln Z + Z must pick the vapour root for a superheated vapour"
    )

    # The consequence, on **both** backends: a superheated vapour is a stable single
    # phase, and the trial seeded on its own side lands on the feed exactly.
    #
    # Both, because this state is the one that tells the two comparisons apart, and
    # each implementation chooses the root independently - a reference that had the
    # wrong form would place the feed on the liquid root here and report `unstable`
    # at about `-7.7`, while the Rust side reported `stable` at zero. A test on the
    # default backend alone would see only one of them.
    with use_backend("python"):
        from_python = stability_test(fluid, T=Q(150.0, "K"), P=Q(1.0, "bar"), z=z)
    with use_backend("rust"):
        from_rust = stability_test(fluid, T=Q(150.0, "K"), P=Q(1.0, "bar"), z=z)

    for backend, result in (("python", from_python), ("rust", from_rust)):
        assert result.verdict is StabilityVerdict.STABLE, (
            f"{backend}: the wrong root calls this unstable: {result.tm}"
        )
        assert result.tm[0] == 0.0, f"{backend}: the vapour-like trial should land on the feed"
        assert result.w[0] == (1.0,), f"{backend}: the feed's own stationary point is the feed"
        assert result.iterations[0] == 1, f"{backend}: it starts at the answer"
        assert result.tm[1] > 0.0, (
            f"{backend}: the liquid-like trial is measured from the vapour root"
        )


def test_a_pure_components_trials_are_both_the_feed_and_one_is_trivial() -> None:
    """With one component both seeds normalise to ``[1.0]``, so both trials sit there.

    The trial claiming the root the feed is actually on has
    ``ln W = ln phi(z) - ln phi(z) = 0`` at its first step: ``tm = 0`` exactly,
    ``w == z``, and it converges in one step. Which trial that is is the *feed's*
    decision, not the seed's - a subcooled liquid sits on the liquid root, so it is
    the second trial that is trivial and the first that reports a positive distance;
    a superheated vapour is the other way round; and above the critical temperature
    there is one admissible root and both are trivial.
    """
    for substance in (METHANE, BUTANE, PROPANE):
        fluid = mixture([substance], kij={})
        for t_c in (200.0, 280.0, 330.0):
            for p_pa in (1e5, 1e6, 5e6):
                result = stability_test(fluid, T=Q(t_c, "K"), P=Q(p_pa, "Pa"), z=[1.0])
                context = f"Tc={substance.Tc}, T={t_c}, P={p_pa}"

                assert result.w == ((1.0,), (1.0,)), f"{context}: both trials are the feed"
                trivial = [i for i, distance in enumerate(result.tm) if distance == 0.0]
                assert trivial, f"{context}: the feed's own root has nothing to iterate"
                for i in trivial:
                    assert result.iterations[i] == 1, f"{context}: trial {i} starts there"
                assert all(distance >= 0.0 for distance in result.tm), (
                    f"{context}: a pure substance would be reported unstable: {result.tm}"
                )
                assert result.verdict is StabilityVerdict.STABLE, f"{context}"

    # And which trial is trivial, on two states that differ only in side: butane's
    # saturation pressure at 330 K is about 6 bar.
    butane = mixture([BUTANE], kij={})
    vapour = stability_test(butane, T=Q(330.0, "K"), P=Q(1e5, "Pa"), z=[1.0])
    assert vapour.tm[0] == 0.0, "a superheated vapour is on the vapour root"
    assert vapour.tm[1] > 0.0

    liquid = stability_test(butane, T=Q(330.0, "K"), P=Q(1e6, "Pa"), z=[1.0])
    assert liquid.tm[1] == 0.0, "a subcooled liquid is on the liquid root"
    assert liquid.tm[0] > 0.0


def test_every_reported_trial_is_a_stationary_point() -> None:
    """Each reported trial is a fixed point of the iteration, and ``tm`` is its value.

    Everything below is recomputed from the public API - the feed's lower-Gibbs root
    included, by the rule the module documents - and then the *returned* ``w`` is fed
    back through the model's own stationarity equation:

    .. code-block:: text

        W_i = exp(d_i - ln phi_i(w))     must reproduce w after normalisation
        tm  = 1 - sum_i W_i              must reproduce the returned distance

    A model that reported a partially-converged iterate, a composition from the wrong
    root, or a distance measured somewhere other than at the ``w`` it returned fails
    this without any expected value to compare against.
    """
    for case in CASES:
        fluid = from_names(list(case["inputs"]["components"]))
        t_c, p_pa = case["inputs"]["T"], case["inputs"]["P"]
        z = list(case["inputs"]["z"])
        result = call(case)

        reduced = reduced_parameters(fluid, t_c, p_pa)
        a_mix, b_mix = mixture_parameters(reduced.a, reduced.b, fluid.kij, z)
        roots = pr_z_factor(a_mix, b_mix)

        # Bound now rather than closed over: a closure defined in a loop and reading
        # the loop's variables is a function whose body depends on when it is called,
        # and `partial` is the same binding made explicit.
        gibbs = partial(gibbs_at_root, reduced, fluid.kij, z)
        candidates = [roots.z_min] if roots.z_min == roots.z_max else [roots.z_min, roots.z_max]
        root = min(candidates, key=gibbs)
        feed = phase_state_at(reduced, fluid.kij, z, root)
        d = [math.log(zi) + lp for zi, lp in zip(z, feed.ln_phi, strict=True)]

        # Trial 0 claims the vapour root and trial 1 the liquid one - the ordering
        # contract - so the fixed-point equation is evaluated on each trial's own
        # side. The root is the one at the *trial's* composition, not the feed's:
        # `phase_state` re-derives it there, exactly as the iteration does.
        for trial, liquid in ((0, False), (1, True)):
            w = list(result.w[trial])
            state = phase_state(reduced, fluid.kij, w, liquid=liquid)
            mole_numbers = [math.exp(di - lp) for di, lp in zip(d, state.ln_phi, strict=True)]
            totals = sum(mole_numbers)
            back = [value / totals for value in mole_numbers]
            for i, (got, want) in enumerate(zip(back, w, strict=True)):
                assert abs(got - want) <= 1e-9, (
                    f"{case['id']}: w[{trial}][{i}] is {want}, but the iteration at that "
                    f"composition gives {got} - it is not a stationary point"
                )
            # Absolute rather than relative: at a trivial trial both sides are of
            # order 1e-16, and a relative comparison between two numbers that size is
            # a comparison of round-off. Every meaningful distance is O(1).
            assert abs((1.0 - totals) - result.tm[trial]) <= 1e-9, (
                f"{case['id']}: tm[{trial}] is reported as {result.tm[trial]}, but the "
                f"iteration at the returned w gives {1.0 - totals}"
            )


def test_a_trial_that_cannot_converge_raises_rather_than_calling_the_feed_stable() -> None:
    """The rule the spec's port notes record as a decision.

    A tangent-plane distance bounds stability only *at* a stationary point, so a
    ``tm`` read from a partial iteration is a fact about the iteration rather than
    about the mixture - and discarding it turns "could not tell" into "stable".
    NeqSim breaks out of its loop and carries on; this raises, as the flash does.

    The state is a feed with an **absent** component, whose reference potential is
    ``-inf`` and whose trial mole number is zero and stays zero, so the rms change in
    ``ln W`` is ``-inf - -inf`` - a ``NaN`` - at every step. The iteration never meets
    a tolerance and runs to the cap: 2000 steps, by design, and an error rather than a
    verdict. The cap is deliberately far above every state measured (the worst on the
    binaries is 507), which is why it takes a degenerate composition to reach it.
    """
    with pytest.raises(SolverNotConvergedError):
        stability_test(methane_butane(), T=Q(330.0, "K"), P=Q(2.5e6, "Pa"), z=[1.0, 0.0])


def test_the_trials_are_reported_vapour_like_first() -> None:
    """``tm``, ``w`` and ``iterations`` are positional, and the order is the contract.

    The spec's first case separates the two trials: the vapour-like one converges to
    the feed (16 steps, ``tm = 0``) and the liquid-like one proves the feed unstable
    (9 steps, ``tm = -0.215``). The rows are the two Wilson seeds' stationary points
    rather than two guesses - at this state methane's K-value is well above one and
    butane's well below, so ``z_i K_i`` is the methane-rich seed.
    """
    result = stability_test(methane_butane(), T=Q(330.0, "K"), P=Q(2.5e6, "Pa"), z=[0.6, 0.4])
    assert result.tm[0] > result.tm[1], (
        f"the vapour-like trial is the one that says nothing: {result.tm}"
    )
    assert result.iterations == (16, 9)
    assert result.w[0][0] > result.w[1][0], "row 0 should be the vapour-like trial"
    assert result.w[1][1] > result.w[0][1], "row 1 should be the liquid-like trial"


def test_every_verdict_is_reachable() -> None:
    """Both members of ``StabilityVerdict`` are reached by some state.

    A value a caller branches on and never sees is worse than an absent one, which is
    why ``RootStructure``'s unreachable variant was removed rather than kept for
    symmetry. This is the same check for this model's enum - and unlike the flash's
    ``Phase``, both members here are easy to reach, so the check is cheap.
    """
    fluids = (methane_butane(), mixture([PROPANE, BUTANE], kij={(0, 1): 0.02}))
    seen = set()
    for fluid in fluids:
        for t_c in (250.0, 280.0, 300.0, 330.0, 350.0, 400.0):
            for p_pa in (1e5, 1e6, 3e6, 1e7, 2e7, 5e7):
                for z in ([0.1, 0.9], [0.5, 0.5], [0.9, 0.1]):
                    try:
                        seen.add(stability_test(fluid, T=Q(t_c, "K"), P=Q(p_pa, "Pa"), z=z).verdict)
                    except SolverNotConvergedError:
                        continue
    assert seen == set(StabilityVerdict), f"unreachable verdict(s): {set(StabilityVerdict) - seen}"


def test_no_distance_is_ever_nan() -> None:
    """A ``NaN`` distance is the one value a caller cannot reason about.

    It satisfies neither ``tm < -1e-8`` nor its negation, so a feed would be called
    ``stable`` by a comparison that never ran - the right answer for the wrong
    reason, and silently. This asserts the distances are numbers wherever a result is
    returned at all.
    """
    fluid = methane_butane()
    checked = 0
    for t_c in (250.0, 280.0, 300.0, 330.0, 350.0, 400.0):
        for p_pa in (1e5, 1e6, 3e6, 1e7, 5e7):
            for z in ([0.1, 0.9], [0.5, 0.5], [0.9, 0.1]):
                try:
                    result = stability_test(fluid, T=Q(t_c, "K"), P=Q(p_pa, "Pa"), z=z)
                except SolverNotConvergedError:
                    continue
                checked += 1
                for i, distance in enumerate(result.tm):
                    assert math.isfinite(distance), f"T={t_c}, P={p_pa}, z={z}: tm[{i}]"
                assert math.isfinite(result.min_t_over_tc) and result.min_t_over_tc > 0.0
    assert checked > 80, f"the sweep only covered {checked} states"


def test_the_verdict_is_the_threshold_applied_to_the_distances() -> None:
    """``unstable`` if and only if some ``tm`` is below ``-1e-8``, and nothing else.

    The identity that ties the two reported things together, asserted on returned
    results rather than recomputed from the inputs, so a result reporting a verdict
    its own distances do not support fails.
    """
    fluid = methane_butane()
    unstable = 0
    for t_c in (250.0, 280.0, 300.0, 330.0, 350.0, 400.0):
        for p_pa in (1e5, 1e6, 3e6, 1e7, 5e7):
            for z in ([0.1, 0.9], [0.5, 0.5], [0.9, 0.1]):
                try:
                    result = stability_test(fluid, T=Q(t_c, "K"), P=Q(p_pa, "Pa"), z=z)
                except SolverNotConvergedError:
                    continue
                expected = StabilityVerdict.STABLE
                if any(distance < -1e-8 for distance in result.tm):
                    expected = StabilityVerdict.UNSTABLE
                    unstable += 1
                assert result.verdict is expected, f"T={t_c}, P={p_pa}, z={z}"
    assert unstable > 0, "the sweep found no unstable feed at all"


def test_a_feed_that_does_not_sum_to_one_is_refused() -> None:
    fluid = methane_butane()
    for z in ([0.6, 0.5], [0.6, 0.3], [0.6, 0.4, 0.0]):
        with pytest.raises(InvalidInputError):
            stability_test(fluid, T=Q(330.0, "K"), P=Q(2.5e6, "Pa"), z=z)


def test_a_negative_mole_fraction_is_refused() -> None:
    with pytest.raises(InvalidInputError):
        stability_test(methane_butane(), T=Q(330.0, "K"), P=Q(2.5e6, "Pa"), z=[1.6, -0.6])


def test_a_non_positive_temperature_or_pressure_is_refused() -> None:
    fluid = methane_butane()
    for t_c, p_pa in ((0.0, 2.5e6), (-1.0, 2.5e6), (330.0, 0.0)):
        with pytest.raises(OutOfRangeError):
            stability_test(fluid, T=Q(t_c, "K"), P=Q(p_pa, "Pa"), z=[0.6, 0.4])


def test_a_malformed_kij_matrix_is_refused() -> None:
    """All three malformations are invisible in the answer, so all three are refused."""
    components = (METHANE, BUTANE)
    with pytest.raises(InvalidInputError):
        Mixture(components=components, kij=((0.0, 0.05), (0.05, 0.0), (0.0, 0.0)))
    with pytest.raises(InvalidInputError):
        Mixture(components=components, kij=((0.1, 0.05), (0.05, 0.0)))
    with pytest.raises(InvalidInputError):
        Mixture(components=components, kij=((0.0, 0.05), (0.07, 0.0)))
    Mixture(components=components, kij=((0.0, 0.05), (0.05, 0.0)))


def test_a_near_critical_component_warns_rather_than_failing() -> None:
    """The bound warns: this model rests on the trials finding *distinct* points."""
    result = stability_test(methane_butane(), T=Q(400.0, "K"), P=Q(2.0e6, "Pa"), z=[0.6, 0.4])
    assert any(w.code == "OUT_OF_VALID_RANGE" for w in result.warnings)
    h.assert_close(result.min_t_over_tc, 400.0 / 425.12, 1e-12, "min_t_over_tc")


def test_the_two_backends_agree_on_every_spec_case() -> None:
    """The cross-language claim, case by case, including the iteration counts.

    ``test_cross_impl.py`` does this for the calc registry; the models are a separate
    registry with a separate bridge table, so they need their own. The iteration
    counts are the sharpest of the layers: a ULP difference that flips one iteration
    shows up here as a number rather than as a mysterious drift.

    ``w`` is compared here for a reason beyond symmetry - it is the field the generated
    Rust table cannot carry, so this is what pins the Rust side's compositions to the
    spec's.
    """
    for case in CASES:
        with use_backend("python"):
            py = call(case)
        with use_backend("rust"):
            rs = call(case)

        assert py.verdict is rs.verdict, f"{case['id']}: verdict"
        assert py.iterations == rs.iterations, f"{case['id']}: iteration counts"
        h.assert_close(py.min_t_over_tc, rs.min_t_over_tc, 1e-12, f"{case['id']}")
        for i, (left, right) in enumerate(zip(py.tm, rs.tm, strict=True)):
            assert_tm(right, left, 1e-12, f"{case['id']} (tm[{i}], rust against python)")
        for i, (left_row, right_row) in enumerate(zip(py.w, rs.w, strict=True)):
            for j, (left, right) in enumerate(zip(left_row, right_row, strict=True)):
                h.assert_close(right, left, 1e-12, f"{case['id']} (w[{i}][{j}])")
        assert [w.code for w in py.warnings] == [w.code for w in rs.warnings]


def test_the_verdict_reaches_both_backends_as_the_same_enum() -> None:
    """The enum must not degrade to a string on one path.

    ``verdict`` crosses the FFI boundary as the spec's spelling and is rebuilt here,
    so ``result.verdict is StabilityVerdict.UNSTABLE`` has to hold whichever backend
    answered. The flash's ``phase`` makes the same promise, and a caller comparing
    against a literal string on one backend and an enum on the other is exactly the
    bug that promise prevents.
    """
    fluid = methane_butane()
    for t_c, p_pa, expected in (
        (330.0, 2.5e6, StabilityVerdict.UNSTABLE),
        (430.0, 6.0e6, StabilityVerdict.STABLE),
    ):
        with use_backend("python"):
            from_python = stability_test(fluid, T=Q(t_c, "K"), P=Q(p_pa, "Pa"), z=[0.6, 0.4])
        with use_backend("rust"):
            from_rust = stability_test(fluid, T=Q(t_c, "K"), P=Q(p_pa, "Pa"), z=[0.6, 0.4])

        for backend, result in (("python", from_python), ("rust", from_rust)):
            assert result.verdict is expected, f"{backend} at T={t_c}, P={p_pa}"
            assert isinstance(result.verdict, StabilityVerdict), (
                f"{backend}: the verdict arrived as {type(result.verdict).__name__}"
            )


def test_the_cases_are_announced_in_the_model_docs() -> None:
    """Every model has a generated page, and the book names it."""
    from pathlib import Path

    root = Path(__file__).resolve().parents[3]
    page = root / "docs" / "src" / "eos" / "stability_test.md"
    assert page.is_file(), f"{page} was not generated"
    summary = (root / "docs" / "src" / "SUMMARY.md").read_text(encoding="utf-8")
    assert "eos/stability_test.md" in summary, "the model is not in the book's contents"
