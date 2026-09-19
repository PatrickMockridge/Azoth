"""Spec-driven tests for the ``eos.pt_flash`` model.

The spec's cases pin the two implementations to each other. These tests are for the
things a case *cannot* say: the identities that hold at any answer, the reductions
that tie the model layer to the registered kernels below it, and the failure modes
that were found by measurement rather than by construction - including the two that
an earlier draft of the model got wrong.
"""

from __future__ import annotations

import math
from typing import Any

import pytest

import _helpers as h
from azoth import _models_gen, ureg, use_backend
from azoth.core.errors import InvalidInputError, OutOfRangeError, SolverNotConvergedError
from azoth.core.result import Phase, PtFlashResult
from azoth.eos import (
    Mixture,
    component,
    from_names,
    mixture,
    pr_alpha_ab,
    pr_departure,
    pr_kappa,
    pr_z_factor,
    pt_flash,
    rachford_rice_binary,
    vdw1f_mix_binary,
)

MODEL_ID = "eos.pt_flash"
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


def call(case: dict[str, Any]) -> PtFlashResult:
    """Run one case declared in the model spec.

    Through :func:`_helpers.model_kwargs`, which is the one place a case's
    declared inputs become arguments: it resolves `components` against the
    databank and hands over the mixture and the ideal-gas model the function
    takes. A hand-built mixture here would be a second fluid, described by the
    case file rather than by NeqSim's tables.
    """
    return pt_flash(**h.model_kwargs(SPEC, case["inputs"]))


@pytest.mark.parametrize("case", CASES, ids=lambda c: c["id"])
def test_spec_case(case: dict[str, Any]) -> None:
    """Run one case declared in the model spec."""
    result = call(case)
    tolerance = case["tolerance"]
    expected = case["expected"]

    for name in ("z_liquid", "z_vapour", "min_t_over_tc"):
        h.assert_close(getattr(result, name), expected[name], tolerance, f"{case['id']} ({name})")

    # `residual` is deliberately *not* pinned as a case value: it is an rms of
    # `ln`-built quantities, so the two implementations differ in its last few ulps
    # and a value near 4.7e-11 has no significant figures to spare. What is
    # reproducible - and what the case does assert - is the iteration count.
    # See the spec's correction 4.
    assert result.residual <= SPEC["algorithm"]["tolerance"], (
        f"{case['id']}: residual {result.residual:e} exceeds the declared tolerance"
    )
    assert result.iterations == expected["iterations"], f"{case['id']}: iteration count"

    assert result.beta is not None, f"{case['id']}: every spec case is a split"
    h.assert_close(result.beta, expected["beta"], tolerance, f"{case['id']} (beta)")

    for name in ("x", "y", "k", "ln_phi_liquid", "ln_phi_vapour"):
        actual = getattr(result, name)
        assert len(actual) == len(expected[name]), f"{case['id']}: {name} length"
        for i, (got, want) in enumerate(zip(actual, expected[name], strict=True)):
            h.assert_close(got, want, tolerance, f"{case['id']} ({name}[{i}])")

    h.assert_consistent(result, case["id"])


def test_every_case_ran() -> None:
    """Guard against a spec edit that silently removes every case."""
    assert len(CASES) >= 3, f"expected several cases, found {len(CASES)}"


def test_the_five_identities_hold_at_every_returned_answer() -> None:
    """Material balance, both normalisations, the K-definition, and the RR residual.

    Recomputed here from the result rather than read from it, so a result reporting
    an identity it had not achieved fails. A plausible-looking wrong answer cannot
    satisfy all five.
    """
    for case in CASES:
        result = call(case)
        z = case["inputs"]["z"]
        beta = result.beta
        assert beta is not None

        for i, zi in enumerate(z):
            material = (1.0 - beta) * result.x[i] + beta * result.y[i]
            assert abs(material - zi) <= 1e-12, f"{case['id']}: material balance at {i}"
            assert abs(result.k[i] - result.y[i] / result.x[i]) <= 1e-12, (
                f"{case['id']}: K[{i}] is not y/x"
            )

        assert abs(sum(result.x) - 1.0) <= 1e-12, f"{case['id']}: sum x"
        assert abs(sum(result.y) - 1.0) <= 1e-12, f"{case['id']}: sum y"

        residual = sum(
            zi * (ki - 1.0) / (1.0 + beta * (ki - 1.0)) for zi, ki in zip(z, result.k, strict=True)
        )
        assert abs(residual) <= 1e-12, f"{case['id']}: Rachford-Rice residual {residual:e}"


def test_the_answer_is_a_gibbs_minimum() -> None:
    """The strongest data-free correctness test: the answer minimises `G / RT`.

    Perturbing `beta` at fixed `K` must raise the total Gibbs energy on both sides,
    and raise it *quadratically* - which a converged-but-wrong answer cannot do.
    """
    fluid = methane_butane()
    z = [0.1, 0.9]
    result = pt_flash(fluid, T=Q(300.0, "K"), P=Q(3.0e6, "Pa"), z=z)
    beta = result.beta
    assert beta is not None

    def g(beta: float) -> float:
        """`G / RT` at a trial vapour fraction, with the K-values held fixed."""
        x = [zi / (1.0 + beta * (ki - 1.0)) for zi, ki in zip(z, result.k, strict=True)]
        y = [ki * xi for ki, xi in zip(result.k, x, strict=True)]
        liquid = sum(
            xi * (math.log(xi) + lp) for xi, lp in zip(x, result.ln_phi_liquid, strict=True)
        )
        vapour = sum(
            yi * (math.log(yi) + lp) for yi, lp in zip(y, result.ln_phi_vapour, strict=True)
        )
        return (1.0 - beta) * liquid + beta * vapour

    base = g(beta)
    for step in (1e-3, 1e-4, 1e-5):
        for offset in (-step, step):
            assert g(beta + offset) - base > 0.0, (
                f"perturbing beta by {offset:e} lowered G/RT - the answer is not a minimum"
            )
    coarse, fine = g(beta + 1e-3) - base, g(beta + 1e-4) - base
    assert abs(coarse / fine - 100.0) < 1.0, (
        f"the rise is not quadratic: {coarse:e} at 1e-3 against {fine:e} at 1e-4"
    )


def test_the_mixture_form_reduces_to_pr_departure_at_one_component() -> None:
    """`N = 1` must reproduce the registered pure-component kernel exactly.

    The mixture fugacity coefficient is the one piece of arithmetic in this layer
    that no registered calculation covers. This is what keeps it honest: at `N = 1`
    the cross-sum factor collapses to 1 and the whole expression becomes
    `eos.pr_departure`'s, which has its own spec, worked example and source.
    """
    from azoth.eos.reference._mixture_state import ReducedParameters
    from azoth.eos.reference._mixture_state import phase_state as _phase_state

    for substance, tr, pr in ((PROPANE, 0.8, 0.25), (BUTANE, 0.7, 0.4), (METHANE, 1.2, 0.9)):
        tc = substance.Tc.to_base_units().magnitude
        pc = substance.Pc.to_base_units().magnitude
        kappa = pr_kappa(substance.omega).kappa
        ab = pr_alpha_ab(kappa, tr, pr)
        sqrt_tr = tr**0.5
        reduced = ReducedParameters(
            a=[ab.a_reduced],
            b=[ab.b_reduced],
            psi=[-kappa * sqrt_tr / (1.0 + kappa * (1.0 - sqrt_tr))],
            # `T*dpsi/dT` at one substance, which is what the mixture's departure heat
            # capacity reduces to. Written out rather than read back from the call under
            # test, so the reduction is checked rather than asserted.
            psi_t=[
                -kappa * (1.0 + kappa) * tr / (2.0 * sqrt_tr * (1.0 + kappa * (1.0 - sqrt_tr)) ** 2)
            ],
            warnings=[],
        )
        state = _phase_state(reduced, ((0.0,),), [1.0], liquid=False)
        z, ln_phi = state.z, state.ln_phi
        kappa = pr_kappa(substance.omega).kappa
        pure = pr_departure(ab.a_reduced, ab.b_reduced, z, kappa, tr)
        h.assert_close(ln_phi[0], pure.ln_phi, 1e-12, f"ln phi for Tc={tc}, Pc={pc}")


def test_the_mixture_parameters_and_vapour_fraction_reduce_to_the_binary_kernels() -> None:
    """`N = 2` must reproduce `eos.vdw1f_mix_binary` and `eos.rachford_rice_binary`.

    A *cross-layer* check no single-language test can replace: the binary kernels
    are registered calculations with their own specs, and the general N-component
    code path has to agree with them where both apply. The mixture parameters match
    to a couple of ulps rather than bit-for-bit, because the registered kernel
    evaluates its three terms longhand while this sums a double loop.
    """
    from azoth.eos.reference._mixture_state import (
        mixture_parameters as _mixture_parameters,
    )
    from azoth.eos.reference._mixture_state import rachford_rice, rachford_rice_bounds

    fluid = methane_butane()
    for t_c, p_pa in ((330.0, 2.5e6), (300.0, 3.0e6), (350.0, 5.0e6)):
        a, b = [], []
        for substance in fluid.components:
            ab = pr_alpha_ab(
                pr_kappa(substance.omega).kappa,
                t_c / substance.Tc.to_base_units().magnitude,
                p_pa / substance.Pc.to_base_units().magnitude,
            )
            a.append(ab.a_reduced)
            b.append(ab.b_reduced)

        for x1 in (0.1, 0.3, 0.6, 0.9):
            x = [x1, 1.0 - x1]
            # The kernel takes the pair's interaction parameter as its last argument;
            # it is read off the fluid rather than written here, so the two terms of
            # the comparison cannot describe different mixtures.
            a_mix, b_mix = _mixture_parameters(a, b, fluid.kij, x)
            mixture_kernel = vdw1f_mix_binary(x1, a[0], a[1], b[0], b[1], fluid.kij[0][1])
            h.assert_close(a_mix, mixture_kernel.a_mix, 1e-12, f"a_mix at x1={x1}")
            h.assert_close(b_mix, mixture_kernel.b_mix, 1e-12, f"b_mix at x1={x1}")

    for t_c, p_pa, z in ((330.0, 2.5e6, [0.6, 0.4]), (300.0, 3.0e6, [0.1, 0.9])):
        result = pt_flash(fluid, T=Q(t_c, "K"), P=Q(p_pa, "Pa"), z=z)
        assert result.beta is not None
        bounds = rachford_rice_bounds(list(result.k))
        assert bounds is not None
        beta = rachford_rice(list(z), list(result.k), bounds, 1e-14, 200)
        closed_form = rachford_rice_binary(z[0], result.k[0], result.k[1])
        h.assert_close(beta, closed_form.beta, 1e-12, f"beta at T={t_c}")

        # And the same reduction through the public API, which is the one a caller
        # would actually be relying on.
        h.assert_close(result.beta, closed_form.beta, 1e-12, f"reported beta at T={t_c}")


def test_a_feed_with_no_rachford_rice_root_is_single_phase() -> None:
    """Every K-value on one side of one means the RR equation has no root.

    Not a numerical failure but a proof: `sum_i y_i = sum_i K_i x_i = 1` with
    `sum_i x_i = 1` is impossible when every `K_i < 1`, and equally when every
    `K_i > 1`. `beta` is therefore absent, and the loop settles on the **first**
    iteration - a diagnosis, not an iteration that failed to converge.
    """
    fluid = methane_butane()
    for t_c, p_pa, expected in (
        (280.0, 100_000.0, Phase.ALL_VAPOUR),
        (330.0, 100_000.0, Phase.ALL_VAPOUR),
        (280.0, 50_000_000.0, Phase.ALL_LIQUID),
        (330.0, 50_000_000.0, Phase.ALL_LIQUID),
    ):
        z = [0.6, 0.4]
        result = pt_flash(fluid, T=Q(t_c, "K"), P=Q(p_pa, "Pa"), z=z)
        assert result.phase is expected, f"T={t_c}, P={p_pa}"
        assert result.beta is None, f"T={t_c}: no root means no vapour fraction"
        assert result.iterations == 1, f"T={t_c}: a no-root diagnosis settles at once"
        assert _is_nan(result.residual), f"T={t_c}: no step completed"
        assert any(w.code == "TRIVIAL_SOLUTION" for w in result.warnings), f"T={t_c}"
        assert list(result.x) == z and list(result.y) == z
        assert all(k != 1.0 for k in result.k)


def _is_nan(value: float) -> bool:
    return value != value


def test_a_negative_flash_reports_its_vapour_fraction() -> None:
    """A converged `beta` outside `[0, 1]` is a reading, not a failure.

    The negative flash is the amount of the absent phase that would have to be
    added to bring the feed to saturation. It is a real solution of the same
    equations - and the material balance holds on it exactly as it does on a split -
    so it is reported with `OUT_OF_VALID_RANGE` rather than suppressed.
    """
    fluid = methane_butane()
    z = [0.6, 0.4]
    result = pt_flash(fluid, T=Q(330.0, "K"), P=Q(1.0e6, "Pa"), z=z)

    assert result.phase is Phase.ALL_VAPOUR
    assert result.beta is not None and result.beta > 1.0
    assert result.iterations == 9, "the iteration ran, so this is not a proof"
    assert any(w.code == "OUT_OF_VALID_RANGE" for w in result.warnings)

    for i, zi in enumerate(z):
        material = (1.0 - result.beta) * result.x[i] + result.beta * result.y[i]
        assert abs(material - zi) <= 1e-12, f"the negative flash breaks the balance at {i}"
        assert result.x[i] > 0.0 and result.y[i] > 0.0


def test_the_phase_label_and_the_vapour_fraction_agree() -> None:
    """`two_phase` means a split, and a split means `two_phase`.

    The spec's `phase` description is a contract about the two fields together: `beta`
    is in `[0, 1]` when the phase is `two_phase`, absent when it is `trivial`, and an
    extrapolation - reported, with a warning - otherwise. Each half is checked
    elsewhere; the implication between them was not, and a model reporting
    `two_phase` with a vapour fraction of 1.7 would satisfy every single-field test.

    Swept rather than spot-checked, and the sweep has to reach every phase value or
    the assertion above it proves nothing for the phases it missed.
    """
    fluid = methane_butane()
    z = [0.6, 0.4]
    seen: set[Phase] = set()

    for temperature in range(240, 460, 10):
        for pressure in (1.0e5, 1.0e6, 2.0e6, 5.0e6, 2.0e7):
            try:
                result = pt_flash(fluid, T=Q(float(temperature), "K"), P=Q(pressure, "Pa"), z=z)
            except (OutOfRangeError, SolverNotConvergedError):
                # A state this model does not cover; the point of the sweep is the
                # states it does.
                continue

            where = f"T={temperature} K, P={pressure} Pa"
            seen.add(result.phase)
            if result.phase is Phase.TWO_PHASE:
                assert result.beta is not None, f"{where}: a split with no vapour fraction"
                assert 0.0 <= result.beta <= 1.0, (
                    f"{where}: `two_phase` with beta = {result.beta}, which is not a "
                    f"split and not a number a caller can use as one"
                )
            elif result.phase is Phase.TRIVIAL:
                assert result.beta is None, f"{where}: a trivial solution with a beta"

    assert seen == set(Phase), f"the sweep reached only {seen}, so it proved little"


def test_a_feed_that_converges_to_the_trivial_solution_says_so() -> None:
    """The model's honest gap, made testable.

    Where the K-values straddle one throughout, the no-root proof never applies and
    a single-phase feed converges to `x = y = z`. There every `K_i` is 1, `g(beta)`
    is identically zero, and the vapour fraction is *indeterminate* rather than out
    of range - measured at -7.7e10 for one feed and -2.2e11 for another, neither
    reproducible. So the contract is that no number is reported and the caller is
    told why, not that a plausible-looking number is.
    """
    fluid = methane_butane()
    z = [0.6, 0.4]
    for t_c in (280.0, 300.0, 330.0):
        result = pt_flash(fluid, T=Q(t_c, "K"), P=Q(20.0e6, "Pa"), z=z)
        assert result.phase is Phase.TRIVIAL, f"T={t_c}"
        assert result.beta is None, f"T={t_c}: a trivial solution has no vapour fraction"
        assert result.iterations > 1, f"T={t_c}: this is an iteration, not a proof"
        assert any(w.code == "TRIVIAL_SOLUTION" for w in result.warnings), f"T={t_c}"
        # Stated exactly, rather than as the last iterate that approached it.
        assert list(result.x) == z and list(result.y) == z
        assert all(k == 1.0 for k in result.k), f"T={t_c}"


def test_every_phase_value_is_reachable() -> None:
    """Every `Phase` member is reached by some state, so none is a dead branch.

    `RootStructure` had a `TWO_ROOTS` variant that could not occur, and it was
    removed rather than kept for symmetry. `ALL_LIQUID` in particular was found by
    sweeping states, not assumed - a caller branching on a value that never arrives
    is worse off than one branching on an absent value.
    """
    fluids = (
        methane_butane(),
        mixture([PROPANE, BUTANE], kij={(0, 1): 0.02}),
    )
    seen = set()
    for fluid in fluids:
        for t_c in (250.0, 280.0, 300.0, 330.0, 350.0, 400.0):
            for p_pa in (1e5, 1e6, 3e6, 1e7, 2e7, 5e7):
                for z in ([0.1, 0.9], [0.5, 0.5], [0.9, 0.1]):
                    try:
                        seen.add(pt_flash(fluid, T=Q(t_c, "K"), P=Q(p_pa, "Pa"), z=z).phase)
                    except (SolverNotConvergedError, InvalidInputError):
                        continue
    assert seen == set(Phase), f"unreachable phase(s): {set(Phase) - seen}"


def test_both_compositions_stay_positive_at_every_returned_answer() -> None:
    """The property the correct Rachford-Rice bracket is *equivalent* to.

    The bracket is exactly the interval on which `1 + beta (K_i - 1) > 0` for every
    `i`, which is what keeps both compositions in `[0, 1]`. So sweeping states and
    asserting that is asserting the bracket - and it is the check that caught the
    second of the two wrong brackets, which returned a vapour mole fraction of 1.162
    for a feed with no root at all.
    """
    fluids = (
        methane_butane(),
        mixture([PROPANE, BUTANE], kij={(0, 1): 0.02}),
    )
    splits = 0
    for fluid in fluids:
        for t_c in (250.0, 280.0, 300.0, 330.0, 350.0):
            for p_pa in (1e5, 1e6, 3e6, 1e7):
                for z in ([0.1, 0.9], [0.5, 0.5], [0.9, 0.1]):
                    try:
                        result = pt_flash(fluid, T=Q(t_c, "K"), P=Q(p_pa, "Pa"), z=z)
                    except (SolverNotConvergedError, InvalidInputError):
                        continue
                    for name, values in (("x", result.x), ("y", result.y)):
                        assert all(0.0 < v <= 1.0 for v in values), (
                            f"T={t_c}, P={p_pa}, z={z}: {name} = {values} is not a composition"
                        )
                    if result.phase is Phase.TWO_PHASE:
                        splits += 1
    assert splits > 20, f"the sweep should find plenty of two-phase states, found {splits}"


def test_a_feed_that_does_not_sum_to_one_is_refused() -> None:
    fluid = methane_butane()
    for z in ([0.6, 0.5], [0.6, 0.3], [0.6, 0.4, 0.0]):
        with pytest.raises(InvalidInputError):
            pt_flash(fluid, T=Q(330.0, "K"), P=Q(2.5e6, "Pa"), z=z)


def test_a_negative_mole_fraction_is_refused() -> None:
    with pytest.raises(InvalidInputError):
        pt_flash(methane_butane(), T=Q(330.0, "K"), P=Q(2.5e6, "Pa"), z=[1.6, -0.6])


def test_a_non_positive_temperature_or_pressure_is_refused() -> None:
    fluid = methane_butane()
    for t_c, p_pa in ((0.0, 2.5e6), (-1.0, 2.5e6), (330.0, 0.0)):
        with pytest.raises(OutOfRangeError):
            pt_flash(fluid, T=Q(t_c, "K"), P=Q(p_pa, "Pa"), z=[0.6, 0.4])


def test_a_malformed_kij_matrix_is_refused() -> None:
    """All three malformations are invisible in the answer, so all three are refused.

    A non-zero diagonal silently rescales that component's attraction; an
    asymmetric matrix has no meaning in a double sum that is symmetric by
    construction; the wrong shape is a different mixture entirely.
    """
    components = (METHANE, BUTANE)
    with pytest.raises(InvalidInputError):
        Mixture(components=components, kij=((0.0, 0.05), (0.05, 0.0), (0.0, 0.0)))
    with pytest.raises(InvalidInputError):
        Mixture(components=components, kij=((0.1, 0.05), (0.05, 0.0)))
    with pytest.raises(InvalidInputError):
        Mixture(components=components, kij=((0.0, 0.05), (0.07, 0.0)))
    Mixture(components=components, kij=((0.0, 0.05), (0.05, 0.0)))


def test_the_kij_helper_rejects_an_out_of_range_pair() -> None:
    with pytest.raises(InvalidInputError):
        mixture([METHANE, BUTANE], kij={(0, 5): 0.1})
    with pytest.raises(InvalidInputError):
        mixture([METHANE, BUTANE], kij={(0, 0): 0.1})


def test_a_near_critical_component_warns_rather_than_failing() -> None:
    """The bound warns: the arithmetic is still defined, it is the accuracy that goes."""
    result = pt_flash(methane_butane(), T=Q(400.0, "K"), P=Q(2.0e6, "Pa"), z=[0.6, 0.4])
    assert any(w.code == "OUT_OF_VALID_RANGE" for w in result.warnings)
    h.assert_close(result.min_t_over_tc, 400.0 / 425.12, 1e-12, "min_t_over_tc")


def test_the_result_is_clean_at_an_ordinary_state() -> None:
    result = pt_flash(methane_butane(), T=Q(330.0, "K"), P=Q(2.5e6, "Pa"), z=[0.6, 0.4])
    assert result.is_clean, f"unexpected warnings: {result.warnings}"


def test_the_reported_roots_are_the_cubics_roots_at_the_reported_compositions() -> None:
    """Ties the model to `eos.pr_z_factor` rather than to its own arithmetic.

    A `z_liquid` that was not the smallest admissible root would make every `ln phi`
    downstream wrong while `beta` stayed entirely plausible.
    """
    from azoth.eos.reference._mixture_state import mixture_parameters as _mixture_parameters

    fluid = methane_butane()
    for t_c, p_pa, z in ((330.0, 2.5e6, [0.6, 0.4]), (300.0, 3.0e6, [0.1, 0.9])):
        result = pt_flash(fluid, T=Q(t_c, "K"), P=Q(p_pa, "Pa"), z=z)
        a, b = [], []
        for substance in fluid.components:
            ab = pr_alpha_ab(
                pr_kappa(substance.omega).kappa,
                t_c / substance.Tc.to_base_units().magnitude,
                p_pa / substance.Pc.to_base_units().magnitude,
            )
            a.append(ab.a_reduced)
            b.append(ab.b_reduced)

        for composition, reported in ((result.x, result.z_liquid), (result.y, result.z_vapour)):
            a_mix, b_mix = _mixture_parameters(a, b, fluid.kij, list(composition))
            roots = pr_z_factor(a_mix, b_mix)
            assert (roots.z_min == reported) or (roots.z_max == reported), (
                f"reported root {reported} is not one of the cubic's"
            )


def test_the_two_backends_agree_on_every_spec_case() -> None:
    """The cross-language claim, case by case, including the iteration count.

    `test_cross_impl.py` does this for the calc registry; the models are a separate
    registry with a separate bridge table, so they need their own. The iteration
    count is the sharpest of the four layers of claim: a ULP difference that flips
    one iteration shows up here as a number rather than as a mysterious drift.
    """
    for case in CASES:
        with use_backend("python"):
            py = call(case)
        with use_backend("rust"):
            rs = call(case)

        assert py.iterations == rs.iterations, f"{case['id']}: iteration count"
        assert py.phase is rs.phase, f"{case['id']}: phase"
        assert py.beta is not None and rs.beta is not None, f"{case['id']}: a split"
        h.assert_close(py.beta, rs.beta, 1e-12, f"{case['id']} (beta)")
        h.assert_close(py.z_liquid, rs.z_liquid, 1e-12, f"{case['id']} (z_liquid)")
        h.assert_close(py.z_vapour, rs.z_vapour, 1e-12, f"{case['id']} (z_vapour)")
        for name in ("x", "y", "k", "ln_phi_liquid", "ln_phi_vapour"):
            for i, (left, right) in enumerate(
                zip(getattr(py, name), getattr(rs, name), strict=True)
            ):
                h.assert_close(left, right, 1e-12, f"{case['id']} ({name}[{i}])")


def test_the_single_phase_diagnosis_reaches_both_backends() -> None:
    """The outcomes that carry no `beta` must agree across languages too.

    A cross-language test that only covered splits would miss the two cases where
    the answer is an *absence*, and an absence is exactly the kind of thing two
    implementations can disagree about silently - one returning `None` and the other
    a sentinel.

    **The iteration count is not asserted here, and at one of these three states it
    differs** (`330 K, 20 MPa`: 48 steps against 37). The count at a trivial solution
    is where the iteration *crossed* the trivial threshold, and it crosses it on a
    slow crawl - so the two kernels, whose `log` and `sqrt` differ in the last bit,
    cross it a few steps apart while agreeing on the answer. The split test above
    still pins the count, and it can: a split converges to a point rather than
    crawling through a threshold.
    """
    fluid = methane_butane()
    for t_c, p_pa, z in (
        (330.0, 1.0e5, [0.6, 0.4]),
        (330.0, 5.0e7, [0.6, 0.4]),
        (330.0, 2.0e7, [0.6, 0.4]),
    ):
        with use_backend("python"):
            py = pt_flash(fluid, T=Q(t_c, "K"), P=Q(p_pa, "Pa"), z=z)
        with use_backend("rust"):
            rs = pt_flash(fluid, T=Q(t_c, "K"), P=Q(p_pa, "Pa"), z=z)

        assert py.beta is None and rs.beta is None, f"T={t_c}, P={p_pa}"
        assert py.phase is rs.phase, f"T={t_c}, P={p_pa}"
        assert py.k == rs.k, f"T={t_c}, P={p_pa}: the trivial K-values"
        assert [w.code for w in py.warnings] == [w.code for w in rs.warnings]


def test_the_cases_are_announced_in_the_model_docs() -> None:
    """Every model has a generated page, and the book names it."""
    from pathlib import Path

    root = Path(__file__).resolve().parents[3]
    page = root / "docs" / "src" / "eos" / "pt_flash.md"
    assert page.is_file(), f"{page} was not generated"
    summary = (root / "docs" / "src" / "SUMMARY.md").read_text(encoding="utf-8")
    assert "eos/pt_flash.md" in summary, "the model is not in the book's contents"


def test_the_second_order_fallback_converges_where_successive_substitution_crawls() -> None:
    """The handover to NeqSim's second-order scheme, and what it buys.

    NeqSim's ``TPflash`` switches to ``SysNewtonRhapsonTPflash`` once its own scheme has
    run ``newtonLimit`` steps. At methane/n-butane 0.6/0.4, 340 K and 120 bar, successive
    substitution needs 218 steps and the second-order scheme needs 36, and the *answer*
    is the same to 1e-7 - the two stop on different measures, so their last digits
    differ.

    **1e-7 rather than 1e-9 because the state at 340 K and 120 bar is marginal**: the
    handover reaches it in 28 steps rather than 36 and the vapour fraction lands 2.1e-9
    away, where the other three states land within 1e-14. Both are functions of where
    the iteration stopped on a crawl.

    The iteration count is the sabotage detector: reversing the sign of the liquid term
    in the Jacobian takes the count to 244, worse than the scheme it replaced, which is
    what a wrong descent direction does.
    """
    fluid = methane_butane()
    for t_c, p_pa, ss_only, hybrid, beta in (
        (340.0, 1.2e7, 218, 36, 0.2099871524789072),
        (340.0, 1.0e7, 53, 26, 0.4902672003063496),
        (350.0, 1.0e7, 68, 26, 0.5851785631591365),
        (320.0, 1.2e7, 88, 28, 0.13061921281724959),
    ):
        r = pt_flash(fluid, T=Q(t_c, "K"), P=Q(p_pa, "Pa"), z=[0.6, 0.4])
        assert r.phase is Phase.TWO_PHASE, f"{t_c}, {p_pa}"
        assert r.beta is not None, f"{t_c}, {p_pa}: a split has a vapour fraction"
        # Two steps of slack: the count is where the Newton's own path stops, and a
        # one-ulp change in the fugacity coefficient moves it by a step or two. The
        # sabotage detector is the distance to the scheme it replaced - 244 against
        # 26 - not the last digit.
        assert r.iterations <= hybrid + 2, (
            f"{t_c}, {p_pa}: {r.iterations} steps, and the fallback should reach it in at "
            f"most {hybrid} (successive substitution alone needs {ss_only})"
        )
        h.assert_close(r.beta, beta, 1e-7, f"{t_c}, {p_pa} (beta)")


def test_the_fallback_converges_a_state_the_outer_scheme_abandons() -> None:
    """A state successive substitution cannot converge inside the cap, which the fallback does.

    It converges to the *trivial* solution rather than to a split, and that is the
    answer rather than a failure: NeqSim reports a single phase at this state too.
    """
    r = pt_flash(methane_butane(), T=Q(355.0, "K"), P=Q(1.2e7, "Pa"), z=[0.6, 0.4])
    assert r.phase is Phase.TRIVIAL
    assert r.beta is None, "a trivial solution has no vapour fraction"
