"""Tests for the ``process.splitter`` unit operation.

The spec's cases pin the answers, and they are the weakest of the checks here.

What the model claims is that **material is conserved** and that **each branch is its
fraction of the feed** - `sum(flows) == n` and `flows[k] == fractions[k] * n`. Both are
arithmetic rather than physics, and both are what a splitter *is*: it is the one unit
operation here with no energy balance and no ideal-gas datum, so there is nothing else it
could be getting wrong.

What is *not* claimed, and has its own test below, is that a fraction set which does not
account for the whole feed is refused rather than normalised. NeqSim normalises; a
splitter's entire output is those numbers, so rescaling them silently would make a
caller's arithmetic error invisible while changing every downstream number.
"""

from __future__ import annotations

from typing import Any

import pytest
from _process import Q, a_mixture, assert_case, call, h

from azoth import _models_gen
from azoth.core.errors import InvalidInputError, OutOfRangeError

MODEL = _models_gen.model("process.splitter")
CASES = MODEL["cases"]


def run(
    *,
    T: float = 300.0,
    P: float = 20.0e5,
    n: float = 10.0,
    fractions: list[float] | None = None,
) -> Any:
    """Call the model with the spec's own mixture. **No ideal-gas datum** - it takes none."""
    import azoth

    return azoth.process.splitter(
        a_mixture(),
        T=Q(T, "K"),
        P=Q(P, "Pa"),
        n=Q(n, "mol/s"),
        z=[0.6, 0.4],
        fractions=[0.3, 0.7] if fractions is None else fractions,
    )


@pytest.mark.parametrize("case", CASES, ids=lambda c: c["id"])
def test_spec_case(case: dict[str, Any]) -> None:
    assert_case(MODEL, case)


def test_every_case_ran() -> None:
    assert CASES, "the spec declares no cases"
    for case in CASES:
        call(MODEL, case)


@pytest.mark.parametrize(
    "fractions",
    [[0.3, 0.7], [1.0], [0.0, 1.0], [0.25, 0.25, 0.5], [0.0, 0.0, 0.0, 1.0]],
)
def test_material_is_conserved(fractions: list[float]) -> None:
    """The branches account for the whole feed, at every fraction set.

    The analogue of the separator's mole balance, and the check that catches a set
    applied twice or not at all. It holds at every answer rather than at one recorded one.
    """
    result = run(n=10.0, fractions=fractions)
    h.assert_close(sum(result.flows), 10.0, 1.0e-12, f"sum of flows with {fractions}")


def test_each_branch_is_its_fraction_of_the_feed() -> None:
    """``flows[k] == fractions[k] * n`` - the definition, and all the model computes."""
    result = run(n=7.0, fractions=[0.2, 0.5, 0.3])
    expected = [1.4, 3.5, 2.1]
    assert len(result.flows) == len(expected)
    for index, (got, want) in enumerate(zip(result.flows, expected, strict=True)):
        h.assert_close(got, want, 1.0e-12, f"flows[{index}]")


def test_the_branches_are_all_at_the_feed_state() -> None:
    """One temperature, one pressure, one phase - because every branch is the same state.

    This is the model's one simplification and it is worth pinning: NeqSim flashes each
    branch because its streams are mutable objects a caller may have written to. Here
    they cannot differ, so one flash answers for all of them.
    """
    result = run(T=300.0, P=20.0e5)

    h.assert_close(result.T.to("K").magnitude, 300.0, 1.0e-12, "branch temperature")
    h.assert_close(result.P.to("Pa").magnitude, 20.0e5, 1.0e-12, "branch pressure")


def test_the_split_is_the_feeds_own_flash() -> None:
    """The phase state reported is ``eos.pt_flash``'s at the feed state, computed here.

    The cross-model check: with the branches indistinguishable, this model's `beta` and
    `phase` must be the flash's for the same state. It is computed from the other model
    rather than recorded, so the two cannot drift together.
    """
    import azoth

    flash = azoth.eos.pt_flash(a_mixture(), T=Q(300.0, "K"), P=Q(20.0e5, "Pa"), z=[0.6, 0.4])
    result = run()

    assert result.phase == flash.phase
    assert result.beta is not None
    assert flash.beta is not None
    h.assert_close(result.beta, flash.beta, 1.0e-12, "vapour fraction against the flash")


@pytest.mark.parametrize("fractions", [[0.3, 0.6], [1.1], [0.5, 0.4, 0.2], [0.0]])
def test_fractions_that_do_not_sum_to_one_are_refused(fractions: list[float]) -> None:
    """Checked rather than corrected. See the module docstring for why."""
    with pytest.raises(InvalidInputError, match="sum"):
        run(fractions=fractions)


def test_a_negative_fraction_is_refused() -> None:
    """A branch that gives back more than it takes is a mixer written backwards."""
    with pytest.raises(InvalidInputError, match="negative"):
        run(fractions=[1.5, -0.5])


def test_no_branches_is_refused() -> None:
    """A splitter with nothing to split into."""
    with pytest.raises(InvalidInputError, match="no branches"):
        run(fractions=[])


def test_a_zero_feed_is_refused() -> None:
    """The fractions would have nothing to act on."""
    with pytest.raises(OutOfRangeError):
        run(n=0.0)


@pytest.mark.requires_rust
def test_the_two_implementations_agree_with_matching_iteration_counts() -> None:
    """Both implementations reach the same split, in the same number of steps."""
    import azoth
    from azoth._dispatch import use_backend

    arguments: dict[str, Any] = {
        "mixture": a_mixture(),
        "T": Q(300.0, "K"),
        "P": Q(20.0e5, "Pa"),
        "n": Q(10.0, "mol/s"),
        "z": [0.6, 0.4],
        "fractions": [0.3, 0.7],
    }
    with use_backend("python"):
        py = azoth.process.splitter(**arguments)
    with use_backend("rust"):
        rs = azoth.process.splitter(**arguments)

    assert rs.phase == py.phase, "phase disagrees"
    assert rs.iterations == py.iterations, (
        f"python took {py.iterations} steps and Rust {rs.iterations}"
    )
    for index, (a, b) in enumerate(zip(rs.flows, py.flows, strict=True)):
        h.assert_close(a, b, 1.0e-12, f"flows[{index}]")
