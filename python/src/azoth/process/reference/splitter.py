"""``process.splitter`` - the splitter's kernel.

Spec: ``specs/models/process/splitter.toml``

The Python twin of ``crates/azoth-process/src/models/splitter.rs`` over
``crates/azoth-process/src/kernels/splitter.rs``, written to mirror them.

# Nothing changes but the flow

A splitter adds no heat and no work and does not flash. Every outlet therefore carries
the feed's composition, pressure, temperature and molar enthalpy, and the fractions are
the whole of the arithmetic. That is also why the outlet enthalpy is not something an
oracle can be asked for directly - see the divergence below.

# NeqSim's ``Splitter.run`` does not conserve enthalpy

It reaches its branch state through ``clone()``, then ``init(0)`` - which leaves a
``SystemPrEos`` holding *two* phase objects, each with the full inventory - then
``addComponent(int, double)``, which mutates every phase in ``phaseArray`` while
``setTotalNumberOfMolesRaw`` counts the change once. The ``TPflash`` after that converges
the active phase and leaves the stale one in place, and ``getEnthalpy()`` sums what it
finds there. Measured on the captured fluid: -20689.86 J/mol in, -93092.87 and -28142.87
out, weighting to -47627.87. This diverges deliberately and the spec says so.
"""

from __future__ import annotations

from dataclasses import dataclass

from azoth.core.errors import InvalidInputError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import SplitterResult
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning
from azoth.eos import components as _components
from azoth.eos.reference.ph_flash import enthalpy_at


def splitter(
    components: list[str],
    feed_n: Q,
    feed_z: list[float],
    feed_p: Q,
    feed_t: Q,
    split_factors: list[float],
) -> SplitterResult:
    """Split a stream into several with the same state, in proportion to ``split_factors``.

    Args:
        components: the fluid's substances, by name.
        feed_n: the inlet's molar flow.
        feed_z: the inlet composition.
        feed_p: the inlet pressure.
        feed_t: the inlet temperature.
        split_factors: the fraction of the feed to each outlet. Normalised, so only the
            ratios matter and the vector's length is the outlet count.

    Returns:
        Every outlet's record: its share of the flow, and the feed's composition,
        pressure, temperature and molar enthalpy.

    Raises:
        InvalidInputError: where there is no outlet, or the factors do not sum to a
            positive value.

    Example:
        >>> import azoth
        >>> q = azoth.ureg.Quantity
        >>> r = splitter(
        ...     ["n-butane", "n-pentane"],
        ...     q(1.0, "mol/s"),
        ...     [0.6, 0.4],
        ...     q(10.0, "bar"),
        ...     q(300.0, "K"),
        ...     [0.3, 0.7],
        ... )
        >>> [round(n.to("mol/s").magnitude, 12) for n in r.products_n]
        [0.3, 0.7]
    """
    spec = _spec()
    checks = checks_for(spec)
    warnings: list[Warning] = []

    n = input_to_si(spec, "feed_n", feed_n)
    t = input_to_si(spec, "feed_t", feed_t)
    p = input_to_si(spec, "feed_p", feed_p)

    apply_checks(checks.on_input, {"feed_t": t}.get, warnings)

    if not split_factors:
        raise InvalidInputError("split_factors", "a splitter needs at least one outlet")

    total = sum(split_factors)
    if total <= 0.0:
        raise InvalidInputError("split_factors", "the split fractions must sum to a positive value")

    states = _route(components, t, p, feed_z, split_factors)
    z = tuple(feed_z)
    return SplitterResult(
        products_n=tuple(from_si(n * fraction, "mol/s") for fraction in states.fractions),
        products_z=tuple(z for _ in states.fractions),
        products_p=tuple(feed_p for _ in states.fractions),
        products_t=tuple(feed_t for _ in states.fractions),
        products_h=tuple(from_si(states.h_in, "J/mol") for _ in states.fractions),
        warnings=tuple(warnings),
    )


@dataclass(frozen=True, slots=True)
class SplitStates:
    """The feed's own state and the share each outlet takes of it, in SI.

    The refusals are not here: they are about the *declared* factors, and
    :func:`splitter` states them before the fluid is resolved. What is here is what the
    split is - one enthalpy, and the normalised fractions that divide the flow.
    """

    h_in: float
    fractions: tuple[float, ...]


def _route(
    components: list[str],
    t: float,
    p: float,
    feed_z: list[float],
    split_factors: list[float],
) -> SplitStates:
    """The splitter's intermediates, in the order it computes them, in SI.

    **One arithmetic, two consumers.** :func:`splitter` builds the result from this and
    :mod:`azoth.process.layers` reports the same numbers against a NeqSim capture, so a
    layer the harness compares is a layer the model actually took. The feed's enthalpy is
    the whole interior: a split is a flow division and nothing else, which is exactly why
    the branch enthalpies the capture carries are the row where the two libraries part.
    """
    mixture, ideal_gas = _components.mixture_of(components, eos="pr")
    # The feed's own state, which every outlet carries: `h` is a function of `(T, P, z)`
    # and none of the three moves.
    h_in, _ = enthalpy_at(mixture, ideal_gas, t, p, feed_z)
    total = sum(split_factors)
    return SplitStates(h_in=h_in, fractions=tuple(factor / total for factor in split_factors))


def _spec() -> dict[str, object]:
    """The generated spec, imported here so the module has no import-time cycle."""
    from azoth._models_gen import model

    return model("process.splitter")
