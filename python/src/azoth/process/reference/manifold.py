"""``process.manifold`` - a mixer and a splitter composed, with a low-flow rule between.

Spec: ``specs/models/process/manifold.toml``

The Python twin of ``crates/azoth-process/src/models/manifold.rs`` over
``crates/azoth-process/src/kernels/manifold.rs``, written to mirror them.

# It is the two kernels and not a third arithmetic

``Manifold.run`` is four statements: ``propagateMinimumFlow()``, ``localmixer.run()``,
``refreshLocalSplitter()`` and ``localsplitter.run()``. So this calls
:func:`azoth.process.reference.mixer` and :func:`azoth.process.reference.splitter` rather
than restating either, and every assumption those two carry holds here - including that
NeqSim's branch enthalpies do not conserve energy and neither of them reproduces them.

# The one rule that is the manifold's own

An inlet whose mass flow is at or below `ProcessEquipmentBaseClass.DEFAULT_MINIMUM_FLOW`
(`1e-20` kg/hr) is not mixed. Measured, a feed stated at zero flow leaves the mixture as
the other feed alone, flow and composition both - the captured third row.
"""

from __future__ import annotations

from azoth.core.errors import InvalidInputError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import ManifoldResult
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning
from azoth.process.reference.mixer import _route as mix_route
from azoth.process.reference.splitter import _route as split_route

#: `ProcessEquipmentBaseClass.DEFAULT_MINIMUM_FLOW`, in kg/hr. The palette declares no way
#: to change it, so this is the whole of the rule here.
DEFAULT_MINIMUM_FLOW_KG_PER_HOUR = 1e-20


def manifold(
    components: list[str],
    feed_n: list[Q],
    feed_z: list[list[float]],
    feed_p: list[Q],
    feed_t: list[Q],
    split_factors: list[float],
) -> ManifoldResult:
    """Join a manifold's feeds and divide the mixture between its outlets.

    Args:
        components: the fluid's substances, by name, shared by every feed.
        feed_n: each feed's molar flow.
        feed_z: each feed's composition, one row per feed.
        feed_p: each feed's pressure.
        feed_t: each feed's temperature.
        split_factors: the fraction of the mixture to each outlet.

    Returns:
        The outlets' records.

    Raises:
        InvalidInputError: where the feeds' shapes disagree, a factor is negative, or no
            feed is above the low-flow threshold.

    Example:
        >>> import azoth
        >>> q = azoth.ureg.Quantity
        >>> r = manifold(
        ...     ["methane", "n-butane"],
        ...     [q(1.0, "mol/s"), q(2.0, "mol/s")],
        ...     [[0.9, 0.1], [0.3, 0.7]],
        ...     [q(30.0, "bar"), q(10.0, "bar")],
        ...     [q(320.0, "K"), q(300.0, "K")],
        ...     [0.25, 0.75],
        ... )
        >>> [round(n.to("mol/s").magnitude, 3) for n in r.products_n]
        [0.75, 2.25]
    """
    spec = _spec()
    checks = checks_for(spec)
    warnings: list[Warning] = []

    n = [input_to_si(spec, "feed_n", value) for value in feed_n]
    p = [input_to_si(spec, "feed_p", value) for value in feed_p]
    t = [input_to_si(spec, "feed_t", value) for value in feed_t]

    apply_checks(
        checks.on_input,
        {
            "feed_t": t[0] if t else None,
            "split_factors": split_factors[0] if split_factors else None,
        }.get,
        warnings,
    )

    # The low-flow filter, in kg/s against a threshold NeqSim states in kg/hr.
    keep = [
        index
        for index in range(len(n))
        if _mass_flow(components, feed_z[index], n[index]) * 3600.0
        > DEFAULT_MINIMUM_FLOW_KG_PER_HOUR
    ]
    if not keep:
        raise InvalidInputError(
            "feeds",
            "a manifold with no feed above the low-flow threshold has nothing to divide",
        )

    mixture = mix_route(
        components,
        [n[index] for index in keep],
        [list(feed_z[index]) for index in keep],
        [p[index] for index in keep],
        [t[index] for index in keep],
        None,
    )
    # The splitter's route takes the mixture as a *state*, which is what the manifold's
    # `refreshLocalSplitter` does by re-attaching its inlet to the mixer's outlet.
    branches = split_route(
        components,
        mixture.product_t.to("K").magnitude,
        mixture.pressure,
        list(mixture.z),
        split_factors,
    )

    return ManifoldResult(
        products_n=tuple(
            from_si(mixture.n_total * fraction, "mol/s") for fraction in branches.fractions
        ),
        products_z=tuple(tuple(mixture.z) for _ in branches.fractions),
        products_p=tuple(from_si(mixture.pressure, "Pa") for _ in branches.fractions),
        products_t=tuple(mixture.product_t for _ in branches.fractions),
        products_h=tuple(from_si(mixture.h_out, "J/mol") for _ in branches.fractions),
        warnings=tuple(warnings),
    )


def _mass_flow(components: list[str], z: list[float], n: float) -> float:
    """A feed's mass flow in kg/s, which the low-flow rule is stated in."""
    from azoth.eos import components as databank

    fluid = databank.mixture_of(components, eos="pr")[0]
    masses = [
        zi * component.molar_mass.to("kg/mol").magnitude  # type: ignore[union-attr]
        for zi, component in zip(z, fluid.components, strict=True)
    ]
    return n * sum(masses)


def _spec() -> dict[str, object]:
    """The generated spec, imported here so the module has no import-time cycle."""
    from azoth._models_gen import model

    return model("process.manifold")
