"""``process.component_splitter`` - a per-component routing, with each outlet flashed.

Spec: ``specs/models/process/component_splitter.toml``

The Python twin of ``crates/azoth-process/src/models/component_splitter.rs`` over
``crates/azoth-process/src/kernels/component_splitter.rs``, written to mirror them.

# The factor is per component and the outlet count is two

``ComponentSplitter.run`` loops ``for i in 0..2`` and reads ``splitFactor[k]`` as the
fraction of component *k* routed to the overhead, so the outlets carry different
*compositions* - which is what separates this from ``process.splitter``.

# Each outlet is flashed, and the class's enthalpies are not

``run`` sets the component moles on an empty fluid, calls ``init(0)`` and runs a ``TPflash``,
so each outlet is a state of its own. The class then *reports* enthalpies that are not those
states: a fresh NeqSim fluid at the first row's bottoms gives ``-20309.24832914077`` where the
class gives ``-14431.31737097011``. This computes the state's, as ``process.splitter`` does
for the same defect.
"""

from __future__ import annotations

from azoth.core.errors import InvalidInputError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import ComponentSplitterResult
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning
from azoth.process.reference.splitter import _route as split_route


def component_splitter(
    components: list[str],
    feed_n: Q,
    feed_z: list[float],
    feed_p: Q,
    feed_t: Q,
    split_factors: list[float],
) -> ComponentSplitterResult:
    """Divide a stream between two outlets component by component.

    Args:
        components: the fluid's substances, by name.
        feed_n: the inlet's molar flow.
        feed_z: the inlet composition.
        feed_p: the inlet pressure, which both outlets flash at.
        feed_t: the inlet temperature, which both outlets flash at.
        split_factors: the fraction of each component to the overhead, one per component.

    Returns:
        Both outlets' records.

    Raises:
        InvalidInputError: where the factors are not one per component, or any is outside
            ``[0, 1]``.

    Example:
        >>> import azoth
        >>> q = azoth.ureg.Quantity
        >>> r = component_splitter(
        ...     ["methane", "n-butane", "n-pentane"],
        ...     q(1.0, "mol/s"),
        ...     [0.5, 0.3, 0.2],
        ...     q(20.0, "bar"),
        ...     q(300.0, "K"),
        ...     [0.98, 0.05, 0.02],
        ... )
        >>> round(r.overhead_n.to("mol/s").magnitude, 3)
        0.509
    """
    spec = _spec()
    checks = checks_for(spec)
    warnings: list[Warning] = []

    n = input_to_si(spec, "feed_n", feed_n)
    p = input_to_si(spec, "feed_p", feed_p)
    t = input_to_si(spec, "feed_t", feed_t)

    apply_checks(
        checks.on_input,
        {"feed_t": t, "split_factors": split_factors[0] if split_factors else None}.get,
        warnings,
    )

    if len(split_factors) != len(feed_z):
        raise InvalidInputError(
            "split_factors",
            f"a component split routes each component, so the factors are one per component: "
            f"{len(split_factors)} factors for {len(feed_z)} components",
        )
    for index, factor in enumerate(split_factors):
        if not 0.0 <= factor <= 1.0:
            raise InvalidInputError(
                "split_factors",
                f"a component's split fraction is a fraction, and component {index}'s is {factor}",
            )

    def outlet(amounts: list[float]) -> tuple[float, list[float]]:
        total = sum(amounts)
        if total <= 0.0:
            return 0.0, list(feed_z)
        return total, [amount / total for amount in amounts]

    overhead_amounts = [n * zi * f for zi, f in zip(feed_z, split_factors, strict=True)]
    bottoms_amounts = [n * zi * (1.0 - f) for zi, f in zip(feed_z, split_factors, strict=True)]
    overhead_n, overhead_z = outlet(overhead_amounts)
    bottoms_n, bottoms_z = outlet(bottoms_amounts)

    # Each outlet is flashed at the feed's own temperature and pressure.
    overhead = split_route(components, t, p, overhead_z, [1.0])
    bottoms = split_route(components, t, p, bottoms_z, [1.0])

    return ComponentSplitterResult(
        overhead_n=from_si(overhead_n, "mol/s"),
        overhead_z=tuple(overhead_z),
        overhead_p=from_si(p, "Pa"),
        overhead_t=from_si(t, "K"),
        overhead_h=from_si(overhead.h_in, "J/mol"),
        bottoms_n=from_si(bottoms_n, "mol/s"),
        bottoms_z=tuple(bottoms_z),
        bottoms_p=from_si(p, "Pa"),
        bottoms_t=from_si(t, "K"),
        bottoms_h=from_si(bottoms.h_in, "J/mol"),
        warnings=tuple(warnings),
    )


def _spec() -> dict[str, object]:
    """The generated spec, imported here so the module has no import-time cycle."""
    from azoth._models_gen import model

    return model("process.component_splitter")
