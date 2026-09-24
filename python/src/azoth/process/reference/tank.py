"""``process.tank`` - the tank's kernel.

Spec: ``specs/models/process/tank.toml``

The Python twin of ``crates/azoth-process/src/models/tank.rs`` over
``crates/azoth-process/src/kernels/tank.rs``, written to mirror them.

# The steady state is a flash the fluid already satisfies

``Tank.run`` clones the inlet and calls
``ops.VUflash(thermoSystem2.getVolume(), thermoSystem2.getInternalEnergy())`` - the
*fluid's* volume, not the vessel's. A fluid at its stable state re-imposing its own ``V``
and ``U`` returns that state, so the answer is the flash the feed already carries:
measured, the captured two-phase row is ``process.separator``'s first row to the last
digit at every field of both outlets.

**The design volume does not reach it.** ``setVolume`` is read by ``getVolume``,
``validateMechanicalDesign``, ``toJson`` and ``displayResult``, and by nothing in ``run``;
the capture's two rows that differ only in ``setVolume`` are identical. The palette entry
declares no ``volume`` for that reason.

**And a tank is not a three-phase machine.** ``run`` never enables ``multiPhaseCheck``, and
its two tests name ``gas`` and ``oil``. On a feed a ``ThreePhaseSeparator`` splits three
ways at the same state, this finds two - the gas leaves carrying a water mole fraction of
``0.234`` against the three-phase vessel's ``0.0016``.
"""

from __future__ import annotations

from azoth.core.range import apply_checks, checks_for
from azoth.core.result import TankResult
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning
from azoth.process.reference.mixer import _route as mix_route
from azoth.process.reference.separator import _route as sep_route


def tank(
    components: list[str],
    feed_n: list[Q],
    feed_z: list[list[float]],
    feed_p: list[Q],
    feed_t: list[Q],
) -> TankResult:
    """Join a tank's inlets and split the result into a gas and a liquid outlet.

    Args:
        components: the fluid's substances, by name, shared by every feed.
        feed_n: each inlet's molar flow.
        feed_z: each inlet's composition, one row per feed.
        feed_p: each inlet's pressure.
        feed_t: each inlet's temperature.

    Returns:
        Both outlets' records.

    Raises:
        InvalidInputError: where the feeds' shapes disagree or carry no flow in total.
        OutOfRangeError: for a state the spec bounds.

    Example:
        >>> import azoth
        >>> q = azoth.ureg.Quantity
        >>> r = tank(
        ...     ["methane", "n-butane"],
        ...     [q(1.0, "mol/s")],
        ...     [[0.7, 0.3]],
        ...     [q(20.0, "bar")],
        ...     [q(300.0, "K")],
        ... )
        >>> round(r.gas_n.to("mol/s").magnitude, 6)
        0.818221
    """
    spec = _spec()
    checks = checks_for(spec)
    warnings: list[Warning] = []

    n = [input_to_si(spec, "feed_n", value) for value in feed_n]
    p = [input_to_si(spec, "feed_p", value) for value in feed_p]
    t = [input_to_si(spec, "feed_t", value) for value in feed_t]
    # The feeds are one port, so the bound the family puts on a temperature resolves the
    # first - `manifold`'s and `component_splitter`'s convention, and this model's one
    # input it can check.
    apply_checks(checks.on_input, {"feed_t": t[0] if t else None}.get, warnings)

    # **The whole kernel is these two calls, and both are the *routes*.** A tank's inlet is
    # a mixer, so the feeds join first; it holds the joined pressure and temperature, so the
    # flash is a plain one with no drop and no duty. Writing either out again here would be a
    # second implementation of a ported kernel. The feeds cross as quantities and the routes
    # take SI, which is `manifold`'s arrangement.
    #
    # **The routes and not the models is the point.** A warning is a model's statement about
    # its own parameter surface - `mixer` reports a skipped `outlet_pressure` bound, which the
    # tank has no input for - so a sub-model's warning is not the tank's. The cross-impl test
    # is what caught this: propagating the mixer's put a skipped bound on an input the tank
    # does not have onto the tank's own result.
    joined = mix_route(components, n, feed_z, p, t, None)
    split = sep_route(
        components,
        joined.n_total,
        list(joined.z),
        joined.pressure,
        joined.product_t,
        0.0,
        0.0,
        None,
    )

    return TankResult(
        gas_n=from_si(split.n_vapour, "mol/s"),
        gas_z=split.y,
        gas_p=from_si(split.p_out, "Pa"),
        gas_t=split.temperature,
        gas_h=from_si(split.vapour_h, "J/mol"),
        liquid_n=from_si(split.n_liquid, "mol/s"),
        liquid_z=split.liquid_z,
        liquid_p=from_si(split.p_out, "Pa"),
        liquid_t=split.temperature,
        liquid_h=from_si(split.liquid_h, "J/mol"),
        warnings=tuple(warnings),
    )


def _spec() -> dict[str, object]:
    """The generated spec, imported here so the module has no import-time cycle."""
    from azoth._models_gen import model

    return model("process.tank")
