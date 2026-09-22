"""``process.mixer`` - the mixer's kernel.

Spec: ``specs/models/process/mixer.toml``

The Python twin of ``crates/azoth-process/src/models/mixer.rs`` over
``crates/azoth-process/src/kernels/mixer.rs``, written to mirror them.

# Mixing is isenthalpic, and the pressure is the lowest feed's

Mixing adds no heat and no work, so the outlet's molar enthalpy is the feeds' weighted by
their flows - and the outlet temperature is the one at which the joined mixture carries it
at the outlet pressure, solved by a ``PHflash``. That is why the outlet is at neither
feed's temperature in general: the feeds may be at different pressures, and the joined
state is at one.

The outlet pressure is the lowest feed pressure unless ``outlet_pressure`` states one, for
the reason a mixer cannot raise its own pressure. NeqSim's ``Mixer.run`` does the same -
it takes the minimum over the active inlets and lets an explicit outlet pressure override
it - and both cases exercise one of the two.
"""

from __future__ import annotations

from azoth.core.errors import InvalidInputError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import MixerResult
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning
from azoth.eos import components as _components
from azoth.eos.reference.ph_flash import enthalpy_at
from azoth.eos.reference.ph_flash import ph_flash as ph_flash_solve


def mixer(
    components: list[str],
    feed_n: list[Q],
    feed_z: list[list[float]],
    feed_p: list[Q],
    feed_t: list[Q],
    outlet_pressure: Q | None = None,
) -> MixerResult:
    """Join several streams into one, conserving molar flow and enthalpy.

    Args:
        components: the fluid's substances, by name, shared by every feed.
        feed_n: each feed's molar flow.
        feed_z: each feed's composition, one row per feed.
        feed_p: each feed's pressure.
        feed_t: each feed's temperature.
        outlet_pressure: a specified outlet pressure, or ``None`` for the lowest feed's.

    Returns:
        The outlet's record: the feeds' flows summed, their compositions weighted by
        flow, the outlet pressure, and the temperature and molar enthalpy the joined
        state implies.

    Raises:
        InvalidInputError: where the feeds' shapes disagree or the outlet pressure is
            not positive.

    Example:
        >>> import azoth
        >>> q = azoth.ureg.Quantity
        >>> r = mixer(
        ...     ["methane", "n-butane"],
        ...     [q(1.0, "mol/s"), q(2.0, "mol/s")],
        ...     [[1.0, 0.0], [1.0, 0.0]],
        ...     [q(10.0, "bar"), q(6.0, "bar")],
        ...     [q(300.0, "K"), q(300.0, "K")],
        ... )
        >>> round(r.product_n.to("mol/s").magnitude, 12)
        3.0
    """
    spec = _spec()
    checks = checks_for(spec)
    warnings: list[Warning] = []

    p_out = (
        None if outlet_pressure is None else input_to_si(spec, "outlet_pressure", outlet_pressure)
    )
    # An omitted optional input leaves its bound unrun, which `apply_checks` reports as a
    # skip rather than a pass - the reason the bound is declared on this input at all.
    apply_checks(checks.on_input, {"outlet_pressure": p_out}.get, warnings)

    if not feed_n:
        raise InvalidInputError("feed_n", "a mixer needs at least one feed")
    for name, got in (("feed_z", len(feed_z)), ("feed_p", len(feed_p)), ("feed_t", len(feed_t))):
        if got != len(feed_n):
            raise InvalidInputError(
                name,
                f"a mixer's feeds are one port, so every field has one entry per feed: "
                f"`feed_n` has {len(feed_n)} and `{name}` has {got}",
            )

    mixture, ideal_gas = _components.mixture_of(components, eos="pr")

    n_total = sum(input_to_si(spec, "feed_n", value) for value in feed_n)
    n_components = len(feed_z[0])
    z = [0.0] * n_components
    h_total = 0.0
    for index, flow in enumerate(feed_n):
        n_i = input_to_si(spec, "feed_n", flow)
        t_i = input_to_si(spec, "feed_t", feed_t[index])
        p_i = input_to_si(spec, "feed_p", feed_p[index])
        h_i, _ = enthalpy_at(mixture, ideal_gas, t_i, p_i, feed_z[index])
        for i, zi in enumerate(feed_z[index]):
            z[i] += n_i * zi
        h_total += n_i * h_i
    z = [zi / n_total for zi in z]

    # Carried as the caller's quantity rather than as a magnitude: `ph_flash` takes
    # quantities, where `enthalpy_at` above takes SI. The two conventions are the existing
    # surface's, not this model's choice.
    if outlet_pressure is None:
        pressure = min(feed_p, key=lambda value: input_to_si(spec, "feed_p", value))
        pressure_si = input_to_si(spec, "feed_p", pressure)
    else:
        pressure = outlet_pressure
        pressure_si = input_to_si(spec, "outlet_pressure", outlet_pressure)
    h_out = h_total / n_total
    solved = ph_flash_solve(mixture, ideal_gas, pressure, from_si(h_out, "J/mol"), z)

    return MixerResult(
        product_n=from_si(n_total, "mol/s"),
        product_z=tuple(z),
        product_p=from_si(pressure_si, "Pa"),
        product_t=solved.T,
        product_h=from_si(h_out, "J/mol"),
        warnings=tuple(warnings),
    )


def _spec() -> dict[str, object]:
    """The generated spec, imported here so the module has no import-time cycle."""
    from azoth._models_gen import model

    return model("process.mixer")
