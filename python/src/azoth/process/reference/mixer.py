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

from dataclasses import dataclass

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

    states = _route(
        components,
        [input_to_si(spec, "feed_n", value) for value in feed_n],
        feed_z,
        [input_to_si(spec, "feed_p", value) for value in feed_p],
        [input_to_si(spec, "feed_t", value) for value in feed_t],
        None if outlet_pressure is None else input_to_si(spec, "outlet_pressure", outlet_pressure),
    )
    return MixerResult(
        product_n=from_si(states.n_total, "mol/s"),
        product_z=states.z,
        product_p=from_si(states.pressure, "Pa"),
        product_t=states.product_t,
        product_h=from_si(states.h_out, "J/mol"),
        warnings=tuple(warnings),
    )


@dataclass(frozen=True, slots=True)
class MixStates:
    """The mixer's intermediates, in the order it computes them, in SI.

    ``h_total`` is the one number here the record does not carry: an outlet's molar
    enthalpy is this over the total flow, and a mean of the right magnitude hides a wrong
    weighting. NeqSim exposes it by name as ``Mixer.calcMixStreamEnthalpy``.
    """

    n_total: float
    z: tuple[float, ...]
    h_total: float
    h_out: float
    pressure: float
    product_t: Q


def _route(
    components: list[str],
    feed_n: list[float],
    feed_z: list[list[float]],
    feed_p: list[float],
    feed_t: list[float],
    outlet_pressure: float | None,
) -> MixStates:
    """The mixer's intermediates, in SI, in the order it computes them.

    **One arithmetic, two consumers.** :func:`mixer` builds the result from this and
    :mod:`azoth.process.layers` reports the same numbers against a NeqSim capture, so a
    layer the harness compares is a layer the model actually took.
    """
    mixture, ideal_gas = _components.mixture_of(components, eos="pr")

    n_total = sum(feed_n)
    # **Feeds that carry nothing are refused rather than divided by.** The outlet
    # composition is the flows' weighted average, so a zero total is a 0/0 and every field
    # of the outlet would be NaN. NeqSim propagates a zero-flow outlet instead, through an
    # `isActive` flag `Stream` has no equivalent of; a caller with a branch that may be idle
    # has to say what it wants, because a NaN state is not an answer. The refusal is here
    # rather than in :func:`mixer` because `kernels/mixer.rs` is where Rust puts it, and a
    # model that reaches the mixer's arithmetic - `tank` and `manifold` do - reaches it too.
    if n_total <= 0.0:
        raise InvalidInputError(
            "feed_n",
            f"a mixer's feeds carry {n_total} mol/s in total, so there is no mixture to report",
        )
    z = [0.0] * len(feed_z[0])
    h_total = 0.0
    for index, n_i in enumerate(feed_n):
        h_i, _ = enthalpy_at(mixture, ideal_gas, feed_t[index], feed_p[index], feed_z[index])
        for i, zi in enumerate(feed_z[index]):
            z[i] += n_i * zi
        h_total += n_i * h_i
    z = [zi / n_total for zi in z]

    h_out = h_total / n_total
    # Carried as a quantity rather than as a magnitude: `ph_flash` takes quantities, where
    # `enthalpy_at` above takes SI. The two conventions are the existing surface's, not this
    # model's choice.
    if outlet_pressure is None:
        pressure = min(feed_p)
        pressure_q: Q = from_si(pressure, "Pa")
    else:
        pressure = outlet_pressure
        pressure_q = from_si(pressure, "Pa")
    solved = ph_flash_solve(mixture, ideal_gas, pressure_q, from_si(h_out, "J/mol"), z)

    return MixStates(
        n_total=n_total,
        z=tuple(z),
        h_total=h_total,
        h_out=h_out,
        pressure=pressure,
        product_t=solved.T,
    )


def _spec() -> dict[str, object]:
    """The generated spec, imported here so the module has no import-time cycle."""
    from azoth._models_gen import model

    return model("process.mixer")
