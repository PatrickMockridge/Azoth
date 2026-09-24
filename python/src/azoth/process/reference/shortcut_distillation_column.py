"""``process.shortcut_distillation_column`` - the FUG column's kernel.

Spec: ``specs/models/process/shortcut_distillation_column.toml``

The Python twin of ``crates/azoth-process/src/models/shortcut_distillation_column.rs`` over
``crates/azoth-process/src/kernels/shortcut_distillation_column.rs``, written to mirror them.

# The class is closed form, and that is why it is ported first

Fenske's minimum stages, Underwood's minimum reflux, Gilliland's actual stages through
Molokanov's fit and Kirkbride's feed tray are all rearrangements of one vector of K-values.
The only loop is Underwood's root, bisected on ``[1 + 1e-6, alpha_LK/HK - 1e-6]`` against the
class's own ``1e-10`` and its own 200 iterations. So the two libraries can be expected to
agree tightly here, and they do - measured to ``3e-12`` on the captured row - which makes a
disagreement in this model a plumbing bug rather than a physics question.

# Three things the class does that a reader would not guess

**The feed quality is the flash's split, and one phase is decided by the phase type.** One
phase gives ``q = 0`` for a gas and ``q = 1`` for anything else; two give the liquid's share
of the moles, which on a molar flash is ``1 - beta``.

**The Wilson fallback constants are the class's.** ``(Pc/P) exp(5.373 (1 + omega) (1 - Tc/T))``
- the class carries ``5.373``, not Wilson's own ``5.37``, and it is that constant that is here.

**The duties are estimates.** ``V = D(R + 1)`` with a hard-coded ``30000`` J/mol average
latent heat for the condenser, and the reboiler takes the condenser's figure with its sign
changed plus one per cent of the feed's **total** enthalpy in W. A condenser duty that does
not depend on what the fluid is is not an energy balance, and it is reported as the class
reports it.

# Why the arithmetic is factored

:func:`_states` is the whole calculation in SI magnitudes and returns every intermediate the
class forms, and :func:`shortcut_distillation_column` is the unit conversion and the result
around it. ``azoth.process.layers`` reads the factored entry point rather than re-deriving
the K-values, so the layer the diff names is the one the kernel took.
"""

from __future__ import annotations

import math
from typing import NamedTuple

from azoth.core.errors import InvalidInputError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import Phase, ShortcutDistillationColumnResult
from azoth.core.units import Q, from_si, input_to_si, quantity
from azoth.core.warnings import Warning
from azoth.eos import components as _components
from azoth.eos.reference._mixture_state import wilson_k
from azoth.eos.reference.ph_flash import enthalpy_at
from azoth.eos.reference.pt_flash import pt_flash

#: The class's iteration cap for Underwood's bisection.
UNDERWOOD_ITERATIONS = 200
#: The class's stopping rule on the function's value.
UNDERWOOD_TOLERANCE = 1.0e-10
#: The class's second stopping rule: a bracket this narrow ends the bisection.
UNDERWOOD_BRACKET = 1.0e-12
#: The class's split-fraction floor, from ``boundedSplitFraction``.
SPLIT_FLOOR = 1.0e-12
#: The class's latent-heat stand-in, J/mol, in ``computeDuties``.
AVERAGE_LATENT_HEAT = 30000.0


class _States(NamedTuple):
    """Everything the class forms, in SI magnitudes.

    The K-values and the relative volatilities are here because the capture has them -
    ``ProcessProbe`` re-runs the class's own flash beside it - and they are the layer that
    localises a divergence: everything downstream of ``alpha`` is rearrangement, and
    everything upstream of it is the flash.
    """

    #: The K-values the split is written in, from the flash or from Wilson.
    k_values: tuple[float, ...]
    #: ``alpha_i = K_i / K_HK``.
    alpha: tuple[float, ...]
    #: ``alpha_LK/HK``.
    relative_volatility: float
    #: The feed quality Underwood's right-hand side takes.
    q: float
    #: Underwood's root.
    theta: float
    #: Fenske's minimum stages.
    minimum_stages: float
    #: The fraction of each component the class sends to the distillate.
    fractions: tuple[float, ...]
    #: Underwood's minimum reflux.
    minimum_reflux_ratio: float
    #: ``minimum_reflux_ratio * reflux_ratio_multiplier``.
    actual_reflux_ratio: float
    #: Molokanov's stage count.
    actual_stages: float
    #: The feed stage from the top.
    feed_tray_number: int
    #: The condenser duty, an estimate.
    condenser_duty: float
    #: The reboiler duty, an estimate.
    reboiler_duty: float
    #: The distillate's molar flow, mol/s.
    distillate_n: float
    #: The distillate's composition.
    distillate_z: tuple[float, ...]
    #: The distillate's molar enthalpy, J/mol.
    distillate_h: float
    #: The bottoms' molar flow, mol/s.
    bottoms_n: float
    #: The bottoms' composition.
    bottoms_z: tuple[float, ...]
    #: The bottoms' molar enthalpy, J/mol.
    bottoms_h: float


def shortcut_distillation_column(
    components: list[str],
    feed_n: Q,
    feed_z: list[float],
    feed_p: Q,
    feed_t: Q,
    light_key: str,
    heavy_key: str,
    light_key_recovery_distillate: float,
    heavy_key_recovery_bottoms: float,
    reflux_ratio_multiplier: float,
    condenser_pressure: Q | None = None,
    reboiler_pressure: Q | None = None,
) -> ShortcutDistillationColumnResult:
    """Solve a shortcut distillation column.

    Args:
        components: the fluid's substances, by name.
        feed_n: the feed's molar flow.
        feed_z: the feed composition, in ``components``' order.
        feed_p: the feed's pressure, and the default for both products'.
        feed_t: the feed's temperature.
        light_key: the light key's name, one of ``components``.
        heavy_key: the heavy key's name, one of ``components``.
        light_key_recovery_distillate: the light key's recovery in the distillate.
        heavy_key_recovery_bottoms: the heavy key's recovery in the bottoms.
        reflux_ratio_multiplier: the actual reflux as a multiple of Underwood's minimum.
        condenser_pressure: the distillate's pressure, or ``None`` for the feed's.
        reboiler_pressure: the bottoms' pressure, or ``None`` for the feed's.

    Returns:
        Both products' records and the eight scalars the class exposes.

    Raises:
        InvalidInputError: for a key that is not a component, a relative volatility at or
            below one, a reflux multiplier at or below one, a split that empties a product,
            or a feed whose flash is a trivial solution.
        OutOfRangeError: for an input the spec bounds.

    Example:
        >>> import azoth
        >>> q = azoth.ureg.Quantity
        >>> r = shortcut_distillation_column(
        ...     ["methane", "n-butane"],
        ...     q(1.0, "mol/s"),
        ...     [0.5, 0.5],
        ...     q(20.0, "bar"),
        ...     q(300.0, "K"),
        ...     "methane",
        ...     "n-butane",
        ...     0.99,
        ...     0.99,
        ...     1.2,
        ... )
        >>> round(r.minimum_stages, 6)
        2.391643
    """
    spec = _spec()
    checks = checks_for(spec)
    warnings: list[Warning] = []

    n = input_to_si(spec, "feed_n", feed_n)
    t = input_to_si(spec, "feed_t", feed_t)
    p = input_to_si(spec, "feed_p", feed_p)

    apply_checks(
        checks.on_input,
        {
            "reflux_ratio_multiplier": reflux_ratio_multiplier,
            "light_key_recovery_distillate": light_key_recovery_distillate,
            "heavy_key_recovery_bottoms": heavy_key_recovery_bottoms,
            "feed_t": t,
        }.get,
        warnings,
    )

    condenser_si = (
        None if condenser_pressure is None else float(condenser_pressure.to("Pa").magnitude)
    )
    reboiler_si = None if reboiler_pressure is None else float(reboiler_pressure.to("Pa").magnitude)

    states = _states(
        components,
        n,
        t,
        p,
        feed_z,
        light_key,
        heavy_key,
        light_key_recovery_distillate,
        heavy_key_recovery_bottoms,
        reflux_ratio_multiplier,
        condenser_si,
        reboiler_si,
    )

    return ShortcutDistillationColumnResult(
        distillate_n=from_si(states.distillate_n, "mol/s"),
        distillate_z=states.distillate_z,
        distillate_p=from_si(
            p if condenser_si is None or condenser_si <= 0.0 else condenser_si, "Pa"
        ),
        distillate_t=feed_t,
        distillate_h=from_si(states.distillate_h, "J/mol"),
        bottoms_n=from_si(states.bottoms_n, "mol/s"),
        bottoms_z=states.bottoms_z,
        bottoms_p=from_si(p if reboiler_si is None or reboiler_si <= 0.0 else reboiler_si, "Pa"),
        bottoms_t=feed_t,
        bottoms_h=from_si(states.bottoms_h, "J/mol"),
        minimum_stages=states.minimum_stages,
        minimum_reflux_ratio=states.minimum_reflux_ratio,
        actual_stages=states.actual_stages,
        actual_reflux_ratio=states.actual_reflux_ratio,
        feed_tray_number=states.feed_tray_number,
        condenser_duty=from_si(states.condenser_duty, "W"),
        reboiler_duty=from_si(states.reboiler_duty, "W"),
        relative_volatility=states.relative_volatility,
        warnings=tuple(warnings),
    )


def _states(
    components: list[str],
    feed_n: float,
    feed_t: float,
    feed_p: float,
    feed_z: list[float],
    light_key: str,
    heavy_key: str,
    light_key_recovery_distillate: float,
    heavy_key_recovery_bottoms: float,
    reflux_ratio_multiplier: float,
    condenser_pressure: float | None,
    reboiler_pressure: float | None,
) -> _States:
    """The class's ``run``, in SI magnitudes: the one arithmetic the kernel and the dump share."""
    lk = _index_of(components, light_key, "light_key")
    hk = _index_of(components, heavy_key, "heavy_key")
    if reflux_ratio_multiplier <= 1.0:
        raise InvalidInputError(
            "reflux_ratio_multiplier",
            f"a reflux multiplier of {reflux_ratio_multiplier} is at or below the minimum: at "
            f"one the class returns an infinite stage count and at less than one it takes its "
            f"silent Y = 0.5 fallback, and neither is a column",
        )

    mixture, ideal_gas = _components.mixture_of(components, eos="pr")
    count = len(feed_z)
    flash = pt_flash(mixture, quantity(feed_t, "K"), quantity(feed_p, "Pa"), list(feed_z))

    if flash.phase == Phase.TWO_PHASE:
        k_values = [
            flash.y[i] / flash.x[i] if flash.x[i] > 1.0e-20 else 1.0e10 for i in range(count)
        ]
        q_quality = 1.0 - (0.0 if flash.beta is None else float(flash.beta))
    elif flash.phase in (Phase.ALL_VAPOUR, Phase.ALL_LIQUID):
        k_values = wilson_k(mixture, feed_t, feed_p)
        q_quality = 0.0 if flash.phase == Phase.ALL_VAPOUR else 1.0
    else:
        raise InvalidInputError(
            "feed_t / feed_p",
            f"the flash at {feed_t} K and {feed_p} Pa is a trivial solution - every K-value "
            f"straddled one, so nothing here proves whether the feed is a gas or a liquid. "
            f"The class reads its phase type, which decides both the Wilson fallback and the "
            f"feed quality, and this model does not decide it",
        )

    k_hk = k_values[hk]
    relative_volatility = k_values[lk] / k_hk
    if relative_volatility <= 1.0:
        raise InvalidInputError(
            "light_key / heavy_key",
            f"{light_key} is not more volatile than {heavy_key}: alpha_LK/HK is "
            f"{relative_volatility}, and the class logs an error and leaves every answer at "
            f"its field initialiser",
        )
    alpha = [k / k_hk for k in k_values]

    # Fenske, on the recoveries clamped exactly as `boundedSplitFraction` clamps them.
    x_lk_d = _bounded(light_key_recovery_distillate)
    x_hk_d = _bounded(1.0 - heavy_key_recovery_bottoms)
    x_lk_b = _bounded(1.0 - light_key_recovery_distillate)
    x_hk_b = _bounded(heavy_key_recovery_bottoms)
    minimum_stages = math.log((x_lk_d / x_hk_d) * (x_hk_b / x_lk_b)) / math.log(relative_volatility)

    theta = _underwood_root(alpha, feed_z, q_quality, lk)
    fractions = _estimate_distillate_fractions(
        alpha,
        lk,
        hk,
        relative_volatility,
        minimum_stages,
        light_key_recovery_distillate,
        heavy_key_recovery_bottoms,
    )

    x_d = [feed_z[i] * fractions[i] for i in range(count)]
    total_d = sum(x_d)
    if total_d > 0.0:
        x_d = [value / total_d for value in x_d]
    r_min_plus_one = 0.0
    for i in range(count):
        if x_d[i] > 1.0e-15 and abs(alpha[i] - theta) > 1.0e-10:
            r_min_plus_one += alpha[i] * x_d[i] / (alpha[i] - theta)
    minimum_reflux_ratio = max(r_min_plus_one - 1.0, 0.0)

    # Gilliland through Molokanov's fit, and the stages it inverts to.
    actual_reflux_ratio = minimum_reflux_ratio * reflux_ratio_multiplier
    x_gilliland = (actual_reflux_ratio - minimum_reflux_ratio) / (actual_reflux_ratio + 1.0)
    if 0.0 <= x_gilliland <= 1.0:
        y_gilliland = 1.0 - math.exp(
            ((1.0 + 54.4 * x_gilliland) / (11.0 + 117.2 * x_gilliland))
            * ((x_gilliland - 1.0) / math.sqrt(x_gilliland))
        )
    else:
        y_gilliland = 0.5
    actual_stages = (y_gilliland + minimum_stages) / (1.0 - y_gilliland)

    # Kirkbride, on the class's own argument - see the spec's assumptions for how it differs
    # from the correlation as usually quoted.
    z_lk_f = feed_z[lk]
    z_hk_f = feed_z[hk]
    b_over_d = _bottoms_to_distillate(feed_z, fractions)
    x_hk_distillate = x_hk_d * z_hk_f
    x_lk_bottoms = x_lk_b * z_lk_f
    if x_lk_bottoms > 1.0e-15 and b_over_d > 1.0e-15:
        kirkbride_ratio = (
            (z_hk_f / z_lk_f) * (x_lk_bottoms / x_hk_distillate) * (b_over_d * b_over_d)
        ) ** 0.206
    else:
        kirkbride_ratio = 0.0
    feed_tray_number = round(actual_stages / (1.0 + kirkbride_ratio)) + 1

    total_distillate_fraction = min(
        max(sum(feed_z[i] * fractions[i] for i in range(count)), 0.001), 0.999
    )
    distillate_flow = feed_n * total_distillate_fraction
    vapour_flow = distillate_flow * (actual_reflux_ratio + 1.0)
    condenser_duty = -vapour_flow * AVERAGE_LATENT_HEAT
    feed_total_enthalpy, _ = enthalpy_at(mixture, ideal_gas, feed_t, feed_p, list(feed_z))
    reboiler_duty = -condenser_duty + (feed_n * feed_total_enthalpy) * 0.01

    distillate_moles = [max(0.0, feed_n * feed_z[i] * fractions[i]) for i in range(count)]
    bottoms_moles = [max(0.0, feed_n * feed_z[i] - distillate_moles[i]) for i in range(count)]
    distillate = _product(mixture, ideal_gas, distillate_moles, condenser_pressure, feed_p, feed_t)
    bottoms = _product(mixture, ideal_gas, bottoms_moles, reboiler_pressure, feed_p, feed_t)

    return _States(
        k_values=tuple(k_values),
        alpha=tuple(alpha),
        relative_volatility=relative_volatility,
        q=q_quality,
        theta=theta,
        minimum_stages=minimum_stages,
        fractions=tuple(fractions),
        minimum_reflux_ratio=minimum_reflux_ratio,
        actual_reflux_ratio=actual_reflux_ratio,
        actual_stages=actual_stages,
        feed_tray_number=feed_tray_number,
        condenser_duty=condenser_duty,
        reboiler_duty=reboiler_duty,
        distillate_n=distillate["n"],
        distillate_z=distillate["z"],
        distillate_h=distillate["h"],
        bottoms_n=bottoms["n"],
        bottoms_z=bottoms["z"],
        bottoms_h=bottoms["h"],
    )


def _product(
    mixture: object,
    ideal_gas: object,
    moles: list[float],
    stated: float | None,
    feed_p: float,
    feed_t: float,
) -> dict[str, object]:
    """One product: the feed's fluid at the split's composition and the stated pressure.

    The class sets the pressure only when the stated one is positive, so an unstated or
    non-positive pressure leaves the product at the feed's.
    """
    total = sum(moles)
    if total <= 0.0:
        raise InvalidInputError(
            "light_key_recovery_distillate / heavy_key_recovery_bottoms",
            "this split leaves one product empty, and a stream with no moles has no composition",
        )
    z = [m / total for m in moles]
    p = feed_p if stated is None or stated <= 0.0 else stated
    h, _ = enthalpy_at(mixture, ideal_gas, feed_t, p, z)  # type: ignore[arg-type]
    return {"n": total, "z": tuple(z), "p": p, "h": h}


def _bounded(split: float) -> float:
    """The class's ``boundedSplitFraction``: a finite fraction pulled inside the interval, and
    a non-finite one replaced by ``0.5``."""
    if not math.isfinite(split):
        return 0.5
    return min(max(split, SPLIT_FLOOR), 1.0 - SPLIT_FLOOR)


def _estimate_distillate_fractions(
    alpha: list[float],
    lk: int,
    hk: int,
    relative_volatility: float,
    minimum_stages: float,
    light_key_recovery_distillate: float,
    heavy_key_recovery_bottoms: float,
) -> list[float]:
    """The fraction of each feed component the class sends to the distillate."""
    fractions = []
    for i in range(len(alpha)):
        if i == lk:
            fraction = light_key_recovery_distillate
        elif i == hk:
            fraction = 1.0 - heavy_key_recovery_bottoms
        elif alpha[i] > relative_volatility:
            fraction = 0.999
        elif alpha[i] < 1.0:
            fraction = 0.001
        else:
            fraction = 1.0 / (
                1.0
                + (alpha[i] / relative_volatility) ** (-minimum_stages)
                * (1.0 - light_key_recovery_distillate)
                / light_key_recovery_distillate
            )
        fractions.append(_bounded(fraction))
    return fractions


def _bottoms_to_distillate(z_feed: list[float], fractions: list[float]) -> float:
    """The class's ``computeBottomsToDistillateRatio``, on the same split fractions."""
    total_d = sum(z_feed[i] * fractions[i] for i in range(len(z_feed)))
    total_b = sum(z_feed[i] * (1.0 - fractions[i]) for i in range(len(z_feed)))
    return total_b / total_d if total_d > 1.0e-15 else 1.0


def _underwood_root(alpha: list[float], z_feed: list[float], q: float, light_key: int) -> float:
    """The class's ``solveUnderwood``: bisection on ``[1 + 1e-6, alpha_LK/HK - 1e-6]``.

    The bracket's upper end is the **light key's** relative volatility and not the array's
    maximum, which is a different number whenever a component is more volatile than the light
    key - and in the captured four-component row two of them are.
    """
    low = 1.0 + 1.0e-6
    high = alpha[light_key] - 1.0e-6
    theta = 0.5 * (low + high)
    for _ in range(UNDERWOOD_ITERATIONS):
        mid = 0.5 * (low + high)
        f_mid = _underwood_function(alpha, z_feed, mid, q)
        if abs(f_mid) < UNDERWOOD_TOLERANCE or (high - low) < UNDERWOOD_BRACKET:
            theta = mid
            break
        f_low = _underwood_function(alpha, z_feed, low, q)
        if f_low * f_mid < 0.0:
            high = mid
        else:
            low = mid
        theta = mid
    return theta


def _underwood_function(alpha: list[float], z_feed: list[float], theta: float, q: float) -> float:
    """``sum_i alpha_i z_i / (alpha_i - theta) - (1 - q)``."""
    total = 0.0
    for i in range(len(alpha)):
        if abs(alpha[i] - theta) > 1.0e-10:
            total += alpha[i] * z_feed[i] / (alpha[i] - theta)
    return total - (1.0 - q)


def _index_of(components: list[str], name: str, which: str) -> int:
    """The position of a named component, or a refusal naming what the caller asked for."""
    if name not in components:
        raise InvalidInputError(which, f"{name} is not one of this fluid's components")
    return components.index(name)


def _spec() -> dict[str, object]:
    """The generated spec, imported here so the module has no import-time cycle."""
    from azoth._models_gen import model

    return model("process.shortcut_distillation_column")
