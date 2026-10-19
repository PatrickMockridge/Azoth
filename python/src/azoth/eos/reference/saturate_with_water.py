"""``eos.saturate_with_water`` - the water a fluid holds when it is water-saturated.

Spec: ``specs/models/eos/saturate_with_water.toml``. NeqSim's ``SaturateWithWater`` floods the
feed with water, flashes, and reads the water dissolved in the non-aqueous phases - a fluid
coexisting with a free aqueous phase is, by definition, saturated - then resets the overall
water to that amount and re-flashes.

This is the pure-Python reference: a second, independent expression of the same operation as
the Rust kernel, held to it by ``test_cross_impl.py``.

**The class mutates a system; this is a function.** The three paths that restore the caller's
composition instead of returning a flooded fluid - the flood flash diverged, no aqueous phase
formed, a residual incipient aqueous phase survived the correction - become "return the feed",
with the reason in the result's warnings.
"""

from __future__ import annotations

import dataclasses
import math

from azoth import _models_gen
from azoth.core.errors import InvalidInputError, OutOfRangeError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import SaturateWithWaterResult
from azoth.core.units import Q, input_to_si, ureg
from azoth.core.warnings import Warning, WarningCode
from azoth.eos.components import mixture_of
from azoth.eos.reference._association import RDF_CHAPMAN
from azoth.eos.reference._mixture_state import reduced_parameters
from azoth.eos.reference.hydrate_inhibitor_wt import AQUEOUS, phase_label
from azoth.eos.reference.tp_multiflash import tp_multiflash

MODEL_ID = "eos.saturate_with_water"

#: The class's `addComponent("water", totalMoles / 100.0)`, as a fraction of the feed.
TRACE = 0.01


def _labels(mixture: object, temperature: float, pressure: float, flash: object) -> list[str]:
    """NeqSim's `PhaseEos.init` label for every phase of a flash."""
    reduced = reduced_parameters(mixture, temperature, pressure)  # type: ignore[arg-type]
    return [
        phase_label(mixture, reduced, list(row), z_factor)  # type: ignore[arg-type]
        for row, z_factor in zip(flash.x, flash.z_factor, strict=True)  # type: ignore[attr-defined]
    ]


def _at_the_feed(
    mixture: object,
    moles: list[float],
    total: float,
    water: int,
    temperature: float,
    pressure: float,
    warnings: list[Warning],
) -> SaturateWithWaterResult:
    """The feed's own state, which is what each recovery path returns."""
    fractions = [value / total for value in moles]
    flash = tp_multiflash(
        mixture,  # type: ignore[arg-type]
        ureg.Quantity(temperature, "K"),
        ureg.Quantity(pressure, "Pa"),
        fractions,
    )
    return SaturateWithWaterResult(
        water_mole_fraction=moles[water] / total,
        phase_count=flash.phase_count,
        phase_fractions=tuple(flash.phase_fractions),
        x=tuple(tuple(row) for row in flash.x),
        z_factor=tuple(flash.z_factor),
        aqueous_present=False,
        warnings=tuple(warnings),
    )


def saturate_with_water(
    components: list[str],
    moles: list[float],
    T: Q,
    P: Q,
    eos: str = "srk",
    associating: bool = False,
) -> SaturateWithWaterResult:
    """The water a fluid holds at saturation, with the phase state at that composition.

    Args:
        components: the fluid's components, by name. Water is added at a hundredth of the
            feed when absent, as the class's own first step does.
        moles: the feed's mole numbers.
        T: absolute temperature the flash runs at.
        P: absolute pressure the flash runs at.
        eos: the cubic the fluid is built on.
        associating: whether to build the fluid as a CPA mixture.

    Returns:
        The overall water mole fraction at saturation, the reported phases and their
        compositions, and whether a free aqueous phase was found.

    Raises:
        InvalidInputError: if ``moles`` is not one entry per component or is not a feed.
        PropertyUnavailableError: if a name is not in the databank.

    See :func:`azoth.eos.saturate_with_water`.
    """
    spec = _models_gen.model(MODEL_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    t_si = input_to_si(spec, "T", T)
    p_si = input_to_si(spec, "P", P)
    apply_checks(checks.on_input, {"T": t_si, "P": p_si}.get, warnings)

    feed = [input_to_si(spec, "moles", value) for value in moles]
    if len(components) != len(feed):
        raise InvalidInputError(
            "moles",
            f"{len(components)} components but {len(feed)} mole numbers; the two must match",
        )
    # **A bound on the input, so it is the class the spec's other bounds use** - which is
    # what the Rust entry raises for the same state.
    if any(not math.isfinite(m) or m < 0.0 for m in feed):
        first = next(m for m in feed if not math.isfinite(m) or m < 0.0)
        raise OutOfRangeError("moles", first, "every mole number must be finite and not negative")
    total = sum(feed)
    if not total > 0.0 or not math.isfinite(total):
        raise OutOfRangeError("moles", total, "the feed's mole numbers do not describe a feed")

    # **The trace of water, and it is not an optimisation to skip.** The class adds it before
    # the flood, so the saturation is reported against a total that includes it.
    names = [name.strip().lower() for name in components]
    if "water" not in names:
        names.append("water")
        feed.append(total * TRACE)
    total = sum(feed)
    input_moles = list(feed)

    mixture, _ = mixture_of(names, eos=eos, associating=associating)
    if associating:
        # **The Statoil phase, which is the generic CPA fluid plus its radial distribution
        # function.** The source's CPA cases are `SystemSrkCPAstatoil`, whose
        # `PhaseSrkCPAs.calc_g` replaces `PhaseCPAInterface`'s Carnahan-Starling with
        # `1/(1 - 1.9 eta)`. Building the Carnahan-Starling fluid instead moves the
        # dissolved water by 15% in both source states, in opposite directions; with this
        # and the temperature-dependent CPA `kij` the two reproduce to 0.2-0.3%.
        mixture = dataclasses.replace(mixture, rdf_form=RDF_CHAPMAN)
    water = names.index("water")

    # The flood: one feed's worth of water.
    flooded_total = total * 2.0
    flooded = list(input_moles)
    flooded[water] += total
    flooded_fractions = [value / flooded_total for value in flooded]

    try:
        flood = tp_multiflash(
            mixture,
            ureg.Quantity(t_si, "K"),
            ureg.Quantity(p_si, "Pa"),
            flooded_fractions,
        )
    # `runFlashSafely` recovers from any failure the flash raises, not from a named set.
    except Exception as error:
        warnings.append(
            Warning(
                WarningCode.FALLBACK_TO_INPUT,
                f"the flood flash failed ({error}); the feed is returned unchanged",
            )
        )
        return _at_the_feed(mixture, input_moles, total, water, t_si, p_si, warnings)

    flood_labels = _labels(mixture, t_si, p_si, flood)
    if AQUEOUS not in flood_labels:
        warnings.append(
            Warning(
                WarningCode.FALLBACK_TO_INPUT,
                "flooding produced no free aqueous phase, so saturation is undefined here and "
                "the feed is returned unchanged",
            )
        )
        return _at_the_feed(mixture, input_moles, total, water, t_si, p_si, warnings)
    aqueous = flood_labels.index(AQUEOUS)

    # The water dissolved in every non-aqueous phase, which is the saturation amount.
    saturated = sum(
        beta * row[water] * flooded_total
        for index, (beta, row) in enumerate(zip(flood.phase_fractions, flood.x, strict=True))
        if index != aqueous
    )

    saturated_moles = list(input_moles)
    saturated_moles[water] = saturated
    saturated_total = sum(saturated_moles)
    if not saturated_total > 0.0 or not math.isfinite(saturated_total):
        raise InvalidInputError(
            "moles",
            "the saturated composition has no moles, which the correction produced",
        )
    saturated_fractions = [value / saturated_total for value in saturated_moles]

    saturated_flash = tp_multiflash(
        mixture,
        ureg.Quantity(t_si, "K"),
        ureg.Quantity(p_si, "Pa"),
        saturated_fractions,
    )
    saturated_labels = _labels(mixture, t_si, p_si, saturated_flash)

    # **A residual incipient aqueous phase is dropped, not flashed away.** The class calls
    # `removePhase` and re-flashes; a flash here is a function of the composition alone.
    kept = [
        index
        for index, label in enumerate(saturated_labels)
        if not (label == AQUEOUS and saturated_flash.phase_count > 1)
    ]
    kept_fraction = sum(saturated_flash.phase_fractions[index] for index in kept)

    return SaturateWithWaterResult(
        water_mole_fraction=saturated / saturated_total,
        phase_count=len(kept),
        phase_fractions=tuple(
            saturated_flash.phase_fractions[index] / kept_fraction for index in kept
        ),
        x=tuple(tuple(saturated_flash.x[index]) for index in kept),
        z_factor=tuple(saturated_flash.z_factor[index] for index in kept),
        aqueous_present=True,
        warnings=tuple(warnings),
    )
