"""``DistillationColumn``'s two capacity-limit families, as the reference computes them.

The Python twin of ``crates/azoth-process/src/column/capacity.rs``. Every step mirrors the Rust
line for line, and the two are compared case by case by ``python/tests/test_cross_impl.py``.

**Two families, and they are not one.** ``getFsFactor`` is ``u*sqrt(rho_g)`` over the **total**
cross-section ``pi*D**2/4``; ``getGasLoadFactor`` is the Souders-Brown
``Ks = u*sqrt(rho_g/(rho_l - rho_g))`` over the same area. Both were measured before either was
ported, in ``validation/neqsim/captures/process_column_capacity.tsv``, and three of that
capture's readings shape this module:

**The two read different densities.** The Fs family takes the *system's* ``getDensity("kg/m3")``
and the K family ``getPhase(0)``'s, so an over-cooled overhead with two phases has them a per
cent apart - ``15.34351809893286`` against ``14.4273646532561`` on the capture's two-phase row.

**The K family's liquid density is usually not the liquid's.** ``getPhase(0)`` of NeqSim's liquid
outlet is the *gas* that outlet carries, gas-first ordering being the class's own, so
``rho_l - rho_g`` falls under its ``10.0`` floor and ``1000.0`` is substituted. On the capture's
solved absorber that is the ordinary path: the fallback's K is ``0.006694885232955643`` against
the real liquid's ``0.00826846105398939``.

**Both velocities are superficial over the total area**, which is what separates them from
``TrayHydraulicsCalculator.getFsFactor()``, whose area is the net one less a downcomer.
"""

from __future__ import annotations

import math
from typing import Any

from azoth.core.units import Q, quantity
from azoth.eos import components as _components
from azoth.eos.reference._mixture_state import reduced_parameters
from azoth.eos.reference.hydrate_inhibitor_wt import GAS, phase_label
from azoth.eos.reference.pr_mass_density import pr_mass_density
from azoth.eos.reference.pr_molar_volume import pr_molar_volume
from azoth.eos.reference.pt_flash import pt_flash
from azoth.process.kernels import Stream

#: ``DistillationColumn``'s own constructor default for ``internalDiameter``.
DEFAULT_INTERNAL_DIAMETER_M = 1.0
#: ``DistillationColumn``'s own default for ``maxAllowableFsFactor``, from its constructor.
DEFAULT_MAX_ALLOWABLE_FS_FACTOR = 2.5
#: ``AbsorptionColumn``'s, set in its constructor - so the absorber pair reads ``3.0`` where the
#: distillation column and ``PackedColumn`` read ``2.5``.
DEFAULT_MAX_ALLOWABLE_FS_FACTOR_ABSORBER = 3.0
#: ``AbsorptionColumn.DEFAULT_MAX_ALLOWABLE_GAS_LOAD_FACTOR``.
DEFAULT_MAX_ALLOWABLE_GAS_LOAD_FACTOR = 0.15
#: ``AbsorptionColumn.MIN_LIQUID_GAS_DENSITY_DIFFERENCE``.
MIN_LIQUID_GAS_DENSITY_DIFFERENCE = 10.0
#: ``AbsorptionColumn.DEFAULT_LIQUID_DENSITY``, which that floor substitutes.
DEFAULT_LIQUID_DENSITY = 1000.0


def _molar_mass(mixture: Any, z: list[float]) -> float:
    """The phase's molar mass, kg/mol - ``PhaseInterface.getMolarMass``."""
    return float(
        sum(
            (component.molar_mass.to("kg/mol").magnitude if component.molar_mass else 0.0)
            * fraction
            for component, fraction in zip(mixture.components, z, strict=True)
        )
    )


def _mass_density(mixture: Any, t: float, p: float, z: list[float], root: float) -> float:
    """``getDensity("kg/m3")``: the cubic's volume at the phase's own root, with the shift."""
    molar_mass = _molar_mass(mixture, z)
    volume = pr_molar_volume(root, quantity(t, "K"), quantity(p, "Pa")).v.to("m**3/mol").magnitude
    shift = mixture.volume_shift(z)
    return float(
        pr_mass_density(quantity(molar_mass, "kg/mol"), quantity(volume - shift, "m**3/mol"))
        .rho.to("kg/m**3")
        .magnitude
    )


def _labelled_phases(
    components: list[str], t: float, p: float, z: list[float]
) -> tuple[Any, list[tuple[str, list[float], float]], list[float]]:
    """The phases, labelled and in the order the class's own system carries them: gas first."""
    mixture, _ideal_gas = _components.mixture_of(components, eos="pr")
    flash = pt_flash(mixture, quantity(t, "K"), quantity(p, "Pa"), list(z))
    reduced = reduced_parameters(mixture, t, p)

    if flash.phase == "all_vapour":
        labelled = [(GAS, list(z), flash.z_vapour)]
    elif flash.phase in ("all_liquid", "trivial"):
        labelled = [
            (phase_label(mixture, reduced, list(z), flash.z_liquid), list(z), flash.z_liquid)
        ]
    else:
        labelled = [
            (GAS, list(flash.y), flash.z_vapour),
            (
                phase_label(mixture, reduced, list(flash.x), flash.z_liquid),
                list(flash.x),
                flash.z_liquid,
            ),
        ]

    shares: list[float] = []
    for label, _phase_z, _root in labelled:
        if len(labelled) == 1:
            shares.append(1.0)
        elif label == GAS:
            shares.append(flash.vapour_fraction if flash.vapour_fraction is not None else 1.0)
        else:
            shares.append(
                1.0 - (flash.vapour_fraction if flash.vapour_fraction is not None else 0.0)
            )
    return mixture, labelled, shares


def _system_per_mole(
    components: list[str], t: float, p: float, z: list[float]
) -> tuple[float, float]:
    """The whole system's mass and volume per mole of it: ``(kg/mol, m**3/mol)``.

    **The sum over the phases, and every capacity-limit reading needs it.** ``getDensity("kg/m3")``
    and ``getFlowRate("m3/sec")`` are *system* quantities in NeqSim, not phase ones, and on a
    two-phase outlet neither is either phase's.
    """
    mixture, labelled, shares = _labelled_phases(components, t, p, z)
    mass = 0.0
    volume = 0.0
    for (_label, phase_z, root), share in zip(labelled, shares, strict=True):
        molar_mass = _molar_mass(mixture, phase_z)
        density = _mass_density(mixture, t, p, phase_z, root)
        mass += share * molar_mass
        volume += share * molar_mass / density
    return mass, volume


def _phase_zero_density(components: list[str], t: float, p: float, z: list[float]) -> float:
    """``getPhase(0)``'s mass density: the gas phase where there is one, the single phase else."""
    mixture, labelled, _shares = _labelled_phases(components, t, p, z)
    _label, phase_z, root = labelled[0]
    return _mass_density(mixture, t, p, phase_z, root)


def resolve_liquid_density(gas_density: float, liquid_density_field: float) -> float:
    """``resolveLiquidDensityForGasLoad``: the near-dry fallback, and its ``NaN``.

    **``NaN`` propagates rather than substituting.** A field the difference cannot be taken
    against leaves the comparison false in both branches, so the caller's ``isnan`` test is what
    refuses it - which is why this returns ``NaN`` rather than a zero the caller would have to
    distinguish.
    """
    if liquid_density_field - gas_density < MIN_LIQUID_GAS_DENSITY_DIFFERENCE:
        density = DEFAULT_LIQUID_DENSITY
    else:
        density = liquid_density_field
    return float("nan") if density <= gas_density else density


def fs_limits(
    gas_out: Stream,
    internal_diameter: Q,
    max_allowable_fs_factor: float,
) -> dict[str, Any]:
    """``getFsFactor`` and its three siblings, at the products of a solved column.

    **The Fs family is the base class's**, so every column has it - and ``PackedColumn``'s
    diameter is the one its own sizing resolved rather than the one a caller stated.
    """
    components = list(gas_out.components)
    diameter = internal_diameter.to("m").magnitude
    area = math.pi * diameter * diameter / 4.0

    gas_t = gas_out.t.to("K").magnitude
    gas_p = gas_out.p.to("Pa").magnitude
    gas_z = list(gas_out.z)
    gas_n = gas_out.n

    # ``getFsFactor``: the area guard comes first and answers zero.
    if area <= 0.0:
        fs_factor = 0.0
    else:
        mass, volume_per_mole = _system_per_mole(components, gas_t, gas_p, gas_z)
        fs_factor = (gas_n * volume_per_mole) / area * math.sqrt(mass / volume_per_mole)

    # ``getMinimumDiameterForFsLimit`` **does not consult the area at all** - its own guard is the
    # limit's - so a column with no diameter still owes this number.
    if max_allowable_fs_factor > 0.0:
        mass, volume_per_mole = _system_per_mole(components, gas_t, gas_p, gas_z)
        minimum_fs = math.sqrt(
            4.0
            * (gas_n * volume_per_mole)
            * math.sqrt(mass / volume_per_mole)
            / (math.pi * max_allowable_fs_factor)
        )
    else:
        minimum_fs = 0.0

    return {
        "fs_factor": fs_factor,
        "fs_factor_utilization": _utilization(fs_factor, max_allowable_fs_factor),
        "fs_factor_within_design_limit": fs_factor <= max_allowable_fs_factor,
        "minimum_diameter_for_fs_limit": minimum_fs,
    }


def gas_load_limits(
    gas_out: Stream,
    liquid_out_vapour: Stream,
    internal_diameter: Q,
    max_allowable_gas_load_factor: float,
) -> dict[str, Any]:
    """``getGasLoadFactor`` and its three siblings, on the absorber pair.

    **``getPhase(0)`` of both outlets, and the second one is not the liquid.** A NeqSim
    ``getLiquidOutStream`` answers a stream whose thermo system is the *tray's*, so its phase 0 -
    the array is gas-first - is the **vapour that tray carries**. ``liquid_out_vapour`` is
    therefore the tray's vapour and not the liquid product: on the capture's stripper the field
    density is ``9.221459471348782`` against a gas of ``11.267630874486859``, so
    ``rho_l - rho_g`` is negative and the ``10.0`` floor is what answers - while the liquid's own
    density would give ``0.0011870273688000244`` against the class's ``0.0009402581113639056``, a
    factor of ``1.26``.
    """
    components = list(gas_out.components)
    diameter = internal_diameter.to("m").magnitude
    area = math.pi * diameter * diameter / 4.0

    gas_t = gas_out.t.to("K").magnitude
    gas_p = gas_out.p.to("Pa").magnitude
    gas_z = list(gas_out.z)
    gas_n = gas_out.n

    # ``getGasSuperficialVelocity``: the same total area, and the same guard as the Fs family's.
    if area <= 0.0:
        velocity = 0.0
    else:
        _mass, volume_per_mole = _system_per_mole(components, gas_t, gas_p, gas_z)
        velocity = (gas_n * volume_per_mole) / area

    if _not_positive(velocity):
        gas_load_factor = 0.0
    else:
        gas_density, liquid_density = _densities(gas_out, liquid_out_vapour)
        if _not_positive(gas_density) or math.isnan(liquid_density):
            gas_load_factor = 0.0
        else:
            gas_load_factor = velocity * math.sqrt(gas_density / (liquid_density - gas_density))

    minimum_gas_load = 0.0
    if max_allowable_gas_load_factor > 0.0:
        gas_density, liquid_density = _densities(gas_out, liquid_out_vapour)
        if not _not_positive(gas_density) and not math.isnan(liquid_density):
            permissible = max_allowable_gas_load_factor * math.sqrt(
                (liquid_density - gas_density) / gas_density
            )
            if permissible > 0.0:
                _mass, volume_per_mole = _system_per_mole(components, gas_t, gas_p, gas_z)
                minimum_gas_load = math.sqrt(
                    4.0 * (gas_n * volume_per_mole) / (math.pi * permissible)
                )

    return {
        "gas_load_factor": gas_load_factor,
        "gas_load_factor_utilization": _utilization(gas_load_factor, max_allowable_gas_load_factor),
        "gas_load_factor_within_design_limit": (gas_load_factor <= max_allowable_gas_load_factor),
        "minimum_diameter_for_gas_load_limit": minimum_gas_load,
    }


def _densities(gas_out: Stream, liquid_out_vapour: Stream) -> tuple[float, float]:
    """The two densities both gas-load getters read: ``getPhase(0)`` of the gas outlet and of the
    liquid outlet's own system, with the near-dry fallback applied to the second."""
    components = list(gas_out.components)
    gas_density = _phase_zero_density(
        components,
        gas_out.t.to("K").magnitude,
        gas_out.p.to("Pa").magnitude,
        list(gas_out.z),
    )
    liquid_density = resolve_liquid_density(
        gas_density,
        _phase_zero_density(
            components,
            liquid_out_vapour.t.to("K").magnitude,
            liquid_out_vapour.p.to("Pa").magnitude,
            list(liquid_out_vapour.z),
        ),
    )
    return gas_density, liquid_density


def _not_positive(value: float) -> bool:
    """The class's own ``!(x > 0.0)`` guard, which a ``NaN`` also fails."""
    return math.isnan(value) or value <= 0.0


def _utilization(value: float, limit: float) -> float:
    """The two utilizations, which share one shape: a limit that is not positive answers zero."""
    return value / limit if limit > 0.0 else 0.0
