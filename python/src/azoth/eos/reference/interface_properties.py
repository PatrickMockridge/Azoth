"""``eos.interface_properties`` - NeqSim's phase-pair dispatch for interfacial tension.

Spec: ``specs/models/eos/interface_properties.toml``, which records
``InterfaceProperties.setInterfacialTensionModel``'s six arms, the string selector and
``getSurfaceTension``'s phase-type dispatch.
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from typing import Any

from azoth.core.errors import InvalidInputError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import InterfacePropertiesResult
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning
from azoth.eos.reference._mixture_state import reduced_parameters
from azoth.eos.reference.cdft_surface_tension import cdft_surface_tension
from azoth.eos.reference.firozabadi_ramley_interface_tension import (
    firozabadi_ramley_interface_tension,
)
from azoth.eos.reference.gt_surface_tension import gt_surface_tension
from azoth.eos.reference.gt_surface_tension_simple import gt_surface_tension_simple
from azoth.eos.reference.hydrate_inhibitor_wt import AQUEOUS, GAS, OIL, phase_label
from azoth.eos.reference.lgt_surface_tension import lgt_surface_tension
from azoth.eos.reference.parachor_mixture_surface_tension import (
    parachor_mixture_surface_tension,
)
from azoth.eos.reference.tp_multiflash import tp_multiflash

CALC_ID = "eos.interface_properties"

#: The class the gas-liquid and gas-aqueous slots install in model 0.
PARACHOR = "ParachorSurfaceTension"
#: The class the liquid-liquid slot installs in model 0, and every slot's Firozabadi arm.
FIROZABADI = "FirozabadiRamleyInterfaceTension"


class _Phase:
    """One phase, with the two volumes the class reads."""

    def __init__(
        self,
        composition: list[float],
        z_factor: float,
        molar_volume: float,
        translated_density: float,
        molar_mass: float,
    ) -> None:
        self.x = composition
        self.z_factor = z_factor
        self.molar_volume = molar_volume
        self.translated_density = translated_density
        self.molar_mass = molar_mass

    def mass_density(self) -> Q:
        return from_si(self.translated_density * self.molar_mass, "kg/m**3")


def interface_properties(
    components: Sequence[str],
    eos: str,
    selection: str,
    influence_parameters: Sequence[float],
    reference_component: float,
    pair_interactions: Sequence[Sequence[float]],
    T: Q,
    P: Q,
    z: Sequence[float],
) -> InterfacePropertiesResult:
    """The interface tension of a flashed phase pair by NeqSim's dispatch.

    The pair decides the slot: GAS/OIL takes the gas-liquid slot, GAS/AQUEOUS the
    gas-aqueous slot, and every other pair the liquid-liquid slot. Models 3 and 5 are
    aliases of 0; the string selector's three calculators can be asked for directly;
    models 1, 2 and 4 install their own calculators.

    The inputs are the union of the arms': the feed, the fluid, the influence
    parameters, the reference component and the interaction matrix. An arm reads only
    what its calculator uses.

    Raises:
        InvalidInputError: for a feed that is not two phases, an unknown selection, a
            component with no parachor or molar mass, or vectors that are not one entry
            per component.
        OutOfRangeError: if an input is outside its declared range.

    The Rust half is ``azoth_eos::interface_properties``; the two run the same dispatch.
    """
    from azoth._models_gen import model
    from azoth.eos import components as databank

    spec: dict[str, Any] = model(CALC_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    count = len(components)
    if (
        len(influence_parameters) != count
        or len(pair_interactions) != count
        or any(len(row) != count for row in pair_interactions)
        or len(z) != count
    ):
        raise InvalidInputError(
            "influence_parameters",
            f"{count} component(s) against {len(influence_parameters)} influence "
            f"parameter(s), {len(pair_interactions)} interaction row(s) and {len(z)} "
            f"mole fraction(s): every one is one entry per component",
        )
    if (
        not math.isfinite(reference_component)
        or reference_component != int(reference_component)
        or int(reference_component) < 0
        or int(reference_component) >= count
    ):
        raise InvalidInputError(
            "reference_component",
            f"`{reference_component}` is not a whole index in 0..{count}",
        )

    temperature = input_to_si(spec, "T", T)
    pressure = input_to_si(spec, "P", P)
    influences = [float(value) for value in influence_parameters]
    apply_checks(
        checks.on_input,
        {
            "T": temperature,
            "P": pressure,
            "influence_parameters": min(influences) if influences else None,
            "z": min(float(value) for value in z) if z else None,
        }.get,
        warnings,
    )

    mixture, _ = databank.mixture_of(list(components), eos=eos)
    flash = tp_multiflash(
        mixture, from_si(temperature, "K"), from_si(pressure, "Pa"), [float(v) for v in z]
    )
    if flash.phase_count != 2:
        raise InvalidInputError(
            "z",
            f"the feed splits into {flash.phase_count} phase(s) at this state, so there is "
            f"no pair to dispatch; the class's `getSurfaceTension` is called with two phases",
        )

    gases = mixture.cubic.name
    if gases == "srk":
        shift = [component.volume_shift_srk for component in mixture.components]
    elif gases == "pr":
        shift = [component.volume_shift for component in mixture.components]
    else:
        shift = [0.0] * count

    gas_constant = 8.3144621
    rt = gas_constant * temperature
    phases: list[_Phase] = []
    for index in range(2):
        composition = [float(value) for value in flash.x[index]]
        molar_mass = 0.0
        for i, fraction in enumerate(composition):
            mass = mixture.components[i].molar_mass
            if mass is None:
                raise InvalidInputError(
                    "components",
                    f"component {i} carries no molar mass, and the mass densities the "
                    f"parachor and Firozabadi slots take need one",
                )
            molar_mass += fraction * float(mass.to("kg/mol").magnitude)
        molar_volume = flash.z_factor[index] * rt / pressure
        translated = sum(fraction * c for fraction, c in zip(composition, shift, strict=True))
        phases.append(
            _Phase(
                composition,
                flash.z_factor[index],
                molar_volume,
                1.0 / (molar_volume - translated),
                molar_mass,
            )
        )
    # NeqSim's system orders its phases lightest first; the flash's order is denser first.
    phases.sort(key=lambda phase: phase.molar_volume, reverse=True)

    reduced = reduced_parameters(mixture, temperature, pressure)
    labels = [phase_label(mixture, reduced, phase.x, phase.z_factor) for phase in phases]
    if labels[0] == GAS and labels[1] == OIL:
        slot = "gas_oil"
    elif labels[0] == GAS and labels[1] == AQUEOUS:
        slot = "gas_aqueous"
    else:
        slot = "liquid_liquid"

    def parachor_arm() -> tuple[float, str, list[Warning]]:
        parachors: list[float] = []
        for index, component in enumerate(mixture.components):
            if component.parachor <= 0.0:
                raise InvalidInputError(
                    "components",
                    f"component {index} carries no parachor, so the slot's sum would count "
                    f"it as zero; the class reads `getParachor()` and would answer a number "
                    f"built on a missing coefficient",
                )
            parachors.append(component.parachor)
        result = parachor_mixture_surface_tension(
            parachors,
            phases[0].mass_density(),
            from_si(phases[0].molar_mass, "kg/mol"),
            phases[0].x,
            phases[1].mass_density(),
            from_si(phases[1].molar_mass, "kg/mol"),
            phases[1].x,
        )
        return (
            result.sigma.to("N/m").magnitude,
            PARACHOR,
            list(result.warnings),
        )

    def firozabadi_arm() -> tuple[float, str, list[Warning]]:
        pseudo = sum(
            fraction * float(mixture.components[i].Tc.to("K").magnitude)
            for i, fraction in enumerate(phases[0].x)
        )
        result = firozabadi_ramley_interface_tension(
            from_si(temperature, "K"),
            from_si(pseudo, "K"),
            phases[0].mass_density(),
            phases[1].mass_density(),
        )
        return (
            result.sigma.to("N/m").magnitude,
            FIROZABADI,
            list(result.warnings),
        )

    def simple_arm() -> tuple[float, str, list[Warning]]:
        result = gt_surface_tension_simple(
            list(components),
            eos,
            influences,
            reference_component,
            pair_interactions,
            from_si(temperature, "K"),
            from_si(pressure, "Pa"),
            [float(v) for v in z],
        )
        return (
            result.sigma.to("N/m").magnitude,
            "GTSurfaceTensionSimple",
            list(result.warnings),
        )

    def lgt_arm() -> tuple[float, str, list[Warning]]:
        result = lgt_surface_tension(
            list(components),
            eos,
            influences,
            reference_component,
            pair_interactions,
            from_si(temperature, "K"),
            from_si(pressure, "Pa"),
            [float(v) for v in z],
        )
        return (result.sigma.to("N/m").magnitude, "LGTSurfaceTension", list(result.warnings))

    if selection in ("model_0", "model_3", "model_5"):
        arm = firozabadi_arm() if slot == "liquid_liquid" else parachor_arm()
    elif selection == "parachor":
        arm = parachor_arm()
    elif selection == "firozabadi_ramley":
        arm = firozabadi_arm()
    elif selection == "model_1":
        if slot == "gas_oil":
            result = gt_surface_tension(
                list(components),
                eos,
                influences,
                from_si(temperature, "K"),
                from_si(pressure, "Pa"),
                [float(v) for v in z],
            )
            arm = (
                result.sigma.to("N/m").magnitude,
                "GTSurfaceTension",
                list(result.warnings),
            )
        else:
            arm = simple_arm()
    elif selection == "model_2":
        arm = lgt_arm()
    elif selection == "model_4":
        if slot == "gas_oil":
            arm = simple_arm()
        elif slot == "gas_aqueous":
            arm = lgt_arm()
        else:
            arm = parachor_arm()
    elif selection == "cdft":
        cdft_result = cdft_surface_tension(
            list(components),
            eos,
            from_si(temperature, "K"),
            phases[0].x,
            phases[1].x,
            from_si(phases[0].molar_volume, "m**3/mol"),
            from_si(phases[1].molar_volume, "m**3/mol"),
        )
        arm = (
            cdft_result.sigma.to("N/m").magnitude,
            "CDFTSurfaceTension",
            list(cdft_result.warnings),
        )
    else:
        raise InvalidInputError(
            "selection",
            f"`{selection}` is not one of `model_0`, `model_3`, `model_5`, `parachor`, "
            f"`firozabadi_ramley`, `model_1`, `model_2`, `model_4` or `cdft`",
        )

    sigma, selected_model, arm_warnings = arm
    warnings.extend(arm_warnings)
    return InterfacePropertiesResult(
        sigma=from_si(sigma, "N/m"),
        selected_model=selected_model,
        warnings=tuple(warnings),
    )
