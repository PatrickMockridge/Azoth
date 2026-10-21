"""``eos.interface_properties`` - NeqSim's phase-pair dispatch for interfacial tension.

Spec: ``specs/models/eos/interface_properties.toml``, which records
``InterfaceProperties.setInterfacialTensionModel``'s model-0 arm and
``getSurfaceTension``'s phase-type dispatch.
"""

from __future__ import annotations

from collections.abc import Sequence
from typing import Any

from azoth.core.errors import InvalidInputError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import InterfacePropertiesResult
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning
from azoth.eos.reference.firozabadi_ramley_interface_tension import (
    firozabadi_ramley_interface_tension,
)
from azoth.eos.reference.parachor_mixture_surface_tension import (
    parachor_mixture_surface_tension,
)

CALC_ID = "eos.interface_properties"

#: The class the gas-liquid and gas-aqueous slots install in model 0.
PARACHOR = "ParachorSurfaceTension"
#: The class the liquid-liquid slot installs in model 0, and every slot's Firozabadi arm.
FIROZABADI = "FirozabadiRamleyInterfaceTension"


def interface_properties(
    interface: str,
    selection: str,
    parachors: Sequence[float],
    tc: Sequence[Q],
    rho_1: Q,
    rho_2: Q,
    molar_mass_1: Q,
    molar_mass_2: Q,
    x_1: Sequence[float],
    x_2: Sequence[float],
    T: Q,
) -> InterfacePropertiesResult:
    """The interface tension of a phase pair by NeqSim's dispatch.

    The pair decides the slot: ``gas_oil`` and ``gas_aqueous`` take Parachor in model 0,
    and every other pair - including a reversed oil/gas one - takes Firozabadi-Ramley.
    Models 3 and 5 are aliases of 0; the string selector's two calculators can be asked
    for directly; models 1, 2 and 4 are refused.

    Args:
        interface: ``gas_oil``, ``gas_aqueous`` or ``liquid_liquid``.
        selection: the class's ``interfacialTensionModel``.
        parachors: each component's parachor.
        tc: each component's critical temperature.
        rho_1: the first phase's mass density.
        rho_2: the second phase's mass density.
        molar_mass_1: the first phase's molar mass.
        molar_mass_2: the second phase's molar mass.
        x_1: the first phase's mole fractions.
        x_2: the second phase's mole fractions.
        T: absolute temperature.

    Returns:
        The tension and the class the selected slot installed.

    Raises:
        InvalidInputError: for an unknown pair or selection, or a refused model.
        OutOfRangeError: if an input is outside its declared range.

    Example:
        >>> import azoth
        >>> q = azoth.ureg.Quantity
        >>> r = interface_properties(
        ...     "gas_oil", "model_0", [77.3, 191.7],
        ...     [q(190.56, "K"), q(425.12, "K")],
        ...     q(19.938328330315997, "kg/m**3"), q(549.3851668858096, "kg/m**3"),
        ...     q(0.022959902730870334, "kg/mol"), q(0.054107691623917306, "kg/mol"),
        ...     [0.8356249351028912, 0.1643750648971088],
        ...     [0.09542082642782022, 0.9045791735721797], q(300.0, "K"))
        >>> round(r.sigma.magnitude, 15), r.selected_model
        (0.009424889197268, 'ParachorSurfaceTension')
    """
    from azoth._models_gen import model

    spec: dict[str, Any] = model(CALC_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    if interface not in ("gas_oil", "gas_aqueous", "liquid_liquid"):
        raise InvalidInputError(
            "interface",
            f"`{interface}` is not one of `gas_oil`, `gas_aqueous` or `liquid_liquid`; the "
            f"pair decides which calculator slot runs",
        )

    if selection in ("model_0", "model_3", "model_5"):
        use_firozabadi = interface == "liquid_liquid"
    elif selection == "parachor":
        use_firozabadi = False
    elif selection == "firozabadi_ramley":
        use_firozabadi = True
    elif selection == "model_1":
        raise InvalidInputError(
            "selection",
            "`model_1` installs `GTSurfaceTension` (`eos.gt_surface_tension`) in the "
            "gas-oil slot and `GTSurfaceTensionSimple` in the other two; all three are "
            "ported, and this dispatch does not carry the influence parameters the gas-oil "
            "calculator takes",
        )
    elif selection == "model_2":
        raise InvalidInputError(
            "selection",
            "`model_2` installs `LGTSurfaceTension` in all three slots. That calculator is "
            "`eos.lgt_surface_tension`, and this dispatch does not carry the influence "
            "parameters, reference component and interaction matrix it takes - so the arm is "
            "refused here rather than assumed",
        )
    elif selection == "model_4":
        raise InvalidInputError(
            "selection",
            "`model_4` installs `GTSurfaceTensionSimple` (`eos.gt_surface_tension_simple`), "
            "`LGTSurfaceTension` (`eos.lgt_surface_tension`) and Parachor in its three slots; "
            "all three are ported, and this dispatch does not carry the influence parameters, "
            "reference component and interaction matrix the two profile calculators take",
        )
    else:
        raise InvalidInputError(
            "selection",
            f"`{selection}` is not one of `model_0`, `model_3`, `model_5`, `parachor`, "
            f"`firozabadi_ramley`, or the refused `model_1`, `model_2`, `model_4`",
        )

    count = len(x_1)
    if len(parachors) != count or len(tc) != count or len(x_2) != count:
        raise InvalidInputError(
            "parachors",
            f"{len(parachors)} parachor(s), {len(tc)} critical temperature(s), {count} and "
            f"{len(x_2)} mole fraction(s): the two phases share one component list, so every "
            f"vector has one entry per component",
        )

    temperature = input_to_si(spec, "T", T)
    density_1 = input_to_si(spec, "rho_1", rho_1)
    density_2 = input_to_si(spec, "rho_2", rho_2)
    mass_1 = input_to_si(spec, "molar_mass_1", molar_mass_1)
    mass_2 = input_to_si(spec, "molar_mass_2", molar_mass_2)
    critical_temperatures = [input_to_si(spec, "tc", value) for value in tc]

    apply_checks(
        checks.on_input,
        {
            "T": temperature,
            "rho_1": density_1,
            "rho_2": density_2,
            "molar_mass_1": mass_1,
            "molar_mass_2": mass_2,
            "x_1": min(float(value) for value in x_1) if x_1 else None,
            "x_2": min(float(value) for value in x_2) if x_2 else None,
        }.get,
        warnings,
    )

    if use_firozabadi:
        pseudo_critical = sum(float(x_1[i]) * critical_temperatures[i] for i in range(count))
        firozabadi = firozabadi_ramley_interface_tension(
            from_si(temperature, "K"),
            from_si(pseudo_critical, "K"),
            from_si(density_1, "kg/m**3"),
            from_si(density_2, "kg/m**3"),
        )
        sigma = firozabadi.sigma
        warnings.extend(firozabadi.warnings)
        selected = FIROZABADI
    else:
        parachor = parachor_mixture_surface_tension(
            parachors,
            from_si(density_1, "kg/m**3"),
            from_si(mass_1, "kg/mol"),
            x_1,
            from_si(density_2, "kg/m**3"),
            from_si(mass_2, "kg/mol"),
            x_2,
        )
        sigma = parachor.sigma
        warnings.extend(parachor.warnings)
        selected = PARACHOR

    return InterfacePropertiesResult(
        sigma=sigma,
        selected_model=selected,
        warnings=tuple(warnings),
    )
