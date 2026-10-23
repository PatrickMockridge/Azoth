"""A :class:`~azoth.eos.mixture.Mixture` rebuilt from a characterization table.

Mirrors ``crates/azoth-characterization/src/assembly.rs``. This is the step NeqSim's class does
with ``addTBPfraction`` + ``profile.applyTo``, and it is a library helper rather than a
calculation: its answer is a fluid, and a fluid is not a result.

A base row is a name: its ``Tc``, ``Pc``, ``omega``, interaction pairs and heat-capacity
polynomial are the databank's, and a row that carries properties of its own has them ignored. A
pseudo row is the opposite and must carry every number a component needs. The Peneloux
translation is derived here from the row's Rackett compressibility, exactly as the databank path
derives it - it is a function of ``racket_z`` and not a field.

``ideal_gas_enthalpy_of_formation`` is read by the blend and has no home on a ``Component`` or an
:class:`~azoth.eos.IdealGasModel`, so it is dropped at assembly.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence
from dataclasses import replace
from typing import Any, cast

import pint

from azoth.core.errors import InvalidInputError, PropertyUnavailableError
from azoth.core.units import Q, ureg
from azoth.eos import (
    Component,
    IdealGasModel,
    mixture,
    pr_peneloux_shift,
    srk_kappa,
    srk_peneloux_shift,
)
from azoth.eos.components import entry, kij_for
from azoth.eos.cubic import CUBICS

from .reference._combiner import CP_D_GETTER_OFFSET

#: The unit each dimensioned row field is stated in when it arrives as a bare number.
UNITS = {
    "molar_mass": "kg/mol",
    "density": "kg/m**3",
    "normal_boiling_point": "K",
    "critical_temperature": "K",
    "critical_pressure": "Pa",
    "critical_volume": "m**3/mol",
    "critical_viscosity": "Pa*s",
    "triple_point_temperature": "K",
    "heat_of_fusion": "J/mol",
    "ideal_gas_enthalpy_of_formation": "J/mol",
    "cp_a": "J/(mol*K)",
    "cp_b": "J/(mol*K**2)",
    "cp_c": "J/(mol*K**3)",
    "cp_d": "J/(mol*K**4)",
}


def _field(row: Any, name: str) -> Any:
    """One field of a row, from a mapping or from a result row's attributes."""
    if isinstance(row, Mapping):
        return row.get(name)
    return getattr(row, name, None)


def _q(row: Any, name: str) -> Q | None:
    """One row field as a quantity, or absent. A bare number is the field's declared unit."""
    value = _field(row, name)
    if value is None:
        return None
    if isinstance(value, pint.Quantity):
        return value
    return cast("Q", ureg.Quantity(value, UNITS[name]))


def _scalar(value: Any) -> float:
    """One row field as a plain number, converting a quantity through its SI magnitude."""
    if isinstance(value, pint.Quantity):
        return float(value.to_base_units().magnitude)
    return float(value)


def _required(row: Any, name: str, index: int) -> Any:
    value = _field(row, name)
    if value is None:
        raise InvalidInputError(
            f"components[{index}].{name}",
            f"a pseudo-component needs its own {name}, and this row states none",
        )
    return value


def _quantity(row: Any, name: str, index: int) -> Q:
    """One required row field as a quantity. A bare number is in the declared unit."""
    value = _required(row, name, index)
    return (
        value if isinstance(value, pint.Quantity) else cast("Q", ureg.Quantity(value, UNITS[name]))
    )


def _number(row: Any, name: str, index: int) -> float:
    """One required row field as an SI magnitude."""
    return float(_quantity(row, name, index).to_base_units().magnitude)


def _pseudo_component(row: Any, index: int, cubic: str) -> Component:
    """One pseudo-component as a ``Component``, its fitted alpha coefficient and shifts."""
    tc = _quantity(row, "critical_temperature", index)
    pc = _quantity(row, "critical_pressure", index)
    omega = _scalar(_required(row, "acentric_factor", index))
    molar_mass = _quantity(row, "molar_mass", index)
    rackett_z = _field(row, "racket_z")
    rackett_z = None if rackett_z is None else _scalar(rackett_z)
    pr_shift = pr_peneloux_shift(omega, tc, pc, rackett_z)
    srk_shift = srk_peneloux_shift(omega, tc, pc, rackett_z)
    parachor = _field(row, "parachor")
    heat_of_fusion = _field(row, "heat_of_fusion")
    triple_point = _field(row, "triple_point_temperature")
    params: tuple[float, ...] = ()
    if cubic == "srk":
        # NeqSim's `addTBPfraction` installs the fitted `m`; a row that states none is a cut
        # whose alpha is Soave's own.
        attractive = _field(row, "attractive_m")
        params = (float(attractive) if attractive is not None else srk_kappa(omega).kappa,)
    return Component(
        Tc=tc,
        Pc=pc,
        omega=omega,
        rackett_z=0.0 if rackett_z is None else rackett_z,
        volume_shift=pr_shift.c.to("m**3/mol").magnitude,
        volume_shift_srk=srk_shift.c.to("m**3/mol").magnitude,
        molar_mass=molar_mass,
        alpha_params=params,
        component_class="hc",
        parachor=0.0 if parachor is None else _scalar(parachor),
        heat_of_fusion=0.0 if heat_of_fusion is None else _scalar(heat_of_fusion),
        triple_point_temperature=0.0 if triple_point is None else _scalar(triple_point),
    )


def _pseudo_cp(row: Any, index: int) -> tuple[float, float, float, float, float]:
    """A pseudo-component's heat-capacity polynomial, all four coefficients required."""
    missing = [name for name in ("cp_a", "cp_b", "cp_c", "cp_d") if _field(row, name) is None]
    if missing:
        raise InvalidInputError(
            f"components[{index}].{missing[0]}",
            f"a pseudo-component needs its own {missing[0]}; a mixture's enthalpy needs the "
            "polynomial rather than a row of zeros",
        )
    cp_d = _number(row, "cp_d", index) - CP_D_GETTER_OFFSET
    # NeqSim's TBP path never sets `CpE`, so the quartic and higher term is absent.
    return (
        _number(row, "cp_a", index),
        _number(row, "cp_b", index),
        _number(row, "cp_c", index),
        cp_d,
        0.0,
    )


def mixture_from_characterization(
    components: Iterable[Any], *, eos: str = "srk"
) -> tuple[Any, IdealGasModel]:
    """A mixture and its ideal-gas model, from a characterization table.

    ``components`` are rows in the order the mixture should carry them - the result of
    :func:`azoth.characterization.combine_reservoir_fluids` is the usual source, and a mapping
    with the same field names is accepted too. A bare number is stated in the unit the spec
    declares for that field; a ``pint`` quantity is converted from whatever unit it carries.

    ``eos`` is ``"srk"``, the class's own fluid, or ``"pr"``; Peng-Robinson has no fitted ``m``
    and ignores a row's ``attractive_m``.

    Raises:
        PropertyUnavailableError: if a base row's name is not in the databank, or its substance
            has no heat-capacity polynomial.
        InvalidInputError: if a pseudo row omits a critical constant, molar mass or
            heat-capacity coefficient, or ``eos`` is not a cubic this builds.
    """
    if eos not in ("pr", "srk"):
        raise InvalidInputError(
            "eos", f"{eos!r} is not a cubic a characterized fluid is built for; use 'pr' or 'srk'"
        )

    rows = list(components)
    built: list[Component] = []
    names: list[str] = []
    base_names: list[str | None] = []
    cp: list[tuple[float, float, float, float, float]] = []
    for index, row in enumerate(rows):
        name = str(_required(row, "name", index))
        names.append(name)
        if bool(_field(row, "pseudo_component")):
            built.append(_pseudo_component(row, index, eos))
            cp.append(_pseudo_cp(row, index))
            base_names.append(None)
        else:
            record = entry(name)
            if record.component_type == "ion":
                raise InvalidInputError(
                    "components",
                    f"{name} is an ion, and a cubic has no notion of one; an ion belongs to "
                    "the electrolyte models",
                )
            if record.cp is None:
                raise PropertyUnavailableError(
                    name,
                    "heat-capacity coefficients",
                    "a mixture's ideal-gas model needs the polynomial, and this substance has "
                    "none - which is what a keycard-added substance carries",
                )
            component = record.component(alpha="srk" if eos == "srk" else None)
            if eos == "srk":
                component = replace(
                    component,
                    alpha_params=(srk_kappa(component.omega).kappa,),
                )
            built.append(component)
            cp.append(record.cp)
            base_names.append(record.name)

    # `kij_for` is keyed by index into the names it is given, so the pairs are remapped onto
    # the mixture's own order - a pseudo row between two base rows moves every index after it.
    known = [name for name in base_names if name is not None]
    at = {name: index for index, name in enumerate(base_names) if name is not None}
    pairs = {
        (at[known[first]], at[known[second]]): value
        for (first, second), value in kij_for(tuple(known), eos=eos).items()
    }
    fluid = mixture(
        built,
        kij=pairs,
        cubic=CUBICS[eos],
        alpha="srk_fitted" if eos == "srk" else "pr",
        names=names,
    )
    ideal_gas = IdealGasModel(
        cp_a=tuple(row[0] for row in cp),
        cp_b=tuple(row[1] for row in cp),
        cp_c=tuple(row[2] for row in cp),
        cp_d=tuple(row[3] for row in cp),
        cp_e=tuple(row[4] for row in cp),
    )
    return fluid, ideal_gas


def oil_assay_to_mixture(
    basis: str,
    cuts: Sequence[Mapping[str, Any]],
    total_assay_mass: Q,
    *,
    eos: str = "srk",
) -> tuple[Any, IdealGasModel]:
    """An oil assay as the mixture its cuts describe - ``OilAssayCharacterisation.apply``.

    ``basis``, ``cuts`` and ``total_assay_mass`` are
    :func:`azoth.characterization.oil_assay_derivation`'s own: each cut is a mapping that states
    its yield and the routes to its molar mass and gravity. The derivation resolves each cut's
    molar mass and mole amount, the cut's card comes from
    :func:`azoth.characterization.tbp_fraction_component`, and the mixture is built by
    :func:`mixture_from_characterization`.

    A cut the assay gives as a standard substance is a base component - its constants are the
    databank's - and no card is built for it. The density a pseudo-component's card is built from
    is the one the derivation resolved, stated or implied by the cut's API gravity.

    **The fluid is constructed rather than mutated**, which is the one place this differs from
    the class's ``addTBPfraction``/``addComponent`` sequence.
    """
    from azoth import _models_gen
    from azoth.characterization.reference.oil_assay_derivation import (
        oil_assay_derivation,
        specific_gravity,
    )
    from azoth.characterization.reference.tbp_fraction_component import (
        tbp_fraction_component as card_of,
    )
    from azoth.core.units import input_records_to_si

    spec = _models_gen.model("characterization.oil_assay_derivation")
    converted: list[dict[str, Any]] = input_records_to_si(spec, "cuts", cuts)
    derived = oil_assay_derivation(basis, cuts, total_assay_mass)

    keys = [
        (
            str(row["standard_component_name"])
            if row["standard_component_name"] is not None
            else f"{row['name']}_PC"
        ).lower()
        for row in converted
    ]
    if len(set(keys)) != len(keys):
        raise InvalidInputError(
            "cuts", "two cuts name one component once the pseudo-component suffix is added"
        )

    rows: list[dict[str, Any]] = []
    for index, name in enumerate(derived.component_name):
        key = str(name).lower()
        if key not in keys:
            raise InvalidInputError(
                "cuts", f"the derivation answered `{name}`, which no cut describes"
            )
        cut = converted[keys.index(key)]
        moles = derived.moles[index]
        if cut["standard_component_name"] is not None:
            rows.append(
                {
                    "fluid": 0.0,
                    "name": str(name),
                    "pseudo_component": False,
                    "moles": moles,
                }
            )
            continue
        gravity = specific_gravity(cut)
        if gravity is None:
            raise InvalidInputError(
                "cuts", f"`{cut['name']}` states neither a density nor an API gravity"
            )
        card = card_of(derived.molar_mass[index], ureg.Quantity(gravity * 1000.0, "kg/m**3"))
        rows.append(
            {
                "fluid": 0.0,
                "name": str(name),
                "pseudo_component": True,
                "moles": moles,
                "molar_mass": derived.molar_mass[index],
                "density": ureg.Quantity(gravity * 1000.0, "kg/m**3"),
                "normal_boiling_point": card.boiling_temperature,
                "critical_temperature": card.critical_temperature,
                "critical_pressure": card.critical_pressure,
                "acentric_factor": card.acentric_factor,
                "critical_volume": card.critical_volume,
                "racket_z": card.racket_z,
                "racket_z_cpa": card.racket_z_cpa,
                "parachor": card.parachor,
                "critical_viscosity": card.critical_viscosity,
                "triple_point_temperature": card.triple_point_temperature,
                "heat_of_fusion": card.heat_of_fusion,
                "ideal_gas_enthalpy_of_formation": card.ideal_gas_enthalpy_of_formation,
                "cp_a": card.cp_a,
                "cp_b": card.cp_b,
                "cp_c": card.cp_c,
                "cp_d": card.cp_d,
                "attractive_m": card.attraction_exponent,
            }
        )
    return mixture_from_characterization(rows, eos=eos)
