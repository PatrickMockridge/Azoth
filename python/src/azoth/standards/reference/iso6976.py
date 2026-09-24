"""``standards.iso6976`` - the calorific values and density of a natural gas.

Spec: ``specs/models/standards/iso6976.toml``

The Python twin of ``crates/azoth-standards/src/iso6976.rs`` over
``crates/azoth-standards/src/iso6976_constants.rs``, written to mirror them.

``M    = sum z_i M_i``, ``Z = 1 - (sum z_i sqrt(b_i))**2``,
``d = (sum z_i M_i) / M_air * Z_air / Z``, ``Hsup = sum z_i Hsup_i``,
``Hinf = sum z_i Hinf_i``, ``rho_ideal = P M / (R T)`` and ``rho_real = rho_ideal / Z``.

# The table, and which revision of it

``data/standards/iso6976.csv``, compiled from NeqSim's ``ISO6976constants.csv`` by
``tools/gen_databank.py``. **NeqSim ships two ISO 6976 tables** - the original and a 2016
revision - and they disagree on the compression factors, the summation factors and the molar
masses. The class ``Stream.LCV()`` builds queries the original, so that is the one compiled;
``Standard_ISO6976_2016``'s table stays vendored and unread, and the manifest says so.

# The two constants

``R`` is ``8.314510`` and the reference pressure ``1.01325`` bar, both as NeqSim's class
declares them rather than as the current CODATA values - the port reproduces the class.

# The reference temperatures cross as kelvin

This library states a temperature in kelvin and the standard states its rows in degrees
Celsius, so the two cross once, here, by the exact ``273.15``. The table's selectors then
match with a tolerance rather than by equality: a caller writes ``288.7`` K for the 60 F
reference and ``288.7 - 273.15`` is ``15.550000000000011``.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from typing import Any

from azoth._data import find
from azoth.core.errors import InvalidInputError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import Iso6976Result
from azoth.core.units import Q, from_si, input_to_si

#: The standard's compiled per-component table.
ISO6976_CSV = "data/standards/iso6976.csv"

#: The class's own molar gas constant, J/(mol*K).
GAS_CONSTANT = 8.314510

#: The standard's reference pressure, Pa: 1.01325 bar.
REFERENCE_PRESSURE = 101_325.0

#: Dry air's molar mass, kg/mol, from `ThermodynamicConstantsInterface.molarMassAir`.
MOLAR_MASS_AIR = 0.02896546

#: 273.15 K, so the standard's degrees Celsius are reachable from a kelvin temperature.
KELVIN_AT_ZERO_CELSIUS = 273.15

#: How close a stated reference temperature has to be to one the table carries.
TOLERANCE = 1.0e-9

#: Air's own compression factor at the standard's volumetric reference temperatures, from
#: the class's field initialisers `Zair0`, `Zair15` and `Zair20`.
AIR_COMPRESSION_FACTOR = ((0.0, 0.99941), (15.0, 0.99958), (15.55, 0.99958), (20.0, 0.99963))


def iso6976(
    components: list[str],
    z: list[float],
    volumetric_reference_temperature: Q,
    energy_reference_temperature: Q,
) -> Iso6976Result:
    """The calorific values and density of a natural gas.

    Args:
        components: the gas's substances, by name. Every one must have a row in the
            standard's table, which is narrower than the component databank.
        z: the composition, in the order the components are listed.
        volumetric_reference_temperature: the temperature the densities and the relative
            density are formed at, **in kelvin**.
        energy_reference_temperature: the temperature the calorific values are stated at, in
            kelvin - the combustion reference.

    Returns:
        The mixture's molar mass, compression factor, relative density, two densities and
        two calorific values.

    Raises:
        InvalidInputError: where the shapes disagree, a component has no row in the
            standard's table, or a reference temperature is not one the table carries.
        OutOfRangeError: for a temperature outside the standard's set.

    Example:
        >>> import azoth
        >>> q = azoth.ureg.Quantity
        >>> r = iso6976(["methane"], [1.0], q(288.15, "K"), q(298.15, "K"))
        >>> round(r.molar_mass.to("kg/mol").magnitude, 6)
        0.016043
    """
    spec = _spec()
    checks = checks_for(spec)

    vol = input_to_si(spec, "volumetric_reference_temperature", volumetric_reference_temperature)
    energy = input_to_si(spec, "energy_reference_temperature", energy_reference_temperature)
    warnings: list[Any] = []
    apply_checks(
        checks.on_input,
        {
            "volumetric_reference_temperature": vol,
            "energy_reference_temperature": energy,
        }.get,
        warnings,
    )

    states = _route(components, z, vol, energy)

    return Iso6976Result(
        molar_mass=from_si(states.molar_mass, "kg/mol"),
        compression_factor=states.compression_factor,
        relative_density=states.relative_density,
        density_ideal=from_si(states.density_ideal, "kg/m**3"),
        density_real=from_si(states.density_real, "kg/m**3"),
        superior_calorific_value=from_si(states.superior_calorific_value, "J/mol"),
        inferior_calorific_value=from_si(states.inferior_calorific_value, "J/mol"),
        warnings=tuple(warnings),
    )


@dataclass(frozen=True, slots=True)
class Iso6976States:
    """The standard's sums, in SI, in the order the class fills them."""

    molar_mass: float
    compression_factor: float
    relative_density: float
    density_ideal: float
    density_real: float
    superior_calorific_value: float
    inferior_calorific_value: float


def _route(
    components: list[str], z: list[float], volumetric_k: float, energy_k: float
) -> Iso6976States:
    """The standard's arithmetic, in SI, in the order it computes it.

    **One arithmetic, two consumers**: :func:`iso6976` builds the result from this, and it is
    the twin of ``crates/azoth-standards/src/iso6976.rs``. The refusals live here for the
    reason Rust puts them in the calculation.
    """
    if len(components) != len(z):
        raise InvalidInputError(
            "z", f"{len(components)} components and {len(z)} mole fractions, which is not a composition"
        )
    if not components:
        raise InvalidInputError(
            "components", "a gas of no components has no molar mass and no calorific value"
        )
    if any(value < 0.0 for value in z):
        raise InvalidInputError("z", "a mole fraction cannot be negative")

    vol_c = volumetric_k - KELVIN_AT_ZERO_CELSIUS
    energy_c = energy_k - KELVIN_AT_ZERO_CELSIUS
    rows = [_row(name) for name in components]

    molar_mass = 0.0
    summation = 0.0
    relative_ideal = 0.0
    superior = 0.0
    inferior = 0.0
    for row, fraction in zip(rows, z, strict=True):
        molar_mass += fraction * float(row["molar_mass_kg_per_mol"])
        summation += fraction * _column(row, "summation_factor", vol_c, "volumetric")
        relative_ideal += fraction * float(row["molar_mass_kg_per_mol"]) / MOLAR_MASS_AIR
        superior += fraction * _column(row, "superior_calorific_value", energy_c, "energy")
        inferior += fraction * _column(row, "inferior_calorific_value", energy_c, "energy")

    compression_factor = 1.0 - summation * summation
    relative_density = relative_ideal * _air_factor(vol_c) / compression_factor
    reference_temperature = vol_c + KELVIN_AT_ZERO_CELSIUS
    density_ideal = REFERENCE_PRESSURE * molar_mass / (GAS_CONSTANT * reference_temperature)
    density_real = density_ideal / compression_factor

    return Iso6976States(
        molar_mass=molar_mass,
        compression_factor=compression_factor,
        relative_density=relative_density,
        density_ideal=density_ideal,
        density_real=density_real,
        superior_calorific_value=superior,
        inferior_calorific_value=inferior,
    )


def _air_factor(celsius: float) -> float:
    """Air's compression factor at a volumetric reference temperature."""
    for value, factor in AIR_COMPRESSION_FACTOR:
        if abs(celsius - value) < TOLERANCE:
            return factor
    raise InvalidInputError(
        "volumetric_reference_temperature",
        f"the standard's air compression factor is tabulated at 0, 15 and 20 C, and {celsius} C "
        f"is not one of them (15.55 C is the 60 F reference and reads the 15 C value)",
    )


def _column(row: dict[str, str], prefix: str, celsius: float, kind: str) -> float:
    """One tabulated value at a reference temperature.

    The columns are named for the reference temperature they belong to, so the selection is a
    name rather than an index: `summation_factor_15c` is the 15 C column, and the 60 F
    reference reads the one its own column has - which for the volumetric quantities is the
    15 C one, the table carrying no 60 F summation factor.
    """
    names: tuple[tuple[float, str], ...]
    if kind == "volumetric":
        names = ((0.0, "0c"), (15.0, "15c"), (15.55, "15c"), (20.0, "20c"))
    else:
        names = ((0.0, "0c"), (15.0, "15c"), (15.55, "60f"), (20.0, "20c"), (25.0, "25c"))
    # The compiled header spells the calorific values with their unit in the name and the
    # summation factor without one, because that is what the two quantities are.
    unit = "" if prefix == "summation_factor" else "_j_per_mol"
    for value, suffix in names:
        if abs(celsius - value) < TOLERANCE:
            return float(row[f"{prefix}{unit}_{suffix}"])
    raise InvalidInputError(
        "energy_reference_temperature" if kind == "energy" else "volumetric_reference_temperature",
        f"the standard's table carries a {prefix.replace('_', ' ')} at "
        f"{', '.join(f'{value:g}' for value, _ in names)} C, and {celsius} C is not one of them",
    )


def _row(name: str) -> dict[str, str]:
    """One component's row from the standard's compiled table."""
    wanted = name.strip().lower()
    for row in _table():
        if row["name"] == wanted:
            return row
    raise InvalidInputError(
        "components",
        f"{name} has no row in ISO 6976's table. The standard defines 56 substances and this "
        f"library compiles the 40 the component databank also carries; a gas naming one of "
        f"the others has no calorific value here.",
    )


def _table() -> list[dict[str, str]]:
    """The compiled table, parsed once per call as the Rust side does."""
    text = find(ISO6976_CSV).read_text(encoding="utf-8")
    return list(csv.DictReader(text.splitlines()))


def _spec() -> dict[str, object]:
    """The generated spec, imported here so the module has no import-time cycle."""
    from azoth._models_gen import model

    return model("standards.iso6976")
