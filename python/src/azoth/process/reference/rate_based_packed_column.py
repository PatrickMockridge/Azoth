"""``process.rate_based_packed_column`` - the segment model, as a reference implementation.

The Python twin of ``crates/azoth-process/src/segment/``. Every step mirrors the Rust line for
line, and the two are compared case by case by ``python/tests/test_cross_impl.py``.

**The packing is inside the equations**, unlike ``process.packed_column``: the wetted area and
the two film coefficients come from ``azoth.hydraulics.packing_hydraulics`` every segment, and
the fluxes are built from them.

# Three mixtures per segment

The gas system, the liquid system, and an *interface* mix - built by cloning the gas and adding
the liquid's positive-mole components, flashed at the interface temperature and the mean
pressure. A fourth, at the segment's own state, is built for the surface tension.

# Two component lists, and they are not the same one

The interface mixture is stated over the **union** of the two systems' components, and its maps
are *reported* over the active **transfer list**. Collapsing the two is a real defect: a
transfer list of ``["CO2"]`` read against the union's indices picks methane's moles for CO2 and
moves the transfer by three orders of magnitude.

# The reference diffusivity is the class's constant

``averageDiffusivity`` averages ``getEffectiveDiffusionCoefficient(i)``, and NeqSim's flash path
never populates that vector - measured, every entry of both phases is ``0.0``, so the class's
``DEFAULT_*`` constant stands in. The pair matrix on the same state *is* real, so the film
model's ratios are computed against a reference that is not. See ``_FALLBACKS``.
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from typing import Any

from azoth.core.errors import InvalidInputError
from azoth.core.result import RateBasedPackedColumnResult
from azoth.core.units import Q, from_si, quantity
from azoth.core.warnings import Warning, WarningCode
from azoth.eos import components as _components
from azoth.eos.reference._mixture_state import reduced_parameters
from azoth.eos.reference.hydrate_inhibitor_wt import AQUEOUS, GAS, OIL, phase_label
from azoth.eos.reference.molar_enthalpy_entropy import molar_enthalpy_entropy
from azoth.eos.reference.parachor_mixture_surface_tension import (
    parachor_mixture_surface_tension,
)
from azoth.eos.reference.ph_flash import enthalpy_at
from azoth.eos.reference.phase_transport import phase_transport
from azoth.eos.reference.pr_mass_density import pr_mass_density
from azoth.eos.reference.pr_molar_volume import pr_molar_volume
from azoth.eos.reference.pt_flash import pt_flash
from azoth.hydraulics.reference.packing_hydraulics import packing_hydraulics

#: The class's own ``DEFAULT_*`` constants: the value each property falls back to.
DEFAULT_GAS_DIFFUSIVITY = 1.5e-5
DEFAULT_LIQUID_DIFFUSIVITY = 1.5e-9
DEFAULT_SURFACE_TENSION = 0.025
DEFAULT_GAS_THERMAL_CONDUCTIVITY = 0.030
DEFAULT_LIQUID_THERMAL_CONDUCTIVITY = 0.60
DEFAULT_GAS_HEAT_CAPACITY = 2200.0
DEFAULT_LIQUID_HEAT_CAPACITY = 4200.0

#: The floors the snapshot lifts the two diffusivities to.
MIN_GAS_DIFFUSIVITY = 1.0e-7
MIN_LIQUID_DIFFUSIVITY = 1.0e-12

#: The vector ``averageDiffusivity`` averages, as NeqSim's own flash path leaves it.
#:
#: **It is empty, and that is the measurement.** The class averages
#: ``getEffectiveDiffusionCoefficient(i)`` over the components whose value is positive, and
#: NeqSim never populates that vector: on the CO2/water absorber's own state every entry of both
#: phases answers ``0.0``, so the sum is empty and the class's ``DEFAULT_*`` constant stands in.
#: What the same state *does* carry is the pair matrix - ``getDiffusionCoefficient(0, 1)`` is
#: ``9.457085848059925e-7`` m2/s for the gas - so the film model's ratios are real against a
#: reference that is not. `azoth.eos.phase_transport` *does* assemble a vector, and it is not the
#: one this class reads.
_NEQSIM_EFFECTIVE_DIFFUSIVITY: tuple[float, ...] = ()

#: ``maxTransferFractionPerSegment``.
MAX_TRANSFER_FRACTION = 0.35

#: ``maxHeatTransferFractionPerSegment``.
MAX_HEAT_FRACTION = 0.50

#: The class's own defaults for the eleven palette parameters.
DEFAULTS: dict[str, Any] = {
    "column_diameter": 1.0,
    "packed_height": 5.0,
    "number_of_segments": 10,
    "packing_type": "Pall-Ring-50",
    "max_iterations": 30,
    "convergence_tolerance": 1.0e-8,
    "mass_transfer_correction": 1.0,
    "heat_transfer_correction": 1.0,
    "mass_transfer_correlation": "onda_1968",
    "film_model": "maxwell_stefan_matrix",
    "heat_transfer_model": "chilton_colburn_analogy",
    "segment_solver": "sequential_explicit",
    "column_solver": "fixed_point_profile",
}

#: The four values the port does not carry, with the class behind each.
REFUSED = {
    "mass_transfer_correlation": {
        "billet_schultes_1999": (
            "`RateBasedPackedColumn.MassTransferCorrelation.BILLET_SCHULTES_1999` is the "
            "class that "
            "would close it, and it is not a correlation: it is a constant multiplier, "
            "`max(0.1, Ch/0.4)` on `kGa` and `max(0.1, Cp)` on `kLa`."
        )
    },
    "segment_solver": {
        "simultaneous_residual": (
            "`RateBasedPackedColumn.SegmentSolver.SIMULTANEOUS_RESIDUAL` is the class that would "
            "close it, and NeqSim disables its own test for the branch: "
            '`@Disabled("TODO: not working per 19.06.2060")`.'
        )
    },
    "column_solver": {
        "equation_oriented": (
            "`RateBasedPackedColumn.ColumnSolver.EQUATION_ORIENTED` is the class that would close "
            "it - a column-wide damped Newton with homotopy continuation."
        )
    },
}

#: The class's own ``validateSetup`` conditions, as (parameter, predicate, message) so the
#: refusal reads once.
_VALIDATED = (
    "column_diameter",
    "packed_height",
    "number_of_segments",
    "max_iterations",
    "convergence_tolerance",
)


def _finite_positive(value: float, fallback: float) -> float:
    """``finitePositive``."""
    return value if math.isfinite(value) and value > 0.0 else fallback


def _non_negative(value: float, fallback: float) -> float:
    """``finiteNonNegative``."""
    return value if math.isfinite(value) and value >= 0.0 else fallback


def _clamp(value: float, low: float, high: float) -> float:
    """``clamp``."""
    return max(low, min(high, value))


def _phase_view(
    components: list[str], pick: str, t: float, p: float, z: list[float], n: float
) -> dict[str, Any]:
    """One phase of a flashed system: the class's ``PhaseInterface``, as a value.

    ``pick`` is ``"gas"`` or ``"liquid"``, and the phase is chosen by ``PhaseEos.init``'s label
    rule rather than by ``beta`` - the same rule ``three_phase_separator`` reads.
    """
    mixture, ideal_gas = _components.mixture_of(components, eos="pr")
    flash = pt_flash(mixture, quantity(t, "K"), quantity(p, "Pa"), z)
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

    if pick == "gas":
        matches = [entry for entry in labelled if entry[0] == GAS]
    else:
        # `getLiquidPhase`'s own ladder: AQUEOUS, then OIL, then the first phase that is not
        # the gas, then phase zero.
        matches = (
            [entry for entry in labelled if entry[0] == AQUEOUS]
            or [entry for entry in labelled if entry[0] == OIL]
            or [entry for entry in labelled if entry[0] != GAS]
        )
    chosen = matches[0] if matches else labelled[0]

    label, phase_z, root = chosen
    if len(labelled) == 1:
        share = 1.0
    elif label == GAS:
        share = flash.beta if flash.beta is not None else 1.0
    else:
        share = 1.0 - (flash.beta if flash.beta is not None else 0.0)

    molar_mass = sum(
        (component.molar_mass.to("kg/mol").magnitude if component.molar_mass else 0.0) * fraction
        for component, fraction in zip(mixture.components, phase_z, strict=True)
    )
    transport = phase_transport(components, label, quantity(t, "K"), quantity(p, "Pa"), phase_z)

    # `getDensity("kg/m3")`: the cubic's volume at the phase's own root, with the Peneloux
    # shift. **The root is the flash's own** - a phase of a split sits on the root the split put
    # it on, and asking the cubic again from its composition can answer a different one.
    volume = pr_molar_volume(root, quantity(t, "K"), quantity(p, "Pa")).v.to("m**3/mol").magnitude
    shift = mixture.volume_shift(phase_z)
    density = (
        pr_mass_density(quantity(molar_mass, "kg/mol"), quantity(volume - shift, "m**3/mol"))
        .rho.to("kg/m**3")
        .magnitude
    )
    cp_mass = (
        molar_enthalpy_entropy(
            mixture, ideal_gas, quantity(t, "K"), quantity(p, "Pa"), phase_z, root
        )
        .cp.to("J/(mol*K)")
        .magnitude
        / molar_mass
    )
    return {
        "z": phase_z,
        "n": n * share,
        "t": t,
        "p": p,
        "molar_mass": molar_mass,
        "density": density,
        "cp_mass": cp_mass,
        "kind": label,
        "transport": transport,
        "components": components,
    }


def _mole_fraction(view: dict[str, Any], name: str) -> float:
    """A component's mole fraction in a phase, or zero where it is absent."""
    components: list[str] = view["components"]
    if name not in components:
        return 0.0
    return max(0.0, float(view["z"][components.index(name)]))


def _average_diffusivity(effective: Sequence[Any], gas_phase: bool) -> tuple[float, bool]:
    """``averageDiffusivity``: the mean of the positive-finite entries, or the class's constant.

    **The vector it averages is empty.** NeqSim's flash path never populates
    ``getEffectiveDiffusionCoefficient``, so the sum is empty on every state and the default
    stands in; the Rust twin carries the same rule with the same measurement.
    """
    values = [
        value.to("m**2/s").magnitude
        for value in effective
        if math.isfinite(value.to("m**2/s").magnitude) and value.to("m**2/s").magnitude > 0.0
    ]
    if values:
        return sum(values) / len(values), False
    return (DEFAULT_GAS_DIFFUSIVITY if gas_phase else DEFAULT_LIQUID_DIFFUSIVITY), True


def _estimate_surface_tension(gas: dict[str, Any], liquid: dict[str, Any]) -> tuple[float, bool]:
    """``estimateSurfaceTension``: a fourth mixture, at the segment's own state.

    **The class asks for a two-phase mixture with a gas in it, and takes the pair's surface
    tension only where it is positive.** Both gates are reached on its own states:

    - the mixture it builds - the gas cloned, the liquid's positive-mole components added,
      flashed at the gas's temperature and pressure - is gas-plus-aqueous, and
      ``InterfaceProperties.getSurfaceTension`` answers ``0.0`` for that pair. Measured across
      dissolved-CO2 loadings of ``0``, ``1e-6``, ``1e-3``, ``1e-2`` and ``0.04``.
    - so ``isFinitePositive(0.0)`` is false and the constant stands in, which is what the
      capture's wetted area is built from: the parachor form answers ``0.05303`` on the same
      mixture, and that moves the wetted area from ``56.9`` to ``32.9``.

    The parachor form is a different branch - an oil pair takes it - and it is not wrong here;
    it is simply not what this class's own path reads.
    """
    names = list(gas["components"])
    for name in liquid["components"]:
        if name not in names:
            names.append(name)
    total = gas["n"] + liquid["n"]
    z = [
        (gas["n"] * _mole_fraction(gas, name) + liquid["n"] * _mole_fraction(liquid, name)) / total
        for name in names
    ]
    try:
        mixed = _phase_view(names, "gas", gas["t"], gas["p"], z, total)
        other = _phase_view(names, "liquid", gas["t"], gas["p"], z, total)
    except (InvalidInputError, ValueError, ZeroDivisionError):
        return DEFAULT_SURFACE_TENSION, True
    if mixed["kind"] != GAS or mixed["n"] <= 0.0 or other["n"] <= 0.0:
        return DEFAULT_SURFACE_TENSION, True
    if other["kind"] == AQUEOUS:
        return DEFAULT_SURFACE_TENSION, True
    try:
        mixture, _ = _components.mixture_of(names, eos="pr")
        parachors = [component.parachor for component in mixture.components]
        if any(value <= 0.0 for value in parachors):
            return DEFAULT_SURFACE_TENSION, True
        sigma = parachor_mixture_surface_tension(
            parachors,
            quantity(mixed["density"], "kg/m**3"),
            quantity(mixed["molar_mass"], "kg/mol"),
            mixed["z"],
            quantity(other["density"], "kg/m**3"),
            quantity(other["molar_mass"], "kg/mol"),
            other["z"],
        )
        value = sigma.sigma.to("N/m").magnitude
        if math.isfinite(value) and value > 0.0:
            return value, False
    except (InvalidInputError, ValueError, ZeroDivisionError):
        pass
    return DEFAULT_SURFACE_TENSION, True


def _film_coefficient(
    view: dict[str, Any],
    index: int,
    base: float,
    reference: float,
    gas_phase: bool,
    matrix_model: bool,
) -> float:
    """``calculateFilmCoefficient``."""
    if not matrix_model:
        return base
    return _maxwell_stefan_film_coefficient(view, index, base, reference, gas_phase)


def _binary_diffusivity(view: dict[str, Any], first: int, second: int, gas_phase: bool) -> float:
    """``binaryDiffusivity``: the pair coefficient, with the class's two-step fallback."""
    if first == second:
        return 0.0
    matrix = float(view["transport"].d_binary[first][second].to("m**2/s").magnitude)
    if math.isfinite(matrix) and matrix > 0.0:
        return matrix
    effective = float(view["transport"].d_effective[first].to("m**2/s").magnitude)
    if math.isfinite(effective) and effective > 0.0:
        return effective
    return DEFAULT_GAS_DIFFUSIVITY if gas_phase else DEFAULT_LIQUID_DIFFUSIVITY


def _scale_film_coefficient(base: float, diffusivity: float, reference: float) -> float:
    """``scaleFilmCoefficient``."""
    reference = _finite_positive(reference, diffusivity)
    scaled = base * _finite_positive(diffusivity, reference) / reference
    return _clamp(scaled, base * 0.02, base * 50.0)


def _mixture_diffusivity(
    view: dict[str, Any], index: int, reference: float, gas_phase: bool
) -> float:
    """``mixtureDiffusivityForComponent``."""
    count = len(view["z"])
    if count <= 1:
        return reference
    resistance = 0.0
    for other in range(count):
        if other != index:
            resistance += max(0.0, view["z"][other]) / _binary_diffusivity(
                view, index, other, gas_phase
            )
    if math.isfinite(resistance) and resistance > 0.0:
        return 1.0 / resistance
    return reference


def _binary_film_coefficient(
    view: dict[str, Any], first: int, second: int, base: float, reference: float, gas_phase: bool
) -> float:
    """``binaryFilmCoefficient``."""
    if first == second:
        return base
    return _scale_film_coefficient(
        base, _binary_diffusivity(view, first, second, gas_phase), reference
    )


def _maxwell_stefan_film_coefficient(
    view: dict[str, Any], index: int, base: float, reference: float, gas_phase: bool
) -> float:
    """``maxwellStefanFilmCoefficient``.

    The resistance matrix is ``(n - 1) x (n - 1)`` and only its diagonal is read. **A phase of
    two components never reaches the matrix**: its one-by-one inverse is the same scalar, and the
    class's own absorber states are two-component.
    """
    if not (math.isfinite(base) and base > 0.0):
        return 0.0
    count = len(view["z"])
    reduced = max(0, count - 1)
    if index >= reduced:
        return _scale_film_coefficient(
            base, _mixture_diffusivity(view, index, reference, gas_phase), reference
        )

    matrix = [[0.0] * reduced for _ in range(reduced)]
    for row in range(reduced):
        row_sum = 0.0
        reference_coefficient = _binary_film_coefficient(
            view, row, reduced, base, reference, gas_phase
        )
        for column in range(count):
            binary = _binary_film_coefficient(view, row, column, base, reference, gas_phase)
            if row != column:
                row_sum += max(0.0, view["z"][column]) / binary
            if column < reduced:
                matrix[row][column] = -max(0.0, view["z"][row]) * (
                    1.0 / binary - 1.0 / reference_coefficient
                )
        matrix[row][row] += row_sum + max(0.0, view["z"][row]) / reference_coefficient

    try:
        # `n` unit right-hand sides - the route JAMA's own `inverse` takes - rather than a
        # second elimination. The twin solves it directly because there is one implementation
        # to compare against, and the comparison is the point.
        inverse = _inverse(matrix)
        coefficient = inverse[index][index]
        if math.isfinite(coefficient) and coefficient > 0.0:
            return _clamp(coefficient, base * 0.02, base * 50.0)
        return base
    except ZeroDivisionError:
        return _scale_film_coefficient(
            base, _mixture_diffusivity(view, index, reference, gas_phase), reference
        )


def _inverse(matrix: list[list[float]]) -> list[list[float]]:
    """Gauss-Jordan inversion, which raises ``ZeroDivisionError`` on a singular matrix - the
    class's ``SingularMatrixException``, and a different fallback from a non-finite answer."""
    size = len(matrix)
    work = [row[:] + [1.0 if i == j else 0.0 for j in range(size)] for i, row in enumerate(matrix)]
    for column in range(size):
        pivot = max(range(column, size), key=lambda row: abs(work[row][column]))
        if work[pivot][column] == 0.0:
            raise ZeroDivisionError("singular")
        work[column], work[pivot] = work[pivot], work[column]
        divisor = work[column][column]
        work[column] = [value / divisor for value in work[column]]
        for row in range(size):
            if row != column and work[row][column] != 0.0:
                factor = work[row][column]
                work[row] = [
                    value - factor * other
                    for value, other in zip(work[row], work[column], strict=True)
                ]
    return [row[size:] for row in work]


def _combine_film_fluxes(gas_flux: float, liquid_flux: float) -> float:
    """``combineFilmFluxes``: the harmonic mean with the class's three guards."""
    if not math.isfinite(gas_flux) or not math.isfinite(liquid_flux):
        return 0.0
    if abs(gas_flux) < 1.0e-30 or abs(liquid_flux) < 1.0e-30:
        return 0.0
    if math.copysign(1.0, gas_flux) != math.copysign(1.0, liquid_flux):
        return 0.0
    magnitude = 1.0 / (1.0 / abs(gas_flux) + 1.0 / abs(liquid_flux))
    return math.copysign(1.0, gas_flux) * magnitude


def _volumetric_heat_transfer_coefficient(
    mass_transfer: float,
    density: float,
    heat_capacity: float,
    viscosity: float,
    diffusivity: float,
    conductivity: float,
    model_none: bool,
    correction: float,
) -> float:
    """``calculateVolumetricHeatTransferCoefficient``: the Chilton-Colburn analogy."""
    values = (mass_transfer, density, heat_capacity, viscosity, diffusivity, conductivity)
    if model_none or any(not (math.isfinite(value) and value > 0.0) for value in values):
        return 0.0
    prandtl = heat_capacity * viscosity / conductivity
    schmidt = viscosity / (density * diffusivity)
    if not (prandtl > 0.0 and math.isfinite(prandtl)) or not (
        schmidt > 0.0 and math.isfinite(schmidt)
    ):
        return 0.0
    return float(
        mass_transfer * density * heat_capacity * (schmidt / prandtl) ** (2.0 / 3.0) * correction
    )


def _combine_heat_transfer_coefficients(gas: float, liquid: float) -> float:
    """``combineHeatTransferCoefficients``."""
    if not (gas > 0.0 and math.isfinite(gas) and liquid > 0.0 and math.isfinite(liquid)):
        return 0.0
    return 1.0 / (1.0 / gas + 1.0 / liquid)


def _interface_temperature(gas_t: float, liquid_t: float, gas: float, liquid: float) -> float:
    """``calculateInterfaceTemperature``."""
    if not (gas > 0.0 and math.isfinite(gas) and liquid > 0.0 and math.isfinite(liquid)):
        return 0.5 * (gas_t + liquid_t)
    return (gas * gas_t + liquid * liquid_t) / (gas + liquid)


def _interface_equilibrium(
    gas: dict[str, Any],
    liquid: dict[str, Any],
    components: list[str],
    temperature: float,
    pressure: float,
) -> dict[str, Any]:
    """``calculateInterfaceEquilibrium``: the third mixture, flashed at the interface state."""
    names = list(gas["components"])
    for name in liquid["components"]:
        if _moles_of(liquid, name) > 0.0 and name not in names:
            names.append(name)
    amounts = [_moles_of(gas, name) + _moles_of(liquid, name) for name in names]
    total = sum(amounts)
    try:
        if total <= 0.0:
            raise InvalidInputError("components", "no moles to mix")
        z = [amount / total for amount in amounts]
        gas_phase = _phase_view(names, "gas", temperature, pressure, z, total)
        liquid_phase = _phase_view(names, "liquid", temperature, pressure, z, total)
    except (InvalidInputError, ValueError, ZeroDivisionError):
        # The class's own branch: the two *unmixed* phases, which keeps a segment alive when the
        # combined mixture will not flash.
        names = list(gas["components"])
        gas_phase = _phase_view(gas["components"], "gas", gas["t"], gas["p"], gas["z"], gas["n"])
        liquid_phase = _phase_view(
            liquid["components"], "liquid", liquid["t"], liquid["p"], liquid["z"], liquid["n"]
        )

    gas_fractions: list[float] = []
    liquid_fractions: list[float] = []
    ratios: list[float] = []
    for component in components:
        x = _mole_fraction(liquid_phase, component)
        y = _mole_fraction(gas_phase, component)
        gas_fractions.append(y)
        liquid_fractions.append(x)
        ratios.append(max(1.0e-12, y / x) if (x > 1.0e-12 and y >= 0.0) else 1.0)
    return {
        "temperature": temperature,
        "components": components,
        "gas": gas_fractions,
        "liquid": liquid_fractions,
        "ratios": ratios,
    }


def _moles_of(stream: dict[str, Any], name: str) -> float:
    """``componentMoles(system, name)``."""
    if name not in stream["components"]:
        return 0.0
    return max(0.0, float(stream["z"][stream["components"].index(name)]) * float(stream["n"]))


def _transport_snapshot(
    gas: dict[str, Any], liquid: dict[str, Any], segment_height: float, settings: dict[str, Any]
) -> tuple[dict[str, Any], dict[str, Any], dict[str, Any], dict[str, bool]]:
    """``calculateTransportSnapshot``."""
    gas_phase = _phase_view(gas["components"], "gas", gas["t"], gas["p"], gas["z"], gas["n"])
    liquid_phase = _phase_view(
        liquid["components"], "liquid", liquid["t"], liquid["p"], liquid["z"], liquid["n"]
    )

    gas_density = _finite_positive(gas_phase["density"], 1.0)
    liquid_density = _finite_positive(liquid_phase["density"], 800.0)
    gas_viscosity = _finite_positive(gas_phase["transport"].mu.to("Pa*s").magnitude, 1.0e-5)
    liquid_viscosity = _finite_positive(liquid_phase["transport"].mu.to("Pa*s").magnitude, 1.0e-3)

    fallbacks = {
        "diffusivity": False,
        "surface_tension": False,
        "thermal_conductivity": False,
        "heat_capacity": False,
    }
    gas_diffusivity, gas_took = _average_diffusivity(_NEQSIM_EFFECTIVE_DIFFUSIVITY, True)
    liquid_diffusivity, liquid_took = _average_diffusivity(_NEQSIM_EFFECTIVE_DIFFUSIVITY, False)
    fallbacks["diffusivity"] = gas_took or liquid_took

    surface_tension, tension_took = _estimate_surface_tension(gas, liquid)
    fallbacks["surface_tension"] = tension_took

    gas_heat_capacity = _finite_positive(gas_phase["cp_mass"], DEFAULT_GAS_HEAT_CAPACITY)
    liquid_heat_capacity = _finite_positive(liquid_phase["cp_mass"], DEFAULT_LIQUID_HEAT_CAPACITY)
    fallbacks["heat_capacity"] = (
        gas_heat_capacity == DEFAULT_GAS_HEAT_CAPACITY
        or liquid_heat_capacity == DEFAULT_LIQUID_HEAT_CAPACITY
    )
    gas_conductivity = _finite_positive(
        gas_phase["transport"].k.to("W/(m*K)").magnitude, DEFAULT_GAS_THERMAL_CONDUCTIVITY
    )
    liquid_conductivity = _finite_positive(
        liquid_phase["transport"].k.to("W/(m*K)").magnitude, DEFAULT_LIQUID_THERMAL_CONDUCTIVITY
    )
    fallbacks["thermal_conductivity"] = (
        gas_conductivity == DEFAULT_GAS_THERMAL_CONDUCTIVITY
        or liquid_conductivity == DEFAULT_LIQUID_THERMAL_CONDUCTIVITY
    )

    hydraulics = packing_hydraulics(
        settings["packing_type"],
        quantity(settings["column_diameter"], "m"),
        quantity(max(segment_height, 1.0e-9), "m"),
        quantity(gas_phase["n"] * gas_phase["molar_mass"], "kg/s"),
        quantity(liquid_phase["n"] * liquid_phase["molar_mass"], "kg/s"),
        quantity(gas_density, "kg/m**3"),
        quantity(liquid_density, "kg/m**3"),
        quantity(gas_viscosity, "Pa*s"),
        quantity(liquid_viscosity, "Pa*s"),
        quantity(surface_tension, "N/m"),
        quantity(gas_diffusivity, "m**2/s"),
        quantity(liquid_diffusivity, "m**2/s"),
        1.0,
    )
    k_ga = _non_negative(hydraulics.k_ga, 0.0) * settings["mass_transfer_correction"]
    k_la = _non_negative(hydraulics.k_la, 0.0) * settings["mass_transfer_correction"]

    none = settings["heat_transfer_model"] == "none"
    gas_heat = _volumetric_heat_transfer_coefficient(
        k_ga,
        gas_density,
        gas_heat_capacity,
        gas_viscosity,
        gas_diffusivity,
        gas_conductivity,
        none,
        settings["heat_transfer_correction"],
    )
    liquid_heat = _volumetric_heat_transfer_coefficient(
        k_la,
        liquid_density,
        liquid_heat_capacity,
        liquid_viscosity,
        liquid_diffusivity,
        liquid_conductivity,
        none,
        settings["heat_transfer_correction"],
    )
    snapshot = {
        "gas_density": gas_density,
        "liquid_density": liquid_density,
        "gas_viscosity": gas_viscosity,
        "liquid_viscosity": liquid_viscosity,
        "gas_diffusivity": max(gas_diffusivity, MIN_GAS_DIFFUSIVITY),
        "liquid_diffusivity": max(liquid_diffusivity, MIN_LIQUID_DIFFUSIVITY),
        "wetted_area": _non_negative(hydraulics.wetted_area, 0.0),
        "k_ga": k_ga,
        "k_la": k_la,
        "gas_heat_capacity": gas_heat_capacity,
        "liquid_heat_capacity": liquid_heat_capacity,
        "gas_heat_transfer_coefficient": gas_heat,
        "liquid_heat_transfer_coefficient": liquid_heat,
        "overall_heat_transfer_coefficient": _combine_heat_transfer_coefficients(
            gas_heat, liquid_heat
        ),
        "interface_temperature": _interface_temperature(
            gas["t"], liquid["t"], gas_heat, liquid_heat
        ),
        "pressure_drop_per_meter": _non_negative(
            hydraulics.pressure_drop_per_meter.to("Pa").magnitude, 0.0
        ),
        "percent_flood": _non_negative(hydraulics.percent_flood, 0.0),
    }
    return snapshot, gas_phase, liquid_phase, fallbacks


def _total_enthalpy(stream: dict[str, Any]) -> float:
    """A system's total enthalpy, in W."""
    return float(stream["h"]) * float(stream["n"])


def _at_temperature(stream: dict[str, Any], temperature: float) -> dict[str, Any]:
    """The same system at a different temperature, with its enthalpy re-derived."""
    mixture, ideal_gas = _components.mixture_of(stream["components"], eos="pr")
    enthalpy, _ = enthalpy_at(mixture, ideal_gas, temperature, stream["p"], list(stream["z"]))
    return {**stream, "t": temperature, "h": enthalpy}


def _apply_component_transfer(
    stream: dict[str, Any], component: str, delta: float
) -> dict[str, Any]:
    """``addComponent``: move moles into or out of a system, at its own T and P."""
    components = list(stream["components"])
    if component not in components:
        components.append(component)
    moles = [_moles_of(stream, name) for name in components]
    index = components.index(component)
    moles[index] += delta
    if moles[index] < 0.0:
        raise InvalidInputError(
            component,
            f"a transfer of {delta} mol/s would leave a negative inventory of {component}",
        )
    total = stream["n"] + delta
    z = [value / total for value in moles] if total > 0.0 else [0.0] * len(components)
    return _at_temperature({**stream, "components": components, "z": z, "n": total}, stream["t"])


def _segment(
    index: int,
    gas_in: dict[str, Any],
    liquid_in: dict[str, Any],
    segment_height: float,
    segment_volume: float,
    transfer_components: list[str],
    settings: dict[str, Any],
) -> dict[str, Any]:
    """``calculateSegment``: one segment's outlets and its record."""
    gas = dict(gas_in)
    liquid = dict(liquid_in)
    inlet_enthalpy = _total_enthalpy(gas) + _total_enthalpy(liquid)

    snapshot, gas_phase, liquid_phase, fallbacks = _transport_snapshot(
        gas, liquid, segment_height, settings
    )
    components = list(transfer_components) or list(
        dict.fromkeys(list(gas["components"]) + list(liquid["components"]))
    )
    interface = _interface_equilibrium(
        gas, liquid, components, snapshot["interface_temperature"], 0.5 * (gas["p"] + liquid["p"])
    )

    transfers: list[tuple[str, float]] = []
    heat_rate = 0.0
    if segment_height > 0.0:
        for component in components:
            transfer = _component_transfer(
                component,
                gas,
                liquid,
                gas_phase,
                liquid_phase,
                snapshot,
                interface,
                segment_volume,
                settings,
            )
            if transfer != 0.0:
                gas = _apply_component_transfer(gas, component, -transfer)
                liquid = _apply_component_transfer(liquid, component, transfer)
                transfers.append((component, transfer))
        gas_rate = (
            _finite_positive(gas_phase["n"] * gas_phase["molar_mass"], 0.0)
            * snapshot["gas_heat_capacity"]
        )
        liquid_rate = (
            _finite_positive(liquid_phase["n"] * liquid_phase["molar_mass"], 0.0)
            * snapshot["liquid_heat_capacity"]
        )
        heat_rate, gas_t, liquid_t = _apply_interphase_heat_transfer(
            gas["t"],
            liquid["t"],
            gas_rate,
            liquid_rate,
            snapshot["overall_heat_transfer_coefficient"],
            segment_volume,
            settings["heat_transfer_model"] == "none",
        )
        gas = _at_temperature(gas, gas_t)
        liquid = _at_temperature(liquid, liquid_t)

    return {
        "gas": gas,
        "liquid": liquid,
        "result": {
            "number": index + 1,
            "height_from_bottom": (index + 0.5) * segment_height,
            "gas_temperature": gas["t"],
            "liquid_temperature": liquid["t"],
            "gas_pressure": gas["p"],
            "liquid_pressure": liquid["p"],
            "gas_molar_flow": gas["n"],
            "liquid_molar_flow": liquid["n"],
            "net_molar_transfer": sum(value for _, value in transfers),
            "heat_transfer_rate": heat_rate,
            "component_transfer": transfers,
            "interface": interface,
            "enthalpy_balance_residual": _total_enthalpy(gas)
            + _total_enthalpy(liquid)
            - inlet_enthalpy,
            "fallbacks": fallbacks,
            **{key: value for key, value in snapshot.items()},
        },
    }


def _component_transfer(
    component: str,
    gas: dict[str, Any],
    liquid: dict[str, Any],
    gas_phase: dict[str, Any],
    liquid_phase: dict[str, Any],
    snapshot: dict[str, Any],
    interface: dict[str, Any],
    segment_volume: float,
    settings: dict[str, Any],
) -> float:
    """``calculateComponentTransfer``: the unbounded film transfer, then the inventory cap."""
    ratios = dict(zip(interface["components"], interface["ratios"], strict=True))
    gas_fractions = dict(zip(interface["components"], interface["gas"], strict=True))
    liquid_fractions = dict(zip(interface["components"], interface["liquid"], strict=True))
    k_value = ratios.get(component, 1.0)
    gas_fraction = _mole_fraction(gas_phase, component)
    liquid_fraction = _mole_fraction(liquid_phase, component)
    gas_interface = gas_fractions.get(component, _clamp(k_value * liquid_fraction, 0.0, 0.999999))
    liquid_interface = liquid_fractions.get(
        component,
        _clamp(gas_interface / k_value, 0.0, 0.999999) if k_value > 1.0e-12 else liquid_fraction,
    )
    gas_driving = gas_fraction - gas_interface
    if abs(gas_driving) < 1.0e-12 or snapshot["k_ga"] <= 0.0 or snapshot["k_la"] <= 0.0:
        return 0.0

    index = gas["components"].index(component)
    matrix = settings["film_model"] == "maxwell_stefan_matrix"
    gas_film = _film_coefficient(
        gas_phase, index, snapshot["k_ga"], snapshot["gas_diffusivity"], True, matrix
    )
    liquid_film = _film_coefficient(
        liquid_phase, index, snapshot["k_la"], snapshot["liquid_diffusivity"], False, matrix
    )
    gas_flux = (
        gas_film
        * (snapshot["gas_density"] / _finite_positive(gas_phase["molar_mass"], 0.020))
        * gas_driving
    )
    liquid_flux = (
        liquid_film
        * (snapshot["liquid_density"] / _finite_positive(liquid_phase["molar_mass"], 0.020))
        * (liquid_interface - liquid_fraction)
    )
    transfer_density = _combine_film_fluxes(gas_flux, liquid_flux)
    if transfer_density == 0.0:
        y_star = _clamp(k_value * liquid_fraction, 0.0, 0.999999)
        overall = 1.0 / (1.0 / gas_film + max(k_value, 1.0e-12) / liquid_film)
        transfer_density = (
            overall
            * (snapshot["gas_density"] / _finite_positive(gas_phase["molar_mass"], 0.020))
            * (gas_fraction - y_star)
        )
    proposed = transfer_density * segment_volume
    if proposed > 0.0:
        return min(proposed, max(0.0, _moles_of(gas, component) * MAX_TRANSFER_FRACTION))
    if proposed < 0.0:
        return -min(-proposed, max(0.0, _moles_of(liquid, component) * MAX_TRANSFER_FRACTION))
    return 0.0


def _apply_interphase_heat_transfer(
    gas_t: float,
    liquid_t: float,
    gas_rate: float,
    liquid_rate: float,
    overall: float,
    segment_volume: float,
    model_none: bool,
) -> tuple[float, float, float]:
    """``applyInterphaseHeatTransfer``."""
    if model_none or not (overall > 0.0 and math.isfinite(overall)):
        return 0.0, gas_t, liquid_t
    difference = gas_t - liquid_t
    if abs(difference) < 1.0e-12:
        return 0.0, gas_t, liquid_t
    if not (gas_rate > 0.0 and math.isfinite(gas_rate)) or not (
        liquid_rate > 0.0 and math.isfinite(liquid_rate)
    ):
        return 0.0, gas_t, liquid_t
    rate = overall * segment_volume * difference
    maximum = min(gas_rate, liquid_rate) * abs(difference) * MAX_HEAT_FRACTION
    rate = math.copysign(1.0, rate) * min(abs(rate), maximum)
    if rate == 0.0:
        return 0.0, gas_t, liquid_t
    return (
        rate,
        max(1.0, gas_t - rate / gas_rate),
        max(1.0, liquid_t + rate / liquid_rate),
    )


def _states(
    gas_components: list[str],
    gas_n: float,
    gas_z: list[float],
    gas_p: float,
    gas_t: float,
    liquid_components: list[str],
    liquid_n: float,
    liquid_z: list[float],
    liquid_p: float,
    liquid_t: float,
    transfer_components: list[str] | None,
    column_diameter: float,
    packed_height: float,
    number_of_segments: int,
    packing_type: str,
    max_iterations: int,
    convergence_tolerance: float,
    mass_transfer_correction: float,
    mass_transfer_correlation: str,
    film_model: str,
    heat_transfer_model: str,
    segment_solver: str,
    column_solver: str,
) -> dict[str, Any]:
    """The profile, factored out so the layer diff is a *read* of the same route.

    ``layers.py`` imports this rather than recomputing anything: a layer computed by a route the
    kernel does not take would be a second implementation.
    """
    for parameter, value in (
        ("mass_transfer_correlation", mass_transfer_correlation),
        ("segment_solver", segment_solver),
        ("column_solver", column_solver),
    ):
        refused = REFUSED.get(parameter, {}).get(value)
        if refused is not None:
            raise InvalidInputError(parameter, f"`{value}` is not ported: {refused}")
    if column_diameter <= 0.0:
        raise InvalidInputError("column_diameter", "Column diameter must be positive")
    if packed_height < 0.0:
        raise InvalidInputError("packed_height", "Packed height can not be negative")
    if number_of_segments < 1:
        raise InvalidInputError("number_of_segments", "At least one segment is required")
    if max_iterations < 1:
        raise InvalidInputError("max_iterations", "the loop runs at least once")
    if convergence_tolerance <= 0.0:
        raise InvalidInputError("convergence_tolerance", "a convergence gate is positive")

    def feed(components: list[str], z: list[float], n: float, p: float, t: float) -> dict[str, Any]:
        mixture, ideal_gas = _components.mixture_of(components, eos="pr")
        enthalpy, _ = enthalpy_at(mixture, ideal_gas, t, p, list(z))
        return {
            "components": list(components),
            "n": n,
            "z": list(z),
            "p": p,
            "t": t,
            "h": enthalpy,
        }

    gas_in = feed(gas_components, gas_z, gas_n, gas_p, gas_t)
    liquid_in = feed(liquid_components, liquid_z, liquid_n, liquid_p, liquid_t)
    whitelist = list(transfer_components) if transfer_components else []
    settings = {
        "column_diameter": column_diameter,
        "packed_height": packed_height,
        "packing_type": packing_type,
        "mass_transfer_correction": mass_transfer_correction,
        "heat_transfer_correction": 1.0,
        "heat_transfer_model": heat_transfer_model,
        "film_model": film_model,
    }
    segment_height = packed_height / number_of_segments
    segment_volume = math.pi * column_diameter * column_diameter / 4.0 * segment_height

    liquid_entering = [liquid_in] * number_of_segments
    previous_gas: dict[str, Any] | None = None
    previous_liquid: dict[str, Any] | None = None
    solution: dict[str, Any] | None = None
    iterations = 0
    residual = math.inf
    converged = False
    for iteration in range(1, max_iterations + 1):
        gas_current = gas_in
        liquid_leaving: list[dict[str, Any]] = []
        segments: list[dict[str, Any]] = []
        for index, entering in enumerate(liquid_entering):
            computation = _segment(
                index, gas_current, entering, segment_height, segment_volume, whitelist, settings
            )
            gas_current = computation["gas"]
            liquid_leaving.append(computation["liquid"])
            segments.append(computation["result"])
        pass_gas, pass_liquid = gas_current, liquid_leaving[0]
        residual = _outlet_residual(previous_gas, pass_gas, previous_liquid, pass_liquid)
        iterations = iteration
        solution = {"gas": pass_gas, "liquid": pass_liquid, "segments": segments}
        if residual <= convergence_tolerance or packed_height == 0.0:
            converged = True
            break
        previous_gas, previous_liquid = pass_gas, pass_liquid
        liquid_entering = [
            liquid_in if index + 1 == number_of_segments else liquid_leaving[index + 1]
            for index in range(number_of_segments)
        ]
    assert solution is not None

    totals: dict[str, float] = {}
    total_absolute = 0.0
    fallbacks = {
        "diffusivity": False,
        "surface_tension": False,
        "thermal_conductivity": False,
        "heat_capacity": False,
    }
    for segment in solution["segments"]:
        for key, value in segment["fallbacks"].items():
            fallbacks[key] = fallbacks[key] or value
        for component, transfer in segment["component_transfer"]:
            total_absolute += abs(transfer)
            totals[component] = totals.get(component, 0.0) + transfer
    return {
        "gas_out": solution["gas"],
        "liquid_out": solution["liquid"],
        "segments": solution["segments"],
        "iterations": iterations,
        "residual": residual,
        "converged": converged,
        "total_absolute": total_absolute,
        "totals": totals,
        "fallbacks": fallbacks,
    }


def _outlet_residual(
    previous_gas: dict[str, Any] | None,
    current_gas: dict[str, Any],
    previous_liquid: dict[str, Any] | None,
    current_liquid: dict[str, Any],
) -> float:
    """``calculateOutletResidual``: each side compared to its own kind, never across."""
    if previous_gas is None or previous_liquid is None:
        return math.inf
    names = list(
        dict.fromkeys(list(current_gas["components"]) + list(current_liquid["components"]))
    )
    residual = 0.0
    for component in names:
        residual = max(
            residual,
            abs(_moles_of(previous_gas, component) - _moles_of(current_gas, component)),
            abs(_moles_of(previous_liquid, component) - _moles_of(current_liquid, component)),
        )
    return residual


def rate_based_packed_column(
    gas_components: list[str],
    gas_n: Q,
    gas_z: list[float],
    gas_p: Q,
    gas_t: Q,
    liquid_components: list[str],
    liquid_n: Q,
    liquid_z: list[float],
    liquid_p: Q,
    liquid_t: Q,
    transfer_components: list[str] | None = None,
    column_diameter: Q | None = None,
    packed_height: Q | None = None,
    number_of_segments: float | None = None,
    packing_type: str | None = None,
    max_iterations: float | None = None,
    convergence_tolerance: Q | None = None,
    mass_transfer_correction: float | None = None,
    mass_transfer_correlation: str | None = None,
    film_model: str | None = None,
    heat_transfer_model: str | None = None,
    segment_solver: str | None = None,
    column_solver: str | None = None,
) -> RateBasedPackedColumnResult:
    """Solve a rate-based packed column.

    Args:
        gas_components: the gas's substances, which enters the bottom segment.
        gas_n: the gas's molar flow.
        gas_z: the gas's composition.
        gas_p: the gas's pressure.
        gas_t: the gas's temperature.
        liquid_components: the solvent's substances, which enters the top segment.
        liquid_n: the solvent's molar flow.
        liquid_z: the solvent's composition.
        liquid_p: the solvent's pressure.
        liquid_t: the solvent's temperature.
        transfer_components: the components the two films walk; omitted is the union.
        column_diameter: the column's internal diameter; the class defaults it to 1.0 m.
        packed_height: the packed height; the class defaults it to 5.0 m. Zero is valid.
        number_of_segments: the axial slices; the class defaults it to ten.
        packing_type: the packing's name; a name nothing carries is ``Pall-Ring-50``.
        max_iterations: the iteration cap; the class defaults it to thirty. A solve that
            reaches it publishes its last iterate.
        convergence_tolerance: the outlet gate, mol/s; the class defaults it to 1e-8.
        mass_transfer_correction: a scale on both film coefficients.
        mass_transfer_correlation: ``onda_1968``; ``billet_schultes_1999`` is refused.
        film_model: ``maxwell_stefan_matrix`` or ``overall_two_resistance``.
        heat_transfer_model: ``chilton_colburn_analogy`` or ``none``, which is exactly zero.
        segment_solver: ``sequential_explicit``; ``simultaneous_residual`` is refused.
        column_solver: ``fixed_point_profile``; ``equation_oriented`` is refused.

    Returns:
        The four ports' records and the per-segment profile.
    """
    # **`x or default` is wrong for a quantity**, and it is wrong in the direction that hides: a
    # zero-valued `Quantity` is falsy, so `packed_height or default` would turn the class's own
    # zero-height state into a five-metre bed. `is None` is the test everywhere below.
    diameter = column_diameter if column_diameter is not None else _default("column_diameter", "m")
    height = packed_height if packed_height is not None else _default("packed_height", "m")
    gate = (
        convergence_tolerance
        if convergence_tolerance is not None
        else _default("convergence_tolerance", "mol/s")
    )
    states = _states(
        list(gas_components),
        float(gas_n.to("mol/s").magnitude),
        list(gas_z),
        float(gas_p.to("Pa").magnitude),
        float(gas_t.to("K").magnitude),
        list(liquid_components),
        float(liquid_n.to("mol/s").magnitude),
        list(liquid_z),
        float(liquid_p.to("Pa").magnitude),
        float(liquid_t.to("K").magnitude),
        list(transfer_components) if transfer_components else None,
        float(diameter.to("m").magnitude),
        float(height.to("m").magnitude),
        int(number_of_segments if number_of_segments is not None else _count("number_of_segments")),
        str(packing_type if packing_type is not None else DEFAULTS["packing_type"]),
        int(max_iterations if max_iterations is not None else _count("max_iterations")),
        float(gate.to("mol/s").magnitude),
        float(
            mass_transfer_correction
            if mass_transfer_correction is not None
            else DEFAULTS["mass_transfer_correction"]
        ),
        str(mass_transfer_correlation or DEFAULTS["mass_transfer_correlation"]),
        str(film_model or DEFAULTS["film_model"]),
        str(heat_transfer_model or DEFAULTS["heat_transfer_model"]),
        str(segment_solver or DEFAULTS["segment_solver"]),
        str(column_solver or DEFAULTS["column_solver"]),
    )
    return _result(states)


def _default(name: str, unit: str) -> Q:
    """The class's own default for a dimensioned parameter, as a quantity."""
    return quantity(float(DEFAULTS[name]), unit)


def _count(name: str) -> int:
    """The class's own default for a count."""
    return int(float(DEFAULTS[name]))


def _result(states: dict[str, Any]) -> RateBasedPackedColumnResult:
    """The flat record, from the profile's own states."""
    segments = states["segments"]
    column = {key: [segment[key] for segment in segments] for key in _SEGMENT_FIELDS}
    transfer_names = list(states["totals"])
    return RateBasedPackedColumnResult(
        gas_out_n=from_si(states["gas_out"]["n"], "mol/s"),
        gas_out_z=tuple(states["gas_out"]["z"]),
        gas_out_p=from_si(states["gas_out"]["p"], "Pa"),
        gas_out_t=from_si(states["gas_out"]["t"], "K"),
        gas_out_h=from_si(states["gas_out"]["h"], "J/mol"),
        liquid_out_n=from_si(states["liquid_out"]["n"], "mol/s"),
        liquid_out_z=tuple(states["liquid_out"]["z"]),
        liquid_out_p=from_si(states["liquid_out"]["p"], "Pa"),
        liquid_out_t=from_si(states["liquid_out"]["t"], "K"),
        liquid_out_h=from_si(states["liquid_out"]["h"], "J/mol"),
        iterations=states["iterations"],
        convergence_residual=from_si(states["residual"], "mol/s"),
        converged=states["converged"],
        total_absolute_molar_transfer=from_si(states["total_absolute"], "mol/s"),
        component_transfer_totals=tuple(
            from_si(states["totals"][name], "mol/s") for name in transfer_names
        ),
        transfer_components=tuple(transfer_names),
        segment_height_from_bottom=tuple(
            from_si(value, "m") for value in column["height_from_bottom"]
        ),
        segment_gas_temperature=tuple(from_si(value, "K") for value in column["gas_temperature"]),
        segment_liquid_temperature=tuple(
            from_si(value, "K") for value in column["liquid_temperature"]
        ),
        segment_gas_pressure=tuple(from_si(value, "Pa") for value in column["gas_pressure"]),
        segment_liquid_pressure=tuple(from_si(value, "Pa") for value in column["liquid_pressure"]),
        segment_gas_molar_flow=tuple(from_si(value, "mol/s") for value in column["gas_molar_flow"]),
        segment_liquid_molar_flow=tuple(
            from_si(value, "mol/s") for value in column["liquid_molar_flow"]
        ),
        segment_gas_density=tuple(from_si(value, "kg/m**3") for value in column["gas_density"]),
        segment_liquid_density=tuple(
            from_si(value, "kg/m**3") for value in column["liquid_density"]
        ),
        segment_gas_viscosity=tuple(from_si(value, "Pa*s") for value in column["gas_viscosity"]),
        segment_liquid_viscosity=tuple(
            from_si(value, "Pa*s") for value in column["liquid_viscosity"]
        ),
        segment_gas_diffusivity=tuple(
            from_si(value, "m**2/s") for value in column["gas_diffusivity"]
        ),
        segment_liquid_diffusivity=tuple(
            from_si(value, "m**2/s") for value in column["liquid_diffusivity"]
        ),
        segment_wetted_area=tuple(float(v) for v in column["wetted_area"]),
        segment_k_ga=tuple(float(v) for v in column["k_ga"]),
        segment_k_la=tuple(float(v) for v in column["k_la"]),
        segment_gas_heat_transfer_coefficient=tuple(
            float(v) for v in column["gas_heat_transfer_coefficient"]
        ),
        segment_liquid_heat_transfer_coefficient=tuple(
            float(v) for v in column["liquid_heat_transfer_coefficient"]
        ),
        segment_overall_heat_transfer_coefficient=tuple(
            float(v) for v in column["overall_heat_transfer_coefficient"]
        ),
        segment_interface_temperature=tuple(
            from_si(value, "K") for value in column["interface_temperature"]
        ),
        segment_heat_transfer_rate=tuple(
            from_si(value, "W") for value in column["heat_transfer_rate"]
        ),
        segment_pressure_drop_per_meter=tuple(
            from_si(value, "Pa") for value in column["pressure_drop_per_meter"]
        ),
        segment_percent_flood=tuple(float(v) for v in column["percent_flood"]),
        segment_net_molar_transfer=tuple(
            from_si(value, "mol/s") for value in column["net_molar_transfer"]
        ),
        segment_enthalpy_balance_residual=tuple(
            from_si(value, "W") for value in column["enthalpy_balance_residual"]
        ),
        warnings=_warnings(states),
    )


#: The per-segment fields the profile carries, which the record flattens into parallel vectors.
_SEGMENT_FIELDS = (
    "height_from_bottom",
    "gas_temperature",
    "liquid_temperature",
    "gas_pressure",
    "liquid_pressure",
    "gas_molar_flow",
    "liquid_molar_flow",
    "gas_density",
    "liquid_density",
    "gas_viscosity",
    "liquid_viscosity",
    "gas_diffusivity",
    "liquid_diffusivity",
    "wetted_area",
    "k_ga",
    "k_la",
    "gas_heat_transfer_coefficient",
    "liquid_heat_transfer_coefficient",
    "overall_heat_transfer_coefficient",
    "interface_temperature",
    "heat_transfer_rate",
    "pressure_drop_per_meter",
    "percent_flood",
    "net_molar_transfer",
    "enthalpy_balance_residual",
)


def _warnings(states: dict[str, Any]) -> tuple[Warning, ...]:
    """One caveat per substituted property, and one for a profile that missed its gate."""
    warnings: list[Warning] = []
    for key, property_name, field in (
        ("diffusivity", "diffusivity", "segment_gas_diffusivity"),
        ("surface_tension", "surface tension", "segment_wetted_area"),
        ("thermal_conductivity", "thermal conductivity", "segment_gas_heat_transfer_coefficient"),
        ("heat_capacity", "heat capacity", "segment_gas_heat_transfer_coefficient"),
    ):
        if states["fallbacks"][key]:
            warnings.append(
                Warning(
                    code=WarningCode.SOLVER_NOT_CONVERGED,
                    message=(
                        f"the {property_name} is the class's own `DEFAULT_*` constant: the "
                        "physical-property model answered nothing, and `RateBasedPackedColumn` "
                        "substitutes rather than refusing"
                    ),
                    field=field,
                )
            )
    if not states["converged"]:
        warnings.append(
            Warning(
                code=WarningCode.SOLVER_NOT_CONVERGED,
                message=(
                    f"the profile reached its {states['iterations']}-iteration cap without "
                    "meeting the gate, and its last iterate is published - which is what the "
                    "class does"
                ),
                field="max_iterations",
            )
        )
    return tuple(warnings)
