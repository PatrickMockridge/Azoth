"""``process.consequence_analysis`` - the QRA composer.

As NeqSim's ``ConsequenceAnalysisEngine`` computes it. Spec:
``specs/models/process/consequence_analysis.toml``

The Python twin of ``crates/azoth-process/src/models/consequence_analysis.rs``: the four
dose formulas over the ported consequence references, the zero-dose answer, the per-outcome
frequency and the summed individual risk.
"""

from __future__ import annotations

import math
from typing import Any

from azoth.core.errors import InvalidInputError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import ConsequenceAnalysisResult, ConsequenceOutcomeResult
from azoth.core.units import quantity
from azoth.core.warnings import Warning
from azoth.process.reference.gaussian_plume import concentration
from azoth.process.reference.jet_fire import jet_fire as jet_fire_model
from azoth.process.reference.pool_fire import pool_fire as pool_fire_model
from azoth.process.reference.probit import probit as probit_model
from azoth.process.reference.vce import vce as vce_model

MODEL_ID = "process.consequence_analysis"

#: The gas constant the class's toxic conversion uses, J/(mol K).
_TOXIC_GAS_CONSTANT = 8.314

_LABELS = {
    "jet_fire": "Jet fire",
    "pool_fire": "Pool fire",
    "vce": "VCE",
    "toxic_dispersion": "Toxic dispersion",
}

_ALLOWED = {
    "jet_fire": (
        "mass_flow",
        "heat_of_combustion",
        "radiative_fraction",
        "relative_humidity",
        "exposure_seconds",
    ),
    "pool_fire": ("pool_diameter", "burn_rate", "air_density", "exposure_seconds"),
    "vce": ("cloud_fuel_mass", "heat_of_combustion", "strength_class", "ambient_pressure"),
    "toxic_dispersion": (
        "emission_rate",
        "effective_stack_height",
        "wind_speed",
        "stability",
        "terrain",
        "molar_mass",
        "temperature",
        "pressure",
        "exposure_seconds",
    ),
}

_OPTIONAL_FIELDS = (
    "mass_flow",
    "heat_of_combustion",
    "radiative_fraction",
    "relative_humidity",
    "pool_diameter",
    "burn_rate",
    "air_density",
    "cloud_fuel_mass",
    "strength_class",
    "ambient_pressure",
    "emission_rate",
    "effective_stack_height",
    "wind_speed",
    "stability",
    "terrain",
    "molar_mass",
    "temperature",
    "pressure",
    "exposure_seconds",
)


def _required(index: int, kind: str, row: dict[str, Any], field: str) -> float:
    value = row.get(field)
    if value is None:
        raise InvalidInputError(f"outcomes[{index}].{field}", f"a {kind} branch requires `{field}`")
    return float(value)


def _require_finite(name: str, value: float) -> None:
    if not math.isfinite(value):
        raise InvalidInputError(
            name, "a non-finite value is not a number this calculation publishes"
        )


def _fatality(a: float, b: float, n: float, dose: float) -> float:
    """The probit of a dose; zero where the class's probit answers `-Infinity`."""
    if dose <= 0.0:
        return 0.0
    return float(probit_model(a, b, n, dose=dose).probability)


def _evaluate_row(
    index: int, row: dict[str, Any], release_frequency: float, receiver_distance: float
) -> tuple[str, str, float, float]:
    """One branch's label, type, frequency and fatality probability.

    The layer dumper calls this too, so the dose arithmetic is written once.
    """
    kind = str(row["outcome_type"])
    if kind not in _ALLOWED:
        raise InvalidInputError(
            f"outcomes[{index}].outcome_type",
            f"an outcome type is one of the class's four branches, and {kind!r} is not",
        )
    probability = float(row["conditional_probability"])
    _require_finite(f"outcomes[{index}].conditional_probability", probability)
    if not 0.0 <= probability <= 1.0:
        raise InvalidInputError(
            f"outcomes[{index}].conditional_probability",
            f"a branch probability is in [0, 1], and {probability} is not",
        )
    a = float(row["probit_a"])
    b = float(row["probit_b"])
    n = float(row["probit_n"])
    for field, value in (("probit_a", a), ("probit_b", b), ("probit_n", n)):
        _require_finite(f"outcomes[{index}].{field}", value)

    for field in _OPTIONAL_FIELDS:
        if row.get(field) is not None and field not in _ALLOWED[kind]:
            raise InvalidInputError(
                f"outcomes[{index}].{field}",
                f"a {kind} branch does not carry `{field}`",
            )

    if kind == "jet_fire":
        mass = _required(index, kind, row, "mass_flow")
        heat = _required(index, kind, row, "heat_of_combustion")
        fraction = _required(index, kind, row, "radiative_fraction")
        exposure = _required(index, kind, row, "exposure_seconds")
        humidity = row.get("relative_humidity")
        jet_result = jet_fire_model(
            quantity(mass, "kg/s"),
            heat,
            fraction,
            relative_humidity=None if humidity is None else float(humidity),
            distance=quantity(receiver_distance, "m"),
        )
        flux = (
            0.0
            if jet_result.incident_heat_flux is None
            else jet_result.incident_heat_flux.to("W/m**2").magnitude
        )
        dose = exposure * math.pow(flux, 4.0 / 3.0) / 1.0e4
    elif kind == "pool_fire":
        diameter = _required(index, kind, row, "pool_diameter")
        burn = _required(index, kind, row, "burn_rate")
        exposure = _required(index, kind, row, "exposure_seconds")
        air = row.get("air_density")
        pool_result = pool_fire_model(
            quantity(diameter, "m"),
            quantity(burn, "kg/(m**2*s)"),
            air_density=None if air is None else float(air),
            distance=quantity(receiver_distance, "m"),
        )
        flux = (
            0.0
            if pool_result.incident_heat_flux is None
            else pool_result.incident_heat_flux.to("W/m**2").magnitude
        )
        dose = exposure * math.pow(flux, 4.0 / 3.0) / 1.0e4
    elif kind == "vce":
        mass = _required(index, kind, row, "cloud_fuel_mass")
        heat = _required(index, kind, row, "heat_of_combustion")
        strength = _required(index, kind, row, "strength_class")
        ambient = row.get("ambient_pressure")
        vce_result = vce_model(
            quantity(mass, "kg"),
            heat,
            strength,
            ambient_pressure=None if ambient is None else quantity(float(ambient), "Pa"),
            distance=quantity(receiver_distance, "m"),
        )
        dose = (
            0.0
            if vce_result.overpressure is None
            else vce_result.overpressure.to("Pa").magnitude
        )
    else:
        emission = _required(index, kind, row, "emission_rate")
        stack = _required(index, kind, row, "effective_stack_height")
        wind = _required(index, kind, row, "wind_speed")
        stability = str(row["stability"])
        terrain = str(row["terrain"])
        molar = _required(index, kind, row, "molar_mass")
        temperature = _required(index, kind, row, "temperature")
        pressure = _required(index, kind, row, "pressure")
        exposure = _required(index, kind, row, "exposure_seconds")
        for field, value in (
            ("emission_rate", emission),
            ("wind_speed", wind),
            ("molar_mass", molar),
            ("temperature", temperature),
            ("pressure", pressure),
        ):
            _require_finite(f"outcomes[{index}].{field}", value)
        if molar <= 0.0:
            raise InvalidInputError(
                f"outcomes[{index}].molar_mass",
                "the ppm conversion divides by the molar mass",
            )
        if temperature <= 0.0 or pressure <= 0.0:
            raise InvalidInputError(
                f"outcomes[{index}].pressure",
                "the ppm conversion divides by the pressure and needs a positive temperature",
            )
        c = concentration(emission, stack, wind, stability, terrain, receiver_distance, 0.0, 0.0)
        ppm = (c / molar) * _TOXIC_GAS_CONSTANT * temperature / pressure * 1.0e6
        dose = ppm * exposure / 60.0

    frequency = release_frequency * probability
    return _LABELS[kind], kind, frequency, _fatality(a, b, n, dose)


def _spec() -> dict[str, Any]:
    """The generated spec, imported here so the module has no import-time cycle."""
    from azoth._models_gen import model

    return model(MODEL_ID)


def consequence_analysis(
    release_frequency_per_year: float,
    receiver_distance: Any,
    outcomes: list[dict[str, object]],
) -> ConsequenceAnalysisResult:
    """The class's four consequence branches over one receiver, and the summed risk.

    Args:
        release_frequency_per_year: the initiating release frequency in 1/year.
        receiver_distance: the distance from the release to the receiver.
        outcomes: the branch rows in caller order, each selecting its type by
            ``outcome_type`` and carrying that type's fields.

    Returns:
        The per-outcome rows and the sum of their fatality frequencies.

    Raises:
        InvalidInputError: for a probability outside ``[0, 1]``, a missing or non-positive
            field its type requires, a field belonging to another type, a non-positive
            receiver distance, a negative release frequency, or a non-finite input.
    """
    from azoth.core.units import input_records_to_si, input_to_si

    spec = _spec()
    checks = checks_for(spec)
    warnings: list[Warning] = []

    frequency = float(release_frequency_per_year)
    distance = input_to_si(spec, "receiver_distance", receiver_distance)
    apply_checks(
        checks.on_input,
        {
            "release_frequency_per_year": frequency,
            "receiver_distance": distance,
        }.get,
        warnings,
    )
    _require_finite("release_frequency_per_year", frequency)
    _require_finite("receiver_distance", distance)

    rows: list[dict[str, Any]] = input_records_to_si(spec, "outcomes", outcomes)
    results = []
    for index, row in enumerate(rows):
        label, kind, outcome_frequency, probability = _evaluate_row(index, row, frequency, distance)
        results.append(
            ConsequenceOutcomeResult(
                outcome=label,
                outcome_type=kind,
                frequency_per_year=outcome_frequency,
                fatality_probability=probability,
                fatality_frequency_per_year=outcome_frequency * probability,
            )
        )

    return ConsequenceAnalysisResult(
        outcomes=tuple(results),
        individual_fatality_risk_per_year=sum(row.fatality_frequency_per_year for row in results),
        warnings=tuple(warnings),
    )
