"""The fluid property provider interface.

Deliberately small: two methods, because the hydraulics slice needs density and
viscosity and nothing else. Speed of sound is *not* here even though the
Darcy-Weisbach spec's Mach assumption would need it - a spec-level assumption
recorded as unchecked is honest, whereas an unused interface method with
half-sourced data behind it is not.

The interface takes a pressure so that a future provider can honour it. The two
built-ins ship data at 1 atm only and say so rather than silently returning a
value for a pressure they were never tabulated at - which, for a gas, would be a
wrong answer of exactly the kind this library exists to prevent.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from itertools import pairwise
from pathlib import Path
from typing import Protocol, runtime_checkable

from azoth._data import find
from azoth.core.errors import InvalidInputError, PropertyUnavailableError
from azoth.core.units import Q, quantity, to_si

#: The pressure the built-in tables were tabulated at.
STANDARD_PRESSURE_PA = 101_325.0

#: Tolerance for treating a supplied pressure as atmospheric. Half a percent is
#: generous enough for "1 atm, roughly" and tight enough that a real process
#: pressure is rejected.
_PRESSURE_TOLERANCE = 0.005


@dataclass(frozen=True, slots=True)
class PropertyPoint:
    """One tabulated row."""

    temperature_c: float
    density_kg_m3: float
    dynamic_viscosity_pa_s: float
    citation: str
    verify_status: str
    source_ref: str | None = None
    source_locator: str | None = None


@runtime_checkable
class PropertyProvider(Protocol):
    """What a fluid property source has to be able to do."""

    @property
    def name(self) -> str:
        """Identifier, e.g. ``water``."""
        ...

    def density(self, temperature: Q, pressure: Q | None = None) -> Q:
        """Density at a temperature, and optionally a pressure."""
        ...

    def dynamic_viscosity(self, temperature: Q, pressure: Q | None = None) -> Q:
        """Dynamic viscosity at a temperature, and optionally a pressure."""
        ...


class TableProvider:
    """A provider backed by a small table of tabulated values.

    Linear interpolation between tabulated points, exact at them. Outside the
    tabulated range it raises rather than extrapolating: water's viscosity
    changes by a factor of six across 0-100 C, so a straight-line extension past
    either end would be a confident wrong number, which is the failure mode this
    whole library is organised against.
    """

    def __init__(self, name: str, path: str) -> None:
        self._name = name
        self._path = path
        self._points = _read_table(find(path))

    @property
    def name(self) -> str:
        """Identifier, e.g. ``water``."""
        return self._name

    @property
    def points(self) -> tuple[PropertyPoint, ...]:
        """The tabulated rows, for diagnostics and tests."""
        return self._points

    @property
    def temperature_range_c(self) -> tuple[float, float]:
        """The tabulated temperature range, in degrees Celsius."""
        return (self._points[0].temperature_c, self._points[-1].temperature_c)

    def _interpolate(self, temperature_c: float, attribute: str) -> float:
        points = self._points
        low, high = self.temperature_range_c
        if not low <= temperature_c <= high:
            raise PropertyUnavailableError(
                self._name,
                attribute,
                f"temperature {temperature_c} C is outside the tabulated range "
                f"{low}-{high} C; this provider does not extrapolate",
            )

        # Exact hits first, so a tabulated point returns the tabulated value
        # rather than an interpolation that happens to land on it.
        for point in points:
            if temperature_c == point.temperature_c:
                return float(getattr(point, attribute))

        for lower, upper in pairwise(points):
            if lower.temperature_c < temperature_c < upper.temperature_c:
                span = upper.temperature_c - lower.temperature_c
                fraction = (temperature_c - lower.temperature_c) / span
                lo = float(getattr(lower, attribute))
                hi = float(getattr(upper, attribute))
                return lo + fraction * (hi - lo)

        raise PropertyUnavailableError(  # pragma: no cover - guarded by the range check
            self._name, attribute, f"no bracketing points for {temperature_c} C"
        )

    def _check_pressure(self, pressure: Q | None) -> None:
        if pressure is None:
            return
        pascals = to_si(pressure, "Pa", "pressure")
        if abs(pascals - STANDARD_PRESSURE_PA) / STANDARD_PRESSURE_PA > _PRESSURE_TOLERANCE:
            raise PropertyUnavailableError(
                self._name,
                "pressure",
                f"the built-in table is tabulated at 1 atm only; {pascals:.0f} Pa "
                f"was requested. Supply your own provider for other pressures.",
            )

    def _temperature_c(self, temperature: Q) -> float:
        kelvin = to_si(temperature, "K", "temperature")
        return kelvin - 273.15

    def density(self, temperature: Q, pressure: Q | None = None) -> Q:
        """Density at a temperature, and optionally a pressure."""
        self._check_pressure(pressure)
        value = self._interpolate(self._temperature_c(temperature), "density_kg_m3")
        return quantity(value, "kg/m**3")

    def dynamic_viscosity(self, temperature: Q, pressure: Q | None = None) -> Q:
        """Dynamic viscosity at a temperature, and optionally a pressure."""
        self._check_pressure(pressure)
        value = self._interpolate(self._temperature_c(temperature), "dynamic_viscosity_pa_s")
        return quantity(value, "Pa*s")


def _read_table(path: Path) -> tuple[PropertyPoint, ...]:
    lines = [
        line
        for line in path.read_text(encoding="utf-8").splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ]
    if len(lines) < 2:  # pragma: no cover - guarded by the test suite
        raise InvalidInputError("fluids", f"{path} has no data rows")

    points: list[PropertyPoint] = []
    for raw in csv.DictReader(lines):
        try:
            points.append(
                PropertyPoint(
                    temperature_c=float(raw["temperature_c"]),
                    density_kg_m3=float(raw["density_kg_m3"]),
                    dynamic_viscosity_pa_s=float(raw["dynamic_viscosity_pa_s"]),
                    citation=raw["citation"],
                    verify_status=raw["verify_status"],
                    # Empty means absent, not an empty string - the same rule the
                    # Rust loader applies, so the two agree on what a row with no
                    # source carries.
                    source_ref=raw.get("source_ref") or None,
                    source_locator=raw.get("source_locator") or None,
                )
            )
        except (KeyError, ValueError) as exc:
            raise InvalidInputError("fluids", f"malformed row in {path}: {exc}") from exc

    points.sort(key=lambda p: p.temperature_c)
    return tuple(points)
