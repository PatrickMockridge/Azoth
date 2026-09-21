"""The Furst short-range ``Wij`` table: the Python reference.

The mirror of ``crates/azoth-eos/src/furst_mixing.rs``. Three sweeps over the components,
each overwriting the last - a table read, a cation sweep and a gas-then-oil sweep - so a
methane/Na+ pair is the gas-ion value and a water/Na+ pair is the correlation's, and neither
is visible from the other's table alone.
"""

from __future__ import annotations

import math
from collections.abc import Callable, Sequence
from dataclasses import dataclass

from azoth.core.errors import InvalidInputError

#: The reference temperature ``Wij(T)`` is written against.
T_REFERENCE = 298.15

#: The reference dielectric constant the predictive correlation is centred on.
EPSILON_WATER_REFERENCE = 78.4


@dataclass(frozen=True, slots=True)
class FurstComponent:
    """One component, as the construction sees it."""

    #: The databank name, lower-cased. The sweeps test it, so it is arithmetic.
    name: str
    #: The ionic charge, in elementary charges.
    charge: float
    #: The **table's** diameter in ångströms, which is the ionic diameter the correlation
    #: uses - not the one derived from the fitted covolume, which the phase's own terms use.
    diameter: float
    #: The dielectric constant at 298.15 K, for the predictive correlation.
    dielectric_at_reference: float

    @property
    def is_cation(self) -> bool:
        """``charge > 0`` exactly, as NeqSim writes it."""
        return self.charge > 0.0

    @property
    def is_anion(self) -> bool:
        """``charge < -0.01``, NeqSim's own threshold."""
        return self.charge < -0.01


@dataclass(frozen=True, slots=True)
class WijTable:
    """``w0``, ``w1`` and ``w2`` for every ordered pair, flattened row-major."""

    w0: tuple[float, ...]
    w1: tuple[float, ...]
    w2: tuple[float, ...]
    n: int

    def wij(self, i: int, j: int, temperature: float) -> float:
        """``Wij(T)`` for one pair.

        **Total by construction.** ``temperature`` is an absolute temperature, refused at the
        model's boundary when it is not positive, so ``1/T`` and ``ln(T/298.15)`` are both
        defined. The ``log`` is NeqSim's own form of the temperature correction and is kept
        rather than rewritten as ``ln T - ln 298.15``: the two differ in the last digits, and
        this is a layer the probe prints.
        """
        at = i * self.n + j
        return (
            self.w0[at]
            + self.w1[at] * (1.0 / temperature - 1.0 / T_REFERENCE)
            + self.w2[at]
            * ((T_REFERENCE - temperature) / temperature + math.log(temperature / T_REFERENCE))
        )

    def wij_dt(self, i: int, j: int, temperature: float) -> float:
        """``dWij/dT`` for one pair."""
        at = i * self.n + j
        return (
            -self.w1[at] / temperature**2
            - self.w2[at] * (T_REFERENCE - temperature) / temperature**2
        )

    def wij_dtdt(self, i: int, j: int, temperature: float) -> float:
        """``d^2Wij/dT^2`` for one pair."""
        at = i * self.n + j
        return (
            2.0 * self.w1[at] / temperature**3
            + self.w2[at] / temperature**2
            + 2.0 * self.w2[at] * (T_REFERENCE - temperature) / temperature**3
        )


@dataclass(frozen=True, slots=True)
class ShortRange:
    """The three sums the short-range term reads."""

    w: float
    w_dt: float
    w_dtdt: float


def _cpa(index: int) -> float:
    """``furstParamsCPA[index]``, a missing coefficient read as zero."""
    return _named("furstParamsCPA", index)


def _t_dep(index: int) -> float:
    """``furstParamsCPA_TDep[index]``."""
    return _named("furstParamsCPA_TDep", index)


def _named(set_name: str, index: int) -> float:
    # Imported here rather than at module scope: `components` imports this module for the
    # Fürst term it resolves, so a top-level import would close a cycle.
    from azoth.eos.components import furst_parameters

    values = furst_parameters(set_name)
    return values[index] if index < len(values) else 0.0


def _glycol_fit(set_name: str, diameter: float, divalent: bool) -> float:
    """One glycol's fit, ``slope * d + intercept`` from the named set.

    The two sets differ at ``[2]``, ``[3]``, ``[6]`` and ``[7]`` - ``4.98e-5`` against
    ``8.0e-5`` at ``[2]`` - so a ``TEG`` pair given the MEG set carries ``2.11x`` the ``Wij``
    its own fit gives. Each glycol reads its own, which is NeqSim's dispatch as of issue 3846's
    fix (#3847).
    """
    if divalent:
        return _named(set_name, 6) * diameter + _named(set_name, 7)
    return _named(set_name, 2) * diameter + _named(set_name, 3)


def _predictive(solvent_dielectric: float, diameter: float, divalent: bool) -> float:
    """``getPredictiveWij``: a linear correction to the water fit in ``1/eps - 1/eps_w``."""
    delta = 1.0 / solvent_dielectric - 1.0 / EPSILON_WATER_REFERENCE
    slope_water = _cpa(6) if divalent else _cpa(2)
    intercept_water = _cpa(7) if divalent else _cpa(3)
    slope_coefficient = -2.5e-3 if divalent else -3.0e-3
    intercept_coefficient = 1.5e-2 if divalent else 1.2e-2
    return (slope_water + slope_coefficient * delta) * diameter + (
        intercept_water + intercept_coefficient * delta
    )


def _cation_anion(cation: float, anion: float, divalent: bool) -> tuple[float, float, float]:
    """``Wij(cation, anion) = p4 (d + d)^4 + p5`` on the monovalent or divalent set."""
    sum4 = (cation + anion) ** 4
    if divalent:
        return (
            _cpa(8) * sum4 + _cpa(9),
            _t_dep(12) * sum4 + _t_dep(13),
            _t_dep(14) * sum4 + _t_dep(15),
        )
    return (
        _cpa(4) * sum4 + _cpa(5),
        _t_dep(4) * sum4 + _t_dep(5),
        _t_dep(6) * sum4 + _t_dep(7),
    )


def _cation_solvent(cation: FurstComponent, solvent: FurstComponent, divalent: bool) -> float:
    """``Wij(cation, neutral)``, by the neutral's name."""
    d = cation.diameter
    name = solvent.name
    if name == "water":
        return _cpa(6) * d + _cpa(7) if divalent else _cpa(2) * d + _cpa(3)
    # **Each glycol reads its own fit.** NeqSim's chain tested `TEG` twice and the first match
    # won, so a TEG pair was given the MEG set - issue 3846, closed upstream in #3847 by the
    # swap this now makes.
    if name in ("teg", "triethylene glycol"):
        return _glycol_fit("furstParamsCPA_TEG", d, divalent)
    if name in ("meg", "ethylene glycol"):
        return _glycol_fit("furstParamsCPA_MEG", d, divalent)
    if name == "mdea":
        return _named("furstParams", 6) if divalent else _named("furstParams", 2)
    if name == "piperazine":
        return _cpa(6) * d + _cpa(7) if divalent else _cpa(2) * d + _cpa(3)
    if name == "methanol":
        index = 6 if divalent else 2
        return _named("furstParamsCPA_MeOH", index)
    if name == "ethanol":
        index = 6 if divalent else 2
        return _named("furstParamsCPA_EtOH", index)
    if name == "mea":
        if divalent:
            return _named("furstParamsCPA_MEA", 6) * d + _named("furstParamsCPA_MEA", 7)
        return _named("furstParamsCPA_MEA", 2) * d + _named("furstParamsCPA_MEA", 3)
    return _predictive(solvent.dielectric_at_reference, d, divalent)


#: The ``furstParamsGasIon`` index pair for a non-polar gas, by name.
_GAS_INDICES: dict[str, tuple[int, int]] = {
    "co2": (0, 1),
    "carbon dioxide": (0, 1),
    "methane": (2, 3),
    "ch4": (2, 3),
    "ethane": (4, 5),
    "c2h6": (4, 5),
    "propane": (6, 7),
    "c3h8": (6, 7),
    "n-butane": (8, 9),
    "i-butane": (8, 9),
    "n-c4": (8, 9),
    "i-c4": (8, 9),
    "butane": (8, 9),
    "nitrogen": (12, 13),
    "n2": (12, 13),
    "h2s": (14, 15),
    "hydrogen sulfide": (14, 15),
    "hydrogensulfide": (14, 15),
    "hydrogen": (16, 17),
    "h2": (16, 17),
}

#: The heavier components the gas pass names as ``C5+``.
_C5_PLUS = (
    "n-pentane",
    "i-pentane",
    "n-c5",
    "i-c5",
    "pentane",
    "n-hexane",
    "hexane",
    "n-c6",
    "n-heptane",
    "heptane",
    "n-c7",
    "n-octane",
    "octane",
    "n-c8",
    "n-nonane",
    "nonane",
    "n-c9",
    "n-decane",
    "decane",
    "n-c10",
)

#: The ``furstParamsOIIon`` index pair for an organic inhibitor, by name.
_OIL_INDICES: dict[str, tuple[int, int]] = {
    "methanol": (0, 1),
    "meoh": (0, 1),
    "meg": (2, 3),
    "ethylene glycol": (2, 3),
    "ethanol": (4, 5),
    "etoh": (4, 5),
}


def gas_indices(name: str) -> tuple[int, int] | None:
    """The gas-ion index pair for a non-polar gas, by name."""
    if name in _C5_PLUS:
        return (10, 11)
    return _GAS_INDICES.get(name)


def oil_indices(name: str) -> tuple[int, int] | None:
    """The oil-ion index pair for an organic inhibitor, by name."""
    return _OIL_INDICES.get(name)


def wij_table(
    components: Sequence[FurstComponent],
    fitted: Callable[[str, str], tuple[float, float, float, bool] | None],
) -> WijTable:
    """Build the pair table for a phase, the way ``ElectrolyteMixRule.calcWij`` builds it.

    Raises:
        InvalidInputError: if ``components`` is empty.
    """
    n = len(components)
    if n == 0:
        raise InvalidInputError(
            "components",
            "the short-range table is a sum over pairs and there are no components",
        )
    w0 = [0.0] * (n * n)
    w1 = [0.0] * (n * n)
    w2 = [0.0] * (n * n)
    calc_wij = [False] * (n * n)
    for i in range(n):
        for j in range(n):
            record = fitted(components[i].name, components[j].name)
            if record is None:
                continue
            w0[i * n + j], w1[i * n + j], w2[i * n + j], calc_wij[i * n + j] = record

    def set_symmetric(i: int, j: int, a: float, b: float, c: float) -> None:
        for x, y in ((i, j), (j, i)):
            w0[x * n + y] = a
            w1[x * n + y] = b
            w2[x * n + y] = c

    def set_w0(i: int, j: int, value: float) -> None:
        """Write the constant term alone, leaving the temperature coefficients.

        The gas and oil sweeps do this, so a pair the cation sweep gave a ``w1`` and ``w2``
        keeps them. Zeroing them there is the natural thing to write and makes ``dW/dT``
        short by 0.02%, which reads as rounding.
        """
        w0[i * n + j] = value
        w0[j * n + i] = value

    # The cation sweep.
    for i in range(n):
        if not components[i].is_cation:
            continue
        divalent = round(components[i].charge) >= 2
        for j in range(n):
            if calc_wij[i * n + j]:
                continue
            if components[j].is_anion:
                set_symmetric(
                    i, j, *_cation_anion(components[i].diameter, components[j].diameter, divalent)
                )
            elif components[j].charge == 0.0:
                value = _cation_solvent(components[i], components[j], divalent)
                d = components[i].diameter
                if components[i].name == "mdea+":
                    coeffs = (0.0, 0.0)
                elif divalent:
                    coeffs = (_t_dep(8) * d + _t_dep(9), _t_dep(10) * d + _t_dep(11))
                else:
                    coeffs = (_t_dep(0) * d + _t_dep(1), _t_dep(2) * d + _t_dep(3))
                set_symmetric(i, j, value, coeffs[0], coeffs[1])

    # The gas sweep, then the oil sweep.
    for i in range(n):
        pair = gas_indices(components[i].name)
        if pair is None:
            continue
        for j in range(n):
            if calc_wij[i * n + j]:
                continue
            if components[j].is_cation:
                set_w0(i, j, _named("furstParamsGasIon", pair[0]))
            elif components[j].is_anion:
                set_w0(i, j, _named("furstParamsGasIon", pair[1]))
    for i in range(n):
        pair = oil_indices(components[i].name)
        if pair is None:
            continue
        for j in range(n):
            if calc_wij[i * n + j]:
                continue
            side = 0 if components[j].is_cation else (1 if components[j].is_anion else None)
            if side is None:
                continue
            value = _named("furstParamsOIIon", pair[side])
            # **Only a non-zero parameter is applied**, so a zero here leaves the cation
            # sweep's value in place rather than overwriting it with nothing.
            if abs(value) > 1.0e-20:
                set_w0(i, j, value)

    return WijTable(tuple(w0), tuple(w1), tuple(w2), n)


def short_range(table: WijTable, mole_numbers: Sequence[float], temperature: float) -> ShortRange:
    """The three sums at a composition.

    Raises:
        InvalidInputError: if the mole numbers do not match the table.
    """
    if len(mole_numbers) != table.n:
        raise InvalidInputError(
            "mole_numbers",
            f"the table is {table.n} x {table.n} and there are {len(mole_numbers)} mole numbers",
        )
    w = w_dt = w_dtdt = 0.0
    for i in range(table.n):
        for j in range(table.n):
            weight = mole_numbers[i] * mole_numbers[j]
            w += weight * table.wij(i, j, temperature)
            w_dt += weight * table.wij_dt(i, j, temperature)
            w_dtdt += weight * table.wij_dtdt(i, j, temperature)
    return ShortRange(w=-w, w_dt=-w_dt, w_dtdt=-w_dtdt)
