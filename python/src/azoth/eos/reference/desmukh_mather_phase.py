"""``eos.desmukh_mather_phase`` - the activity coefficients of a Desmukh-Mather phase.

Spec: ``specs/models/eos/desmukh_mather_phase.toml``. A *direct* model: no iteration, so no
algorithm block.

This is the pure-Python reference: a second, independent expression of the same physics as
the Rust kernel, held to it by ``test_cross_impl.py``. What is here is the branch structure,
which is where the two models of this pair differ: the solvent is chosen by
``REFERENCESTATETYPE`` and the Poynting correction is chosen by the *name* ``water``, so
they are two questions and not one.
"""

from __future__ import annotations

import math
from collections.abc import Sequence

from azoth import _models_gen
from azoth.core.errors import InvalidInputError, PropertyUnavailableError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import DesmukhMatherPhaseResult
from azoth.core.units import Q, input_to_si
from azoth.core.warnings import Warning
from azoth.eos import components as _components
from azoth.eos.reference import _henry
from azoth.eos.reference.antoine_vapor_pressure import saturation_pressure

MODEL_ID = "eos.desmukh_mather_phase"

#: The Debye-Huckel `A` of the expression, dimensionless.
A = 1.174

#: The `B` that scales the ionic diameter, in reciprocal metres.
B = 3.32384e9

#: The gas constant the Poynting correction uses, J/(mol K). NeqSim's ``Component.R``.
R = 8.314462618

#: The ``aij`` of the ``MDEA``/``MDEA`` diagonal, which is **not in the interaction table**.
MDEA_DIAGONAL = -0.0828487

#: The fugacity coefficient an ion gets. The smallest of the tranche's three constants.
INSOLUBLE_ION = 1.0e-15


def desmukh_mather_phase(
    components: Sequence[str],
    T: Q,
    P: Q,
    x: Sequence[float],
) -> DesmukhMatherPhaseResult:
    """The activity coefficients of a phase whose non-ideality is Desmukh and Mather's.

    ``ln gamma*_i = -A z_i^2 sqrt(I)/(1 + B d_i 1e-10 sqrt(I)) + 2 sum_j beta_ij m_j`` and
    ``gamma_i = m_i M_solvent exp(ln gamma*_i) / x_i``, with ``beta_ij = aij + bij T`` read
    from the interaction table. The fugacity coefficient is NeqSim's
    ``(gamma / gamma^infinity) H / P`` branch.

    Args:
        components: the substances the phase is made of, by name.
        T: absolute temperature.
        P: absolute pressure.
        x: the phase's mole fractions; non-negative and summing to one, with a
            ``solvent``-reference component among them.

    Returns:
        The activity coefficients, the molalities, the ionic strength, the mean solvent
        molar mass and the fugacity coefficients.

    Raises:
        InvalidInputError: if a component is not in the databank, if ``x`` is not a
            composition, or if the mixture carries no ``solvent``-reference component.
        PropertyUnavailableError: if a ``solvent`` component carries no Antoine correlation
            or a neutral solute one carries no Henry row.
        OutOfRangeError: if ``T`` or ``P`` is not positive.

    Example:
        >>> import azoth
        >>> q = azoth.ureg.Quantity
        >>> r = desmukh_mather_phase(
        ...     ["water", "na+", "cl-", "co2"], q(313.15, "K"), q(500000.0, "Pa"),
        ...     [0.89, 0.04, 0.04, 0.03],
        ... )
        >>> round(r.ionic_strength, 9)
        2.494799901
    """
    spec = _models_gen.model(MODEL_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    t_si = input_to_si(spec, "T", T)
    p_si = input_to_si(spec, "P", P)
    apply_checks(checks.on_input, {"T": t_si, "P": p_si}.get, warnings)

    n = len(components)
    if len(x) != n:
        raise InvalidInputError("x", f"a mixture of {n} components has {len(x)} mole fractions")
    bad = next((i for i, value in enumerate(x) if value < 0.0), None)
    if bad is not None:
        raise InvalidInputError("x", f"x[{bad}] is {x[bad]} but a mole fraction cannot be negative")
    total = sum(x)
    if abs(total - 1.0) > 1.0e-9:
        raise InvalidInputError(
            "x",
            f"the mole fractions sum to {total}, not to one. Renormalising them here "
            "would make a composition error invisible in every number downstream, so it "
            "is refused instead",
        )

    entries = [_components.entry(name) for name in components]
    names = [entry.name for entry in entries]
    charge = [entry.ionic_charge for entry in entries]
    diameter = [entry.deshmukh_mather_diameter for entry in entries]

    # **The solvent is whatever `REFERENCESTATETYPE` says it is**, which is NeqSim's rule
    # here and not `PhasePitzer`'s. An amine tagged `solvent` joins water in the sum.
    solvent = [i for i, entry in enumerate(entries) if entry.reference_state == _components.SOLVENT]
    solvent_weight = 0.0
    for i in solvent:
        mass = entries[i].molar_mass
        if mass is not None:
            solvent_weight += x[i] * mass.to_base_units().magnitude
    solvent_moles = sum(x[i] for i in solvent)
    if not solvent_weight > 0.0 or not solvent_moles > 0.0:
        raise InvalidInputError(
            "components",
            "the mixture carries no `solvent`-reference component, so its molalities would "
            "all be zero. NeqSim divides by a zero solvent weight and returns a phase of "
            "nothing",
        )
    solvent_molar_mass = solvent_weight / solvent_moles

    molality = [xi / solvent_weight for xi in x]
    ionic_strength = 0.5 * sum(m * z * z for m, z in zip(molality, charge, strict=True))

    ln_gamma: list[float] = []
    for i in range(n):
        star = _ln_gamma_star(names, molality, charge, diameter, t_si, i)
        value = max(molality[i], 0.0) * solvent_molar_mass * math.exp(star) / max(x[i], 1.0e-30)
        # **NeqSim's own fallback**: `gamma` is set to `exp(ln gamma)` when the conversion
        # is not finite or not positive, which never happens for a present component.
        ln_gamma.append(math.log(value) if value > 0.0 and math.isfinite(value) else star)

    ln_phi: list[float] = []
    for i, entry in enumerate(entries):
        is_solvent = entry.reference_state == _components.SOLVENT
        neutral = entry.ionic_charge == 0.0
        # **Four branches, and the first is a name test.** `ComponentDesmukhMather.fugcoef`
        # asks for `water` *by name* before it asks what the reference state is, so a
        # solvent-reference component that is not water gets no Poynting correction.
        if entry.name.lower() == "water":
            p0 = saturation_pressure(entry, t_si, warnings)
            mass = entry.molar_mass
            molar_volume = 1.0e-3 * (mass.to_base_units().magnitude if mass is not None else 0.0)
            poynting = math.exp(molar_volume / (R * t_si) * (p_si - p0))
            coefficient = math.exp(ln_gamma[i]) * p0 / p_si * poynting
        elif neutral and is_solvent:
            p0 = saturation_pressure(entry, t_si, warnings)
            coefficient = math.exp(ln_gamma[i]) * p0 / p_si
        elif neutral and entry.reference_state == _components.SOLUTE:
            infinite = _ln_gamma_infinite_dilution(
                entry.name, entry.ionic_charge, diameter[i], t_si
            )
            coefficient = (
                _henry.effective_coefficient(entry, t_si)
                * 1.0e5
                / p_si
                * math.exp(ln_gamma[i] - infinite)
            )
        else:
            coefficient = INSOLUBLE_ION
        if not coefficient > 0.0 or not math.isfinite(coefficient):
            raise PropertyUnavailableError(
                entry.name,
                "fugacity coefficient",
                f"the branch its reference state selects gives {coefficient}, which is not a "
                f"coefficient. NeqSim returns it without comment",
            )
        ln_phi.append(math.log(coefficient))

    return DesmukhMatherPhaseResult(
        gamma=tuple(math.exp(value) for value in ln_gamma),
        ln_gamma=tuple(ln_gamma),
        molality=tuple(molality),
        ionic_strength=ionic_strength,
        solvent_molar_mass=solvent_molar_mass,
        ln_phi=tuple(ln_phi),
        warnings=tuple(warnings),
    )


def _ln_gamma_star(
    names: Sequence[str],
    molality: Sequence[float],
    charge: Sequence[float],
    diameter: Sequence[float],
    temperature_k: float,
    component: int,
) -> float:
    """``ln gamma*``, the expression before the mole-fraction conversion."""
    ionic_strength = 0.5 * sum(m * z * z for m, z in zip(molality, charge, strict=True))
    sqrt_i = math.sqrt(ionic_strength)
    debye_huckel = (
        -A * charge[component] ** 2 * sqrt_i / (1.0 + B * diameter[component] * 1.0e-10 * sqrt_i)
    )
    pairs = 0.0
    for j, m in enumerate(molality):
        # **NeqSim's loop skips neither the diagonal nor water**: it sums every component
        # and relies on the diagonal's `beta` being zero, and `ComponentDesmukhMather`
        # tests for water by name before adding its term.
        if j == component or names[j].lower() == "water":
            continue
        pairs += 2.0 * _beta(names[component], names[j], temperature_k) * m
    return debye_huckel + pairs


def _beta(first: str, second: str, temperature_k: float) -> float:
    """``beta_ij = aij + bij T``, either order round."""
    if first.lower() == "mdea" and second.lower() == "mdea":
        return MDEA_DIAGONAL
    record = _components.desmukh_mather_pair(first, second)
    if record is None:
        return 0.0
    aij, bij = record
    return aij + bij * temperature_k


def _ln_gamma_infinite_dilution(
    name: str, charge: float, diameter: float, temperature_k: float
) -> float:
    """``gamma^infinity``: the same expression at ``Phase.initRefPhases``'s reference state.

    That phase is the solute at ``1e-10`` mol and **water at 10 mol**, so its solute mole
    fraction is ``1e-11`` and its solvent weight is ``10 M_water``.
    """
    water = _components.WATER_MOLAR_MASS
    weight = 10.0 * water
    names = [name, "water"]
    molality = [1.0e-10 / weight, 10.0 / weight]
    charges = [charge, 0.0]
    diameters = [diameter, 0.0]
    star = _ln_gamma_star(names, molality, charges, diameters, temperature_k, 0)
    reference_fraction = 1.0e-10 / (10.0 + 1.0e-10)
    value = molality[0] * water * math.exp(star) / reference_fraction
    return math.log(value) if value > 0.0 and math.isfinite(value) else star
