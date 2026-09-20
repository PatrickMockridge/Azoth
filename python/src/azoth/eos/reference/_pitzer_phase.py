"""Pitzer's activity-coefficient arithmetic, as NeqSim's `ComponentGePitzer` carries it.

The Python mirror of ``crates/azoth-eos/src/pitzer_phase.rs``, held to the same constants
by the tests on both sides. The catalogue and the dataset selection are
``_pitzer_catalog``'s; this is the arithmetic that reads them.

**The ion branch is one long expression and the three pieces that are easy to lose are
the ones the flags gate.** ``fBprime`` is not accumulated when PHREEQC's common-ion terms
are active, because :func:`common_ion_contribution` computes its own over every pair;
``beta2`` applies only for a 2:2 pair or when the dataset activated it elsewhere; and
``E_theta`` enters twice, once through ``F`` and once through the theta sum.

``ln(gamma_M) = z_M^2 F + sum_a m_a (2 B_Ma + Z C_Ma)`` with
``F = -A_phi [ sqrt(I)/(1 + b sqrt(I)) + (2/b) ln(1 + b sqrt(I)) ] + sum_c sum_a m_c m_a B'_ca``,
NeqSim's own comment.
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from typing import Any

from azoth.eos.reference import _electrolyte as electrolyte
from azoth.eos.reference import _pitzer_electrostatic as electrostatic
from azoth.eos.reference._pitzer_catalog import find

#: The gas constant is not used here; the Debye-Huckel parameter carries its own arithmetic.


def debye_huckel_a_phi(temperature_k: float) -> float:
    """The Debye-Huckel parameter ``A_phi(T)``, in ``(mol/kg)^-1/2``.

    ``ComponentGePitzer.debyeHuckelAphi``, divided by the three its callers divide by.
    Water density from Kell (1975) with the IAPWS-style branch above 100 C and a 700 kg/m3
    floor; the dielectric constant from Archer and Wang (1990) floored at 20; then
    ``1.4006e6 sqrt(rho_g_per_cm3) / (eps T)^1.5``.

    **There is no clamp below 0 C** - the Kell polynomial is extrapolated - and the
    density floor is what keeps the cold end from going negative under the square root.
    """
    celsius = temperature_k - 273.15
    density = (
        999.83
        + 5.0948e-2 * celsius
        - 7.5722e-3 * celsius**2
        + 3.8907e-5 * celsius**3
        - 1.2e-7 * celsius**4
    )
    if celsius > 100.0:
        excess = celsius - 100.0
        density = 958.0 - 1.08 * excess - 0.0028 * excess**2
    density = max(density, 700.0)
    per_cm3 = density / 1000.0

    dielectric = 87.740 - 0.40008 * celsius + 9.398e-4 * celsius**2 - 1.410e-6 * celsius**3
    dielectric = max(dielectric, 20.0)
    # **`math.pow` and not `**`, because `**` is typed `Any` and is not.** typeshed
    # returns `Any` from `float ** float` for the reason `docs/src/calculus/numerics.md`
    # gives: a negative base with a fractional exponent promotes to a *complex* in Python.
    # The base here is `dielectric * T` with the dielectric floored at 20 and a positive
    # temperature, so it is positive by construction and `math.pow` is total - which is
    # also the form whose return type is `float`, so the type system can see it.
    return 1.4006e6 * math.sqrt(per_cm3) / math.pow(dielectric * temperature_k, 1.5)


def _g(x: float) -> float:
    """``2 (1 - (1 + x) exp(-x)) / x^2``, zero below the floor as NeqSim has it."""
    if x <= DERIVATIVE_FLOOR:
        return 0.0
    return 2.0 * (1.0 - (1.0 + x) * math.exp(-x)) / (x * x)


def _g_prime(x: float) -> float:
    """The derivative of :func:`_g`, zero below the same floor."""
    if x <= DERIVATIVE_FLOOR:
        return 0.0
    return -2.0 * (1.0 - (1.0 + x + x * x / 2.0) * math.exp(-x)) / (x * x)


#: The temperature both of NeqSim's Pitzer forms state their value at, in K.
#:
#: ``PitzerParameterDatasets.PHREEQC_REFERENCE_TEMPERATURE_K`` and the ``298.15`` the
#: Silvester form divides by are the same number, which is why one constant serves both.
REFERENCE_TEMPERATURE_K = 298.15

#: Within this many kelvin of the reference, the catalogue form returns its constant term
#: unchanged rather than evaluating the polynomial (``PitzerTemperatureFunction``).
REFERENCE_TOLERANCE_K = 1.0e-3


def catalog_value(a: Sequence[float], temperature_k: float) -> float:
    """`PitzerTemperatureFunction.valueAt`, the PHREEQC six-coefficient form.

    ``a0 + a1(1/T - 1/Tr) + a2 ln(T/Tr) + a3(T - Tr) + a4(T^2 - Tr^2) + a5(1/T^2 - 1/Tr^2)``

    **The reference window is a branch and not a rounding guard.** Within
    :data:`REFERENCE_TOLERANCE_K` of the reference the constant term is returned exactly,
    so ``value_at(298.1501)`` is ``a0`` and not ``a0`` plus about ``1e-8``.
    """
    if abs(temperature_k - REFERENCE_TEMPERATURE_K) < REFERENCE_TOLERANCE_K:
        return a[0]
    inverse = 1.0 / temperature_k
    inverse_reference = 1.0 / REFERENCE_TEMPERATURE_K
    return (
        a[0]
        + a[1] * (inverse - inverse_reference)
        + a[2] * math.log(temperature_k / REFERENCE_TEMPERATURE_K)
        + a[3] * (temperature_k - REFERENCE_TEMPERATURE_K)
        + a[4] * (temperature_k**2 - REFERENCE_TEMPERATURE_K**2)
        + a[5] * (inverse**2 - inverse_reference**2)
    )


def silvester_value(at_25: float, t1: float, t2: float, temperature_k: float) -> float:
    """Silvester and Pitzer's form, from `PitzerParameters.csv`.

    ``value_25 + t1 (1/T - 1/Tr) + t2 ln(T/Tr)``, and **zero `t1` and `t2` is the flat
    case**, where NeqSim returns ``value_25`` rather than evaluating the sum. The two
    differ only in the last bit, but they differ, so the branch is reproduced.
    """
    if abs(t1) < 1.0e-20 and abs(t2) < 1.0e-20:
        return at_25
    return (
        at_25
        + t1 * (1.0 / temperature_k - 1.0 / REFERENCE_TEMPERATURE_K)
        + t2 * math.log(temperature_k / REFERENCE_TEMPERATURE_K)
    )


def alpha1(first_charge: float, second_charge: float) -> float:
    """``PhasePitzer.getPitzerAlpha1``'s **bounded** form: 1.4 for a 2:2 pair, else 2.0.

    **Three call sites in NeqSim use an unbounded form instead** - ``getGamma``,
    ``getWaterGamma`` and ``phreeqcBinaryBprime`` all test ``|z| >= 1.5`` with no upper
    bound, so a 3-valent ion gets 1.4 there and 2.0 here. The two agree for every charge
    the vendored data carries (the highest is 2), so this is a latent inconsistency.
    :func:`alpha1_as_used` is what the activity coefficient is built from.
    """

    def bounded(z: float) -> bool:
        """NeqSim's own two-clause test, upper bound included."""
        return 1.5 <= abs(z) < 2.5

    return 1.4 if bounded(first_charge) and bounded(second_charge) else 2.0


def alpha1_as_used(first_charge: float, second_charge: float) -> float:
    """The ``alpha1`` the activity coefficient is actually built from."""
    if abs(first_charge) >= 1.5 and abs(second_charge) >= 1.5:
        return 1.4
    return 2.0


def alpha2(first_charge: float, second_charge: float) -> float:
    """``PhasePitzer.getPitzerAlpha2``: 12.0 when monovalent **or** 2:2, else 50.0.

    The 2:2 clause is not redundant with the monovalent one and is why ``CaCl2``'s
    qualified ``B2`` row survives: that row is a 2:1 term with ``alpha2 = 12``, and a
    2:2-only branch would have discarded it.
    """
    one, two = abs(first_charge), abs(second_charge)
    monovalent = one < 1.5 or two < 1.5
    two_two = 1.5 <= one < 2.5 and 1.5 <= two < 2.5
    return 12.0 if monovalent or two_two else 50.0


#: `b`, the Debye-Huckel denominator constant, the same `1.2` in every Pitzer expression.
B = 1.2

#: Below these an `x`, an ionic strength or a fitted `beta2` is treated as absent.
DERIVATIVE_FLOOR = 1.0e-12
BETA2_FLOOR = 1.0e-20

#: Below this a component's charge makes it a neutral rather than an ion.
ION_CHARGE = 0.5


def ln_gamma(
    molality: Sequence[float],
    charge: Sequence[float],
    temperature_k: float,
    a_phi: float,
    parameters: Any,
    component: int,
    *,
    common_ion_terms: bool,
    non_two_two_beta2: bool,
    unequal_charge_same_sign: bool,
) -> float:
    """``ln gamma_M`` for one ion: ``ComponentGePitzer.getGamma``'s ion branch.

    NeqSim's own comment gives the shape:

    .. code-block:: text

        ln(gamma_M) = z_M^2 F + sum_a m_a (2 B_Ma + Z C_Ma)
        F = -A_phi [ sqrt(I)/(1 + b sqrt(I)) + (2/b) ln(1 + b sqrt(I)) ]
              + sum_c sum_a m_c m_a B'_ca

    and the three pieces that are easy to lose are the ones the flags gate: ``fBprime`` is
    *not* accumulated when the common-ion terms are active, because
    :func:`common_ion_contribution` computes its own; ``beta2`` applies only for a 2:2 pair
    or when the dataset activated it elsewhere; and ``E_theta`` enters twice, once through
    ``F`` and once through the theta sum.

    ``parameters`` answers the six per-pair getters, so the caller supplies whichever
    dataset the selection rule chose.
    """
    z = charge[component]
    ionic_strength = electrolyte.ionic_strength(molality, charge)
    sqrt_i = math.sqrt(ionic_strength)

    # f^phi = -A_phi [ sqrt(I)/(1 + b sqrt(I)) + (2/b) ln(1 + b sqrt(I)) ]
    debye_huckel = -a_phi * (sqrt_i / (1.0 + B * sqrt_i) + (2.0 / B) * math.log(1.0 + B * sqrt_i))

    z_sum = sum(m * abs(x) for m, x in zip(molality, charge, strict=False) if abs(x) > ION_CHARGE)
    common_ion = (
        common_ion_contribution(molality, charge, temperature_k, parameters, z)
        if common_ion_terms
        else 0.0
    )

    m_this = molality[component]
    total = 0.0
    b_prime_sum = 0.0

    for j, other in enumerate(charge):
        if j == component or other * z >= 0.0:
            continue
        m_j = molality[j]
        beta0 = parameters.beta0(component, j, temperature_k)
        beta1 = parameters.beta1(component, j, temperature_k)
        cphi = parameters.cphi(component, j, temperature_k)

        alpha_one = alpha1_as_used(z, other)
        two_two = abs(z) >= 1.5 and abs(other) >= 1.5
        x1 = alpha_one * sqrt_i
        b_value = beta0 + beta1 * _g(x1)
        b_derivative = (
            beta1 * _g_prime(x1) / ionic_strength if ionic_strength > DERIVATIVE_FLOOR else 0.0
        )

        if two_two or non_two_two_beta2:
            beta2 = parameters.beta2(component, j, temperature_k)
            if abs(beta2) > BETA2_FLOOR:
                x2 = alpha2(z, other) * sqrt_i
                b_value += beta2 * _g(x2)
                if ionic_strength > DERIVATIVE_FLOOR:
                    b_derivative += beta2 * _g_prime(x2) / ionic_strength

        c_value = cphi / (2.0 * math.sqrt(abs(z * other)))
        total += m_j * (2.0 * b_value + z_sum * c_value)
        if not common_ion_terms:
            b_prime_sum += m_this * m_j * b_derivative

    e_theta_prime = 0.0
    if unequal_charge_same_sign:
        for first in range(len(molality)):
            if abs(charge[first]) < ION_CHARGE:
                continue
            for second in range(first + 1, len(molality)):
                other = charge[second]
                if charge[first] * other <= 0.0 or abs(charge[first] - other) < DERIVATIVE_FLOOR:
                    continue
                _, derivative = electrostatic.calculate(charge[first], other, ionic_strength, a_phi)
                e_theta_prime += molality[first] * molality[second] * derivative

    for j, other in enumerate(charge):
        if j == component or abs(other) < ION_CHARGE or other * z <= 0.0:
            continue
        m_j = molality[j]
        theta = parameters.theta(component, j, temperature_k)
        e_theta = 0.0
        if unequal_charge_same_sign and abs(z - other) >= DERIVATIVE_FLOOR:
            e_theta, _ = electrostatic.calculate(z, other, ionic_strength, a_phi)
        total += m_j * 2.0 * (theta + e_theta)

        for k, third in enumerate(charge):
            if third * z >= 0.0 or abs(third) < ION_CHARGE:
                continue
            total += m_j * molality[k] * parameters.psi(component, j, k, temperature_k)

    return z * z * (debye_huckel + b_prime_sum + e_theta_prime) + total + common_ion


def common_ion_contribution(
    molality: Sequence[float],
    charge: Sequence[float],
    temperature_k: float,
    parameters: Any,
    target_charge: float,
) -> float:
    """PHREEQC's common B-prime and C0 contributions to one ion's ``ln gamma``.

    ``ComponentGePitzer.phreeqcCommonIonContribution``: it sums over **every**
    cation-anion pair rather than over this ion's pairs, and adds a ``C0`` term the main
    loop has no counterpart for.
    """
    ionic_strength = electrolyte.ionic_strength(molality, charge)
    sqrt_i = math.sqrt(ionic_strength)
    b_prime = 0.0
    cphi_sum = 0.0
    for cation, cation_charge in enumerate(charge):
        if cation_charge <= 0.0:
            continue
        for anion, anion_charge in enumerate(charge):
            if anion_charge >= 0.0:
                continue
            first, second = molality[cation], molality[anion]
            b_prime += (
                first
                * second
                * _binary_b_prime(
                    charge, temperature_k, parameters, cation, anion, ionic_strength, sqrt_i
                )
            )
            cphi_sum += (
                first
                * second
                * parameters.cphi(cation, anion, temperature_k)
                / (2.0 * math.sqrt(abs(cation_charge * anion_charge)))
            )
    return target_charge**2 * b_prime + abs(target_charge) * cphi_sum


def _binary_b_prime(
    charge: Sequence[float],
    temperature_k: float,
    parameters: Any,
    first: int,
    second: int,
    ionic_strength: float,
    sqrt_i: float,
) -> float:
    """One pair's ``B'`` as the common-ion sum computes it, which applies ``beta2``
    wherever the dataset carries one rather than gating on the 2:2 test."""
    x1 = alpha1_as_used(charge[first], charge[second]) * sqrt_i
    derivative = 0.0
    if x1 > DERIVATIVE_FLOOR and ionic_strength > DERIVATIVE_FLOOR:
        derivative = parameters.beta1(first, second, temperature_k) * _g_prime(x1) / ionic_strength
    beta2 = parameters.beta2(first, second, temperature_k)
    x2 = alpha2(charge[first], charge[second]) * sqrt_i
    if abs(beta2) > BETA2_FLOOR and x2 > DERIVATIVE_FLOOR and ionic_strength > DERIVATIVE_FLOOR:
        derivative += beta2 * _g_prime(x2) / ionic_strength
    return derivative


#: The water molar mass NeqSim hard-codes in the osmotic expression, in kg/mol.
#:
#: ``getWaterGamma`` writes ``double Mw = 18.015`` in g/mol and divides by 1000, and the
#: databank's water row is ``0.018015 kg/mol`` - the two agree exactly.
WATER_MOLAR_MASS = 0.018015


def _electrostatic_phi(
    first_charge: float,
    second_charge: float,
    ionic_strength: float,
    a_phi: float,
    unequal_charge_same_sign: bool,
) -> float:
    """``E_theta + I dE_theta/dI``, the combination the *osmotic* route uses.

    The ion branch takes the value and the derivative separately - through the theta sum
    and through ``F`` - so this sum appears here and nowhere else.
    """
    if not unequal_charge_same_sign or abs(first_charge - second_charge) < DERIVATIVE_FLOOR:
        return 0.0
    value, derivative = electrostatic.calculate(first_charge, second_charge, ionic_strength, a_phi)
    return value + ionic_strength * derivative


def ln_gamma_water(
    molality: Sequence[float],
    charge: Sequence[float],
    temperature_k: float,
    a_phi: float,
    parameters: Any,
    solvent: int,
    x_water: float,
    neutral_osmotic: float,
    *,
    non_two_two_beta2: bool,
    unequal_charge_same_sign: bool,
    neutral_interactions_active: bool,
) -> float:
    """``ln gamma_w``, as ``ComponentGePitzer.getWaterGamma`` computes it.

    **This route does not go through the ion expression at all.** The solvent's activity
    coefficient comes from the Pitzer *osmotic* coefficient, and that is built from a
    different binary function:

    .. code-block:: text

        phi - 1 = (2 / sum m) [ -A_phi I^1.5/(1 + b sqrt(I))
                                + sum_c sum_a m_c m_a (B^phi_ca + Z C_ca) + theta + psi ]
        B^phi_ca = beta0 + beta1 exp(-alpha1 sqrt(I))
        ln a_w = -phi M_w sum m
        ln gamma_w = ln a_w - ln x_w

    where the ion branch's ``B`` is ``beta0 + beta1 g(alpha1 sqrt(I))``. **``exp(-alpha
    sqrt(I))`` and ``g(alpha sqrt(I))`` are different functions**, so a port that shared
    them would be wrong at every state.

    NeqSim returns the ideal ``gamma = 1`` when the water mole fraction falls below
    ``1e-10``, and that is reproduced here. :func:`_electrolyte.water_activity_coefficient`
    refuses the same state, and the two differ deliberately.
    """
    ionic_strength = electrolyte.ionic_strength(molality, charge)
    sqrt_i = math.sqrt(ionic_strength)
    n = len(molality)

    # Charged components, plus the neutral solutes when that layer is active - and
    # **`charge != 0` exactly**, a different test from the `|z| >= 0.5` the ion branch uses.
    sum_molalities = sum(
        molality[k]
        for k in range(n)
        if charge[k] != 0.0 or (neutral_interactions_active and k != solvent)
    )
    if sum_molalities < 1.0e-10:
        return 0.0

    # f^phi = -A_phi I^1.5 / (1 + b sqrt(I)); the base is non-negative by construction.
    f_phi = -a_phi * ionic_strength**1.5 / (1.0 + B * sqrt_i)

    # Z = sum |z_i| m_i over **every** component, water included; its charge is zero.
    z = sum(abs(x) * m for m, x in zip(molality, charge, strict=False))

    binary = 0.0
    for cation, cation_charge in enumerate(charge):
        if cation_charge <= 0.0:
            continue
        for anion, anion_charge in enumerate(charge):
            if anion_charge >= 0.0:
                continue
            beta0 = parameters.beta0(cation, anion, temperature_k)
            beta1 = parameters.beta1(cation, anion, temperature_k)
            cphi = parameters.cphi(cation, anion, temperature_k)
            alpha_one = alpha1_as_used(cation_charge, anion_charge)
            b_phi = beta0 + beta1 * math.exp(-alpha_one * sqrt_i)
            if (abs(cation_charge) >= 1.5 and abs(anion_charge) >= 1.5) or non_two_two_beta2:
                beta2 = parameters.beta2(cation, anion, temperature_k)
                if abs(beta2) > BETA2_FLOOR:
                    b_phi += beta2 * math.exp(-alpha2(cation_charge, anion_charge) * sqrt_i)
            product = abs(cation_charge * anion_charge)
            c_value = cphi / (2.0 * math.sqrt(product)) if product > 0.0 else 0.0
            binary += molality[cation] * molality[anion] * (b_phi + z * c_value)

    theta_psi = 0.0
    for first, first_charge in enumerate(charge):
        for second in range(first + 1, n):
            second_charge = charge[second]
            same_sign = (first_charge > 0.0 and second_charge > 0.0) or (
                first_charge < 0.0 and second_charge < 0.0
            )
            if not same_sign:
                continue
            m_first, m_second = molality[first], molality[second]
            theta = parameters.theta(first, second, temperature_k)
            theta_psi += (
                m_first
                * m_second
                * (
                    theta
                    + _electrostatic_phi(
                        first_charge,
                        second_charge,
                        ionic_strength,
                        a_phi,
                        unequal_charge_same_sign,
                    )
                )
            )
            # The third ion is the opposite sign to the pair.
            for third, third_charge in enumerate(charge):
                if third_charge * first_charge >= 0.0:
                    continue
                theta_psi += (
                    m_first
                    * m_second
                    * molality[third]
                    * parameters.psi(first, second, third, temperature_k)
                )

    phi = 1.0 + (2.0 / sum_molalities) * (f_phi + binary + theta_psi + neutral_osmotic)
    ln_a_w = electrolyte.ln_water_activity(phi, sum_molalities, WATER_MOLAR_MASS)
    if x_water < 1.0e-10:
        return 0.0
    return ln_a_w - math.log(x_water)


#: The neutral-solute families: a neutral with a neutral, an ion, or a pair of ions.
LAMBDA, ZETA, MU, ETA = "LAMBDA", "ZETA", "MU", "ETA"


class NeutralInteraction:
    """One neutral-family interaction over a tuple of components.

    **The coefficients come from the tuple's repetition structure**, not from the file:

    =========================  =========================  =====================
    tuple                      ``log_gamma_coefficients``  ``osmotic_coefficient``
    =========================  =========================  =====================
    two equal indexes          ``[1, 1]``                 ``0.5``
    two different              ``[2, 2]``                 ``1.0``
    three, family ``MU``       ``[m, m, m]``              ``m`` for ``m`` in {1, 3, 6}
    three, otherwise           ``[1, 1, 1]``              ``1.0``
    =========================  =========================  =====================

    with ``m`` counting the tuple's distinct permutations. **The repeated-species case is
    the one to get right** - PHREEQC differentiates both slots before accumulating them,
    which is why it is ``[1, 1]`` and a half rather than ``[2, 2]`` and one.
    """

    __slots__ = ("family", "form", "indexes", "log_gamma_coefficients", "osmotic_coefficient")

    def __init__(self, family: str, indexes: Sequence[int], form: Sequence[float]) -> None:
        self.indexes = sorted(indexes)
        self.family = family
        self.form = form
        if len(self.indexes) == 2:
            if self.indexes[0] == self.indexes[1]:
                self.log_gamma_coefficients, self.osmotic_coefficient = [1.0, 1.0], 0.5
            else:
                self.log_gamma_coefficients, self.osmotic_coefficient = [2.0, 2.0], 1.0
        elif len(self.indexes) == 3 and family == MU:
            if self.indexes[0] == self.indexes[2]:
                multiplicity = 1.0
            elif self.indexes[0] == self.indexes[1] or self.indexes[1] == self.indexes[2]:
                multiplicity = 3.0
            else:
                multiplicity = 6.0
            self.log_gamma_coefficients = [multiplicity] * 3
            self.osmotic_coefficient = multiplicity
        else:
            self.log_gamma_coefficients = [1.0] * len(self.indexes)
            self.osmotic_coefficient = 1.0

    def log_gamma_contribution(
        self, molality: Sequence[float], component: int, temperature: float
    ) -> float:
        """This interaction's contribution to one component's ``ln gamma``.

        The component contributes once per *position* it occupies, and the molality product
        is over the tuple's other members - so a component appearing twice contributes
        twice from one tuple.
        """
        parameter = catalog_value(self.form, temperature)
        contribution = 0.0
        for position, index in enumerate(self.indexes):
            if index != component:
                continue
            product = 1.0
            for other, other_index in enumerate(self.indexes):
                if other != position:
                    product *= molality[other_index]
            contribution += self.log_gamma_coefficients[position] * product * parameter
        return contribution

    def osmotic_contribution(self, molality: Sequence[float], temperature: float) -> float:
        """This interaction's contribution to PHREEQC's osmotic sum."""
        product = 1.0
        for index in self.indexes:
            product *= molality[index]
        return self.osmotic_coefficient * product * catalog_value(self.form, temperature)


def ln_gamma_neutral(
    interactions: Sequence[NeutralInteraction],
    molality: Sequence[float],
    component: int,
    temperature: float,
) -> float:
    """Every neutral-family contribution to one component's ``ln gamma``."""
    return sum(i.log_gamma_contribution(molality, component, temperature) for i in interactions)


def osmotic_neutral(
    interactions: Sequence[NeutralInteraction],
    molality: Sequence[float],
    temperature: float,
) -> float:
    """Every neutral-family contribution to the osmotic sum."""
    return sum(i.osmotic_contribution(molality, temperature) for i in interactions)


def catalogue_interactions(
    species: Sequence[str],
    charge: Sequence[float],
    ions: Sequence[int],
    neutrals: Sequence[int],
) -> list[NeutralInteraction] | None:
    """The neutral interactions the PHREEQC catalogue imposes on a topology.

    ``PitzerParameterDatasets.applyCatalogNeutralRows``: for each neutral, a ``LAMBDA`` with
    every neutral **including itself** and with every ion, and a ``ZETA`` with every
    cation-anion pair. **The sign filter on that pair is not optional** - without it the
    lookup asks for ``(CO2, Na+, Na+)``, which the catalogue does not carry, and the whole
    layer is abandoned.

    ``None`` for a tuple the catalogue does not carry, which is the coverage rule's answer
    and not a defect: the whole dataset is abandoned rather than half-applied.
    """
    out: list[NeutralInteraction] = []
    for position, neutral in enumerate(neutrals):
        for second in neutrals[position:]:
            form = find(LAMBDA, [species[neutral], species[second]])
            if form is None:
                return None
            out.append(NeutralInteraction(LAMBDA, [neutral, second], form))
        for ion in ions:
            form = find(LAMBDA, [species[neutral], species[ion]])
            if form is None:
                return None
            out.append(NeutralInteraction(LAMBDA, [neutral, ion], form))
        for cation in ions:
            if charge[cation] <= 0.0:
                continue
            for anion in ions:
                if charge[anion] >= 0.0:
                    continue
                form = find(ZETA, [species[neutral], species[cation], species[anion]])
                if form is None:
                    return None
                out.append(NeutralInteraction(ZETA, [neutral, cation, anion], form))
    return out
