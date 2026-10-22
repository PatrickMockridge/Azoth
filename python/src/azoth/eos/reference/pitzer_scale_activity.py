"""``eos.pitzer_scale_activity`` - NeqSim's ``PitzerScaleActivityModel``.

A binary Pitzer model over the PHREEQC ``pitzer.dat`` parameter sets: mean activity
coefficients from the six-term temperature expression, and trace ions mapped at an
equivalent ionic strength by the MacInnes convention.

Spec: ``specs/calcs/eos/pitzer_scale_activity.toml``.
"""

from __future__ import annotations

import math
from typing import Any

from azoth.core.errors import InvalidInputError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import PitzerScaleActivityResult
from azoth.core.units import Q, input_to_si
from azoth.core.warnings import Warning
from azoth.eos.reference._pitzer_catalog import find

CALC_ID = "eos.pitzer_scale_activity"

PITZER_B = 1.2
MIN_IONIC_STRENGTH = 1.0e-12
REFERENCE_TEMPERATURE_K = 298.15

IONS = (
    "sodium",
    "chloride",
    "calcium",
    "barium",
    "strontium",
    "ferrous",
    "sulphate",
    "carbonate",
    "bicarbonate",
)

#: cation, anion, cation stoichiometry, anion stoichiometry, charges, alpha1, alpha2.
SALTS: dict[str, tuple[str, str, float, float, float, float, float, float]] = {
    "nacl": ("Na+", "Cl-", 1.0, 1.0, 1.0, 1.0, 2.0, 0.0),
    "cacl2": ("Ca++", "Cl-", 1.0, 2.0, 2.0, 1.0, 2.0, 12.0),
    "bacl2": ("Ba++", "Cl-", 1.0, 2.0, 2.0, 1.0, 2.0, 0.0),
    "srcl2": ("Sr++", "Cl-", 1.0, 2.0, 2.0, 1.0, 2.0, 0.0),
    "fecl2": ("Fe++", "Cl-", 1.0, 2.0, 2.0, 1.0, 2.0, 0.0),
    "na2so4": ("Na+", "SO4--", 2.0, 1.0, 1.0, 2.0, 2.0, 0.0),
    "na2co3": ("Na+", "CO3--", 2.0, 1.0, 1.0, 2.0, 2.0, 0.0),
    "nahco3": ("Na+", "HCO3-", 1.0, 1.0, 1.0, 1.0, 2.0, 12.0),
}


def _family(
    salt: tuple[str, str, float, float, float, float, float, float], name: str
) -> tuple[float, ...]:
    values = find(name, (salt[0], salt[1]))
    return values if values is not None else (0.0,) * 6


def _temperature_parameter(coefficients: tuple[float, ...], temperature_k: float) -> float:
    tr = REFERENCE_TEMPERATURE_K
    return (
        coefficients[0]
        + coefficients[1] * (1.0 / temperature_k - 1.0 / tr)
        + coefficients[2] * math.log(temperature_k / tr)
        + coefficients[3] * (temperature_k - tr)
        + coefficients[4] * (temperature_k * temperature_k - tr * tr)
        + coefficients[5] * (1.0 / (temperature_k * temperature_k) - 1.0 / (tr * tr))
    )


def _beta_gamma_term(beta: float, alpha: float, ionic_strength: float) -> float:
    if alpha <= 0.0 or beta == 0.0:
        return 0.0
    strength = max(ionic_strength, MIN_IONIC_STRENGTH)
    x = alpha * math.sqrt(strength)
    numerator = 1.0 - (1.0 + x - 0.5 * x * x) * math.exp(-x)
    return 2.0 * beta * numerator / (alpha * alpha * strength)


def _mean(
    salt: tuple[str, str, float, float, float, float, float, float],
    molality: float,
    temperature_k: float,
    debye_huckel_a: float,
) -> float:
    if molality == 0.0:
        return 1.0
    _, _, ns, na, zs, za, alpha1, alpha2 = salt
    ionic_strength = 0.5 * molality * (ns * zs * zs + na * za * za)
    sqrt_i = math.sqrt(max(ionic_strength, MIN_IONIC_STRENGTH))
    a_phi = debye_huckel_a * math.log(10.0) / 3.0
    f_gamma = -a_phi * (
        sqrt_i / (1.0 + PITZER_B * sqrt_i) + 2.0 / PITZER_B * math.log(1.0 + PITZER_B * sqrt_i)
    )
    beta0 = _temperature_parameter(_family(salt, "B0"), temperature_k)
    beta1 = _temperature_parameter(_family(salt, "B1"), temperature_k)
    beta2 = _temperature_parameter(_family(salt, "B2"), temperature_k)
    b_gamma = 2.0 * beta0 + _beta_gamma_term(beta1, alpha1, ionic_strength)
    if alpha2 > 0.0 and abs(beta2) > 0.0:
        b_gamma += _beta_gamma_term(beta2, alpha2, ionic_strength)
    nu = ns + na
    b_factor = 2.0 * ns * na / nu
    c_factor = 3.0 * math.pow(ns * na, 1.5) / nu
    c_phi = _temperature_parameter(_family(salt, "C0"), temperature_k)
    ln_gamma = (
        abs(zs * za) * f_gamma
        + molality * b_factor * b_gamma
        + molality * molality * c_factor * c_phi
    )
    return math.exp(ln_gamma)


def _mean_at_equivalent_strength(
    name: str, ionic_strength: float, temperature_k: float, debye_huckel_a: float
) -> float:
    salt = SALTS[name]
    _, _, ns, na, zs, za, _, _ = salt
    per_molality = 0.5 * (ns * zs * zs + na * za * za)
    return _mean(salt, ionic_strength / per_molality, temperature_k, debye_huckel_a)


def pitzer_scale_activity(
    temperature: Q,
    ionic_strength_molal: Q,
    debye_huckel_a: float,
    salt: str,
    salt_molality: Q,
    ion: str,
    cation: str,
    anion: str,
) -> PitzerScaleActivityResult:
    """The binary Pitzer activity coefficients of a scale-forming brine.

    Args:
        temperature: absolute temperature; the expressions are about 298.15 K.
        ionic_strength_molal: the strength that configures the trace-ion mapping.
        debye_huckel_a: the Debye-Huckel A in log10 units.
        salt: the binary salt whose mean coefficient is wanted.
        salt_molality: the salt's molality.
        ion: the ion whose mapped coefficient is reported.
        cation: the product's cation.
        anion: the product's anion.

    Raises:
        InvalidInputError: for a salt or ion name outside the model's set.
        OutOfRangeError: if a temperature or strength is outside the spec's range.

    Example:
        >>> import azoth
        >>> q = azoth.ureg.Quantity
        >>> r = pitzer_scale_activity(
        ...     q(298.15, "K"), q(1.0, "mol/kg"), 0.5107140231833494, "nacl",
        ...     q(1.0, "mol/kg"), "sulphate", "barium", "sulphate")
        >>> abs(r.mean_activity_coefficient - 0.6566052005348894) < 1e-12
        True
    """
    from azoth._registry_gen import spec as _spec_for

    spec: dict[str, Any] = _spec_for(CALC_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    temperature_k = input_to_si(spec, "temperature", temperature)
    strength_si = input_to_si(spec, "ionic_strength_molal", ionic_strength_molal)
    molality_si = input_to_si(spec, "salt_molality", salt_molality)
    apply_checks(
        checks.on_input,
        {
            "temperature": temperature_k,
            "ionic_strength_molal": strength_si,
            "debye_huckel_a": debye_huckel_a,
            "salt_molality": molality_si,
        }.get,
        warnings,
    )

    if salt not in SALTS:
        raise InvalidInputError("salt", f"`{salt}` is not one of the model's eight binary salts")
    for name, role in ((ion, "ion"), (cation, "cation"), (anion, "anion")):
        if name not in IONS:
            raise InvalidInputError(role, f"`{name}` is not one of the model's nine ions: {IONS}")

    selected = SALTS[salt]
    mean = _mean(selected, molality_si, temperature_k, debye_huckel_a)

    strength = max(strength_si, MIN_IONIC_STRENGTH)
    gamma_nacl = _mean_at_equivalent_strength("nacl", strength, temperature_k, debye_huckel_a)

    def divalent(name: str) -> float:
        value = _mean_at_equivalent_strength(name, strength, temperature_k, debye_huckel_a)
        return math.pow(value, 3.0) / (gamma_nacl * gamma_nacl)

    gamma_nahco3 = _mean_at_equivalent_strength("nahco3", strength, temperature_k, debye_huckel_a)
    gammas = (
        gamma_nacl,
        gamma_nacl,
        divalent("cacl2"),
        divalent("bacl2"),
        divalent("srcl2"),
        divalent("fecl2"),
        divalent("na2so4"),
        divalent("na2co3"),
        gamma_nahco3 * gamma_nahco3 / gamma_nacl,
    )

    def gamma_of(name: str) -> float:
        return gammas[IONS.index(name)]

    return PitzerScaleActivityResult(
        mean_activity_coefficient=mean,
        ion_activity_coefficient=gamma_of(ion),
        activity_coefficient_product=gamma_of(cation) * gamma_of(anion),
        gamma_sodium=gammas[0],
        gamma_chloride=gammas[1],
        gamma_calcium=gammas[2],
        gamma_barium=gammas[3],
        gamma_strontium=gammas[4],
        gamma_ferrous=gammas[5],
        gamma_sulphate=gammas[6],
        gamma_carbonate=gammas[7],
        gamma_bicarbonate=gammas[8],
        warnings=tuple(warnings),
    )
