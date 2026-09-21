"""``eos.scale_saturation_ratio`` - one salt's saturation ratio in a brine.

```text
m_i = x_i / (x_water M_water)                    [mol/kg water]
Ksp = exp(A/T + B + C ln T + D T + E/T**2) exp(-dV (P - P0)/(R T))
SR  = (gamma1 m1)**stoc1 (gamma2 m2)**stoc2 a_water**waterstoc / Ksp
```

Spec: ``specs/calcs/eos/scale_saturation_ratio.toml``, which carries the three name
overrides, the clamps and the measurement that pins both worked cases.

NeqSim's ``CheckScalePotential``, one row of its ``compsalt`` walk. **`Vdelta` is in cm3/mol
and the pressure term's `R` is `83.1446`** - the term is dimensionless as written, and the SI
constant against a cm3 volume is `1e4` out.
"""

from __future__ import annotations

from azoth._registry_gen import spec as _spec_for
from azoth.core.errors import InvalidInputError, OutOfRangeError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import ScaleSaturationRatioResult
from azoth.core.units import Q, input_to_si
from azoth.core.warnings import Warning

CALC_ID = "eos.scale_saturation_ratio"

#: The gas constant in the units `Vdelta` and the pressure term are stated in: cm3 bar/(mol K).
R_CM3_BAR = 83.1446

#: The pressure the volume term is referred to, in bar - one atmosphere.
REFERENCE_PRESSURE_BAR = 1.01325

#: Below this molality the ion activity product is taken as zero.
MIN_MOLALITY = 1.0e-30

#: The clamp on `ln SR`, NeqSim's own: `exp(69)` is about `1e30`.
LN_SR_CLAMP = 69.0


def scale_saturation_ratio(
    salt: str,
    x1: float,
    x2: float,
    x_water: float,
    gamma1: float,
    gamma2: float,
    water_activity: float,
    T: Q,
    P: Q,
    h3o_molality: Q | None = None,
) -> ScaleSaturationRatioResult:
    """One salt's saturation ratio in an aqueous phase.

    Args:
        salt: the name as ``compsalt`` spells it, which selects the row and the three
            overrides.
        x1: the first ion's mole fraction in the aqueous phase.
        x2: the second ion's mole fraction.
        x_water: water's mole fraction in the same phase.
        gamma1: the first ion's activity coefficient, from the phase's own model.
        gamma2: the second ion's.
        water_activity: water's activity, read only where the row states a non-zero
            ``waterstoc``.
        T: absolute temperature.
        P: absolute pressure.
        h3o_molality: the hydrogen ion's molality, **read for ``FeS`` alone**.

    Raises:
        InvalidInputError: if the salt is not in the table, or ``FeS`` is asked for without a
            hydrogen molality.
        OutOfRangeError: if ``T`` or ``P`` is not positive, or the solubility product has no
            value.
    """
    from azoth.eos import components as databank

    spec = _spec_for(CALC_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []
    t_si = input_to_si(spec, "T", T)
    p_si = input_to_si(spec, "P", P)
    apply_checks(checks.on_input, {"T": t_si, "P": p_si}.get, warnings)

    record = databank.salt(salt)
    if record is None:
        raise InvalidInputError(
            "salt", f"`{salt}` is not a row of `compsalt`, so it has no solubility product"
        )
    t = t_si
    p_bar = p_si / 1.0e5

    # The row's own correlation, which three names then override.
    import math

    ksp = math.exp(
        record.ksp[0] / t
        + record.ksp[1]
        + record.ksp[2] * math.log(t)
        + record.ksp[3] * t
        + record.ksp[4] / (t * t)
    )
    if salt == "NaCl":
        ksp = 92.78 - 0.407 * t + 0.000747 * t * t
    elif salt == "CaCO3":
        # Plummer & Busenberg (1982), for calcite.
        log10_ksp = -171.9065 - 0.077993 * t + 2839.319 / t + 71.595 * math.log10(t)
        ksp = 10.0**log10_ksp
    elif salt == "FeCO3":
        # Greenberg & Tomson (1992).
        log10_ksp = -59.3498 - 0.041377 * t - 2.1963 / t + 24.5724 * math.log10(t)
        ksp = 10.0**log10_ksp
    elif salt == "FeS":
        # **The one salt whose product carries an ion molality**, and NeqSim skips the salt
        # outright where the phase has no `H3O+` - which a scalar calc cannot do, so the
        # caller says which it is by stating the molality or by not calling.
        if h3o_molality is None:
            raise InvalidInputError(
                "h3o_molality",
                "`FeS`'s solubility product is multiplied by the hydrogen-ion molality, and "
                "NeqSim skips the salt altogether where the phase has no `H3O+`; this calc "
                "has no phase to look in, so the caller states it",
            )
        ksp *= input_to_si(spec, "h3o_molality", h3o_molality)

    # The pressure term, guarded exactly as NeqSim guards it.
    if abs(record.volume_delta) > 1.0e-10 and p_bar > 1.013:
        correction = -record.volume_delta * (p_bar - REFERENCE_PRESSURE_BAR) / (R_CM3_BAR * t)
        ksp *= math.exp(min(max(correction, -LN_SR_CLAMP), LN_SR_CLAMP))
    if not ksp > 0.0 or not math.isfinite(ksp):
        raise OutOfRangeError(
            "Ksp", ksp, f"the solubility product of {salt} has no value at {t} K and {p_bar} bar"
        )

    water = databank.entry("water")
    if water.molar_mass is None:
        raise InvalidInputError(
            "components", "water carries no molar mass, and a molality is per kilogram of it"
        )
    molar_mass_water = float(water.molar_mass.to("kg/mol").magnitude)
    m1 = x1 / (x_water * molar_mass_water)
    m2 = x2 / (x_water * molar_mass_water)

    # **A ratio too far below saturation to matter is a zero**, and the guard is before the
    # logarithm rather than after it.
    saturation_ratio = 0.0
    iap = 0.0
    if m1 >= MIN_MOLALITY and m2 >= MIN_MOLALITY and gamma1 > 0.0 and gamma2 > 0.0:
        ln_iap = record.cation_stoichiometry * math.log(
            gamma1 * m1
        ) + record.anion_stoichiometry * math.log(gamma2 * m2)
        iap = (gamma1 * m1) ** record.cation_stoichiometry * (
            gamma2 * m2
        ) ** record.anion_stoichiometry
        if record.water_stoichiometry > 0.0:
            if not water_activity > 0.0 or not math.isfinite(water_activity):
                return ScaleSaturationRatioResult(
                    saturation_ratio=0.0,
                    ion_activity_product=0.0,
                    solubility_product=ksp,
                    warnings=tuple(warnings),
                )
            ln_iap += record.water_stoichiometry * math.log(water_activity)
            iap *= water_activity**record.water_stoichiometry
        ln_sr = ln_iap - math.log(ksp)
        if ln_sr < -LN_SR_CLAMP:
            saturation_ratio = 0.0
        elif ln_sr > LN_SR_CLAMP:
            saturation_ratio = math.exp(LN_SR_CLAMP)
        else:
            saturation_ratio = math.exp(ln_sr)

    return ScaleSaturationRatioResult(
        saturation_ratio=saturation_ratio,
        ion_activity_product=iap,
        solubility_product=ksp,
        warnings=tuple(warnings),
    )


__all__ = ["scale_saturation_ratio"]
