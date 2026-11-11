"""``eos.production_chemical_scale_scenario`` - NeqSim's ``ProductionChemicalScaleScenario``.

A produced water before and after chemical dosing: pH adjusters as acid/base equivalents on
a closed carbonate balance, H2S scavengers by capacity, scale inhibitors as kinetic controls,
and the coupled mineral equilibrium at both states.

Spec: ``specs/models/eos/production_chemical_scale_scenario.toml``.
"""

from __future__ import annotations

import math
from typing import Any, cast

from azoth.core.range import apply_checks, checks_for
from azoth.core.record_outputs_gen import (
    ChemicalWarningRow,
    ScenarioChangeRow,
    ScenarioIndexRow,
)
from azoth.core.result import ProductionChemicalScaleScenarioResult
from azoth.core.units import Q, from_si, input_records_to_si, input_to_si
from azoth.core.warnings import Warning
from azoth.eos.reference.multi_mineral_scale_equilibrium import (
    multi_mineral_scale_equilibrium,
)

CALC_ID = "eos.production_chemical_scale_scenario"

#: The class's constants: bicarbonate and H2S molar masses, and the 25 C water ion product.
MW_HCO3 = 61.016
MW_H2S = 34.08
WATER_KW_25C = 1.0e-14


def _key(active: str) -> str:
    return active.lower().replace("-", "").replace(" ", "")


def _base_equivalents(active: str, active_mg_l: float) -> float | None:
    key = _key(active)
    grams_per_litre = active_mg_l / 1000.0
    if "naoh" in key or "sodiumhydroxide" in key or "caustic" in key:
        return grams_per_litre / 39.997
    if "na2co3" in key or "sodiumcarbonate" in key or "sodaash" in key:
        return 2.0 * grams_per_litre / 105.9888
    if "mdea" in key or "methyldiethanolamine" in key:
        return grams_per_litre / 119.163
    if key == "mea" or "monoethanolamine" in key:
        return grams_per_litre / 61.084
    if "hcl" in key or "hydrochloricacid" in key:
        return -grams_per_litre / 36.461
    return None


def _scavenger_capacity_mg_l(active: str, active_mg_l: float) -> float | None:
    key = _key(active)
    if "meatriazine" in key:
        molar_mass, per_mol = 189.26, 1.0
    elif "mmatriazine" in key:
        molar_mass, per_mol = 129.21, 1.0
    elif "iron" in key or key == "fe":
        molar_mass, per_mol = 55.845, 1.0
    elif "glyoxal" in key:
        molar_mass, per_mol = 58.04, 1.0 / 1.5
    else:
        return None
    return active_mg_l / molar_mass * per_mol * MW_H2S


def _carbonate_k1(temperature_k: float) -> float:
    log_k1 = (
        -356.3094
        - 0.06091964 * temperature_k
        + 21834.37 / temperature_k
        + 126.8339 * math.log10(temperature_k)
        - 1684915.0 / (temperature_k * temperature_k)
    )
    return math.pow(10.0, log_k1)


def _carbonate_k2(temperature_k: float) -> float:
    log_k2 = (
        -107.8871
        - 0.03252849 * temperature_k
        + 5151.79 / temperature_k
        + 38.92561 * math.log10(temperature_k)
        - 563713.9 / (temperature_k * temperature_k)
    )
    return math.pow(10.0, log_k2)


def _carbonate_fractions(ph: float, k1: float, k2: float) -> tuple[float, float, float]:
    h = math.pow(10.0, -ph)
    denominator = h * h + k1 * h + k1 * k2
    return h * h / denominator, k1 * h / denominator, k1 * k2 / denominator


def _alkalinity(total_carbon: float, ph: float, k1: float, k2: float) -> float:
    alpha = _carbonate_fractions(ph, k1, k2)
    h = math.pow(10.0, -ph)
    oh = WATER_KW_25C / h
    return total_carbon * (alpha[1] + 2.0 * alpha[2]) + oh - h


def _equilibrium(
    ions_mg_l: dict[str, float],
    tds_mg_l: float,
    ph: float,
    temperature_k: float,
    pressure_pa: float,
    activity_model: str,
    pitzer_ionic_strength_molal: float | None,
) -> Any:
    def mass(name: str) -> Q:
        return from_si(ions_mg_l[name] * 1e-3, "kg/m**3")

    return multi_mineral_scale_equilibrium(
        mass("calcium"),
        mass("barium"),
        mass("strontium"),
        mass("iron"),
        mass("magnesium"),
        mass("sodium"),
        mass("bicarbonate"),
        mass("sulphate"),
        from_si(tds_mg_l * 1e-3, "kg/m**3"),
        from_si(temperature_k, "K"),
        from_si(pressure_pa, "Pa"),
        from_si(0.0, "Pa"),
        ph,
        False,
        False,
        activity_model,
        None,
        None,
        pitzer_ionic_strength_molal,
    )


def production_chemical_scale_scenario(
    calcium: Q,
    barium: Q,
    strontium: Q,
    iron: Q,
    magnesium: Q,
    sodium: Q,
    sulphate: Q,
    bicarbonate: Q,
    total_dissolved_solids: Q,
    ph: float,
    dissolved_h2s: Q,
    temperature: Q,
    pressure: Q,
    activity_model: str,
    chemicals: list[dict[str, Any]],
    pitzer_ionic_strength_molal: float | None = None,
) -> ProductionChemicalScaleScenarioResult:
    """The treatment-aware scale scenario of a produced water.

    Args:
        calcium: calcium ion concentration.
        barium: barium ion concentration.
        strontium: strontium ion concentration.
        iron: ferrous iron concentration.
        magnesium: magnesium ion concentration.
        sodium: sodium ion concentration.
        sulphate: sulphate ion concentration.
        bicarbonate: bicarbonate concentration as HCO3-.
        total_dissolved_solids: total dissolved solids.
        ph: the untreated aqueous pH.
        dissolved_h2s: dissolved H2S as H2S.
        temperature: absolute temperature.
        pressure: absolute pressure.
        activity_model: ``davies``, ``bdot`` or ``pitzer_binary``.
        chemicals: the dosed chemicals, as the class's ``ProductionChemical`` records.
        pitzer_ionic_strength_molal: the Pitzer strength; absent uses the estimate.

    Returns:
        The treated state, the two index sets, their change and the class's warnings.

    Raises:
        OutOfRangeError: if an input is outside the spec's hard range.

    Example:
        >>> import azoth
        >>> q = azoth.ureg.Quantity
        >>> row = {"chemical_type": "ph_adjuster", "name": "caustic",
        ...        "active_ingredient": "NaOH", "active_wt_pct": 50.0,
        ...        "dosage": q(0.1, "kg/m**3")}
        >>> r = production_chemical_scale_scenario(
        ...     q(2.0, "kg/m**3"), q(0.05, "kg/m**3"), q(0.1, "kg/m**3"),
        ...     q(0.005, "kg/m**3"), q(0.5, "kg/m**3"), q(20.0, "kg/m**3"),
        ...     q(0.5, "kg/m**3"), q(0.8, "kg/m**3"), q(80.0, "kg/m**3"), 6.5,
        ...     q(0.02, "kg/m**3"), q(333.15, "K"), q(1.0e7, "Pa"), "davies",
        ...     [row], None)
        >>> abs(r.treated_ph - 6.612286718254484) < 1e-12
        True
    """
    from azoth._models_gen import model as _model_for

    spec: dict[str, Any] = _model_for(CALC_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    values: dict[str, float] = {
        name: input_to_si(spec, name, value)
        for name, value in {
            "calcium": calcium,
            "barium": barium,
            "strontium": strontium,
            "iron": iron,
            "magnesium": magnesium,
            "sodium": sodium,
            "sulphate": sulphate,
            "bicarbonate": bicarbonate,
            "total_dissolved_solids": total_dissolved_solids,
            "dissolved_h2s": dissolved_h2s,
            "temperature": temperature,
            "pressure": pressure,
        }.items()
    }
    values["ph"] = ph
    apply_checks(checks.on_input, values.get, warnings)

    rows = cast("list[dict[str, Any]]", input_records_to_si(spec, "chemicals", chemicals))

    ions_mg_l = {
        "calcium": values["calcium"] * 1000.0,
        "barium": values["barium"] * 1000.0,
        "strontium": values["strontium"] * 1000.0,
        "iron": values["iron"] * 1000.0,
        "magnesium": values["magnesium"] * 1000.0,
        "sodium": values["sodium"] * 1000.0,
        "sulphate": values["sulphate"] * 1000.0,
        "bicarbonate": values["bicarbonate"] * 1000.0,
    }
    tds_mg_l = values["total_dissolved_solids"] * 1000.0
    dissolved_h2s_mg_l = values["dissolved_h2s"] * 1000.0

    equivalents_mol_l = 0.0
    capacity_mg_l = 0.0
    scale_inhibitor = False
    model_warnings: list[ChemicalWarningRow] = []

    for row in rows:
        active_mg_l = max(0.0, row["dosage"]) * 1000.0 * max(0.0, row["active_wt_pct"]) / 100.0
        chemical_type = str(row["chemical_type"])
        if chemical_type in ("acid", "ph_adjuster"):
            capacity = row.get("alkalinity_capacity")
            equivalents = (
                capacity * active_mg_l * 1.0e-6
                if capacity is not None
                else _base_equivalents(str(row["active_ingredient"]), active_mg_l)
            )
            if equivalents is None:
                model_warnings.append(
                    ChemicalWarningRow(
                        code=f"unsupported_ph_adjuster_{row['name']}",
                        message=(
                            "No acid/base-equivalent model for active ingredient "
                            f"{row['active_ingredient']}"
                        ),
                    )
                )
            else:
                equivalents_mol_l += equivalents
        elif chemical_type == "h2s_scavenger":
            capacity = row.get("h2s_capacity")
            scavenger = (
                capacity * active_mg_l
                if capacity is not None
                else _scavenger_capacity_mg_l(str(row["active_ingredient"]), active_mg_l)
            )
            if scavenger is None:
                model_warnings.append(
                    ChemicalWarningRow(
                        code=f"unsupported_h2s_scavenger_{row['name']}",
                        message=(
                            "No H2S stoichiometric model for active ingredient "
                            f"{row['active_ingredient']}"
                        ),
                    )
                )
            else:
                capacity_mg_l += scavenger
                if "triazine" in _key(str(row["active_ingredient"])):
                    model_warnings.append(
                        ChemicalWarningRow(
                            code="triazine_spent_product",
                            message=(
                                "Check dithiazine/amorphous spent-product deposition; "
                                "scavenger removal is not mineral equilibrium"
                            ),
                        )
                    )
        elif chemical_type == "scale_inhibitor":
            scale_inhibitor = True
        else:
            model_warnings.append(
                ChemicalWarningRow(
                    code=f"no_equilibrium_effect_{row['name']}",
                    message=(
                        "Chemical retained in audit trail; no supported stoichiometric "
                        "effect on this aqueous scale scenario"
                    ),
                )
            )

    treated_ph, treated_bicarbonate_mol_l = _apply_base_dose(
        ph, ions_mg_l["bicarbonate"], equivalents_mol_l, values["temperature"], model_warnings
    )
    treated_bicarbonate_mg_l = treated_bicarbonate_mol_l * MW_HCO3 * 1000.0
    residual_h2s_mg_l = max(0.0, dissolved_h2s_mg_l - capacity_mg_l)

    baseline = _equilibrium(
        ions_mg_l,
        tds_mg_l,
        ph,
        values["temperature"],
        values["pressure"],
        activity_model,
        pitzer_ionic_strength_molal,
    )
    treated_ions = dict(ions_mg_l)
    treated_ions["bicarbonate"] = treated_bicarbonate_mg_l
    treated = _equilibrium(
        treated_ions,
        tds_mg_l,
        treated_ph,
        values["temperature"],
        values["pressure"],
        activity_model,
        pitzer_ionic_strength_molal,
    )

    baseline_si = tuple(
        ScenarioIndexRow(mineral=m.name, saturation_index=m.initial_saturation_index)
        for m in baseline.minerals
    )
    treated_si = tuple(
        ScenarioIndexRow(mineral=m.name, saturation_index=m.initial_saturation_index)
        for m in treated.minerals
    )
    treated_by_name = {row.mineral: row.saturation_index for row in treated_si}
    changes: list[ScenarioChangeRow] = []
    for index_row in baseline_si:
        after = treated_by_name.get(index_row.mineral)
        change = (
            None
            if index_row.saturation_index is None or after is None
            else after - index_row.saturation_index
        )
        changes.append(ScenarioChangeRow(mineral=index_row.mineral, change=change))
    si_change = tuple(changes)

    delta_calcite = next((row.change for row in si_change if row.mineral == "CaCO3"), None)
    if scale_inhibitor:
        model_warnings.append(
            ChemicalWarningRow(
                code="inhibitor_is_kinetic",
                message=(
                    "Scale inhibitor does not change thermodynamic SI; evaluate MIC and "
                    "growth control separately with ScaleControlAssessor"
                ),
            )
        )
    if delta_calcite is not None and delta_calcite > 0.25:
        model_warnings.append(
            ChemicalWarningRow(
                code="chemical_induced_calcite_risk",
                message=f"Base-equivalent dosing increased calcite SI by {delta_calcite:.2f}",
            )
        )
    if capacity_mg_l > 0.0 and residual_h2s_mg_l > 0.0:
        model_warnings.append(
            ChemicalWarningRow(
                code="h2s_scavenger_under_capacity",
                message=(
                    "Calculated dissolved-H2S load exceeds the practical active scavenger capacity"
                ),
            )
        )

    return ProductionChemicalScaleScenarioResult(
        treated_ph=treated_ph,
        treated_bicarbonate=from_si(treated_bicarbonate_mg_l * 1e-3, "kg/m**3"),
        base_dose_equivalents=from_si(equivalents_mol_l * 1000.0, "mol/m**3"),
        scavenger_capacity=from_si(capacity_mg_l * 1e-3, "kg/m**3"),
        residual_dissolved_h2s=from_si(residual_h2s_mg_l * 1e-3, "kg/m**3"),
        baseline_si=baseline_si,
        treated_si=treated_si,
        si_change=si_change,
        model_warnings=tuple(model_warnings),
        warnings=tuple(warnings),
    )


def _apply_base_dose(
    ph: float,
    bicarbonate_mg_l: float,
    equivalents_mol_l: float,
    temperature_k: float,
    warnings: list[ChemicalWarningRow],
) -> tuple[float, float]:
    if bicarbonate_mg_l <= 0.0 or equivalents_mol_l == 0.0:
        return ph, max(0.0, bicarbonate_mg_l) / (MW_HCO3 * 1000.0)
    k1 = _carbonate_k1(temperature_k)
    k2 = _carbonate_k2(temperature_k)
    initial_hco3 = bicarbonate_mg_l / (MW_HCO3 * 1000.0)
    initial_alpha1 = _carbonate_fractions(ph, k1, k2)[1]
    if initial_alpha1 < 1.0e-12:
        warnings.append(
            ChemicalWarningRow(
                code="carbonate_balance_ill_conditioned",
                message="Initial pH/HCO3 combination cannot define total inorganic carbon",
            )
        )
        return ph, initial_hco3
    total_carbon = initial_hco3 / initial_alpha1
    target = _alkalinity(total_carbon, ph, k1, k2) + equivalents_mol_l

    low, high = 2.0, 14.0
    for _ in range(100):
        mid = 0.5 * (low + high)
        if _alkalinity(total_carbon, mid, k1, k2) < target:
            low = mid
        else:
            high = mid
    new_ph = 0.5 * (low + high)
    if new_ph > 13.99:
        warnings.append(
            ChemicalWarningRow(
                code="ph_solver_upper_bound",
                message="Dose exceeds the closed-carbonate model pH 14 search bound",
            )
        )
    return new_ph, total_carbon * _carbonate_fractions(new_ph, k1, k2)[1]
