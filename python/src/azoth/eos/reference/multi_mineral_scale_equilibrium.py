"""``eos.multi_mineral_scale_equilibrium`` - NeqSim's ``MultiMineralScaleEquilibrium``.

The coupled solid/liquid equilibrium of five minerals over one brine, by the class's greedy
coordinate descent: the single most-violating mineral is advanced by a relaxed step until no
extent exceeds the threshold, with under-saturated solids allowed to redissolve.

Spec: ``specs/models/eos/multi_mineral_scale_equilibrium.toml``.
"""

from __future__ import annotations

import math
from typing import Any

from azoth.core.errors import InvalidInputError
from azoth.core.range import apply_checks, checks_for
from azoth.core.record_outputs_gen import ActivityCoefficientRow, MineralEquilibriumRow
from azoth.core.result import MultiMineralScaleEquilibriumResult
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning, WarningCode
from azoth.eos.reference.pitzer_scale_activity import ion_coefficients
from azoth.eos.reference.scale_prediction import scale_prediction as _predictor

CALC_ID = "eos.multi_mineral_scale_equilibrium"

EPS = 1.0e-18
DH_B = 0.3283

MW_BASO4 = 233.39
MW_SRSO4 = 183.68
MW_CASO4 = 136.14
MW_CACO3 = 100.09
MW_FECO3 = 115.85

CA, BA, SR, FE, SO4, CO3 = range(6)

#: The nine-ion index of each free-concentration position, for the Pitzer products.
PITZER_INDEX = {CA: 2, BA: 3, SR: 4, FE: 5, SO4: 6, CO3: 7}

ACTIVITY_MODELS = ("davies", "bdot", "pitzer_binary")


def _mol_per_l(quantity: Q) -> float:
    return quantity.to("mol/m**3").magnitude / 1000.0


def multi_mineral_scale_equilibrium(
    calcium: Q,
    barium: Q,
    strontium: Q,
    iron: Q,
    magnesium: Q,
    sodium: Q,
    bicarbonate: Q,
    sulphate: Q,
    total_dissolved_solids: Q,
    temperature: Q,
    pressure: Q,
    co2_partial_pressure: Q,
    ph: float,
    auto_ph: bool,
    second_order_pressure_correction: bool,
    activity_model: str,
    bdot_ion_size: float | None = None,
    bdot_b: float | None = None,
    pitzer_ionic_strength_molal: float | None = None,
    water_volume: Q | None = None,
    max_iterations: float | None = None,
    relaxation: float | None = None,
    convergence_threshold: Q | None = None,
) -> MultiMineralScaleEquilibriumResult:
    """The coupled multi-mineral solid/liquid equilibrium of a water analysis.

    Args:
        calcium: calcium ion concentration.
        barium: barium ion concentration.
        strontium: strontium ion concentration.
        iron: ferrous iron concentration.
        magnesium: magnesium ion concentration.
        sodium: sodium ion concentration.
        bicarbonate: bicarbonate concentration.
        sulphate: sulphate concentration.
        total_dissolved_solids: TDS, which sets the predictor's ionic strength.
        temperature: absolute temperature.
        pressure: absolute pressure.
        co2_partial_pressure: CO2 partial pressure, read by the auto-pH path.
        ph: the water's pH, ignored when ``auto_ph`` is on and its gate holds.
        auto_ph: whether to estimate the pH with the class's own Henry constant.
        second_order_pressure_correction: whether the Ksp correction carries it.
        activity_model: ``davies``, ``bdot`` or ``pitzer_binary``.
        bdot_ion_size: the B-dot ion size; ``None`` takes 5 Angstrom.
        bdot_b: the B-dot coefficient; ``None`` takes 0.041.
        pitzer_ionic_strength_molal: the Pitzer strength; ``None`` takes the fallback.
        water_volume: the volume the grams are over; ``None`` takes one litre.
        max_iterations: the coordinate-descent cap; ``None`` takes 200000.
        relaxation: the under-relaxation; ``None`` takes 0.5.
        convergence_threshold: the stopping extent; ``None`` takes 1e-12 mol/m3.

    Returns:
        The assemblage, the residuals and the coefficients the solve used.

    Raises:
        InvalidInputError: for an activity-model name outside the three.
        OutOfRangeError: if an input is outside its declared range.

    Example:
        >>> import azoth
        >>> q = azoth.ureg.Quantity
        >>> r = multi_mineral_scale_equilibrium(
        ...     q(0.0, "kg/m**3"), q(1.0, "kg/m**3"), q(0.0, "kg/m**3"),
        ...     q(0.0, "kg/m**3"), q(0.0, "kg/m**3"), q(0.0, "kg/m**3"),
        ...     q(0.0, "kg/m**3"), q(0.7, "kg/m**3"), q(0.02, "kg/m**3"),
        ...     q(298.15, "K"), q(1.013e5, "Pa"), q(0.0, "Pa"), 7.0, False, False,
        ...     "davies", water_volume=q(0.001, "m**3"))
        >>> abs(r.total_mass.magnitude - 1.6973864119588732) < 1e-9
        True
    """
    from azoth._models_gen import model as _model_for

    spec: dict[str, Any] = _model_for(CALC_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    values = {
        name: input_to_si(spec, name, value)
        for name, value in (
            ("calcium", calcium),
            ("barium", barium),
            ("strontium", strontium),
            ("iron", iron),
            ("magnesium", magnesium),
            ("sodium", sodium),
            ("bicarbonate", bicarbonate),
            ("sulphate", sulphate),
            ("total_dissolved_solids", total_dissolved_solids),
            ("temperature", temperature),
            ("pressure", pressure),
            ("co2_partial_pressure", co2_partial_pressure),
            ("water_volume", water_volume),
            ("convergence_threshold", convergence_threshold),
            ("bdot_ion_size", bdot_ion_size),
            ("bdot_b", bdot_b),
            ("pitzer_ionic_strength_molal", pitzer_ionic_strength_molal),
        )
        if value is not None
    }
    values["ph"] = ph
    if relaxation is not None:
        values["relaxation"] = relaxation
    if convergence_threshold is not None:
        threshold_mol_l = values["convergence_threshold"] * 1000.0
    else:
        threshold_mol_l = 1.0e-15
    apply_checks(checks.on_input, values.get, warnings)

    if activity_model not in ACTIVITY_MODELS:
        raise InvalidInputError(
            "activity_model", f"`{activity_model}` is not one of {', '.join(ACTIVITY_MODELS)}"
        )

    ion_size = values.get("bdot_ion_size", 5.0)
    b_dot = values.get("bdot_b", 0.041)
    iterations_cap = 200000 if max_iterations is None else int(max_iterations)
    relax = 0.5 if relaxation is None else relaxation
    litres = 1.0 if water_volume is None else values["water_volume"] * 1000.0
    temperature_k = values["temperature"]

    state = _predictor(
        calcium,
        barium,
        strontium,
        iron,
        magnesium,
        sodium,
        bicarbonate,
        sulphate,
        total_dissolved_solids,
        temperature,
        pressure,
        co2_partial_pressure,
        ph,
        auto_ph,
        second_order_pressure_correction,
    )

    gammas: tuple[float, ...] = ()
    if activity_model == "pitzer_binary":
        sodium_mol_l = _mol_per_l(state.total_sodium)
        strength = values.get(
            "pitzer_ionic_strength_molal",
            sodium_mol_l
            if sodium_mol_l > 0.0
            else max(0.0, values["total_dissolved_solids"] * 1000.0 / 58440.0),
        )
        gammas = ion_coefficients(temperature_k, state.debye_huckel_a, strength)
        gamma2 = gammas[2]
        coefficients = [
            ActivityCoefficientRow(ion=name, activity_coefficient=gammas[index])
            for name, index in (
                ("Na+", 0),
                ("Cl-", 1),
                ("Ca++", 2),
                ("Ba++", 3),
                ("Sr++", 4),
                ("Fe++", 5),
                ("SO4--", 6),
                ("CO3--", 7),
                ("HCO3-", 8),
            )
        ]
    elif activity_model == "bdot":
        strength = max(_mol_per_l(state.ionic_strength), 0.0)
        sqrt_i = math.sqrt(strength)
        log_gamma = (
            -state.debye_huckel_a * 4.0 * sqrt_i / (1.0 + ion_size * DH_B * sqrt_i)
            + b_dot * strength
        )
        gamma2 = math.pow(10.0, log_gamma)
        coefficients = [ActivityCoefficientRow(ion="divalent", activity_coefficient=gamma2)]
    else:
        gamma2 = state.divalent_activity_coefficient
        coefficients = [ActivityCoefficientRow(ion="divalent", activity_coefficient=gamma2)]

    def product(cation: int, anion: int) -> float:
        if activity_model != "pitzer_binary":
            return gamma2 * gamma2
        return float(gammas[PITZER_INDEX[cation]]) * float(gammas[PITZER_INDEX[anion]])

    free = [0.0] * 6
    free[CA] = max(_mol_per_l(state.free_calcium), 0.0)
    free[BA] = max(_mol_per_l(state.total_barium), 0.0)
    free[SR] = max(_mol_per_l(state.total_strontium), 0.0)
    free[FE] = max(_mol_per_l(state.total_iron), 0.0)
    free[SO4] = max(_mol_per_l(state.free_sulphate), 0.0)
    free[CO3] = max(_mol_per_l(state.free_carbonate), 0.0)
    initial_free = list(free)

    minerals: list[dict[str, Any]] = []

    def register(
        cation: int, anion: int, ksp: float, molar_mass: float, name: str, common: str
    ) -> None:
        if free[cation] <= EPS or free[anion] <= EPS or not ksp > 0.0:
            return
        minerals.append(
            {
                "cation": cation,
                "anion": anion,
                "ksp": ksp,
                "molar_mass": molar_mass,
                "name": name,
                "common": common,
                "product": product(cation, anion),
                "initial_si": 0.0,
                "mol": 0.0,
            }
        )

    register(BA, SO4, state.ksp_barite, MW_BASO4, "BaSO4", "barite")
    register(SR, SO4, state.ksp_celestite, MW_SRSO4, "SrSO4", "celestite")
    register(CA, SO4, state.ksp_anhydrite, MW_CASO4, "CaSO4", "anhydrite")
    register(CA, CO3, state.ksp_calcite, MW_CACO3, "CaCO3", "calcite")
    register(FE, CO3, state.ksp_siderite, MW_FECO3, "FeCO3", "siderite")

    for mineral in minerals:
        mineral["initial_si"] = _saturation_index(
            free[mineral["cation"]], free[mineral["anion"]], mineral["product"], mineral["ksp"]
        )

    iterations = 0
    stopped_before_limit = False
    for _ in range(iterations_cap):
        best: dict[str, Any] | None = None
        best_abs = 0.0
        best_dxi = 0.0
        for mineral in minerals:
            cation = free[mineral["cation"]]
            anion = free[mineral["anion"]]
            k_prime = mineral["ksp"] / mineral["product"]
            b = cation + anion
            c = cation * anion - k_prime
            disc = max(b * b - 4.0 * c, 0.0)
            xi_eq = 0.5 * (b - math.sqrt(disc))
            if xi_eq > 0.0:
                cap = min(cation, anion) - EPS
                dxi = max(min(xi_eq, cap), 0.0)
            elif xi_eq < 0.0 and mineral["mol"] > EPS:
                dxi = max(xi_eq, -mineral["mol"])
            else:
                dxi = 0.0
            if abs(dxi) > best_abs:
                best_abs = abs(dxi)
                best_dxi = dxi
                best = mineral
        if best is None or best_abs < threshold_mol_l:
            stopped_before_limit = True
            break
        step = relax * best_dxi
        free[best["cation"]] -= step
        free[best["anion"]] -= step
        best["mol"] += step
        if best["mol"] < 0.0:
            best["mol"] = 0.0
        if free[best["cation"]] < 0.0:
            free[best["cation"]] = 0.0
        if free[best["anion"]] < 0.0:
            free[best["anion"]] = 0.0
        iterations += 1
    limit_reached = not stopped_before_limit and iterations >= iterations_cap
    if limit_reached:
        warnings.append(
            Warning(
                code=WarningCode.SOLVER_NOT_CONVERGED,
                message=(
                    "the coordinate descent reached its iteration cap; the answer is the "
                    "last iterate"
                ),
                field=None,
            )
        )

    rows: list[MineralEquilibriumRow] = []
    total_mass_mg_l = 0.0
    violation = 0.0
    for mineral in minerals:
        mass_mg_l = mineral["mol"] * mineral["molar_mass"] * 1000.0
        final_si = _saturation_index(
            free[mineral["cation"]], free[mineral["anion"]], mineral["product"], mineral["ksp"]
        )
        total_mass_mg_l += mass_mg_l
        limiting = (
            free[mineral["cation"]] <= 10.0 * EPS or free[mineral["anion"]] <= 10.0 * EPS
        ) and mineral["mol"] > EPS
        rows.append(
            MineralEquilibriumRow(
                name=mineral["name"],
                common_name=mineral["common"],
                initial_saturation_index=mineral["initial_si"],
                final_saturation_index=final_si,
                activity_coefficient_product=mineral["product"],
                precipitated_amount=from_si(mineral["mol"] * 1000.0, "mol/m**3"),
                precipitated_mass=from_si(mass_mg_l * 1.0e-3, "kg/m**3"),
                total_precipitated_mass=from_si(mass_mg_l * 1.0e-6 * litres, "kg"),
                limiting_ion_exhausted=limiting,
            )
        )
        residual = abs(final_si) if mineral["mol"] > EPS else max(final_si, 0.0)
        if not math.isfinite(residual):
            residual = math.inf
        violation = max(violation, residual)

    precipitated_by_ion = [0.0] * 6
    for mineral in minerals:
        precipitated_by_ion[mineral["cation"]] += mineral["mol"]
        precipitated_by_ion[mineral["anion"]] += mineral["mol"]
    ion_balance = max(
        abs(initial_free[ion] - free[ion] - precipitated_by_ion[ion]) for ion in range(6)
    )

    return MultiMineralScaleEquilibriumResult(
        minerals=tuple(rows),
        total_mass=from_si(total_mass_mg_l * 1.0e-3, "kg/m**3"),
        activity_coefficients=tuple(coefficients),
        divalent_activity_coefficient=gamma2,
        iterations=iterations,
        iteration_limit_reached=limit_reached,
        maximum_complementarity_violation=violation,
        maximum_ion_balance_residual=from_si(ion_balance * 1000.0, "mol/m**3"),
        residual_calcium=from_si(free[CA] * 1000.0, "mol/m**3"),
        residual_barium=from_si(free[BA] * 1000.0, "mol/m**3"),
        residual_strontium=from_si(free[SR] * 1000.0, "mol/m**3"),
        residual_ferrous=from_si(free[FE] * 1000.0, "mol/m**3"),
        residual_sulphate=from_si(free[SO4] * 1000.0, "mol/m**3"),
        residual_carbonate=from_si(free[CO3] * 1000.0, "mol/m**3"),
        warnings=tuple(warnings),
    )


def _saturation_index(c_cation: float, c_anion: float, product: float, ksp: float) -> float:
    iap = product * max(c_cation, EPS) * max(c_anion, EPS)
    if iap <= 0.0 or ksp <= 0.0:
        return math.nan
    return math.log10(iap / ksp)
