"""``eos.hydrate_fraction``: the Python reference.

The mirror of ``crates/azoth-eos/src/hydrate_fraction.rs``:

```text
f_i = z_i - beta w_i                   the fluid the hydrate leaves, as mole numbers
solve  beta = min_i z_i / w_i          the fraction at which a component is exhausted
```

NeqSim's ``TPHydrateFlash`` computes the same quantity and **its state is not one**: the
hydrate it leaves has mole fractions summing to ``1.1201``, and the phases hold ``+0.0874``
water and ``-0.0734`` methane against the feed. The divergence is recorded upstream and in
the spec; this is the solve that closes.

What decides between no hydrate and all of the water is the same objective
:mod:`azoth.eos.reference.hydrate_formation_temperature` drives to zero, read at the **feed**:
the objective is monotone in the fraction, so a feed above its formation temperature has no
hydrate at all, and one below it takes the whole of the water there is.
"""

from __future__ import annotations

import math

from azoth.core.errors import InvalidInputError, SolverNotConvergedError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import HydrateFractionResult, HydrateStructure, PtFlashResult
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning
from azoth.eos.mixture import Mixture
from azoth.eos.reference import _hydrate
from azoth.eos.reference._mixture_state import phase_state, reduced_parameters
from azoth.eos.reference.pt_flash import pt_flash

MODEL_ID = "eos.hydrate_fraction"


def _reference_water_fugacity(eos: str, t: float, p: float) -> float:
    """The fugacity of pure water at a state on the cubic the fluid runs, in Pa."""
    from azoth.eos import components as databank

    water, _ = databank.mixture_of(["water"], eos=eos)
    reduced = reduced_parameters(water, t, p)
    state = phase_state(reduced, water.kij, [1.0], liquid=True)
    return math.exp(state.ln_phi[0]) * p


def _guest_fugacities(flash: PtFlashResult, fluid: list[float], p: float) -> list[float]:
    """The guests' fugacities at a flash, in Pa, from the phase ``setFug`` reads."""
    composition = list(flash.y) if flash.phase == "two_phase" else list(fluid)
    coefficients = flash.ln_phi_liquid if flash.phase == "all_liquid" else flash.ln_phi_vapour
    return [
        fraction * math.exp(coefficient) * p
        for fraction, coefficient in zip(composition, coefficients, strict=True)
    ]


def _fluid_at(z: list[float], in_hydrate: list[float], beta: float) -> list[float]:
    """The fluid the hydrate leaves at one fraction, as mole fractions.

    ``f_i = z_i - beta w_i`` is that fluid's own mole numbers, summing to ``1 - beta``. A
    component the hydrate has taken all of comes out at zero rather than negative, which is
    the state the bound is *defined* by.
    """
    fluid = [max(total - beta * taken, 0.0) for total, taken in zip(z, in_hydrate, strict=True)]
    total = sum(fluid)
    return [fraction / total for fraction in fluid]


def _maximum_fraction(z: list[float], in_hydrate: list[float]) -> float:
    """The largest fraction a feed's own composition allows: ``min_i z_i / w_i``."""
    ratios = [total / taken for total, taken in zip(z, in_hydrate, strict=True) if taken > 1.0e-12]
    return min(ratios) if ratios else math.inf


def hydrate_fraction(
    components: list[str], T: Q, P: Q, z: list[float], eos: str = "srk"
) -> HydrateFractionResult:
    """The fraction of a feed that is hydrate at a temperature and pressure.

    **The components cross by name**, because the hydrate's guest tables are keyed by name
    and a :class:`~azoth.eos.mixture.Component` carries none.

    Args:
        components: the substances, by name. At least one must be water and one a former.
        T: absolute temperature; it must be below the formation temperature, or the answer
            is zero.
        P: absolute pressure.
        z: the overall mole fractions. Checked rather than renormalised.
        eos: the cubic the fluid runs, which is also the one the hydrate's reference water
            phase is built from.

    Raises:
        InvalidInputError: if the fluid has no water, or nothing in it occupies a cage.
        OutOfRangeError: if ``T``, ``P`` or a trial's cavity sum has no value.
        SolverNotConvergedError: if the fraction and the composition do not settle.
    """
    from azoth import _models_gen
    from azoth.eos import components as databank

    spec = _models_gen.model(MODEL_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []
    t_si = input_to_si(spec, "T", T)
    p_si = input_to_si(spec, "P", P)
    apply_checks(checks.on_input, {"T": t_si, "P": p_si}.get, warnings)

    hydration = _hydrate.hydration_for(components)
    if hydration.water_index is None:
        raise InvalidInputError("components", "the hydration carries no water index")
    water_index = hydration.water_index
    mixture, _ = databank.mixture_of(components, eos=eos)

    algorithm = spec["algorithm"]
    tolerance = float(algorithm["tolerance"])
    max_iterations = int(algorithm["max_iterations"])

    def cages(fluid: list[float]) -> tuple[list[float], int, float, list[float]]:
        """One trial: flash the fluid, fill the cages by what the flash leaves."""
        flash = pt_flash(mixture, from_si(t_si, "K"), from_si(p_si, "Pa"), fluid)
        fugacities = _guest_fugacities(flash, fluid, p_si)
        reference = _reference_water_fugacity(eos, t_si, p_si)
        structure, coefficient = _hydrate.stable_structure(
            hydration.guests, fugacities, t_si, p_si, reference
        )
        return (
            _hydrate.composition(hydration.guests, fugacities, structure, t_si, water_index),
            structure,
            coefficient,
            fugacities,
        )

    # **The feed's own state decides whether there is any hydrate at all.** The objective is
    # only read here: at the bound the fluid holds no water, so the ratio has no value at it,
    # and the fixed point below needs the cages rather than the residual.
    iterations = 1
    composition, structure, coefficient, fugacities = cages(list(z))
    residual = math.log(coefficient * p_si / fugacities[water_index])
    if residual > 0.0:
        return _finish(
            mixture, t_si, p_si, z, composition, structure, 0.0, residual, iterations, warnings
        )

    # **Below it the hydrate takes the water there is**, and the fraction and the cages are
    # one fixed point: the fraction is ``min_i z_i/w_i`` from the cages' own composition, and
    # the cages are filled from the fluid that fraction leaves.
    beta = _maximum_fraction(z, composition)
    for step in range(2, max_iterations + 1):
        if not math.isfinite(beta):
            break
        fluid = _fluid_at(z, composition, beta)
        composition, structure, _coefficient, _fugacities = cages(fluid)
        iterations = step
        following = _maximum_fraction(z, composition)
        if abs(following - beta) < tolerance:
            return _finish(
                mixture,
                t_si,
                p_si,
                z,
                composition,
                structure,
                following,
                residual,
                iterations,
                warnings,
            )
        beta = following

    raise SolverNotConvergedError(iterations=iterations, residual=beta, tolerance=tolerance)


def _finish(
    mixture: Mixture,
    t_si: float,
    p_si: float,
    z: list[float],
    composition: list[float],
    structure: int,
    beta: float,
    residual: float,
    iterations: int,
    warnings: list[Warning],
) -> HydrateFractionResult:
    """The answer, with the material balance it exists to keep measured and reported."""
    fluid = _fluid_at(z, composition, beta)
    flash = pt_flash(mixture, from_si(t_si, "K"), from_si(p_si, "Pa"), fluid)
    split = flash.beta if flash.beta is not None else 0.0
    if flash.phase == "two_phase":
        x, y = flash.x, flash.y
    else:
        x = y = tuple(fluid)
    balance_error = 0.0
    for i, total in enumerate(z):
        left = ((1.0 - split) * x[i] + split * y[i]) * (1.0 - beta)
        balance_error = max(balance_error, abs(beta * composition[i] + left - total))

    return HydrateFractionResult(
        beta=beta,
        structure=(
            HydrateStructure.STRUCTURE_I if structure == 0 else HydrateStructure.STRUCTURE_II
        ),
        balance_error=balance_error,
        iterations=iterations + 1,
        residual=residual,
        warnings=tuple(warnings),
    )


__all__ = ["hydrate_fraction"]
