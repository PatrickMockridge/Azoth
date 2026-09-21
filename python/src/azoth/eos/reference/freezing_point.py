"""``eos.freezing_point`` - the temperature at which a pure substance's solid appears.

The Python reference, the mirror of ``crates/azoth-eos/src/freezing_point.rs``:

```text
solve  (g_fluid(T, P) - g_solid(T, P)) / (R T) = 0   for T
```

Two things the Rust module's docstring states and a reader of *this* file needs in the same
breath, because they are where a port goes wrong silently:

* **the solid is calibrated, not raw.** ``eos.parahydrogen_solid_phase`` is the thesis
  equation as it stands and its reference is fixed elsewhere - by NeqSim's
  ``SystemLeachmanEos``, which shifts it to the para-Leachman liquid at the triple point.
  Unshifted it misses the triple point's own definition, the two Gibbs energies equal there,
  by ``3.12`` J/mol.
* **the fluid's root is NeqSim's rule**: a gas below the triple-point pressure and a liquid at
  and above it, which are different states at the same temperature and pressure.
"""

from __future__ import annotations

import math
from collections.abc import Callable
from typing import Any

from azoth.core.errors import (
    AzothError,
    InvalidInputError,
    OutOfRangeError,
    SolverNotConvergedError,
)
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import FreezingPointResult, Phase
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning
from azoth.eos.reference import hydrogen_phase as fluid
from azoth.eos.reference.parahydrogen_solid_phase import parahydrogen_solid_phase
from azoth.eos.reference.pt_flash import pt_flash
from azoth.eos.reference.tp_solid_flash import METHANE_NEVER_FREEZES, _tabulated_solid_fugacity

MODEL_ID = "eos.freezing_point"

#: NeqSim's ``ThermodynamicConstantsInterface.R``, which its residual divides by. **Not** the
#: equation's own: the two reference equations carry ``8.31451`` and ``8.3144621``, each a
#: rounded value of its own era, and this is the third.
R = 8.3144621

#: The tolerance on the dimensionless residual and on the bracket.
RESIDUAL_TOLERANCE = 1.0e-10
BRACKET_TOLERANCE = 1.0e-12

#: How far the search may walk, in kelvin, and how many steps it may take.
MINIMUM_TEMPERATURE = 0.5
MAXIMUM_TEMPERATURE = 5000.0
SPAN_GROWTH = 1.8
MAXIMUM_BRACKET_STEPS = 50
MAXIMUM_SOLVER_STEPS = 100

#: The triple point, and the enthalpy of fusion there: the coordinates the calibration is
#: measured at. Para-hydrogen's, from NeqSim's ``ParaHydrogenSolidHelmholtzEquation``.
TRIPLE_POINT_TEMPERATURE = 13.8033
TRIPLE_POINT_PRESSURE = 7042.0
TRIPLE_POINT_ENTHALPY_OF_FUSION = 118.0

#: The substances this has a solid equation for. One name, because
#: ``thermo/util/solid/`` holds two equations and the other is argon's.
SUBSTANCES = ("para-hydrogen",)


def _solid_state(t: float, p: float) -> tuple[float, float]:
    """The raw solid's molar Gibbs energy and entropy at a state, in J/mol and J/(mol*K)."""
    state = parahydrogen_solid_phase(T=from_si(t, "K"), P=from_si(p, "Pa"))
    return state.g.to("J/mol").magnitude, state.s.to("J/(mol*K)").magnitude


def _fluid_state(t: float, p: float, root: str) -> tuple[float, float]:
    """The fluid's molar Gibbs energy and entropy, on the root named, at a state."""
    density = fluid._solve_density(t, p, "para") if root == "gas" else _dense_density(t, p)
    properties = fluid._properties(t, density, "para")
    return properties["g"], properties["s"]


def _dense_density(t: float, p: float) -> float:
    density = fluid._solve_density_dense(t, p, "para")
    if density is None:
        raise OutOfRangeError(
            "P",
            p,
            "the freezing search asks for the fluid's dense root, and this pressure has none "
            "at this temperature within the range the equation is fitted to",
        )
    return density


def calibration() -> tuple[float, float]:
    """The two constants NeqSim's ``SystemLeachmanEos`` computes at the triple point.

    ``(gibbs_shift, entropy_shift)``, in J/mol and J/(mol*K). Measurements of the two models
    rather than tunables: the first comes out ``-3.122683``, which is the offset the capture
    measures between NeqSim's solid and the raw one.
    """
    g_liquid, s_liquid = _fluid_state(TRIPLE_POINT_TEMPERATURE, TRIPLE_POINT_PRESSURE, "liquid")
    g_raw, s_raw = _solid_state(TRIPLE_POINT_TEMPERATURE, TRIPLE_POINT_PRESSURE)
    return (
        g_liquid - g_raw,
        s_liquid - s_raw - TRIPLE_POINT_ENTHALPY_OF_FUSION / TRIPLE_POINT_TEMPERATURE,
    )


def calibrated_solid_gibbs(t: float, p: float, shifts: tuple[float, float]) -> float:
    """The solid's molar Gibbs energy at a state, on the reference the liquid is on."""
    gibbs_shift, entropy_shift = shifts
    g_raw, _ = _solid_state(t, p)
    # A shift of `-S (T - T_tp) + G` in the Helmholtz energy is the same shift in the Gibbs
    # energy: it carries no volume, so `g = a + P v` passes it through.
    return g_raw - entropy_shift * (t - TRIPLE_POINT_TEMPERATURE) + gibbs_shift


def residual(t: float, p: float, root: str, shifts: tuple[float, float]) -> float:
    """The dimensionless residual the solve drives to zero."""
    g_fluid, _ = _fluid_state(t, p, root)
    return (g_fluid - calibrated_solid_gibbs(t, p, shifts)) / (R * t)


def _solve(
    residual_fn: Callable[[float], float], start: float, algorithm: Any
) -> tuple[float, int, float]:
    """The temperature a residual crosses zero at, by NeqSim's own walk.

    The search starts where the algorithm says and expands a span either side of it: the answer
    is a property of the substance and the pressure, so the start is a path and not a state. A
    trial the equation cannot evaluate is skipped rather than fatal, which is NeqSim's
    ``evaluateResidualIfValid`` - near the triple point one side of a span can ask for a dense
    root that does not exist.

    Raises:
        SolverNotConvergedError: if no span brackets a sign change, or the bracket collapses
            without reaching the tolerance.
    """
    start_residual = residual_fn(start)
    if abs(start_residual) < RESIDUAL_TOLERANCE:
        return start, 1, start_residual

    span = max(start * 0.01, 0.1)
    low, low_residual = start, start_residual
    # The bracket's upper side is only read through its sign; it starts at the same point
    # because a bracket needs two sides and either can move first.
    high, high_residual = start, start_residual
    iterations = 1
    bracketed = False
    for _ in range(MAXIMUM_BRACKET_STEPS):
        candidate_low = max(start - span, MINIMUM_TEMPERATURE)
        candidate_high = min(start + span, MAXIMUM_TEMPERATURE)
        # **Any error, not just an out-of-range one**, because that is what the Rust kernel
        # catches: its `if let Ok(value) = residual(...)` discards every `Err`. The tabulated
        # route's residual is a whole flash, so it can fail as `SolverNotConverged` too.
        try:
            low, low_residual = candidate_low, residual_fn(candidate_low)
            iterations += 1
        except AzothError:
            pass
        try:
            high, high_residual = candidate_high, residual_fn(candidate_high)
            iterations += 1
        except AzothError:
            pass
        if low_residual * high_residual <= 0.0:
            bracketed = True
            break
        span *= SPAN_GROWTH

    if not bracketed:
        raise _not_converged(iterations, start_residual, algorithm)

    for _ in range(MAXIMUM_SOLVER_STEPS):
        mid = 0.5 * (low + high)
        value = residual_fn(mid)
        iterations += 1
        if abs(value) < RESIDUAL_TOLERANCE or high - low < BRACKET_TOLERANCE:
            return mid, iterations, value
        if low_residual * value <= 0.0:
            high = mid
        else:
            low, low_residual = mid, value

    raise _not_converged(iterations, 0.5 * (low + high), algorithm)


def _log(x: float) -> float:
    """``log``, giving NaN outside its domain as Java and Rust do."""
    try:
        return math.log(x)
    except ValueError:
        return math.nan


def _tabulated_residual(
    mixture: Any, z: list[float], index: int, solid: str, eos: str, t: float, p: float
) -> float:
    """NeqSim's ``calculateLogEquilibriumResidual`` outside its Helmholtz branch.

    ``residual = ln z_k - ln( sum_p beta_p phi_solid_k / phi_kp )``, the multiphase appearance
    condition over the phases the flash found. The largest logarithm is taken out before
    exponentiating: the terms are ratios of fugacity coefficients spanning twenty orders of
    magnitude and a direct sum underflows.
    """
    if not z[index] > 0.0:
        raise InvalidInputError(
            "solid",
            f"`{solid}` has no overall mole fraction, and a solid cannot appear from a "
            f"substance the feed does not have",
        )
    flash = pt_flash(mixture, from_si(t, "K"), from_si(p, "Pa"), z)
    # **Methane never freezes in NeqSim**: `ComponentSolid.fugcoef` returns `1e30` for it before
    # any arithmetic, so its residual has no sign change and the search refuses.
    if solid.strip().lower() == "methane":
        phi_solid = METHANE_NEVER_FREEZES
    else:
        phi_solid = _tabulated_solid_fugacity(
            mixture, index, solid, t, p, from_si(t, "K"), from_si(p, "Pa"), eos
        )
    # **Not `math.log`: Python raises where Java and Rust give NaN.** NeqSim's
    # `ComponentSolid.fugcoef` can return a negative coefficient - the class makes one for
    # water by dividing two signed expressions - and Rust's `f64::ln` makes that a `NaN` the
    # search then skips or refuses on, while `math.log` aborts the twin at the same state. The
    # same helper, and the same reason, as `saft_vr_mie_phase._log`.
    ln_phi_solid = _log(phi_solid)

    # The fluid's phases, as `(amount, ln phi of this component)`. **The match is on the
    # phase and not on whether a beta came back**, which is the Rust kernel's order and not an
    # equivalent: a flash can report `all_liquid` with a beta still set, and the two orders
    # then read a different phase's coefficient.
    terms: list[float] = []
    if flash.phase is Phase.TWO_PHASE:
        raw = flash.beta if flash.beta is not None else 0.0
        terms.append(_log(1.0 - raw) + ln_phi_solid - flash.ln_phi_liquid[index])
        terms.append(_log(raw) + ln_phi_solid - flash.ln_phi_vapour[index])
    elif flash.phase is Phase.ALL_LIQUID:
        terms.append(ln_phi_solid - flash.ln_phi_liquid[index])
    else:
        terms.append(ln_phi_solid - flash.ln_phi_vapour[index])
    maximum = max(terms)
    if not math.isfinite(maximum):
        raise OutOfRangeError(
            "P",
            p,
            "no phase of the flashed fluid has a finite fugacity contribution at this state",
        )
    scaled = sum(math.exp(term - maximum) for term in terms)
    return _log(z[index]) - maximum - _log(scaled)


def freezing_point(components: list[str], z: list[float], solid: str, P: Q) -> FreezingPointResult:
    """The temperature at which a fluid's solid-forming substance freezes.

    Two routes, as NeqSim has two: ``para-hydrogen`` against its calibrated solid Helmholtz
    equation, every other substance against the tabulated solid whose residual is the multiphase
    appearance condition.

    Raises:
        InvalidInputError: if ``z`` does not match ``components``, if ``solid`` is not one of
            them, if ``P`` is not positive, or if the tabulated route's candidate carries no
            melt data.
        OutOfRangeError: from a trial's state.
        SolverNotConvergedError: if no candidate's residual brackets a sign change.
    """
    from azoth import _models_gen

    spec = _models_gen.model(MODEL_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []
    p_si = input_to_si(spec, "P", P)
    apply_checks(checks.on_input, {"P": p_si}.get, warnings)

    if len(components) != len(z):
        raise InvalidInputError(
            "z",
            f"{len(components)} components and {len(z)} mole fractions; a candidate's own "
            f"fraction is what decides whether it can appear",
        )
    names = [name.strip().lower() for name in components]
    wanted = solid.strip().lower()
    if wanted not in names:
        raise InvalidInputError(
            "solid",
            f"{solid!r} is not one of the fluid's components, so it is not a substance whose "
            f"freezing point this can solve for",
        )
    index = names.index(wanted)

    algorithm = spec["algorithm"]
    start = float(algorithm.get("initial_temperature", 14.0))

    if wanted == "para-hydrogen":
        # NeqSim's own rule, from `initializeCoexistingFluidPhase`: the fluid is a gas below
        # the triple-point pressure and a liquid at and above it.
        root = "gas" if p_si < TRIPLE_POINT_PRESSURE else "liquid"
        shifts = calibration()
        temperature, iterations, value = _solve(
            lambda t: residual(t, p_si, root, shifts), start, algorithm
        )
    else:
        from azoth.eos import components as databank

        mixture, _ = databank.mixture_of(components, eos="srk")
        temperature, iterations, value = _solve(
            lambda t: _tabulated_residual(mixture, z, index, solid, "srk", t, p_si),
            start,
            algorithm,
        )

    return FreezingPointResult(
        temperature=from_si(temperature, "K"),
        component=components[index],
        iterations=iterations,
        residual=value,
        warnings=tuple(warnings),
    )


def _not_converged(iterations: int, value: float, algorithm: dict[str, object]) -> Exception:
    return SolverNotConvergedError(
        iterations=iterations,
        residual=value,
        tolerance=float(str(algorithm["tolerance"])),
    )


__all__ = ["freezing_point"]
