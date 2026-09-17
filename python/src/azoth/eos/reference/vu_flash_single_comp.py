"""``eos.vu_flash_single_comp`` - the volume/internal-energy state of a pure component.

Spec: ``specs/models/eos/vu_flash_single_comp.toml``. NeqSim's ``VUflashSingleComp``.

A pure component has no two-phase split for a flash to find. One K-value is either above
one or below it, so :func:`azoth.eos.reference.pt_flash.pt_flash` reports a single phase
at every temperature and pressure and ``eos.vu_flash``'s iteration - a Newton on ``(P, T)``
around a flash - has nothing to converge on. It refuses, and it refuses for every pure
feed in the two-phase region, which is a large part of what a depressurisation
calculation asks for.

The state is not unknown there, it is *constrained*: a pure component at a pressure is on
its saturation line, so the temperature is the saturation temperature and the split is a
lever rule on the two saturated internal energies. That is this model, and it is exact
for a pure component where a flash would be an approximation of it.
"""

from __future__ import annotations

from typing import Any

from azoth import _models_gen
from azoth.core.errors import InvalidInputError, OutOfRangeError, SolverNotConvergedError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import Phase, VuFlashSingleCompResult
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning, WarningCode
from azoth.eos.mixture import Mixture
from azoth.eos.reference._mixture_state import phase_state, reduced_parameters
from azoth.eos.reference.molar_enthalpy_entropy import IdealGasModel, molar_enthalpy_entropy
from azoth.eos.reference.pr_molar_volume import MOLAR_GAS_CONSTANT
from azoth.eos.reference.pure_saturation import pure_saturation

MODEL_ID = "eos.vu_flash_single_comp"


def _saturated(
    mixture: Mixture, ideal_gas: IdealGasModel, t: float, p: float
) -> tuple[tuple[float, float], tuple[float, float]]:
    """The saturated liquid's and vapour's `(u, v)` at a state.

    Both are the *feed's* own composition - there is only one - so the two differ in the
    cubic root alone, which is what "saturated" means for a pure component.
    """
    reduced = reduced_parameters(mixture, t, p)
    x = [1.0]
    liquid = phase_state(reduced, mixture.kij, x, liquid=True)
    vapour = phase_state(reduced, mixture.kij, x, liquid=False)
    out = []
    for z in (liquid.z, vapour.z):
        state = molar_enthalpy_entropy(mixture, ideal_gas, from_si(t, "K"), from_si(p, "Pa"), x, z)
        v = z * MOLAR_GAS_CONSTANT * t / p
        out.append((float(state.h.to("J/mol").magnitude) - p * v, v))
    return out[0], out[1]


def _saturation_temperature(mixture: Mixture, p: float, algorithm: dict[str, Any]) -> float:
    """The temperature at which a pure component's saturation pressure is ``p``."""
    component = mixture.components[0]
    tc = component.Tc.to_base_units().magnitude
    if p >= component.Pc.to_base_units().magnitude:
        raise OutOfRangeError(
            "P",
            p,
            f"the critical pressure is {component.Pc.to('Pa').magnitude:.6e} Pa, so {p:.6e} Pa "
            f"is above the saturation line at every temperature and a pure component there "
            f"has no two-phase state to find - NeqSim's own guard is the same comparison",
        )
    bracket = algorithm["bracket"]
    lo = float(bracket["lower"]) * tc
    hi = float(bracket["upper"]) * tc

    # Whether the saturation pressure at `t` is at most `p`, which is the question the
    # bisection asks. `None` is a temperature where `eos.pure_saturation` cannot be
    # evaluated - only reachable at the top of the bracket, where the roots have
    # coalesced - and there the answer is *below*: the saturation pressure at such a
    # temperature is at or above the critical pressure, and the critical pressure has
    # already been excluded above. So `None` reads as "above", which is one-sided and
    # needs no temperature at which the calculation is assumed to work.
    def at_most(t: float) -> bool | None:
        try:
            result = pure_saturation(component.Tc, component.Pc, component.omega, from_si(t, "K"))
        except SolverNotConvergedError:
            return None
        return float(result.p_sat.to_base_units().magnitude) <= p

    for _ in range(int(algorithm["max_iterations"])):
        mid = 0.5 * (lo + hi)
        if (hi - lo) / mid <= float(algorithm["tolerance"]):
            break
        if at_most(mid) is True:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def vu_flash_single_comp(
    mixture: Mixture, ideal_gas: IdealGasModel, P: Q, V: Q, U: Q
) -> VuFlashSingleCompResult:
    """The temperature and split of a pure component at a pressure and internal energy.

    Raises:
        InvalidInputError: if the mixture is not a single component.
        OutOfRangeError: if ``P`` is at or above the critical pressure, or ``U`` lies
            outside the two saturated internal energies at it - both are states this
            model has no split for rather than ones it failed to find.

    Example:
        >>> import azoth
        >>> from azoth.eos.components import mixture_of
        >>> q = azoth.ureg.Quantity
        >>> mix, ig = mixture_of(["propane"])
        >>> r = vu_flash_single_comp(mix, ig, q(1.0e6, "Pa"),
        ...                          q(0.001059848052516, "m**3/mol"),
        ...                          q(-7797.318485008, "J/mol"))
        >>> round(r.beta, 6)
        0.5
    """
    spec = _models_gen.model(MODEL_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    if len(mixture) != 1:
        raise InvalidInputError(
            "components",
            f"a pure component's saturation state needs one component and this mixture has "
            f"{len(mixture)}. A mixture's two-phase split is `eos.vu_flash`",
        )
    p_si = input_to_si(spec, "P", P)
    apply_checks(checks.on_input, {"P": p_si}.get, warnings)

    algorithm = spec["algorithm"]
    temperature = _saturation_temperature(mixture, p_si, algorithm)
    liquid, vapour = _saturated(mixture, ideal_gas, temperature, p_si)
    (u_liq, v_liq), (u_vap, v_vap) = liquid, vapour
    u_si = U.to("J/mol").magnitude
    span = u_vap - u_liq
    if span <= 0.0 or u_si < u_liq or u_si > u_vap:
        raise OutOfRangeError(
            "U",
            u_si,
            f"the saturated liquid's internal energy at {temperature:.4f} K and {p_si:.6e} Pa "
            f"is {u_liq:.6e} J/mol and the saturated vapour's is {u_vap:.6e}, so a feed of "
            f"{u_si:.6e} J/mol is outside the two-phase region and has no split",
        )

    beta = (u_si - u_liq) / span
    implied = (1.0 - beta) * v_liq + beta * v_vap
    given = V.to("m**3/mol").magnitude
    if abs(implied - given) > 1.0e-6 * abs(implied):
        warnings.append(
            Warning(
                code=WarningCode.OUT_OF_VALID_RANGE,
                message=(
                    f"the volume asked for, {given:.9e} m**3/mol, is not the volume this "
                    f"split implies, {implied:.9e}. The pressure and the internal energy fix "
                    f"the volume for a pure component - the two phases are saturated - so the "
                    f"volume is reported from them rather than solved for, and the caller's "
                    f"is inconsistent with the state they asked about"
                ),
            )
        )

    return VuFlashSingleCompResult(
        T=from_si(temperature, "K"),
        beta=beta,
        V=from_si(implied, "m**3/mol"),
        phase=Phase.TWO_PHASE,
        warnings=tuple(warnings),
    )
