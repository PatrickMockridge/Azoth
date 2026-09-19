"""``eos.tp_flash_saft`` - the SAFT-VR-Mie flash, the pure-Python reference.

The second, independent expression of the physics the Rust ``azoth_eos::tp_flash_saft``
computes.

**It is a different algorithm from every other flash here, and that is the model.** NeqSim
dispatches a ``SystemSAFTVRMie`` to ``TPflashSAFT`` by ``instanceof``, and that class reaches
its K-values by successive substitution - Wilson K-values to seed, Rachford-Rice inside for
the vapour fraction, and each phase's fugacity coefficients from a separate single-phase solve
at that trial composition - rather than by the Newton scheme the generic ``TPflash`` uses.

**And NeqSim's own dispatch cannot split anything.** ``TPflashSAFT.run()`` reads
``getz()`` before its own ``init(0)``, and on a fresh system that is zero for every
component, so the Rachford-Rice solves a zero feed and the vapour fraction walks to zero. The
numbers this is checked against are the class's answer on a system that *was* initialised,
from ``validation/neqsim/SaftVrMieFlashProbe.java``; the finding is written up at
``~/Desktop/neqsim-tpflashsaft-reads-a-zero-feed.md``.

The layers are :mod:`azoth.eos.reference.saft_vr_mie_phase`'s - the same module the phase
model uses - because a flash's only new arithmetic is the loop.
"""

from __future__ import annotations

import math
from collections.abc import Sequence

from azoth.core.errors import InvalidInputError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import Phase, TpFlashSaftResult
from azoth.core.units import Q, input_to_si
from azoth.core.warnings import Warning
from azoth.eos.components import entry
from azoth.eos.reference.saft_vr_mie_phase import (
    _Component,
    _ln_phi,
    _molar_volume,
)

MODEL_ID = "eos.tp_flash_saft"

#: NeqSim's ``MAX_SS_ITER``.
_MAX_SS_ITER = 50

#: The relative change in one K-value the loop stops at. **NeqSim's ``K_TOL`` is ``1e-6``
#: and this is ``1e-5``**, because the K-values carry this model's ``eta`` differences and
#: their own noise floor is about ``2e-6``: the threshold decides which branch the flash
#: takes, so one below that floor is a different answer rather than a stricter test.
_K_TOLERANCE = 1.0e-5

#: NeqSim's ``BETA_MIN``: outside ``(BETA_MIN, 1 - BETA_MIN)`` the flash reports one phase.
_BETA_MIN = 1.0e-10


def wilson_k(
    pc: Sequence[float], omega: Sequence[float], tc: Sequence[float], t: float, p: float
) -> list[float]:
    """``exp(ln(Pc/P) + 5.373 (1 + omega)(1 - Tc/T))``, which is ``Component.init``'s.

    A ratio of pressures, so it is unit-free as long as the two agree - the databank's
    pascals and NeqSim's bar give the same seed.
    """
    return [
        math.exp(math.log(pc[i] / p) + 5.373 * (1.0 + omega[i]) * (1.0 - tc[i] / t))
        for i in range(len(pc))
    ]


def _rachford_rice(z: Sequence[float], k: Sequence[float]) -> float:
    """The physical root of the Rachford-Rice equation, by bisection on its bracket.

    The bracket is the interval free of poles - ``1 + beta (K_i - 1) > 0`` for every ``i``
    - which holds the same root NeqSim's bisection on ``[0, 1]`` reaches whenever one
    exists. **The bracket is wider than ``[0, 1]`` for an ordinary seed** (at
    ``K = [5.58, 0.0138]`` it is ``(-0.218, 1.014)``), so the *root* is what is checked
    against the flash's bound and not the bracket.
    """
    lo = -math.inf
    hi = math.inf
    for value in k:
        if value > 1.0:
            lo = max(lo, 1.0 / (1.0 - value))
        elif value < 1.0:
            hi = min(hi, 1.0 / (1.0 - value))
        else:
            raise InvalidInputError(
                "k", "a K-value of exactly one puts a pole at infinity and makes the sum degenerate"
            )
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        total = sum(zi * (ki - 1.0) / (1.0 + mid * (ki - 1.0)) for zi, ki in zip(z, k, strict=True))
        if total > 0.0:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def _compositions(
    z: Sequence[float], k: Sequence[float], beta: float
) -> tuple[list[float], list[float]]:
    """The two phases' compositions at a vapour fraction."""
    x = [zi / (1.0 + beta * (ki - 1.0)) for zi, ki in zip(z, k, strict=True)]
    y = [ki * xi for ki, xi in zip(k, x, strict=True)]
    return x, y


def tp_flash_saft(
    components: list[str],
    T: Q,
    P: Q,
    z: list[float],
) -> TpFlashSaftResult:
    """The SAFT-VR-Mie flash at a temperature and a pressure.

    Args:
        components: the substance names, resolved against the databank with their
            SAFT-VR-Mie set and the cubic constants the Wilson seed is built from.
        T: absolute temperature.
        P: absolute pressure.
        z: the feed's mole fractions, checked rather than renormalised.

    Raises:
        InvalidInputError: if a component has no SAFT-VR-Mie set, or ``z`` is not a
            composition of the right length.
        OutOfRangeError: if ``T`` or ``P`` is not positive, or a trial phase has no root.

    See :func:`azoth.eos.tp_flash_saft`.
    """
    from azoth import _models_gen

    spec = _models_gen.model(MODEL_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    t_si = input_to_si(spec, "T", T)
    p_si = input_to_si(spec, "P", P)
    apply_checks(checks.on_input, {"T": t_si, "P": p_si}.get, warnings)

    if len(components) != len(z):
        raise InvalidInputError(
            "z", f"{len(components)} components but {len(z)} fractions; the two must match"
        )
    total = sum(z)
    if abs(total - 1.0) > 1.0e-9:
        raise InvalidInputError(
            "z",
            f"the mole fractions sum to {total}, not to one. Renormalising them here would "
            f"make a composition error invisible in every number downstream, so it is "
            f"refused instead",
        )

    resolved = [_Component.of(name) for name in components]
    records = [entry(name) for name in components]
    seed = wilson_k(
        [r.Pc.to("pascal").magnitude for r in records],
        [r.omega for r in records],
        [r.Tc.to("kelvin").magnitude for r in records],
        t_si,
        p_si,
    )

    def coefficients(comp: Sequence[float], side: str) -> tuple[list[float], float]:
        v, z_factor, _ = _molar_volume(resolved, comp, t_si, p_si, side)
        return _ln_phi(resolved, comp, t_si, v), z_factor

    k = list(seed)
    beta = 0.5
    iterations = 0
    residual = math.inf
    converged = False
    liquid: list[float] = []
    vapour: list[float] = []
    liquid_coefficients: list[float] = []
    vapour_coefficients: list[float] = []
    z_liquid = z_vapour = 0.0

    for step in range(_MAX_SS_ITER):
        iterations = step + 1
        beta = _rachford_rice(z, k)
        if not (_BETA_MIN <= beta <= 1.0 - _BETA_MIN):
            break
        x, y = _compositions(z, k, beta)

        liquid_coefficients, liquid_z = coefficients(x, "liquid")
        vapour_coefficients, vapour_z = coefficients(y, "vapour")

        largest = 0.0
        for i in range(len(k)):
            updated = math.exp(liquid_coefficients[i] - vapour_coefficients[i])
            relative = abs(updated - k[i]) / max(abs(k[i]), 1.0e-10)
            largest = max(largest, relative)
            k[i] = updated
        residual = largest
        liquid, vapour = x, y
        z_liquid, z_vapour = liquid_z, vapour_z
        if largest < _K_TOLERANCE:
            converged = True
            break

    if converged:
        return TpFlashSaftResult(
            beta=beta,
            x=tuple(liquid),
            y=tuple(vapour),
            k=tuple(k),
            ln_phi_liquid=tuple(liquid_coefficients),
            ln_phi_vapour=tuple(vapour_coefficients),
            z_liquid=z_liquid,
            z_vapour=z_vapour,
            phase=Phase.TWO_PHASE,
            iterations=iterations,
            residual=residual,
            warnings=tuple(warnings),
        )

    # **One phase, and which one is a Gibbs comparison.** The ideal part is identical in
    # both - same temperature, pressure and composition - so `sum_i x_i ln phi_i` decides it.
    gas_coefficients, gas_z = coefficients(z, "vapour")
    liquid_coefficients, liquid_z = coefficients(z, "liquid")
    gas_gibbs = sum(zi * p for zi, p in zip(z, gas_coefficients, strict=True))
    liquid_gibbs = sum(zi * p for zi, p in zip(z, liquid_coefficients, strict=True))
    return TpFlashSaftResult(
        beta=None,
        x=tuple(z),
        y=tuple(z),
        k=tuple(k),
        ln_phi_liquid=tuple(liquid_coefficients),
        ln_phi_vapour=(),
        z_liquid=liquid_z,
        z_vapour=gas_z,
        phase=Phase.ALL_VAPOUR if gas_gibbs <= liquid_gibbs else Phase.ALL_LIQUID,
        iterations=iterations,
        residual=residual,
        warnings=tuple(warnings),
    )


__all__ = ["tp_flash_saft"]
