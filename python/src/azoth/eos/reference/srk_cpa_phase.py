"""The Soave-Redlich-Kwong CPA phase state, the pure-Python reference.

The second, independent expression of the physics the Rust `azoth_eos::srk_cpa_phase`
computes: the cubic's attraction and covolume replaced by each component's fitted
`aCPA_SRK`/`bCPA_SRK`, mixed with the `cpakij_SRK` column, and the Wertheim association
contribution added to the residual Helmholtz energy.

The arithmetic is the family's alone in one place - see `_cpa_phase` - and what is here is
the family: `eos="srk"`, which selects the fitted set and the interaction column this model
reads. **The fluid is resolved from names here, and again in Rust**, which is what makes the
two-kernel comparison cover the resolution as well as the arithmetic: an associating mixture
mixes with `cpakij_SRK` and a classical one with `KIJPR`, and on water/methanol those differ
by a factor of two.
"""

from __future__ import annotations

from azoth.core.result import SrkCpaPhaseResult
from azoth.core.units import Q, ureg
from azoth.eos.reference._cpa_phase import cpa_phase

MODEL_ID = "eos.srk_cpa_phase"


def srk_cpa_phase(
    components: list[str],
    T: Q,
    P: Q,
    z: list[float],
    compressed_phase: str,
) -> SrkCpaPhaseResult:
    """One SRK-CPA phase's state at a temperature, pressure and composition.

    See :func:`azoth.eos.srk_cpa_phase`; the arithmetic is `_cpa_phase`'s and what is here
    is the family.
    """
    state = cpa_phase(MODEL_ID, "srk", components, T, P, z, compressed_phase)
    return SrkCpaPhaseResult(
        z_factor=state.z_factor,
        ln_phi=state.ln_phi,
        h_res=ureg.Quantity(state.h_res, "J/mol"),
        s_res=ureg.Quantity(state.s_res, "J/(mol*K)"),
        warnings=state.warnings,
    )


__all__ = ["srk_cpa_phase"]
