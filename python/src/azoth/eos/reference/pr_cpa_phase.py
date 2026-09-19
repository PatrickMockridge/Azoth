"""The Peng-Robinson CPA phase state, the pure-Python reference.

The second, independent expression of the physics the Rust `azoth_eos::pr_cpa_phase`
computes: the cubic's attraction and covolume replaced by each component's fitted
`aCPA_PR`/`bCPA_PR`, mixed with the `cpakij_PR` column, and the Wertheim association
contribution added to the residual Helmholtz energy.

**The first model to read the PR family's CPA columns**, and the reason the two families are
separate selections rather than one converted into the other: water's `kappa_AB` is 0.0692
for SRK against 0.046473789 for PR, and its fitted covolume is 1.4515 against 1.456360879 -
so a model that read the wrong set would be a different fluid.

**There is no NeqSim oracle for this model, and that is a finding rather than an omission.**
Against the pinned 3.20.0 jar, `SystemPrCPA` builds `ComponentSrkCPA` components that carry
their sites, and its phase never sums them: `PhaseSrkCPA` does that in its init path and
`PhasePrCPA` has the field and the setter and no block that sets it, so water reports four
sites and the phase's own total is zero. Every association term is then computed over
nothing and the flash is a Peng-Robinson run wearing the name.
`validation/neqsim/PrCpaFlash.java` prints both counts for exactly that reason. The finding
is recorded with the tranche's other two; this model's case carries its own numbers and says
so.
"""

from __future__ import annotations

from azoth.core.result import PrCpaPhaseResult
from azoth.core.units import Q, ureg
from azoth.eos.reference._cpa_phase import cpa_phase

MODEL_ID = "eos.pr_cpa_phase"


def pr_cpa_phase(
    components: list[str],
    T: Q,
    P: Q,
    z: list[float],
    compressed_phase: str,
) -> PrCpaPhaseResult:
    """One PR-CPA phase's state at a temperature, pressure and composition.

    See :func:`azoth.eos.pr_cpa_phase`; the arithmetic is `_cpa_phase`'s and what is here
    is the family.
    """
    state = cpa_phase(MODEL_ID, "pr", components, T, P, z, compressed_phase)
    return PrCpaPhaseResult(
        z_factor=state.z_factor,
        ln_phi=state.ln_phi,
        h_res=ureg.Quantity(state.h_res, "J/mol"),
        s_res=ureg.Quantity(state.s_res, "J/(mol*K)"),
        warnings=state.warnings,
    )


__all__ = ["pr_cpa_phase"]
