"""The pure-Python reference implementation.

Always importable, with or without the compiled Rust extension. Two reasons it is
a public module rather than a private one:

* It is the *reference* the Rust implementation is checked against. A reference
  that can only be reached through the thing it validates is not much of a check.
* It is the fallback when the extension is not built, and a fallback that is
  harder to reach than the fast path is a fallback nobody tests.

The bodies here are written to mirror the Rust line for line, so a reviewer can
read the two side by side against the published equation. Where they diverge, it
is a bug in one of them - and the cross-implementation tests are what catch it.
"""

from __future__ import annotations

from azoth.hydraulics.reference.choked_flow_area import choked_flow_area
from azoth.hydraulics.reference.control_valve_cv import CV_TO_SI, control_valve_cv
from azoth.hydraulics.reference.crane_k_factors import crane_k_factors
from azoth.hydraulics.reference.darcy_weisbach import add_fitting_loss, darcy_weisbach
from azoth.hydraulics.reference.fittings import (
    REGISTRY_PATH,
    Fitting,
    VerifyStatus,
    find_fitting,
    registry,
)
from azoth.hydraulics.reference.friction_factor_colebrook import (
    friction_factor_colebrook,
    fully_rough_limit,
)
from azoth.hydraulics.reference.friction_factor_haaland import (
    friction_factor_haaland,
)
from azoth.hydraulics.reference.friction_factor_swamee_jain import (
    friction_factor_swamee_jain,
)
from azoth.hydraulics.reference.orifice_flow import orifice_flow
from azoth.hydraulics.reference.pump_power import STANDARD_GRAVITY_M_S2, pump_power
from azoth.hydraulics.reference.reynolds_number import reynolds_number
from azoth.hydraulics.reference.solver import (
    Convergence,
    SolverOutcome,
    fixed_point,
    require_converged,
)

__all__ = [
    "CV_TO_SI",
    "REGISTRY_PATH",
    "STANDARD_GRAVITY_M_S2",
    "Convergence",
    "Fitting",
    "SolverOutcome",
    "VerifyStatus",
    "add_fitting_loss",
    "choked_flow_area",
    "control_valve_cv",
    "crane_k_factors",
    "darcy_weisbach",
    "find_fitting",
    "fixed_point",
    "friction_factor_colebrook",
    "friction_factor_haaland",
    "friction_factor_swamee_jain",
    "fully_rough_limit",
    "orifice_flow",
    "pump_power",
    "registry",
    "require_converged",
    "reynolds_number",
]
