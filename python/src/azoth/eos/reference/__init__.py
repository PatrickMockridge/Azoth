"""Pure-Python reference implementations of the equations-of-state calculations.

Every calculation in azoth has two implementations, and this is the Python one. It
is not a fallback for a missing extension: it is the second opinion the Rust core
is checked against, and the test suite runs both against every case in the specs.
See ``README.md`` for why a single implementation was rejected.

The bodies here are written to mirror the Rust line for line, so a reviewer can
read the two side by side against the published equation. Where they diverge, it is
a bug in one of them - and the cross-implementation tests are what catch it.
"""

from __future__ import annotations

from azoth.eos.reference.pr_alpha_ab import OMEGA_A, OMEGA_B, pr_alpha_ab
from azoth.eos.reference.pr_departure import pr_departure
from azoth.eos.reference.pr_kappa import pr_kappa
from azoth.eos.reference.pr_mass_density import pr_mass_density
from azoth.eos.reference.pr_molar_volume import MOLAR_GAS_CONSTANT, pr_molar_volume
from azoth.eos.reference.pr_z_factor import pr_z_factor
from azoth.eos.reference.prsv_kappa import prsv_kappa
from azoth.eos.reference.pure_saturation import pure_saturation
from azoth.eos.reference.rachford_rice_binary import rachford_rice_binary
from azoth.eos.reference.vdw1f_mix_binary import vdw1f_mix_binary

__all__ = [
    "MOLAR_GAS_CONSTANT",
    "OMEGA_A",
    "OMEGA_B",
    "pr_alpha_ab",
    "pr_departure",
    "pr_kappa",
    "pr_mass_density",
    "pr_molar_volume",
    "pr_z_factor",
    "prsv_kappa",
    "pure_saturation",
    "rachford_rice_binary",
    "vdw1f_mix_binary",
]
