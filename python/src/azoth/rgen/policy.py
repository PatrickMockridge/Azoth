"""Which result fields two implementations may legitimately disagree about.

These are the two sets ``python/tests/_helpers.py`` uses when it compares the kernels.
They live here as well because a wheel cannot import the tests, and the bundle's agreement
check must apply the same rule or it will call a solver's iteration count a defect. A test
asserts the two copies are set-equal, so a change on either side fails rather than drifts.

``DIAGNOSTIC_FIELDS`` are residuals: they are evidence that a loop converged, meaningful
only down to the tolerance the loop declares, so they are compared on that absolute bound
rather than relatively. ``UNCOMPARED_FIELDS`` are path quantities - an iteration count, a
loop's last error, its verdict - whose value is where the arithmetic stopped rather than
what it computed. Both implementations' own tests pin them against the capture.
"""

from __future__ import annotations

from typing import Final

__all__ = ["DIAGNOSTIC_FIELDS", "UNCOMPARED_FIELDS"]

#: Result fields compared on the solver's declared absolute bound, not relatively.
DIAGNOSTIC_FIELDS: Final[frozenset[str]] = frozenset(
    {
        "balance_error",
        "bottom_specification_residual",
        "charge_moles",
        "charge_residual",
        "chemical_deviation",
        "convergence_residual",
        "element_balance_difference",
        "element_balance_residual",
        "element_residual",
        "final_error",
        "max_element_residual",
        "max_log_fugacity_residual",
        "max_material_balance_residual",
        "max_reaction_log_residual",
        "maximum_complementarity_violation",
        "maximum_ion_balance_residual",
        "net_charge_moles",
        "reaction_log_residual",
        "residual",
        "tm",
        "top_specification_residual",
    }
)

#: Result fields that are recorded but not judged: a path quantity is not an answer.
UNCOMPARED_FIELDS: Final[frozenset[str]] = frozenset(
    {
        "condition_number_history",
        "converged",
        "energy_residual",
        "error",
        "iterations",
        "mass_residual",
        "outer_iterations",
        "passes",
        "temperature_residual",
        "total_inner_iterations",
        "total_iterations",
    }
)
