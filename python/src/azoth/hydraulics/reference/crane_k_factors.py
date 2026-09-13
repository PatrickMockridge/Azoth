"""``hydraulics.crane_k_factors`` - fitting losses by the equivalent-length
method.

```text
K = f_t * sum(n_ld for each fitting)
```

Spec: ``specs/calcs/hydraulics/crane_k_factors.yaml``

# This calc cannot validate its own inputs

The *method* is standard: a fitting's resistance coefficient is its equivalent
length ratio times the friction factor. The **coefficients** come from
``data/fittings/crane_k_factors.csv``, where every row is currently an estimated
dummy value - a placeholder of plausible magnitude, not from Crane TP-410 or any
other standard.

The consequence is unusual and worth stating plainly: because the coefficients
are placeholders, there is no correct value for this calc to be checked against,
and **no test in this repository can detect a wrong coefficient**. The tests
validate the arithmetic and the registry lookup. A pressure drop computed from
this data can be wrong by a factor of two and still look entirely reasonable, so
the placeholder status is recorded in the data file and in this calc's spec
rather than announced by a warning on the result.
"""

from __future__ import annotations

from collections.abc import Sequence

from azoth._registry_gen import spec as _spec_for
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import KComponent, KFactorsResult
from azoth.core.warnings import Warning
from azoth.hydraulics.reference.fittings import find_fitting

CALC_ID = "hydraulics.crane_k_factors"


def crane_k_factors(fittings: Sequence[str], f_t: float) -> KFactorsResult:
    """Total resistance coefficient for a list of fittings.

    Args:
        fittings: fitting ids, resolved against the registry. An unknown id
            raises rather than being skipped.
        f_t: the Darcy friction factor used as the basis. Crane specifies ``f_T``,
            the fully turbulent friction factor at the nominal fitting size;
            using the actual friction factor at the flow Reynolds number instead
            is a common and slightly more accurate variant, and is what the
            ``azoth pipe`` CLI does. Both are defensible and they give
            different answers, which is why this is an explicit input rather than
            something the function guesses.

    Raises:
        UnknownFittingError: for an id not in the registry.
        OutOfRangeError: if ``f_t`` is not positive, or if the list is empty - an
            empty list would silently return zero loss, which reads as "no
            fittings" when the caller may have meant to supply some.

    Example:
        >>> r = crane_k_factors(["90_elbow", "gate_valve_open"], 0.018)
        >>> round(r.k_total, 12)
        0.684
        >>> len(r.components)
        2
    """
    spec = _spec_for(CALC_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    apply_checks(checks.on_input, lambda name: f_t if name == "f_t" else None, warnings)

    count = float(len(fittings))
    apply_checks(checks.derived, lambda name: count if name == "n_fittings" else None, warnings)

    components: list[KComponent] = []
    for fitting_id in fittings:
        row = find_fitting(fitting_id)
        components.append(KComponent(fitting_id=row.id, n_ld=row.n_ld, k=f_t * row.n_ld))

    return KFactorsResult(
        k_total=sum(c.k for c in components),
        f_t=f_t,
        components=tuple(components),
        warnings=tuple(warnings),
    )
