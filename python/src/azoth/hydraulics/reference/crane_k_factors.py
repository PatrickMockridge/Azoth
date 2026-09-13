"""``hydraulics.crane_k_factors`` - fitting losses by the equivalent-length
method.

```text
K = f_t * sum(n_ld for each fitting)
```

Spec: ``specs/calcs/hydraulics/crane_k_factors.yaml``

# This calc cannot validate its own inputs, and says so

The *method* is standard: a fitting's resistance coefficient is its equivalent
length ratio times the friction factor. The **coefficients** come from
``data/fittings/crane_k_factors.csv``, where every row is currently an estimated
dummy value - a placeholder of plausible magnitude, not from Crane TP-410 or any
other standard.

The consequence is unusual and worth stating plainly: because the coefficients
are placeholders, there is no correct value for this calc to be checked against,
and **no test in this repository can detect a wrong coefficient**. The tests
validate the arithmetic and the registry lookup. That is why every result built
from estimated rows carries an ``ESTIMATED_DATA`` warning: a pressure drop
computed from this data can be wrong by a factor of two and still look entirely
reasonable.
"""

from __future__ import annotations

from collections.abc import Sequence

from azoth._registry_gen import spec as _spec_for
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import KComponent, KFactorsResult
from azoth.core.warnings import Warning, WarningCode
from azoth.hydraulics.reference.fittings import VerifyStatus, find_fitting

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
        >>> r.has_warning(azoth.core.WarningCode.ESTIMATED_DATA)
        True
    """
    spec = _spec_for(CALC_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    apply_checks(checks.on_input, lambda name: f_t if name == "f_t" else None, warnings)

    count = float(len(fittings))
    apply_checks(checks.derived, lambda name: count if name == "n_fittings" else None, warnings)

    components: list[KComponent] = []
    estimated: list[str] = []
    unverified: list[str] = []
    for fitting_id in fittings:
        row = find_fitting(fitting_id)
        if row.status is VerifyStatus.ESTIMATED_DUMMY:
            estimated.append(row.id)
        elif row.status is VerifyStatus.UNVERIFIED:
            unverified.append(row.id)
        components.append(KComponent(fitting_id=row.id, n_ld=row.n_ld, k=f_t * row.n_ld))

    # Provenance warnings, driven by the data rather than hardcoded: promote a row
    # and the warning for it stops firing with no code change.
    #
    # Two levels, because there are two different things to say. A placeholder is
    # not engineering data at all; a cited-but-unconfirmed value is a real
    # published figure that nobody has checked against an authoritative copy.
    # Both deserve a warning, and collapsing them into one would either overstate
    # the first or understate the second.
    if estimated:
        warnings.append(
            Warning(
                code=WarningCode.ESTIMATED_DATA,
                message=(
                    f"{len(estimated)} of {len(fittings)} fitting(s) use ESTIMATED DUMMY "
                    f"coefficients that are not engineering data "
                    f"({', '.join(estimated)}). This resistance coefficient is a "
                    f"placeholder and must not be used to size equipment. Populate "
                    f"data/fittings/crane_k_factors.csv from the primary standard and "
                    f"set verify_status=verified."
                ),
            )
        )
    if unverified:
        warnings.append(
            Warning(
                code=WarningCode.UNVERIFIED_SOURCE,
                message=(
                    f"{len(unverified)} of {len(fittings)} fitting(s) use coefficients "
                    f"that are cited but NOT CONFIRMED by a named verifier against an "
                    f"authoritative copy of the source ({', '.join(unverified)}). They "
                    f"may be correct; nobody has checked. Treat this resistance "
                    f"coefficient as provisional and confirm it before sizing equipment."
                ),
            )
        )

    return KFactorsResult(
        k_total=sum(c.k for c in components),
        f_t=f_t,
        components=tuple(components),
        warnings=tuple(warnings),
    )
