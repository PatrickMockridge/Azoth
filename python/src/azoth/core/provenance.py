"""What a number rests on, carried by the number.

A result's values say what was computed; its ``warnings`` say which declared checks
could not be evaluated. Neither says *what code produced it*, and that is the gap this
closes: a reader holding a number - a person or an agent - should not have to go and
find `provenance.json`, the spec, and the two kernels to say where it came from.

# Static half and dynamic half

The static half is a fact about the *calculation* and never changes between calls: its
id, the spec and the two implementations with a hash of each, the source its equation
is attributed to, and how far the answer is externally checked. `tools/gen_registry.py`
embeds it, so it travels with an installed wheel rather than only with a checkout.

The dynamic half is a fact about *this call*: which declared checks were skipped, and
which warnings were raised. `RANGE_CHECK_SKIPPED` is the one that matters most, because
it is the difference between a bound that was tested and a bound that was never
reached - and that difference is invisible in the value.

# Why the status is derived rather than declared

A spec refuses a verification-status field, and this does not add one back. The status
is measured from the evidence the tree already holds - the tests a spec ships and the
external cases under `validation/` - so it cannot disagree with what it summarises. A
declared status is a status that can be wrong, which is the failure this whole module
exists to make impossible.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from typing import Any, Final, cast

from azoth.core import serialise
from azoth.core.warnings import Warning, WarningCode

#: The document key a result's block is written under.
PROVENANCE_KEY: Final[str] = "provenance"

#: Keys a result's JSON document carries that are not dataclass fields.
#:
#: Exported rather than repeated in the test that asserts the document's shape, so the
#: test cannot hold its own copy of a list the library owns.
FRAMEWORK_DOCUMENT_KEYS: Final[tuple[str, ...]] = (PROVENANCE_KEY,)

#: Warning codes that mean *a check did not run*, as opposed to one that ran and found
#: something. The distinction is the library's oldest rule: "checked and fine" and
#: "never checked" are different states and must not read the same.
SKIPPED_CHECK_CODES: Final[frozenset[WarningCode]] = frozenset({WarningCode.RANGE_CHECK_SKIPPED})


@dataclass(frozen=True, slots=True)
class Provenance:
    """The static half and the dynamic half of one result, together.

    Frozen, and a dataclass, for the same reason a result is: :mod:`azoth.core.serialise`
    walks ``dataclasses.fields``, so this block is written by the same one writer as
    everything else and needs no codec of its own.

    **Flat rather than nested**, and deliberately: the block travels inside every result
    document and inside an MCP envelope carrying many stream records at once, so a
    ``{path, sha256}`` object per file would roughly triple its size to carry a ``present``
    flag that has no state - a shipped kernel that is missing is a defect the generator
    refuses, not a case a reader handles.
    """

    #: The spec id, e.g. ``hydraulics.darcy_weisbach``.
    calc_id: str
    #: The spec's own name for the calculation.
    name: str
    #: The spec file the bounds and the equation come from, and its hash.
    spec_path: str
    spec_sha256: str
    #: Both implementations, each with its hash. Neither is a wrapper for the other.
    python_path: str
    python_sha256: str
    rust_path: str
    rust_sha256: str
    #: What the spec attributes its equation to, e.g. ``Crane TP-410``.
    source: str
    #: ``verified``, ``partially_verified``, ``source_needed`` or ``unverified``.
    verification: str
    #: How many external cases exist for this id, and the spec's own test tally.
    validation_cases: int
    tests_active: int
    tests_skipped: int
    #: Declared checks this call could not evaluate, by the field they are about.
    skipped_checks: tuple[str, ...]
    #: Every warning raised, by code, in the warning vocabulary's own order.
    warning_codes: tuple[str, ...]
    #: True when the call carried no warnings at all.
    clean: bool

    @classmethod
    def of(cls, static: Mapping[str, Any], warnings: Sequence[Warning]) -> Provenance:
        """Merge the generated static half with what this call actually reported.

        The dynamic half is read from ``warnings`` rather than passed separately, because
        a caller able to report "which checks were skipped" apart from the warnings that
        said so would be a second answer to a question the result has already answered.
        """
        # `field` is optional on a warning, and a skipped check that names no field is
        # one this tuple cannot carry - it is a tuple of *fields*, so an unnamed entry
        # would put a `None` in it rather than a name.
        skipped = {
            warning.field
            for warning in warnings
            if warning.code in SKIPPED_CHECK_CODES and warning.field is not None
        }
        raised = {warning.code for warning in warnings}
        return cls(
            calc_id=static["calc_id"],
            name=static["name"],
            spec_path=static["spec"]["path"],
            spec_sha256=static["spec"]["sha256"],
            python_path=static["code"][0]["path"],
            python_sha256=static["code"][0]["sha256"],
            rust_path=static["code"][1]["path"],
            rust_sha256=static["code"][1]["sha256"],
            source=static["source"],
            verification=static["verification"],
            validation_cases=static["validation_cases"],
            tests_active=static["tests_active"],
            tests_skipped=static["tests_skipped"],
            skipped_checks=tuple(sorted(skipped)),
            # Ordered by the vocabulary rather than by the call, so two runs raising the
            # same warnings in a different order produce the same block - which is what
            # lets the cross-language test compare them.
            warning_codes=tuple(code.value for code in WarningCode if code in raised),
            clean=not warnings,
        )

    def to_dict(self) -> dict[str, Any]:
        """This block as JSON-shaped data, written by the walker every result uses."""
        return cast("dict[str, Any]", serialise.to_dict(self))
