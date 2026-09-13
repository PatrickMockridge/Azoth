"""The fittings registry: equivalent-length coefficients read from
``data/fittings/crane_k_factors.csv``.

The Rust side embeds the same file with ``include_str!``, so both languages read
one file rather than two copies that are supposed to match.

Reading the same bytes is not the same as being checked against each other, and
the parsed coefficients are **not** compared row by row. This side cannot reach
the Rust parser - the compiled extension exposes the calcs but not this registry -
so there is nothing here to compare against, and the gap is recorded rather than
described away.

**Every coefficient in that file is currently an estimated dummy value, not
engineering data.** See the file's header. :attr:`Fitting.is_estimated` reports
it, and ``crane_k_factors`` turns it into a warning on every affected result so a
computed pressure drop cannot be mistaken for a design-grade one.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from enum import StrEnum
from functools import lru_cache

from azoth._data import find
from azoth.core.errors import InvalidInputError, UnknownFittingError

#: Repo-relative path. The Rust crate embeds the same file.
REGISTRY_PATH = "data/fittings/crane_k_factors.csv"

#: Column order as declared in the CSV header.
COLUMNS = (
    "fitting_id",
    "family",
    "name",
    "n_ld",
    "f_t_basis",
    "citation",
    "verify_status",
    "source_ref",
    "source_locator",
)


def _optional(raw: str) -> str | None:
    """An empty CSV field means absent, not an empty string."""
    stripped = raw.strip()
    return stripped or None


class VerifyStatus(StrEnum):
    """Provenance of a coefficient."""

    #: Placeholder for software testing. Not from any source. Not engineering
    #: data.
    ESTIMATED_DUMMY = "estimated_dummy"
    #: Read from a secondary public reference, not checked against the primary
    #: standard.
    UNVERIFIED = "unverified"
    #: Read from the primary standard by a named engineer.
    VERIFIED = "verified"


@dataclass(frozen=True, slots=True)
class Fitting:
    """One row of the registry."""

    #: Stable identifier used by callers and specs.
    id: str
    #: Coarse category: ``bend``, ``valve``, and so on.
    family: str
    #: Human-readable name.
    name: str
    #: Equivalent length ratio, ``L_eq / D``, for fully turbulent flow.
    n_ld: float
    #: What the coefficient is multiplied by, conventionally ``f_t``.
    f_t_basis: str
    #: Where the value came from, or a statement that it came from nowhere.
    citation: str
    #: How far the value can be trusted.
    status: VerifyStatus
    #: The document the value was read from, in a fetchable form.
    #:
    #: ``arweave:<txid>`` is preferred: an Arweave transaction ID is the hash of
    #: its content, so the document is immutable, independently timestamped, and
    #: fetchable byte-for-byte by anyone. That makes a single number's provenance
    #: auditable rather than a matter of trusting whoever typed it.
    source_ref: str | None
    #: Where inside that document to look, e.g. "Table 2, 90 deg elbow".
    source_locator: str | None

    @property
    def is_estimated(self) -> bool:
        """True when this coefficient is a placeholder rather than a measurement."""
        return self.status is VerifyStatus.ESTIMATED_DUMMY


def _read_rows(lines: list[str]) -> list[Fitting]:
    reader = csv.DictReader(lines)
    if reader.fieldnames is None:  # pragma: no cover - guarded by the test suite
        raise InvalidInputError("fittings", f"{REGISTRY_PATH} has no header row")
    missing = [c for c in COLUMNS if c not in reader.fieldnames]
    if missing:
        raise InvalidInputError("fittings", f"{REGISTRY_PATH} is missing columns {missing}")

    rows: list[Fitting] = []
    for raw in reader:
        try:
            n_ld = float(raw["n_ld"])
        except ValueError as exc:
            raise InvalidInputError("n_ld", f"not a number in row {raw['fitting_id']!r}") from exc
        try:
            status = VerifyStatus(raw["verify_status"])
        except ValueError as exc:
            raise InvalidInputError(
                "verify_status",
                f"unknown value {raw['verify_status']!r} in row {raw['fitting_id']!r}",
            ) from exc
        rows.append(
            Fitting(
                id=raw["fitting_id"],
                family=raw["family"],
                name=raw["name"],
                n_ld=n_ld,
                f_t_basis=raw["f_t_basis"],
                citation=raw["citation"],
                status=status,
                source_ref=_optional(raw["source_ref"]),
                source_locator=_optional(raw["source_locator"]),
            )
        )
    return rows


@lru_cache(maxsize=1)
def registry() -> tuple[Fitting, ...]:
    """The registry, parsed once.

    The header block is stripped before the CSV parser sees it. It is
    documentation and provenance rather than decoration, so the reader has to
    tolerate it rather than the file having to be machine-only.
    """
    path = find(REGISTRY_PATH)
    lines = [
        line
        for line in path.read_text(encoding="utf-8").splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ]
    if len(lines) < 2:  # pragma: no cover - guarded by the test suite
        raise InvalidInputError("fittings", f"{REGISTRY_PATH} has no data rows")
    return tuple(_read_rows(lines))


def find_fitting(fitting_id: str) -> Fitting:
    """Look up one fitting by id.

    Raises:
        UnknownFittingError: if the id is not in the registry. An error rather
            than a skip: silently treating an unknown fitting as zero loss would
            under-report pressure drop, which is the dangerous direction to be
            wrong in.
    """
    for row in registry():
        if row.id == fitting_id:
            return row
    raise UnknownFittingError(fitting_id)
