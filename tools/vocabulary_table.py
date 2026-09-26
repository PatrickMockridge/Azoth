#!/usr/bin/env python3
"""The vocabulary table, as the generators that need a dimension read it.

`specs/vocabulary/vocabulary.toml` maps each canonical unit to the dimension it carries, and two
generators need that mapping for the same reason: a spec states an input's **unit** and the
implementation takes a **typed quantity**, so the dimension is the only thing the two can be
compared by.

**One reader, because two would drift.** `gen_model_inputs.py` derives each model input's
dimension this way and `gen_registry.py` derives each calculation input's, and the assertions
that hold them to the specs compare against this mapping — so a second copy would be a second
answer to the question the whole correspondence rests on.
"""

from __future__ import annotations

import sys
import tomllib
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
VOCABULARY = ROOT / "specs" / "vocabulary" / "vocabulary.toml"


def load_vocabulary() -> dict[str, str]:
    """Each canonical unit's dimension id, in the table's own order."""
    table = tomllib.loads(VOCABULARY.read_text(encoding="utf-8"))
    units = {unit["id"]: unit["dimension"] for unit in table["units"]}
    if not units:
        sys.exit(f"vocabulary_table: no units under {VOCABULARY}")
    return units


def dimension_of(unit: str | None, vocabulary: dict[str, str], where: str) -> str | None:
    """The dimension an input's declared unit carries, or `None` where there is no unit.

    **A unit the vocabulary does not carry is refused rather than defaulted.** `spec_lint`
    already holds every declared unit to the schema's enum, so this cannot happen - and a
    default here would be a dimension absent for a reason nobody can see, which is the state
    this mapping exists to remove.
    """
    if unit is None:
        return None
    if unit not in vocabulary:
        sys.exit(f"vocabulary_table: {where} declares unit {unit!r}, which is not in {VOCABULARY}")
    return vocabulary[unit]
