"""``process.event_tree`` - the event-tree enumeration.

As NeqSim's ``EventTreeAnalyzer`` computes it. Spec: ``specs/models/process/event_tree.toml``

The Python twin of ``crates/azoth-process/src/models/event_tree.rs``: the full ``2^n``
enumeration in bit order, the path spelling and the products.
"""

from __future__ import annotations

import math
from typing import Any

from azoth.core.errors import InvalidInputError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import EventTreeOutcome, EventTreeResult
from azoth.core.units import input_records_to_si
from azoth.core.warnings import Warning

MODEL_ID = "process.event_tree"


def _spec() -> dict[str, Any]:
    """The generated spec, imported here so the module has no import-time cycle."""
    from azoth._models_gen import model

    return model(MODEL_ID)


def event_tree(
    initiating_frequency_per_year: float,
    branches: list[dict[str, object]],
) -> EventTreeResult:
    """The full ``2^n`` event tree in bit order.

    Args:
        initiating_frequency_per_year: the initiating event's frequency in 1/year.
        branches: the pivotal events in branch order, each ``{name, probability_of_yes}``;
            an empty list is a tree with no branches.

    Returns:
        One row per leaf with its path, probability and frequency.

    Raises:
        InvalidInputError: for a negative initiating frequency, a branch probability
            outside ``[0, 1]`` or a NaN, or more than the algorithm's branch cap.
    """
    spec = _spec()
    checks = checks_for(spec)
    warnings: list[Warning] = []

    frequency = float(initiating_frequency_per_year)
    apply_checks(checks.on_input, {"initiating_frequency_per_year": frequency}.get, warnings)
    if not math.isfinite(frequency):
        raise InvalidInputError(
            "initiating_frequency_per_year",
            "a non-finite frequency is not a number this calculation publishes",
        )

    rows = input_records_to_si(spec, "branches", branches)
    cap = int(spec["algorithm"]["max_iterations"])
    if len(rows) > cap:
        raise InvalidInputError(
            "branches",
            f"the tree enumerates 2^n leaves, and {len(rows)} branches exceeds the {cap} "
            "the procedure carries",
        )
    probabilities: list[float] = []
    for index, row in enumerate(rows):
        probability = float(row["probability_of_yes"])
        if math.isnan(probability) or not 0.0 <= probability <= 1.0:
            raise InvalidInputError(
                f"branches[{index}].probability_of_yes",
                f"a branch probability is in [0, 1], and {probability} is not",
            )
        probabilities.append(probability)

    outcomes = []
    for leaf in range(1 << len(rows)):
        probability = 1.0
        parts: list[str] = []
        for bit, row in enumerate(rows):
            yes = (leaf >> bit) & 1 == 1
            probability *= probabilities[bit] if yes else 1.0 - probabilities[bit]
            parts.append(f"{row['name']}={'Y' if yes else 'N'}")
        outcomes.append(
            EventTreeOutcome(
                path="; ".join(parts),
                probability=probability,
                frequency_per_year=probability * frequency,
            )
        )

    return EventTreeResult(outcomes=tuple(outcomes), warnings=tuple(warnings))
