"""``characterization.tbp_grouping`` - a phase's components binned by normal boiling point.

Spec: ``specs/models/characterization/tbp_grouping.toml``. Oracle:
``validation/neqsim/captures/tbp_grouping_probe.tsv``.

Mirrors ``crates/azoth-characterization/src/tbp_grouping.rs`` line for line.

# One name, three implementations

``groupTBPfractions`` is declared once and implemented three times. ``PlusCharacterize``'s
returns ``true`` and computes nothing. ``TBPCharacterize``'s does the binning and is reachable
only from inside its own package, because ``Characterise.TBPCharacterise`` has no accessor.
``Phase``'s is the one ported: public on a public class, though **not** on the
``PhaseInterface`` that ``SystemInterface.getPhase`` returns, so a caller needs a cast.

Nothing in NeqSim calls any of the three, and nothing in its tests either.

# The bins are boiling points, not carbon numbers

Fourteen thresholds in degrees Celsius, the lowest at 69.2, each opening a bin indexed from six.
Bins 0 to 5 exist in the returned array and no component can reach them, a component below 69.2
reaches no bin at all, and two components between one pair of thresholds sum into one bin - so
the twenty entries do not add up to one and are not meant to.

The comparison is on ``getNormalBoilingPoint("C")``, which is the stored value less 273.15. This
takes the stored kelvin and subtracts the same 273.15, so the comparison sees the same number.
"""

from __future__ import annotations

from azoth import _models_gen
from azoth.core.errors import InvalidInputError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import TbpGroupingResult
from azoth.core.units import Q, input_to_si
from azoth.core.warnings import Warning

MODEL_ID = "characterization.tbp_grouping"

#: The bins the answer always carries, of which only the last fourteen can be filled.
BIN_COUNT = 20

#: The fourteen thresholds, in degrees Celsius, paired with the bin each opens. Descending,
#: because that is the order the class tests them in.
THRESHOLDS = (
    (331.0, 19),
    (317.0, 18),
    (303.0, 17),
    (287.0, 16),
    (271.1, 15),
    (253.9, 14),
    (235.9, 13),
    (216.8, 12),
    (196.4, 11),
    (174.6, 10),
    (151.3, 9),
    (126.1, 8),
    (98.9, 7),
    (69.2, 6),
)

#: The kelvin-to-Celsius step ``getNormalBoilingPoint("C")`` makes.
ABSOLUTE_ZERO_CELSIUS = 273.15


def tbp_grouping(
    boiling_point: list[Q],
    mole_fraction: list[float],
) -> TbpGroupingResult:
    """A phase's components grouped into boiling-point bins.

    ``boiling_point`` is the **stored kelvin**, which is the number the class subtracts 273.15
    from before comparing against its thresholds in degrees Celsius. A component below the lowest
    threshold is grouped nowhere rather than into bin zero, and the twenty entries do not sum to
    one because of it.

    Raises:
        InvalidInputError: if the two vectors are not the same length.
    """
    spec = _models_gen.model(MODEL_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    temperatures = [input_to_si(spec, "boiling_point", value) for value in boiling_point]
    fractions = [float(value) for value in mole_fraction]
    # The spec declares no range on either input, which is the house rule for a vector.
    apply_checks(checks.on_input, lambda _name: None, warnings)

    if len(temperatures) != len(fractions):
        raise InvalidInputError(
            "mole_fraction",
            f"the two vectors are one component list and must be the same length: "
            f"{len(temperatures)} boiling points against {len(fractions)} mole fractions",
        )

    group_fraction = [0.0] * BIN_COUNT
    for temperature, fraction in zip(temperatures, fractions, strict=True):
        celsius = temperature - ABSOLUTE_ZERO_CELSIUS
        for threshold, bin_index in THRESHOLDS:
            if celsius >= threshold:
                group_fraction[bin_index] += fraction
                break

    apply_checks(checks.derived, lambda _name: None, warnings)

    return TbpGroupingResult(
        group_fraction=tuple(group_fraction),
        warnings=tuple(warnings),
    )
