#!/usr/bin/env python3
"""Check that an *installed* azoth can find its own data.

Run this with a wheel installed, from a directory outside the repository:

    cd /tmp && /path/to/clean-venv/bin/python /path/to/repo/tools/check_wheel_data.py

# Why this is separate from the test suite

The suite runs against the source tree, where ``azoth/_data.py`` walks up and finds
``data/`` at the repository root. That path works no matter what the wheel contains,
so every test in the suite would pass on a wheel that shipped no data at all. The one
thing a source-tree run cannot exercise is the thing that was broken: resolving the
tables from an *installed* distribution.

# What it asserts

1. ``azoth`` is the installed copy, not the source tree. Otherwise this proves nothing
   about the artifact - and it would still pass.
2. A calculation whose data comes from the fittings registry works, and reports the
   estimated-data warning its placeholder coefficients require.
3. A fluid property lookup works.
4. Neither resolved data file lives inside the repository. If one does, the source
   tree was reachable and assertion 1 was passed by accident.

Exit status is non-zero if any of those fails.
"""

from __future__ import annotations

import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent


def fail(message: str) -> None:
    print(f"check_wheel_data: {message}", file=sys.stderr)
    raise SystemExit(1)


def main() -> int:
    import azoth
    from azoth._data import find
    from azoth.core.warnings import WarningCode
    from azoth.hydraulics import crane_k_factors
    from azoth.properties import provider_for

    installed = Path(azoth.__file__).resolve()
    print(f"  azoth imported from {installed}")

    if installed.is_relative_to(REPO_ROOT):
        fail(
            f"azoth resolved to the source tree ({installed}) rather than an installed "
            f"wheel. Run this from outside {REPO_ROOT} with the wheel's interpreter, "
            f"or it proves nothing about the artifact."
        )

    # 1. The fittings registry, which is read from a CSV the wheel must carry.
    result = crane_k_factors(["90_elbow", "gate_valve_open"], 0.018)
    if result.k_total <= 0:
        fail(f"crane_k_factors returned a non-positive coefficient: {result.k_total}")
    if not result.has_warning(WarningCode.ESTIMATED_DATA):
        fail(
            "the placeholder coefficients did not produce an ESTIMATED_DATA warning. "
            "Either the data file is no longer the placeholder registry, or the "
            "warning has stopped firing - both are worth stopping for."
        )
    print(f"  crane_k_factors: k_total={result.k_total} (with the expected warning)")

    # 2. A fluid table, which is a separate file and a separate parser.
    water = provider_for("water")
    density = water.density(azoth.ureg.Quantity(20.0, "degC"))
    if not 990.0 < density.magnitude < 1005.0:
        fail(f"water density at 20 C came back as {density}, which is not water")
    print(f"  water at 20 C: {density:.2f}")

    # 3. Neither file may have come from the repository.
    for relative in ("data/fittings/crane_k_factors.csv", "data/fluids/water.csv"):
        resolved = find(relative).resolve()
        if resolved.is_relative_to(REPO_ROOT):
            fail(
                f"{relative} resolved inside the repository ({resolved}), so the "
                f"source tree was reachable after all and this run did not exercise "
                f"the packaged path."
            )
        print(f"  {relative} -> {resolved}")

    print("check_wheel_data: OK")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
