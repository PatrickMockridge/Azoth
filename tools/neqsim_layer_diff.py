#!/usr/bin/env python3
"""Name the first layer a model and its NeqSim capture disagree at.

    python tools/neqsim_layer_diff.py <case-id> [--tolerance T]
    python tools/neqsim_layer_diff.py --all

A case compares a *total*, and a total that is 4% out says nothing about where. This
compares the intermediates instead: `azoth.eos.layers` and `azoth.process.layers` dump the
reference kernel's own layers under the capture's own key names, so the two sides meet
without a mapping table, and the first key that is not where it is declared to be is printed
with both numbers beside it.

**Two tiers, two kinds of case id.** The `eos` cases are `gen_neqsim_cases.CASES`, named by
their own id. A process case is named `model::case` - `process.pump::the_work_divides_by_the
_efficiency` - which is exactly the id `python/tests/test_process_layer_diff.py` parametrises
by, so a failing test's name pastes straight in here.

The capture is read through that generator's readers rather than through a second
parser, for the reason the two exist at all - a second parser is a second set of
answers about what a probe printed. The process captures are read by
`azoth.process.layers` for the same reason.

`--all` walks every case of both tiers and is what the two layer tests run, so a divergence
fails the build naming its layer. Nothing here needs a JDK or the NeqSim jar: the captures
are committed and the layers are azoth's.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent.parent
TOOLS = ROOT / "tools"
sys.path.insert(0, str(TOOLS))

import gen_neqsim_cases as gen  # noqa: E402  (the path above is what makes it importable)


def case_inputs(case_id: str) -> tuple[dict[str, Any], float]:
    """A case's inputs and tolerance, from the file the case *is*.

    Read from `validation/eos/` rather than from the generator's table so that the
    tolerance a divergence is judged against is the one the case itself compares with -
    the table and the file are kept equal by `--check`, and this takes the file.
    """
    path = ROOT / "validation" / "eos" / f"{case_id}.json"
    if not path.is_file():
        raise SystemExit(f"neqsim_layer_diff: no {path.relative_to(ROOT)}")
    case: dict[str, Any] = json.loads(path.read_text(encoding="utf-8"))
    return case["inputs"], float(case["tolerance"])


def comparable(case: Any, inputs: dict[str, Any]) -> list[tuple[str, float, float]]:
    """Every layer both sides have, as `(key, neqsim, azoth)` in the dump's order.

    The order is the dumper's and not the capture's, because "the first divergence" is
    a statement about the build order of the model - the same reason the probes print
    their rows in it.
    """
    from azoth.eos import layers as azoth_layers

    text = (gen.CAPTURES / case.capture).read_text(encoding="utf-8")
    record = gen.SHAPES[case.shape](text)[case.state]
    dumped = azoth_layers.rows(case.calc, inputs)
    out = []
    for key, mine in dumped:
        if key not in record:
            continue
        try:
            theirs = float(record[key].strip().strip("[]"))
        except ValueError:
            # A capture key that is not a number - an enum name, say. Nothing to
            # difference, and nothing a dumper should have offered.
            continue
        out.append((key, theirs, mine))
    return out


def relative(mine: float, theirs: float) -> float:
    """The relative difference, with an absolute one where the oracle is zero.

    A layer that is exactly zero has no relative difference to take, and the layers
    where that happens are not incidental - an ideal-mixing `ln gamma` or a vanishing
    departure is a number a port is most likely to get wrong in a way a ratio cannot
    see.
    """
    if theirs == 0.0:
        return abs(mine)
    return abs(mine / theirs - 1.0)


def first_divergence(
    rows: list[tuple[str, float, float]], tolerance: float
) -> tuple[str, float, float, float] | None:
    """The first key past the tolerance, or `None`."""
    for key, theirs, mine in rows:
        gap = relative(mine, theirs)
        if gap > tolerance:
            return key, theirs, mine, gap
    return None


def is_process(case_id: str) -> bool:
    """Whether an id names a `process.*` model's case rather than an `eos` one."""
    return case_id.startswith("process.") or case_id.startswith("process::")


def process_case(case_id: str) -> Any:
    """The `LayerCase` a `model::case` id names."""
    from azoth.process import layers as process_layers

    model, _, case = case_id.partition("::")
    for entry in process_layers.LAYER_CASES:
        if entry.model == model and entry.case == case:
            return entry
    raise SystemExit(
        f"neqsim_layer_diff: no process case {case_id!r}; "
        f"`--all` carries {', '.join(f'{c.model}::{c.case}' for c in process_layers.LAYER_CASES)}"
    )


def report_process(case_id: str, tolerance: float | None, quiet: bool) -> bool:
    """One process case. `True` when every layer is where it is declared to be.

    **The rules are `azoth.process.layers.compare`'s**, the same object the pytest gate
    applies, so this cannot report a divergence the build does not fail on - or the reverse.
    """
    from azoth.process import layers as process_layers

    entry = process_case(case_id)
    inputs, case_tolerance = process_layers.case_inputs(entry.model, entry.case)
    bound = case_tolerance if tolerance is None else tolerance

    diffs = process_layers.compare(entry, inputs, tolerance=tolerance)
    if not diffs:
        raise SystemExit(
            f"neqsim_layer_diff: {case_id} has no layer in common with {entry.capture}, "
            f"so there is nothing to compare and nothing to report"
        )
    missed = [diff for diff in diffs if not diff.ok]
    if not missed:
        if not quiet:
            print(f"neqsim_layer_diff: {case_id}: {len(diffs)} layer(s) agree to {bound:g}")
        return True

    first = missed[0]
    print(
        f"neqsim_layer_diff: {case_id}: first divergence is `{first.key}`"
        f"{'' if first.by == 'case' else f' (a {first.by})'}",
        file=sys.stderr,
    )
    print(f"    neqsim  {first.capture!r}", file=sys.stderr)
    print(f"    azoth   {first.azoth!r}", file=sys.stderr)
    print(f"    {first.describe().lstrip('.')}", file=sys.stderr)
    for later in missed[1:]:
        print(f"    also    {later.describe().lstrip('.')}", file=sys.stderr)
    return False


def report(case_id: str, tolerance: float | None, quiet: bool) -> bool:
    """One case. `True` when nothing diverged."""
    from azoth.eos import layers as azoth_layers

    case = next((c for c in gen.CASES if c.id == case_id), None)
    if case is None:
        raise SystemExit(f"neqsim_layer_diff: {case_id} is not a recorded case")
    if case.calc not in azoth_layers.DUMPERS:
        # Named outright this is a caller's mistake; `--all` filters first.
        raise SystemExit(f"neqsim_layer_diff: {case_id} has no layer dumper")
    inputs, case_tolerance = case_inputs(case_id)
    bound = case_tolerance if tolerance is None else tolerance

    rows = comparable(case, inputs)
    if not rows:
        raise SystemExit(
            f"neqsim_layer_diff: {case_id} has no layer in common with "
            f"{case.capture}, so there is nothing to compare and nothing to report"
        )
    found = first_divergence(rows, bound)
    if found is None:
        if not quiet:
            print(f"neqsim_layer_diff: {case_id}: {len(rows)} layer(s) agree to {bound:g}")
        return True

    key, theirs, mine, gap = found
    print(
        f"neqsim_layer_diff: {case_id}: first divergence is `{key}` at "
        f"{gap:.3e} relative, over {bound:g}",
        file=sys.stderr,
    )
    print(f"    neqsim  {theirs!r}", file=sys.stderr)
    print(f"    azoth   {mine!r}", file=sys.stderr)
    for later in rows[rows.index((key, theirs, mine)) + 1 :]:
        gap = relative(later[2], later[1])
        if gap > bound:
            print(
                f"    also    `{later[0]}` at {gap:.3e}: {later[2]!r} against {later[1]!r}",
                file=sys.stderr,
            )
    return False


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("case", nargs="?", help="a recorded case id")
    parser.add_argument("--all", action="store_true", help="every recorded case")
    parser.add_argument("--tolerance", type=float, default=None, help="override the case's")
    args = parser.parse_args(argv)

    if args.all == (args.case is not None):
        parser.error("give a case id or --all, not both and not neither")

    from azoth.process import layers as process_layers

    if args.all:
        from azoth.eos import layers as azoth_layers

        dumped = [case.id for case in gen.CASES if case.calc in azoth_layers.DUMPERS]
        skipped = [case.id for case in gen.CASES if case.calc not in azoth_layers.DUMPERS]
        if skipped:
            print(f"neqsim_layer_diff: no dumper for {', '.join(skipped)}")
        if process_layers.NO_INTERIOR:
            # Named on separate lines, because the joined list is the one string here that
            # grows with the tier and an f-string cannot be wrapped by the formatter.
            interiorless = ",\n    ".join(sorted(process_layers.NO_INTERIOR))
            print(f"neqsim_layer_diff: no interior for {interiorless}")
        ids = dumped + [f"{entry.model}::{entry.case}" for entry in process_layers.LAYER_CASES]
    else:
        ids = [args.case]
    quiet = args.all
    ok = True
    for case_id in ids:
        if is_process(case_id):
            ok = report_process(case_id, args.tolerance, quiet) and ok
        else:
            ok = report(case_id, args.tolerance, quiet) and ok
    if args.all:
        print(f"neqsim_layer_diff: {len(ids)} case(s), {'no divergence' if ok else 'FAILED'}")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
