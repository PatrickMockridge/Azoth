#!/usr/bin/env python3
"""Name the first layer a model and its NeqSim capture disagree at.

    python tools/neqsim_layer_diff.py <case-id> [--tolerance T]
    python tools/neqsim_layer_diff.py --all

A case compares a *total*, and a total that is 4% out says nothing about where. This
compares the intermediates instead: `azoth.eos.layers` dumps the reference kernel's
own layers under the capture's own key names, so the two sides meet without a mapping
table, and the first key whose relative difference exceeds the case's tolerance is
printed with both numbers beside it.

**Every case in `gen_neqsim_cases.CASES` has one of these**, and a divergence is meant
to be read in this order: the layers first, then the case.

The capture is read through that generator's readers rather than through a second
parser, for the reason the two exist at all - a second parser is a second set of
answers about what a probe printed.

`--all` walks every case and is what `python/tests/test_layer_diff.py` runs, so a
divergence fails the build naming its layer. Nothing here needs a JDK or the NeqSim
jar: the captures are committed and the layers are azoth's.
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

    if args.all:
        from azoth.eos import layers as azoth_layers

        dumped = [case.id for case in gen.CASES if case.calc in azoth_layers.DUMPERS]
        skipped = [case.id for case in gen.CASES if case.calc not in azoth_layers.DUMPERS]
        if skipped:
            print(f"neqsim_layer_diff: no dumper for {', '.join(skipped)}")
        ids = dumped
    else:
        ids = [args.case]
    quiet = args.all
    ok = True
    for case_id in ids:
        ok = report(case_id, args.tolerance, quiet) and ok
    if args.all:
        print(f"neqsim_layer_diff: {len(ids)} case(s), {'no divergence' if ok else 'FAILED'}")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
