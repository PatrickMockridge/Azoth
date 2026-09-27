#!/usr/bin/env python3
"""Capture the documents the front end's mirrors are held to.

`ui/test/fixtures/` holds three documents the library wrote, and every shape
`ui/src/wire/types.ts` declares is checked against them. **They are a *copy* of the library's
output, and a copy goes stale**: the envelope gained a `provenance` record, a spec moved a
parameter, and the fixtures went on describing a library that no longer existed until
`ui/test/currency.test.ts` failed - which it did, twice, on other people's commits.

That test is the *gate*: it asks the module a build just produced for the same three documents
and compares, so it runs in the `ui` job with the module that job already builds and costs
nothing. This tool is the *recapture*, because a gate can say "the fixture is stale" and cannot
say what to write in its place. They are the same three commands in two languages.

    python tools/gen_ui_fixtures.py            # write
    python tools/gen_ui_fixtures.py --check    # fail if a fixture is not what the CLI writes

**The commands are the ones `ui/test/fixtures.test.ts` names**, and they are the whole pin: a
fixture captured some other way asserts something about a document nobody can reproduce.
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
FIXTURES = ROOT / "ui" / "test" / "fixtures"
DEMO = "specs/flowsheets/demo.toml"

#: The position one fixture's edit sets. A **position** and nothing else, so the assertions
#: that read it can say one gesture is in the document *and* in the graph, and the run it
#: reports is the shipped demo's own.
SET_POSITION = '{"command":"set_position","node":"instance:sep1","x":1234,"y":56}'

#: And the removal the refused fixture is. The heater goes, two ports are left underfed, and
#: the document is one the checker refuses - which is a state a canvas draws rather than an
#: error, which is why this capture's exit code is not zero.
REMOVE_HEATER = '{"command":"remove_instance","id":"hx1"}'

#: The fixture, and the CLI arguments that write it.
CASES: tuple[tuple[str, list[str]], ...] = (
    ("catalogue.json", ["forms", "--tools"]),
    (
        "envelope.json",
        ["edit", "--flowsheet", DEMO, "--command", SET_POSITION, "--run", "--json"],
    ),
    ("broken.json", ["edit", "--flowsheet", DEMO, "--command", REMOVE_HEATER, "--json"]),
)

#: The exit codes a capture may answer with and still have written its document.
#:
#: **`2` is the document being refused, not a failure of the command.** `azoth edit` answers
#: `2` when the edit left a document the checker refuses - `ok: false` with the diagnostics,
#: which is `broken.json`'s whole content - and treating that as an error would make this tool
#: unable to capture the one fixture whose point is that a document can be broken.
ACCEPTED = (0, 2)


def capture(arguments: list[str]) -> str:
    """Run the CLI, and answer with what it wrote to stdout."""
    proc = subprocess.run(
        ["cargo", "run", "-q", "-p", "azoth-cli", "--", *arguments],
        cwd=ROOT,
        capture_output=True,
        text=True,
        check=False,
    )
    if proc.returncode not in ACCEPTED:
        sys.exit(
            f"gen_ui_fixtures: `azoth {' '.join(arguments)}` exited {proc.returncode}, which is "
            f"neither a document nor a refused one:\n{proc.stderr}"
        )
    if not proc.stdout.strip():
        # A capture that wrote nothing would make every fixture empty and every mirror agree
        # with it, which is the vacuous pass this tool and the test beside it exist to refuse.
        sys.exit(f"gen_ui_fixtures: `azoth {' '.join(arguments)}` wrote nothing")
    return proc.stdout


def main() -> None:
    check = "--check" in sys.argv
    stale: list[str] = []
    for name, arguments in CASES:
        written = capture(arguments)
        path = FIXTURES / name
        if check:
            if not path.exists() or path.read_text(encoding="utf-8") != written:
                stale.append(name)
            continue
        path.write_text(written, encoding="utf-8")
        print(f"gen_ui_fixtures: wrote ui/test/fixtures/{name}")
    if stale:
        print(
            f"gen_ui_fixtures: {', '.join(stale)} is out of date against the library this tree "
            f"builds. Run `python tools/gen_ui_fixtures.py`.",
            file=sys.stderr,
        )
        raise SystemExit(1)
    if check:
        print(f"gen_ui_fixtures: {len(CASES)} fixture(s) up to date")


if __name__ == "__main__":
    main()
