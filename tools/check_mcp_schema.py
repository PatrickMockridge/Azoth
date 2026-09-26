#!/usr/bin/env python3
"""The vendored MCP schema is the revision it says it is, and upstream still publishes it.

`azoth mcp` is a hand-rolled implementation of a published protocol, and the one thing such an
implementation cannot check about itself is whether it writes what the specification says. That is
what the vendored schema under `python/tests/validation/vendor/mcp/` is for, and
`python/tests/validation/test_mcp_conformance.py` holds every message both transports write to it.

**The pin is the URL and the digest together, and neither was checked by anything.** A file that
drifted from what was vendored would keep passing every conformance test, because those tests read
the same file. And upstream moving is invisible: the schema carries no `$id`, so nothing in it names
the revision it is, and a new revision is not a change to this tree at all.

So there are two checks here, with two different failure modes, which is why they are not one
command:

    python tools/check_mcp_schema.py              # the digests, against `NOTICE`
    python tools/check_mcp_schema.py --upstream   # and upstream, against each URL in `NOTICE`

The first belongs in CI on every push. The second does not, and the reason is the exit code: a
network call in a push gate makes every build depend on GitHub being reachable, and a fetch that
failed would then be indistinguishable from a revision that moved. An unreachable upstream is
therefore exit 2 - *the check could not decide* - and drift is exit 1, which is what lets the
scheduled run say which of the two happened.

**`NOTICE` is read rather than restated.** A digest written down a second time is a second pin, and
the two would disagree the day one of them moved.
"""

from __future__ import annotations

import argparse
import hashlib
import sys
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent

#: Where the schemas are vendored, and the file that pins them.
VENDOR = Path("python/tests/validation/vendor/mcp")
NOTICE = "NOTICE"

#: What a fetch is given before it is called unreachable.
TIMEOUT = 30

#: What a fetch can fail with, none of which is a revision having moved.
UNREACHABLE = (OSError, TimeoutError)


def sha256(path: Path) -> str:
    """The digest of a file's bytes, which is what a pin is."""
    return hashlib.sha256(path.read_bytes()).hexdigest()


def pin(notice: str) -> dict[str, dict[str, str]]:
    """Every schema `NOTICE` records, with the digest and the URL it records for each.

    Two line shapes - `<name> sha256 <digest>` and `<name> <url>` - told apart by the token after
    the name rather than by their position, so a record's lines can be split across a continuation
    without this having to know how the block is laid out.
    """
    found: dict[str, dict[str, str]] = {}
    for line in notice.splitlines():
        parts = line.split()
        # **The first line of a record carries its label and the continuations do not** - `Digest:`
        # names two files and only the first is on the line with the word - so the label is dropped
        # rather than required, which is what lets both shapes parse the same way.
        if parts and parts[0].endswith(":"):
            parts = parts[1:]
        if not parts or not parts[0].startswith("schema-"):
            continue
        name = parts[0]
        if len(parts) == 3 and parts[1] == "sha256":
            found.setdefault(name, {})["digest"] = parts[2]
        elif len(parts) == 2 and parts[1].startswith("https://"):
            found.setdefault(name, {})["url"] = parts[1]
    return found


def declaration(root: Path, records: dict[str, dict[str, str]]) -> list[str]:
    """Every way the pin is incomplete, or does not describe this tree.

    **These are not drift, and the distinction is the exit code.** A file with no digest has nothing
    to be held to and a record with no file is a digest of something that is not here: in both cases
    the check has learned nothing, which is a different report from *these bytes moved*. It is the
    same split `check_doc_claims.py` makes between a broken declaration and a page that drifted.

    **Both directions.** A vendored file `NOTICE` records nothing for is a pin with a hole in it,
    and a record naming no file is what a file deleted without its record looks like - and it would
    otherwise pass, because the loop over the files would never reach it.
    """
    found: list[str] = []
    vendor = root / VENDOR
    if not vendor.is_dir():
        return [f"{VENDOR} is not there, so there is nothing to hold to a digest"]

    files = sorted(vendor.glob("schema-*.json"))
    if not files:
        found.append(f"no `schema-*.json` under {VENDOR}, so this check decided nothing")
    names = {path.name for path in files}

    for name in sorted(names):
        record = records.get(name)
        if record is None:
            found.append(f"{VENDOR}/{name} is vendored and {NOTICE} records nothing for it")
            continue
        if "digest" not in record:
            found.append(f"{NOTICE} records no digest for {name}")
        if "url" not in record:
            found.append(
                f"{NOTICE} records no URL for {name}, so nothing can re-fetch it and a new "
                f"revision would be invisible"
            )
    for name in sorted(set(records) - names):
        found.append(f"{NOTICE} records {name}, and nothing under {VENDOR} is that file")
    return found


def drift(root: Path, records: dict[str, dict[str, str]]) -> list[str]:
    """Every vendored file whose bytes are not what `NOTICE` pins.

    A file the declaration half has already reported is skipped rather than reported twice: where
    there is no digest there is nothing to compare, and the same sentence twice is noise.
    """
    found: list[str] = []
    for path in sorted((root / VENDOR).glob("schema-*.json")):
        record = records.get(path.name)
        pinned = None if record is None else record.get("digest")
        if pinned is None:
            continue
        if (actual := sha256(path)) != pinned:
            found.append(
                f"{VENDOR}/{path.name} is not what was vendored: {NOTICE} pins {pinned}, the file "
                f"is {actual}"
            )
    return found


def drifted(records: dict[str, dict[str, str]]) -> list[str]:
    """Each URL fetched and hashed, against the digest `NOTICE` pins for it.

    # Errors
    Whatever the fetch raised. That is the caller's to report as *undecided* rather than as drift:
    a check that could not reach upstream has not learned anything about the revision.
    """
    found: list[str] = []
    for name, record in sorted(records.items()):
        url = record["url"]
        with urllib.request.urlopen(url, timeout=TIMEOUT) as response:
            body = response.read()
        actual = hashlib.sha256(body).hexdigest()
        if actual != record["digest"]:
            found.append(
                f"{url} is not what {NOTICE} pins for {name}: pinned {record['digest']}, upstream "
                f"is {actual} - a revision whose bytes moved is one somebody has to read before "
                f"the port can claim to follow it"
            )
    return found


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "--upstream",
        action="store_true",
        help="also re-fetch each URL in NOTICE and compare; needs the network",
    )
    args = parser.parse_args()

    records = pin((ROOT / NOTICE).read_text(encoding="utf-8"))

    # A pin that does not describe this tree is not a disagreement about a revision: the check has
    # nothing to compare, so it says so rather than reporting a pass it did not earn.
    broken = declaration(ROOT, records)
    if broken:
        for problem in broken:
            print(f"check_mcp_schema: {problem}", file=sys.stderr)
        raise SystemExit(2)

    stale = drift(ROOT, records)
    if stale:
        for problem in stale:
            print(f"check_mcp_schema: {problem}", file=sys.stderr)
        raise SystemExit(1)

    if not args.upstream:
        print(f"check_mcp_schema: OK ({len(records)} schema(s), each the digest NOTICE pins)")
        return

    try:
        moved = drifted(records)
    except UNREACHABLE as error:
        print(f"check_mcp_schema: upstream could not be read: {error}", file=sys.stderr)
        print("check_mcp_schema: this is not a claim that the revision moved", file=sys.stderr)
        raise SystemExit(2) from None
    if moved:
        for problem in moved:
            print(f"check_mcp_schema: {problem}", file=sys.stderr)
        raise SystemExit(1)
    print(f"check_mcp_schema: OK ({len(records)} schema(s), upstream unchanged since vendoring)")


if __name__ == "__main__":
    main()
