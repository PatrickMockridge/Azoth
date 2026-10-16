#!/usr/bin/env python3
"""Is the nightly still running, and still green?

    python tools/check_nightly.py --check

**Why this exists.** The whole scoped-push design rests on one thing the tree cannot see:
that the nightly actually runs. GitHub disables a scheduled workflow after roughly sixty days
of repository inactivity, and it does so without telling anyone - no failure, no email, just a
gate that stops existing. `release.yml` now runs the full surface at a tag, so a release
cannot ship on a tree nothing verified. But between releases, a schedule that stopped would
mean every push was verified against nothing, and nothing would say so.

Exit codes follow the tree's convention: **1 is a fact about the nightly, 2 is could-not-
decide**. They are separate because "the API was unreachable" and "the nightly is 90 hours
old" are not the same finding, and a check that conflated them would either cry wolf or go
quiet exactly when it mattered.

`evaluate` is a pure function of the runs, so the three findings can be tested from fixtures
without a network - see `python/tests/test_check_nightly.py`.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import urllib.error
import urllib.request
from datetime import UTC, datetime

#: The workflow this watches. Named rather than passed, because watching a different one is a
#: different check and would need a different set of thresholds.
WORKFLOW = "nightly.yml"

#: How old the newest run may be before the schedule counts as stopped. The cycle is 24 hours,
#: so 36 leaves room for a late queue and a slow run without hiding a missed day.
MAX_AGE_HOURS = 36.0

#: How long the nightly may go without a green run. Two cycles: one red is a finding to
#: triage, two is a gate that has stopped being looked at.
MAX_RED_HOURS = 48.0


def _when(run: dict) -> datetime:
    """A run's timestamp, as an aware datetime."""
    return datetime.fromisoformat(run["created_at"].replace("Z", "+00:00"))


def evaluate(
    runs: list[dict], now: datetime, max_age_hours: float, max_red_hours: float
) -> tuple[int, str]:
    """`(exit code, message)` for a list of workflow runs, newest-agnostic.

    Three findings, each with its own message so a red is readable without the API payload in
    hand: never run, stopped, and running-but-not-green.
    """
    if not runs:
        return 1, (
            f"the {WORKFLOW} workflow has never run. Until it does, a green push is a claim "
            "about a diff and nothing anywhere is checking the tree."
        )

    newest = max(runs, key=_when)
    age_hours = (now - _when(newest)).total_seconds() / 3600.0
    if age_hours > max_age_hours:
        return 1, (
            f"the newest {WORKFLOW} run is {age_hours:.0f} hours old, over the {max_age_hours:.0f}"
            "-hour ceiling. GitHub disables a schedule after about sixty days without activity, "
            "and it does not say so. Re-enable the workflow, then dispatch it once."
        )

    green = [run for run in runs if run.get("conclusion") == "success"]
    if not green:
        return 1, (
            f"the {WORKFLOW} has not succeeded in the {len(runs)} run(s) on record, the newest "
            f"being {newest.get('conclusion')!r}. A red nightly is a finding to triage, not "
            "scenery - it is the only thing running the tests no diff claimed."
        )

    newest_green = max(green, key=_when)
    red_hours = (now - _when(newest_green)).total_seconds() / 3600.0
    if red_hours > max_red_hours:
        return 1, (
            f"the {WORKFLOW} has not been green for {red_hours:.0f} hours, over the "
            f"{max_red_hours:.0f}-hour ceiling. Triage it before the next tranche; a red one "
            "left standing makes every later red unreadable."
        )

    return 0, (
        f"{WORKFLOW} last succeeded {red_hours:.0f} hours ago, and its newest run is "
        f"{age_hours:.0f} hours old."
    )


def fetch(repo: str, token: str, workflow: str, limit: int = 30) -> list[dict]:
    """The workflow's recent runs, newest first, from the GitHub API. Raises on any failure."""
    url = f"https://api.github.com/repos/{repo}/actions/workflows/{workflow}/runs?per_page={limit}"
    request = urllib.request.Request(
        url,
        headers={
            "Accept": "application/vnd.github+json",
            "Authorization": f"Bearer {token}",
            "X-GitHub-Api-Version": "2022-11-28",
            "User-Agent": "azoth-check-nightly",
        },
    )
    with urllib.request.urlopen(request, timeout=30) as response:
        payload = json.load(response)
    runs = payload.get("workflow_runs")
    if not isinstance(runs, list):
        raise ValueError(f"the API returned no `workflow_runs` list: {sorted(payload)}")
    return runs


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--check", action="store_true", help="the only mode there is")
    parser.add_argument("--repo", default=os.environ.get("GITHUB_REPOSITORY", ""))
    parser.add_argument("--workflow", default=WORKFLOW)
    parser.add_argument("--max-age-hours", type=float, default=MAX_AGE_HOURS)
    parser.add_argument("--max-red-hours", type=float, default=MAX_RED_HOURS)
    args = parser.parse_args()

    if not args.repo:
        print("check_nightly: could not decide - no --repo and no $GITHUB_REPOSITORY")
        return 2
    token = os.environ.get("GITHUB_TOKEN", "")
    if not token:
        print("check_nightly: could not decide - no $GITHUB_TOKEN to read the API with")
        return 2
    try:
        runs = fetch(args.repo, token, args.workflow)
    except (urllib.error.URLError, ValueError, OSError) as exc:
        # **A check that cannot read the API has not found the nightly healthy.** Returning 0
        # here would make an outage look like a pass, which is the one thing this whole
        # machinery is built not to do.
        print(f"check_nightly: could not decide - {exc}")
        return 2

    code, message = evaluate(runs, datetime.now(UTC), args.max_age_hours, args.max_red_hours)
    print(f"check_nightly: {message}")
    return code


if __name__ == "__main__":
    sys.exit(main())
