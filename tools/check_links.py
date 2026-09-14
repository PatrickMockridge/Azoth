#!/usr/bin/env python3
"""Check that every internal link in the documentation resolves.

# Why this exists rather than mdbook-linkcheck

`mdbook-linkcheck` 0.7.7 requires mdBook ^0.4 and this project uses 0.5, and it
is no longer maintained. `mdbook-linkcheck2` is maintained and would work, but it
costs a multi-minute Rust build in every CI run to check a few dozen relative
links across roughly ten generated pages.

This script is stdlib-only, runs in milliseconds, and checks the case that
actually breaks here: a link to a page that does not exist. mdBook does not fail
on those - it renders a link the reader discovers is dead - so without this the
docs-build job would pass while shipping broken navigation.

# What it does not do

External links are not fetched. Checking them makes a docs build depend on
whatever those servers are doing, which is a flaky test dressed up as a
thorough one; the pages here are reachable and stable, and the cost of an
occasional dead external link is lower than the cost of a build that fails for
reasons the author cannot fix.

Usage:
    python tools/check_links.py
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DOCS_SRC = ROOT / "docs" / "src"

#: Markdown outside the book. The README's links are as capable of being wrong
#: as the book's, and a dead one there is as broken as a dead one in the book.
ROOT_PAGES = (
    "README.md",
    "CODE_OF_CONDUCT.md",
    "SECURITY.md",
    "SPEC.md",
    "validation/README.md",
    "databank/README.md",
)

#: Markdown inline links and images: `[text](target)` and `![alt](target)`.
LINK = re.compile(r"!?\[[^\]]*\]\(([^)]+)\)")

#: A target with a scheme (http:, https:, mailto:) or protocol-relative (`//`)
#: is external and not fetched.
EXTERNAL = re.compile(r"^([a-zA-Z][a-zA-Z0-9+.-]*:|//)")


def strip_anchor(target: str) -> str:
    """Drop any `#fragment` - only the file part is checked."""
    return target.split("#", 1)[0]


def check_file(path: Path, errors: list[str]) -> int:
    """Check one markdown file. Returns the number of links examined."""
    checked = 0
    for lineno, line in enumerate(path.read_text(encoding="utf-8").splitlines(), start=1):
        for match in LINK.finditer(line):
            target = match.group(1).strip()
            if not target or EXTERNAL.match(target):
                continue
            checked += 1
            file_part = strip_anchor(target)
            if not file_part:
                # A bare anchor: links within the same page.
                continue
            resolved = (path.parent / file_part).resolve()
            if resolved.is_dir():
                resolved = resolved / "index.md"
            if not resolved.exists():
                errors.append(
                    f"{path.relative_to(ROOT)}:{lineno}: link target {target!r} does not exist"
                )
    return checked


def check_summary_paths(errors: list[str]) -> int:
    """Every page must appear in SUMMARY.md, or mdBook will not build it.

    A markdown file that exists but is unreferenced is silently dropped from the
    book, so a page can be generated and never appear.
    """
    summary = DOCS_SRC / "SUMMARY.md"
    if not summary.exists():  # pragma: no cover - guarded by the test suite
        errors.append("docs/src/SUMMARY.md does not exist; run tools/gen_docs.py")
        return 0

    referenced = {
        (DOCS_SRC / strip_anchor(m.group(1).strip())).resolve()
        for m in LINK.finditer(summary.read_text(encoding="utf-8"))
    }
    checked = 0
    for page in sorted(DOCS_SRC.rglob("*.md")):
        if page.name == "SUMMARY.md":
            continue
        checked += 1
        if page.resolve() not in referenced:
            errors.append(
                f"{page.relative_to(ROOT)} is not listed in SUMMARY.md, so mdBook "
                f"will not include it in the book"
            )
    return checked


def main() -> int:
    if not DOCS_SRC.exists():
        sys.exit(f"check_links: {DOCS_SRC} does not exist")

    errors: list[str] = []
    pages = sorted(p for p in DOCS_SRC.rglob("*.md") if p.name != "SUMMARY.md")
    if not pages:
        sys.exit(f"check_links: no pages found under {DOCS_SRC}")

    links = sum(check_file(page, errors) for page in pages)
    links += check_summary_paths(errors)

    # Root-level markdown, checked for broken links but not for SUMMARY
    # membership, which only applies inside the book.
    for name in ROOT_PAGES:
        page = ROOT / name
        if not page.exists():
            errors.append(f"{name} is listed as a root page but does not exist")
            continue
        links += check_file(page, errors)
        pages.append(page)

    if errors:
        for error in errors:
            print(f"  ERROR  {error}", file=sys.stderr)
        print(f"\ncheck_links: FAILED with {len(errors)} broken link(s)", file=sys.stderr)
        return 1

    print(f"check_links: OK ({links} link(s) across {len(pages)} page(s))")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
