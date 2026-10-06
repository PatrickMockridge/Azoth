"""Every claim this tree makes is of a kind `docs/src/architecture/claims.md` names.

Two vocabularies answer one question - *what carries this claim?* - from two sides: the calculus
pages write a `*Status:` out of three words, and `lean/guards.toml` writes an `owner` out of five.
`claims.md` is where the two are related, and this is what holds it: a seventh status word, or a
sixth owner, fails here rather than being read as emphasis.

**The failure this exists for is a word that looks like a status and is not one.** `**probably**`,
`**mostly**`, `**believed**` - each reads as a kind of assurance and names no artifact, which is the
shape this repository refuses everywhere else: a claim that says something is covered without saying
what covers it.
"""

from __future__ import annotations

import re
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
CLAIMS = REPO_ROOT / "docs" / "src" / "architecture" / "claims.md"
GUARDS = REPO_ROOT / "lean" / "guards.toml"

#: The trees whose `*Status:*` lines are read. The same pair `check_doc_claims.swept_pages`
#: sweeps, and for the same reason: those two are where a claim about the tree is stated.
TREES = ("calculus", "architecture")

#: A status line, and the words it holds in bold. Read only from `*Status:` because that is what
#: makes a word a status: `index.md`'s "**one of them is proved**" is prose about a claim, not a
#: claim's own status.
# **The block, not the line.** A page is free to wrap an italic span and `rho.md` does, so a
# line-anchored pattern reads only its first line and its `**proved**` and `**specified**`
# were compared by nothing at all.
_STATUS = re.compile(r"^\*Status: (.*?)\*\s*$", re.M | re.S)
_BOLD = re.compile(r"\*\*([a-z][a-z-]*)\*\*")

#: The kinds `claims.md` names, and nothing else is one.
KINDS = frozenset({"proved", "checked", "data", "unheld", "specified", "characterised"})

#: The owners `lean/guards.toml` may name. A closed vocabulary the generator and the checker
#: already enforce; held again here so the two lists cannot drift apart silently.
OWNERS = frozenset({"lean", "check", "construction", "port", "unformalised"})


def _status_words() -> dict[str, list[str]]:
    """**The word each status block opens with**, by the page it is on.

    The *first* bold word and not every one: a status block states its status first and then
    argues, so the bold words after it are the prose's emphasis - `session.md` says
    `**signature**` and `view.md` says `**two**` - and reading those as statuses would report
    two pages for writing ordinary English.
    """
    out: dict[str, list[str]] = {}
    for tree in TREES:
        for page in sorted((REPO_ROOT / "docs" / "src" / tree).glob("*.md")):
            words = []
            for block in _STATUS.findall(page.read_text(encoding="utf-8")):
                bold = _BOLD.findall(block)
                if bold:
                    words.append(bold[0])
            if words:
                out[page.relative_to(REPO_ROOT).as_posix()] = words
    return out


def _owners() -> set[str]:
    """Every `owner` a guard or a class entry names."""
    text = GUARDS.read_text(encoding="utf-8")
    return set(re.findall(r'^owner = "([a-z]+)"', text, re.M))


def test_no_page_invents_a_status_word() -> None:
    """**A word that names no artifact is the failure, not the word's spelling.**

    Every status a page writes has to be one `claims.md` names - so a claim's kind is something a
    reader can look up and a sweep can hold, rather than an assurance a writer chose.
    """
    pages = _status_words()
    assert pages, (
        "no `*Status:*` line was found in either tree - this is looking in the wrong place"
    )

    unknown = {
        page: sorted(set(words) - KINDS) for page, words in pages.items() if set(words) - KINDS
    }
    assert not unknown, (
        f"{unknown}: these pages use a status word `docs/src/architecture/claims.md` does not "
        f"name, so nothing says what carries the claim. Either it is one of "
        f"{sorted(KINDS)}, or the kind it is needs naming on that page and in the table."
    )


def test_the_kinds_the_table_names_are_the_ones_in_use() -> None:
    """Both directions, because either alone is satisfiable by accident.

    A table naming a kind no page uses is a row nobody can check; a page using a kind the table
    does not name is what the test above catches. Together they say the vocabulary is closed and
    exact - which is the only thing that makes `claims.md` a rule rather than a description.
    """
    used = {word for words in _status_words().values() for word in words}
    assert used, "no status word was read"
    unused = sorted(KINDS - {"data", "unheld"} - used)
    assert not unused, (
        f"{unused} are named as kinds in `claims.md` and no page states one. `data` and `unheld` "
        f"are the two that never appear as a page's own status - they are the kinds the databank "
        f"and the guard manifest carry - so a third unused one is a row nobody can check."
    )


def test_every_guard_owner_is_in_the_vocabulary() -> None:
    """The other vocabulary, held to the same table.

    `tools/guards.py` refuses an owner outside its five, so this cannot fire while that runs - and
    it is here anyway, because the two vocabularies are what `claims.md` exists to relate and one
    of them silently widening would leave the relation wrong.
    """
    owners = _owners()
    assert owners, (
        "no `owner` was read from lean/guards.toml - the parse is looking for the wrong key"
    )
    assert owners <= OWNERS, f"{sorted(owners - OWNERS)}: not an owner this tree knows"
    assert owners, "no owner is named at all, so nothing holds anything"
