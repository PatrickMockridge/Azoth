"""One card, two readers, and the numbers they resolve it to.

A keycard is read by `python/src/azoth/keycard.py` and by `azoth_eos::card`. The
comparison is over the overlay, the part both readers resolve rather than carry:
`components` and `kij`, name by name and pair by pair. The card is the shipped
template, for the reason `test_data_agreement.py` reads the real tables: a comparison
over a fixture written for the occasion agrees with itself.
"""

from __future__ import annotations

import tomllib
from pathlib import Path
from typing import Any

import pytest

from azoth import _core, keycard
from azoth.core.errors import InvalidInputError, KeycardError
from azoth.eos import components

pytestmark = pytest.mark.requires_rust

REPO_ROOT = Path(__file__).resolve().parents[2]
TEMPLATE = REPO_ROOT / "keycard.example.toml"


def rust_card(text: str) -> Any:
    """The card as Rust reads it, as the overlay it resolves to."""
    return _core.card_overlay(text)


def python_card(text: str, path: Path | None = None) -> keycard.Keycard:
    """The card as Python reads it, from the same text."""
    return keycard.use(tomllib.loads(text), path=path)


def test_both_readers_take_the_shipped_template() -> None:
    """The template is a card *both* implementations accept.

    The guard that makes the comparisons below readable: a comparison that failed
    because one reader refused the document would say nothing about the values, and
    the failure would look like a disagreement about a number.
    """
    text = TEMPLATE.read_text(encoding="utf-8")

    assert python_card(text).component("methane") is not None
    assert _core.overlay_component_rows(rust_card(text)), "the Rust reader resolved nothing"


def test_the_two_readers_resolve_the_same_components() -> None:
    """Every parameter of every name the card states, compared value by value.

    Read through the databank on both sides - `components.entry` on one,
    `overlay_entry_row` on the other - so what is compared is the *effective* record a
    calculation would use, not the card's fragment. A card naming only `omega` keeping
    the shipped `Tc` and `Pc` is the rule this checks in both implementations at once.
    """
    text = TEMPLATE.read_text(encoding="utf-8")
    card = python_card(text, TEMPLATE)
    overlay = rust_card(text)

    compared = 0
    for name in sorted(card.components):
        mine = components.entry(name, card=card)
        theirs = _core.overlay_entry_row(name, overlay)
        assert theirs.name == name
        assert theirs.tc_k == pytest.approx(mine.Tc.to("K").magnitude), name
        assert theirs.pc_pa == pytest.approx(mine.Pc.to("Pa").magnitude), name
        assert theirs.acentric_factor == pytest.approx(mine.omega), name
        compared += 1

    assert compared, "the card names no components, so nothing was compared"
    assert compared == len(_core.overlay_component_rows(overlay)), (
        "the two readers disagree about how many substances the card names"
    )


def test_the_two_readers_resolve_the_same_pairs() -> None:
    """Every interaction pair, compared as the value a mixing rule would read."""
    text = TEMPLATE.read_text(encoding="utf-8")
    card = python_card(text, TEMPLATE)

    pairs = _core.overlay_kij_rows(rust_card(text))
    assert pairs, "the card states no pairs, so nothing was compared"

    for first, second, theirs in pairs:
        mine = components.kij_for((first, second), card=card)[(0, 1)]
        assert theirs == pytest.approx(mine), f"{first}/{second}"


def test_the_comparison_notices_a_value_that_moved() -> None:
    """The comparison above is not comparing two sides of one reader.

    The card's `Tc` is changed by one kelvin in the *text*, and the Rust answer for the
    changed card must differ from the Python answer for the unchanged one. Without this,
    a test that had accidentally read the same side twice - or a card with no values in
    it at all - would pass every assertion in this file.
    """
    text = TEMPLATE.read_text(encoding="utf-8")
    moved = text.replace("value = 190.0", "value = 191.0", 1)
    assert moved != text, "the template's methane Tc is no longer 190.0 kelvin"

    mine = components.entry("methane", card=python_card(text, TEMPLATE))
    theirs = _core.overlay_entry_row("methane", rust_card(moved))

    assert theirs.tc_k != pytest.approx(mine.Tc.to("K").magnitude)


def test_the_two_readers_refuse_the_same_unit() -> None:
    """A unit neither vocabulary carries is refused by both.

    The refusal is the same *rule* implemented twice, so the two are held to it together:
    a card one reader took and the other refused would be a card that works until the
    language changes.
    """
    text = 'schema_version = 2\n[components.methane.Tc]\nvalue = 1.0\nunit = "kelvin"\n'

    with pytest.raises(KeycardError, match="not in the vocabulary"):
        python_card(text)
    with pytest.raises(InvalidInputError, match="kelvin"):
        rust_card(text)


def test_the_two_readers_agree_about_a_self_pair() -> None:
    """A substance paired with itself is refused on both sides, and named as the same
    defect: it would silently rescale that substance's attraction."""
    text = (
        'schema_version = 2\n[[kij]]\ncomponent_a = "methane"\n'
        'component_b = "methane"\nvalue = 0.1\n'
    )

    with pytest.raises(KeycardError, match="does not interact with itself"):
        python_card(text)
    with pytest.raises(InvalidInputError, match="does not interact with itself"):
        rust_card(text)
