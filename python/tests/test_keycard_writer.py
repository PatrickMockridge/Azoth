"""A card goes out the way it came in, and reaches a flash through the file.

`azoth.keycard.to_toml` and `dump` are the writer half of a format whose two *readers* are held to
each other by `test_card_agreement`. A writer is the third thing that knows the format, so what it
has to be held to is not itself: every test here writes a card and reads it back through **both**
readers, and the last one carries the card all the way to a flash.

**What the writer loses, and why that is not a defect to fix.** `Keycard` keeps the parameters, the
compiled sections and the three class entries. It does not keep citations, comments, the author's
section order or a model's `alpha_parameters` - so the text is what the card *resolves to* rather
than what its author typed. The tests assert the value, never the bytes of the input.
"""

from __future__ import annotations

import dataclasses
import tomllib
from pathlib import Path

import pytest

from azoth import _core, characterization, eos, keycard
from azoth.core.units import ureg

pytestmark = pytest.mark.requires_rust

REPO_ROOT = Path(__file__).resolve().parents[2]
EXAMPLE = REPO_ROOT / "keycard.example.toml"
BASELINE = REPO_ROOT / "databank" / "keycard.toml"

q = ureg.Quantity

#: One cut, the same one `test_characterization_front_door` works.
CUT_MOLAR_MASS = 0.5
CUT_DENSITY = 880.0


def written(card: keycard.Keycard) -> keycard.Keycard:
    """The card as it comes back through the writer."""
    return keycard.use(tomllib.loads(keycard.to_toml(card)))


def test_the_example_card_round_trips_as_a_value() -> None:
    """The template's own card, out and back, is the same card.

    `path` is left out of the comparison deliberately: it is where the card was *read from* rather
    than anything the card says, and a card written to a string has no such place.
    """
    card = keycard.load(EXAMPLE)
    again = written(card)
    assert dataclasses.replace(again, path=card.path) == card
    assert again.path is None, "a card read from a string names no file"


def test_writing_twice_gives_the_same_bytes() -> None:
    """**The order is the writer's and not the reader's.** A card built from a mapping in another
    order writes the same text, so a diff of two cards says what changed rather than how they were
    assembled, and writing is a fixed point.
    """
    card = keycard.load(EXAMPLE)
    once = keycard.to_toml(card)
    assert keycard.to_toml(written(card)) == once


def test_the_shipped_baseline_survives_the_writer() -> None:
    """**The generator's own output through the writer, and then through both readers.**

    `databank/keycard.toml` is emitted from the CSVs by `tools/gen_keycard.py` and its `--check`
    gate is a byte comparison of that emission, so nothing has ever read it back. This is the
    stronger claim: the 348 components it states survive a round trip as values, and the Rust
    reader resolves the written text to the same rows as the Python one.
    """
    card = keycard.load(BASELINE)
    text = keycard.to_toml(card)
    assert dataclasses.replace(written(card), path=card.path) == card

    overlay = _core.card_overlay(text)
    rows = _core.overlay_component_rows(overlay)
    assert len(rows) == len(card.components), "the written card resolved to a different row count"
    # And the Python reader resolves the same text to the same names.
    assert sorted(row.name for row in rows) == sorted(written(card).components)


def test_a_characterisation_becomes_a_card_and_reaches_a_flash(tmp_path: Path) -> None:
    """**The point of the writer, and the one test that cannot be satisfied by agreeing with
    itself.**

    A cut is a molar mass and a density. `card_from_cuts` runs `tbp_cut_properties` over it, the
    card carries the answer, `dump` writes it, and the written file is read back and flashed. Two
    independent things are asserted about that flash:

    * it is the same as **the card held in memory**, so the file route is not a second route; and
    * it **differs from the same cut flashed with `alpha="srk"`**, so the exponent the card states
      was read rather than fallen back on. A writer that dropped `alpha_param_1` would make the
      written card agree with Soave, and no writer can satisfy both.
    """
    card = characterization.card_from_cuts(
        ["c19"], [CUT_MOLAR_MASS * q("kg/mol")], [CUT_DENSITY * q("kg/m**3")]
    )
    stated = card.components["c19"].parameters
    assert set(stated) == {
        "Tc",
        "Pc",
        "omega",
        "molar_mass",
        "normal_liquid_density",
        "alpha_param_1",
    }

    path = tmp_path / "c19_keycard.toml"
    keycard.dump(card, path)
    from_file = keycard.load(path)

    T, P, z = q(400.0, "K"), q(2_000_000.0, "Pa"), [1.0]
    in_memory = eos.pt_flash(
        eos.from_names(["c19"], card=card, eos="srk", alpha="srk_fitted"), T, P, z
    )
    through_the_file = eos.pt_flash(
        eos.from_names(["c19"], card=from_file, eos="srk", alpha="srk_fitted"), T, P, z
    )
    soave = eos.pt_flash(eos.from_names(["c19"], card=from_file, eos="srk", alpha="srk"), T, P, z)
    assert through_the_file.z_vapour == pytest.approx(in_memory.z_vapour, rel=1e-12)
    assert through_the_file.z_liquid == pytest.approx(in_memory.z_liquid, rel=1e-12)
    assert through_the_file.vapour_fraction == pytest.approx(in_memory.vapour_fraction, rel=1e-12)
    assert through_the_file.z_vapour != pytest.approx(soave.z_vapour, rel=1e-9), (
        "a fitted exponent and Soave's give the same Z through the file, so the written card's "
        "value is not reaching the alpha function"
    )
    path.unlink()
