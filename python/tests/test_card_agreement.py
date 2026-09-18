"""One card, two readers, and the numbers they resolve it to.

A keycard is read by `python/src/azoth/keycard.py` and by `azoth_eos::card`. The
comparison is over the overlay, the part both readers resolve rather than carry:
`components`, `associations` and `kij`, name by name and pair by pair. The card is the
shipped template, for the reason `test_data_agreement.py` reads the real tables: a
comparison over a fixture written for the occasion agrees with itself.
"""

from __future__ import annotations

import re
import tomllib
from pathlib import Path
from typing import Any

import pytest

from azoth import _core, keycard
from azoth.core._units_gen import CANONICAL_UNITS
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
    for name in sorted(card.names()):
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


def test_the_two_readers_resolve_the_same_association() -> None:
    """A card's association, resolved through the databank on both sides.

    The association is the parameter set a card could state last, and the two readers
    hold it in different scales: the card states SI and the shipped table carries
    NeqSim's internal one. So this compares the *resolved* record on both sides rather
    than the numbers either reader was handed - a crossing applied twice, or not at all,
    is a factor of ten thousand million and would otherwise be found by a flash.
    """
    text = TRANSPORTED_ASSOCIATION_CARD
    card = python_card(text)
    overlay = rust_card(text)

    stated = card.association_for("water")
    assert stated is not None, "the card states water's association"
    assert stated.scheme == "2B"
    assert stated.parameters["a_srk"].magnitude == pytest.approx(0.12277)

    mine = components.entry("water", card=card).association
    theirs = _core.overlay_entry_row("water", overlay)

    assert mine is not None
    assert theirs.association_scheme == mine.scheme
    assert theirs.association_sites == mine.sites
    assert theirs.association_energy == pytest.approx(mine.energy)
    assert theirs.association_a_srk == pytest.approx(mine.a_srk)
    assert theirs.association_b_srk == pytest.approx(mine.b_srk)
    assert theirs.association_m_srk == pytest.approx(mine.m_srk)
    assert theirs.association_a_pr == pytest.approx(mine.a_pr)
    assert theirs.association_b_pr == pytest.approx(mine.b_pr)
    assert theirs.association_m_pr == pytest.approx(mine.m_pr)
    # The table's own scale, so the card's SI value of 0.12277 arrives as 12277.
    assert theirs.association_a_srk == pytest.approx(12277.0)
    # Water ships `4C` with a count of four; the card states `2B`, so the count is the new
    # scheme's. Both readers, because a count kept from the other scheme is the silent
    # resolution the record's `sites` field exists to prevent.
    assert theirs.association_sites == 2 == mine.sites


#: A card stating an association, with the shipped water's numbers so both sides have a
#: table row to merge against. The scheme is *changed* from the shipped `4C` to `2B`, so a
#: reader that ignored the scheme and kept the table's would be visible.
TRANSPORTED_ASSOCIATION_CARD = (
    "schema_version = 2\n"
    '[associations.water]\nscheme = "2B"\n'
    '[associations.water.energy]\nvalue = 16655.0\nunit = "J/mol"\n'
    '[associations.water.a_srk]\nvalue = 0.12277\nunit = "Pa*m**6/mol**2"\n'
    '[associations.water.b_srk]\nvalue = 1.4515e-5\nunit = "m**3/mol"\n'
)


def test_the_two_readers_refuse_the_same_association_parameter() -> None:
    """A parameter name and a scheme name neither reader knows, refused by both.

    Held to the same rule for the reason the unit test below is: a card one reader takes
    and the other refuses works until the language changes.
    """
    for body, match in (
        ('[associations.water]\nscheme = "3B"\n', "3B"),
        (
            '[associations.water]\nscheme = "4C"\n'
            '[associations.water.epsilon]\nvalue = 1.0\nunit = "J/mol"\n',
            "epsilon",
        ),
    ):
        text = f"schema_version = 2\n{body}"

        with pytest.raises(KeycardError, match=match):
            python_card(text)
        with pytest.raises(InvalidInputError, match=match):
            rust_card(text)


def test_the_two_readers_refuse_the_same_association_unit() -> None:
    """An attraction stated as a ratio: the dimension check, on both sides."""
    text = (
        'schema_version = 2\n[associations.water]\nscheme = "4C"\n'
        '[associations.water.a_srk]\nvalue = 0.12277\nunit = "dimensionless"\n'
    )

    with pytest.raises(KeycardError, match="cannot be read as"):
        python_card(text)
    with pytest.raises(InvalidInputError, match=re.escape("Pa*m**6/mol**2")):
        rust_card(text)


#: One pair, three columns. `water/methanol` ships `-0.0789` classical against `-0.153`
#: associating, so a reader that read one column for the other is visible; and the two
#: associating columns are overridden to different values, so a reader that keyed them by
#: pair alone is visible too.
TRANSPORTED_KIJ_CARD = (
    "schema_version = 2\n"
    "[[kij]]\n"
    'component_a = "water"\ncomponent_b = "methanol"\nvalue = -0.0789\n'
    "cpa_value_srk = -0.08\ncpa_value_pr = -0.31\n"
)


def test_the_two_readers_resolve_the_same_associating_pairs() -> None:
    """The associating interaction columns, resolved on both sides.

    A separate comparison from the classical one because they are separate columns of the
    same table, and because the family is what selects between them: the two are fits
    rather than one converted, so a card may state one and not the other.
    """
    card = python_card(TRANSPORTED_KIJ_CARD)
    overlay = rust_card(TRANSPORTED_KIJ_CARD)

    # In the order the Rust pairs come in - the lower-sorting name first - so the pair the
    # comparison names is the pair `cpa_kij_for` indexes by.
    names = ("methanol", "water")
    for family in components.CPA_FAMILIES:
        theirs = _core.overlay_cpa_kij_rows(overlay, family)
        assert theirs, f"the card states no pair, so nothing was compared for {family}"

        mine = components.cpa_kij_for(names, family, card=card)
        for first, second, value in theirs:
            index = (names.index(first), names.index(second))
            assert index in mine, f"{family} {first}/{second} resolved to nothing"
            assert value == pytest.approx(mine[index]), f"{family} {first}/{second}"

    assert card.cpa_kij_for("methanol", "water", "srk") == pytest.approx(-0.08)
    assert card.cpa_kij_for("methanol", "water", "pr") == pytest.approx(-0.31)


def test_a_card_states_a_cpa_pair_without_the_classical_one_only_by_stating_both() -> None:
    """`value` is required, and both readers say so.

    An absent classical value would leave the shipped one in force, so a row carrying only
    a `cpa_value_*` would look like it had overridden a column it had not.
    """
    text = (
        'schema_version = 2\n[[kij]]\ncomponent_a = "water"\n'
        'component_b = "methanol"\ncpa_value_srk = -0.08\n'
    )

    with pytest.raises(KeycardError, match="value"):
        python_card(text)
    with pytest.raises(InvalidInputError, match="value"):
        rust_card(text)


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


A_MATRIX_CARD = """\
schema_version = 2
keyholder.name = "Example Engineering Ltd"

# Deliberately not square: a square matrix cannot tell a transposition from the truth,
# so a comparison over one would pass whichever way round the two readers built it.
[coefficients."eos.uniquac_activity_coefficients".aij]
value = [[0.0, -71.0, 5.5], [209.0, 0.0, 7.5]]
unit = "K"
citation = "a representative UNIQUAC matrix; NeqSim carries no such table"

[coefficients."hydraulics.orifice_flow".Cd]
value = 0.61
unit = "dimensionless"

[coefficients."hydraulics.choked_flow_area".d]
value = 50.0
unit = "mm"
"""


def _python_coefficients(card: keycard.Keycard) -> dict[tuple[str, str], tuple[str, list[Any]]]:
    """Each coefficient as `(unit, rows)`, a scalar and a vector both being rows.

    The shape is the point: a matrix flattened to a list of numbers would compare equal
    to a vector of the same numbers, which is the one mistake this comparison is here
    to catch.
    """
    out: dict[tuple[str, str], tuple[str, list[Any]]] = {}
    for calc_id, arguments in card.coefficients.items():
        for name, value in arguments.items():
            # A vector and a matrix are lists of quantities, one per entry - the
            # convention every declared vector and matrix input follows. `rows` is the
            # shape the Rust side reports: a matrix's own rows, a vector as one row, a
            # number as a one-entry row.
            nested = _magnitudes(value)
            if not isinstance(nested, list):
                rows: list[Any] = [nested]
            elif nested and not isinstance(nested[0], list):
                rows = [nested]
            else:
                rows = nested
            out[(calc_id, name)] = (_spelling(value), rows)
    return out


def _spelling(value: Any) -> str:
    """The unit as the *format* spells it, which is what the Rust reader reports.

    `pint` names the unit `kelvin` and the vocabulary names it `K`; both readers were
    handed `K` and validated it against that vocabulary, so the comparison is of the
    spelling the format uses rather than of either library's own name for it.
    """
    name = str(_leaves(value)[0].units)
    for spelled, canonical in CANONICAL_UNITS.items():
        if canonical == name:
            return spelled
    raise AssertionError(f"{name} is not a spelling this vocabulary carries")


def _magnitudes(value: Any) -> Any:
    """A coefficient's SI base magnitudes, in the shape it was declared in.

    Base units rather than as written, because that is the number a calculation uses -
    the same choice `overlay_kij_rows` makes when it compares the *effective*
    interaction parameter rather than the stated one. A `50 mm` coefficient is `0.05`,
    and a comparison over the written numbers would call two readers that disagree
    about the conversion equal.
    """
    if isinstance(value, list):
        return [_magnitudes(entry) for entry in value]
    return float(value.to_base_units().magnitude)


def _leaves(value: Any) -> list[Any]:
    """Every quantity in a coefficient, in order, so its unit can be read once."""
    if isinstance(value, list):
        return [leaf for entry in value for leaf in _leaves(entry)]
    return [value]


def test_the_two_readers_resolve_the_same_coefficients() -> None:
    """The card format's third section, compared value by value and shape by shape.

    `components` and `kij` had a comparison; `coefficients` had none, so a value one
    reader accepted and the other refused - or one read as a matrix and the other as a
    vector - was checked by nothing. The card the format's *own* schema admits is now
    read by both and the two answers compared.
    """
    card = python_card(A_MATRIX_CARD)
    mine = _python_coefficients(card)
    theirs = {
        (calc_id, name): (unit, _rows(rows, cols, values))
        for calc_id, name, unit, rows, cols, values in _core.card_coefficients(A_MATRIX_CARD)
    }

    assert set(mine) == set(theirs), "the readers name different coefficients"
    assert mine, "the card states no coefficients, so nothing was compared"
    for key in sorted(mine):
        unit, rows = mine[key]
        assert unit == theirs[key][0], f"{key}: unit {unit} vs {theirs[key][0]}"
        assert rows == theirs[key][1], f"{key}: {rows} vs {theirs[key][1]}"


def _rows(rows: int, cols: int, values: list[float]) -> list[Any]:
    """The Rust side's `(rows, cols, values)` as the nested list Python builds."""
    if rows == 1 and cols == 1:
        return [values[0]]
    if rows == 1:
        return [list(values)]
    return [list(values[i * cols : (i + 1) * cols]) for i in range(rows)]


def test_both_readers_refuse_a_ragged_coefficient() -> None:
    """A shape neither reader will read, refused by both rather than by one.

    Serde accepts `[[1, 2], [3]]` - it is a list of lists of numbers - and `pint` will
    build a quantity from it, so a shape check that lived on one side only would let the
    other through with a matrix the calculation cannot index.
    """
    ragged = (
        'schema_version = 2\n[coefficients."eos.uniquac_activity_coefficients".aij]\n'
        'value = [[0.0, -71.0], [209.0]]\nunit = "K"\n'
    )
    with pytest.raises(KeycardError, match="rectangular"):
        python_card(ragged)
    with pytest.raises(InvalidInputError, match="rectangular"):
        _core.card_coefficients(ragged)


def test_both_readers_refuse_an_empty_coefficient() -> None:
    """A value with no numbers is not a default, on either side."""
    empty = (
        'schema_version = 2\n[coefficients."eos.uniquac_activity_coefficients".aij]\n'
        'value = []\nunit = "K"\n'
    )
    with pytest.raises(KeycardError, match="empty"):
        python_card(empty)
    with pytest.raises(InvalidInputError, match="rectangular"):
        _core.card_coefficients(empty)


def test_a_card_stating_a_zero_pair_keeps_it_on_both_sides() -> None:
    """A card's `value = 0.0` means "reset this pair to ideal mixing", not "no value".

    The databank's own zero *is* an absence - it is the ideal-mixture default - so the
    two are distinguished by where the value came from, not by the number. Read as a
    number the distinction is lost, and a caller overriding a fitted pair back to ideal
    mixing would be silently given the fitted one back.
    """
    text = (
        'schema_version = 2\n[[kij]]\ncomponent_a = "methane"\n'
        'component_b = "n-butane"\nvalue = 0.0\n'
    )
    card = python_card(text)
    mine = components.kij_for(("methane", "n-butane"), card=card)
    assert mine == {(0, 1): 0.0}, "the card's zero must survive the Python reader"

    theirs = _core.overlay_kij_rows(rust_card(text))
    assert [(a, b, value) for a, b, value in theirs] == [("methane", "n-butane", 0.0)], (
        "and the Rust reader's, so the two agree"
    )

    # The pair is genuinely one the table fits, so a zero is a change rather than a
    # restatement of what is already there.
    assert components.kij_for(("methane", "n-butane")) != {(0, 1): 0.0}
