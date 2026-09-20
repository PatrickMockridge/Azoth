"""The PHREEQC Pitzer catalogue and the rule that chooses it, as Python carries them.

The Python mirror of the Rust tests in `crates/azoth-eos/src/pitzer_catalog.rs`, asserting
the same constants. The two measured cases - a brine the catalogue covers and one it does
not - are reproduced without a `SystemPitzer`, because the rule is about a phase's
*topology* rather than its state.
"""

from __future__ import annotations

import pytest

from azoth.core.errors import InvalidInputError
from azoth.eos.reference import _pitzer_catalog as catalog
from azoth.eos.reference import _pitzer_phase as phase


def test_the_catalogue_is_read() -> None:
    """The probe's own reading: `b0(Na+,Cl-) = 0.07534` under PHREEQC."""
    sodium_chloride = catalog.find("B0", ["Na+", "Cl-"])
    assert sodium_chloride is not None
    assert sodium_chloride[0] == pytest.approx(0.07534, abs=1.0e-12)

    # **Either order round**, because the key is sorted - the catalogue writes `Cl- Na+`
    # and a caller holds the cation first.
    assert catalog.find("B0", ["Cl-", "Na+"]) == sodium_chloride

    # **And PHREEQC's spelling of a charge is canonicalised**, so a caller holding NeqSim's
    # `Ba++` finds the row the catalogue wrote as `Ba+2`.
    assert catalog.find("B0", ["Ba++", "Cl-"]) is not None
    assert catalog.canonical_species("Ba+2") == "Ba++"
    assert catalog.canonical_species("SO4-2") == "SO4--"
    assert catalog.canonical_species("Fe+3") == "Fe+++"
    # `+3` is tested before `+2`, or `Fe+3` would come back `Fe+` with a stray `3`.
    assert catalog.canonical_species("B(OH)4-") == "B(OH)4-"

    # The families the enum declares and the file does not carry are empty, which is a
    # different answer from a family that does not exist.
    for family in ("MU", "ETA", "ALPHAS"):
        assert catalog.family_rows(family) == ()
    assert len(catalog.family_rows("B0")) == 54
    assert sum(len(catalog.family_rows(f)) for f in catalog.FAMILY_SPECIES) == 268


def _species(
    ions: list[tuple[str, float]], neutrals: list[tuple[str, str]]
) -> list[catalog.Species]:
    out = [catalog.Species("water", 0.9, 0.0, "H2O", False)]
    out += [catalog.Species(name, 0.01, 0.0, formula, False) for name, formula in neutrals]
    out += [catalog.Species(name, 0.05, charge, "", False) for name, charge in ions]
    return out


def test_an_ion_free_phase_falls_back_before_the_catalogue_is_consulted() -> None:
    """**Why `water + CO2` takes the CSV** although the catalogue carries `CO2|CO2`.

    `tryApplyCompletePhreeqcPitzerCatalog` returns `false` on an empty ion list before it
    looks at the catalogue at all, so this is not a coverage question.
    """
    assert catalog.find("LAMBDA", ["CO2", "CO2"]) is not None
    assert catalog.select_dataset(_species([], [("CO2", "CO2")])) == ("legacy", "NoIons")


def test_the_coverage_rule_decides_the_dataset() -> None:
    """A brine the catalogue covers takes it; one it does not falls back."""
    assert catalog.select_dataset(_species([("Na+", 1.0), ("Cl-", -1.0)], [])) == ("phreeqc", None)

    # **Hydrogen carbonate is the measured case**: `B0` and `B1` exist, `C0` does not, and
    # the whole dataset is abandoned rather than completed from the other one.
    assert catalog.find("B0", ["Na+", "HCO3-"]) is not None
    assert catalog.find("B1", ["Na+", "HCO3-"]) is not None
    assert catalog.find("C0", ["Na+", "HCO3-"]) is None
    assert catalog.select_dataset(_species([("Na+", 1.0), ("HCO3-", -1.0)], [])) == (
        "legacy",
        "C0 for Na+, HCO3-",
    )


def test_an_optional_family_does_not_gate_the_selection() -> None:
    """`B2` is looked up rather than required, so its presence is not a gate."""
    assert catalog.find("B2", ["Mg++", "SO4--"]) is not None
    assert catalog.select_dataset(_species([("Mg++", 2.0), ("SO4--", -2.0)], [])) == (
        "phreeqc",
        None,
    )


def test_a_hydrocarbon_is_excluded_from_the_neutral_topology() -> None:
    """The formula test is what the method exists for: methane keeps the type `normal`."""
    methane = catalog.Species("methane", 0.1, 0.0, "CH4", False)
    assert catalog.is_hydrocarbon(methane)
    # `CO2` has an oxygen, so it is an active neutral - and the catalogue covers it.
    assert not catalog.is_hydrocarbon(catalog.Species("CO2", 0.1, 0.0, "CO2", False))
    assert catalog.is_hydrocarbon(catalog.Species("default", 0.1, 0.0, "", True))


def test_a_trace_component_does_not_enter_the_topology() -> None:
    """Lithium has no catalogue rows, so if the trace entry counted this would fall back."""
    mixture = _species([("Na+", 1.0), ("Cl-", -1.0)], [])
    mixture.append(catalog.Species("Li+", 1.0e-25, 1.0, "", False))
    assert catalog.select_dataset(mixture) == ("phreeqc", None)


def test_both_temperature_forms_match_neqsim() -> None:
    """The two forms against the built phase, from the probe's own columns.

    `validation/neqsim/PitzerArithmetic.java` prints each across temperature from a live
    `SystemPitzer`. The catalogue pair is `Na+/Cl-` (six coefficients) and the CSV pair is
    `Na+/HCO3-`, whose two coefficients are zero and which is therefore flat.
    """
    pair = catalog.find("B0", ["Na+", "Cl-"])
    assert pair is not None

    # At the reference the value is `a0`, and it stays `a0` inside the 1e-3 K window.
    assert phase.catalog_value(pair, 298.15) == pytest.approx(0.07534, abs=1e-15)
    assert phase.catalog_value(pair, 298.1501) == pytest.approx(0.07534, abs=1e-15)

    for temperature, expected in (
        (273.15, 0.0493895679117),
        (323.15, 0.0892387746618),
        (373.15, 0.100154205617),
    ):
        assert phase.catalog_value(pair, temperature) == pytest.approx(expected, abs=1e-12)

    # The CSV's flat case: zero `t1`/`t2` gives the same number at every temperature.
    for temperature in (273.15, 298.15, 323.15, 373.15):
        assert phase.silvester_value(0.0277, 0.0, 0.0, temperature) == 0.0277


def test_the_alpha_coefficients_follow_the_charge() -> None:
    """The oracle prints `Na+/Cl-` as `alpha1 = 2.0`, `alpha2 = 12.0`."""
    assert phase.alpha1_as_used(1.0, -1.0) == 2.0
    assert phase.alpha2(1.0, -1.0) == 12.0
    assert phase.alpha1_as_used(2.0, -2.0) == 1.4
    assert phase.alpha2(2.0, -2.0) == 12.0
    assert phase.alpha1_as_used(2.0, -1.0) == 2.0
    assert phase.alpha2(2.0, -1.0) == 12.0
    assert phase.alpha2(3.0, -1.0) == 12.0, "the chloride is monovalent"

    # **The one charge where the two alpha1 definitions part**, unreachable with vendored
    # data because the highest charge it carries is 2.
    assert phase.alpha1_as_used(3.0, -3.0) == 1.4
    assert phase.alpha1(3.0, -3.0) == 2.0


def _brine(entries: list[tuple[str, float, float]]) -> tuple[list[catalog.Species], float]:
    """One mole of mixture's worth, and the solvent mass the audit's threshold scales by."""
    water = 0.9
    out = [catalog.Species("water", water, 0.0, "H2O", False)]
    out += [catalog.Species(name, moles, charge, "", False) for name, charge, moles in entries]
    return out, water * 0.018015


def test_a_mixed_catalogue_brine_is_complete() -> None:
    """The catalogue covers a mixed brine, so nothing is missing in any family."""
    mixture, mass = _brine(
        [("Na+", 1.0, 0.01), ("K+", 1.0, 0.01), ("Cl-", -1.0, 0.01), ("SO4--", -2.0, 0.01)]
    )
    selection = catalog.select_dataset(mixture)
    assert selection == ("phreeqc", None)
    audit = catalog.coverage(mixture, selection, mass)
    assert audit.is_complete(), audit
    assert audit.active_cations == ("K+", "Na+")
    assert audit.active_anions == ("Cl-", "SO4--")
    assert audit.dataset_id == catalog.PHREEQC_DATASET_ID
    catalog.require_complete(audit)


def test_a_pair_with_no_parameters_is_refused() -> None:
    """**A pair with no parameters is refused, not evaluated at zero.**

    The sabotage check: `NH4+`/`Cl-` appears in neither dataset, and the audit must say so
    rather than let the model compute a substance as an ideal solution.
    """
    mixture, mass = _brine([("NH4+", 1.0, 0.01), ("Cl-", -1.0, 0.01)])
    selection = catalog.select_dataset(mixture)
    assert selection == ("legacy", "B0 for NH4+, Cl-"), "no ammonium rows in the catalogue"
    assert catalog.find("B0", ["NH4+", "Cl-"]) is None

    audit = catalog.coverage(mixture, selection, mass)
    assert not audit.is_complete()
    assert audit.missing_binary == ("Cl-|NH4+",)
    assert audit.missing_theta == ()
    assert audit.missing_psi == ()

    # **The refusal, and its exact wording** - the string is NeqSim's, because a reader
    # comparing a port's failure against NeqSim's exception must not have to tell two
    # formats apart.
    with pytest.raises(InvalidInputError) as raised:
        catalog.require_complete(audit)
    assert raised.value.reason == audit.diagnostic()

    # **And NeqSim would still initialize it.** A single cation and a single anion is not a
    # mixed topology, and `validateParameterCoverageOncePerState` only enforces the audit
    # when it is - so this brine evaluates the pair at zero.
    assert not catalog.has_mixed_primary_salt_topology(mixture, mass)


def test_a_mixed_legacy_brine_is_incomplete() -> None:
    """**Measured: NeqSim throws from `init(1)` here.** The legacy CSV has no same-sign rows."""
    mixture, mass = _brine(
        [
            ("Na+", 1.0, 0.01),
            ("K+", 1.0, 0.01),
            ("Cl-", -1.0, 0.01),
            ("HCO3-", -1.0, 0.01),
        ]
    )
    selection = catalog.select_dataset(mixture)
    assert selection[0] == "legacy", "hydrogen carbonate is in no catalogue family"
    audit = catalog.coverage(mixture, selection, mass)

    # **`HCO3-` is absent from `active_anions`** although it is the larger of the two
    # anions, because the primary-salt audit excludes the reaction species.
    assert audit.active_cations == ("K+", "Na+")
    assert audit.active_anions == ("Cl-",)
    assert audit.missing_binary == ()
    assert audit.missing_theta == ("K+|Na+",)
    assert audit.missing_psi == ("K+|Na+|Cl-",)
    assert catalog.has_mixed_primary_salt_topology(mixture, mass)
    assert audit.diagnostic() == (
        "Pitzer parameter coverage incomplete for dataset "
        "'neqsim-legacy-pitzer-parameters-v1': activeCations=[K+, Na+], "
        "activeAnions=[Cl-], missingBinary=[], missingTheta=[K+|Na+], "
        "missingPsi=[K+|Na+|Cl-]"
    )


def test_the_reaction_audit_keeps_the_species_the_primary_audit_drops() -> None:
    """The reaction variant is a different observable, which is why both are exposed."""
    mixture, mass = _brine([("Na+", 1.0, 0.01), ("HCO3-", -1.0, 0.01), ("CO3--", -2.0, 0.01)])
    selection = catalog.select_dataset(mixture)

    primary = catalog.coverage(mixture, selection, mass)
    assert primary.active_anions == ()
    assert primary.is_complete()

    reaction = catalog.reaction_coverage(mixture, selection, mass)
    assert reaction.active_anions == ("CO3--", "HCO3-")
    assert not reaction.is_complete()
    assert reaction.missing_theta == ("CO3--|HCO3-",)
    assert reaction.missing_psi == ("CO3--|HCO3-|Na+",)


def test_a_trace_ion_leaves_the_audited_topology() -> None:
    """**The audit's threshold is a molality and the selection's is a mole count.**

    A trace ion can have chosen the dataset and still not be audited.
    """
    mixture, mass = _brine([("Na+", 1.0, 0.01), ("Cl-", -1.0, 0.01), ("K+", 1.0, 1.0e-11)])
    molality = 1.0e-11 / mass
    assert catalog.ACTIVE_MOLES < molality < catalog.ACTIVE_ION_MOLALITY

    selection = catalog.select_dataset(mixture)
    assert selection == ("phreeqc", None)
    audit = catalog.coverage(mixture, selection, mass)
    assert audit.active_cations == ("Na+",)
    assert not catalog.has_mixed_primary_salt_topology(mixture, mass)
