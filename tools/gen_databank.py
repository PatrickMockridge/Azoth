#!/usr/bin/env python3
"""Generate `data/components/` from NeqSim's component and interaction tables.

    python tools/gen_databank.py --check                     # against the vendored copies
    python tools/gen_databank.py /tmp/neqsim-check/neqsim    # against a fresh checkout

The default argument is `databank/sources/neqsim/`, which holds NeqSim's `COMP.csv`
and `INTER.csv` verbatim; a checkout refreshes those two when NeqSim is bumped, and
every other run reads what is committed. `NOTICE` carries the attribution, and each
row's `citation` names the source once because that is the column the loader reads.

# The units

Read out of NeqSim's own loader, `thermo/component/Component.java` around line 305,
rather than guessed from the column names:

    molarmass  divided by 1000        -> kg/mol      (the file is g/mol)
    TC         plus 273.15            -> K           (the file is degrees Celsius)
    PC         unchanged              -> bar         (see below)
    critvol    unchanged              -> cm**3/mol
    liqdens    unchanged              -> g/cm**3

Each is applied once, here, and the result is asserted against values known
independently of NeqSim.

`PC` is read out rather than copied from a comment. NeqSim's
`criticalCompressibilityFactor` is `Pc * Vc / R / Tc / 10`, and methane's
`Z_c = 0.2874` falls out of that arithmetic with `Pc` in bar and `Vc` in cm**3/mol -
which fixes those two and rules out Pa, kPa and m**3/mol.

# What is carried across

`COMPONENT_COLUMNS` below is the list. The reasons are not here: every one of
`COMP.csv`'s 170 columns and `INTER.csv`'s 39 is dispositioned in
`databank/manifest.toml`, and this tool asserts at startup that the two agree.
`COMP_EXT.csv` (86 MB, 76,705 rows) is not vendored at all, and the manifest says why.
"""

from __future__ import annotations

import argparse
import csv
import io
import math
import re
import sys
from collections.abc import Callable
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import manifest as manifest_module

ROOT = Path(__file__).resolve().parent.parent
SOURCES = ROOT / "databank" / "sources" / "neqsim"
OUT_DIR = ROOT / "data" / "components"

#: The small tables carried verbatim, re-rendered like the component and interaction
#: tables so their compiled headers match the manifest's `as` names: group-contribution
#: parameters and the MBWR-32 coefficients need no unit conversion and no row filtering.
#: Each entry is `(manifest id, source file, compiled name)`.
VENDORED_FILES = (
    ("neqsim/UNIFACcomp.csv", "UNIFACcomp.csv", "UNIFACcomp.csv"),
    ("neqsim/UNIFACGroupParam.csv", "UNIFACGroupParam.csv", "UNIFACGroupParam.csv"),
    ("neqsim/UNIFACInterParam.csv", "UNIFACInterParam.csv", "UNIFACInterParam.csv"),
    ("neqsim/UNIFACInterParamB.csv", "UNIFACInterParamB.csv", "UNIFACInterParamB.csv"),
    ("neqsim/UNIFACInterParamC.csv", "UNIFACInterParamC.csv", "UNIFACInterParamC.csv"),
    ("neqsim/UNIFACcompUMRPRU.csv", "UNIFACcompUMRPRU.csv", "UNIFACcompUMRPRU.csv"),
    ("neqsim/UNIFACInterParamA_UMR.csv", "UNIFACInterParamA_UMR.csv", "UNIFACInterParamA_UMR.csv"),
    (
        "neqsim/UNIFACInterParamA_UMRMC.csv",
        "UNIFACInterParamA_UMRMC.csv",
        "UNIFACInterParamA_UMRMC.csv",
    ),
    ("neqsim/UNIFACInterParamB_UMR.csv", "UNIFACInterParamB_UMR.csv", "UNIFACInterParamB_UMR.csv"),
    (
        "neqsim/UNIFACInterParamB_UMRMC.csv",
        "UNIFACInterParamB_UMRMC.csv",
        "UNIFACInterParamB_UMRMC.csv",
    ),
    ("neqsim/UNIFACInterParamC_UMR.csv", "UNIFACInterParamC_UMR.csv", "UNIFACInterParamC_UMR.csv"),
    (
        "neqsim/UNIFACInterParamC_UMRMC.csv",
        "UNIFACInterParamC_UMRMC.csv",
        "UNIFACInterParamC_UMRMC.csv",
    ),
    ("neqsim/MBWR32param.csv", "MBWR32param.csv", "mbwr32.csv"),
    # The electrolyte pair parameters. Same shape as the rest of this list - no unit
    # conversion and no row filtering - and compiled rather than read from
    # `databank/sources/` for the same reason: `sources/` is the vendored upstream copy,
    # kept for reproducibility and not shipped, and a wheel carries `data/`.
    ("neqsim/PitzerParameters.csv", "PitzerParameters.csv", "PitzerParameters.csv"),
    ("neqsim/COMPSALT.csv", "COMPSALT.csv", "COMPSALT.csv"),
)


#: The PHREEQC Pitzer catalogue NeqSim bundles, and the families it carries.
#:
#: **A `.dat` and not a CSV**, so it has its own parser rather than a `VENDORED_FILES`
#: entry. The format is PHREEQC's: a `-FAMILY` line opens a section, and each row below it
#: is `species... coefficients...`, whitespace-separated, up to six coefficients.
#:
#: **The species names are PHREEQC's and the lookup key is NeqSim's.**
#: `PhreeqcPitzerParameterCatalog` rewrites a numeric suffix - `Ba+2` becomes `Ba++`,
#: `SO4-2` becomes `SO4--` - and then sorts the names and joins them with `|`, so a row is
#: found whichever order its species are written in. That canonical form is what the
#: compiled table carries, because it is what the model looks up by.
PHREEQC_CATALOG = (
    "neqsim/thermo/phase/phreeqc/pitzer-b0b3be767158ccc3322d2c816625cf470045e67e.dat",
    "phreeqc/pitzer-b0b3be767158ccc3322d2c816625cf470045e67e.dat",
    "PitzerPhreeqc.csv",
)

#: The Fürst electrolyte's fitted parameters, which are **Java rather than data**.
#:
#: `FurstElectrolyteConstants` is 937 lines of hardcoded `static double[]` and does no file
#: I/O, so there is no upstream CSV to carry - the model's parameters exist only as source
#: code, and `ComponentModifiedFurstElectrolyteEos`'s constructor reads `furstParams[0]` and
#: `[1]` to build an ion's covolume. **This is the first vendored source that is not
#: tabular**, and it is vendored as bytes rather than transcribed for the reason every other
#: source is: a transcription that drifts from NeqSim is undetectable, and a vendored file
#: can be diffed against the pinned commit.
#:
#: `PHREEQC_CATALOG`'s shape - a source, a vendored path, an output - because it has the
#: same property: a format of its own and therefore its own parser.
FURST_CONSTANTS = (
    "neqsim/src/main/java/neqsim/thermo/util/constants/FurstElectrolyteConstants.java",
    "FurstElectrolyteConstants.java",
    "furst_parameters.csv",
)

#: ISO 6976's per-component gas-quality constants, compiled to `data/standards/`.
#:
#: **The older of NeqSim's two tables, because `Stream.LCV()` reads it.** NeqSim ships
#: `ISO6976constants.csv` and `ISO6976constants2016.csv` and has a class for each -
#: `Standard_ISO6976` queries `ISO6976constants`, and `Standard_ISO6976_2016 extends
#: Standard_ISO6976` queries the revision. **They disagree**: on methane the older file has
#: `Z15 = 0.998000` against `0.998020`, `srtb15 = 0.044700` against `0.044520`, a molar mass
#: of `16.043000` against `16.042460` g/mol and a superior calorific value of `890.630`
#: against `890.580` kJ/mol at 25 °C - a revision of the constants moves them, which is what
#: a revision is for. `unit_ops.flare`'s duty is
#: `Stream.LCV()` and `Stream.LCV()` builds `new Standard_ISO6976(fluid, 0, 15.55,
#: "volume")`, so the table the flare reads is this one. Compiling the 2016 file instead
#: would be porting the sibling class, and the port would then disagree with the machine it
#: is porting.
ISO6976_CONSTANTS = ("ISO6976constants.csv", "iso6976.csv")

#: Where the compiled standard constants go.
STANDARDS_DIR = ROOT / "data" / "standards"

#: `GibbsReactor`'s own species database: the element vectors and the per-species
#: Gibbs/enthalpy polynomials, compiled to `data/reactors/`.
#:
#: **The two files are one table and are compiled together**, because neither is a row on
#: its own: `GibbsReactDatabase.csv` carries the element vector and the formation
#: properties, `DatabaseGibbsFreeEnergyCoeff.csv` the order-5 Gibbs and enthalpy
#: polynomials, and `GibbsReactor.loadGibbsDatabase` joins them on the lowercased molecule
#: name. A species in the first and not the second takes the class's fallback branch.
#:
#: **This is not the databank's formation properties.** `GIBBSENERGYOFFORMATION` is already
#: `used` from P10 and its values are not these; the class reads its own table, and a port
#: composed from the databank would answer with different numbers than the machine it ports.
GIBBS_REACTOR = (
    "GibbsReactDatabase/GibbsReactDatabase.csv",
    "gibbs_reactor.csv",
)

#: The per-species order-5 Gibbs and enthalpy polynomials, compiled to its own file.
#:
#: **Two compiled files rather than one joined table, because that is what the class
#: holds**: `loadGibbsDatabase` builds `extraCoeffMap` from this file and `componentMap`
#: from the other, and joins them at lookup time by molecule name, storing all-NaN where
#: this file has no row. Joining them here would put that fallback in the data instead of
#: in the reader, and a species this file omits would then have to be written as an empty
#: cell rather than being absent.
GIBBS_COEFFS = (
    "GibbsReactDatabase/DatabaseGibbsFreeEnergyCoeff.csv",
    "gibbs_reactor_coeffs.csv",
)

#: Where the compiled reactor data goes.
REACTORS_DIR = ROOT / "data" / "reactors"

#: The compiled table's columns, in order.
#:
#: **All eight element counts are carried, under the CSV's own names.** The class declares
#: `elementNames = {"O","N","C","H","S","Ar","Z"}` - seven - and reads indices 0..6 of an
#: array the loader fills with eight, so its `Ar` is this file's `Na` and its `Z` is this
#: file's `Ar`, and this file's `Z` is read by nothing. Compiling the eight under their own
#: names keeps that a fact about the kernel's indexing rather than about the data, which is
#: where it has to be for the capture to pin it with a number.
#:
#: **`a`..`d` are not carried.** The class parses them, stores them, clones them in a getter
#: and reads them nowhere: the fallback Gibbs branch calls
#: `calculateCorrectedHeatCapacityCoeffs`, which reads `system.getComponent(i).getCpA()`
#: from NeqSim's own component database. Four columns nothing reads are four columns the
#: manifest dispositions as `not-a-value` rather than four this file copies.
GIBBS_REACTOR_HEADER: tuple[str, ...] = (
    "name",
    "o",
    "n",
    "c",
    "h",
    "s",
    "na",
    "ar",
    "z",
    "hf298",
    "gf298",
    "sf298",
    "citation",
)

#: The coefficient table's columns, in the upstream file's own order.
GIBBS_COEFF_HEADER: tuple[str, ...] = (
    "name",
    "ag",
    "bg",
    "cg",
    "dg",
    "eg",
    "fg",
    "ah",
    "bh",
    "ch",
    "dh",
    "eh",
    "gh",
    "citation",
)

#: How many species each family's rows name, from `PhreeqcPitzerParameterCatalog.Family`.
#: A family the file does not carry compiles to no rows, which is how `MU`, `ETA` and
#: `ALPHAS` come out: the enum declares them and the shipped catalogue has none.
PHREEQC_FAMILIES: dict[str, int] = {
    "B0": 2,
    "B1": 2,
    "B2": 2,
    "C0": 2,
    "THETA": 2,
    "PSI": 3,
    "LAMBDA": 2,
    "ZETA": 3,
    "MU": 3,
    "ETA": 3,
    "ALPHAS": 2,
}


#: The NeqSim release this was generated from, for the `citation` column and for
#: `NOTICE`. A version and a commit, because a databank is only reproducible against
#: one revision of the file it came from.
#:
#: **Read from `databank/manifest.toml` rather than restated.** What a row's `citation`
#: names and what the manifest records were two copies of one fact, and they had come
#: apart: the manifest said `805cf0f` while this said `dedba873`, so the citation in every
#: compiled row named a revision the vendored files do not match. The manifest is the
#: record of what was taken - the same reason `COMPONENT_COLUMNS` is read from it.
def _upstream(key: str) -> str:
    """One field of the NeqSim upstream, or a failure that says which."""
    found, problems = manifest_module.read()
    if problems:
        raise SystemExit("gen_databank: " + "\n  ".join(problems))
    for upstream in found.upstreams:
        if upstream.id == "neqsim":
            return str(getattr(upstream, key))
    raise SystemExit("gen_databank: the manifest declares no `neqsim` upstream")


NEQSIM_VERSION = _upstream("version")
NEQSIM_COMMIT = _upstream("commit")
#: **The revision, named by the commit rather than the version.** The manifest's `version`
#: is NeqSim's own `pom.xml` property, which reads the same at the release tag and on master
#: - fifty commits apart - so a citation built from it names two different trees. The commit
#: is the pin, and `master` is the only readable name that resolves to the bytes vendored
#: here. The version is kept beside it because a reader recognises it.
#: The ISO 6976 table's own citation: the same upstream pin, the file that answered.
ISO6976_CITATION = (
    f"NeqSim master ({NEQSIM_COMMIT[:7]}) ISO6976constants.csv (Equinor/NTNU), Apache-2.0, "
    f"retrieved {_upstream('retrieved')}"
)

#: The reactor database's own citation: one pin and two files, because the table is the
#: join of them.
GIBBS_REACTOR_CITATION = (
    f"NeqSim master ({NEQSIM_COMMIT[:7]}) GibbsReactDatabase.csv + "
    f"DatabaseGibbsFreeEnergyCoeff.csv (Equinor/NTNU), Apache-2.0, "
    f"retrieved {_upstream('retrieved')}"
)

CITATION = (
    f"NeqSim master ({NEQSIM_COMMIT[:7]}) COMP.csv (Equinor/NTNU), Apache-2.0, "
    f"retrieved {_upstream('retrieved')}"
)

#: Output column -> (NeqSim column, conversion). The conversions are the ones
#: documented above, each with the reason it is that and not another.
#: How an upstream value becomes the unit the manifest records for that column.
#:
#: Every entry is a place NeqSim's file and this library disagree about a unit, and the
#: comment names both sides. A column absent from this table is carried exactly as
#: NeqSim stores it, which the manifest records as `neqsim-internal` - not a guess, and
#: countable, so the columns whose unit nobody has established are visible.
CONVERSIONS: dict[str, Callable[[str], object]] = {
    "MOLARMASS": lambda v: float(v) / 1000.0,  # g/mol -> kg/mol
    "TC": lambda v: float(v) + 273.15,  # degC -> K
    "PC": lambda v: float(v) * 1.0e5,  # bar -> Pa
    "CRITVOL": lambda v: float(v) * 1.0e-6,  # cm3/mol -> m3/mol
    "LIQDENS": lambda v: float(v) * 1000.0,  # g/cm3 -> kg/m3
    "NORMBOIL": lambda v: float(v) + 273.15,  # degC -> K
    "sigmaSAFT": lambda v: float(v) / 1.0e10,  # Angstrom -> m
    "sigmaSAFTVRMie": lambda v: float(v) / 1.0e10,  # Angstrom -> m
    # A defect in the upstream file, not a unit: six rows of this column are written
    # with a comma decimal separator - `propane/CO2 = 0,1241` and five more alkane/CO2
    # pairs. NeqSim's own reader is `Double.parseDouble`, which throws on every one, so
    # the file has never been read along that path. The values rise monotonically with
    # carbon number, which is what a Soreide-Whitson kij does, so the comma is a decimal
    # point and the fix is stated here rather than applied silently to the source.
    "KIJWhitsonSoriede": lambda v: float(v.replace(",", ".")),
}


def read_rows(path: Path) -> list[dict[str, str]]:
    """Rows of a NeqSim resource CSV, decoded strictly.

    Strictly rather than with `errors="replace"`, which is what this did. A replaced
    byte is silent: a name or a formula arrives with U+FFFD in it and the row parses
    anyway, so the corruption ships. The unit conversions are guarded by `ROUND_TRIP`
    and the encoding was guarded by nothing.
    """
    try:
        with path.open(newline="", encoding="utf-8") as handle:
            return [row for row in csv.DictReader(handle) if row.get("NAME") or row.get("COMP1")]
    except UnicodeDecodeError as error:
        raise ValueError(
            f"{path}: not valid UTF-8 (byte {error.start}). Read strictly because a "
            f"replaced byte is an undetectable corruption in a shipped component row."
        ) from error


def _text_columns(source: Path) -> frozenset[str]:
    """The columns whose values are not all numbers, read from the file itself.

    **Read rather than listed**, because a list is a second copy of the file that can be
    wrong: `associationscheme` holds `4C`, `HydrateFormer` holds `yes`, and a column
    coerced to a number is a row silently dropped - which is how the first run of this
    lost two components and three hundred interaction pairs without saying so.
    """
    text: set[str] = set()
    for row in read_rows(source):
        for name, value in row.items():
            if name in text:
                continue
            stripped = (value or "").strip()
            if not stripped:
                continue
            try:
                float(stripped)
            except ValueError:
                text.add(name)
    return frozenset(text)


#: **Both sources, not just the component table.** This was `COMP.csv` alone, and the
#: consequence was that `INTER.csv`'s `HVTYPE` and `WSTYPE` - the selectors saying which
#: pairs a mixing rule's fitted columns belong to - were sent to `float` and the
#: generator refused the file. Reading both is the same rule applied to the other file,
#: rather than a second list: a column is text because its values are.
TEXT_COLUMNS: frozenset[str] = _text_columns(SOURCES / "COMP.csv") | _text_columns(
    SOURCES / "INTER.csv"
)


def _converter(column: str) -> Callable[[str], object]:
    """The function that turns one upstream value into the unit the manifest states.

    `CONVERSIONS` is the record of every column where the two units differ, so a column
    absent from it is carried as NeqSim stores it and the manifest's `unit` says what
    that is - a stated unit where it is known, `neqsim-internal` where it is not. The
    two are told apart in the manifest rather than here, because the difference between
    "dimensionless" and "nobody has established this" is what a reader needs and a
    conversion table cannot carry it.

    **`CONVERSIONS` is asked first, and the order is load-bearing.**
    `KIJWhitsonSoriede` is the one column in both tables: six of its rows are written
    with a comma decimal separator, which is what makes it look like text, and the
    conversion is what repairs them. Asking `TEXT_COLUMNS` first would carry the comma
    through as a string and quietly lose the fix.
    """
    if column in CONVERSIONS:
        return CONVERSIONS[column]
    if column in TEXT_COLUMNS:
        # `NAME` is lower-cased as well as trimmed: it is the key every other table and
        # the keycard join on, and `INTER.csv` holds the same names in lower case. A
        # text column carried with `str.strip` alone matches nothing.
        return (lambda v: v.strip().lower()) if column == "NAME" else str.strip
    return float


def _component_columns() -> tuple[tuple[str, str, Callable[[str], object]], ...]:
    """The columns to carry, **read from the manifest** rather than listed again.

    This used to be a tuple written out here, next to a manifest that listed the same
    columns, held together by a disagreement check. The check was the symptom: two
    lists that had to be kept equal by hand. There is one list now, and the manifest is
    it - so a column is carried because the manifest says so, and no edit can make the
    two disagree.
    """
    found, problems = manifest_module.read()
    if problems:
        raise SystemExit("gen_databank: " + "\n  ".join(problems))
    return tuple(
        (column.as_field or "", column.name, _converter(column.name))
        for column in found.file("neqsim/COMP.csv").columns
        if column.disposition in manifest_module.CARRIED
    )


COMPONENT_COLUMNS: tuple[tuple[str, str, Callable[[str], object]], ...] = _component_columns()

#: Which of NeqSim's component types a cubic can describe, and therefore which rows of
#: `COMP.csv` the compiled table keeps.
#:
#: Defined in `tools/manifest.py` and imported here, because the `empty-upstream` claims
#: in the manifest are about exactly these rows: the checker that decides them and the
#: generator that selects them have to be answering the same question.
from manifest import KEEP_TYPES  # noqa: E402


#: The interaction parameters, keyed by ordered pair. `KIJPR` is Peng-Robinson's,
#: which is the only cubic this library implements today.
def _inter_columns() -> tuple[tuple[str, str, Callable[[str], object]], ...]:
    """The interaction columns to carry, read from the manifest for the same reason
    `_component_columns` is: one list, and the manifest is it."""
    found, problems = manifest_module.read()
    if problems:
        raise SystemExit("gen_databank: " + "\n  ".join(problems))
    return tuple(
        (column.as_field or "", column.name, _converter(column.name))
        for column in found.file("neqsim/INTER.csv").columns
        if column.disposition in manifest_module.CARRIED
    )


INTER_COLUMNS: tuple[tuple[str, str, Callable[[str], object]], ...] = _inter_columns()

#: The header of the interaction file, derived from the column spec above so the two
#: cannot disagree about the order.
KIJ_HEADER: tuple[str, ...] = tuple(name for name, _, _ in INTER_COLUMNS)

#: The header of the component file, derived from the column spec above so the two
#: cannot disagree about the order.
COMPONENT_HEADER = (*(name for name, _, _ in COMPONENT_COLUMNS), "citation")

#: Values known independently of NeqSim, which is the point: if the conversions above
#: are wrong, these catch it rather than agreeing with the file they were read from.
#: Methane's critical point is a published constant, not a NeqSim output.
ROUND_TRIP = {
    "methane": {
        "tc_k": 190.564,
        "pc_pa": 4.5992e6,
        "acentric_factor": 0.0115,
        "molar_mass_kg_per_mol": 0.016043,
    },
}


def resource_dir(checkout: Path) -> Path:
    """Where NeqSim's data resources actually live.

    Its loader asks the classpath for `data/COMP.csv`. On disk that is
    `src/main/resources/data/` in a source checkout and `data/` in a packaged one,
    so both are tried rather than one assumed. The path itself is tried too, so the
    vendored copies under `databank/sources/neqsim/` - which hold the two files
    directly - work as an argument.
    """
    for candidate in (
        checkout / "src" / "main" / "resources" / "data",
        checkout / "data",
        checkout,
    ):
        if (candidate / "COMP.csv").is_file():
            return candidate
    raise FileNotFoundError(f"no COMP.csv under {checkout}")


def build_components(source: Path) -> list[dict[str, str]]:
    """`data/components/components.csv` as text rows."""
    out: list[dict[str, str]] = []
    for row in read_rows(source / "COMP.csv"):
        name = (row.get("NAME") or "").strip()
        if not name:
            continue
        if (row.get("COMPTYPE") or "").strip() not in KEEP_TYPES:
            # Not a substance a cubic describes. See KEEP_TYPES.
            continue
        try:
            values = {}
            for out_name, src, convert in COMPONENT_COLUMNS:
                raw = (row.get(src) or "").strip()
                # An empty cell is a column this component has no value for - an
                # acid has no CPA association volume - so it is carried empty. It is
                # not zero, and it is not the same thing as a value that will not
                # parse. Treating the two alike is what dropped `hno3` and `h2so4`
                # from the compiled table without a word.
                values[out_name] = "" if not raw else str(convert(raw))
        except (KeyError, ValueError):
            # A non-empty value NeqSim could not parse either. Skipped rather than
            # shipped as a partial record: a component missing its critical pressure is
            # not a component, and half a row would look like one.
            continue
        values["citation"] = CITATION
        out.append(values)
    return out


def build_kij(source: Path, known: set[str]) -> list[dict[str, str]]:
    """`data/components/kij.csv`: every carried column, for pairs this databank has both of.

    The row filter is the one this file always had - both components present in the
    compiled component table, and a Peng-Robinson `KIJPR` to carry - so the row set is
    unchanged. The columns are now every one the manifest carries rather than the
    interaction parameter alone: a pair table that carried one model's parameter and
    dropped the other twenty-seven would be the slice this change exists to remove.
    """
    out: list[dict[str, str]] = []
    for row in read_rows(source / "INTER.csv"):
        a = (row.get("COMP1") or "").strip().lower()
        b = (row.get("COMP2") or "").strip().lower()
        if a not in known or b not in known:
            continue
        raw = (row.get("KIJPR") or "").strip()
        if not raw:
            continue
        try:
            float(raw)
        except ValueError:
            continue
        carried: dict[str, str] = {}
        for out_name, source_name, convert in INTER_COLUMNS:
            text = (row.get(source_name) or "").strip()
            if out_name in ("component_a", "component_b"):
                carried[out_name] = a if out_name == "component_a" else b
            elif not text:
                carried[out_name] = ""
            else:
                # `repr` for a number and `str` for a string. `repr` is here because it
                # round-trips a float exactly, which `str` does not in general - and it
                # is wrong for a text column, where it writes the quotes into the file:
                # `HVTYPE` came out as `'HV'`, which parses as a float nowhere and
                # matches `"HV"` nowhere either.
                value = convert(text)
                carried[out_name] = repr(value) if isinstance(value, float) else str(value)
        out.append(carried)
    return out


def _unifac_columns(file_id: str) -> tuple[tuple[str, str], ...]:
    """The carried columns of one UNIFAC table, as `(compiled, upstream)` pairs."""
    found, problems = manifest_module.read()
    if problems:
        raise SystemExit("gen_databank: " + "\n  ".join(problems))
    return tuple(
        (column.as_field or "", column.name)
        for column in found.file(file_id).columns
        if column.disposition in manifest_module.CARRIED
    )


def build_unifac(source: Path, file_id: str) -> tuple[tuple[str, ...], list[dict[str, str]]]:
    """One UNIFAC table, re-rendered with the manifest's `as` header, values verbatim."""
    columns = _unifac_columns(file_id)
    header = tuple(as_field for as_field, _ in columns)
    rows: list[dict[str, str]] = []
    with (source).open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            rows.append(
                {as_field: (row.get(upstream) or "").strip() for as_field, upstream in columns}
            )
    return header, rows


#: The compiled ISO 6976 table's columns. **Named for what they are, not for NeqSim's
#: spellings**: `srtb15` is the standard's summation factor at 15 °C, and the two names it
#: would suggest - a compression factor and a square root - are both wrong, since the
#: *square of the sum* is what the compression factor subtracts.
ISO6976_HEADER = (
    "name",
    "molar_mass_kg_per_mol",
    "compression_factor_0c",
    "compression_factor_15c",
    "compression_factor_20c",
    "summation_factor_0c",
    "summation_factor_15c",
    "summation_factor_20c",
    "superior_calorific_value_j_per_mol_0c",
    "superior_calorific_value_j_per_mol_15c",
    "superior_calorific_value_j_per_mol_20c",
    "superior_calorific_value_j_per_mol_25c",
    "superior_calorific_value_j_per_mol_60f",
    "inferior_calorific_value_j_per_mol_0c",
    "inferior_calorific_value_j_per_mol_15c",
    "inferior_calorific_value_j_per_mol_20c",
    "inferior_calorific_value_j_per_mol_25c",
    "inferior_calorific_value_j_per_mol_60f",
    "carbon_count",
    "citation",
)

#: NeqSim column -> compiled column, for the fields that cross unchanged in value.
ISO6976_PASSTHROUGH = {
    "Z0": "compression_factor_0c",
    "Z15": "compression_factor_15c",
    "Z20": "compression_factor_20c",
    "srtb0": "summation_factor_0c",
    "srtb15": "summation_factor_15c",
    "srtb20": "summation_factor_20c",
}

#: The calorific values, whose NeqSim columns are in kJ/mol and whose compiled columns are
#: in J/mol - the unit `azoth-core` carries, and the one the specification states.
ISO6976_CALORIFIC = {
    "Hsupmolar0": "superior_calorific_value_j_per_mol_0c",
    "Hsupmolar15": "superior_calorific_value_j_per_mol_15c",
    "Hsupmolar20": "superior_calorific_value_j_per_mol_20c",
    "Hsupmolar25": "superior_calorific_value_j_per_mol_25c",
    "Hsupmolar60F": "superior_calorific_value_j_per_mol_60f",
    "Hinfmolar0": "inferior_calorific_value_j_per_mol_0c",
    "Hinfmolar15": "inferior_calorific_value_j_per_mol_15c",
    "Hinfmolar20": "inferior_calorific_value_j_per_mol_20c",
    "Hinfmolar25": "inferior_calorific_value_j_per_mol_25c",
    "Hinfmolar60F": "inferior_calorific_value_j_per_mol_60f",
}


def build_iso6976(
    source: Path, components: set[str]
) -> tuple[tuple[str, ...], list[dict[str, str]]]:
    """The standard's gas-quality constants, for the components the databank carries.

    **The two unit conversions are NeqSim's own, read out of its use rather than guessed.**
    `MolarMass` crosses as g/mol: `Standard_ISO6976`'s `relDensIdeal` divides it by
    `molarMassAir`, `28.96546` g/mol, and the class's mass-basis branch divides `Mmix` by
    `1000` on its way to a value per kilogram - so the file is g/mol and the compiled column
    is kg/mol, as the component table's is. `Hsupmolar*` and `Hinfmolar*` cross as kJ/mol:
    methane's `Hsupmolar25` is `890.630`, which is the standard's `890.63` kJ/mol, and
    `Stream.LCV()` multiplies the returned value by `1.0e3` to reach joules.

    **Rows are the components the databank carries, matched by name.** The file's 56 rows
    include sixteen the component table does not have - propylene, acetylene, carbon
    monoxide and the rest - and a composition is built from databank names, so those rows
    are unreachable rather than merely unused.
    """
    rows: list[dict[str, str]] = []
    skipped: list[str] = []
    reader = csv.DictReader(source.open(encoding="utf-8"))
    for row in reader:
        name = (row["ComponentName"] or "").strip().lower()
        if name not in components:
            skipped.append(row["ComponentName"])
            continue
        values = {"name": name, "citation": ISO6976_CITATION}
        for upstream, compiled in ISO6976_PASSTHROUGH.items():
            values[compiled] = repr(float(row[upstream]))
        for upstream, compiled in ISO6976_CALORIFIC.items():
            values[compiled] = repr(float(row[upstream]) * 1000.0)
        values["molar_mass_kg_per_mol"] = repr(float(row["MolarMass"]) / 1000.0)
        values["carbon_count"] = str(int(row["numberOfCarbon"]))
        rows.append(values)
    if not rows:
        raise ValueError(f"{source.name}: no row matched a component the databank carries")
    rows.sort(key=lambda row: row["name"])
    return ISO6976_HEADER, rows


#: The compiled columns that carry an element count or charge, and the twelve that carry
#: the order-5 Gibbs and enthalpy polynomials. Both are slices of a header rather than
#: second lists, so they cannot drift from it.
GIBBS_ELEMENT_COLUMNS: tuple[str, ...] = GIBBS_REACTOR_HEADER[1:9]
GIBBS_POLY_COLUMNS: tuple[str, ...] = GIBBS_COEFF_HEADER[1:13]

#: Species whose formation properties are known independently of NeqSim, for the check
#: below. Values are CODATA's standard formation Gibbs energy and enthalpy at 298.15 K in
#: kJ/mol and the standard molar entropy in J/(mol·K), from the same tables the file's own
#: column headings name. Agreement is to the file's three significant figures, because the
#: file is a three-significant-figure table.
GIBBS_ROUND_TRIP: dict[str, dict[str, float]] = {
    "CO2": {"gf298": -394.4, "hf298": -393.5, "sf298": 213.8},
    "water": {"gf298": -228.6, "hf298": -241.8, "sf298": 188.8},
    "methane": {"gf298": -50.5, "hf298": -74.9, "sf298": 186.3},
    "ammonia": {"gf298": -16.4, "hf298": -45.9, "sf298": 192.8},
}


def _gibbs_rows(path: Path) -> list[list[str]]:
    """One of the two reactor files, read as `GibbsReactor.loadGibbsDatabase` reads it.

    Semicolon-separated with comma decimals - `-393,51` - and the class parses every field
    as `Double.parseDouble(parts[i].trim().replace(",", "."))`. Both files carry a UTF-8
    BOM, which the class never sees because it reads through a `Scanner`; decoded as
    `utf-8-sig` it does not become part of the first column's name.

    A row too short for the caller is **refused rather than skipped**. The class logs a
    warning and drops it, which for a shipped table is a species silently missing.
    """
    lines = path.read_text(encoding="utf-8-sig").splitlines()
    return [
        [field.strip().replace(",", ".") for field in line.split(";")]
        for line in lines[1:]
        if line.strip() and not line.strip().startswith("#")
    ]


def build_gibbs_reactor(
    directory: Path, components: set[str]
) -> tuple[tuple[str, ...], list[dict[str, str]]]:
    """`GibbsReactDatabase.csv`: the element vectors and the 298 K formation properties.

    **Rows are the species a fluid can name.** The lookup is
    `componentMap.get(system.getComponent(i).getComponentName().toLowerCase())`, so a row
    whose name is not a component azoth's table carries is unreachable. Of the file's 79
    rows, 49 name such a component; the other 30 are dispositioned in the manifest rather
    than compiled.

    **The first eight columns are counts and the ninth is a charge.** They are read
    positionally from a `;`-split line, so a non-whole value means the offsets are wrong -
    the failure no downstream number would show. Seven of the eight are counts of atoms and
    cannot be negative; the eighth is a charge, which is whole on every row and negative on
    the anions (`OH-` is -1, `SO4--` is -2).

    **The `A`..`D` columns are not compiled.** The class parses them, stores them, clones
    them in a getter and reads them nowhere: the fallback Gibbs branch calls
    `calculateCorrectedHeatCapacityCoeffs`, which reads `system.getComponent(i).getCpA()`
    from NeqSim's own component database. Four columns nothing reads are four the manifest
    dispositions as `not-a-value` rather than four this file copies.
    """
    rows: list[dict[str, str]] = []
    for parts in _gibbs_rows(directory / GIBBS_REACTOR[0]):
        if len(parts) < 16:
            raise ValueError(
                f"{GIBBS_REACTOR[0]}: the row for {parts[0]!r} has {len(parts)} columns, "
                f"where the loader reads at least 16 - a molecule, eight element counts, "
                f"four heat-capacity coefficients and three formation properties"
            )
        name = parts[0]
        if name.lower() not in components:
            continue
        elements = [float(value) for value in parts[1:9]]
        for column, value in zip(GIBBS_ELEMENT_COLUMNS, elements, strict=True):
            floor = -8 if column == "z" else 0
            if value < floor or value != int(value):
                raise ValueError(
                    f"{GIBBS_REACTOR[0]}: {name!r} has {value!r} for {column}, and this "
                    f"column is a whole number no smaller than {floor}"
                )
        values = {
            "name": name,
            "hf298": repr(float(parts[13])),
            "gf298": repr(float(parts[14])),
            "sf298": repr(float(parts[15])),
            "citation": GIBBS_REACTOR_CITATION,
        }
        for column, value in zip(GIBBS_ELEMENT_COLUMNS, elements, strict=True):
            values[column] = str(int(value))
        rows.append(values)

    if not rows:
        raise ValueError(f"{GIBBS_REACTOR[0]}: no row named a component azoth carries")
    return GIBBS_REACTOR_HEADER, rows


def build_gibbs_coeffs(
    directory: Path, species: set[str]
) -> tuple[tuple[str, ...], list[dict[str, str]]]:
    """`DatabaseGibbsFreeEnergyCoeff.csv`: the order-5 Gibbs and enthalpy polynomials.

    **A separate file, joined at lookup time**, because that is what `loadGibbsDatabase`
    does: `componentMap` and `extraCoeffMap` are built independently and joined by molecule
    name, and a species with no row here takes the class's fallback branch. Joining them at
    compile time would move that fallback out of the reader.

    **Rows are the species the reactor table compiled**, not the file's own 18: the three it
    omits - `HNO2`, `NO` and `SO3` - are species azoth's component table does not carry, so
    no fluid can name them and the join can never reach their coefficients.

    **The join is on the lowercased name, and this lowercases both sides** because the class
    does: `componentMap` is keyed by `molecule.toLowerCase()` and `extraCoeffMap` by
    `parts[0].trim().toLowerCase()`. The two files do not agree on case - this one writes
    `co2` where that one writes `CO2` - so a case-sensitive join keeps 8 of the 15 rows and
    says nothing.
    """
    wanted = {name.lower() for name in species}
    rows: list[dict[str, str]] = []
    for parts in _gibbs_rows(directory / GIBBS_COEFFS[0]):
        if len(parts) != 13:
            raise ValueError(
                f"{GIBBS_COEFFS[0]}: a row has {len(parts)} columns where the loader "
                f"reads 13 - a molecule and six Gibbs and six enthalpy coefficients"
            )
        name = parts[0]
        if name.lower() not in wanted:
            continue
        values = {"name": name, "citation": GIBBS_REACTOR_CITATION}
        for index, column in enumerate(GIBBS_POLY_COLUMNS):
            values[column] = repr(float(parts[index + 1]))
        rows.append(values)
    if not rows:
        raise ValueError(f"{GIBBS_COEFFS[0]}: no row named a compiled reactor species")
    return GIBBS_COEFF_HEADER, rows


def check_gibbs_reactor(rows: list[dict[str, str]]) -> list[str]:
    """The formation properties against values known independently of NeqSim.

    **The defence against a uniform misalignment.** Every column here is read positionally
    from a `;`-split line, so swapping two of them moves all 49 rows together and every row
    still looks like a row. The file is a three-significant-figure table, so the check is to
    that precision and no tighter - the class's own polynomials reproduce these same values
    only to a few kJ/mol, which is why they are not the check.
    """
    by_name = {row["name"]: row for row in rows}
    problems: list[str] = []
    for name, expected in GIBBS_ROUND_TRIP.items():
        row = by_name.get(name)
        if row is None:
            problems.append(f"{name}: not in the compiled reactor table")
            continue
        for field, want in expected.items():
            got = float(row[field])
            if abs(got - want) > 0.5:
                problems.append(f"{name}.{field}: got {got!r}, expected about {want!r}")
    return problems


def build_phreeqc(source: Path) -> tuple[tuple[str, ...], list[dict[str, str]]]:
    """The PHREEQC Pitzer catalogue, re-rendered so both kernels read it the same way.

    Reproduces `PhreeqcPitzerParameterCatalog.load`, which is short and worth reading
    beside this:

    * A `#` starts a comment and the rest of the line is dropped; a blank line and the
      `PITZER` banner are skipped.
    * A line beginning `-` opens a family, uppercased, with `LAMDA` accepted for
      `LAMBDA` - PHREEQC's own spelling in older databases. **A family the enum does not
      name closes the section** rather than erroring, so rows under an unknown heading
      are skipped.
    * A row is whitespace-split: the first `speciesCount` tokens are the species and the
      next at most six are the coefficients. **Fewer than six is not an error** - the
      rest are zero, which is how a PHREEQC row with one coefficient compiles.
    * The species are **canonicalised** (`Ba+2` -> `Ba++`, `SO4-2` -> `SO4--`) and then
      sorted, and the key is the sorted names joined with `|`, so a row is found
      whichever order its species are written in.
    * **A duplicate key is an error**, as it is upstream.

    The six coefficients are the PHREEQC temperature form's own, written `a0`..`a5` in
    the compiled header because that is how `PitzerTemperatureFunction` indexes them:

    ```text
    phi(T) = a0 + a1 (1/T - 1/Tr) + a2 ln(T/Tr) + a3 (T - Tr)
                + a4 (T^2 - Tr^2) + a5 (1/T^2 - 1/Tr^2)
    ```

    with `Tr` a reference temperature. **Their units are relative to the family's own** -
    `a1` carries `U K`, `a3` carries `U / K` - so no single unit per column exists, which
    the manifest records rather than guesses. A different polynomial from the four
    coefficients `PitzerParameters.csv` carries.
    """
    # `a0`..`a5`, which is the form's own indexing (`PitzerTemperatureFunction.valueAt`
    # reads `coefficients[0]` as the constant term) rather than the file's 1st..6th. An
    # off-by-one between a table's column names and the formula that reads them is the
    # kind of confusion that survives a review.
    header = ("family", "species_key", *[f"a{i}" for i in range(6)])
    rows: list[dict[str, str]] = []
    seen: set[tuple[str, str]] = set()
    family: str | None = None
    for line in source.read_text(encoding="utf-8").splitlines():
        data = line.split("#", 1)[0].strip()
        if not data or data.upper() == "PITZER":
            continue
        if data.startswith("-"):
            section = data[1:].strip().upper()
            section = "LAMBDA" if section == "LAMDA" else section
            family = section if section in PHREEQC_FAMILIES else None
            continue
        if family is None:
            continue
        tokens = data.split()
        count = PHREEQC_FAMILIES[family]
        if len(tokens) <= count:
            raise ValueError(
                f"{source.name}: `{data}` names {len(tokens)} token(s) and the {family} "
                f"family needs {count} species and at least one coefficient"
            )
        species = [canonical_species(token) for token in tokens[:count]]
        coefficients = [0.0] * 6
        for index, token in enumerate(tokens[count : count + 6]):
            coefficients[index] = float(token)
        # Folded, because the catalogue spells a species as NeqSim's `getComponentName()`
        # does - `Na+`, `HCO3-` - and this library's component table spells it lower case.
        # The lookup is a resolution between two namespaces, and case is the only thing
        # that differs between them.
        key = "|".join(sorted(name.lower() for name in species))
        if (family, key) in seen:
            raise ValueError(
                f"{source.name}: duplicate {family} row for {key}, which upstream refuses"
            )
        seen.add((family, key))
        rows.append(
            {
                "family": family,
                "species_key": key,
                **{f"a{i}": repr(value) for i, value in enumerate(coefficients)},
            }
        )
    return header, rows


def build_furst(source: Path) -> tuple[tuple[str, ...], list[dict[str, str]]]:
    """`FurstElectrolyteConstants`' arrays, re-rendered so a reader can find one.

    **Long form - `set`, `index`, `value` - because the arrays are four different
    lengths.** `furstParams` carries six coefficients, the electrolyte-CPA sets ten,
    `furstParamsCPA_TDep` sixteen and `furstParamsGasIon` eighteen, so a fixed
    `param1..param6` header would have to pad or truncate, and either is a silent lie about
    what the source holds.

    Comments are stripped before the match, because the class keeps every superseded fit as
    a commented-out array beside the live one - eight of them for `furstParams` alone. A
    parser that matched those would compile a parameter set NeqSim stopped using, and the
    numbers would still look plausible.

    A duplicate set name is refused, as a missing value is. The class's own `furstParams`
    is a `public static` field that `setFurstParams` reassigns, so a file that declared it
    twice would leave the compiled table's meaning to whichever line came last.
    """
    header = ("set", "index", "value")
    text = source.read_text(encoding="utf-8")
    # Block comments first, then line comments: a `//` inside a `/* */` block would
    # otherwise end the strip early and expose an array the block was hiding.
    text = re.sub(r"/\*.*?\*/", "", text, flags=re.DOTALL)
    text = re.sub(r"//[^\n]*", "", text)

    rows: list[dict[str, str]] = []
    seen: set[str] = set()
    for match in re.finditer(
        r"public static (?:final )?double\[\]\s+(\w+)\s*=\s*\{(.*?)\};", text, re.DOTALL
    ):
        name, body = match.group(1), match.group(2)
        if name in seen:
            raise ValueError(f"{source.name}: `{name}` is declared twice")
        seen.add(name)
        values = [token.strip() for token in body.split(",")]
        values = [value for value in values if value]
        if not values:
            raise ValueError(f"{source.name}: `{name}` is an empty array")
        for index, value in enumerate(values):
            try:
                number = float(value)
            except ValueError:
                raise ValueError(
                    f"{source.name}: `{name}[{index}]` is {value!r}, which is not a number"
                ) from None
            if not math.isfinite(number):
                raise ValueError(f"{source.name}: `{name}[{index}]` is {value!r}")
            rows.append({"set": name, "index": str(index), "value": repr(number)})
    if not rows:
        raise ValueError(f"{source.name}: no `public static double[]` array was found")
    return header, rows


def canonical_species(name: str) -> str:
    """One species name in NeqSim's spelling: `PhreeqcPitzerParameterCatalog`.

    A numeric charge suffix is rewritten to NeqSim's repeated-sign form and nothing else
    changes - `Ba+2` to `Ba++`, `SO4-2` to `SO4--`, `B(OH)4-` left alone. So this is a
    rewrite and not a table, which is what makes it safe: a species the catalogue carries
    and NeqSim's databank does not still compiles, and the model refuses it later.
    """
    for suffix, sign in (("+3", "+++"), ("-3", "---"), ("+2", "++"), ("-2", "--")):
        if name.endswith(suffix):
            return name[: -len(suffix)] + sign
    return name


def render(header: tuple[str, ...], rows: list[dict[str, str]]) -> str:
    """A file's whole contents, as it will be written."""
    buffer = io.StringIO()
    writer = csv.DictWriter(buffer, fieldnames=list(header), lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    return buffer.getvalue()


def check_round_trip(rows: list[dict[str, str]]) -> list[str]:
    """Assert the conversions against values known independently of NeqSim.

    This is the only defence against the failure this file is most likely to have: a
    conversion that is off by a factor, applied uniformly, producing 258 rows that
    all agree with each other and none of which are right.
    """
    by_name = {row["name"]: row for row in rows}
    problems: list[str] = []
    for name, expected in ROUND_TRIP.items():
        row = by_name.get(name)
        if row is None:
            problems.append(f"{name}: not in the generated databank")
            continue
        for field, want in expected.items():
            got = float(row[field])
            if abs(got - want) / max(abs(want), 1e-300) > 1e-4:
                problems.append(f"{name}.{field}: got {got!r}, expected {want!r}")
    return problems


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "neqsim",
        type=Path,
        nargs="?",
        default=SOURCES,
        help=(
            "path to a NeqSim checkout, or to a directory holding COMP.csv and "
            f"INTER.csv directly. Defaults to {SOURCES.relative_to(ROOT)}, the vendored "
            f"copies, which is what makes --check runnable without a checkout."
        ),
    )
    parser.add_argument("--check", action="store_true", help="report drift; write nothing")
    args = parser.parse_args(argv)

    try:
        resources = resource_dir(args.neqsim)
    except FileNotFoundError:
        print(
            f"gen_databank: no COMP.csv under {args.neqsim}. Point this at a NeqSim "
            f"checkout, or at {SOURCES.relative_to(ROOT)}.",
            file=sys.stderr,
        )
        return 1

    try:
        components = build_components(resources)
        # The catalogue is a `.dat` and its own parser, so its refusals - a short row, a
        # duplicate key - arrive here rather than from the shared path.
        phreeqc = build_phreeqc(resources / PHREEQC_CATALOG[1])
        # Java rather than data, so its own parser too. See `FURST_CONSTANTS`.
        furst = build_furst(SOURCES / FURST_CONSTANTS[1])
    except ValueError as error:
        print(f"gen_databank: {error}", file=sys.stderr)
        return 1

    problems = check_round_trip(components)
    if problems:
        print(
            "gen_databank: the conversions do not reproduce independently known "
            "values. This is a bug in COMPONENT_COLUMNS, not in NeqSim:",
            file=sys.stderr,
        )
        for problem in problems:
            print(f"  {problem}", file=sys.stderr)
        return 1

    kij = build_kij(resources, {row["name"] for row in components})
    iso6976 = build_iso6976(SOURCES / ISO6976_CONSTANTS[0], {row["name"] for row in components})
    reactor = build_gibbs_reactor(SOURCES, {row["name"] for row in components})
    coeffs = build_gibbs_coeffs(SOURCES, {row["name"] for row in reactor[1]})

    problems = check_gibbs_reactor(reactor[1])
    if problems:
        print(
            "gen_databank: the compiled reactor table does not reproduce values known "
            "independently of NeqSim. This is a bug in build_gibbs_reactor's column "
            "offsets, not in NeqSim:",
            file=sys.stderr,
        )
        for problem in problems:
            print(f"  {problem}", file=sys.stderr)
        return 1

    outputs = [
        (OUT_DIR / "components.csv", COMPONENT_HEADER, components),
        (OUT_DIR / "kij.csv", KIJ_HEADER, kij),
        *[
            (OUT_DIR / compiled, *build_unifac(resources / source, file_id))
            for file_id, source, compiled in VENDORED_FILES
        ],
        (OUT_DIR / PHREEQC_CATALOG[2], *phreeqc),
        (OUT_DIR / FURST_CONSTANTS[2], *furst),
        (STANDARDS_DIR / ISO6976_CONSTANTS[1], *iso6976),
        (REACTORS_DIR / GIBBS_REACTOR[1], *reactor),
        (REACTORS_DIR / GIBBS_COEFFS[1], *coeffs),
    ]

    if args.check:
        stale = [
            path
            for path, header, rows in outputs
            if (path.read_text(encoding="utf-8") if path.is_file() else "") != render(header, rows)
        ]
        for path in stale:
            print(f"gen_databank: {path.relative_to(ROOT)} is out of date")
        if stale:
            return 1
        print(
            f"gen_databank: {len(components)} component(s), {len(kij)} kij row(s), "
            f"{len(phreeqc[1])} PHREEQC Pitzer row(s), {len(furst[1])} Furst "
            f"constant(s), {len(reactor[1])} Gibbs reactor species and "
            f"{len(coeffs[1])} coefficient row(s) up to date"
        )
        return 0

    for path, header, rows in outputs:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(render(header, rows), encoding="utf-8")

    print(
        f"gen_databank: {len(components)} component(s) and {len(kij)} kij row(s) "
        f"from NeqSim {NEQSIM_VERSION}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
