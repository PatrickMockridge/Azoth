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
CITATION = (
    f"NeqSim v{NEQSIM_VERSION} COMP.csv (Equinor/NTNU), Apache-2.0, "
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

    outputs = [
        (OUT_DIR / "components.csv", COMPONENT_HEADER, components),
        (OUT_DIR / "kij.csv", KIJ_HEADER, kij),
        *[
            (OUT_DIR / compiled, *build_unifac(resources / source, file_id))
            for file_id, source, compiled in VENDORED_FILES
        ],
        (OUT_DIR / PHREEQC_CATALOG[2], *phreeqc),
        (OUT_DIR / FURST_CONSTANTS[2], *furst),
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
            f"constant(s) up to date"
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
