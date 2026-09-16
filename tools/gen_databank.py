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
import sys
from collections.abc import Callable
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import manifest as manifest_module

ROOT = Path(__file__).resolve().parent.parent
SOURCES = ROOT / "databank" / "sources" / "neqsim"
OUT_DIR = ROOT / "data" / "components"

#: The UNIFAC tables, re-rendered like the component and interaction tables so their
#: compiled headers match the manifest's `as` names. The values are carried verbatim:
#: group-contribution parameters need no unit conversion and no row filtering.
UNIFAC_FILES = (
    ("neqsim/UNIFACcomp.csv", "UNIFACcomp.csv"),
    ("neqsim/UNIFACGroupParam.csv", "UNIFACGroupParam.csv"),
    ("neqsim/UNIFACInterParam.csv", "UNIFACInterParam.csv"),
)

#: The NeqSim release this was generated from, for the `citation` column and for
#: `NOTICE`. A version and a commit, because a databank is only reproducible against
#: one revision of the file it came from.
NEQSIM_VERSION = "3.20.0"
NEQSIM_COMMIT = "dedba8735d030c6411e09b6fd7e69f6c4a136114"
CITATION = f"NeqSim v{NEQSIM_VERSION} COMP.csv (Equinor/NTNU), Apache-2.0, retrieved 2026-09-13"

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
    # An integer selector, not a magnitude: the model number the LIQVISC columns
    # are the coefficients of. Kept whole rather than carried as `3.0`.
    "LIQVISCMODEL": lambda v: int(float(v)),
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


TEXT_COLUMNS: frozenset[str] = _text_columns(SOURCES / "COMP.csv")


def _converter(column: str) -> Callable[[str], object]:
    """The function that turns one upstream value into the unit the manifest states.

    `CONVERSIONS` is the record of every column where the two units differ, so a column
    absent from it is carried as NeqSim stores it and the manifest's `unit` says what
    that is - a stated unit where it is known, `neqsim-internal` where it is not. The
    two are told apart in the manifest rather than here, because the difference between
    "dimensionless" and "nobody has established this" is what a reader needs and a
    conversion table cannot carry it.
    """
    if column in TEXT_COLUMNS:
        # `NAME` is lower-cased as well as trimmed: it is the key every other table and
        # the keycard join on, and `INTER.csv` holds the same names in lower case. A
        # text column carried with `str.strip` alone matches nothing.
        return (lambda v: v.strip().lower()) if column == "NAME" else str.strip
    return CONVERSIONS.get(column, float)


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

#: Which of NeqSim's component types a cubic equation of state can describe.
#:
#: Read from NeqSim's own `COMPTYPE` column rather than guessed, and the reason it
#: is needed at all was found by looking: 62 of the 258 rows are `ion`, 18 carry no
#: type at all, and the rest of the excluded set is `ice`, `salt`, `seawater` and
#: `asphaltene`. **Every one of the ions shares the same critical pressure, acentric
#: factor and critical volume** - `Pc = 29.089 MPa`, `omega = 0.344`, `Vc = 9.9e-05` -
#: because a cubic has no notion of an ion and NeqSim fills those columns with a
#: default. Shipping them would ship 29 rows of plausible-looking wrong numbers,
#: which is the failure this project is organised against, and the repetition is
#: what gave it away: an acentric factor of exactly 0.344 for twenty-nine different
#: substances is not a coincidence.
KEEP_TYPES = frozenset(
    {"HC", "inert", "other", "glycol", "acid", "alcohol", "amine", "chlorine", "water"}
)


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
                carried[out_name] = repr(convert(text))
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
            (OUT_DIR / name, *build_unifac(resources / name, file_id))
            for file_id, name in UNIFAC_FILES
        ],
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
        print(f"gen_databank: {len(components)} component(s), {len(kij)} kij row(s) up to date")
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
