#!/usr/bin/env python3
"""Generate `data/components/` from a NeqSim checkout.

    python tools/gen_databank.py /tmp/neqsim-check/neqsim
    python tools/gen_databank.py /tmp/neqsim-check/neqsim --check

# What this is

azoth ships no component data of its own, and every calculation in the `eos`
namespace takes `Tc`, `Pc` and `omega` as caller arguments because of it. NeqSim is
Apache-2.0 and has the databank; this compiles a slice of it into the canonical CSVs
both implementations read, exactly as `gen_registry.py` compiles specs and
`gen_user_data.py` compiles a user's file.

`NOTICE` at the repository root carries the attribution. It is not repeated in the
generated file's rows beyond the `citation` column, which names the source once per
row because that is the column the loader reads.

# The units, which are the whole risk

Read out of NeqSim's own loader - `thermo/component/Component.java` around line 305
- rather than guessed from the column names:

    molarmass  divided by 1000        -> kg/mol      (the file is g/mol)
    TC         plus 273.15            -> K           (the file is degrees Celsius)
    PC         unchanged              -> bar         (see the check below)
    critvol    unchanged              -> cm**3/mol
    liqdens    unchanged              -> g/cm**3

Every one of those is applied once, here, and the result is asserted against values
known independently of NeqSim. Guessing any of them would be the `mm` bug again: a
factor that is silently wrong and looks entirely reasonable.

`PC` is the one read out rather than copied from a comment. NeqSim's
`criticalCompressibilityFactor` is `Pc * Vc / R / Tc / 10`, and methane's
`Z_c = 0.2874` falls out of that arithmetic with `Pc` in bar and `Vc` in cm**3/mol -
which is what fixes those two and rules out Pa, kPa and m**3/mol.

# What is kept, and what is not

Only the columns the calculations use: identity, critical properties, acentric
factor, and the binary interaction parameters. `COMP.csv` has 170 columns and most
of them are hydrate, wax, electrolyte and SAFT data for physics this library does not
have; shipping them would be shipping data nothing reads. The association columns
the CPA and PC-SAFT models will need are a later addition, and they are left out
until those models exist and their units can be read the same way these were.

`COMP_EXT.csv` (86 MB, 76,705 rows) is not vendored at all. Nearly all of it is heavy
and characterised fluids, which azoth cannot use: it has no TBP or plus-fraction
characterisation.
"""

from __future__ import annotations

import argparse
import csv
import io
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
OUT_DIR = ROOT / "data" / "components"

#: The NeqSim release this was generated from, for the `citation` column and for
#: `NOTICE`. A version and a commit, because a databank is only reproducible against
#: one revision of the file it came from.
NEQSIM_VERSION = "3.20.0"
NEQSIM_COMMIT = "dedba8735d030c6411e09b6fd7e69f6c4a136114"
CITATION = f"NeqSim v{NEQSIM_VERSION} COMP.csv (Equinor/NTNU), Apache-2.0, retrieved 2026-09-13"

#: Output column -> (NeqSim column, conversion). The conversions are the ones
#: documented above, each with the reason it is that and not another.
COMPONENT_COLUMNS: tuple[tuple[str, str, object], ...] = (
    ("name", "NAME", lambda v: v.strip().lower()),
    ("cas", "CASnumber", lambda v: v.strip()),
    ("formula", "FORMULA", lambda v: v.strip()),
    ("molar_mass_kg_per_mol", "MOLARMASS", lambda v: float(v) / 1000.0),
    ("tc_k", "TC", lambda v: float(v) + 273.15),
    ("pc_pa", "PC", lambda v: float(v) * 1.0e5),
    ("acentric_factor", "ACSFACT", float),
    ("critical_volume_m3_per_mol", "CRITVOL", lambda v: float(v) * 1.0e-6),
    ("liquid_density_kg_per_m3", "LIQDENS", lambda v: float(v) * 1000.0),
)

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
KIJ_HEADER = ("component_a", "component_b", "kij_pr")

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
    so both are tried rather than one assumed.
    """
    for candidate in (checkout / "src" / "main" / "resources" / "data", checkout / "data"):
        if (candidate / "COMP.csv").is_file():
            return candidate
    raise FileNotFoundError(f"no COMP.csv under {checkout}")


def read_rows(path: Path) -> list[dict[str, str]]:
    """Rows of a NeqSim resource CSV."""
    with path.open(newline="", encoding="utf-8", errors="replace") as handle:
        return [row for row in csv.DictReader(handle) if row.get("NAME") or row.get("COMP1")]


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
            values = {out_name: convert(row[src]) for out_name, src, convert in COMPONENT_COLUMNS}
        except (KeyError, ValueError):
            # A row NeqSim cannot parse either. Skipped rather than shipped as a
            # partial record: a component missing its critical pressure is not a
            # component, and half a row would look like one.
            continue
        values["citation"] = CITATION
        out.append(values)
    return out


def build_kij(source: Path, known: set[str]) -> list[dict[str, str]]:
    """`data/components/kij.csv`, restricted to pairs this databank has both of."""
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
            value = float(raw)
        except ValueError:
            continue
        out.append({"component_a": a, "component_b": b, "kij_pr": repr(value)})
    return out


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
    parser.add_argument("neqsim", type=Path, help="path to a NeqSim checkout")
    parser.add_argument("--check", action="store_true", help="report drift; write nothing")
    args = parser.parse_args(argv)

    try:
        resources = resource_dir(args.neqsim)
    except FileNotFoundError:
        print(
            f"gen_databank: no COMP.csv under {args.neqsim}. Point this at a NeqSim "
            f"checkout - the repository does not vendor one.",
            file=sys.stderr,
        )
        return 1

    components = build_components(resources)
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
