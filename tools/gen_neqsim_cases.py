#!/usr/bin/env python3
"""Generate the recorded NeqSim validation cases from the committed probe captures.

    python tools/gen_neqsim_cases.py [--check]

A **capture** is a probe's standard output, committed verbatim under
``validation/neqsim/captures/``. The probes are run by hand, against the pinned
``neqsim-f0c7436.jar``, and the capture is what they printed. This reads a capture and emits
``validation/eos/<name>_against_neqsim.json``, the shape
``python/tests/validation/test_validation_cases.py`` already consumes.

**It reads the capture rather than running the JVM**, which is why ``--check`` needs neither
a JDK nor the jar - the jar is gitignored, so a gate that needed it could not run in CI.
The consequence is that a capture is a pinned artifact: regenerating one means running the
probe again and committing its output, and the case moves with it.

# The table is the provenance

``CASES`` names, for each case, the capture, the state in it, and **the capture key each
azoth output is taken from**. That is not decoration. ``CpaSweep``'s own header carries the
same table for the same reason: a quantity whose provenance is not written down beside it
gets misread, and this project has already spent a session on an ``a`` that was read as a
pure component's attraction when it was an interaction row sum.

A case's ``derivation`` is generated from the capture - the command, the row, and the
numbers the probe printed - so what the case says the oracle reported is what the oracle
reported. The ``source`` block is authored, and says why the case is worth having and what
it does not cover.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from dataclasses import dataclass
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
CAPTURES = ROOT / "validation" / "neqsim" / "captures"
OUT = ROOT / "validation" / "eos"

#: `key = value`, where the key is everything since the last assignment on the line and the
#: value is either a bracketed vector or one bare token. `phase[0] lnPhi = [1, 2]` and
#: `aHS=0.1166 gHS=1.0740` are both one match per key; a component line's leading prose
#: (`component[0] methane m=1.0`) ends up keyed by the whole run before the `=`, which is
#: why no case maps a key that only ever appears there.
_PAIR = re.compile(r"([^\s=][^=]*?)\s*=\s*(\[[^\]]*\]|\S+)")


def pairs(text: str, first: bool = False) -> dict[str, str]:
    """Every `key = value` in a run of text.

    **A repeated key keeps the last by default and the first on request**, and which one is
    right is a property of the capture rather than a preference. `CpaSweep`'s `case` rows
    and `PcsaftProbe`'s blocks are flat, so it never comes up. `UmrCpaProbe` prints the
    single-phase state's `Z` and `lnPhi` and then the two flashed phases' own under the same
    bare names, and it is the first that the case is about; `SaftVrMieFlashProbe` prints the
    same figure for each of the three ways it runs the class, and it is the last that
    belongs to the initialised one.
    """
    out: dict[str, str] = {}
    for line in text.splitlines():
        stripped = line.strip()
        if not stripped or stripped.startswith("#") or stripped.startswith("--"):
            continue
        if "=" in line:
            found = [(key.strip(), value) for key, value in _PAIR.findall(line)]
        else:
            # `UmrCpaProbe`'s own shape: a bare key, whitespace, a number. Only a line whose
            # second token is a number is taken, so a probe's prose cannot set a key.
            parts = line.split()
            try:
                float(parts[1])
            except (IndexError, ValueError):
                continue
            found = [(parts[0], " ".join(parts[1:]))]
        for key, value in found:
            if first and key in out:
                continue
            out[key] = value
    return out


def vector(value: str) -> list[float]:
    """A capture value as a list: a bracketed one entry by entry, a bare one as itself.

    **Both shapes occur for the same output.** `ln_phi` is one component per capture key in
    `CpaSweep` and a bracketed composition in `SaftVrMieFlashProbe`, and the case says which
    entry of the *output* vector each mapping fills rather than how the capture spells it.
    """
    stripped = value.strip()
    if not stripped.startswith("["):
        return [float(stripped)]
    inner = stripped.removeprefix("[").removesuffix("]")
    return [float(part) for part in inner.split(",")]


# --- the four capture shapes ---------------------------------------------------------


def line_records(text: str) -> list[dict[str, str]]:
    """One record per data line: `CpaSweep`, `PrCpaFlash`."""
    out = []
    for line in text.splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        out.append(pairs(line))
    return out


def block_records(text: str, first: bool = False) -> list[dict[str, str]]:
    """One record per state block: `PcsaftProbe`, `SaftVrMieProbe`, `UmrCpaProbe`.

    A comment line opens a block, and a probe's banner is several comment lines before the
    first state - so the banner opens blocks of its own that carry no keys, and those are
    dropped rather than counted. **Indexing the states by position after that**, rather than
    by their header text: the header text is what a case quotes, and a case that both
    selected and quoted by it would be quoting itself.
    """
    blocks: list[list[str]] = []
    for line in text.splitlines():
        if line.startswith("#"):
            blocks.append([line])
        elif blocks:
            blocks[-1].append(line)
    records = [pairs("\n".join(block), first) for block in blocks]
    return [record for record in records if record]


#: How `SaftVrMieFlashProbe` labels each way it runs the class: the marker on the
#: `phases=...` line, and the marker it puts on the composition lines beneath it.
_FLASH_WAYS = (("no init(0) first", "no-init"), ("init(0) first", "init(0)"))


def flash_records(text: str) -> list[dict[str, str]]:
    """One record per state, from the probe that runs the class three ways.

    `SaftVrMieFlashProbe`'s output is prose - a section per way of running the class, with
    the compositions on lines that carry no key of their own - so the two other shapes do
    not fit it. **The `init(0)`-first answer is the one a case takes**, for the reason
    `~/Desktop/neqsim-tpflashsaft-reads-a-zero-feed.md` records: `run()` reads `getz()`
    before its own `system.init(0)`, so on a fresh system the Rachford-Rice solves a zero
    feed and never finds the split.
    """
    out: list[dict[str, str]] = []
    for block in text.split("\n\n"):
        # `phases=` and not the banner, which also names `TPflashSAFT`.
        if "phases=" not in block:
            continue
        record = pairs(block)
        for reported, tagged in _FLASH_WAYS:
            for phase in ("GAS", "OIL"):
                found = re.search(rf"{re.escape(tagged)} type={phase} x = (\[[^\]]*\])", block)
                if found:
                    record[f"{reported} {phase} x"] = found.group(1)
            # Anchored to the line's own content: `"no init(0) first"` *contains*
            # `"init(0) first"`, so an unanchored search takes the no-init line's vapour
            # fraction for the initialised one's.
            beta = re.search(
                rf"(?:^|\n)(?:TPflashSAFT, )?{re.escape(reported)}\s+phases=\d+\s+beta=(\S+)",
                block,
            )
            if beta:
                record[f"{reported} beta"] = beta.group(1)
        k = re.search(r"^\s*K = (\[[^\]]*\])", block, re.MULTILINE)
        if k:
            record["K"] = k.group(1)
        out.append(record)
    return out


SHAPES = {
    "line": line_records,
    "block": block_records,
    "block-first": lambda text: block_records(text, first=True),
    "flash": flash_records,
}


# --- the cases -----------------------------------------------------------------------


@dataclass(frozen=True)
class Mapping:
    """One azoth output field, and the capture key it is read from."""

    field: str
    key: str
    #: For a vector output, which entry this mapping fills; `None` means a scalar output.
    index: int | None = None
    #: The output field, when the output is a vector and `field` names the same one for
    #: every entry. `None` means the capture key's own name is the output field.
    under: str | None = None


@dataclass(frozen=True)
class Case:
    capture: str
    shape: str
    state: int
    id: str
    calc: str
    tol: float
    #: The probe's command line, as it was run.
    command: str
    #: What the state's own comment header in the capture says, so the derivation can quote
    #: the line rather than search for it.
    header: str
    attribution: str
    notes: str
    #: What the case checks and what it does not, in the author's words.
    note: str
    inputs: dict[str, object]
    eos: str | None
    mappings: tuple[Mapping, ...]


CASES: tuple[Case, ...] = (
    Case(
        capture="cpa_sweep.tsv",
        shape="line",
        state=0,
        id="water_methanol_srk_cpa_against_neqsim",
        calc="eos.srk_cpa_phase",
        tol=1.0e-10,
        command="java -cp .:neqsim-f0c7436.jar CpaSweep 300 100 0.6",
        header="case 0",
        attribution=(
            "NeqSim master, SystemSrkCPA with setMixingRule(10), run from "
            "validation/neqsim/CpaSweep.java"
        ),
        notes=(
            "The single-phase liquid state the probe's `init(1)` leaves, at the state the "
            "SRK-CPA kernel was first verified at. `lnPhi` is "
            "`ComponentInterface.getFugacityCoefficient()` on the same phase."
        ),
        note=(
            "The whole association model in one comparison: the fitted `aCPA`/`bCPA` "
            "substitution, the `cpakij_srk` interaction column, the association's pressure, "
            "the site-fraction solve and the fugacity coefficient built on it."
        ),
        inputs={
            "components": ["water", "methanol"],
            "T": 300.0,
            "P": 10000000.0,
            "z": [0.6, 0.4],
            "compressed_phase": "liquid",
        },
        eos=None,
        mappings=(
            Mapping("z_factor", "Z"),
            Mapping("ln_phi", "lnPhi[0]", index=0, under="ln_phi"),
            Mapping("ln_phi", "lnPhi[1]", index=1, under="ln_phi"),
        ),
    ),
    Case(
        capture="pcsaft_probe.tsv",
        shape="block",
        state=1,
        id="methane_butane_pcsaft_against_neqsim",
        calc="eos.pcsaft_rahmat_phase",
        # `1e-7` and not `1e-8`: see the note. Upstream reworked both PCSAFT phases after
        # 3.20.0, which is the revision this conversion is a port of.
        tol=1.0e-07,
        command="java -cp .:neqsim-f0c7436.jar PcsaftProbe 350 30 methane 0.6 n-butane 0.4",
        header="# methane/n-butane at T=350",
        attribution=(
            "NeqSim master, SystemPCSAFT - which builds PhasePCSAFTRahmat - run from "
            "validation/neqsim/PcsaftProbe.java"
        ),
        notes=(
            "The probe's `volumeSAFT` is the molar volume its own solve converged to, and "
            "its `v` line is the class's `getMolarVolume()`, which returns `Infinity` here "
            "and is not recorded. **The mixture is deliberate**: methane alone has `m = 1`, "
            "so `m_bar - 1` and every term carrying it vanish."
        ),
        note=(
            "This pair's `kij` is `0.022` in the `KIJPCSAFT` column, which no shipped model "
            "read before this one. **The departures are not recorded, and the tolerance is "
            "`1e-7` rather than the `1e-8` its siblings carry to say so**: `PhasePCSAFT."
            "getdDSAFTdT` carried a spurious `3 d_i**2` and its `Hres` was 3% out at 3.20.0, "
            "which is the revision this conversion is a port of. Upstream reworked both PCSAFT "
            "phases afterwards - dropping the override and adding a finite-difference "
            "consistency test - so these numbers are `3.8e-8` from master's where they were "
            "`4.6e-11` from 3.20.0's. **Where the two differ is localised and worth "
            "recording**: the volume and `Z` are `1.2e-9` apart and n-butane's `ln phi` "
            "`2.4e-10`, while **methane's is `3.8e-8`** - the `m = 1` component, where "
            "`m_bar - 1` and every term carrying it vanish. So it is a *derivative* and not "
            "a state, and the two libraries' corrections are not the same correction. "
            "Re-porting is the next tranche's call; see "
            "`~/Desktop/neqsim-pcsaft-hard-chain-temperature-derivative.md`."
        ),
        inputs={
            "components": ["methane", "n-butane"],
            "T": 350.0,
            "P": 3000000.0,
            "z": [0.6, 0.4],
            "compressed_phase": "vapour",
        },
        eos=None,
        mappings=(
            Mapping("z_factor", "Z"),
            Mapping("v", "volumeSAFT"),
            Mapping("ln_phi", "lnPhi[0]", index=0, under="ln_phi"),
            Mapping("ln_phi", "lnPhi[1]", index=1, under="ln_phi"),
        ),
    ),
    Case(
        capture="saft_vr_mie_probe.tsv",
        shape="block",
        state=1,
        id="methane_butane_saft_vr_mie_against_neqsim",
        calc="eos.saft_vr_mie_phase",
        tol=1.0e-04,
        command="java -cp .:neqsim-f0c7436.jar SaftVrMieProbe 350 30 methane 0.6 n-butane 0.4",
        header="# methane/n-butane at T=350",
        attribution=(
            "NeqSim master, SystemSAFTVRMie - which builds PhaseSAFTVRMie - run from "
            "validation/neqsim/SaftVrMieProbe.java"
        ),
        notes=(
            "The probe's `v` line is `getVolume()/n` in the class's internal scale, not a "
            "molar volume, so the volume is not recorded here; the compressibility and the "
            "fugacity coefficients are. **The tolerance is `1e-4`** because this model's "
            "last digits are an arithmetic rather than a model: NeqSim takes the `eta` "
            "derivatives of `gHS` and of the dispersion by central difference at a relative "
            "`1e-5` and the crate reproduces that step, so the two sides are differencing "
            "the same nearly-equal energies."
        ),
        note=(
            "Every layer reproduces one at a time in `crates/azoth-eos/tests/saft_vr_mie.rs` "
            "- the diameter, the chain contact value, `g1`, `g2`, the dispersion, the "
            "Helmholtz energy and the pressure. The departures are not recorded: "
            "`PhaseSAFTVRMie.dF_HC_SAFTdT` leaves out the chain contact value's own "
            "temperature dependence; see "
            "`~/Desktop/neqsim-saft-vr-mie-chain-contact-value-temperature.md`."
        ),
        inputs={
            "components": ["methane", "n-butane"],
            "T": 350.0,
            "P": 3000000.0,
            "z": [0.6, 0.4],
            "compressed_phase": "vapour",
        },
        eos=None,
        mappings=(
            Mapping("z_factor", "Z"),
            Mapping("ln_phi", "lnPhi[0]", index=0, under="ln_phi"),
            Mapping("ln_phi", "lnPhi[1]", index=1, under="ln_phi"),
        ),
    ),
    Case(
        capture="saft_vr_mie_flash.tsv",
        shape="flash",
        state=1,
        id="methane_butane_saft_flash_against_neqsim",
        calc="eos.tp_flash_saft",
        tol=1.0e-04,
        command="java -cp .:neqsim-f0c7436.jar SaftVrMieFlashProbe",
        header="# methane/n-butane at T=250",
        attribution=(
            "NeqSim master, SystemSAFTVRMie with TPflashSAFT, run from "
            "validation/neqsim/SaftVrMieFlashProbe.java"
        ),
        notes=(
            "**The `init(0)`-first answer, and only that one.** `TPflashSAFT.run()` reads "
            "`getz()` before its own `system.init(0)`, so on a fresh system it solves a zero "
            "feed and never finds the split; the probe runs it three ways and they disagree. "
            "`~/Desktop/neqsim-tpflashsaft-reads-a-zero-feed.md` is the write-up, and the "
            "probe's own two other reports are left out of this case rather than "
            "reproduced. The probe does not print the two compressibility factors, so "
            "`z_liquid` and `z_vapour` are not checked. **The tolerance is the loop's own "
            "stopping rule**, not a measurement of the model."
        ),
        note=(
            "The compositions and the vapour fraction are one quantity here, and this is "
            "the tranche's only two-phase case."
        ),
        inputs={
            "components": ["methane", "n-butane"],
            "T": 250.0,
            "P": 3000000.0,
            "z": [0.6, 0.4],
        },
        eos=None,
        mappings=(
            Mapping("beta", "init(0) first beta"),
            Mapping("x", "init(0) first OIL x", index=0, under="x"),
            Mapping("x", "init(0) first OIL x", index=1, under="x"),
            Mapping("y", "init(0) first GAS x", index=0, under="y"),
            Mapping("y", "init(0) first GAS x", index=1, under="y"),
            Mapping("k", "K", index=0, under="k"),
            Mapping("k", "K", index=1, under="k"),
        ),
    ),
    Case(
        capture="umr_cpa_probe.tsv",
        shape="block-first",
        state=4,
        id="methane_water_umr_cpa_against_neqsim",
        calc="eos.umr_cpa_phase",
        tol=1.0e-09,
        command="java -cp .:neqsim-f0c7436.jar UmrCpaProbe",
        header="# methane/water at T=298.150 K, P=70",
        attribution=(
            "NeqSim master, SystemUMRCPAEoS - which builds PhaseUMRCPA - run from "
            "validation/neqsim/UmrCpaProbe.java"
        ),
        notes=(
            "**The state the lifecycle test measures its water-in-gas envelope at**, 298.15 K "
            "and 70 bar, and the last of the probe's nine. `HresTP` and `SresTP` are "
            "NeqSim's own and are recorded: UMR-CPA is a Peng-Robinson with an excess-Gibbs "
            "mixing rule, not a Helmholtz model, so the two defects the PC-SAFT and "
            "SAFT-VR-Mie departures carry do not reach it."
        ),
        note=(
            "This is the case the UMR-CPA port was built against, and it is what caught the "
            "two shared-code defects: the association's departure carrying a spurious "
            "`A/(RT)`, and an excess-Gibbs rule's departure not being the classical "
            "`psi_bar` form."
        ),
        inputs={
            "components": ["methane", "water"],
            "T": 298.15,
            "P": 7000000.0,
            "z": [0.98, 0.02],
            "compressed_phase": "vapour",
        },
        eos=None,
        mappings=(
            Mapping("z_factor", "Z"),
            Mapping("ln_phi", "lnPhi[0]", index=0, under="ln_phi"),
            Mapping("ln_phi", "lnPhi[1]", index=1, under="ln_phi"),
            Mapping("h_res", "HresTP_J_per_mol"),
            Mapping("s_res", "SresTP_J_per_molK"),
        ),
    ),
    Case(
        capture="furst_probe.tsv",
        shape="block",
        state=2,
        id="methane_water_sodium_chloride_furst_against_neqsim",
        calc="eos.furst_electrolyte_phase",
        tol=1.0e-10,
        command="java -cp .:neqsim-f0c7436.jar FurstProbe",
        header="# the shipped test: methane water Na+ Cl- phase 1",
        attribution=(
            "NeqSim master, SystemFurstElectrolyteEos - which builds "
            "PhaseModifiedFurstElectrolyteEos - run from validation/neqsim/FurstProbe.java"
        ),
        notes=(
            "**The aqueous phase of the shipped test, which is the state the model was built "
            "against**, and the layers are why this case exists: the shielding solve, the "
            "Born radius, the solvent dielectric and the short-range table all reach "
            "`ln phi` as one number, and this is what says which of them moved. `Z` and "
            "`lnPhi` are the state `specs/cases/eos/furst_electrolyte_phase.toml` pins."
        ),
        note=(
            "`W = -3.25444241835884e-07` is the value the gas-ion pass gives the methane/Na+ "
            "pair; the predictive correlation alone gives `-1.92e-07`, so the layer is live. "
            "`eps_dT = -0.359218709298880` is the solvent dielectric's temperature "
            "derivative, which the Mod2004 revision zeroes."
        ),
        inputs={
            "components": ["methane", "water", "Na+", "Cl-"],
            "T": 298.15,
            "P": 1001325.0,
            "x": [
                0.000225745660581355,
                0.997778050427449,
                0.000998101955985164,
                0.000998101955985164,
            ],
            "compressed_phase": "liquid",
        },
        eos=None,
        mappings=(
            Mapping("z_factor", "Z"),
            Mapping("ln_phi", "lnPhi[0]", index=0, under="ln_phi"),
            Mapping("ln_phi", "lnPhi[1]", index=1, under="ln_phi"),
            Mapping("ln_phi", "lnPhi[2]", index=2, under="ln_phi"),
            Mapping("ln_phi", "lnPhi[3]", index=3, under="ln_phi"),
        ),
    ),
)


def expected_of(case: Case, record: dict[str, str]) -> dict[str, object]:
    """The case's `expected` block, read out of one capture record."""
    out: dict[str, object] = {}
    vectors: dict[str, list[tuple[int, float]]] = {}
    for mapping in case.mappings:
        if mapping.key not in record:
            raise SystemExit(
                f"gen_neqsim_cases: {case.capture} state {case.state} has no `{mapping.key}`"
            )
        raw = record[mapping.key].strip()
        if mapping.under is None:
            out[mapping.field] = float(raw)
            continue
        # A bracketed capture value carries the whole vector and `index` picks from it; a
        # bare one *is* the entry, and `index` only says where in the output it goes.
        assert mapping.index is not None
        taken = vector(raw)[mapping.index] if raw.startswith("[") else float(raw)
        vectors.setdefault(mapping.under, []).append((mapping.index, taken))
    for field, indexed in vectors.items():
        out[field] = [value for _, value in sorted(indexed)]
    return out


#: A bracketed run of JSON scalars, however `json.dumps` wrapped it.
_ARRAY = re.compile(r"\[\s*([^\[\]{}]*?)\s*\]", re.DOTALL)


def compact(text: str) -> str:
    """Put every array of scalars back on one line, as the rest of the directory has them.

    `json.dumps(indent=2)` gives a composition one number per line, which is nine lines for
    a binary. The other cases in `validation/eos/` carry a vector on the line its key is on,
    and a reader comparing two of them should not have to hold two layouts in mind.
    """
    return _ARRAY.sub(
        lambda m: "[" + ", ".join(p.strip() for p in m.group(1).split(",")) + "]", text
    )


def render_case(case: Case, record: dict[str, str], header: str) -> str:
    """One case's JSON, in the order the directory's other files use."""
    case_json: dict[str, object] = {
        "id": case.id,
        "calc": case.calc,
        "source": {
            "attribution": case.attribution,
            "verification": "verified",
            "notes": case.notes,
        },
        "derivation": (
            f"{case.note}\n\n"
            f"The probe was run as `{case.command}`, reading `{case.capture}`; this is its "
            f"state {case.state}, whose own header line is:\n\n"
            f"    {header.strip()}\n\n"
            f"The case compares against the values that header's block prints, which are "
            f"in the capture verbatim."
        ),
        "inputs": case.inputs,
        "expected": expected_of(case, record),
        "tolerance": case.tol,
    }
    if case.eos is not None:
        case_json["eos"] = case.eos
    return compact(json.dumps(case_json, indent=2)) + "\n"


def render() -> dict[Path, str]:
    """Every case file, keyed by where it goes."""
    out: dict[Path, str] = {}
    for case in CASES:
        text = (CAPTURES / case.capture).read_text(encoding="utf-8")
        records = SHAPES[case.shape](text)
        if case.state >= len(records):
            raise SystemExit(
                f"gen_neqsim_cases: {case.capture} has {len(records)} states, not {case.state + 1}"
            )
        found = next((line for line in text.splitlines() if case.header in line), case.header)
        # The line itself, trimmed: a `line`-shaped capture's record is one enormous row and
        # quoting it whole would bury the derivation in keys the case does not read.
        header = found if len(found) <= 200 else found[:200] + " ..."
        out[OUT / f"{case.id}.json"] = render_case(case, records[case.state], header)
    return out


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--check", action="store_true", help="report drift; write nothing")
    args = parser.parse_args(argv)

    files = render()
    stale = []
    for path, text in sorted(files.items()):
        current = path.read_text(encoding="utf-8") if path.is_file() else ""
        if current != text:
            stale.append(path)
        elif args.check:
            print(f"gen_neqsim_cases: {path.relative_to(ROOT)} up to date")

    if args.check:
        for path in stale:
            print(f"gen_neqsim_cases: {path.relative_to(ROOT)} is out of date", file=sys.stderr)
        if stale:
            return 1
        return 0

    for path in stale:
        path.write_text(files[path], encoding="utf-8")
        print(f"gen_neqsim_cases: wrote {path.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
