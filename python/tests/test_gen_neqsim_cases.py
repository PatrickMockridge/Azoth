"""The recorded NeqSim cases, and the generator that derives them from the captures.

`tools/gen_neqsim_cases.py` reads a probe's committed stdout and writes the
`validation/eos/*.json` a case is. Two things about it can be wrong and neither shows
up as a failing validation test:

* the generated files can drift from the captures they claim to come from, which is what
  `--check` is for and what the drift job runs; and
* the *reading* can be wrong, in which case the case records a number the probe never
  printed - and it passes, because both sides of the comparison moved together.

So this exercises the readers as well as the gate: which of two repeated keys a shape
keeps, and that a probe's prose cannot set a key.
"""

from __future__ import annotations

import importlib
import subprocess
import sys
from pathlib import Path
from types import ModuleType
from typing import Any

REPO_ROOT = Path(__file__).resolve().parents[2]
TOOLS = REPO_ROOT / "tools"


def generator() -> ModuleType:
    """`tools/gen_neqsim_cases.py`, imported by name.

    The tools are not a package and sit outside mypy's `files`, so the import is dynamic -
    the arrangement `test_ci_scripts.py` and `test_prose_lint.py` both use.
    """
    sys.path.insert(0, str(TOOLS))
    try:
        return importlib.import_module("gen_neqsim_cases")
    finally:
        sys.path.pop(0)


def test_the_committed_cases_are_current() -> None:
    """The gate the drift job runs, run here so it fails in the suite too."""
    result = subprocess.run(
        [sys.executable, str(TOOLS / "gen_neqsim_cases.py"), "--check"],
        capture_output=True,
        text=True,
        check=False,
    )

    assert result.returncode == 0, (
        f"the recorded NeqSim cases are out of date - run "
        f"`python tools/gen_neqsim_cases.py`\n{result.stdout}{result.stderr}"
    )


def test_a_stale_capture_is_reported(tmp_path: Path) -> None:
    """`--check` is a comparison against the captures, not against nothing.

    The captures are redirected to a copy with one number moved, and the emit has to
    stop matching the committed file. Without this the gate could pass on a generator
    that returned the committed text unconditionally.
    """
    tool: Any = generator()
    original = tool.CAPTURES
    # The value the first case records for `z_factor`, moved in the last digit.
    moved = "0.105050962879418"
    try:
        tool.CAPTURES = tmp_path
        for case in tool.CASES:
            text = (original / case.capture).read_text(encoding="utf-8")
            if case is tool.CASES[0]:
                text = text.replace(moved, "0.105050962879419", 1)
            (tmp_path / case.capture).write_text(text, encoding="utf-8")
        rendered = tool.render()
    finally:
        tool.CAPTURES = original

    stale = [path for path in rendered if rendered[path] != path.read_text(encoding="utf-8")]
    assert stale, "a capture with a moved number must produce a different case"


def test_every_case_reads_the_keys_it_names() -> None:
    """Every mapping resolves against the state it selects.

    This is the check that a rename in a probe's output cannot pass: the generator raises
    rather than writing a case with the key missing, and a case with a *wrong* key would
    be a number from somewhere else in the capture.
    """
    tool: Any = generator()
    for case in tool.CASES:
        text = (tool.CAPTURES / case.capture).read_text(encoding="utf-8")
        records = tool.SHAPES[case.shape](text)
        assert case.state < len(records), f"{case.id}: {case.capture} has {len(records)} states"
        found = tool.expected_of(case, records[case.state])
        assert found, f"{case.id} recorded nothing"


def test_a_repeated_key_keeps_the_last_unless_asked_for_the_first() -> None:
    """Which occurrence a shape keeps is a property of the capture, not a preference.

    `UmrCpaProbe` prints the single-phase state's `Z` and then the two flashed phases' own
    under the same bare name; the case is about the first. `SaftVrMieFlashProbe` prints the
    same figure for each of the three ways it runs the class, and the initialised one is
    the last. Both are read from the real captures rather than from a fixture, because the
    disagreement between the two answers is the whole reason the distinction exists.
    """
    tool: Any = generator()

    umr = (tool.CAPTURES / "umr_cpa_probe.tsv").read_text(encoding="utf-8")
    first = tool.SHAPES["block-first"](umr)[4]
    last = tool.SHAPES["block"](umr)[4]
    assert first["Z"] != last["Z"], (
        "the two readings must differ, or this test is checking nothing: "
        f"{first['Z']} against {last['Z']}"
    )
    # The single-phase state's `Z` is the one the case records, and it is followed in the
    # capture by the flashed gas and aqueous phases'.
    assert first["Z"] == "0.860124667999392", first["Z"]

    flash = (tool.CAPTURES / "saft_vr_mie_flash.tsv").read_text(encoding="utf-8")
    assert tool.SHAPES["flash"](flash)[1]["init(0) first beta"] == "0.535808206580083"


def test_a_bare_key_line_must_carry_a_number() -> None:
    """`UmrCpaProbe` writes `key   value` with no `=`, so a whitespace branch reads it.

    **That branch must not read prose.** A probe that printed a sentence would otherwise
    set a key to a word, and the case would fail to parse rather than to match - so the
    second token has to be a number before anything is taken.
    """
    tool: Any = generator()
    assert tool.pairs("alpha_T   0.839630150810119")["alpha_T"] == "0.839630150810119"
    assert tool.pairs("phase                      PhaseUMRCPA") == {}
    assert tool.pairs("some sentence with no number in it") == {}


def test_every_mapped_field_is_an_output_of_the_model() -> None:
    """A case checks fields its model declares, and nothing else.

    `test_validation_cases.py` reports an unknown field as a failure rather than skipping
    it, so a typo in the table is caught there - but only once the case runs. This catches
    it here, against the spec, which is where the name is actually declared.

    It is also the check that a case records only what its capture carries: `PcsaftProbe`
    and `SaftVrMieProbe` print no `HresTP`, so those cases check the state and not the
    departures. That is not the same as the models having none - both have them, and both
    are held to a difference of NeqSim's own energy rather than to a printed number.
    """
    from azoth import _models_gen

    tool: Any = generator()
    for case in tool.CASES:
        spec = _models_gen.model(case.calc)
        assert spec is not None, f"{case.calc} is not a model"
        declared = set(spec["outputs"])
        mapped = {mapping.under or mapping.field for mapping in case.mappings}
        assert mapped <= declared, f"{case.id}: {sorted(mapped - declared)} are not outputs"
        assert "cp_res" not in declared, (
            "no shipped model reports a residual heat capacity, and nothing here should claim one"
        )
