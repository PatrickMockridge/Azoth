"""The simple wrappers are generated, and the hand-written ones are capped.

`tools/gen_python_wrappers.py` emits `crates/azoth-python/src/wrappers_gen.rs` - one
`#[pyfunction]` per id whose kernel parameter names are the spec's declared input names, which is
the plan's own "the signature *is* the declared input list". `docs-drift` regenerates and diffs
it; the `--check` below runs the same gate in the suite.

The cap is the other half. What this file cannot derive is the mixture-expanded models - the `eos`
ids whose kernels take a `Mixture` the spec does not declare - plus the handful whose kernel
signature carries a conversion this file has no rule for, a slice of records being the newest.
Those stay hand-written in the namespace modules, and this file holds their number, so a new one
is a decision someone made rather than a line that drifted in.
"""

from __future__ import annotations

import importlib
import re
import subprocess
import sys
from pathlib import Path
from types import ModuleType

REPO_ROOT = Path(__file__).resolve().parents[2]
BINDING = REPO_ROOT / "crates" / "azoth-python" / "src"

#: The hand-written `pub fn`s the namespace modules may carry. Measured when the generator
#: landed; lowering it is the point, and raising it needs a reason in the diff.
#:
#: 47 -> 46 when `process.packed_column` moved: slice 5c's last change made the generator able to
#: derive it, and the committed `wrappers_gen.rs` was not regenerated, so for four commits the id
#: was both covered and hand-written - which the disjointness below also catches.
#:
#: 46 -> 45 when `eos.wilson_activity_coefficients` moved, and the cause was the bridge rather than
#: the wrapper: this file decides coverage from the bridge's own call, and the *hand-written*
#: bridge passed eight arguments where the `#[pyfunction]` takes eleven. Generating the bridge made
#: the call the extension's own, the arity matched, and the emitted wrapper was byte-identical to
#: the one in `eos.rs` - bar the `build_mixture_with_mass` every other mixture id already uses.
#:
#: 45 -> 52 when `characterization` joined the scanned modules, and the change is to the guard
#: rather than to the tree: the tuple below listed six of the seven binding modules, so every
#: wrapper in `characterization.rs` was outside the cap. Seven were there, all of them registered
#: ids, and `characterize_to_reference` had just become the eighth thing a hand-written wrapper
#: would have been added to a module nothing counted - which is exactly what this number exists to
#: notice. No wrapper was written or deleted by the change.
#:
#: 52 -> 50 when the two record wrappers moved, and this time the count came *down*: `convert`
#: learned the slice-of-records rule (`&[Row]`, where the row is a spec `row_type`), and
#: `tools/gen_python_records.py` emits the `FromPyObject` structs the rule builds each row into.
#: `assay_mass_fractions` and `characterize_to_reference` are the two ids whose kernels take one.
#:
#: 50 -> 51 when `eos.water_dew_point_temperature` landed, and this one is genuinely underivable
#: rather than a guard being widened again. Its kernel takes `(components, z, p, cubic)`, and the
#: generator emits a wrapper whose parameters are the spec's declared *inputs* - `P` and `eos`
#: here - so the two names disagree and the id has no generated wrapper. `tp_multiflash_wax` and
#: `hydrate_fraction` are the same shape for the same reason: the fluid crosses as names and this
#: side resolves them, because the fixed point starts from each `water`-or-`MEG` component's own
#: `melting_point_temperature` and that is a databank column rather than a caller's argument.
#: 51 -> 52 when `eos.water_dew_point_multiphase` landed, the same shape as the sibling above
#: and for the same reason - a name-taking kernel with an `eos` enum the spec declares as `eos`
#: where the kernel calls it `cubic`, so the generated wrapper's parameter list is not the
#: kernel's and the generator refuses rather than guessing at the mapping.
#: 52 -> 53 when `eos.water_dew_point_line` landed, again the same shape: its kernel takes
#: `(components, z, minimum_pressure, maximum_pressure, eos)` while the spec declares `eos` and
#: the two pressures, so the generated wrapper's parameters are not the kernel's.
#: 53 -> 54 when `eos.saturate_with_water` landed, and its refusal is a different one: the
#: kernel's parameter names *are* the spec's inputs, and the generator stops at `associating`,
#: which the spec declares optional and the kernel takes as a `bool`. An optional flag crosses
#: the boundary as its default rather than as an absence, and no emitted signature can say
#: which; the wrapper spells the default once and the spec states the input.
#: 54 -> 55 when `reactions.reactive_co2_brine_phase_equilibrium` landed: the kernel names the
#: two state variables `temperature` and `pressure` where the spec declares `T` and `P`, the
#: same name mismatch the first two carry.
#: 55 -> 56 when `eos.lgt_surface_tension` landed: its kernel takes a `Mixture` built from the
#: `components` names and a `Vec<Vec<f64>>` interaction matrix. The first is the record/mixture
#: shape the generator refuses by design; the second is a matrix whose row length is a
#: per-component quantity the spec states in prose, not in a column count.
#: 56 -> 57 when `eos.gt_surface_tension_simple` landed, the same shape as LGT and for the
#: same two reasons.
#: 57 -> 58 when `eos.gt_surface_tension` landed, again the `components` shape.
#: 58 -> 59 when `eos.cdft_surface_tension` landed, the same `components` shape with two
#: dimensionless vector inputs the generator's row rules cannot size.
#: 59 -> 60 when `eos.interface_properties` moved: its kernel now takes a `Mixture` and
#: the interaction matrix, the same two reasons as LGT, so the generated wrapper it had
#: could no longer be derived.
#: 60 -> 59 when the generator learned an `interval = true` input:
#: `thermal.conduction_plane_wall`'s `dT` is a `TemperatureInterval` it had refused, so
#: that binding - and the hand-written `thermal.rs` module with it - could go.
#: 59 -> 60 when `process.pump`'s `curve` became an optional `record_list`: the generator
#: refuses an `Option<&[Row]>` (a record slice it can derive, an absent one it cannot), so
#: the binding moved into the process module where the row is built by hand.
#: 60 -> 61 when `process.filter` gained two optional `record_list`s for the same reason.
KNOWN_HAND_WRITTEN = 61


def _tools_module(name: str) -> ModuleType:
    sys.path.insert(0, str(REPO_ROOT / "tools"))
    try:
        return importlib.import_module(name)
    finally:
        sys.path.pop(0)


def _hand_written() -> list[str]:
    """Every registered id whose wrapper is still written out in a namespace module.

    **Every binding module is read**, which is the point: a module left out of this tuple is a
    module whose hand-written wrappers the cap cannot see, and the cap is the only thing standing
    between "the generator cannot derive this" and "somebody wrote one and moved on".
    """
    rust_index = _tools_module("rust_index")
    registered = {function for _, _, function in rust_index.implementations()}
    out = []
    for source in sorted(BINDING.glob("*.rs")):
        # The generated modules emit the same shape, so they are the one thing that has to be
        # named: a `wrappers_gen.rs` counted here would put every derived wrapper inside a cap on
        # the underived ones.
        if source.name.endswith("_gen.rs"):
            continue
        text = rust_index._blank_comments(source.read_text(encoding="utf-8"))
        out.extend(re.findall(r"^pub fn (\w+)\(\s*py: Python<'_>", text, re.M))
    return [name for name in out if name in registered]


def test_the_committed_wrappers_are_current() -> None:
    """The drift gate `docs-drift` runs, run here so it fails in the suite too."""
    result = subprocess.run(
        [sys.executable, str(REPO_ROOT / "tools" / "gen_python_wrappers.py"), "--check"],
        capture_output=True,
        text=True,
        check=False,
    )

    assert result.returncode == 0, (
        f"the wrappers are out of date - run `python tools/gen_python_wrappers.py`\n"
        f"{result.stdout}{result.stderr}"
    )


def test_the_hand_written_wrappers_are_capped() -> None:
    """Every wrapper left in a namespace module is one the generator cannot derive."""
    hand = _hand_written()

    assert len(hand) == KNOWN_HAND_WRITTEN, (
        f"the namespace modules carry {len(hand)} hand-written wrappers against the capped "
        f"{KNOWN_HAND_WRITTEN}: {sorted(set(hand))}. If the new one is genuinely underivable, "
        f"raise the cap here and say why; if its kernel's parameter names are the spec's declared "
        f"inputs, it belongs in the generated file instead."
    )


def test_no_id_is_both_generated_and_hand_written() -> None:
    """The two sets are disjoint. A wrapper emitted twice would not link, and one removed from
    the modules while the generator also refused it would be an id nothing registers."""
    wrappers = _tools_module("gen_python_wrappers")
    generated = {function for _, function in wrappers.covered_ids()}

    assert generated & set(_hand_written()) == set(), "an id is both generated and hand-written"


def test_the_boundary_reader_sees_both_halves_of_the_bridge() -> None:
    """**A file the reader cannot see is a check that has been turned off.**

    `bridge_call` and `bridge_arguments` decide the emitted wrapper's arity and parameter order,
    and they read the bridge for the call. When the bridge's adapters were generated, most of the
    calls moved into `_rust_bridge_gen.py`; a reader that read only `_rust_bridge.py` then found no
    call for those ids, which *passed* both checks rather than failing them - and admitted five ids
    the generator had refused on purpose (`phase_transport`, `absorption_column`, `ejector`,
    `heat_exchanger`, `stripping_column`) while the committed `wrappers_gen.rs` still lacked them,
    so `--check` was red for the whole of the next slice. Every generated adapter makes exactly one
    `_core.<name>(` call, so the reader has to find one for each.
    """
    wrappers = _tools_module("gen_python_wrappers")
    bridge = (REPO_ROOT / "python" / "src" / "azoth" / "_rust_bridge_gen.py").read_text("utf-8")
    block = bridge.partition("__all__ = [")[2].partition("]")[0]
    names = re.findall(r'^    "(\w+)",$', block, re.M)

    assert names, "the generated bridge declares no adapters"
    unseen = sorted(name for name in names if wrappers.bridge_call(name) is None)
    assert not unseen, (
        f"{len(unseen)} generated adapter(s) whose `_core` call the boundary reader cannot see, so "
        f"its arity and order checks are silently off for them: {unseen}"
    )
