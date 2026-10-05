"""The Rust backend's adapters are generated, and this is the ratchet that says so.

`tools/gen_python_bridge.py` emits `python/src/azoth/_rust_bridge_gen.py` for every id whose
signature is the public wrapper's and whose boundary is the spec's own inputs; `_rust_bridge.py`
re-exports those and keeps the rest. `docs-drift` regenerates and diffs it, and the `--check`
below runs the same gate in the suite.

The hand-written cap is the other half. What the generator cannot derive is the ids whose boundary
is a `Mixture` or a record the spec does not name - and the count of those is a decision someone
made rather than a line that drifted in, so it is held here.
"""

from __future__ import annotations

import ast
import importlib
import re
import subprocess
import sys
from pathlib import Path
from types import ModuleType

REPO_ROOT = Path(__file__).resolve().parents[2]
BRIDGE = REPO_ROOT / "python" / "src" / "azoth" / "_rust_bridge.py"

#: The adapters `_rust_bridge.py` may still carry by hand, and the cause each one is left for:
#: three that take a keycard the public wrapper resolves into an argument before the bridge is
#: reached (`ge_uniquac_phase`, `uniquac_activity_coefficients`, `orifice_flow`); one whose record
#: is per-component *arrays* rather than fields (`bwrs_phase`); and four single ids -
#: `capillary_dew_point` and `rate_based_packed_column`, whose optional inputs the extension takes
#: non-optionally with no default to send, `hydrogen_phase`, whose wrapper takes a mode the spec
#: does not declare, and `pure_saturation`, which takes the critical constants the spec resolves
#: from a component name. Lowering this is the point, and raising it needs a reason in the diff.
KNOWN_HAND_WRITTEN = 8


def _tools_module(name: str) -> ModuleType:
    sys.path.insert(0, str(REPO_ROOT / "tools"))
    try:
        return importlib.import_module(name)
    finally:
        sys.path.pop(0)


def _hand_written() -> list[str]:
    """Every registered id `_rust_bridge.py` still writes out itself."""
    rust_index = _tools_module("rust_index")
    registered = {function for _, _, function in rust_index.implementations()}
    text = rust_index._blank_comments(BRIDGE.read_text(encoding="utf-8"))
    found = re.findall(r"^def (\w+)\(", text, re.M)
    return [name for name in found if name in registered]


def test_the_committed_adapters_are_current() -> None:
    """The drift gate `docs-drift` runs, run here so it fails in the suite too."""
    result = subprocess.run(
        [sys.executable, str(REPO_ROOT / "tools" / "gen_python_bridge.py"), "--check"],
        capture_output=True,
        text=True,
        check=False,
    )

    assert result.returncode == 0, (
        f"the generated adapters are out of date - run `python tools/gen_python_bridge.py`\n"
        f"{result.stdout}{result.stderr}"
    )


def test_the_hand_written_adapters_are_capped() -> None:
    """Every adapter left in `_rust_bridge.py` is one the generator cannot derive."""
    hand = _hand_written()

    assert len(hand) == KNOWN_HAND_WRITTEN, (
        f"`_rust_bridge.py` carries {len(hand)} hand-written adapters against the capped "
        f"{KNOWN_HAND_WRITTEN}: {sorted(set(hand))}. If the new one is genuinely underivable, "
        f"raise the cap here and say why; if its signature is the public wrapper's and its inputs "
        f"are the spec's, it belongs in the generated file instead."
    )


def test_no_id_is_both_generated_and_hand_written() -> None:
    """The two sets are disjoint, and together they are the registry.

    A name defined twice would shadow rather than fail, which is the quiet half of this: the
    generated one would win and the hand-written one would be dead code nobody noticed.
    """
    generator = _tools_module("gen_python_bridge")
    generated = {calc_id.rpartition(".")[2] for calc_id, _, _, _ in generator.covered()[0]}
    registered = {function for _, _, function in _tools_module("rust_index").implementations()}

    both = generated & set(_hand_written())
    assert not both, f"an id is both generated and hand-written: {sorted(both)}"
    assert generated | set(_hand_written()) == registered, (
        "the generated and hand-written sets do not account for the registry: "
        f"{sorted(registered - generated - set(_hand_written()))}"
    )


def test_a_scalar_the_wrapper_calls_a_float_is_not_refused_by_the_bridge() -> None:
    """**A parameter the caller may leave out has to be a shape the bridge accepts.**

    A public wrapper's default is what reaches the extension when the argument is not stated, and
    `input_to_si` refuses a bare number where the spec declares a unit. So a dimensioned input the
    signature defaults to a bare number is a call the Rust backend cannot make:
    `process.distillation_column`'s `temperature_tolerance: float = 1.0e-6` was one, and the three
    side-draw tests that omit it were red for four commits because nothing read the pair. The
    annotation decides, exactly as it does for a vector - a parameter annotated `Q` alone is one
    the caller must state as a quantity, and `_si` is the lax half for everything else - so this
    asserts the emitted call matches the signature rather than restating what `argument_expr` does.
    """
    generator = _tools_module("gen_python_bridge")
    wrappers = generator.public_signatures()
    specs = generator.specs()
    # **The committed file, not the generator's own output.** Reading `covered()` would compare
    # `argument_expr` with itself and pass whatever it emitted, which is the one thing this has to
    # be able to fail on.
    text = (REPO_ROOT / "python" / "src" / "azoth" / "_rust_bridge_gen.py").read_text("utf-8")
    bodies = {
        node.name: ast.get_source_segment(text, node) or ""
        for node in ast.parse(text).body
        if isinstance(node, ast.FunctionDef)
    }

    checked, refused = 0, []
    for function_name, source in bodies.items():
        calc_id = f"{generator.function_ids().get(function_name, '')}"
        function = wrappers.get(function_name)
        if function is None:
            continue
        for parameter in function.args.posonlyargs + function.args.args + function.args.kwonlyargs:
            declaration = specs[calc_id].get("inputs", {}).get(parameter.arg)
            if declaration is None:
                continue
            unit = declaration.get("unit")
            if unit is None or unit == "dimensionless":
                continue
            annotation = ast.unparse(parameter.annotation) if parameter.annotation else ""
            # **`Q` and not `float`.** Measured over the six namespaces, the annotations a wrapper
            # writes for a dimensioned scalar are `Q`, `float`, `Q | None`, `float | None` and
            # `float | Q`, and only the first two say a quantity is the only shape accepted.
            if "Q" in annotation and "float" not in annotation:
                continue
            checked += 1
            if f'input_to_si(spec, "{parameter.arg}", ' in source:
                refused.append(f"{calc_id}.{parameter.arg} ({annotation or 'unannotated'})")

    assert checked, "no adapter takes a dimensioned scalar its wrapper does not call a quantity"
    assert not refused, (
        f"{refused} are dimensioned inputs the public wrapper does not require as quantities, and "
        f"the generated call refuses the bare number that is then its only shape"
    )


def test_every_molar_mass_boundary_finds_the_sentence_it_carries() -> None:
    """**A guard whose words are read out of another file has to find them.**

    `_molar_masses` raises with the sentence the id's own reference implementation already uses,
    looked for in that file and then in the public wrappers it calls. A reader that stops looking
    leaves the id hand-written, which is safe - and indistinguishable from an id that genuinely has
    no guard, which is exactly the shape that hid five ids from the wrapper generator. So the ids
    whose kernel takes a `molar_mass` behind a `mixture` are compared against the ones the reader
    found a sentence for, rather than against the cap.
    """
    generator = _tools_module("gen_python_bridge")
    reasons = generator.molar_mass_reasons()
    wrappers = generator.public_signatures()
    wanted = sorted(
        function
        for function, types in generator.wrapper_types().items()
        if "molar_mass" in types
        and function in wrappers
        and "mixture" in {p.arg for p in wrappers[function].args.args}
    )

    assert wanted, "no kernel takes a molar mass behind a mixture - nothing is being checked here"
    missing = [name for name in wanted if not any(k.endswith(f".{name}") for k in reasons)]
    assert not missing, (
        f"{missing} take a molar mass and the reference reader found no sentence to raise with, so "
        f"their guard is silently absent rather than the id being left hand-written on purpose"
    )
