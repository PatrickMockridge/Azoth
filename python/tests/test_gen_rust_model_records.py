"""The process-model record generator, and the ratchet that keeps it the only writer.

`crates/azoth-process/src/models/*.rs` used to carry one `impl CalcResult` block each. Every line
of one is a function of two declarations that already exist - the spec's id and the struct's own
`pub` fields - so `tools/gen_rust_model_records.py` emits them into `model_records_gen.rs`.

Three things could go wrong quietly, and this file holds each:

* a spec edited without regenerating, which the drift job catches but only in CI;
* a model added without a block, which would fail the registry contract late and cryptically;
* a block hand-written back into a model file. Two impls of one trait for one type do not compile,
  so Rust catches it - but it catches it as a duplicate-impl error somewhere in the crate rather
  than as "model `mixer` grew back the code the generator owns", which is what a reader needs.
"""

from __future__ import annotations

import importlib
import re
import subprocess
import sys
from pathlib import Path
from types import ModuleType
from typing import Any

REPO_ROOT = Path(__file__).resolve().parents[2]
MODELS = REPO_ROOT / "crates" / "azoth-process" / "src" / "models"
MODULE_DECL = re.compile(r"^pub mod ([a-z_][A-Za-z0-9_]*);", re.MULTILINE)


def _tools_module(name: str) -> ModuleType:
    sys.path.insert(0, str(REPO_ROOT / "tools"))
    try:
        return importlib.import_module(name)
    finally:
        sys.path.pop(0)


def test_the_committed_records_are_what_the_generator_emits() -> None:
    """The gate the drift job runs, run here too so a stale tree fails in the suite."""
    result = subprocess.run(
        [sys.executable, str(REPO_ROOT / "tools" / "gen_rust_model_records.py"), "--check"],
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 0, (
        "model_records_gen.rs is out of date - run "
        f"`python tools/gen_rust_model_records.py`\n{result.stdout}{result.stderr}"
    )


def test_no_model_file_hand_writes_its_calc_result_block() -> None:
    """The ratchet: the generator is the only writer of these blocks."""
    offenders = sorted(
        path.name for path in MODELS.glob("*.rs") if "impl CalcResult for" in path.read_text()
    )
    assert offenders == [], (
        "these model files hand-write a block the generator emits: " + ", ".join(offenders)
    )


def test_every_model_module_has_a_generated_block() -> None:
    """A new model file with no block would be a silent hole in the registry."""
    modules = MODULE_DECL.findall((MODELS / "mod.rs").read_text())
    assert modules, "no `pub mod` lines found in models/mod.rs - the reader is wrong"
    gen = _tools_module("gen_rust_model_records")
    covered = {result.calc_id for result in gen.process_results()[0]}
    missing = sorted(name for name in modules if f"process.{name}" not in covered)
    assert missing == [], f"no generated record for: {', '.join(missing)}"


def test_the_emitter_renders_a_result_the_tree_does_not_have() -> None:
    """A generator only ever run over the tree's own shape has untested branches."""
    rust_index = _tools_module("rust_index")
    synthetic: Any = rust_index.ResultType(
        calc_id="process.synthetic",
        rust_name="SyntheticResult",
        item_path="azoth_process::models::synthetic::SyntheticResult",
        module_file=MODELS / "synthetic.rs",
        fields=(("outlet_n", "outlet_n"), ("warnings", "warnings")),
    )
    block = _tools_module("gen_rust_model_records").emit_block(synthetic)
    assert "impl CalcResult for SyntheticResult {" in block
    assert "const CALC_ID: &'static str = \"process.synthetic\";" in block
    assert '"outlet_n",' in block
    assert "&self.warnings" in block
