# azoth

azoth is a chemical-engineering calculation library (thermodynamics, hydraulics,
heat transfer) with an agentic layer on top.

- **Skills** (90) live in `skills/`, each a `SKILL.md`. To load them into Claude Code,
  run `python tools/export_skills.py --target claude` (writes `.claude/skills/`).
- **Orchestration** is in `agents/README.md`. The one agent is a HAZOP team of four
  roles, defined in `agents/hazop/README.md` and exported to `.claude/agents/` and
  `.opencode/agents/`; run `/hazop` to start a study.
- **Runtime**: `pip install '.[agent]'`, from a checkout - `azoth-engine` is not on PyPI
  yet - installs the DeepSeek Harness SDK for the team orchestration; the skills
  themselves need no runtime.

## Running tests and builds

This box has 47 GiB RAM and 2 GiB swap, with no OOM guard (`systemd-oomd` inactive,
`earlyoom` absent), so a runaway thrashes the desktop instead of failing one process.

- **`python tools/test_scope.py --run` is the gate.** It reads the diff and runs what it
  reaches: the crates whose sources changed *and every crate that depends on them*, the
  namespaces and model cases on both backends, the generator checks whose output could have
  moved. A path it cannot place is an error rather than a shrug, so nothing escapes the map
  in silence. `--explain` says why each command is there; `--lane iterate` is the same set
  with the slow global sweeps deselected, for the inner loop.
- **The whole suite runs nightly**, on a 24-hour cycle, in
  `.github/workflows/nightly.yml` — on a runner with the memory this box lacks. It is the
  only thing that runs the tests no diff claimed, so it is not a thing to type, and a red
  one is a finding to triage rather than scenery. `python tools/test_scope.py --lane
  nightly` prints exactly what that job runs.
- **The trade this buys, stated plainly:** a push is verified against what it reaches, not
  against the tree, so something the nightly would have caught can be red for up to a day.
  That is what it costs to stop paying an hour per change. A release tag should require a
  green nightly.
- **Wrap every heavy command** in `tools/gated.sh`: it runs the command under
  `systemd-run --user --scope -p MemoryMax=12G -p MemorySwapMax=0` and caps the JVM
  heap at 2 GiB. Heavy = `cargo build|test|clippy|run|bench`, `pytest`, `maturin`,
  `lake`, `wasm-pack`, `mdbook`, any `java`/`javac`. `tools/oracle_sweep.py` spawns one
  `javac` plus ~124 JVMs, each defaulting to an 11.8 GiB heap.
- `-j 10` on cargo; `--test-threads 10` after `--` on tests.
- Never run two heavy commands at once.
- `.claude/hooks/guard_bash.py` (a PreToolUse hook) enforces the wrapper and refuses a
  repo-wide test run. `AZOTH_FULL_GATE=1`, visible in the transcript, is the only way past
  that refusal — for the rare local full run, not for a push.

