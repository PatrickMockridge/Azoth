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

- **Test what the change touches.** A change to one model cannot break another.
  The full gate set runs **once per tranche, at close-out, before the push** — never
  while iterating, never per edit.
- **Wrap every heavy command** in `tools/gated.sh`: it runs the command under
  `systemd-run --user --scope -p MemoryMax=12G -p MemorySwapMax=0` and caps the JVM
  heap at 2 GiB. Heavy = `cargo build|test|clippy|run|bench`, `pytest`, `maturin`,
  `lake`, `wasm-pack`, `mdbook`, any `java`/`javac`. `tools/oracle_sweep.py` spawns one
  `javac` plus ~124 JVMs, each defaulting to an 11.8 GiB heap.
- `-j 10` on cargo; `--test-threads 10` after `--` on tests.
- Never run two heavy commands at once.
- `.claude/hooks/guard_bash.py` (a PreToolUse hook) enforces the wrapper and the
  close-out gate. `AZOTH_FULL_GATE=1`, visible in the transcript, is the only way past
  the gate.

