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
  reaches: the crates whose sources changed *and every crate that depends on them*, and the
  Python mirror of each of those crates — a kernel change runs `test_cross_impl.py`, because
  that agreement test is the whole reason every calculation is written twice. A path it
  cannot place is an error rather than a shrug, and so is a path a rule *claims and then
  places nowhere*. `--explain` says why each command is there; `--lane iterate` deselects the
  `slow_sweep`-marked global sweeps, for the inner loop, and both lanes name the ones they
  reach. **It also runs what CI runs whatever the diff says** — the five unconditional jobs, from
  `CLOSE_OUT`, so the local run is the push rather than the push's diff-shaped half — and the
  generator checks run first, because a stale generated file is minutes to find and an hour to be
  told about by a test.
- **Three checks say whether the scoped set is worth trusting.** `--audit` walks every
  tracked path through the table, so the day a directory is added and no rule knows it the
  gate says so; `--verify-selection` collects every emitted target and refuses one that
  selects nothing, because an empty collection is a green run over nothing; and the scoped
  Python run is instrumented (`AZOTH_COVERAGE=1`), with `tools/check_coverage.py` refusing a
  changed file no line of which ran — a target can collect forty tests and touch none of the
  code the diff changed. It asserts one fact per file, never a percentage: a threshold is a
  number to game, and the file nobody touched keeps its zero.
- **Mutation testing is not run, deliberately.** The instrument this repository already has is
  an independent second implementation, which catches the class a mutation score proxies for —
  a wrong operator in one language — at the cost of one comparison rather than N mutant
  builds. That is a stated limit, not an oversight.
- **On a push, four jobs run whatever the diff says** — `spec-validate`, `rust-lint`
  (`cargo fmt` only; the workspace clippy is in the nightly), `python-lint` and
  `provenance`, all flat-cost and build-free — and **the other eight run only if the diff
  reaches them**. `tools/test_scope.py --plan` makes that decision once, in the `scope` job,
  and every job reads its own line, so a job that does not run is a job the summary said
  `false` about rather than an absence nobody explained.
- **The whole suite runs nightly**, 04:17 UTC, in `.github/workflows/nightly.yml` — which
  *calls* `ci.yml` with `mode: full` rather than repeating it, so "the nightly runs
  everything" is a fact about one list. It is the only thing that runs the tests no diff
  claimed, so it is not a thing to type, and a red one is a finding to triage rather than
  scenery.
- **The trade this buys, stated plainly:** a push is verified against what it reaches, not
  against the tree, so something the nightly would have caught can be red for up to a day.
  That is what it costs to stop paying an hour per change. A release tag requires the same
  full surface, and `release.yml` runs it before building a wheel.
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

