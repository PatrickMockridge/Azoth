# HAZOP orchestrator

The first agent: four roles chaining the `safety` and `process` skills through the
HAZOP pipeline `guideword → deviation → cause → consequence → safeguard`, specified
in [`docs/src/calculus/orchestration.md`](../../docs/src/calculus/orchestration.md).

Roles, and the skills each chains:

- **chair** — runs the node-and-guideword pass; sequences the other roles and
  carries no method of its own.
- **process** — `azoth-relief-load-screening`, `azoth-depressurization-screening`,
  `azoth-psv-orifice-screening`, and the rest of the `process/` catalog.
- **safety** — `azoth-safety-function-coverage-screening`,
  `azoth-vacuum-collapse-screening`, `azoth-flare-radiation-screening`, and the
  rest of the `safety/` catalog.
- **scribe** — `azoth-final-report-writing-style` and `azoth-prose-english`, to
  render the worksheet.

The worksheet is the serialisation of the run — `@hazop` — so a run re-entered
with `*` is the same run.

Runnable via `python/azoth/agents` (the DeepSeek Harness boundary). The runtime is
a dev preview; this file fixes the agent, and the runtime executes it.

Setup: [`agents/README.md`](../README.md).
