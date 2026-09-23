# Agents

The orchestration layer: skills combined into agents, and agents composed into teams.

## What this is

- A **skill** is a self-contained specialism an agent loads when a task matches. Each
  is one `SKILL.md` — the same convention Claude Code and OpenCode use — plus runnable
  examples and tests. There are 90 under [`skills/`](../skills/README.md); machine
  metadata lives in `skills.toml`, validated by `tools/validate_skills.py`.
- An **agent** is an orchestrator that chains skills and carries no engineering method
  of its own; the method lives in a skill.
- A **team** is a composition of role-processes over message channels. The formal
  statement is [`docs/src/agentic/orchestration.md`](../docs/src/agentic/orchestration.md):
  a role is a named process, a team its parallel composition.

The one agent today is **HAZOP**: four roles — chair, process, safety, scribe — over
the pipeline `guideword → deviation → cause → consequence → safeguard`. Its definition
is [`agents/hazop/README.md`](hazop/README.md).

## Use the skills

The skills are plain `SKILL.md` files, so Claude Code or OpenCode can load them
directly — point the tool at `skills/`, or copy the ones you want. No export is needed
for that path.

## Set up the team runtime (DeepSeek Harness)

1. **Prerequisites** — Python 3.12, and azoth (`pip install azoth`, or clone this
   repository).
2. **Install the runtime** — `pip install 'azoth[agent]'` pulls
   `deepseek-harness-sdk` and its bundled `dsh` runtime binary.
3. **Export the skills** — `python tools/export_skills.py --target dsh` writes the
   catalog into `.dsh/skills/`, which DeepSeek Harness discovers at project rank.

## Run the HAZOP team

`python/azoth/agents` is the boundary to that runtime — the one place that imports
`deepseek_harness_sdk`, so the dev-preview churn stays contained. It is a shim today:
the orchestration API lands with the runtime, and the runtime is not built. The agent
definition is [`agents/hazop/README.md`](hazop/README.md); the formal spec is
[`docs/src/agentic/orchestration.md`](../docs/src/agentic/orchestration.md).

## Caveats

- DeepSeek Harness is a v0.1 developer preview: compatibility-breaking changes are
  promised, and it accepts no pull requests yet. The runtime boundary is thin so the
  churn stays in one place.
- The team's formal logic is the rho-calculus
  ([`docs/src/calculus/rho.md`](../docs/src/calculus/rho.md)). Barbed *congruence* — the
  context closure of barbed bisimulation — is stated in
  [`docs/src/calculus/barbs.md`](../docs/src/calculus/barbs.md) and not yet proved.
