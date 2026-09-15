# The agentic layer

azoth is a library, and an agent that uses it has to know how — which units, which
calculation ids, what a result's warnings mean, when to reach for a keycard. That
knowledge is written down as **skills**: self-contained specialisms an agent loads when a
task matches, each a `SKILL.md` with runnable examples and tests. [Skills](./skills.md)
is their shape; [the roadmap](./roadmap.md) maps them to the port tranches.

The layer is native to the repository rather than an external project. Skills live under
[`skills/`](../../../skills/README.md), their machine-readable metadata in
[`skills.toml`](../../../skills.toml), and
[`tools/validate_skills.py`](../../../tools/validate_skills.py) holds every skill to the
shape.

**The one agent is HAZOP** — a team of four roles (chair, process, safety, scribe) chained
through `guideword → deviation → cause → consequence → safeguard`. Its definition is
[`agents/hazop/README.md`](../../../agents/hazop/README.md), exported to `.claude/agents/`
and `.opencode/agents/`; the formal statement of a team as a composition of role-processes
is [Orchestration](./orchestration.md).

## Run it

```bash
pip install 'azoth[agent]'          # DeepSeek Harness SDK + runtime
python tools/export_skills.py       # write the skills into each tool's directory
```

`python/azoth/agents` is the one place that imports the DeepSeek Harness SDK; the skills
themselves need no runtime. [`agents/README.md`](../../../agents/README.md) has the full
setup and the caveats.
