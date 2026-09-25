# azoth

azoth is a chemical-engineering calculation library (thermodynamics, hydraulics,
heat transfer) with an agentic layer on top.

- **Skills** (90) live in `skills/`, each a `SKILL.md`. To load them into Claude Code,
  run `python tools/export_skills.py --target claude` (writes `.claude/skills/`).
- **Orchestration** is in `agents/README.md`. The one agent is a HAZOP team of four
  roles, defined in `agents/hazop/README.md` and exported to `.claude/agents/` and
  `.opencode/agents/`; run `/hazop` to start a study.
- **Runtime**: `pip install 'azoth-engine[agent]'` installs the DeepSeek Harness SDK for the
  team orchestration; the skills themselves need no runtime.
