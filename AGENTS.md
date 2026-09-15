# Working with azoth

azoth is a chemical-engineering calculation library with an agentic layer on top.

- **Skills** (90) live in [`skills/`](skills/README.md): self-contained specialisms,
  each a `SKILL.md`. To load them into Claude Code or OpenCode, run
  `python tools/export_skills.py` (writes `.claude/skills/` and `.opencode/skills/`).
- **Orchestration** is in [`agents/README.md`](agents/README.md). The one agent is a
  HAZOP team of four roles, defined in `agents/hazop/README.md` and exported to
  `.claude/agents/` and `.opencode/agents/`; run `/hazop` to start a study.

Setup:

```bash
pip install 'azoth[agent]'          # DeepSeek Harness SDK + runtime binary (team runtime)
python tools/export_skills.py       # write the skills into each tool's directory
```

See [`agents/README.md`](agents/README.md) for the full setup and caveats.
