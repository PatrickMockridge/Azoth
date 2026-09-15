# Working with azoth

azoth is a chemical-engineering calculation library with an agentic layer on top.

- **Skills** (90) live in [`skills/`](skills/README.md): self-contained specialisms,
  each a `SKILL.md`, loadable by Claude Code or OpenCode. Catalog in `skills.toml`.
- **Orchestration** (agents and teams) is documented in
  [`agents/README.md`](agents/README.md). The one agent is a HAZOP team of four roles.

To run the team orchestration:

```bash
pip install 'azoth[agent]'          # DeepSeek Harness SDK + runtime binary
python tools/export_dsh_skills.py   # write the skills into .dsh/skills/
```

See [`agents/README.md`](agents/README.md) for the full setup and caveats.
