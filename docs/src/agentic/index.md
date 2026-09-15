# The agentic layer

azoth is a library, and an agent that uses it has to know how — which units, which
calculation ids, what a result's warnings mean, when to reach for a keycard. That
knowledge is written down here as **skills**: self-contained specialisms an agent
loads when a task matches, each a `SKILL.md` with runnable examples and tests.

The layer is native to the repository rather than an external project. Skills live
under [`skills/`](../../../skills/README.md), their machine-readable metadata in
[`skills.toml`](../../../skills.toml), and
[`tools/validate_skills.py`](../../../tools/validate_skills.py) holds every skill to
the shape set out on [Skills](./skills.md).

What exists today is **specialisms** — one skill and the standard it follows.
Orchestration, and the agentic engineering teams built on it, come later under
[`agents/`](../../../agents/README.md), where an agent chains skills and carries no
method of its own.
