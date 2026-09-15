# Skills

A skill is a self-contained specialism an agent loads when a task matches. It is
guidance plus runnable examples and tests, never the engineering method itself: the
method is azoth's, and a skill teaches how to call it.

## The shape

Each skill is one folder under `skills/<domain>/<skill>/`:

```text
run-calculation/
├── SKILL.md
├── examples/
└── tests/
```

`SKILL.md` is the guidance an agent reads. `examples/` are copy-pasteable and
runnable. `tests/` hold the API surface the skill documents, so a drift between the
skill and the library fails rather than rots. `src/` is optional and omitted when the
skill has no logic of its own — the normal case, because the method is azoth's.

## Naming

- Folder and id use kebab-case: `run-calculation`.
- The catalog name is `azoth-` plus the folder: `azoth-run-calculation`.
- A domain folder is added only when an existing one is unsuitable; the first is
  `library`, for skills about using azoth itself.

## The manifest

Machine-readable metadata lives in [`skills.toml`](../../../skills.toml) — the single
source of truth, and the only place a skill's name, version, description and basis
are written. A `SKILL.md` carries no frontmatter and repeats nothing
machine-readable, so there is no second copy to drift.

Each entry declares, and [`tools/validate_skills.py`](../../../tools/validate_skills.py)
enforces:

| Field | Rule |
|---|---|
| `name` | `azoth-[a-z0-9]+(-[a-z0-9]+)*` |
| `version` | semantic version `x.y.z` |
| `description` | must contain `USE WHEN:`, so an agent can decide when to load it |
| `calculation_basis` | `azoth`, `screening`, `advisory`, `data-retrieval` or `hybrid` |
| `path` | a real `SKILL.md` under `skills/` |
| `tags` | a non-empty list |

`calculation_basis` says how a skill's numbers are produced, so an agent can weigh
rigour programmatically: `azoth` drives the validated library (the analogue of
NeqSim's `neqsim-java`); `screening` is an educational placeholder that references
azoth; `advisory` is methodology with no numeric output; `data-retrieval` returns
data and performs no calculation; `hybrid` mixes screening with at least one azoth
step.

## Required sections

Every `SKILL.md` carries these headings, in this order:

`When to Use`, `Inputs`, `Outputs`, `How a calculation runs`, `Python usage pattern`,
`The keycard`, `Validation checklist`, `Common mistakes`, `Limitations`,
`Related Azoth functionality`, `References`.

## Skill vs agent

A **skill** is a reusable method or body of knowledge. An **agent** is an
orchestrator that chains skills and carries no method of its own. Agents do not exist
yet; when they do they live under [`agents/`](../../../agents/README.md).

## The first skill

[`azoth-run-calculation`](../../../skills/library/run-calculation/SKILL.md) is the
one skill today, and the pattern the rest will follow: teach an agent to call azoth,
read a result's warnings, and supply a keycard.
