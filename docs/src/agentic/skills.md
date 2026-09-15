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

## Domains

Skills are grouped under `skills/<domain>/`:

- `library/` — cross-cutting skills about using and extending azoth itself.
- `eos/`, `hydraulics/`, `thermal/` — `azoth`-basis skills for the calculations the
  library implements today.
- NeqSim's ten engineering domains — `process/`, `safety/`, `flow-assurance/`,
  `pvt/`, `subsea/`, `field-development/`, `environment/`, `subsurface/`,
  `reporting/`, `engineering-data/` — where a `screening`/`advisory`/`data-retrieval`
  skill is a placeholder until the port reaches the tranche that backs it.

See [The roadmap](./roadmap.md) for which tranche backs each domain.

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
| `tranche` | optional; the `P0`–`P12` tranche that will back a `screening` skill |
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
orchestrator that chains skills and carries no method of its own. The one agent today is
the HAZOP team under [`agents/hazop/`](../../../agents/hazop/README.md); its formal
statement is [Orchestration](./orchestration.md).

## The skills today

Eleven `azoth`-basis skills cover what the library implements: four under
`library/`, four under `eos/`, two under `hydraulics/` and one under `thermal/`.
[`azoth-run-calculation`](../../../skills/library/run-calculation/SKILL.md) is the
first and the pattern the rest follow. The remainder of NeqSim's catalog is mapped,
tranche by tranche, on [The roadmap](./roadmap.md).
