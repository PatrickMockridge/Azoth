# Orchestration

The calculus states a calculation — a unit operation, a keycard, a process. This
page states an engineering team in the same language: a composition of
role-processes communicating over message channels, and a HAZOP as one such
process. It is the same machinery as the flowsheet, applied to agents rather than
to unit operations.

## A team is a process

A **role** is a named process: a procedure quoted as `@role` and re-entered as
`*role` ([Reflection](../calculus/rho.md)). A **team** is the parallel composition of its
roles:

```
team  ::=  chair | scribe | process | safety | ...
```

Each role reads the channels it owns and writes the channels downstream roles
read, as a unit operation reads its inlets and writes its outlets
([Processes](../calculus/process.md)). The channels carry **findings**, not streams: a
deviation, a cause, a safeguard.

**Claim (a team is a composition).** A team is barbed-bisimilar to the parallel
composition of its roles, each role a process whose only free names are its
declared channels.

*Status: **specified**. There is no orchestration runtime, so
`Azoth.Orchestration` does not exist; the claim is what the first agent is checked
against.*

## A HAZOP is a process

A HAZOP passes each node and guideword through a fixed pipeline:

```
guideword → deviation → cause → consequence → safeguard
```

as a process `hazop`, each stage a role reading the previous stage's channel and
writing the next. The worksheet is the serialisation of that process: `@hazop`
turns the procedure into data, and `*` re-enters it.

**Claim (a HAZOP worksheet round-trips).** Quoting a HAZOP and dropping it back is
the same HAZOP.

*Status: **specified**, and it is not independent of [rho.md](../calculus/rho.md)'s round
trip: it is that theorem at the process `hazop`. It is named here because a reader
checking the orchestration layer wants the instance, not the general statement.*

## The first agent

[`agents/hazop/`](../../../agents/hazop/README.md) is the first agent: the roles
above, the skills each chains, and the pipeline, in the terms this page fixes. It
is runnable through `python/azoth/agents`, the DeepSeek Harness boundary. The tool
schema a model is handed — `@tool` and `*tool`, the same reflection surface a GUI
side panel or a notebook offers — is [the middleware](../architecture/middleware.md),
deferred until P12.

## What is not here

There is deliberately no orchestration runtime. A team is a process the calculus
states; what runs it is the boundary in `python/azoth/agents`, and the claims above
are proved when that runtime exists — the same "specified until the layer exists"
discipline as the process and reflection layers.
