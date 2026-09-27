# The session

[Processes and channels](./process.md) states a unit operation as a process, and [the
rendering](./view.md) states what an editor draws from it. Between them is the layer
`docs/src/architecture/middleware.md` calls the session: the document, the graph a canvas draws
from it, and the doors a front-end reaches it through. This page fixes what that layer's claims
are, and `Azoth/Graph.lean` and `Azoth/Session.lean` prove the ones that are about functions.

## The graph, and what an edge may join

A connection is not a wire. A port declares a **channel type** — a finite record of fields, each
at a dimension — and a connection is legal when the two records agree: the same field names, each
at the same dimension and shape.

**Claim (an accepted edge supplies what the consumer declares).** If the checker accepts a
producer/consumer pair, then every field the consumer's port declares is a field the producer's
port declares, at the same dimension and shape.

*Status: **proved**, as `Azoth.Graph.check_sound`, with `Azoth.Graph.check_refl` for the case the
documents are made of — a record agrees with itself, which is what makes the judgment total on
what a palette can produce. Two witnesses keep it from being vacuous:
`Azoth.Graph.check_refuses_a_missing_field` shows the judgment can refuse, so soundness is not a
sentence about a function that always accepts; and
`Azoth.Graph.compatibility_is_not_equality` shows two records that declare the same fields in a
**different order** passing the judgment and not being equal, so the claim is about the record
and not about the way it was written down.*

**Measured: the rule cannot fire on the shipped palette.** All 77 of its port field-sets declare
the same five dimensions — `molar_flow`, `dimensionless`, `pressure`,
`thermodynamic_temperature`, `molar_energy`, five declarations each — so every pair of its ports
is compatible and no connection in `specs/` reaches the refusal. What reaches it is a palette a
caller loaded, which `Workspace::open` takes; the checker's `Mismatch` names the field, so a
front end marks the field rather than the edge when one exists to mark.

## The door is a name

`crates/azoth-cli/src/session.rs` holds one `Session`, and five doors call it: `/rpc` and `/mcp`
over HTTP, `azoth mcp` over stdio, the wasm module in a page, and `azoth.process.Session` in a
notebook. The claim is that the door changes nothing about what a call does.

**Claim (every door applies the same transition).** What differs between the doors is how a
refusal is *expressed* — an HTTP status, a JSON-RPC error code, an exception — and nothing about
the edit an accepted call makes. So an edit through one door is visible through the other.

*Status: **proved**, as `Azoth.Session.every_door_uses_the_shared_transition`, which reads a
table of the five doors and closes by `decide` — the shape `Azoth/Gate.lean` has, so a door added
with a transition of its own fails to close. `Azoth.Session.the_doors_do_not_agree_on_a_status`
is the witness beside it: the statuses differ, so the lemma is about a table whose rows are
distinguishable rather than one row written five times. What the *statuses* are is each
transport's own business and is not this layer's to fix; `crates/azoth-cli/tests/mcp_http.rs`
reads an edit made through `/mcp` back through `/rpc`, which is the claim at the boundary.*

## What `dirty` is for

A session holds the document, the values a run left **and the document those values were computed
from**, and a flag saying whether a widget should treat them as current.

**Claim (where there are values, they are the document's).** An edit keeps the values and changes
the document, so it breaks that; a run restores it; and a run that failed leaves no values at all
rather than the previous answer.

*Status: **proved**, as `Azoth.Session.a_run_restores_it` and, for the other direction,
`Azoth.Session.an_edit_leaves_the_values_behind` — which is what makes the flag a *warning* rather
than bookkeeping: the state it marks is reachable and the invariant really is broken there.
`Azoth.Session.a_failed_run_keeps_none` is the code's own sentence, "they described a document
that no longer exists, and showing them beside the current one is the mistake this flag exists to
prevent", stated as the thing the operation returns.*

## A resource is a read

**Claim (a read is not an edit).** A resource a client may read answers with the session it was
given, unchanged, for every resource — so a resource surface cannot become a second way to change
a document.

*Status: **specified**. `Azoth/Session.lean` carries the operator — `Azoth.Session.read` over a
closed set of resources, and `Azoth.Session.reading_leaves_the_state_alone` — and the *thing it
would operate on* is not built: `azoth mcp` and `azoth serve` serve tools, and no resource is
published. That is the reading [Reflection and feedback](./rho.md) takes of its own unbuilt half,
and the tranche that publishes the surface is the one that flips this line to proved.*

*Enforcement: construction — `crates/azoth-process/src/check.rs` decides a connection's legality and names the field when it refuses, `crates/azoth-cli/src/session.rs` is the one session every transport calls, and `crates/azoth-process/src/middleware/session.rs` is where an edit sets the flag the claim above is about.*
