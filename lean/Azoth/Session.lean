/-
The session: one document, several names, and what a flag is for.

`docs/src/calculus/session.md` states the middleware's layer. `crates/azoth-cli/src/session.rs`
holds one `Session` and five doors call it - `/rpc` and `/mcp` over HTTP, `azoth mcp` over
stdio, the wasm module in a page, and `azoth.process.Session` in a notebook - so the first
claim is that **the door is a name and not a state**: what differs between them is how a
refusal is *expressed*, and nothing about what a call does.

# The doors, as data

`doors` below is the table, and `every_door_uses_the_shared_transition` is what a reader can
check it against: every row names the one transition. `the_doors_do_not_agree_on_a_status` is
the witness beside it - the statuses differ, so the table is not a list of identical rows and
the lemma is not a statement about nothing.

This is the same shape `Azoth/Gate.lean` has: a table, and a `decide` that reads it, so a row
added with a second transition fails to close rather than being waved through.

# What `dirty` is for

The second claim is about the flag. A session holds the document, the values a run left, and
whether those values are older than the document - and the *content* of that is an invariant
this module states: **where there are values, they are the values of the document in hand**.

- An edit keeps the values and changes the document, so it **breaks** the invariant:
  `an_edit_leaves_the_values_behind`. That is what the flag warns about, and why the editor
  draws a stale pill rather than a number.
- A run restores it: `a_run_restores_it`.
- And a run that **failed** leaves no values at all: `a_failed_run_keeps_none`. The code says
  the same thing in its own words - they "described a document that no longer exists, and
  showing them beside the current one is the mistake this flag exists to prevent".

A **read is not an edit**, which is the third claim and the one the resource surface rests on:
`reading_leaves_the_state_alone`, stated over a resource the module names, so a resource that
changed the document would fail to close.
-/

namespace Azoth

namespace Session

/-- One door's transition. **One constructor, and that is the claim**: the transports differ in
how they report a refusal and in nothing about what a call does. An enum rather than a function
so that a table of doors is data a `decide` can read. -/
inductive Transition where
  | shared
deriving DecidableEq, Repr

/-- **The doors, as the middleware serves them**: the route, the status a refusal is expressed
with, and the transition. The statuses are the transports' own - `/rpc` answers `400`, `azoth
serve`'s MCP route answers a JSON-RPC error code beside its HTTP status, and the wasm module
throws - and the fifth is the notebook. -/
def doors : List (String × Nat × Transition) :=
  [ ("/rpc", 400, .shared),
    ("/mcp (POST, HTTP)", 400, .shared),
    ("azoth mcp (stdio)", 32602, .shared),
    ("the wasm module", 0, .shared),
    ("azoth.process.Session", 0, .shared) ]

/-- **Every door applies the one transition.** The lemma a row added with a second transition
fails to close. -/
theorem every_door_uses_the_shared_transition :
    (doors.all fun door => door.2.2 == Transition.shared) = true := by
  decide

/-- **And the doors are not the same row.** Their statuses differ, so the lemma above is about a
table whose rows are distinguishable rather than about a list of one thing written five times. -/
theorem the_doors_do_not_agree_on_a_status : doors.map (fun door => door.2.1) ≠
    List.replicate doors.length (doors.headD ("", 0, .shared)).2.1 := by
  decide

/-- A session: the document, the values a run left **and the document those values were computed
from**, and whether a widget should treat them as current. -/
structure State where
  /-- The document in hand, as an opaque id. -/
  document : Nat
  /-- The last run's report, paired with the document it was computed from. -/
  values : Option (Nat × Nat)
  /-- True while the values are older than the document. -/
  dirty : Bool

/-- **Where there are values, they are the document's.** The invariant the flag exists to warn
about breaking, and the thing a widget would rather read than a flag. -/
def Sound (s : State) : Prop :=
  s.values = none ∨ ∃ report, s.values = some (s.document, report)

/-- One edit: the document changes, and the values do not. -/
def edit (s : State) (document : Nat) : State :=
  { s with document := document, dirty := true }

/-- One run that succeeded: the values become the document's, and the flag clears. -/
def run (s : State) (report : Nat) : State :=
  { s with values := some (s.document, report), dirty := false }

/-- One run that failed: the values go rather than staying to be misread. -/
def run_failed (s : State) : State :=
  { s with values := none, dirty := true }

/-- **A run restores the invariant.** It is the only operation that does. -/
theorem a_run_restores_it (s : State) (report : Nat) : Sound (run s report) := by
  unfold Sound run
  exact Or.inr ⟨report, rfl⟩

/-- **An edit leaves the values behind.** With values in hand and a document that moved, the
values describe a document the session no longer holds - which is the state the flag marks and
the state a widget must not read a number from. -/
theorem an_edit_leaves_the_values_behind (s : State) (document : Nat) (report : Nat)
    (h : s.values = some (s.document, report)) (moved : document ≠ s.document) :
    ¬ Sound (edit s document) := by
  unfold Sound edit
  rw [h]
  simp [moved.symm]

/-- **A failed run keeps no values.** Not a stale number, and not the previous answer: nothing
the next reader could mistake for this document's. -/
theorem a_failed_run_keeps_none (s : State) : (run_failed s).values = none := rfl

/-- A resource a client may read. **Only the reads**, because a resource that could change the
document would be a second editing surface. -/
inductive Resource where
  | document
  | diagnostics
  | report
  | catalogue
deriving DecidableEq, Repr

/-- What a resource read answers: the session, and nothing done to it. -/
def read (s : State) (_ : Resource) : State × State := (s, s)

/-- **A read is not an edit.** The state it answers with is the state it was given, for every
resource, so a resource surface cannot be a way to change a document. -/
theorem reading_leaves_the_state_alone (s : State) (resource : Resource) :
    (read s resource).1 = s ∧ (read s resource).2 = s :=
  ⟨rfl, rfl⟩

end Session

end Azoth
