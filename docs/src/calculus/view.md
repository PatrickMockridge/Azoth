# The rendering

React and xyflow are functional rendering: a component is a function of the state it is handed,
and a canvas is a graph drawn from the document. So the editor is a *lambda* layer over the
process calculus — and the claims it makes are claims about functions, which is what this page
fixes and what `Azoth/View.lean` proves.

The layer is stated in prose at [the middleware](../architecture/middleware.md): the editor "holds
no copy of the flowsheet", the drag is "the one piece of local state", and a unit set "is a reading
of a run, not a second run". Each of those is a theorem below, and the reason to prove them is that
each is the kind of statement a refactor breaks silently — a second `useState` changes nothing a
test in `ui/` would notice unless the test is about that.

## The state is the document

`ui/src/App.tsx` holds one `Envelope` and no reducer over the graph, so a panel is a function of
the envelope and two panels given the same envelope draw the same thing. That is a fact about the
shape of the code rather than a theorem, and the interesting consequence is the *asynchrony*: the
doors answer with promises, so two edits in flight can come back the other way round, and the
editor keeps the answer whose ticket is the newest.

**Claim (the answer displayed is the last call's).** Let each call take a ticket and answer with
one; let the editor absorb answers in the order they arrive, keeping the newer ticket. Then the
answer displayed is the answer of the call with the greatest ticket issued, whatever order the
answers arrive in, and it is one of the answers rather than a combination of several.

*Status: **proved**, as `Azoth.View.foldr_keep_ticket` and `Azoth.View.foldr_keep_mem`, over the
`keep` of `Azoth.View.keep_ticket` and `Azoth.View.keep_mem`. The two are what make the claim
whole: the first says the ticket of the answer on screen is the maximum of the tickets that
arrived, and the second says that answer is one of them — without the second, a rule that
*computed* a document from the tickets would satisfy the first. `Azoth.View.keep_keeps_the_newer`
is the witness that `keep` is not a projection: an editor that displayed its first argument would
agree with `keep_mem` and fail this.*

## The drag, and the one piece of local state

xyflow applies a position change on every frame of a drag and the library is told once, on
release. So there is a display the canvas shows while the gesture is in flight and a command the
document receives when it settles, and the difference between them is the whole of the editor's
local state.

**Claim (a release settles on the document).** After a release the canvas draws the position the
document carries for that node; every other node draws what it drew before; and while the gesture
is in flight the canvas draws a position the document does not carry, which is why the release has
to send one.

*Status: **proved**, as `Azoth.View.drawn_after_release`, `Azoth.View.a_drag_moves_one_node` and
`Azoth.View.a_drag_in_flight_is_not_the_document`. The second is the one that carries "one piece":
`Azoth.View.setPosition` writes one node and leaves the function elsewhere unchanged, so the
gesture is a reading with an exception rather than a second copy of the canvas.*

## The display unit

A reader picks a set of units; a document keeps the unit its spec declares; and a field shows and
accepts the one the set names. So a conversion happens in both directions, and for a degree
Celsius it is not a scale — `si = (v + offset) * factor` — which is why `ui/src/state/units.ts`
holds a formula over the catalogue's own `factor` and `offset` rather than a table.

**Claim (a display conversion round-trips).** For any unit whose factor is not zero, converting a
value to SI and back returns the value, and converting an SI value into the display unit and back
returns that.

*Status: **proved**, as `Azoth.View.Shift.fromSi_toSi` and `Azoth.View.Shift.toSi_fromSi`. The
factor is a hypothesis and not a decoration: a unit worth nothing has no inverse, and the theorem
says which units those are rather than excluding them by a type nobody checks.*

**Claim (two points determine a unit).** Two units that agree at one value need not be the same
unit, and two that agree at two different values are.

*Status: **proved**, as `Azoth.View.one_point_does_not_determine_a_shift` and
`Azoth.View.two_points_determine_a_shift`. This is the reason
`python/tests/test_units_cross_library.py` asks `pint` for a unit's own answer at **two** values
rather than one: at a single value a scale and a shifted scale are indistinguishable, so a table
that named the wrong one would pass a one-point check. The witness pair is `⟨2, 0⟩` and `⟨1, 1⟩`,
which both answer `2` at `1`.*

## What this layer does not claim

**The projection never refuses a document**, and that is deliberate rather than proved here: a
canvas has to draw a broken flowsheet, because drawing it is how a reader fixes it. An instance
naming a unit operation the palette does not carry gets a node with no ports, and an edge naming an
instance nobody declared gets the node id it implies and no node to attach to — which
`crates/azoth-process/src/middleware/graph.rs` states as a rule and its tests hold it to. It is a
**characterisation**: nothing in this repository could violate it, because it is how the projection
is written.

**A component's purity is not a theorem either.** The layer's own state is one envelope, so a
component that read anything else would be reading something the editor does not hold; what keeps
that true is `tools/check_ui.py`'s third rule, which refuses a part of the document held in local
state, and `ui/test/currency.test.ts`, which holds the wire's mirrors to the library that emits
them.

*Enforcement: construction — `ui/src/App.tsx` holds one envelope and one ticket counter, so a panel has nothing else to read, and `ui/src/wire/decode.ts` is the parse that turns a field which moved into a refusal at the boundary rather than an `undefined` three files from it.*
