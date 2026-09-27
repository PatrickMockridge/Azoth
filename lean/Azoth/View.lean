/-
The rendering is a function of the document.

`docs/src/calculus/view.md` states the editor's layer, and this module is its proof. Three
things `docs/src/architecture/middleware.md` claims in prose are claims about functions, and
each is a theorem below.

**The editor holds no copy of the flowsheet.** The state is one envelope and every panel is a
reading of it. What that buys is the *ticket discipline* in `ui/src/App.tsx`: a call takes a
number, the door answers with a promise, and an answer is displayed only if its ticket is the
newest issued. Two edits in flight can come back the other way round, and the property that
makes that harmless is not obvious from the code - so it is the first thing proved here. The
tickets of answers in flight are distinct (one is issued per call), and under that hypothesis
the answer displayed is the one the last call left, **whatever order the answers arrive in**.

**The drag is the one piece of local state, and it does not survive its release.** xyflow
applies a position change every frame and the library is told once, on release. So there is a
function for the gesture in flight and another for the document, they agree on the node the
gesture is about, and neither touches any other node - which is the whole of "one piece of
local state".

**A display unit's conversion is a bijection, and one point does not determine it.** A degree
Celsius is affine - `si = (v + offset) * factor` - so `factor` alone cannot express it, which
is why `degC` and `degF` are a table of their own in the vocabulary. Two theorems here: the
round trip in both directions for any `factor ≠ 0`, and the reason `python/tests/
test_units_cross_library.py` asks `pint` at **two** values rather than one - two units that
agree at a single value need not be the same unit, and two values do determine one.

# What is not here

**No per-unit rows.** The display table's two units would each be an instance of
`fromSi_toSi`, which is quantified over every factor and offset and so already covers them; a
row per unit would be the same lemma written twice with the numbers substituted, and the
numbers are the generator's to state. The generic theorem is the stronger claim, and the
gate's value comes from every entry in it being one a change could break.
-/

import Mathlib.Data.Rat.Defs
import Mathlib.Tactic

namespace Azoth

namespace View

/-! ## The answer displayed is the last call's

One call, one ticket, one answer. The editor keeps the newest and drops the rest. -/

/-- One call's answer: the ticket it took, and the document it left.

`Doc` is the document's type rather than a shape here - what these theorems are about is the
*ticket*, and a model that fixed a document's representation would be about the wrong thing. -/
structure Answer (Doc : Type) where
  /-- The number `App.tsx`'s `issued` counter gave the call. Distinct per call. -/
  ticket : Nat
  /-- What that call left, which is the whole envelope. -/
  document : Doc

/-- **Keep the newer answer.** The ticket decides, and the arrival order does not. -/
def keep {Doc : Type} (a b : Answer Doc) : Answer Doc :=
  if a.ticket ≤ b.ticket then b else a

/-- The ticket of what is kept is the larger of the two. -/
theorem keep_ticket {Doc : Type} (a b : Answer Doc) :
    (keep a b).ticket = max a.ticket b.ticket := by
  unfold keep
  split_ifs with h
  · exact (max_eq_right h).symm
  · exact (max_eq_left (Nat.le_of_not_le h)).symm

/-- **What is kept is one of the answers.** The editor never invents a document. -/
theorem keep_mem {Doc : Type} (a b : Answer Doc) : keep a b = a ∨ keep a b = b := by
  unfold keep
  split_ifs with h
  · exact Or.inr rfl
  · exact Or.inl rfl

/-- Absorbing an answer into a kept one keeps the larger ticket. -/
theorem keep_keeps_newest {Doc : Type} (a b : Answer Doc) (h : a.ticket ≤ b.ticket) :
    keep a b = b := by
  unfold keep
  exact if_pos h

/-- **The displayed answer is the newest one issued, whatever order the answers arrive in.**

This is the induction `App.tsx`'s ticket counter needs: draining a list of answers in any
order reaches an answer whose ticket is the maximum of the list's, so the document displayed
is the one the last call left. `l.foldr keep x` is "absorb each answer as it lands", with `x`
the answer already on screen. -/
theorem foldr_keep_ticket {Doc : Type} (l : List (Answer Doc)) (x : Answer Doc) :
    (l.foldr keep x).ticket = l.foldr (fun a b => max a.ticket b) x.ticket := by
  induction l with
  | nil => rfl
  | cons a rest ih => rw [List.foldr_cons, keep_ticket, ih, List.foldr_cons]

/-- **And it is one of the answers that arrived** - the fold over a non-empty list is a member
of it. Together with `foldr_keep_ticket` this is the claim: the envelope on screen is the
envelope of the last call, and not a recombination of several. -/
theorem foldr_keep_mem {Doc : Type} (l : List (Answer Doc)) (x : Answer Doc) :
    l.foldr keep x = x ∨ l.foldr keep x ∈ l := by
  induction l with
  | nil => exact Or.inl rfl
  | cons a rest ih =>
    rw [List.foldr_cons]
    rcases keep_mem a (rest.foldr keep x) with h | h
    · rw [h]
      exact Or.inr (by simp)
    · rcases ih with h' | h'
      · rw [h, h']
        exact Or.inl rfl
      · rw [h]
        exact Or.inr (by simp [h'])

/-- **A witness that `keep` is not a projection.** `keep a b = b` where a `keep` that returned
its first argument would answer `a`, so the theorems above are about a real choice rather than
one whose two sides are definitionally equal. -/
theorem keep_keeps_the_newer : keep ⟨1, "first"⟩ ⟨2, "second"⟩ = ⟨2, "second"⟩ := by
  simp [keep]

/-! ## The gesture and the document

One node moves, the canvas shows the gesture while it is in flight, and the release is what
the document is told. -/

/-- A placed node's position on the canvas. -/
structure Pos where
  x : Int
  y : Int
deriving DecidableEq, Repr

/-- The document, as the canvas reads it: a position per node id. -/
abbrev Canvas := String → Pos

/-- The gesture in flight: which node, and where it has been dragged to. -/
structure Gesture where
  node : String
  target : Pos

/-- **What the canvas draws for one node.** A reading, and not a second state. -/
def draw (c : Canvas) (node : String) : Pos := c node

/-- The command a release sends: one node's position, and nothing else. -/
def setPosition (c : Canvas) (node : String) (p : Pos) : Canvas :=
  fun i => if i = node then p else c i

/-- The canvas while the drag is in flight: the gesture's own position for its node. -/
def dragging (c : Canvas) (g : Gesture) : Canvas := setPosition c g.node g.target

/-- **The release is the document being told once**, and it is the same function as the
gesture's own reading - so what a drag costs is a *display*, and what it settles to is a
command. -/
def released (c : Canvas) (g : Gesture) : Canvas := setPosition c g.node g.target

/-- **The position the canvas draws after a release is the position the document carries.**
The claim the middleware page makes about the one piece of local state. -/
theorem drawn_after_release (c : Canvas) (g : Gesture) : draw (released c g) g.node = g.target := by
  unfold draw released setPosition
  exact if_pos rfl

/-- **A drag touches one node.** Every other node draws what it drew before, which is what
"the one piece of local state" means: the gesture is not a second copy of the graph. -/
theorem a_drag_moves_one_node (c : Canvas) (g : Gesture) (i : String) (h : i ≠ g.node) :
    draw (released c g) i = draw c i := by
  unfold draw released setPosition
  exact if_neg h

/-- **And it is a display rather than a document**: while the gesture is in flight the canvas
draws a position the document does not carry, which is why the release has to send one. -/
theorem a_drag_in_flight_is_not_the_document (c : Canvas) (g : Gesture)
    (h : c g.node ≠ g.target) : draw (dragging c g) g.node ≠ draw c g.node := by
  unfold draw dragging setPosition
  rw [if_pos rfl]
  exact h.symm

/-! ## The display unit

`si = (v + offset) * factor`, and `v = si / factor - offset`. -/

/-- What a display unit is to the library: a factor, and an offset for a shifted scale.

**The offset is the whole reason `degC` and `degF` are a table of their own.** A scale is this
with an offset of zero, so one structure covers both and the conversion is one formula. -/
structure Shift where
  factor : ℚ
  offset : ℚ
deriving DecidableEq, Repr

/-- A value in the display unit, in SI. The front end's `value / factor - offset` backwards. -/
def Shift.toSi (s : Shift) (v : ℚ) : ℚ := (v + s.offset) * s.factor

/-- An SI value, in the display unit. -/
def Shift.fromSi (s : Shift) (x : ℚ) : ℚ := x / s.factor - s.offset

/-- **The round trip, in the direction a field shows.** A reader typing a value in the unit a
set names gets back the value they typed. -/
theorem Shift.fromSi_toSi (s : Shift) (h : s.factor ≠ 0) (v : ℚ) :
    s.fromSi (s.toSi v) = v := by
  unfold Shift.fromSi Shift.toSi
  field_simp
  ring

/-- **And in the direction a run's value is displayed.** A quantity the library computed in SI
survives being shown and being read back. -/
theorem Shift.toSi_fromSi (s : Shift) (h : s.factor ≠ 0) (x : ℚ) :
    s.toSi (s.fromSi x) = x := by
  unfold Shift.fromSi Shift.toSi
  field_simp
  ring

/-- The degree Celsius, as the vocabulary's display table carries it: the kelvin's own scale,
shifted by `273.15`. -/
def celsius : Shift := ⟨1, 27315 / 100⟩

/-- The kelvin, which is the same scale with no shift. -/
def kelvin : Shift := ⟨1, 0⟩

/-- The two are different units, which is not a fact about a number. -/
theorem celsius_is_not_kelvin : celsius ≠ kelvin := by
  norm_num [celsius, kelvin]

/-- **One point cannot tell a scale from a shifted one.**

Two units that agree on a single value need not be the same unit - here `⟨2, 0⟩` and `⟨1, 1⟩`
both answer `2` at `1`, and answer `4` and `3` at `2`. This is why
`python/tests/test_units_cross_library.py` asks `pint` for a unit's value at two points: at
one, a scale and a shifted scale are indistinguishable, and a table that named the wrong one
would pass. -/
theorem one_point_does_not_determine_a_shift :
    ∃ s t : Shift, s ≠ t ∧ s.toSi 1 = t.toSi 1 ∧ s.toSi 2 ≠ t.toSi 2 := by
  refine ⟨⟨2, 0⟩, ⟨1, 1⟩, ?_, ?_, ?_⟩ <;> norm_num [Shift.toSi]

/-- **Two points do.** If two shifts agree at two different values then they are one shift,
which is the other half of the two-point check: asking twice is not diligence, it is what
determines the answer. -/
theorem two_points_determine_a_shift (s t : Shift) (hs : s.factor ≠ 0) (v w : ℚ)
    (hvw : v ≠ w) (hv : s.toSi v = t.toSi v) (hw : s.toSi w = t.toSi w) : s = t := by
  have key : (w - v) * (s.factor - t.factor) = 0 := by
    unfold Shift.toSi at hv hw
    nlinarith [hv, hw]
  have hfactor : s.factor = t.factor := by
    have hwv : w - v ≠ 0 := sub_ne_zero.mpr (Ne.symm hvw)
    exact sub_eq_zero.mp ((mul_eq_zero.mp key).resolve_left hwv)
  have hoffset : s.offset = t.offset := by
    have hcancel : (v + s.offset) * s.factor = (v + t.offset) * s.factor := by
      unfold Shift.toSi at hv
      rw [← hfactor] at hv
      exact hv
    have := mul_right_cancel₀ hs hcancel
    linarith
  cases s
  cases t
  simp_all

end View

end Azoth
