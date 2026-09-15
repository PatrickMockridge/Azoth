/-
The reflection fragment of the rho-calculus, as `docs/src/calculus/rho.md` states it.

Rho adds two operators to the process calculus - `@P`, quoting a process as a name,
and `*x`, dropping a name back into a process - and the claim it exists for is the
round trip: `*@P ≅ P` and `@*x ≅ x`, so a process carried on a channel and re-entered
is the process that was carried. That is what makes an interoperation surface
possible, and the reason this module exists.

This formalises the *binder-free* fragment of reflection: the name and process terms,
and the round-trip reduction `drop (quote P) → P`. The input, output and restriction
prefixes of the pi-calculus, and therefore the barbed-congruence reading of
`reflection_is_a_bijection` in `rho.md`, belong to `Azoth.Process`, which does not
exist yet. The two injectivity theorems below are the part of "reflection is a
bijection" that a fragment with no communication can already state: quoting and
dropping do not collapse two distinct things into one.
-/

namespace Azoth

namespace Rho

mutual
  /-- A name is a free index, or a process quoted as data (`@P`). -/
  inductive Name where
    | free : Nat → Name
    | quote : Process → Name

  /-- A process is the inactive process, a parallel composition, or a name dropped
  back into a process (`*x`). -/
  inductive Process where
    | nil : Process
    | par : Process → Process → Process
    | drop : Name → Process
end

/-- One-step reduction: the round trip, closed under parallel composition.

`par_left` and `par_right` say the two sides of a parallel composition reduce
independently, which is what makes the calculus concurrent rather than sequential.
-/
inductive Reduces : Process → Process → Prop where
  | drop_quote (P : Process) : Reduces (Process.drop (Name.quote P)) P
  | par_left {P P' Q : Process} : Reduces P P' → Reduces (Process.par P Q) (Process.par P' Q)
  | par_right {P Q Q' : Process} : Reduces Q Q' → Reduces (Process.par P Q) (Process.par P Q')

/-- Quoting is injective: two processes that quote to the same name are equal. -/
theorem quote_injective {P Q : Process} (h : Name.quote P = Name.quote Q) : P = Q := by
  injection h

/-- Dropping is injective: two names that drop to the same process are equal. -/
theorem drop_injective {x y : Name} (h : Process.drop x = Process.drop y) : x = y := by
  injection h

/-- The round trip: a process quoted and dropped reduces back to itself. -/
theorem round_trip (P : Process) : Reduces (Process.drop (Name.quote P)) P :=
  Reduces.drop_quote P

/-- The size of a process, the measure the reflection reduction is bounded by. -/
def size : Process → Nat
  | .nil => 1
  | .par P Q => 1 + size P + size Q
  | .drop (.free _) => 2
  | .drop (.quote P) => 2 + size P

/-- The reflection reduction strictly decreases size, so the round trip terminates. -/
theorem reduction_decreases_size {P Q : Process} (h : Reduces P Q) : size Q < size P := by
  induction h with
  | drop_quote P =>
      simp [size]
  | par_left _ ih =>
      simp [size]
      exact ih
  | par_right _ ih =>
      simp [size]
      exact ih

end Rho

end Azoth
