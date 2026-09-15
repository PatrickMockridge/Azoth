/-
The process calculus of the rho layer, as `docs/src/calculus/rho.md` states it.

The calculus has two kinds of term, mutually defined: a `Name` is a de Bruijn index
or a process quoted as data (`@P`), and a `Process` is the inactive process, a
parallel composition, a name sent or received on a channel, or a name dropped back
into a process (`*x`). Reflection is the round trip between the two - `drop (quote P)`
reduces to `P` - and communication is the higher-order send/receive, realised by the
capture-avoiding substitution below. The barbed-congruence reading of
`reflection_is_a_bijection` is built on this in `Azoth.Barb`.
-/

namespace Azoth

namespace Rho

mutual
  /-- A name is a de Bruijn index, or a process quoted as data (`@P`). -/
  inductive Name where
    | free : Nat → Name
    | quote : Process → Name

  /-- A process is the inactive process, a parallel composition, a name sent or
  received on a channel, or a name dropped back into a process (`*x`).

  `out a b P` sends the name `b` on channel `a` and continues as `P`; `in a P`
  receives on channel `a` and continues as `P` with the received name bound as the
  de Bruijn index 0. -/
  inductive Process where
    | nil : Process
    | par : Process → Process → Process
    | out : Name → Name → Process → Process
    | in : Name → Process → Process
    | drop : Name → Process
end

mutual
  /-- Shift the free indices at or above `c` up by one, making room for a new binder. -/
  def liftName (c : Nat) : Name → Name
    | .free n => .free (if c ≤ n then n + 1 else n)
    | .quote P => .quote (lift c P)

  def lift (c : Nat) : Process → Process
    | .nil => .nil
    | .par P Q => .par (lift c P) (lift c Q)
    | .out a b P => .out (liftName c a) (liftName c b) (lift c P)
    | .in a P => .in (liftName c a) (lift (c + 1) P)
    | .drop x => .drop (liftName c x)
end

/-- The lift of a substitution: index 0 stays 0 (the binder's own name), and index
`n + 1` is `σ n` shifted up by one. -/
def upSubst (σ : Nat → Name) : Nat → Name
  | 0 => .free 0
  | n + 1 => liftName 0 (σ n)

mutual
  /-- Apply a simultaneous substitution to a name or a process. -/
  def substName (σ : Nat → Name) : Name → Name
    | .free n => σ n
    | .quote P => .quote (subst σ P)

  def subst (σ : Nat → Name) : Process → Process
    | .nil => .nil
    | .par P Q => .par (subst σ P) (subst σ Q)
    | .out a b P => .out (substName σ a) (substName σ b) (subst σ P)
    | .in a P => .in (substName σ a) (subst (upSubst σ) P)
    | .drop x => .drop (substName σ x)
end

/-- The substitution that puts `b` in place of index 0 and decrements the rest. -/
def singleSubst (b : Name) : Nat → Name
  | 0 => b
  | n + 1 => .free n

/-- Substitute the single name `b` for de Bruijn index 0. -/
def subst0 (b : Name) (P : Process) : Process := subst (singleSubst b) P

/-- Substitute the single name `b` for de Bruijn index 0 in a name. -/
def subst0Name (b : Name) (x : Name) : Name := substName (singleSubst b) x

/-- One-step reduction: the round trip and communication, closed under parallel
composition.

`par_left` and `par_right` say the two sides of a parallel composition reduce
independently, which is what makes the calculus concurrent rather than sequential.
-/
inductive Reduces : Process → Process → Prop where
  | drop_quote (P : Process) : Reduces (Process.drop (Name.quote P)) P
  | comm (a b : Name) (P Q : Process) :
      Reduces (Process.par (Process.out a b P) (Process.in a Q))
              (Process.par P (subst0 b Q))
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

/-- The process carried by a send is determined by the send. -/
theorem out_injective {a b : Name} {P Q : Process}
    (h : Process.out a b P = Process.out a b Q) : P = Q := by
  injection h

/-- The continuation of a receive is determined by the receive. -/
theorem in_injective {a : Name} {P Q : Process}
    (h : Process.in a P = Process.in a Q) : P = Q := by
  injection h

/-- Substituting for index 0 puts the name in place of the bound name. -/
theorem subst0_var_zero (b : Name) : subst0Name b (.free 0) = b := by
  rfl

/-- Substituting decrements the indices above 0. -/
theorem subst0_var_succ (b : Name) (n : Nat) : subst0Name b (.free (n + 1)) = .free n := by
  rfl

/-- Under a binder the outer index 0 becomes the lifted name, so the substitution
does not capture the binder's own name. -/
theorem subst0_under_binder (b c : Name) :
    subst0 b (.in c (.drop (.free 1))) = .in (subst0Name b c) (.drop (liftName 0 b)) := by
  rfl

/-- A communication: send a name, receive it, and drop it. -/
theorem comm_reduction (a b : Name) :
    Reduces (.par (.out a b .nil) (.in a (.drop (.free 0)))) (.par .nil (.drop b)) := by
  simpa using (Reduces.comm a b .nil (.drop (.free 0)))

end Rho

end Azoth
