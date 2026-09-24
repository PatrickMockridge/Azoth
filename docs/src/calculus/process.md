# Processes and channels

A unit operation is not a tuple of arguments. It is a process with inlets and
outlets, and this page fixes what that means. The claim it exists to support is
that **a unit operation is a pure function** — not as a style rule, but as a
theorem about the process its channels describe.

## Channels, and the type of a channel

A **channel** is a name. A **channel type** is a finite record — a set of field
names, each with a dimension — together with a polarity:

```
Δ  ::=  { f₁ : δ₁, …, fₙ : δₙ }
T  ::=  Δ in   |   Δ out
```

A **stream** is a value of a channel type: for each field, a quantity. Nothing
else crosses a channel. A channel carries its record and no hidden state, which is
what makes the next section work.

**Polarity is an annotation and not an intrinsic.** The pi-calculus has no
direction — a name is a name, and `x!(v)` and `x?(y)` are both uses of it — so a
direction has to be said once, at the port, and then respected by everything
downstream. `separator` has one inlet and two outlets; `mixer` has *many* inlets
and one. That arity is a property of the process, declared per port, and the
declaration is checked against the fields the port names.

**A flux is a flow with a witnessed area.** A molar flow is mol/s; a molar flux is
mol/(m²·s), and the difference is a division by an area. The calculus carries the
area as a field of a channel the process already holds, so a flux cannot be
constructed without the area that makes it one:

```
flux :  n / A        where A is a field of a channel in scope
```

That is the whole content of the distinction, and it is why a flux is a derived
field rather than a unit the vocabulary has to carry separately: it is a quotient
of two things the process has.

## Conservation is linearity

The discipline is that **each inlet channel is consumed exactly once and each
outlet produced exactly once**. Write `⊢ P ok` for a process used linearly.

Under that discipline the balances are not extra checks, they are lemmas:

**Claim (the balances close).** For a unit operation `U` with inlets `c₁ … cₘ` and
outlets `c'₁ … c'ₙ`, with `n(c)` the molar-flow field of a channel:

```
⊢ U ok   ⟹   Σᵢ n(cᵢ)  =  Σⱼ n(c'ⱼ)          moles
             Σᵢ ṅ(cᵢ)  =  Σⱼ ṅ(c'ⱼ)          mass, through the molar masses
             Σᵢ ḣ(cᵢ)  =  Σⱼ ḣ(c'ⱼ) + duty    enthalpy
```

*Status: **specified**. The `ports` layer exists: `specs/unit_ops/` declares 29 unit
operations, each with typed channels carrying a polarity and a multiplicity, and
`azoth_process::validate` holds a flowsheet to the discipline above in Rust — a `one` port
consumed exactly once, a `many` port at least once, every feed and product used once, and
every loop passing through a declared recycle. What does not exist is the theorem:
`lean/Azoth/` has no `Process.lean`, so the balance is *checked* and not *proved*, and
`Azoth.Process.balances_close_under_linearity` is for the tranche that states it.*

The reason this belongs here rather than in a runtime check is that a balance
asserted at runtime is a check that can be skipped, and a balance that follows from
the typing cannot be stated wrongly in the first place. A process whose channels
are not used linearly is a process that does not typecheck, and the balance is what
falls out of that judgement.

## The adequacy claim

**Claim (a unit operation is a pure function).** Let `U` be a process whose only
free names are its declared channels, whose body is a total function `f` of its
inlets, and which performs no other interaction. Then `U` is barbed-congruent to
the pure evaluation of `f`:

```
U  ≅  f(x₁, …, xₘ)
```

*Status: **specified**. The tier this claim is about half-exists: `crates/azoth-process`
carries the channel types, the stream record, the palette loader and the checker,
`specs/unit_ops/` declares 29 unit operations, and thirteen of them carry kernels. So the claim
is no longer that there is nothing to check. It is that **nothing checks it**:
`lean/Azoth/` has no `Process.lean`, and an implementation agreeing with a port declaration
is exactly what this claim asserts and nothing tests.
`Azoth.Process.unit_op_is_extensional` is for the tranche that proves it.*

Three hypotheses, and each is doing work: *only its declared channels*, so `U`
cannot read anything the caller did not supply; *total*, so it cannot fail on a
value its type admits; *no other interaction*, so it cannot observe or change
anything outside. Remove any one and the claim is false, which is the point — the
claim is not that unit operations happen to be pure, it is that **these three
properties are what purity is**, and a unit operation that has them is one.

That is why the specification of a unit operation is its port declaration. The
arity, the fields, and their dimensions are the whole of what an implementation
has to agree with, and the theorem says an implementation that agrees with them is
a pure function of its inlets — in both languages, whatever the code inside looks
like.

## Feedback is restriction

A recycle is not a special construct. It is pi's restriction: a fresh name bound in
a scope and used twice, once written and once read.

```
recycle  ::=  (ν c)( U | V | c̄⟨stream⟩ | c(x). rest )
```

The name `c` is not visible outside the bracket, so the loop cannot be opened from
outside, and the two uses of it are the loop's entry and exit. Nothing else is
needed to express a tear stream, which is a consequence of the calculus rather
than a feature added to it.

**Claim (a recycle has a fixed point).** If the loop's transfer function is a
contraction on the states its channel admits, the recycle process is
barbed-congruent to the unique solution of its loop equations.

*Status: **specified**. A recycle needs restriction and a fixed point, and the
process layer has neither a way to write one nor anything to run one with.
`Azoth.Process.recycle_is_a_fixed_point` is for the tranche that gives it both.*

Uniqueness is the part that matters and the part that is easy to leave out. A
converged loop is only *the* answer when the fixed point is unique; where it is
not, iteration finds one of them and the specification has to say which, which is
the same division `eos.pt_flash` makes when it distinguishes a diagnosed single
phase from a `trivial` convergence it can only admit to. The contraction
hypothesis is what buys uniqueness here.

## What is not in this layer

**There is no interpreter and no flowsheet engine.** Pi and rho are the language
the process layer is *stated* in: a model spec declares its ports, the Lean
development proves the claims above about that declaration once the declaration
exists, and the Rust and Python implementations remain the hand-written kernels
they already are. Building
an executor here would put a second flowsheet layer beside the port, and
[Reflection](./rho.md) is where the reason that is a later question rather than
this one is set out.
