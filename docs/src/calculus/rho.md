# Reflection and feedback

[Processes](./process.md) states a unit operation and a recycle in the
pi-calculus. Two things it cannot state are: a process that is *named* rather than
written out, and a loop whose body is not known until it is passed in. Both are
what the rho-calculus adds, and both are what an interoperation surface needs.

## Reflection and drop

Rho adds two operators and two reduction rules:

```
@P        quote: a process, as a name
*x        drop:  a name, as a process

x!(P) | x?(y).Q   →   Q{P/y}          higher-order output
*x    | @P        →   P | x           drop and quote cancel
```

The first rule is what makes the calculus higher-order: what travels on a channel
can be a *process*, not only a value. The second is what makes the two operators
inverses.

**Claim (serialisation round-trips).** `*@P ≅ P` and `@*x ≅ x`, structurally —
so quoting a process and dropping it back is the identity, and a process carried
on a channel and re-entered is the process that was carried.

*Status: the reflection and communication fragment, and the round trip **as
structural congruence**, are **proved**; the barbed-congruence reading is still
**specified**. `Azoth/Rho.lean` formalises the operators — `@P`, `*x`, and the
communication prefixes `out`/`in` — and proves the round-trip *reduction*
`drop (quote P) → P` (`Azoth.Rho.round_trip`), the injectivity of quoting, dropping,
sending and receiving (`Azoth.Rho.quote_injective`, `drop_injective`,
`out_injective`, `in_injective`), the capture-avoiding substitution and its
round-trip (`Azoth.Rho.subst0_var_zero`, `subst0_var_succ`, `subst0_under_binder`),
the communication reduction (`Azoth.Rho.comm_reduction`), and the round trip as `≅`
(`Azoth.Rho.round_trip_congr` — the `*@P ≅ P` half of `reflection_is_a_bijection`;
`@*x ≅ x` is the same theorem at `P = drop x`). What remains *specified* is the
stronger barbed-bisimulation reading, a coinductive layer.*

This is the claim that makes an interoperation surface possible rather than
aspirational. A flowsheet handed to another tool, or stored in a file, or drawn in
an editor, is a process turned into data and turned back — and the theorem is that
this round trip loses nothing. Without it, "a GUI describes a plant" is a claim
about a serialiser; with it, it is a claim about the calculus.

## Why feedback wants reflection

A recycle is expressible in pi alone, as
[Processes](./process.md) shows: bind a fresh name, write it once, read it once.
What pi cannot express is a loop whose *body* is a parameter — a tear stream whose
unit operations are supplied later, or a module applied to a loop. That is
reflection's job, because a process that arrives on a channel is a name and a name
can be re-entered as a process:

```
module(x)  ::=  *x | @( recycle-in-x )
```

So a module is a process that takes a loop and returns a loop, and the loop it
returns is the one it was given with its own body wrapped around it. That is what
makes a flowsheet composable at the level of *modules* rather than only at the
level of streams.

**Claim (a reflected loop has the same fixed point).** A recycle whose loop channel
carries its body by reflection is barbed-congruent to the recycle written out
directly, so reflection changes how a loop is *named* and not what it computes.

*Status: **specified**, for the same reason: it is about a loop expressed by
reflection, and there is neither.*

## What this layer is for

Reflection is where three things that look unrelated turn out to be one:

- **Interoperation.** A process as a value, so a description of a plant is data.
- **Feedback.** A recycle with a parameterised body, above.
- **Extension.** A caller supplying a unit operation this library does not
  implement, as a process rather than as a new kernel in two languages.

**The first is built, and it is this operator made concrete.** A flowsheet is a
process; quoting it and dropping it back is the identity; and that round trip is
what an editor, a file and a notebook each do to the same document. There is one
interoperation layer — `crates/azoth-process`'s `middleware` module — and three
bindings of it, for a browser, a notebook and a shell, none of which holds a
document of its own. **The executor is built as well, and it is a use of this
rather than a layer beside it**: [the middleware](../architecture/middleware.md)
states how, and this page does not restate it.

**The other two are not built**, and what the calculus contributes is that they are
the *same* problem, with the operator for it already written. A loop whose body
arrives as a process, and a unit operation the caller supplies as a process, are
both **specified**: the layer's operators exist and the things they would operate
on do not. What is not built is named at the foot of
[the middleware](../architecture/middleware.md) rather than copied here.

The round trip is what makes a front-end a *reading* of one artifact instead of a
translation of it, and it is why the three bindings were built on this operator
rather than each inventing a serialiser. The claim this page exists for is the
round trip, and the round trip is a theorem before it is a feature.
