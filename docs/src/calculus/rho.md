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

*Status: not yet proved. Will be `Azoth.Rho.reflection_is_a_bijection`.*

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

*Status: not yet proved. Will be `Azoth.Rho.reflected_feedback_is_congruent`.*

## What this layer is for, and what it is not yet

Reflection is where three things that look unrelated turn out to be one:

- **Feedback.** A recycle with a parameterised body, above.
- **Interoperation.** A process as a value, so a description of a plant is data.
  Neither the current code nor any tool this library talks to has this yet.
- **Extension.** A caller supplying a unit operation this library does not
  implement, as a process rather than as a new kernel in two languages.

None of the three is implemented, and this page does not claim otherwise. What it
fixes is that they are the *same* problem, and that the calculus already has the
operator for it — so a future flowsheet layer is a use of this, rather than a
parallel language invented beside it.

**There is deliberately no executor here.** Turning a process into something a
machine runs is a flowsheet engine, which is tranche P12 and has to arrive after
the unit operations it would compose. The claim this page exists for is the round
trip, and the round trip is a theorem before it is a feature.
