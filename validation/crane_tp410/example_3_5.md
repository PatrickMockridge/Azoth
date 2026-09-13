# crane_tp410/example_3_5

A validation case for `hydraulics.darcy_weisbach`.

**This case does not validate against a verified source.** Read
[the attribution](#the-attribution-problem) before treating it as evidence of
anything.

## The case

| Input | Value | Unit |
|---|---|---|
| `f` | 0.02 | dimensionless |
| `L` | 100 | m |
| `D` | 0.1 | m |
| `rho` | 998 | kg/m³ |
| `v` | 1.5 | m/s |

Expected `dp` = **22455.0 Pa**, tolerance `1e-9` relative.

## The derivation

$$
\Delta P = f \cdot \frac{L}{D} \cdot \frac{\rho v^2}{2}
$$

```
L/D        = 100 / 0.1          = 1000
rho*v**2/2 = 998 * 1.5**2 / 2   = 1122.75
dP         = 0.02 * 1000 * 1122.75
           = 22455.0 Pa
```

Checkable by hand, which is the point: the arithmetic is ours and needs no
source.

## The attribution problem

The project brief attributed this case to **Crane TP-410 Example 3-5** and
supplied `expected: {dP: 2245.5}`.

Neither holds up.

**The value is wrong by exactly a factor of ten.** The brief's own equation, applied
to the brief's own inputs, gives 22455.0 Pa. And 2245.5 is not a rounding of that
- it is exactly `rho * v**2`, the dynamic pressure, which appears nowhere in the
Darcy-Weisbach equation. (The velocity head `rho*v**2/2` is 1122.75; the pressure
drop is that times `f*L/D` = 20.) A value that equals an unrelated expression in
the same problem looks like a transcription slip rather than a different
convention - somebody dropped the `f * L/D` factor and doubled what remained.

**The source could not be confirmed.** The example could not be located, and
Crane TP-410's publicly documented numbering for this relation is Eq. 1-6 / 1-7
(and the resistance-of-bends equation 2-20), not 3-2 as the brief also stated.

So the case is recorded as `source_needed`: the inputs are kept, because they are
a perfectly good case, and the expected value is **derived from the equation**
rather than taken from the brief.

## Why this is not simply marked "verified"

There is no path to `verified` here that does not involve a person with a
legitimate copy of Crane working the example themselves. Reproducing a worked
example out of the standard is exactly what the project's copyright rule forbids,
and inventing a citation for a value we computed ourselves would be worse than
recording the gap.

## What this case does establish

That the implementation computes what the equation says, for a set of inputs that
look like a real water line. The arithmetic is independently checkable, and the
same inputs appear as the spec's own worked example, so a change that broke the
equation would fail in two places.

## What it does not establish

Anything about Crane TP-410, or that these inputs are inside the range the
equation is validated over. At `rho = 998`, `v = 1.5`, `D = 0.1` and
`mu = 1.002e-3` the Reynolds number is about 149400 - turbulent, which is fine -
but this case does not supply a viscosity, so nothing here checks the regime.
