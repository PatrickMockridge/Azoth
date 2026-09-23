# The vocabulary table

The calculus above is hand-written and proved. The vocabulary it is written
against is neither: it is one table of data, compiled into four artefacts. This
page is the specification of that table — the one part of this section whose
authority is a generator rather than a proof.

Unlike the rest of this section, everything on this page is implemented and tested
in this tree today.

## The table

`specs/vocabulary/vocabulary.toml`, validated by
`specs/schema/vocabulary.schema.json`. Three sections:

```toml
slots = ["L", "M", "T", "I", "Th", "N", "J"]

[[dimensions]]
id = "molar_entropy"
exponents = [2, 1, -2, 0, -1, -1, 0]

[[units]]
id = "mm"
dimension = "length"
pint = "millimeter"
uom = "length::millimeter"
rust_ctor = "millimeters"

[[units]]
id = "mol/s"
dimension = "molar_flow"
pint = "mole/second"
```

`slots` is the basis of [the dimension group](./dimensions.md), in the order every
exponent vector is written in — `uom::si::ISQ`'s order, and the order is not
cosmetic: everything generated writes its exponents this way, so a different order
permutes all seven coordinates at once.

A unit declares three things and nothing else: **its dimension**, **its name in
`pint`**, and **its conversion path in `uom`**. When the path is present it also
names the hand-written constructor in `crates/azoth-core/src/units.rs` that the
generated conversion calls; the generator never invents a constructor name, and a
row that claims a path and names none is refused. A unit with **no `uom` key** is one
whose dimension `uom` does not carry, and its conversion is the identity.

## The rule: no magnitude is written down

**No conversion factor appears in the table, or in any generated file, or anywhere
else in this repository.**

A factor is a number `uom` and `pint` each already know. A number written in a
third place is a number that can disagree with both of them, and it is not
hypothetical here: a version of `crates/azoth-core/src/units.rs` paired each
canonical unit with a hand-typed SI magnitude, and `mm`'s said `1.0e-3` while the
surrounding machinery treated it as an SI base unit — a factor of a thousand, on
one side only, squared by the bore diameter of an orifice. `batch/_core.py`
records the same lesson for the Python side.

So the table declares what a unit **is** and never what it is **worth**. The
generated Rust carries the conversion a calculation runs and no expected value
beside it.

## `uom: null`, and why it is a value rather than a gap

Some dimensions have no `uom` quantity, and the table says so with a null rather
than with a comment. This is structural, not an omission: `uom`'s `quantity!` macro
must run inside `uom::si`, and the `Units` trait that makes `.new::<u>()` compile
is assembled by `system!` from a closed list. **No crate outside `uom` can add a
quantity to `uom::si`**, so a dimension `uom` does not carry is one azoth has to
name itself.

Rows with a null path, and the generated conversion for each is the identity
because the unit is already its own SI base unit:

| Unit | Dimension | Why `uom` has no quantity |
|---|---|---|
| `mol/s` | `N·T⁻¹` | `MolarFlux` is mol/(m²·s) and `MolarConcentration` is mol/m³; neither is this |
| `J/(mol*K**2)` … `J/(mol*K**5)` | `M·L²·T⁻²·Θ⁻ⁿ·N⁻¹` | the Cp polynomial coefficients; `uom` carries no fractional power of the temperature unit |
| `1/K` … `1/K**3` | `Θ⁻¹` … `Θ⁻³` | the dielectric-constant coefficients `ε(T) = d₀ + d₁/T + d₂T + d₃T² + d₄T³`; `uom` carries a temperature and not its reciprocal |

A null is therefore the *list of dimensions azoth names itself*, stated rather
than left for a reader to discover — and a row that claims a path for one of them
is refused at generation rather than failing later as an obscure type error.

## What is generated, and what is not

`tools/gen_vocabulary.py` compiles the table into four artefacts, and `--check`
fails if any of them is stale:

| Artefact | Is |
|---|---|
| `specs/schema/unit.schema.json` | the closed enum a spec's `unit:` must name |
| `crates/azoth-core/src/unit_vocab_gen.rs` | `UNIT_NAMES`, the dimensions, the conversions |
| `python/src/azoth/core/_units_gen.py` | the same, as the `pint` map and the slot list |
| `azoth.keycard.UNIT_VOCABULARY` | what the keycard loader admits |

**The split is the same one the port already makes.** Two kernels are hand-written
and only two, and everything that merely *names* them is generated. Here the
hand-written things are the Lean development and the boundary constructors in
`units.rs` — which are prose, and stay prose — and the generated thing is the
vocabulary that names them.

The **`uom` quantity type is derived from the exponents**, through a map the
generator holds and refuses to be incomplete. That is what makes the generated
assertion bite: `crates/azoth-core/src/unit_vocab_gen.rs` contains, for each unit,
a line saying the constructor produces the quantity *this unit's exponents* name —
so declaring `mm` with the wrong dimension fails to compile rather than being waved
through by a name that happens to match.

## The checks, and what each one can see

Four, and they are complementary rather than redundant. Each was verified by
breaking the thing it guards.

| Check | Where | Catches |
|---|---|---|
| Generator consistency | `gen_vocabulary.py` | a dimension and a `uom` path that disagree with each other |
| Compile-time dimension | `unit_vocab_gen.rs`, `the_table_agrees_with_uom` | a constructor that produces a different quantity than the exponents name |
| `pint` dimensionality | `python/tests/test_units_cross_library.py` | table exponents that `pint` disagrees with |
| Factor agreement | the same file | the two units libraries disagreeing about a factor |
| Lean dimension | `lean/Azoth/Vocabulary.lean`, one theorem per unit | a unit whose exponents name a different dimension than `lean-units` does |

The factor check is the one that was missing. It compares `pint`'s own answer for a
unit name against `azoth._core.unit_si_factor`, which runs the very conversion a
calculation runs — so **neither side is a literal**. Pointing `mm` at `pint: meter`
reports *"pint says one of it is 1.0 in SI base, and the Rust conversion this
library actually runs gives 0.001"*, which is the defect above, caught by a test
rather than by a reviewer.

## The Lean leg, and why its right-hand side is written by hand

`lean/Azoth/Vocabulary.lean` carries one theorem per unit:

```lean
theorem u_mm_dimension :
    dimOf "mm" = some (Dimension.Length) := by
  simp only [dimOf, units, ...] <;> simp <;> module
```

The left is the table's exponents, read through `Azoth.Dim.slots`. The right is a
dimension with a name **`lean-units` wrote** — and `LEAN_DIMENSIONS` in
`tools/gen_vocabulary.py` is what decides which name each unit gets. That map is
keyed by **unit**, not by dimension, and that is the whole reason the theorem is a
check:

- Keyed by dimension it would prove nothing. A table entry whose exponents were
  wrong would select the expression for those wrong exponents, and the theorem would
  hold however wrong the table was.
- Keyed by unit it is a hand-written claim about what each unit *is* — the
  counterpart of `rust_ctor`, which is what makes the Rust assertion bite — and the
  theorem forces the table to agree with it.

So declaring `mm` an area fails to prove, with a goal of `2 = 1`: the table's second
exponent against the first one `Length` has. A unit with no entry in the map is
refused at generation rather than skipped, so the map stays total as the table
grows. And because the map is data, `Axioms.lean`'s sibling `Gate.lean` is generated
from it too — a hand-maintained list of thirty-seven names would go stale the first
time a unit was added, and go stale silently, with the theorem proved and nothing
gating it.

## Adding a unit

Add a row to the table and run `python tools/gen_vocabulary.py`. If the dimension
is new, the generator stops and asks whether `uom` carries it — and "it does not"
is an answer, spelled `null`. If `uom` does carry it, add a constructor to
`crates/azoth-core/src/units.rs` and name it in the row.

Nothing else changes: the schema, the Rust list, the Python map and the keycard
loader are all compiled from the one table, and `gen_vocabulary.py --check` runs in
both the spec-validation and the drift jobs, so an edit that was not regenerated
fails a build rather than a calculation.
