# The databank

Everything azoth ships as *input*: the upstream data it vendors, the record of what it
took from each and why it took no more, and in time the baseline keycard both derive
into and the compiled files that come out of it.

Specified in [the specification](../docs/src/architecture/specification.md#the-keycard-is-where-responsibility-sits).
This page is the operational half: what is here now, what is not, and how to check it.

## The stages, and which of them exist

```
databank/sources/          upstream files, whole, at a named revision      EXISTS
        |  compile             tools/gen_databank.py
        v
data/components/           the files both languages read                   EXISTS
data/reactors/             GibbsReactor's species database                 EXISTS
        |  derive              tools/gen_keycard.py
        v
databank/keycard.toml      the baseline card: the subset a user overrides  EXISTS
```

The order is one-way: the compiled files come from the sources, and the baseline card
comes from the compiled files. The card is the **subset a user may override** - the
cubic's `Tc`, `Pc` and `omega` and the Peng-Robinson `kij` - stated in the format a
user's card takes, so the library's data and a user's data are one object. It is not
the source of the compiled files, which carry every column the manifest dispositions;
`databank/compiled/` does not exist because the compiled files live at `data/`, where <!-- doc-claims-ok: the sentence says this path is absent -->
the two languages read them at fixed paths - Rust embeds them with `include_str!`, and
`python/src/azoth/_data.py` finds them by walking up from the package.

## `sources/` holds upstream files whole

`databank/sources/neqsim/COMP.csv` and `INTER.csv` are NeqSim's own bytes. Not a
slice, not a transcription: byte-identical to the release named in the manifest, which
is what lets `tools/gen_databank.py --check` reproduce the shipped files without a
NeqSim checkout. Refreshing them is what a checkout is for.

NeqSim is Apache-2.0; the attribution is in [`NOTICE`](../NOTICE).

## `manifest.toml` is the record of what was taken

Every resource file NeqSim ships is declared once, and every column of every one of
them whose columns are enumerated — 33 files, 1,526 columns — is listed with what was
done with it and a reason; the other six are vendored whole, three of them because
they are not tabular at all. `tools/check_manifest.py` prints the tally:

```
check_manifest: OK (37 vendored file(s), 1526 column(s), 1502 carried of which 1185 read, 0 not-vendored entr(ies))
  317  carried, nothing reads it yet
  1090  carried with no unit NeqSim states (neqsim-internal)
  276  not-ported
   19  not-a-value
   10  empty-upstream
    1  superseded-by
    7  unreachable-upstream
   14  uncalled-upstream
   11  unread-upstream
```

**The porting backlog is "carried, nothing reads it yet" — 317 columns — and it is
the number that matters.** NeqSim is the target, not a reference: each carried column is a
physical property whose model NeqSim implements and azoth has not ported, and the
`not-ported` reason names the class that would close it — **276 columns carry it**, carried
and dropped together, with a handful `unreachable-upstream`, fourteen `uncalled-upstream`, eleven
`unread-upstream` and the rest `not-a-value` or `empty-upstream` — `PhaseHydrate`,
`CPAMixingRuleHandler`, `SolidFlash1`, `PhasePCSAFTa`, `ParachorSurfaceTension` and the rest.
The check refuses a `not-ported` reason with no NeqSim name in it, so the list cannot drift
back into being somewhere to put a column.

Nothing here is "out of scope". That word was in an earlier draft of this vocabulary
and it was wrong: a file it labelled out of scope was work not yet done, and filing it
as a decision made an incomplete port look like a boundary.

`not-yet` is the smaller and nearer list: the model exists here and the data does not
reach it. Those entries name a registered id that would read the column once it does.

`unreachable-upstream` is the one that is not work at all. The class named would read
the column, and **nothing in NeqSim constructs it** — measured with a pattern that admits
a fully qualified `new pkg.Class(`, because the bare `new X(` grep is not a liveness test
(`ThermodynamicOperations.java:190` writes one that way). `SystemPrMathiasCopeman` is the
example: the only class that mentions the Mathias-Copeman alpha, never constructed, and
`mcpr1` appears in no source file at all.

`uncalled-upstream` is the same finding one level down, and it needs its own word because
the measurement is different: the class **is** live, and the member that would read the
column is not. `ComponentSolid.fugcoef(PhaseInterface)` has its solid-vapour-pressure branch
commented out and returns `fugcoef2`; the overload that reads the column,
`ComponentSolid.fugcoef(double, double)`, is called by nothing; and the one other member that
reaches the same columns, `ComponentHydrate.getEmptyHydrateStructureVapourPressure`'s
`type == -1` branch, is settable from nowhere in `src/main`. A reachability check on the
class answers yes, which is why the class-level test does not cover this case: the ten
columns that carry the kind are read by a route that exists and cannot be walked.

`unread-upstream` is the third of the family and the emptiest: there is no reader to test.
Not a class nobody constructs and not a member nobody calls — **no reader at all**.
`ReactionKSPdata.csv` is the example, and upstream says it in as many words, at
`NeqSimDataBase.java:597`, on the line above the `updateTable` call that loads it:
`// Table ReactionKSPdata is not in use anywhere`. Its only other reference is
`DataCatalogRunner.addTable`, which lists it for the MCP data catalogue and reads none of
it. Two of the family's three names would be false here — there is no class to find
unconstructed and no member to find uncalled — which is why it has its own word rather than
a widened one. A reachability check cannot see this kind either: it searches for a reader,
and the search that comes back empty is the answer.

A reason is a quoted flow mapping, so `grep 'reason: not-ported'` matches nothing.
Read the tally the check prints; there is no grep for it.

A reader asking "is this databank thin because we decided, or because nobody looked?"
should be able to answer it from this file alone. That is the whole of its purpose.

## Checking it

```bash
python tools/check_manifest.py                     # the manifest against the files
python tools/gen_databank.py --check               # the component files against the sources
python tools/gen_reaction_data.py --check          # the reaction files against the sources
```

All three run in CI. Between them: every vendored file is declared and present; the
manifest's columns and the source's header agree in both directions; the `used` columns
and the compiled header agree in both directions; the row counts are what the manifest
says; and each generator's own column list cannot disagree with the manifest without
failing before it writes anything.

**Two generators rather than one** because the compiled files answer to two crates:
`data/components/` is the equation-of-state data and `data/reactions/` is the reaction
data, and a tool that wrote both would be the one place a namespace's data could be
generated into the wrong tree.

## What none of this can tell you

**Whether the vendored slice is current.** `manifest.toml` records the NeqSim commit it
was last checked against (`f0c7436c6923766b1e22957b7075f650600457a7`, master), and
nothing here can tell you a newer NeqSim exists. The `not_vendored` list was built by
expanding NeqSim's resource directory by hand, at that commit.

**The `version` field does not identify the revision.** It names the `revision` property in
NeqSim's own `pom.xml`, which reads `3.21.0` both at the `v3.21.0` tag and on master — 50
commits apart — so the same version names two different trees. **The `commit` is the pin**,
and a document wanting a readable name for it says `NeqSim master`, because that is the only
name that resolves to the bytes vendored here.

The alternative - fetching from GitHub in CI - would make every build depend on a
network service the project does not control, to answer a question that only matters
when someone deliberately decides to re-vendor. So it is a human step, and the
staleness is visible in the file rather than hidden.

## Moving the pin

Moving it is a deliberate act with a procedure, because the two things it invalidates are
invisible from here: the vendored columns, and the oracle numbers hard-coded into
`crates/*/tests/*.rs`. Those tests are the only thing that would notice a moved number, and
their bars are wide enough that a move of `1e-9` passes them - which is how two whole test
files stayed on the previous revision through the last refresh.

1. **Build the jar** from the checkout, as `.gitignore` says, and name it for its commit.
2. **`tools/pin_impact.py --from <old> --to <new>`** - which NeqSim classes changed, which of
   them this tree cites, and which of those moved a *numeric literal* rather than a brace.
   It needs no jar and no probe run, and it names the files to re-measure. The classes it
   reports differ from a plain `git diff` because a formatting sweep reaches nothing.
3. **`tools/oracle_sweep.py`** - re-runs every probe, rebuilds every capture, and reports
   both what moved and any crate-test literal still holding the *previous* value. It exits
   non-zero on a stale literal. Run it before regenerating the captures: the committed
   capture *is* the previous pin, and that is what makes the check decidable with one jar.
4. **Re-vendor and regenerate**, as the sections above describe, and re-capture.

**The models with no probe are the residual gap.** `crates/*/tests/` holds reference-equation
tests - ammonia, Vega, Span-Wagner, Leachman, argon and para-hydrogen solid - whose NeqSim
classes nothing in `validation/neqsim/` drives, so a move there is caught by step 2 and not
by step 3. In the last window six of those classes changed and exactly one of them,
`Ammonia2023`, mattered: it gained two methods and the five others changed only in
`64efef1 Spotless`, which step 2 says by reporting their literals as identical.
