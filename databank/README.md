# The databank

Everything azoth ships as *input*: the upstream data it vendors, the record of what it
took from each and why it took no more, and in time the baseline keycard both derive
into and the compiled files that come out of it.

Specified in [the specification, S5](../docs/src/architecture/specification.md#s5-it-ships-data-and-the-keycard-extends-it).
This page is the operational half: what is here now, what is not, and how to check it.

## The stages, and which of them exist

```
databank/sources/          upstream files, whole, at a named revision      EXISTS
        |  compile             tools/gen_databank.py
        v
data/components/           the files both languages read                   EXISTS
        |  derive              tools/gen_keycard.py
        v
databank/keycard.toml      the baseline card: the subset a user overrides  EXISTS
```

The order is one-way: the compiled files come from the sources, and the baseline card
comes from the compiled files. The card is the **subset a user may override** - the
cubic's `Tc`, `Pc` and `omega` and the Peng-Robinson `kij` - stated in the format a
user's card takes, so the library's data and a user's data are one object. It is not
the source of the compiled files, which carry every column the manifest dispositions;
`databank/compiled/` does not exist because the compiled files live at `data/`, where
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
them (35 files, 1,496 columns) is listed with what was done with it and a reason.
`tools/check_manifest.py` prints the tally:

```
check_manifest: OK (35 vendored file(s), 1496 column(s), 1459 carried of which 14 read, 0 not-vendored entr(ies))
  1445  carried, nothing reads it yet
  1317  carried with no unit NeqSim states (neqsim-internal)
   10  not-a-value
   22  empty-upstream
    5  superseded-by
```

**The porting backlog is "carried, nothing reads it yet" — 1,445 columns — and it is the
number that matters.** NeqSim is the target, not a reference: each carried column is a
physical property whose model NeqSim implements and azoth has not ported, and the
`not-ported` reason names the class that would close it — `PhaseHydrate`,
`CPAMixingRuleHandler`, `SolidFlash1`, `PhasePCSAFTa`, `ParachorSurfaceTension` and the
rest. The check refuses a `not-ported` reason with no NeqSim name in it, so the list
cannot drift back into being somewhere to put a column.

Nothing here is "out of scope". That word was in an earlier draft of this vocabulary
and it was wrong: a file it labelled out of scope was work not yet done, and filing it
as a decision made an incomplete port look like a boundary.

`not-yet` is the smaller and nearer list: the model exists here and the data does not
reach it. Those entries name a registered id that would read the column once it does.

A reason is a quoted flow mapping, so `grep 'reason: not-ported'` matches nothing.
Read the tally the check prints; there is no grep for it.

A reader asking "is this databank thin because we decided, or because nobody looked?"
should be able to answer it from this file alone. That is the whole of its purpose.

## Checking it

```bash
python tools/check_manifest.py                     # the manifest against the files
python tools/gen_databank.py --check               # the shipped files against the sources
```

Both run in CI. Between them: every vendored file is declared and present; the manifest's
columns and the source's header agree in both directions; the `used` columns and the
compiled header agree in both directions; the row counts are what the manifest says;
and `tools/gen_databank.py`'s own column list cannot disagree with the manifest without
failing before it writes anything.

## What none of this can tell you

**Whether the vendored slice is current.** `manifest.toml` records the NeqSim commit it
was last checked against (`dedba8735d030c6411e09b6fd7e69f6c4a136114`, v3.20.0), and
nothing here can tell you a newer NeqSim exists. The `not_vendored` list was built by
expanding NeqSim's resource directory by hand, at that commit.

The alternative - fetching from GitHub in CI - would make every build depend on a
network service the project does not control, to answer a question that only matters
when someone deliberately decides to re-vendor. So it is a human step, and the
staleness is visible in the file rather than hidden.
