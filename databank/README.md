# The databank

Everything azoth ships as *input*: the upstream data it vendors, the record of what it
took from each and why it took no more, and in time the baseline keycard both derive
into and the compiled files that come out of it.

Specified in [the specification, S5](../docs/src/spec.md#s5-it-ships-data-and-the-keycard-extends-it).
This page is the operational half: what is here now, what is not, and how to check it.

## The three stages, and which of them exist

```
databank/sources/          upstream files, whole, at a named revision      EXISTS
        |
        |  derive
        v
databank/keycard.toml      the baseline keycard: what azoth ships           NOT BUILT
        |
        |  compile
        v
databank/compiled/         the files both languages actually read          NOT BUILT
```

The order is one-way. A source produces a card, and a card produces the compiled files.

**Only the first stage is built.** The compiled files still live at `data/`, where the
two languages read them at fixed paths - Rust embeds them with `include_str!`, and
`python/src/azoth/_data.py` finds them by walking up from the package. Nothing here
has moved, so nothing here is broken by the move taking a while.

## `sources/` holds upstream files whole

`databank/sources/neqsim/COMP.csv` and `INTER.csv` are NeqSim's own bytes. Not a
slice, not a transcription: byte-identical to the release named in the manifest, which
is what lets `tools/gen_databank.py --check` reproduce the shipped files without a
NeqSim checkout. Refreshing them is what a checkout is for.

NeqSim is Apache-2.0; the attribution is in [`NOTICE`](../NOTICE).

## `manifest.yaml` is the record of what was taken

Every column of `COMP.csv` (170) and `INTER.csv` (39) is listed once, with what was
done with it and a reason. `tools/check_manifest.py` prints the tally:

```
check_manifest: OK (2 vendored file(s), 209 column(s), 12 used, 17 not-vendored entr(ies))
  139  not-ported
   21  not-yet
   10  not-a-value
   22  empty-upstream
    5  superseded-by
```

**`not-ported` is the porting backlog, and it is the number that matters.** NeqSim is
the target, not a reference: 139 of these columns are a physical property whose model
NeqSim implements and azoth has not ported, and each entry names the class that would
close it — `PhaseHydrate`, `CPAMixingRuleHandler`, `SolidFlash1`, `PhasePCSAFTa`,
`ParachorSurfaceTension` and the rest. The check refuses a `not-ported` reason with no
NeqSim name in it, so the list cannot drift back into being somewhere to put a column.

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

**Whether the vendored slice is current.** `manifest.yaml` records the NeqSim commit it
was last checked against (`dedba8735d030c6411e09b6fd7e69f6c4a136114`, v3.20.0), and
nothing here can tell you a newer NeqSim exists. The `not_vendored` list was built by
expanding NeqSim's resource directory by hand, at that commit.

The alternative - fetching from GitHub in CI - would make every build depend on a
network service the project does not control, to answer a question that only matters
when someone deliberately decides to re-vendor. So it is a human step, and the
staleness is visible in the file rather than hidden.
