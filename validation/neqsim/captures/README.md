# Probe captures

A **capture** is a probe's standard output, committed verbatim. The probes live one
directory up, are compiled and run by hand against the pinned `neqsim-f0c7436.jar`, and what
they print is what lands here.

```bash
cd validation/neqsim
javac -proc:none -cp neqsim-f0c7436.jar CpaSweep.java
java -cp .:neqsim-f0c7436.jar CpaSweep > captures/cpa_sweep.tsv
```

**A probe that declares NeqSim's own package compiles with `-d .` and runs under that
package's name**, because what it needs is package-private; `WaxReferenceProbe` and
`PcsaftCompositionProbe` are the two, and each says in its header why.

```bash
javac -proc:none -cp neqsim-f0c7436.jar -d . PcsaftCompositionProbe.java
java -cp .:neqsim-f0c7436.jar neqsim.thermo.component.PcsaftCompositionProbe \
  > captures/pcsaft_composition_probe.tsv
```

**A state is passed as `T_K P_bara name z name z ...`** - the name and its mole fraction are
two arguments, not one `name:z` token. Seven probes take that form, and each refuses an odd
argument count rather than silently dropping a name whose fraction is missing. The other
three that read arguments take `T_K P_bara n_water` triplets.

**A state that a probe only prints when it is invoked for one needs its own recipe here**,
because this file is the registry `tools/oracle_sweep.py` reads: a run with no recipe has
no committed capture to compare against, so the next pin move cannot say whether it moved.
`PcsaftProbe`'s three aside from its two defaults are the ones a test asserts.

```bash
java -cp .:neqsim-f0c7436.jar PcsaftProbe 150 50 methane 1 > captures/pcsaft_probe_methane_150_50.tsv
java -cp .:neqsim-f0c7436.jar PcsaftProbe 300 100 propane 1 > captures/pcsaft_probe_propane_300_100.tsv
java -cp .:neqsim-f0c7436.jar PcsaftProbe 400 50 methane 1 > captures/pcsaft_probe_methane_400_50.tsv
```

**`ProcessProbe` takes a unit operation's name and prints that unit operation's inlet and
outlet on the palette's record**, so one recipe per unit operation and one capture each. It
is the only probe here that drives `neqsim.process` rather than `neqsim.thermo`.

```bash
java -cp .:neqsim-f0c7436.jar ProcessProbe pump > captures/process_pump.tsv
java -cp .:neqsim-f0c7436.jar ProcessProbe splitter > captures/process_splitter.tsv
java -cp .:neqsim-f0c7436.jar ProcessProbe mixer > captures/process_mixer.tsv
java -cp .:neqsim-f0c7436.jar ProcessProbe separator > captures/process_separator.tsv
java -cp .:neqsim-f0c7436.jar ProcessProbe throttling_valve > captures/process_throttling_valve.tsv
java -cp .:neqsim-f0c7436.jar ProcessProbe heat_exchanger > captures/process_heat_exchanger.tsv
java -cp .:neqsim-f0c7436.jar ProcessProbe stream > captures/process_stream_properties.tsv
```

They are committed rather than regenerated in CI because the jar is gitignored, so a gate
that ran the JVM could not run on a runner that has no NeqSim checkout. Committing the
output is what lets `tools/gen_neqsim_cases.py --check` be a build gate: it reads a capture
and emits `validation/eos/*.json`, so a case and the probe output it came from move
together or the gate fails.

**`tools/oracle_sweep.py` runs them all at once and is the step a pin move starts with.**
It re-runs every probe and every recipe below, compares each against the capture committed
here, and reports the keys that moved - and then any number in a Rust test that still holds
the *previous* value, which is what a commit regenerating these cannot see on its own. It
also names the probes with no capture and the captures no run reproduces; both are gaps of
the same kind, one on each side.

**A capture is a pinned artifact, not a live one.** Regenerating one means running the probe
again, and the cases derived from it move with it — which is the point, and also the reason
they are here rather than fetched.

**One capture cannot be regenerated at all.** `multi_salt_probe.tsv` is the multi-mineral
complementarity loop's only oracle — a brine of `Ca++`, `SO4--`, `CO3--` and `HCO3-` whose
three minerals take three rounds — and it was written before `PitzerParameterCoverage`
existed. Both pinned revisions refuse that brine, because the legacy dataset carries no theta
for `Ca++|Na+` or `Cl-|SO4--`, so the capture is a *record* of a state no current NeqSim will
run rather than something to re-measure. `tools/oracle_sweep.py` reports it as a capture with
no probe, which is the honest form of that.

Each capture's own header says what it is and, where the probe emits a key table, what each
key came from.

## Which captures a layer diff can read

`tools/neqsim_layer_diff.py` compares a model's intermediates against a capture key by key,
and it can only read `key = value` rows. **Two of the electrolyte captures are not in that
shape at all**: `pitzer_arithmetic.tsv` and `soreide_alpha_derivatives.tsv` print aligned
columns of numbers with no keys - `T`, `alpha`, `dalpha/dT`, and the Debye-Hückel parameter
against temperature - so they are read by hand, and the models built against them are pinned
by their cases rather than by a layer diff.

`pitzer_probe.tsv`, `soreide_whitson_probe.tsv` and `ge_electrolyte_probe.tsv` *are* keyable,
and still have none, because what they print is **which parameters a model selected** - the
dataset id, the active ions, the missing binary terms, the `kij` table - rather than a
model's own arithmetic. A layer diff exists to say which layer of a calculation moved, and
that is not a question those captures answer.

`furst_probe.tsv` was written in the readable shape and in the order the phase builds its
layers, for exactly this reason, and `python/src/azoth/eos/layers.py` carries the dumper
that meets it.
