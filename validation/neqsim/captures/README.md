# Probe captures

A **capture** is a probe's standard output, committed verbatim. The probes live one
directory up, are compiled and run by hand against the pinned `neqsim-3.20.0.jar`, and what
they print is what lands here.

```bash
cd validation/neqsim
javac -proc:none -cp neqsim-3.20.0.jar CpaSweep.java
java -cp .:neqsim-3.20.0.jar CpaSweep > captures/cpa_sweep.tsv
```

They are committed rather than regenerated in CI because the jar is gitignored, so a gate
that ran the JVM could not run on a runner that has no NeqSim checkout. Committing the
output is what lets `tools/gen_neqsim_cases.py --check` be a build gate: it reads a capture
and emits `validation/eos/*.json`, so a case and the probe output it came from move
together or the gate fails.

**A capture is a pinned artifact, not a live one.** Regenerating one means running the probe
again, and the cases derived from it move with it — which is the point, and also the reason
they are here rather than fetched.

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
