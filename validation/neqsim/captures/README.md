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
