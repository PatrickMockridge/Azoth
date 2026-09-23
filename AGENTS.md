# Working with azoth

azoth is a chemical-engineering calculation library with an agentic layer on top.

- **Skills** (90) live in [`skills/`](skills/README.md): self-contained specialisms,
  each a `SKILL.md`. To load them into a tool, run `python tools/export_skills.py
  --target <dsh|claude|opencode>`; the default, `all`, writes all three.
- **Orchestration** is in [`agents/README.md`](agents/README.md). The one agent is a
  HAZOP team of four roles, defined in `agents/hazop/README.md` and exported to
  `.claude/agents/` and `.opencode/agents/`; run `/hazop` to start a study.

Setup:

```bash
pip install 'azoth[agent]'          # DeepSeek Harness SDK + runtime binary (team runtime)
python tools/export_skills.py       # write the skills into each tool's directory
```

See [`agents/README.md`](agents/README.md) for the full setup and caveats.

## The NeqSim oracle

azoth is a port of NeqSim, so NeqSim is the differential oracle: a disagreement with
it is the interesting failure. The 60 MB jar is **not vendored** (only NeqSim's
*data* files are, under `databank/sources/neqsim/`), and the drivers beside the
captures in `validation/neqsim/` are compiled by hand — see
[`validation/README.md`](validation/README.md). This is how to get the jar back on a
machine that has none; it needs a JDK and network, and nothing else. Measured on
2026-09-23 with JDK 21 — NeqSim's Maven wrapper fetches its own Maven.

```bash
git clone https://github.com/equinor/neqsim /tmp/neqsim-check/neqsim
cd /tmp/neqsim-check/neqsim
git checkout f0c7436c69                 # the commit the jar is named after
./mvnw -q -DskipTests package           # → target/neqsim-3.21.0.jar
cp target/neqsim-3.21.0.jar /path/to/azoth/validation/neqsim/neqsim-f0c7436.jar
```

**The name is the commit.** Every driver header says `-cp neqsim-f0c7436.jar`, which is
`neqsim-3.21.0.jar` built at `f0c7436c69` and renamed for the commit it came from — so
the captures belong to that revision, and moving to a newer NeqSim is a rename that
rewrites every driver header with it. The checkout lives in `/tmp` and is disposable:
rebuilding it is this recipe.

```bash
cd validation/neqsim
javac -proc:none -cp neqsim-f0c7436.jar HybridEosGeReactiveProbe.java
java -cp .:neqsim-f0c7436.jar HybridEosGeReactiveProbe > captures/hybrid_eos_ge_reactive_probe.tsv
```

The captures are committed and CI reads them, so a new one is a capture file plus the
numbers in the case beside it, never a live JVM in the build.
