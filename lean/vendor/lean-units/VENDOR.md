# lean-units, vendored

Taken from <https://github.com/ecyrbe/lean-units>, MIT, at commit
`5322fc77dd4e3ec8d31af2c757a47abbeb94d500`.

Vendored rather than required from a package index so that the definitions the
calculus is written against are the ones in this tree: `ecyrbe/lean-units` is a
young library, and a dimension system whose meaning can change under a build is
not one a proof should rest on.

Every file under `LeanUnits/` is unmodified, as are `LeanUnits.lean`,
`lakefile.toml`, `lean-toolchain` and `lake-manifest.json` - the manifest is what
pins Mathlib. `Examples/` and `.github/` are not taken, because nothing here
builds them. `LICENSE` travels with the copy.

`NOTICE` at the repository root records the attribution.

Its formal-proof side is a work in progress upstream, so this development treats
it as a *definitions* library: the dimension and unit structures, and the
`auto_dim`/`auto_equiv` tactics. Nothing here is proved by appealing to one of its
theorems, and `Azoth/Axioms.lean` enforces that automatically, because it reports
the transitive axiom set of every theorem this project claims.
