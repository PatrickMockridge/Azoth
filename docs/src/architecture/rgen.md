# The r-gen bridge

r-gen models engineering work as a claim, a validation and the K-Asset the pair becomes.
Azoth computes the work and already carries its evidence: a spec per calculation, two
implementations, and a provenance block on every result naming the spec, both kernels and
the external cases.

`azoth.rgen` turns a calculation run into that evidence in r-gen's shape. It is a
**producer of documents**, not an r-gen client: it has no dependency on r-gen, computes no
physics of its own beyond running the registered calculation on both implementations, and
writes a directory another program signs and deploys.

## The mapping

| azoth | r-gen |
|---|---|
| the spec and both implementations' hashes | the claim's `definitionHash` |
| the run's inputs, backend and result document | the claim's `statementHash` |
| the external validation cases and the spec's own | `documents/cases.json` |
| the two implementations' agreement | `documents/agreement.json` |
| the four scores and the counts they came from | `documents/quality.json` |
| the verdict and the comment | the validation payload |
| the ordered documents | the asset's provenance anchors |

The quantity is an integer count of validated cases, because r-gen's wire has no floats.
The computed values live in the anchored result document; the chain carries how many cases
stood behind them, not the numbers themselves.

## Documents and hashes

A bundle is a directory with `payloads/`, `documents/`, `manifest.json` and `inputs.json`.
Seven documents form a DAG, each naming its parent by hash:

```
definition -> build -> statement -> cases -> agreement -> quality -> validation-evidence
```

The bytes of each file under `documents/` are exactly the bytes its hash names, so a
verifier re-hashes the file rather than re-serialising it. The hash is BLAKE2b-256, the
digest r-gen's own demo seeding uses; azoth's `provenance.json` stays SHA-256, and the two
are named rather than compared.

`build.json` is the stable source identity: the git commit, a dirty flag read from tracked
files only, the package and core versions and `provenance.json`'s hash. The platform and
extension record is in `manifest.json` and is outside every hashed document, because a
compiled extension is platform-specific and the claim is about this build's arithmetic.

## The four scores

The validation carries exactly the four keys `kasset.conventions/v1` names. The derivation
and its counts are in `documents/quality.json`, so a consumer can recompute the judgement.

- **correctness** is measured against each case's expected values - the external oracle
  cases where the calc has them. Two implementations agreeing is not correctness: they are
  the same spec written twice. It is completeness.
- **completeness** is the fraction of the declared case set that was replayed.
- **traceability** counts five provenance items: the spec hash, the Python kernel hash, the
  Rust kernel hash, a `verified` or `partially_verified` status with at least one external
  case, and a clean build identity. A dirty tracked tree caps it at 3.
- **fitness** falls with warnings: none is 5, one is 4, two or more is 3.

A run whose inputs match a declared case is validated against that case. A run that matches
none is exploratory: it is recorded, it is not counted, and correctness cannot rise above 3.
A case set with no passing case is refused, because a claim with nothing validated has no
quantity the conventions allow.

## The split

`python -m azoth.rgen bundle` writes the bundle. The r-gen side checks it and signs the
deploys: a claim by the claimant's key, a validation by a different validator key and
credential, the K-Asset, and the seven anchors. The two identities are deployer keys, which
the contract can check; the validator's run is the mirrored implementation, which is an
implementation-level check rather than an organisational one. The independent evidence is
the external oracle case set, and a production validator runs the bundle on its own
checkout.

The bridge has no Rust twin because it computes no physics. Its correctness is the
determinism of its bytes, the shared comparison policy and the independent verifier on the
r-gen side.

## Refusals, named

- A calculation with a **`record_list` input** is refused with the calc id and the field.
  `thermal.zone_rating` and `thermal.incremental_zone_analysis` are the two today.
- **Models** are not calculations in this sense; a model bundle is a separate unit.
- A **validation needs two available implementations**; a claim may be built with one, but a
  validation with one would not be a validation.
- A **declared case set with no passing case** is refused rather than given a zero quantity.
- An installed wheel carries no external `validation/` tree, so its case set is the spec's
  own; the provenance block still records the derived status and the external case count.

## Command

```sh
python -m azoth.rgen bundle \
  --calc eos.antoine_vapor_pressure \
  --inputs inputs.json \
  --out artifacts/antoine \
  --project azoth-demo --jurisdiction US-TX \
  --claim-backend rust --validation-backend python

python -m azoth.rgen verify --dir artifacts/antoine
```

`deadlineHeight` and `credentialUri` are placeholders in the written payloads. The contract
refuses a zero deadline and an empty credential, so a bundle deployed without substituting
them fails loudly rather than filing a claim with no clock.

*Enforcement: check — `python -m azoth.rgen verify` re-hashes every document and payload, the r-gen harness re-hashes them independently, and `python/tests/rgen/` holds the document graph, the integer-only payloads and the shared comparison policy to their tests.*
