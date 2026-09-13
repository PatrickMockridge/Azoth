# Verifying azoth

This library asks to be trusted with engineering numbers, so it should be
possible to check rather than believe. This document says what you can verify,
how to do it, and - just as importantly - what none of it proves.

The trust comes from git, signatures, and reproducible hashes, which are boring
and work. azoth operates no ledger, issues no token, and contains no smart
contract or consensus code.

One external service does appear: **source documents are pinned by Arweave
transaction ID**. That is a deliberate and narrow use - a tx ID is the hash of
its content, so pinning a standard there makes a citation checkable byte-for-byte
by anyone, forever, without trusting this repository. The property being used is
*content addressing plus an independent timestamp*, not decentralisation for its
own sake. See [How a source is referenced](#how-a-source-is-referenced).

## What is and is not established

| Question | Answered by | Limits |
|---|---|---|
| Did this code come from this repository, unmodified? | Signed tags, `provenance.json` | Only as strong as the key behind the tag |
| Was this wheel built from that commit, in CI? | cosign keyless signature | Proves the workflow identity, not a person |
| Is this number what the equation gives? | Worked examples, cross-language tests | Only for the cases tested - see each calc's page |
| Is this equation the right one for my situation? | **You.** The docs give source, range and assumptions | No signature helps here |

The last row is the one that matters. A verified artifact means the code is what
it claims to be. It does not mean the correlation is valid for your fluid, your
roughness, or your Reynolds number. Each calculation's page in the documentation
states the range it was validated over and the assumptions it makes, and those
are the things to check against your problem.

## Verify a release tag

Tags are signed with an SSH key. You need the maintainer's public key.

```bash
# The key that signs this repository's tags:
#   ssh-ed25519 AAAAC3NzaC1lZDI1NTE5AAAAIG/xPE6bujR5o37Fw2wU4r5zZZ8su3tQ3rbRoeAfZ/mL

git clone https://github.com/PatrickMockridge/Azoth
cd azoth

# Tell git which keys may sign for which identity. The principal must be the
# committer email, not the key comment.
echo 'patrickmockridge@gmail.com ssh-ed25519 AAAAC3NzaC1lZDI1NTE5AAAAIG/xPE6bujR5o37Fw2wU4r5zZZ8su3tQ3rbRoeAfZ/mL' \
  > /tmp/azoth-allowed-signers
git -c gpg.ssh.allowedSignersFile=/tmp/azoth-allowed-signers tag -v v0.1.0
```

A good result looks like:

```
Good "git" signature for patrickmockridge@gmail.com with ED25519 key SHA256:2RgnJwvzml7F6w384VZdS6aSPmlQYSwE4SVVaZqjRqY
```

That fingerprint is the thing to compare against a source you trust
independently - a conference talk, a signed email, the maintainer in person. A
signature you fetched from the same place as the tag proves much less.

**What this does not prove.** That the key was not stolen, and that the person
behind it is someone you should trust about fluid mechanics. It proves the tag
and the commit it points at came from the holder of that key.

## Verify a wheel

Release wheels are signed with [cosign](https://docs.sigstore.dev/) using
**keyless** signing. There is no long-lived key to leak: cosign uses the CI job's
short-lived OIDC identity, and the signature records which repository workflow
produced it.

```bash
# Install cosign, then download the wheel and its bundle from the release.
WHEEL=$(ls azoth-*.whl)
cosign verify-blob "$WHEEL" \
  --bundle "${WHEEL}.sigstore.json" \
  --certificate-identity-regexp '(?i)^https://github\.com/patrickmockridge/azoth/\.github/workflows/release\.yml@refs/tags/v.*$' \
  --certificate-oidc-issuer https://token.actions.githubusercontent.com
```

A good result reports `Verified OK` with the certificate identity and issuer
echoed back.

**What this does not prove.** It proves the artifact was produced by *that GitHub
Actions workflow*, in *that repository*, at *that tag*. It does not prove a human
reviewed it, and it depends on the Sigstore public-good infrastructure being
honest. It is a stronger claim than "downloaded from a URL", and a weaker one
than "I built it myself".

If you want the strongest available position, build from source at a verified
tag - the next section is how you check the result.

## Reproduce a calculation

Every calculation has a worked example whose inputs, expected outputs and
derivation are in its documentation page. Reproducing one needs no trust in this
project at all: the equation is printed, the numbers are printed, and the
arithmetic is checked by hand.

For the Darcy-Weisbach worked example:

```
dP = f * (L / D) * (rho * v**2 / 2)
   = 0.02 * (100 / 0.1) * (998 * 1.5**2 / 2)
   = 0.02 * 1000 * 1122.75
   = 22455.0 Pa
```

Then check the code agrees, in either language:

```bash
cargo run -p azoth-cli -- pipe --fluid water --flow 10 --diameter 0.1 --length 100

python -c "
import azoth
q = azoth.ureg.Quantity
r = azoth.hydraulics.darcy_weisbach(
    0.02, q(100.0,'m'), q(0.1,'m'), q(998.0,'kg/m**3'), q(1.5,'m/s'))
print(r.dp)
"
```

The two implementations are independent - one in Python, one in Rust - and the
test suite runs every spec case through both and compares. If they ever disagree,
that is a bug in one of them, and CI fails.

### How closely they agree depends on the arithmetic, and the claim is layered

An earlier version of this page said "they agree bit-for-bit on this example",
without qualification. That is true of `darcy_weisbach`, whose arithmetic is
`+ - * /` only - IEEE-754 requires those to be correctly rounded, so identical
source order gives identical bits. It is **not** true of everything here, and
saying it unqualified was the kind of claim this project exists to avoid.

| Claim | Where it holds | How it is checked |
|---|---|---|
| **Bit-identical** | Kernels using only `+ - * /` and `sqrt`. Covers `pr_kappa`, `pr_alpha_ab`, `pr_molar_volume`, `rachford_rice_binary`, `vdw1f_mix_binary`, and the flash's Rachford-Rice bisection | `to_bits()` equality, one CI runner |
| **Within a declared tolerance** | Anything using `ln`, `acos` or `cos` - `pr_departure`, the cubic's root locator - and any platform other than the CI runner | The tolerance the spec case declares |
| **The same algorithm** | Every solver and every model | **Identical iteration counts** on every case. The sharpest cheap check: a difference of one ulp that flips an iteration shows up as `13 vs 14` rather than as a mystifying `1e-11` drift |
| **The right answer** | The flash, saturation, bubble/dew | External `validation/` cases plus the identities that need no data |

Two consequences worth knowing as a caller. A count that agrees is strong
evidence the two ran the same loop; a count that differs means one of them
stopped at a different step even if the answers look close. And test cases for
iterated calculations have to be chosen *away* from a convergence boundary, since
a one-ulp difference there flips the count - if you add one, pick a state where
the residual at the stopping step is comfortably below the tolerance rather than
just under it.

**Read the warnings.** A result that used out-of-range inputs, or skipped a check
because an optional input was missing, says so. `r.warnings` is not decoration.

**Read `pt_flash`'s `phase`, not its `beta`.** The vapour fraction is `None` when
there is genuinely none to report - a feed with no two-phase solution, or an
iteration that converged to `x = y = z`. In the second case the model cannot tell
you *which* single phase the feed is: that needs a stability analysis it does not
perform. A `beta` outside `[0, 1]` but present is the negative flash, which is a
real reading rather than a failure. All of this is in the result's docstring and
in the model's spec page.

## Verify `provenance.json`

Each release attaches a `provenance.json` recording the git commit and tag, and a
SHA-256 for every spec, code file, test, data table and lock file that went into
it. Verifying it is a single command:

```bash
git checkout v0.1.0
python tools/provenance.py --verify provenance.json
```

A good result is `provenance: OK (N file(s) match provenance.json)`. Any file that
does not match is named individually:

```
  MISMATCH  hydraulics.darcy_weisbach.spec: specs/calcs/hydraulics/darcy_weisbach.yaml
            hashes to f35073abd2c29530..., recorded 0000000000000000...
```

Hashes are of the files' **committed bytes**, not of any canonicalised form. That
is deliberate: a canonicalisation scheme would be a second definition of the
content, and the two would eventually disagree.

**What this does not prove.** That the record itself is honest. `provenance.json`
is signed with the same cosign machinery as the wheels - verify both, and check
the record's `git.commit` matches the tag you checked in the first step.

### Reading the record

- `calcs[]` - one entry per calculation: the spec, the code in both languages,
  and the tests, each hashed. A file that does not exist is recorded as
  `"present": false` rather than omitted, so a missing test is visible.
- `shared[]` - files that affect **every** result: the range-check machinery, the
  solver, the fittings loader. Per-calc hashes are that calc's own files, not a
  transitive closure; read the two sections together.
- `data[]` - the fitting and fluid tables. These change results, so they are
  hashed.
- `git.dirty` - true means the record describes a working tree with uncommitted
  changes, which nobody else can reproduce. A release is never dirty.

## How a source is referenced

A data row may carry a `citation`: a string naming where its values came from, written
for a person to read. **Nothing requires one, and nothing validates one.**

There is deliberately no machine-readable `source_ref` or `source_locator` beside it.
An earlier version required both, with a rule that they agree — and that was a form to
fill in rather than a fact, because no tool can check whether a person read the document
a reference points at. A field that cannot be checked, required anyway, teaches people
to fill it in rather than to know the answer.

What the library owes instead is disclosure, and it is specific: `NOTICE` says what
ships and where it came from; the shipped data records that provenance in its own
citation column; and the engineer who supplied a keycard is accountable for their own
values and their right to use them. See
[Specification, S6](docs/src/spec.md#s6-provenance-is-the-engineers-job-not-the-librarys).

## What is deliberately not here

**No blockchain, no token, no smart contract.** Signatures and hashes answer
"is this what it claims to be" without a distributed ledger, and a ledger would
add a dependency without adding an answer.

**No long-lived signing keys in the repository.** The wheel signing is keyless by
design, so there is no key to leak and no secret to rotate.

## Future work, not present in v1

Two things would strengthen the trust story and are **not implemented**. They are
listed so nobody assumes more than exists.

**OpenTimestamps anchoring of releases.** A signed tag proves a commit existed,
but not *when* - a key holder can backdate a tag, and a compromised key can forge
history. An OpenTimestamps proof would anchor a commit hash into Bitcoin, giving
an independent lower bound on when it existed. It needs no new trust assumption
beyond Bitcoin's existence, and no token.

Note that this is a different question from the provenance of a *value*. Recording
which document a value was read from says where the number came from. Anchoring a
release commit would establish when *this code* existed. Neither is implemented for
values beyond the citation string, and the second is not implemented at all.

**A multi-validator registry.** Today, a result is validated by whoever wrote its
worked example. A registry where independent parties reproduce a calc's worked
example and sign their agreement would let a user see *how many* independent
checks a calculation has passed, rather than only that one exists. That is a
social mechanism with a data format attached, not a cryptographic one, which is
why "count the signatures" is the right design and a consensus protocol is not.

Neither is a substitute for the thing that matters most: a competent person
checking that the equation and its source are right for the situation at hand.
