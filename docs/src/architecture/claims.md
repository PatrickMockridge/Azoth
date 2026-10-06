# What a claim is worth

`docs/src/calculus/index.md` states what a *calculus* claim can be. This page states what a claim in
**this tree** can be, which is a larger set: the calculus is one of the things here that makes a
claim, and the databank, the card and the kernels are others.

**A claim is one of these, and each names what carries it. A claim that names nothing is `unheld`,
and has to say why** — which is the same rule `lean/guards.toml` makes for a guard whose domain
nobody can bound, and the same rule `docs/src/calculus/numerics.md` makes for a partial function.

| kind | what it is | carried by | held to it by |
|---|---|---|---|
| **proved** | a statement about the algebra: dimensions, identities, domain conditions, round-trips | a gated theorem in `lean/Azoth/` | `tools/check_lean_axioms.py`, `tools/audit_lean_claims.py` |
| **checked** | a statement about the implementation: this gate refuses, this checker fails | a runtime refusal or a CI checker | `tools/check_doc_claims.py`'s enforcement sweep |
| **data** | a **magnitude**: `Tc = 190.564 K`, a fitted `kij`, a Crane coefficient | a databank row with a `citation`, and a `verify_status` where one is recorded | `tools/provenance.py`, `tools/check_user_data.py` |
| **unheld** | a claim nothing carries | an `owner = "unformalised"` entry with a reason | `tools/check_guards.py` |
| **specified** | the layer does not exist yet, and the claim is what will make its tranche checkable | nothing yet, by definition | the tranche that builds the layer |
| **characterised** | the layer exists and nothing this repository can do would violate the claim | nothing, deliberately | nothing — proving it would add a theorem nobody consults |

**The last two are the calculus's own statuses and not kinds of their own.** A specified claim is one
that will be *proved* or *checked* when its layer lands, and a characterised one is a claim whose
falsification is impossible from inside this tree. Both are **unheld** in the strict sense, and the
`nothing` enforcement marker is the sweep's name for that — a page may claim it **only** where no
status says *proved*, which is the one thing `tools/check_doc_claims.py` refuses.

**And `data` is the kind no tree can promote.** `Tc = 190.564 K` is a measurement somebody made. It
is never *proved*, and `docs/src/calculus/values.md` says what *can* be done about it instead: a
domain declared, a bound enforced, and the value's truth left where it belongs. That is why the
databank's `citation` records where a value came from and never that it is right, and why
`tools/check_user_data.py` checks provenance rather than values.

**Two vocabularies, one idea.** `lean/guards.toml` names an `owner` for each edge case out of a
closed five — `lean`, `check`, `construction`, `port`, `unformalised` — and the calculus pages name a
`Status` out of three. They answer the same question from two sides: *what carries this?* The
mapping is the table above. A guard owned by `lean` is **proved**; by `check` or `construction`, it
is **checked**; by `port`, it is a fidelity claim and is **checked** by the oracle; and
`unformalised` is **unheld**, which is the only one of the five that has to say why.

The vocabulary is closed, and that is what `python/tests/test_claim_kinds.py` holds: every bold word
after a `*Status:` in the two normative trees, and every `owner` in `lean/guards.toml`, has to be one
the table above names.

*Status: **checked**, and the claim is narrower than it looks. The tree's claims are of the kinds above and of no others, and a seventh word fails a test rather than being read as emphasis.*

*Enforcement: check — `python/tests/test_claim_kinds.py` reads every status line in `docs/src/calculus/` and `docs/src/architecture/` and every `owner` in `lean/guards.toml`, and refuses a word this page does not name; `tools/check_doc_claims.py` holds each page's enforcement marker to a path that resolves, and `tools/check_guards.py` holds each `unformalised` owner to its reason.*
