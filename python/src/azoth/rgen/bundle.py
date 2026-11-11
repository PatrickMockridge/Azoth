"""Assemble a calculation run into the documents and payloads r-gen anchors.

The order of the hashes is the order of the evidence: the definition (what the calculation
is), the build (which source answered), the statement (this run's result), the cases (what
the calc's own evidence says), the agreement (the two implementations compared), the
quality (the four scores and their counts), and the validation evidence (the verdict and
the comment). Each document names its parent by hash, so an anchor whose parent is absent
is refused by the contract rather than read as a root.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Final

from azoth import _dispatch
from azoth._dispatch import Backend
from azoth._registry_gen import spec

from . import canonical, compare, inputs, payloads, quality, runner
from . import cases as case_module

__all__ = ["Bundle", "BundleRefused", "build"]

SCHEMA_BUNDLE: Final[str] = "azoth.rgen-bundle/v1"
SCHEMA_DEFINITION: Final[str] = "azoth.definition/v1"
SCHEMA_BUILD: Final[str] = "azoth.build/v1"
SCHEMA_STATEMENT: Final[str] = "azoth.statement/v1"
SCHEMA_CASES: Final[str] = "azoth.cases/v1"
SCHEMA_AGREEMENT: Final[str] = "azoth.agreement/v1"
SCHEMA_QUALITY: Final[str] = "azoth.quality/v1"
SCHEMA_VALIDATION: Final[str] = "azoth.validation-evidence/v1"


class BundleRefused(ValueError):
    """A bundle that cannot be built, with the reason named."""


@dataclass(frozen=True)
class Bundle:
    """The documents, their hashes, the payloads and the manifest."""

    documents: dict[str, dict[str, Any]]
    hashes: dict[str, str]
    claim: dict[str, Any]
    validation: dict[str, Any]
    anchors: list[dict[str, str]]
    manifest: dict[str, Any]

    def write(self, directory: Path) -> None:
        """Write the exact canonical bytes each hash names."""
        documents_dir = directory / "documents"
        payloads_dir = directory / "payloads"
        documents_dir.mkdir(parents=True, exist_ok=True)
        payloads_dir.mkdir(parents=True, exist_ok=True)

        for name, document in self.documents.items():
            data = canonical.dumps(document)
            if canonical.hash_bytes(data) != self.hashes[name]:
                raise BundleRefused(f"documents/{name}.json hashes differently from its record")
            (documents_dir / f"{name}.json").write_bytes(data)

        (payloads_dir / "claim.json").write_bytes(canonical.dumps(self.claim))
        (payloads_dir / "validation.json").write_bytes(canonical.dumps(self.validation))
        (payloads_dir / "anchors.json").write_bytes(canonical.dumps(self.anchors))
        (directory / "manifest.json").write_bytes(canonical.dumps(self.manifest))
        (directory / "inputs.json").write_bytes(
            canonical.dumps(self.documents["statement"]["inputs"])
        )


def _traceability(provenance: Mapping[str, Any], build: Mapping[str, Any]) -> int:
    items = 0
    if provenance.get("spec_sha256"):
        items += 1
    if provenance.get("python_sha256"):
        items += 1
    if provenance.get("rust_sha256"):
        items += 1
    cases = provenance.get("validation_cases", 0)
    if provenance.get("verification") in ("verified", "partially_verified") and int(cases) > 0:
        items += 1
    if build.get("git_commit") and build.get("tracked_dirty") is False:
        items += 1
    return items


def build(
    calc_id: str,
    kwargs: Mapping[str, Any],
    *,
    project: str,
    jurisdiction: str,
    claim_backend: Backend = "rust",
    validation_backend: Backend = "python",
    discipline: str = "process",
    ctype: str = "Validator",
    role: str = "validator",
    level: int = 1,
    deadline_height: int = 0,
    comment_until_height: int | None = None,
) -> Bundle:
    """Run a calculation on both implementations and assemble the bundle."""
    available = _dispatch.available()
    for backend in (claim_backend, validation_backend):
        if backend not in available:
            raise BundleRefused(
                f"the {backend!r} backend is not available; have {sorted(available)}"
            )
    if claim_backend == validation_backend:
        raise BundleRefused(
            "a validation needs a different implementation from the claim; "
            f"both are {claim_backend!r}"
        )

    declared = case_module.declared_cases(calc_id)
    if len(declared) == 0:
        raise BundleRefused(
            f"{calc_id!r} declares no value-carrying case, so a claim would have no validated case"
        )

    encoded = inputs.encode(calc_id, kwargs)
    run_kwargs = inputs.decode(calc_id, encoded)
    diagnostic_bound = runner.solver_tolerance(calc_id)
    outputs = spec(calc_id)["outputs"]

    claim_run = runner.run(calc_id, run_kwargs, claim_backend)
    validation_run = runner.run(calc_id, run_kwargs, validation_backend)
    matched = runner.match_case(calc_id, run_kwargs, declared)

    if matched is not None:
        tolerance = float(matched["tolerance"])
    else:
        strictest = case_module.strictest_tolerance(declared)
        if strictest is None:
            raise BundleRefused("no declared case states a tolerance")
        tolerance = strictest

    mirrored = compare.between_runs(
        claim_run["result"],
        validation_run["result"],
        tolerance=tolerance,
        diagnostic_bound=diagnostic_bound,
    )
    run_expected = None
    if matched is not None:
        run_expected = compare.against_expected(
            claim_run["result"],
            matched["expected"],
            outputs,
            tolerance=float(matched["tolerance"]),
            diagnostic_bound=diagnostic_bound,
        )
    run_disagree = bool(mirrored["disagreements"]) or (
        run_expected is not None and bool(run_expected["disagreements"])
    )

    reports = runner.replay(
        calc_id,
        declared,
        [claim_backend, validation_backend],
        diagnostic_bound=diagnostic_bound,
    )
    passed = sum(1 for report in reports if report["passed"])
    if passed == 0:
        raise BundleRefused(f"{calc_id!r}: no declared case passed, so the quantity is zero")

    build_document = {
        "schema": SCHEMA_BUILD,
        **runner.build_identity(),
    }
    build_hash = canonical.hash_value(build_document)

    provenance = claim_run["result"].get("provenance", {})
    evidence = {
        "declared": len(declared),
        "replayed": len(reports),
        "passed": passed,
        "warnings": len(claim_run["result"].get("warnings", [])),
        "run_matched": matched is not None,
        "run_disagree": run_disagree,
        "traceability_items": _traceability(provenance, build_document),
        "tracked_dirty": build_document.get("tracked_dirty"),
    }
    scored = quality.score(evidence)

    definition = {
        "schema": SCHEMA_DEFINITION,
        "calc_id": calc_id,
        "spec": spec(calc_id),
        "provenance": provenance,
    }
    statement = {
        "schema": SCHEMA_STATEMENT,
        "calc_id": calc_id,
        "inputs": encoded,
        "backend": claim_backend,
        "result": claim_run["result"],
        "build": build_hash,
    }
    cases_document = {
        "schema": SCHEMA_CASES,
        "calc_id": calc_id,
        "declared": declared.as_dict(),
        "replay": reports,
        "run_matched": None if matched is None else matched["id"],
        "tolerance": tolerance,
    }
    agreement = {
        "schema": SCHEMA_AGREEMENT,
        "calc_id": calc_id,
        "claimant_backend": claim_backend,
        "validator_backend": validation_backend,
        "validator_result": validation_run["result"],
        "mirrored": mirrored,
        "run_expected": run_expected,
        "run_matched": None if matched is None else matched["id"],
    }
    quality_document = {
        "schema": SCHEMA_QUALITY,
        "calc_id": calc_id,
        "quality": scored["quality"],
        "evidence": scored["evidence"],
        "verdict": scored["verdict"],
    }
    comment = payloads.comment_text(calc_id, str(scored["verdict"]), passed)
    validation_evidence = {
        "schema": SCHEMA_VALIDATION,
        "calc_id": calc_id,
        "verdict": scored["verdict"],
        "quality": scored["quality"],
        "comment": comment,
        "statement_hash": canonical.hash_value(statement),
        "cases_hash": canonical.hash_value(cases_document),
        "agreement_hash": canonical.hash_value(agreement),
        "quality_hash": canonical.hash_value(quality_document),
    }

    documents = {
        "definition": definition,
        "build": build_document,
        "statement": statement,
        "cases": cases_document,
        "agreement": agreement,
        "quality": quality_document,
        "validation-evidence": validation_evidence,
    }
    hashes = {name: canonical.hash_value(document) for name, document in documents.items()}

    claim = payloads.claim_request(
        calc_id=calc_id,
        definition_hash=hashes["definition"],
        build_hash=hashes["build"],
        statement_hash=hashes["statement"],
        quantity=passed,
        project=project,
        jurisdiction=jurisdiction,
        discipline=discipline,
        ctype=ctype,
        role=role,
        level=level,
        deadline_height=deadline_height,
        comment_until_height=comment_until_height,
    )
    validation = payloads.validation_request(
        verdict=str(scored["verdict"]),
        quality=scored["quality"],
        characteristics_hash=hashes["validation-evidence"],
        agreement_hash=hashes["agreement"],
        quality_hash=hashes["quality"],
        cases_hash=hashes["cases"],
        comment=comment,
    )
    anchor_list = payloads.anchors(hashes)

    manifest = {
        "schema": SCHEMA_BUNDLE,
        "calc_id": calc_id,
        "documents": {
            name: {
                "path": f"documents/{name}.json",
                "algorithm": canonical.HASH_ALGORITHM,
                "hash": digest,
            }
            for name, digest in hashes.items()
        },
        "payloads": {
            "claim": canonical.hash_value(claim),
            "validation": canonical.hash_value(validation),
            "anchors": canonical.hash_value(anchor_list),
        },
        "quantity": {"unit": "validated-cases", "value": passed},
        "quality": scored["quality"],
        "verdict": scored["verdict"],
        "environment": runner.environment(),
    }

    return Bundle(
        documents=documents,
        hashes=hashes,
        claim=claim,
        validation=validation,
        anchors=anchor_list,
        manifest=manifest,
    )
