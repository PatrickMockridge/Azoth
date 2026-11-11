"""``python -m azoth.rgen``: emit a bundle, or verify one already written."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any, Literal, cast

from . import canonical
from .bundle import BundleRefused, build
from .inputs import InputError, UnknownCalc, UnsupportedCalc

Backend = Literal["python", "rust"]


def _backend(value: str) -> Backend:
    if value not in ("python", "rust"):
        raise argparse.ArgumentTypeError(f"{value!r} is not a backend; choose python or rust")
    return cast(Backend, value)


def _inputs(path: Path) -> dict[str, Any]:
    document = canonical.load(path)
    if not isinstance(document, dict):
        raise InputError(f"{path}: the inputs file must be a JSON object")
    return cast("dict[str, Any]", document)


def _bundle(args: argparse.Namespace) -> int:
    try:
        result = build(
            args.calc,
            _inputs(args.inputs),
            project=args.project,
            jurisdiction=args.jurisdiction,
            claim_backend=args.claim_backend,
            validation_backend=args.validation_backend,
            discipline=args.discipline,
            ctype=args.ctype,
            role=args.role,
            level=args.level,
            deadline_height=args.deadline_height,
            comment_until_height=args.comment_until_height,
        )
    except (
        InputError,
        BundleRefused,
        canonical.PayloadError,
        UnknownCalc,
        UnsupportedCalc,
    ) as error:
        print(f"azoth.rgen: {error}", file=sys.stderr)
        return 2
    result.write(args.out)
    if args.json:
        print(json.dumps(result.manifest, indent=2, sort_keys=True))
    else:
        quantity = result.manifest["quantity"]["value"]
        quality = result.manifest["quality"]
        print(
            f"azoth.rgen: {args.calc} -> {args.out} "
            f"({quantity} validated case(s), verdict {result.manifest['verdict']}, "
            f"quality {quality['completeness']}/{quality['correctness']}/"
            f"{quality['traceability']}/{quality['fitness']})"
        )
    return 0


def _verify(args: argparse.Namespace) -> int:
    directory = args.directory
    manifest = canonical.load(directory / "manifest.json")
    if not isinstance(manifest, dict):
        print("azoth.rgen: manifest.json is not an object", file=sys.stderr)
        return 1
    problems: list[str] = []

    documents = manifest.get("documents", {})
    if not isinstance(documents, dict):
        print("azoth.rgen: manifest documents are not an object", file=sys.stderr)
        return 1
    hashes: dict[str, str] = {}
    for name, record in documents.items():
        if not isinstance(record, dict):
            problems.append(f"{name}: record is not an object")
            continue
        path = directory / str(record.get("path"))
        if not path.is_file():
            problems.append(f"{name}: {path} is absent")
            continue
        actual = canonical.hash_bytes(path.read_bytes())
        if actual != record.get("hash"):
            problems.append(f"{name}: hashes {actual}, recorded {record.get('hash')}")
        hashes[str(name)] = str(record.get("hash"))

    for payload_name in ("claim", "validation", "anchors"):
        path = directory / "payloads" / f"{payload_name}.json"
        if not path.is_file():
            problems.append(f"payloads/{payload_name}.json is absent")
            continue
        actual = canonical.hash_value(canonical.load(path))
        recorded = manifest.get("payloads", {}).get(payload_name)
        if actual != recorded:
            problems.append(f"payloads/{payload_name}.json hashes {actual}, recorded {recorded}")

    anchors_path = directory / "payloads" / "anchors.json"
    if anchors_path.is_file():
        anchors = canonical.load(anchors_path)
        if isinstance(anchors, list):
            seen: set[str] = set()
            for anchor in anchors:
                if not isinstance(anchor, dict):
                    problems.append("an anchor is not an object")
                    continue
                parent = anchor.get("parentHash")
                if parent not in (None, "") and parent not in seen:
                    problems.append(
                        f"anchor {anchor.get('kind')}: parent {parent} is not yet anchored"
                    )
                seen.add(str(anchor.get("hash")))

    if problems:
        for problem in problems:
            print(f"azoth.rgen: {problem}", file=sys.stderr)
        return 1
    print(f"azoth.rgen: verified {len(hashes)} document(s) and 3 payload(s)")
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="python -m azoth.rgen")
    sub = parser.add_subparsers(dest="command", required=True)

    emitted = sub.add_parser("bundle", help="run a calculation and write its bundle")
    emitted.add_argument("--calc", required=True, help="a registered calculation id")
    emitted.add_argument("--inputs", type=Path, required=True, help="a JSON object of inputs")
    emitted.add_argument("--out", type=Path, required=True, help="the bundle directory")
    emitted.add_argument("--project", required=True)
    emitted.add_argument("--jurisdiction", required=True)
    emitted.add_argument("--claim-backend", type=_backend, default="rust")
    emitted.add_argument("--validation-backend", type=_backend, default="python")
    emitted.add_argument("--discipline", default="process")
    emitted.add_argument("--ctype", default="Validator")
    emitted.add_argument("--role", default="validator")
    emitted.add_argument("--level", type=int, default=1)
    emitted.add_argument("--deadline-height", type=int, default=0)
    emitted.add_argument("--comment-until-height", type=int, default=None)
    emitted.add_argument("--json", action="store_true")

    checked = sub.add_parser("verify", help="re-hash an existing bundle")
    checked.add_argument("--dir", dest="directory", type=Path, required=True)

    args = parser.parse_args(argv)
    if args.command == "bundle":
        return _bundle(args)
    return _verify(args)


if __name__ == "__main__":
    raise SystemExit(main())
