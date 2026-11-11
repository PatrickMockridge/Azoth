"""The bytes an r-gen anchor names, and the hash of those bytes.

A document is written once, hashed once, and the file keeps exactly the bytes that were
hashed. Re-serialising a loaded document would be a second writer, and the two would
eventually disagree about key order or a float's spelling.

The hash is BLAKE2b-256 in lowercase hex. That is the digest r-gen's own
``harness/demo_seed.py`` computes for the documents it anchors, so an azoth bundle and a
hand-seeded r-gen document are hashed the same way. Azoth's ``provenance.json`` stays
SHA-256; the two algorithms are named wherever a hash is written and are never compared
with each other.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any, Final

__all__ = ["HASH_ALGORITHM", "PayloadError", "dumps", "hash_bytes", "hash_value", "load"]

#: The algorithm name written beside every hash the bundle produces.
HASH_ALGORITHM: Final[str] = "blake2b-256"

_DIGEST_SIZE: Final[int] = 32


class PayloadError(ValueError):
    """A value that cannot be written to a bundle or sent to r-gen."""


def dumps(value: Any) -> bytes:
    """Canonical JSON bytes.

    Sorted keys, no insignificant whitespace, UTF-8, and no NaN: a non-finite magnitude is
    written as ``null`` by the result walk before it reaches here, and a raw NaN reaching
    this function is a defect rather than a value to encode.
    """
    try:
        text = json.dumps(
            value,
            sort_keys=True,
            separators=(",", ":"),
            ensure_ascii=False,
            allow_nan=False,
        )
    except (TypeError, ValueError) as error:
        raise PayloadError(f"not canonical JSON: {error}") from error
    return text.encode("utf-8")


def hash_bytes(data: bytes) -> str:
    """BLAKE2b-256 of bytes already written."""
    return hashlib.blake2b(data, digest_size=_DIGEST_SIZE).hexdigest()


def hash_value(value: Any) -> str:
    """BLAKE2b-256 of a value's canonical bytes."""
    return hash_bytes(dumps(value))


def _reject_constant(token: str) -> Any:
    raise PayloadError(f"{token} is not JSON; a non-finite magnitude has no document form")


def load(path: Path) -> Any:
    """Read a JSON file, refusing the non-finite tokens ``json`` accepts by default."""
    return json.loads(path.read_text(encoding="utf-8"), parse_constant=_reject_constant)
