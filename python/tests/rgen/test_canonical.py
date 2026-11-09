"""The bundle's bytes and hashes are stable and refuse what JSON cannot carry."""

from __future__ import annotations

from pathlib import Path

import pytest

from azoth.rgen import canonical


def test_blake2b_vector_matches_r_gen() -> None:
    assert canonical.hash_bytes(b"abc") == (
        "bddd813c634239723171ef3fee98579b94964e3bb1cb3e427262c8c068d52319"
    )


def test_key_order_does_not_change_the_bytes() -> None:
    left = canonical.dumps({"b": 2, "a": 1})
    right = canonical.dumps({"a": 1, "b": 2})
    assert left == right == b'{"a":1,"b":2}'


def test_a_float_is_written_at_its_shortest_round_trip() -> None:
    assert canonical.dumps({"v": 0.1}) == b'{"v":0.1}'


def test_nan_is_refused() -> None:
    with pytest.raises(canonical.PayloadError):
        canonical.dumps({"v": float("nan")})


def test_loading_refuses_the_non_finite_tokens(tmp_path: Path) -> None:
    path = tmp_path / "payload.json"
    path.write_text('{"v": NaN}', encoding="utf-8")
    with pytest.raises(canonical.PayloadError):
        canonical.load(path)


def test_hash_value_is_the_hash_of_its_bytes() -> None:
    value = {"a": [1, "x", None]}
    assert canonical.hash_value(value) == canonical.hash_bytes(canonical.dumps(value))
