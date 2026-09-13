"""Shared test configuration.

The one job here is to make "the extension is missing" behave differently
depending on why it is missing.

Ordinarily a test needing the Rust extension skips when it has not been built -
that is correct for a developer working only on Python. But in CI, where the
cross-language agreement is a promise the project makes, a silent skip means the
promise is checked by nothing and the suite stays green anyway. Setting
``AZOTH_REQUIRE_RUST=1`` turns those skips into failures.
"""

from __future__ import annotations

import os

import pytest

from azoth import available


def pytest_collection_modifyitems(items: list[pytest.Item]) -> None:
    """Skip Rust-only tests, or fail them when Rust is required."""
    if "rust" in available():
        return
    require_rust = os.environ.get("AZOTH_REQUIRE_RUST") == "1"
    for item in items:
        if "requires_rust" not in item.keywords:
            continue
        if require_rust:
            item.add_marker(
                pytest.mark.xfail(
                    reason=(
                        "AZOTH_REQUIRE_RUST=1 but azoth._core is not built. The "
                        "cross-implementation guarantee must not quietly degrade into "
                        "a skip - build the extension with `maturin develop`, or unset "
                        "the variable if you are working on Python only."
                    ),
                    strict=True,
                    run=True,
                )
            )
        else:
            item.add_marker(
                pytest.mark.skip(reason="azoth._core is not built; run `maturin develop`")
            )
