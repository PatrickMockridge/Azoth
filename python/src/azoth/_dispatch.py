"""Choosing between the pure-Python and Rust implementations.

Both implementations are always reachable, and the choice is never silent: the
backend in use can be queried, and ``AZOTH_REQUIRE_RUST=1`` turns "the
extension is missing" from a quiet fallback into a hard failure.

That last setting matters more than it looks. Without it, a broken wheel build
makes every test pass on the Python path while the cross-language agreement the
project promises is exercised by nothing - the suite stays green and the
guarantee evaporates. CI runs the cross-implementation job with it set.
"""

from __future__ import annotations

import importlib
import os
from collections.abc import Callable, Iterator
from contextlib import contextmanager
from types import ModuleType
from typing import Any, Literal

Backend = Literal["python", "rust"]

#: Env var that forces a particular backend. ``auto`` (the default) prefers Rust.
_ENV_BACKEND = "AZOTH_BACKEND"

#: Env var that makes a missing extension fatal.
_ENV_REQUIRE_RUST = "AZOTH_REQUIRE_RUST"


class RustBackendUnavailableError(RuntimeError):
    """The Rust extension was requested but cannot be used."""


def _load_extension() -> ModuleType | None:
    """Import ``azoth._core`` if it was built.

    Only ``ImportError`` is swallowed. A module that exists but fails to import
    for some other reason - a missing shared library, an ABI mismatch - is a real
    problem and must not be mistaken for "not built".
    """
    try:
        return importlib.import_module("azoth._core")
    except ImportError:
        return None


_EXTENSION: ModuleType | None = _load_extension()

#: Set by `use_backend` to override selection within a block.
_override: Backend | None = None


def available() -> frozenset[Backend]:
    """Which backends can be used right now."""
    if _EXTENSION is None:
        return frozenset({"python"})
    return frozenset({"python", "rust"})


def select() -> Backend:
    """Which backend calculations will use."""
    if _override is not None:
        return _override

    env = os.environ.get(_ENV_BACKEND, "auto").strip().lower()
    if env in ("python", "rust"):
        if env == "rust" and _EXTENSION is None:
            raise RustBackendUnavailableError(
                f"{_ENV_BACKEND}=rust but azoth._core is not built. "
                f"Build it with `maturin develop`."
            )
        return env  # type: ignore[return-value]

    if os.environ.get(_ENV_REQUIRE_RUST) == "1" and _EXTENSION is None:
        raise RustBackendUnavailableError(
            f"{_ENV_REQUIRE_RUST}=1 but azoth._core is not built. "
            f"Cross-implementation tests must not silently fall back to Python - "
            f"that is the whole reason this setting exists."
        )

    return "rust" if _EXTENSION is not None else "python"


@contextmanager
def use_backend(name: Backend) -> Iterator[None]:
    """Force a backend for the duration of a block.

    Used by the cross-implementation tests, which need to call both sides in the
    same process.

    Raises:
        RustBackendUnavailableError: if `name` is not available.
    """
    global _override
    if name not in available():
        raise RustBackendUnavailableError(
            f"backend {name!r} is not available; have {sorted(available())}"
        )
    previous, _override = _override, name
    try:
        yield
    finally:
        _override = previous


def extension() -> ModuleType | None:
    """The compiled extension module, or ``None`` if it was not built."""
    return _EXTENSION


def resolve(calc_id: str) -> Callable[..., Any]:
    """The callable implementing `calc_id` on the selected backend.

    Always returns something that produces the Python result dataclass, so a
    caller never has to know which backend answered - the Rust binding's own
    result objects are adapted before they leave this module.
    """
    # Both segments come from the id. The namespace used to be hardcoded as
    # `hydraulics` while only the function name was derived, so a calc in any other
    # namespace resolved its reference implementation from `azoth.hydraulics.
    # reference.<name>` and failed with an ImportError *at call time* - not at
    # build time, not at import time, but the first time a caller used it. That is
    # the latest possible moment for a wiring mistake to surface, which is what made
    # it worth deriving rather than listing.
    namespace, _, function_name = calc_id.rpartition(".")
    # Annotated rather than inferred: `getattr` returns Any, which would make the
    # return below an implicit Any and fail `--strict` at the boundary of exactly
    # the function whose job is to keep the two implementations interchangeable.
    reference: Callable[..., Any] = getattr(
        importlib.import_module(f"azoth.{namespace}.reference.{function_name}"),
        function_name,
    )

    if select() == "python":
        return reference

    try:
        bridge = importlib.import_module("azoth._rust_bridge")
    except ImportError as exc:
        raise RustBackendUnavailableError(
            "azoth._core is built but azoth._rust_bridge is not, so Rust "
            "results cannot be adapted to the Python result types. This is a "
            "packaging bug, not a missing feature."
        ) from exc
    resolved: Callable[..., Any] = bridge.resolve(calc_id)
    return resolved


def describe() -> dict[str, Any]:
    """What this installation can do, for diagnostics and provenance."""
    return {
        "available": sorted(available()),
        "selected": select(),
        "extension_loaded": _EXTENSION is not None,
    }
