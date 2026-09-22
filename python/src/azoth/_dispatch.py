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
from collections.abc import Callable, Iterator, Mapping
from contextlib import contextmanager
from types import ModuleType
from typing import Any, Literal, get_type_hints

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


def reference_for(calc_id: str) -> Callable[..., Any]:
    """The pure-Python implementation of `calc_id`, found by its id.

    A calculation's id **is** its address: `hydraulics.darcy_weisbach` names
    `azoth.hydraulics.reference.darcy_weisbach.darcy_weisbach`, with nothing in
    between. So this is a lookup rather than a table, and adding a calculation adds
    no entry anywhere.

    Deriving the whole path rather than only the function name is what keeps a calc in
    any namespace from resolving from `azoth.hydraulics.reference.<name>`, where it
    would fail with an `ImportError` *at call time* - the latest possible moment for a
    wiring mistake to surface.
    """
    namespace, _, function_name = calc_id.rpartition(".")
    # Annotated rather than inferred: `getattr` returns Any, which would make the
    # return an implicit Any and fail `--strict` at the boundary of exactly the
    # function whose job is to keep the two implementations interchangeable.
    resolved: Callable[..., Any] = getattr(
        importlib.import_module(f"azoth.{namespace}.reference.{function_name}"),
        function_name,
    )
    return resolved


#: The input a model names its one fluid by when it names only one.
MODEL_COMPONENTS_INPUT = "components"

#: The suffix a model names a *second* fluid by: `hot_components` gives the function
#: `hot_mixture` and `hot_ideal_gas`.
MODEL_COMPONENTS_SUFFIX = "_components"


def fluid_inputs(model: Mapping[str, Any]) -> list[tuple[str, str]]:
    """Every fluid a model names, as `(prefix, declared input name)`.

    A model declares a fluid with an input of ``type = "components"``, and the boundary
    resolves it into the objects the implementation actually takes: a `Mixture` and an
    `IdealGasModel`, built from the databank so a case and a caller reach the same fluid.
    Which argument names those land on is **derived from the input's own name** -
    `components` gives `mixture`/`ideal_gas`, `hot_components` gives
    `hot_mixture`/`hot_ideal_gas` - because a second convention beside the first would be
    a second thing to keep in step.

    **A model that names two fluids is the reason this exists.** `process.heat_exchanger`
    is the first: its two ports carry different fluids, so one `components` input could not
    describe it, and everything downstream - the case runner, the bridge, the stub
    generator - had assumed exactly one. The derivation is written for the general case so
    that a third would need no edit here.
    """
    out: list[tuple[str, str]] = []
    for name, declaration in model["inputs"].items():
        if declaration.get("type") != "components":
            continue
        prefix = (
            "" if name == MODEL_COMPONENTS_INPUT else name[: -len(MODEL_COMPONENTS_SUFFIX)] + "_"
        )
        out.append((prefix, name))
    return out


def result_type(calc_id: str) -> type[Any]:
    """The result dataclass a calculation produces, read from its own annotation.

    Derived rather than listed, and the annotation is the right place to read it
    from: `def friction_factor_colebrook(...) -> ColebrookResult` already says what
    the function returns, and a table beside it is a second answer to a question the
    code has answered. The table this replaces could not have been derived by naming
    convention either - the calc is `friction_factor_colebrook` and its result is
    `ColebrookResult` - which is exactly why it was hand-written, and exactly why it
    could drift.

    Raises:
        KeyError: if the function carries no return annotation. Every calculation
            does, and one that does not is a defect rather than a case to handle.
    """
    hints = get_type_hints(reference_for(calc_id))
    if "return" not in hints:
        raise KeyError(
            f"{calc_id!r} has no return annotation, so its result type cannot be "
            f"read. Every calculation in the registry annotates what it returns."
        )
    resolved: type[Any] = hints["return"]
    return resolved


def result_types() -> dict[str, type[Any]]:
    """Every calculation and model id, mapped to the result type it produces."""
    from azoth._models_gen import MODELS
    from azoth._registry_gen import CALCS

    return {entry["id"]: result_type(entry["id"]) for entry in [*CALCS, *MODELS]}


def resolve(calc_id: str) -> Callable[..., Any]:
    """The callable implementing `calc_id` on the selected backend.

    Always returns something that produces the Python result dataclass, so a
    caller never has to know which backend answered - the Rust binding's own
    result objects are adapted before they leave this module.
    """
    reference = reference_for(calc_id)

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
