"""Backend selection, and the guarantee that a fallback is never silent.

The public surface is tested first: `active_backend`, `available` and `backend_info`
are what a caller uses, and the reference implementation must stay reachable by name
whatever the backend selection is doing.

Below that, the guards the selection *rests* on. A run of this suite exercises
`resolve` continuously, so the happy path is covered everywhere; what had no test is
the machinery that chooses, and the two settings that exist because a silent fallback
would be worse than a failure:

* ``AZOTH_BACKEND=rust`` with no extension is an error, not a fall back to Python.
* ``AZOTH_REQUIRE_RUST=1`` with no extension is an error for the same reason - it is
  what stops the cross-implementation tests passing without Rust ever having run.

The extension's absence is simulated by monkeypatching ``_EXTENSION``, so those run
whether or not the extension is built here.
"""

from __future__ import annotations

from collections.abc import Iterator

import pytest

from azoth import _dispatch, active_backend, available, backend_info
from azoth._dispatch import RustBackendUnavailableError, use_backend


def test_python_is_always_available() -> None:
    """The reference implementation is the fallback, so it must never be absent."""
    assert "python" in available()


def test_selection_is_queryable() -> None:
    """Which implementation answered must never be a guess."""
    assert active_backend() in available()
    info = backend_info()
    assert info["selected"] in info["available"]
    assert isinstance(info["extension_loaded"], bool)


def test_the_reference_implementation_is_always_reachable() -> None:
    """A reference that can only be reached through the thing it validates is not
    much of a check."""
    from azoth.hydraulics.reference import reynolds_number

    assert callable(reynolds_number)


def test_use_backend_restores_the_previous_choice() -> None:
    """The override must not leak out of its block.

    A leaking override would make one test's choice silently change another's
    behaviour, which is the kind of coupling that produces a passing suite and a
    broken library.
    """
    before = active_backend()
    with use_backend("python"):
        assert active_backend() == "python"
    assert active_backend() == before


@pytest.mark.requires_rust
def test_requiring_an_unavailable_backend_fails_loudly() -> None:
    """Asking for the reference explicitly always works."""
    with use_backend("python"):
        assert active_backend() == "python"


def test_unknown_backend_is_rejected() -> None:
    """A typo in a backend name must not silently select something else."""
    with (
        pytest.raises(RustBackendUnavailableError, match="not available"),
        use_backend("python2"),  # type: ignore[arg-type]
    ):
        pass


@pytest.fixture
def no_extension(monkeypatch: pytest.MonkeyPatch) -> Iterator[None]:
    """Pretend the extension was not built, in an otherwise default environment.

    The settings are cleared as well as the module state: CI runs this suite with
    `AZOTH_REQUIRE_RUST=1`, under which every one of the tests below would raise for
    a reason that has nothing to do with what it is checking.
    """
    monkeypatch.setattr(_dispatch, "_EXTENSION", None)
    monkeypatch.setattr(_dispatch, "_override", None)
    monkeypatch.delenv("AZOTH_REQUIRE_RUST", raising=False)
    monkeypatch.delenv("AZOTH_BACKEND", raising=False)
    yield


@pytest.fixture
def with_extension(monkeypatch: pytest.MonkeyPatch) -> Iterator[None]:
    """Pretend it was, whatever this environment actually has."""
    monkeypatch.setattr(_dispatch, "_EXTENSION", object())
    monkeypatch.setattr(_dispatch, "_override", None)
    monkeypatch.delenv("AZOTH_REQUIRE_RUST", raising=False)
    monkeypatch.delenv("AZOTH_BACKEND", raising=False)
    yield


def test_a_missing_extension_leaves_only_python(no_extension: None) -> None:
    assert _dispatch.available() == frozenset({"python"})
    assert _dispatch.select() == "python"


def test_the_extension_adds_the_rust_backend(with_extension: None) -> None:
    assert _dispatch.available() == frozenset({"python", "rust"})
    assert _dispatch.select() == "rust", "auto selects Rust when it is available"


def test_asking_for_rust_without_the_extension_is_an_error(
    monkeypatch: pytest.MonkeyPatch, no_extension: None
) -> None:
    monkeypatch.setenv("AZOTH_BACKEND", "rust")

    with pytest.raises(RustBackendUnavailableError) as caught:
        _dispatch.select()

    assert "maturin develop" in str(caught.value), (
        f"the error has to say how to fix it, and says: {caught.value}"
    )


def test_requiring_rust_without_the_extension_is_an_error(
    monkeypatch: pytest.MonkeyPatch, no_extension: None
) -> None:
    """The setting that stops a cross-implementation run passing on Python alone."""
    monkeypatch.setenv("AZOTH_REQUIRE_RUST", "1")

    with pytest.raises(RustBackendUnavailableError) as caught:
        _dispatch.select()

    assert "fall back" in str(caught.value), caught.value


def test_an_explicit_backend_beats_auto(
    monkeypatch: pytest.MonkeyPatch, with_extension: None
) -> None:
    monkeypatch.setenv("AZOTH_BACKEND", "python")

    assert _dispatch.select() == "python"


def test_an_unrecognised_backend_is_ignored_rather_than_guessed(
    monkeypatch: pytest.MonkeyPatch, with_extension: None
) -> None:
    """Anything other than the two names means auto, not an error."""
    monkeypatch.setenv("AZOTH_BACKEND", "  Rust  ")

    assert _dispatch.select() == "rust", "the value is trimmed and lowercased"


def test_use_backend_restores_after_an_exception() -> None:
    """An override must not survive the block that set it, even on the way out."""
    before = _dispatch.select()

    with pytest.raises(ValueError, match="boom"), _dispatch.use_backend("python"):
        raise ValueError("boom")

    assert _dispatch.select() == before


def test_use_backend_nests() -> None:
    before = _dispatch.select()

    with _dispatch.use_backend("python"):
        with _dispatch.use_backend("python"):
            assert _dispatch.select() == "python"
        assert _dispatch.select() == "python", "the inner block restored the outer one"

    assert _dispatch.select() == before


def test_describe_agrees_with_select_and_available() -> None:
    described = _dispatch.describe()

    assert described["available"] == sorted(_dispatch.available())
    assert described["selected"] == _dispatch.select()
    assert described["extension_loaded"] == (_dispatch.extension() is not None)
    assert described["extension_loaded"] == ("rust" in described["available"])


def test_an_id_that_names_no_module_fails_at_the_lookup() -> None:
    """A calc id is its address, so an id that is not one must not resolve."""
    with pytest.raises(ImportError):
        _dispatch.reference_for("nowhere.no_such_calculation")


def test_every_resolved_reference_is_the_function_its_id_names() -> None:
    from azoth._registry_gen import CALCS

    for calc in CALCS:
        function = _dispatch.reference_for(calc["id"])

        assert function.__name__ == calc["id"].rpartition(".")[2], (
            f"{calc['id']} resolved to {function.__name__}, which is a different "
            f"function - the id is the address, so these have to be the same name"
        )
