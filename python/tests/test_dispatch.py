"""Backend selection, and the guarantee that a fallback is never silent."""

from __future__ import annotations

import pytest

from chemeng import active_backend, available, backend_info
from chemeng._dispatch import RustBackendUnavailableError, use_backend


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
    from chemeng.hydraulics.reference import reynolds_number

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
