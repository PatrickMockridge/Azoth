"""What the Rust backend's adapters need that is not an adapter.

`python/src/azoth/_rust_bridge.py` holds one function per registered id and a handful of helpers
they all call. The functions are generated now - `tools/gen_python_bridge.py` emits
`azoth/_rust_bridge_gen.py` - and a generator cannot emit a helper that is itself logic, so the
helpers live here and both modules import them.

**The move is what makes the generated module importable.** A generated module that imported
`_warnings` from `_rust_bridge` would be importing the file that re-exports it, which is the cycle
this file exists to break.
"""

from __future__ import annotations

from collections.abc import Sequence
from typing import Any

from azoth import _core
from azoth.core.units import Q, input_to_si
from azoth.core.warnings import Warning, WarningCode

__all__ = ["_association_spec", "_si", "_warnings"]


def _warnings(raw: Sequence[_core.Warning]) -> tuple[Warning, ...]:
    """Convert transported warnings to the real ones.

    Codes become the enum rather than staying strings, so a caller can compare
    ``warning.code == WarningCode.TRANSITIONAL_FLOW`` regardless of which backend
    produced it. A code the Python side does not know is impossible - the two
    sets are asserted equal by a test - but the ValueError would surface here if
    it ever happened, which is better than a string quietly failing to match.
    """
    return tuple(Warning(WarningCode(w.code), w.message, w.field) for w in raw)


def _si(spec: dict[str, object], name: str, value: float | Q) -> float:
    """A declared input as its SI magnitude, quantity or not.

    **The boundary is mixed by design.** A vector whose spec declares a unit arrives as a
    pint quantity and has to be converted; a dimensionless one - ``chem_ref``,
    ``log_activity``, the element matrix - arrives as the bare number it is, and
    ``input_to_si`` refuses that rather than passing it through.
    """
    # A quantity or a bare number. `Q` is a type alias rather than a class, so the check
    # is made the other way round: the numeric branch is the one `isinstance` can name.
    if isinstance(value, int | float):
        return float(value)
    return input_to_si(spec, name, value)


def _association_spec(mixture: Any) -> Any:
    """The mixture's association, in the form the Rust boundary takes.

    **Every model whose Python side takes a ``Mixture`` crosses this**, and none of them
    may default it away. A mixture whose association does not cross is a *different
    fluid* that converges: ``eos.pt_flash`` sent nine arguments and no association, so
    the Rust backend ran a classical SRK flash on a fluid carrying the CPA interaction
    column and returned ``all_liquid`` where the associating model splits at
    ``beta = 0.208383589``. The twelve numbers below are the fields of
    :class:`AssociationParameters` after the scheme, in the order Rust's
    ``AssociationRecord`` declares them, in the internal scale the table states them in.
    """
    schemes: list[str] = []
    values: list[list[float]] = []
    for component in mixture.components:
        record = component.association
        if record is None:
            schemes.append("")
            values.append([0.0] * 12)
            continue
        schemes.append(record.scheme)
        values.append(
            [
                float(record.sites),
                record.energy,
                record.volume_srk,
                record.a_srk,
                record.b_srk,
                record.m_srk,
                record.volume_pr,
                record.a_pr,
                record.b_pr,
                record.m_pr,
                record.racket_z,
                record.volume_correction,
            ]
        )
    return _core.AssociationSpec(bool(mixture.associating), schemes, values)
