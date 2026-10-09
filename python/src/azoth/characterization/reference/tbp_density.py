"""``characterization.tbp_density`` - a cut's specific gravity from its boiling point and molar
mass.

Spec: ``specs/models/characterization/tbp_density.toml``. Oracle:
``validation/neqsim/captures/tbp_closure_probe.tsv``.

Mirrors ``crates/azoth-characterization/src/tbp_density.rs`` line for line.

# One member of four

Only the 1980 Riazi-Daubert pair inverts for specific gravity. The other three are refused by
name, each an ``[[unported]]`` row so ``tools/check_unported.py`` holds this half and the Rust
one to the same three keys.

# The pair is not symmetric

``calcMolarMass`` divides by 1000 to reach kg/mol and this multiplies by 1000 to reach g/mol, so
the two directions cannot share a conversion; sharing one is out by ``10^6`` on a density.
"""

from __future__ import annotations

from azoth import _models_gen
from azoth.characterization.reference import _unported
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import TbpDensityResult
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning

MODEL_ID = "characterization.tbp_density"

DEFAULT_CLOSURE = "riazi_daubert_1980"


def _riazi_daubert_1980(boiling_point: float, molar_mass: float) -> float:
    """The 1980 pair, rearranged for specific gravity. Returns g/cm3.

    NeqSim writes the exponent as ``-1.0 / -1.0164``; the two negatives cancel and it is written
    here as the one positive power it is.
    """
    molar_mass_gmol = molar_mass * 1000.0
    rankine = boiling_point * 1.8
    return (4.5673e-5 * rankine**2.1962 / molar_mass_gmol) ** (1.0 / 1.0164)


def tbp_density(
    boiling_point: Q,
    molar_mass: Q,
    closure: str | None = None,
) -> TbpDensityResult:
    """A cut's normal liquid density, from its normal boiling point and its molar mass.

    ``closure`` is one of the four ``TbpClosure`` members. Absent takes
    ``riazi_daubert_1980``, the only member that supports this direction; the other three are
    refused, each naming the class that would close it.

    Raises:
        InvalidInputError: for the three closures that cannot be inverted.
        OutOfRangeError: if `boiling_point` or `molar_mass` is not positive.
    """
    spec = _models_gen.model(MODEL_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    values = {
        "boiling_point": input_to_si(spec, "boiling_point", boiling_point),
        "molar_mass": input_to_si(spec, "molar_mass", molar_mass),
    }
    apply_checks(checks.on_input, values.get, warnings)

    selected = closure if closure is not None else DEFAULT_CLOSURE
    # One literal per row, so the checker reads the key without following control flow - a key a
    # reader has to compute is one `tools/check_unported.py` cannot compare against the spec.
    if selected == "riazi_daubert_1987":
        raise _unported.refuse("closure=riazi_daubert_1987")
    if selected == "soreide":
        raise _unported.refuse("closure=soreide")
    if selected == "tbp_model":
        raise _unported.refuse("closure=tbp_model")

    # The correlation answers in g/cm3 and this declares kg/m3.
    density = _riazi_daubert_1980(values["boiling_point"], values["molar_mass"])

    apply_checks(checks.derived, lambda _name: None, warnings)

    return TbpDensityResult(
        density=from_si(density * 1000.0, "kg/m**3"),
        warnings=tuple(warnings),
    )
