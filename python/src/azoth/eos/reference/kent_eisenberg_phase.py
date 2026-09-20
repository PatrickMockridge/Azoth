"""``eos.kent_eisenberg_phase`` - the fugacity coefficients of a Kent-Eisenberg phase.

Spec: ``specs/models/eos/kent_eisenberg_phase.toml``. A *direct* model: no iteration, so no
algorithm block.

This is the pure-Python reference: a second, independent expression of the same physics as
the Rust kernel, held to it by ``test_cross_impl.py``. What is here is the branch dispatch -
a component's ``REFERENCESTATETYPE`` and charge decide which of three expressions its
fugacity coefficient is, and the activity coefficient is one for all of them.
"""

from __future__ import annotations

import math
from collections.abc import Sequence

from azoth import _models_gen
from azoth.core.errors import InvalidInputError, PropertyUnavailableError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import KentEisenbergPhaseResult
from azoth.core.units import Q, input_to_si
from azoth.core.warnings import Warning
from azoth.eos import components as _components
from azoth.eos.reference import _henry
from azoth.eos.reference.antoine_vapor_pressure import antoine_vapor_pressure

MODEL_ID = "eos.kent_eisenberg_phase"

#: The fugacity coefficient an ion gets: ``ComponentKentEisenberg.fugcoef``'s constant.
#:
#: **This model's own**, not a shared one: NeqSim gives `1e12` in `ComponentGePitzer` and
#: `1e-15` in `ComponentDesmukhMather`.
INSOLUBLE_ION = 1.0e8


def kent_eisenberg_phase(
    components: Sequence[str],
    T: Q,
    P: Q,
    x: Sequence[float],
) -> KentEisenbergPhaseResult:
    """The fugacity coefficients of a phase whose activity coefficients are one.

    ``PhaseKentEisenberg`` overrides ``getActivityCoefficient`` to return ``1.0`` for every
    component, so this model's whole content is the branch ``ComponentKentEisenberg.fugcoef``
    takes: a ``solvent`` reference state gets ``P0_i(T)/P``, a neutral with any other
    reference state gets ``H_i(T)/P``, and an ion gets :data:`INSOLUBLE_ION`.

    Args:
        components: the substances the phase is made of, by name.
        T: absolute temperature.
        P: absolute pressure.
        x: the phase's mole fractions; non-negative and summing to one. They enter no
            branch and are checked because a caller supplying a composition means it.

    Returns:
        The activity coefficients - all one - and the fugacity coefficients.

    Raises:
        InvalidInputError: if a component is not in the databank, if ``x`` is not a
            composition, or if a ``solvent``-reference component carries no Antoine
            correlation.
        PropertyUnavailableError: if a component's Henry row is the all-zero filler.
        OutOfRangeError: if ``T`` or ``P`` is not positive.

    Example:
        >>> import azoth
        >>> q = azoth.ureg.Quantity
        >>> r = kent_eisenberg_phase(
        ...     ["water", "na+", "cl-", "co2"], q(313.15, "K"), q(500000.0, "Pa"),
        ...     [0.89, 0.04, 0.04, 0.03],
        ... )
        >>> r.gamma
        (1.0, 1.0, 1.0, 1.0)
        >>> round(r.ln_phi[1], 9)
        18.420680744
    """
    spec = _models_gen.model(MODEL_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    t_si = input_to_si(spec, "T", T)
    p_si = input_to_si(spec, "P", P)
    apply_checks(checks.on_input, {"T": t_si, "P": p_si}.get, warnings)

    n = len(components)
    if len(x) != n:
        raise InvalidInputError("x", f"a mixture of {n} components has {len(x)} mole fractions")
    bad = next((i for i, value in enumerate(x) if value < 0.0), None)
    if bad is not None:
        raise InvalidInputError("x", f"x[{bad}] is {x[bad]} but a mole fraction cannot be negative")
    total = sum(x)
    if abs(total - 1.0) > 1.0e-9:
        raise InvalidInputError(
            "x",
            f"the mole fractions sum to {total}, not to one. Renormalising them here "
            "would make a composition error invisible in every number downstream, so it "
            "is refused instead",
        )

    ln_phi: list[float] = []
    for name in components:
        entry = _components.entry(name)
        if entry.reference_state == _components.SOLVENT:
            coefficient = _vapour_pressure(entry, t_si, warnings) / p_si
        elif entry.ionic_charge == 0.0:
            # `_henry` reports in bar and this phase's `P` is in pascals, so the
            # conversion is here rather than in the correlation.
            coefficient = _henry.effective_coefficient(entry, t_si) * 1.0e5 / p_si
        else:
            coefficient = INSOLUBLE_ION
        if not coefficient > 0.0 or not math.isfinite(coefficient):
            raise PropertyUnavailableError(
                entry.name,
                "fugacity coefficient",
                f"the branch its reference state selects gives {coefficient}, which is not "
                f"a coefficient. NeqSim returns it without comment",
            )
        ln_phi.append(math.log(coefficient))

    return KentEisenbergPhaseResult(
        # **Identically one**, which is the model's defining property rather than a
        # computed column - `PhaseKentEisenberg.getActivityCoefficient` returns it.
        gamma=tuple([1.0] * n),
        ln_gamma=tuple([0.0] * n),
        ln_phi=tuple(ln_phi),
        warnings=tuple(warnings),
    )


def _vapour_pressure(
    entry: _components.DatabankEntry, temperature_k: float, warnings: list[Warning]
) -> float:
    """``P0_i(T)`` from the component's own Antoine row, in pascals."""
    if entry.antoine is None:
        raise PropertyUnavailableError(
            entry.name,
            "Antoine vapour-pressure coefficients",
            "its reference state is `solvent`, so `ComponentKentEisenberg.fugcoef` gives it "
            "`P0_i(T) / P` and there is no correlation to evaluate. A keycard supplies the "
            "parameters a cubic reads and not these",
        )
    form = _components.form_from_type(entry.antoine_type, entry.antoine[4])
    if not form:
        raise PropertyUnavailableError(
            entry.name,
            "Antoine vapour-pressure coefficients",
            f"its row is marked `{entry.antoine_type}`, which upstream retracted - the rows "
            f"that carried one repeated filler tuple. NeqSim evaluates the filler and "
            f"returns a number; this library refuses it",
        )
    result = antoine_vapor_pressure(*entry.antoine, form, entry.Tc, entry.Pc, _q(temperature_k))
    warnings.extend(result.warnings)
    return result.p_sat.to_base_units().magnitude


def _q(value: float) -> Q:
    """A bare kelvin quantity, for the correlation's own argument."""
    import azoth

    return azoth.ureg.Quantity(value, "K")
