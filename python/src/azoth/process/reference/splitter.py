"""``process.splitter`` - one feed divided into branches.

Spec: ``specs/models/process/splitter.yaml``, which carries the provenance, the reason
there is no ideal-gas block, and the divergence over normalising the fractions.

Every branch leaves at the feed's own temperature, pressure and composition, so one
isothermal flash describes all of them. The model runs it **once** and reports one
``phase`` and one ``beta``, with the branch flows as a vector - not a shortcut, but the
same answer stated in a form that cannot disagree with itself.
"""

from __future__ import annotations

from azoth import _models_gen
from azoth.core.errors import InvalidInputError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import SplitterResult
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning
from azoth.eos.mixture import Mixture
from azoth.eos.reference.pt_flash import pt_flash

MODEL_ID = "process.splitter"

#: How far the fractions may be from summing to one before the caller is told. Not zero:
#: a caller writing thirds as decimals cannot sum them exactly, and refusing that would
#: be refusing arithmetic rather than an error.
FRACTION_TOLERANCE = 1.0e-9


def splitter(
    mixture: Mixture,
    T: Q,
    P: Q,
    n: Q,
    z: list[float],
    fractions: list[float],
) -> SplitterResult:
    """One feed divided into branches at the feed's own temperature and pressure.

    Args:
        mixture: the components and their interaction parameters.
        T: the feed's absolute temperature, and every branch's.
        P: the feed's absolute pressure, and every branch's.
        n: the feed's molar flow rate, which the fractions divide.
        z: the feed's mole fractions, and every branch's.
        fractions: the share each branch takes, summing to one.

    Returns:
        The branch flows, and the phase state they all share.

    Raises:
        InvalidInputError: if ``fractions`` is empty, holds a negative entry, or does not
            sum to one to within :data:`FRACTION_TOLERANCE`.
        OutOfRangeError: if the feed's state is non-positive, or its flow is not positive.

    Note:
        There is **no ideal-gas argument**, and that is a statement about the model rather
        than an omission: a splitter does no energy balance, so it needs no datum. It is
        the only unit operation here that takes none.
    """
    spec = _models_gen.model(MODEL_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    t_si = input_to_si(spec, "T", T)
    p_si = input_to_si(spec, "P", P)
    n_si = input_to_si(spec, "n", n)
    apply_checks(checks.on_input, {"T": t_si, "P": p_si, "n": n_si}.get, warnings)

    if not fractions:
        raise InvalidInputError(
            "fractions", "a splitter with no branches has nothing to split into"
        )
    negative = next((f for f in fractions if f < 0.0), None)
    if negative is not None:
        raise InvalidInputError(
            "fractions",
            f"a negative share ({negative}) asks a branch to give back more than it "
            f"takes, which is a mixer written backwards",
        )
    total = sum(fractions)
    if abs(total - 1.0) > FRACTION_TOLERANCE:
        raise InvalidInputError(
            "fractions",
            f"the shares must account for the whole feed and sum to {total}, not 1. A "
            f"splitter that loses or invents material is not a splitter, and normalising "
            f"here would hide the arithmetic error that produced this",
        )

    # One flash: every branch is at the same state by construction.
    flash = pt_flash(mixture, from_si(t_si, "K"), from_si(p_si, "Pa"), list(z))
    warnings.extend(flash.warnings)

    return SplitterResult(
        T=from_si(t_si, "K"),
        P=from_si(p_si, "Pa"),
        phase=flash.phase,
        beta=flash.beta,
        flows=tuple(fraction * n_si for fraction in fractions),
        iterations=flash.iterations,
        warnings=tuple(warnings),
    )
