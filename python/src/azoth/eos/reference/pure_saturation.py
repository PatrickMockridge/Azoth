"""``eos.pure_saturation`` - pure-component saturation pressure.

The pressure at which a pure component's vapour and liquid roots have equal fugacity,
found by bisection on reduced pressure. A *model* rather than a calculation: what the
spec pins down is the procedure, and this module reads the procedure from
``azoth._models_gen`` rather than choosing it.

Spec: ``specs/models/eos/pure_saturation.yaml``, which carries the procedure - the
bracketing rule, the tolerance and the cap - since a procedure that differs between two
implementations reaches a slightly different answer.

Every number is a registered calc's: :func:`azoth.eos.pr_kappa`,
:func:`azoth.eos.pr_alpha_ab`, :func:`azoth.eos.pr_z_factor` and
:func:`azoth.eos.pr_departure`. Only the search is here.
"""

from __future__ import annotations

from azoth import _models_gen
from azoth.core.errors import OutOfRangeError, SolverNotConvergedError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import PureSaturationResult, RootStructure
from azoth.core.solver import Convergence
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning
from azoth.eos.reference.pr_alpha_ab import pr_alpha_ab
from azoth.eos.reference.pr_departure import pr_departure
from azoth.eos.reference.pr_kappa import pr_kappa
from azoth.eos.reference.pr_z_factor import pr_z_factor

MODEL_ID = "eos.pure_saturation"


def pure_saturation(Tc: Q, Pc: Q, omega: float, T: Q) -> PureSaturationResult:
    """The saturation pressure of a pure component at a temperature.

    Args:
        Tc: critical temperature.
        Pc: critical pressure.
        omega: acentric factor. All three are the caller's; `azoth.eos.components.entry`
            resolves them from the databank by name, and this function is the scalar
            kernel underneath.
        T: absolute temperature. Must be below ``Tc``: above the critical temperature
            a pure component has no saturation pressure, and this refuses rather than
            returning a plausible-looking extrapolation.

    Returns:
        The saturation pressure, the common ``ln_phi`` at it, and the search's own
        report.

    Raises:
        OutOfRangeError: if ``T``, ``Tc`` or ``Pc`` is not positive, or if
            ``T >= Tc``. The bound is on the ratio ``T/Tc``, so the error names what
            is actually wrong rather than picking one of the two inputs arbitrarily.
        SolverNotConvergedError: if the bisection hits its cap.

    Example:
        >>> import azoth
        >>> q = azoth.ureg.Quantity
        >>> r = pure_saturation(q(369.83, "K"), q(4_248_000.0, "Pa"), 0.1523, q(300.0, "K"))
        >>> round(r.p_sat.magnitude / 1e5, 4)
        9.9791
    """
    spec = _models_gen.model(MODEL_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    # `Tc` and `Pc` are converted here rather than through `input_to_si`: that reads the
    # unit out of the spec's declaration of the input, and the spec no longer declares
    # them - a component's constants come from the databank by name, and by the time a
    # caller holds them they are already quantities in the units their constructor took.
    # `T` is still a declared input, so it still goes through the spec.
    values = {
        "Tc": Tc.to("K").magnitude,
        "Pc": Pc.to("Pa").magnitude,
        "omega": omega,
        "T": input_to_si(spec, "T", T),
    }

    # The positivity of `Tc` and `Pc` is refused here rather than by a bound in the
    # spec. They used to be declared `valid_range` entries: the model took the
    # constants as arguments, so a bound on them was a bound on an input. The model
    # names a *component* now and the constants come from the databank, which leaves
    # no input for a bound to guard. `Component.__post_init__` refuses a non-positive
    # pair on the mixture path; this is the same refusal on the scalar one, and it is
    # what keeps `Tc` from being a divisor by zero.
    if not values["Tc"] > 0.0:
        raise OutOfRangeError(
            "Tc",
            values["Tc"],
            "a critical temperature is absolute and positive by definition; it is also "
            "a divisor in the reduced temperature",
        )
    if not values["Pc"] > 0.0:
        raise OutOfRangeError(
            "Pc",
            values["Pc"],
            "a critical pressure is positive by definition, and it converts the reduced "
            "answer back to pascals",
        )

    apply_checks(checks.on_input, values.get, warnings)

    reduced_temperature = values["T"] / values["Tc"]
    apply_checks(
        checks.derived,
        lambda name: reduced_temperature if name == "t_over_tc" else None,
        warnings,
    )

    algorithm = spec["algorithm"]
    bracket = algorithm["bracket"]
    kappa = pr_kappa(omega).kappa

    def pair_at(pr: float) -> tuple[float, float] | None:
        """The two extreme fugacities at a trial reduced pressure, or `None`.

        `None` when the cubic has no liquid branch there - which is exactly what
        `pr_z_factor`'s `root_structure` reports, so this does not second-guess the
        kernel's admissibility rule.
        """
        ab = pr_alpha_ab(kappa, reduced_temperature, pr)
        roots = pr_z_factor(ab.a_reduced, ab.b_reduced)
        if roots.root_structure is not RootStructure.THREE_ROOTS:
            return None
        liquid = pr_departure(ab.a_reduced, ab.b_reduced, roots.z_min, kappa, reduced_temperature)
        vapour = pr_departure(ab.a_reduced, ab.b_reduced, roots.z_max, kappa, reduced_temperature)
        return (liquid.ln_phi, vapour.ln_phi)

    # 1. Bracket. Linear and with the spec's step count, so both implementations land
    #    on the same point - the scan's resolution is the one thing that moves the
    #    answer, and it is in the spec for that reason.
    upper: float | None = None
    for step in range(bracket["steps"]):
        fraction = step / (bracket["steps"] - 1)
        pr = bracket["lower"] + (bracket["upper"] - bracket["lower"]) * fraction
        if pair_at(pr) is not None:
            upper = pr
    if upper is None:  # pragma: no cover - guarded by the t_over_tc bound
        raise OutOfRangeError(
            "t_over_tc",
            reduced_temperature,
            f"no pressure between Pr = {bracket['lower']} and Pr = {bracket['upper']} "
            f"has a liquid branch, so there is nothing to equate. A below-critical "
            f"temperature should always have one, so this means the scan's bounds are "
            f"wrong rather than the state",
        )
    lo, hi = bracket["lower"], upper
    convergence = Convergence.parse(algorithm["convergence"])

    # 2. Bisect on the bracket's width.
    iterations = 0
    for step in range(1, algorithm["max_iterations"] + 1):
        iterations = step
        mid = 0.5 * (lo + hi)
        width = hi - lo
        if convergence is Convergence.ABSOLUTE:
            close_enough = width <= algorithm["tolerance"]
        else:
            close_enough = width <= algorithm["tolerance"] * abs(mid)
        if close_enough:
            break
        pair = pair_at(mid)
        if pair is None:  # pragma: no cover - the bracket's top is a three-root point
            break
        if pair[0] - pair[1] > 0.0:
            lo = mid
        else:
            hi = mid

    reduced_pressure = 0.5 * (lo + hi)
    residual = (hi - lo) / (2.0 * abs(reduced_pressure))
    final = pair_at(reduced_pressure)
    if final is None:  # pragma: no cover - the same
        raise SolverNotConvergedError(iterations, residual, algorithm["tolerance"])

    return PureSaturationResult(
        p_sat=from_si(reduced_pressure * values["Pc"], "Pa"),
        ln_phi=final[0],
        iterations=iterations,
        residual=residual,
        warnings=tuple(warnings),
    )
