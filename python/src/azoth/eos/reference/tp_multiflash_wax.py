"""``eos.tp_multiflash_wax`` - how much of a feed is wax at a state.

```text
seed   the two-phase flash, plus a wax phase
solve  Q(beta) = sum_k beta_k - sum_i z_i ln E_i      E_i = sum_k beta_k / phi_ik
keep   the wax only where the solve converges and it is above the floor
```

Spec: ``specs/models/eos/tp_multiflash_wax.toml``, which carries the seeding, the floor's
rule and the states it is checked at.

NeqSim's ``TPmultiflashWAX``. **It adds one phase to
:mod:`azoth.eos.reference.tp_multiflash` and reuses the solve** rather than restating it: the
fraction Newton is Michelsen's on the whole vector, and a phase whose coefficients are known
is a phase it can carry whatever their source.

**A substance that is not a wax former is excluded by a number**: NeqSim's
``ComponentWax.fugcoef`` returns ``1e50`` for one, and ``x_i = z_i/(E_i phi_i)`` is then zero
to any precision the answer is read at.
"""

from __future__ import annotations

from azoth import _models_gen
from azoth.core.errors import InvalidInputError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import TpMultiflashWaxResult
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning
from azoth.eos.cubic import PR
from azoth.eos.mixture import Mixture
from azoth.eos.reference._mixture_state import reduced_parameters
from azoth.eos.reference.pt_flash import pt_flash
from azoth.eos.reference.tp_multiflash import _Phase, _phases_of, _solve_phase_fractions
from azoth.eos.reference.wax_solid_fugacity import wax_solid_fugacity

MODEL_ID = "eos.tp_multiflash_wax"

#: NeqSim's marker for a component that cannot be in a wax phase.
NOT_A_WAX_FORMER = 1.0e50

#: Below this a phase is not there, the same floor ``eos.tp_multiflash``'s merge uses.
WAX_FLOOR = 1.1e-12


def _wax_coefficients(mixture: Mixture, T: Q, P: Q) -> list[float]:
    """The wax solid's fugacity coefficient for every component at a state.

    Raises:
        InvalidInputError: if a wax former carries no molar mass, no heat of fusion or no
            triple-point temperature, each of which the model reads rather than defaults.
    """
    eos = "pr" if mixture.cubic is PR else "srk"
    out: list[float] = []
    for index, component in enumerate(mixture.components):
        if not component.wax_former:
            out.append(NOT_A_WAX_FORMER)
            continue
        if (
            component.molar_mass is None
            or component.heat_of_fusion <= 0.0
            or component.triple_point_temperature <= 0.0
        ):
            raise InvalidInputError(
                "components",
                f"component {index} is a wax former and carries no melt data; the wax model "
                "reads it rather than defaulting, and `eos.tbp_fraction_properties` is what "
                "gives a cut one",
            )
        out.append(
            wax_solid_fugacity(
                # `Tc`, `Pc` and the molar mass are already quantities; the melt data is a
                # float out of the table and is given its unit here.
                molar_mass=component.molar_mass,
                tc=component.Tc,
                pc=component.Pc,
                omega=component.omega,
                heat_of_fusion=from_si(component.heat_of_fusion, "J/mol"),
                triple_point_temperature=from_si(component.triple_point_temperature, "K"),
                T=T,
                P=P,
                eos=eos,
            ).fugacity_coefficient
        )
    return out


def tp_multiflash_wax(
    components: list[str], T: Q, P: Q, z: list[float], eos: str = "srk"
) -> TpMultiflashWaxResult:
    """The fraction of a feed that is wax at a temperature and pressure.

    **The components cross by name**, because the wax flag and the melt data come from the
    databank's own columns and a :class:`~azoth.eos.mixture.Component` carries no name.

    Args:
        components: the substances, by name. At least one must be a wax former.
        T: absolute temperature.
        P: absolute pressure.
        z: the overall mole fractions. Checked rather than renormalised.
        eos: the cubic the fluid runs, which is also the one each wax former's reference
            liquid is built from.

    Returns:
        The wax fraction, the split, and ``converged`` - which is False where the answer is
        the two-phase flash's rather than the three-phase solve's.

    Raises:
        InvalidInputError: if a wax former carries no melt data.
        OutOfRangeError: if ``T`` or ``P`` is not positive.
    """
    from azoth.eos import components as databank

    spec = _models_gen.model(MODEL_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []
    t_si = input_to_si(spec, "T", T)
    p_si = input_to_si(spec, "P", P)
    apply_checks(checks.on_input, {"T": t_si, "P": p_si}.get, warnings)

    mixture, _ = databank.mixture_of(components, eos=eos)
    reduced = reduced_parameters(mixture, t_si, p_si)
    algorithm = spec["algorithm"]
    tolerance = float(algorithm["tolerance"])
    cap = int(algorithm["max_iterations"])

    # **The two-phase flash seeds it**, and is also the answer if the wax does not survive.
    flash = pt_flash(mixture, T, P, z)
    split_of_flash = flash.beta if flash.beta is not None else 0.0

    def two_phase(iterations: int, residual: float, converged: bool) -> TpMultiflashWaxResult:
        return TpMultiflashWaxResult(
            wax_fraction=0.0,
            phase_count=2,
            beta=(1.0 - split_of_flash, split_of_flash),
            x=(tuple(flash.x), tuple(flash.y)),
            iterations=iterations,
            residual=residual,
            converged=converged,
            warnings=tuple(warnings),
        )

    # A feed with no wax former in it has no solid to find, and the two-phase answer is the
    # whole of it.
    if not any(component.wax_former for component in mixture.components):
        return two_phase(0, 0.0, True)

    coefficients = _wax_coefficients(mixture, T, P)
    phases = [
        *_phases_of(flash, reduced, mixture.kij, z),
        # The wax's starting share is the floor: there is no trial composition to make,
        # because the phase's coefficients do not depend on one.
        _Phase(WAX_FLOOR, list(z), liquid=True, wax=True),
    ]
    solved, iterations, residual, converged = _solve_phase_fractions(
        reduced, mixture.kij, z, phases, tolerance, cap, wax_coefficients=coefficients
    )
    wax = solved[2].fraction
    if not converged or wax <= WAX_FLOOR:
        return two_phase(iterations, residual, converged)

    return TpMultiflashWaxResult(
        wax_fraction=wax,
        phase_count=len(solved),
        beta=tuple(phase.fraction for phase in solved),
        x=tuple(tuple(phase.composition) for phase in solved),
        iterations=iterations,
        residual=residual,
        converged=converged,
        warnings=tuple(warnings),
    )


__all__ = ["tp_multiflash_wax"]
