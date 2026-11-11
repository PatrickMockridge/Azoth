"""``characterization.tbp_fraction_component`` - what ``addTBPfraction`` makes of a cut.

Spec: ``specs/models/characterization/tbp_fraction_component.toml``. Oracle:
``validation/neqsim/captures/tbp_fraction_component_probe.tsv``.

Mirrors ``crates/azoth-characterization/src/tbp_fraction_component.rs``.

# Two acentric factors, and they are not the same number

The TBP model's Edmister form and its fitted ``m`` do not imply each other. ``addTBPfraction``
sets the model's value, then ``AttractiveTermSrk.setm`` overwrites it with the root of
``-0.176 w**2 + 1.574 w + (0.48 - m)`` - so the component carries the *root*, while the critical
volume, which recomputes the model's own, is built from the Edmister value.

# The Rackett compressibility is a flash

``calcRacketZ`` reads the molar volume of a one-component SRK phase at 15 C and one atmosphere,
built from the cut's constants and fitted ``m``, and turns that volume and the Peneloux
coefficient into ``Z_RA = 0.29441 - c/(0.40768 R Tc/Pc)``. NeqSim's gas constant is
``8.3144621``; the formula is a difference of two volumes, so the constant is carried here rather
than taken from the library.
"""

from __future__ import annotations

from azoth import _models_gen
from azoth.characterization.reference import _combiner
from azoth.characterization.reference.tbp_cut_properties import tbp_cut_properties
from azoth.core.errors import PropertyUnavailableError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import TbpFractionComponentResult
from azoth.core.units import Q, from_si, ureg
from azoth.core.warnings import Warning, WarningCode
from azoth.eos import Component, mixture, pt_flash, srk_kappa
from azoth.eos.cubic import SRK

MODEL_ID = "characterization.tbp_fraction_component"

#: NeqSim's `ThermodynamicConstantsInterface.R`, in J/(mol K).
NEQSIM_R = 8.3144621

#: The reference system's temperature and pressure, `273.15 + 15` K and 1.01325 bar.
REFERENCE_TEMPERATURE = 288.15
REFERENCE_PRESSURE_BAR = 1.01325

#: `PedersenTBPModelSRK`'s `TPBracketcoefs`: the Rackett intercept and the Peneloux coefficient.
RACKETT_INTERCEPT = 0.29441
PENELOUX_COEFFICIENT = 0.40768


def _racket_z(
    molar_mass_grams: float,
    gravity: float,
    tc: Q,
    pc_bar: float,
    attraction_exponent: float,
    acentric: float,
) -> float:
    """``calcRacketZ``: the reference SRK phase's molar volume, as a Peneloux constant."""
    component = Component(
        Tc=tc,
        Pc=ureg.Quantity(pc_bar * 1.0e5, "Pa"),
        omega=acentric,
        alpha_params=(attraction_exponent,),
    )
    reference = mixture([component], cubic=SRK, alpha="srk_fitted")
    flash = pt_flash(
        reference,
        ureg.Quantity(REFERENCE_TEMPERATURE, "K"),
        ureg.Quantity(REFERENCE_PRESSURE_BAR * 1.0e5, "Pa"),
        [1.0],
    )
    volume = flash.z_liquid * NEQSIM_R * REFERENCE_TEMPERATURE / REFERENCE_PRESSURE_BAR
    peneloux = volume - molar_mass_grams / (gravity * 10.0)
    return RACKETT_INTERCEPT - peneloux / (
        PENELOUX_COEFFICIENT * NEQSIM_R * tc.to("K").magnitude / pc_bar
    )


def tbp_fraction_component(
    molar_mass: Q,
    density: Q,
) -> TbpFractionComponentResult:
    """The component ``addTBPfraction`` builds from a cut, by the class's default TBP model.

    ``molar_mass`` is the cut's molar mass and ``density`` its normal liquid density at 15 C and
    1 atm. The answer is the card a cubic reads: the critical constants, the fitted alpha
    coefficient, the Peneloux compressibility, the critical volume and viscosity, the parachor,
    the fusion and formation data, the ideal-gas heat-capacity polynomial and the Lennard-Jones
    pair.

    **The acentric factor is the one the fitted ``m`` implies**, which is not the model's Edmister
    value; the critical volume is built from the model's own, as the class does.
    """
    spec = _models_gen.model(MODEL_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    mass = float(molar_mass.to("kg/mol").magnitude)
    density_kg = float(density.to("kg/m**3").magnitude)
    apply_checks(checks.on_input, {"molar_mass": mass, "density": density_kg}.get, warnings)

    cut = tbp_cut_properties(molar_mass, density)
    # The inner model's own warnings are the cut's state; a skipped optional bound is its
    # bookkeeping and not this id's.
    warnings.extend(
        warning for warning in cut.warnings if warning.code is not WarningCode.RANGE_CHECK_SKIPPED
    )
    attraction_exponent = cut.attraction_exponent
    if attraction_exponent is None:
        raise PropertyUnavailableError(
            "m",
            "attraction exponent",
            "the default TBP model always states one, and this path is only that model",
        )
    molar_mass_grams = mass * 1000.0
    gravity = density_kg / 1000.0
    tc = cut.tc
    tc_kelvin = float(tc.to("K").magnitude)
    pc_bar = float(cut.pc.to("Pa").magnitude) / 1.0e5
    acentric = _combiner.implied_acentric_factor(attraction_exponent)
    # `init` re-derives the attractive term's coefficient from the acentric factor, so the
    # component carries Soave's form on the root rather than the model's fitted `m`.
    attraction_exponent = srk_kappa(acentric).kappa
    model_acentric = cut.acentric_factor

    racket = _racket_z(molar_mass_grams, gravity, tc, pc_bar, attraction_exponent, acentric)

    critical_volume_cm3 = (0.2918 - 0.0928 * model_acentric) * NEQSIM_R * tc_kelvin / pc_bar * 10.0
    if critical_volume_cm3 < 0.0:
        critical_volume_cm3 = (0.2918 - 0.0928) * NEQSIM_R * tc_kelvin / pc_bar * 10.0
    watson = float(cut.boiling_temperature.to("K").magnitude * 1.8) ** (1.0 / 3.0) / gravity
    correction = ((12.8 - watson) * (10.0 - watson) / (10.0 * acentric)) ** 2
    cp_a = (
        (-0.33886 + 0.02827 * watson - 0.26105 * correction + 0.59332 * acentric * correction)
        * 4.18682
        * mass
        * 1.0e3
    )
    cp_b = (
        (
            -(0.9291 - 1.1543 * watson + 0.0368 * watson * watson) * 1.0e-4
            + correction * (4.56 - 9.48 * acentric) * 1.0e-4
        )
        * 4.18682
        * mass
        * 1.8
        * 1.0e3
    )
    cp_c = (
        (-1.6658e-7 + correction * (0.536 - 0.6828 * acentric) * 1.0e-7)
        * 4.18682
        * mass
        * 1.8
        * 1.8
        * 1.0e3
    )
    triple_point = 374.5 + 0.02617 * molar_mass_grams - 20172.0 / molar_mass_grams
    pc_atmosphere = pc_bar / REFERENCE_PRESSURE_BAR

    apply_checks(checks.derived, lambda _name: None, warnings)

    return TbpFractionComponentResult(
        critical_temperature=tc,
        critical_pressure=cut.pc,
        boiling_temperature=cut.boiling_temperature,
        acentric_factor=acentric,
        attraction_exponent=attraction_exponent,
        critical_volume=from_si(critical_volume_cm3 / 1.0e6, "m**3/mol"),
        racket_z=racket,
        racket_z_cpa=racket,
        parachor=59.3 + 2.34 * molar_mass_grams,
        critical_viscosity=from_si(
            7.94830
            * molar_mass_grams**0.5
            * pc_bar ** (2.0 / 3.0)
            / tc_kelvin ** (1.0 / 6.0)
            * 1.0e-7,
            "Pa*s",
        ),
        triple_point_temperature=from_si(triple_point, "K"),
        heat_of_fusion=from_si(0.1426 / 0.238845 * molar_mass_grams * triple_point, "J/mol"),
        ideal_gas_enthalpy_of_formation=from_si(-1462600.0 * mass - 47566.0, "J/mol"),
        cp_a=from_si(cp_a, "J/(mol*K)"),
        cp_b=from_si(cp_b, "J/(mol*K**2)"),
        cp_c=from_si(cp_c, "J/(mol*K**3)"),
        # NeqSim's TBP path never sets `CpE`; `cp_d`'s getter adds `1e-10`.
        cp_d=from_si(1.0e-10, "J/(mol*K**4)"),
        lennard_jones_diameter=from_si(
            (2.3551 - 0.087 * acentric) * (tc_kelvin / pc_atmosphere) ** (1.0 / 3.0) * 1.0e-10,
            "angstrom",
        ),
        lennard_jones_energy=from_si(tc_kelvin * (0.7915 + 0.1693 * acentric * acentric), "K"),
        warnings=tuple(warnings),
    )
