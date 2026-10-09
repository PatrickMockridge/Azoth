"""``characterization.tbp_cut_properties`` - a TBP cut's critical properties.

Spec: ``specs/calcs/characterization/tbp_cut_properties.toml``

Mirrors ``crates/azoth-characterization/src/tbp_cut_properties.rs`` line for line, so a
reviewer can read the two side by side. The oracle is
``validation/neqsim/captures/characterization_probe.tsv``.

# The unit boundary

NeqSim's correlations take a molar mass in **g/mol** and a specific gravity in **g/cm3**. The
declared inputs are kg/mol and kg/m3 and this converts once, at the top.

# Two upstream behaviours reproduced rather than repaired

``PedersenSRKHeavyOil`` is a **no-op**: its class re-declares the coefficient fields, shadowing
its parent's, and the parent's ``calcTC`` reads the parent's. ``TwuModel`` computes a molar mass
with ``solveMW`` and never reads it back; the call is kept and its result dropped.
"""

from __future__ import annotations

import math

from azoth._registry_gen import spec as _spec_for
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import TbpCutPropertiesResult
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning

CALC_ID = "characterization.tbp_cut_properties"

#: The reference pressure the acentric-factor correlations divide by, in bar.
REFERENCE_PRESSURE_BAR = 1.01325

#: The model's default when the input is absent, which is what ``getModel`` does with an
#: unrecognised name too.
DEFAULT_MODEL = "pedersen_srk"

#: ``[c0, c1, c2, c3, c4]`` per row - row 0 for ``tc``, row 1 for ``pc`` (``c4`` the density
#: exponent) and row 2 for the alpha exponent ``m``.
SRK_OIL = (
    (163.12, 86.052, 0.43475, -1877.4, 0.0),
    (-0.13408, 2.5019, 208.46, -3987.2, 1.0),
    (0.7431, 0.0048122, 0.0096707, -3.7184e-6, 0.0),
)
SRK_HEAVY_OIL = (
    (8.3063e2, 1.75228e1, 4.55911e-2, -1.13484e4, 0.0),
    (8.02988e-1, 1.78396, 1.56740e2, -6.96559e3, 0.25),
    (-4.7268e-2, 6.02931e-2, 1.21051, -5.76676e-3, 0.0),
)
PR_OIL = (
    (73.4043, 97.3562, 0.618744, -2059.32, 0.0),
    (0.0728462, 2.18811, 163.91, -4043.23, 0.25),
    (0.373765, 0.00549269, 0.0117934, -4.93049e-6, 0.0),
)
PR_HEAVY_OIL = (
    (9.13222e2, 1.01134e1, 4.54194e-2, -1.3587e4, 0.0),
    (1.28155, 1.26838, 1.67106e2, -8.10164e3, 0.25),
    (-2.3838e-1, 6.10147e-2, 1.32349, -6.52067e-3, 0.0),
)

PEDERSEN_HEAVY_SWITCH = 1120.0
BOILING_POINT_SWITCH = 540.0
RIAZI_DAUBERT_SWITCH = 300.0

_PEDERSEN = (
    "pedersen_srk",
    "pedersen_srk_heavy_oil",
    "pedersen_pr",
    "pedersen_pr2",
    "pedersen_pr_heavy_oil",
)


def _pedersen_coefs(model: str, molar_mass: float):
    """The coefficient set a model uses at this molar mass, in g/mol."""
    heavy = molar_mass >= PEDERSEN_HEAVY_SWITCH
    if model in ("pedersen_srk", "pedersen_srk_heavy_oil"):
        # `pedersen_srk_heavy_oil` shadows its parent's fields and changes nothing.
        return SRK_HEAVY_OIL if heavy else SRK_OIL
    if model == "pedersen_pr_heavy_oil":
        # Its constructor assigns the heavy set to *both* fields.
        return PR_HEAVY_OIL
    return PR_HEAVY_OIL if heavy else PR_OIL


def _srk_tb(molar_mass: float, density: float) -> float:
    """The boiling point the SRK family correlates, in K."""
    if molar_mass < BOILING_POINT_SWITCH:
        return 2.0e-6 * molar_mass**3 - 0.0035 * molar_mass**2 + 2.4003 * molar_mass + 171.74
    return 97.58 * molar_mass**0.3323 * density**0.04609


def _base_tb(molar_mass: float, density: float) -> float:
    """``TBPBaseModel.calcTB``'s fallback, in K."""
    return (molar_mass / 5.805e-5 * density**0.9371) ** (1.0 / 2.3776)


def _acentric_edmister(tc: float, tb: float, pc_bar: float) -> float:
    return 3.0 / 7.0 * math.log10(pc_bar / REFERENCE_PRESSURE_BAR) / (tc / tb - 1.0) - 1.0


def _acentric_kesler_lee(tc: float, tb: float, pc_bar: float, density: float) -> float:
    tbr = tb / tc
    pbr = REFERENCE_PRESSURE_BAR / pc_bar
    if tbr < 0.8:
        return (
            math.log(pbr) - 5.92714 + 6.09649 / tbr + 1.28862 * math.log(tbr) - 0.169347 * tbr**6
        ) / (15.2518 - 15.6875 / tbr - 13.4721 * math.log(tbr) + 0.43577 * tbr**6)
    kw = tb ** (1.0 / 3.0) / density
    return -7.904 + 0.1352 * kw - 0.007465 * kw * kw + 8.359 * tbr + (1.408 - 0.01063 * kw) / tbr


def _twu_tfunc(mw: float, tb: float) -> float:
    phi = math.log(mw)
    return (
        math.exp(
            5.1264 + 2.71579 * phi - 0.28659 * phi * phi - 39.8544 / phi - 0.122488 / phi / phi
        )
        - 13.7512 * phi
        + 19.6197 * phi * phi
        - tb
    )


def _twu_solve_mw(tb: float) -> float:
    """``TwuModel.solveMW``: a damped Newton with a central-difference gradient.

    The loop condition is NeqSim's own, ``abs(error) > 1e-6 and iter < 1000 or iter < 3``.
    """
    mw = tb / (5.8 - 0.0052 * tb)
    error = 1.0
    iteration = 0
    while True:
        iteration += 1
        previous = mw
        gradient = (_twu_tfunc(mw + 1.0, tb) - _twu_tfunc(mw - 1.0, tb)) / 2.0
        mw -= 0.5 * _twu_tfunc(mw, tb) / gradient
        error = abs(mw - previous)
        if not ((abs(error) > 1e-6 and iteration < 1000) or iteration < 3):
            break
    return mw


def _twu_tc_pc(density: float, tb: float) -> tuple[float, float]:
    """``TwuModel``'s critical temperature and pressure, ``(K, bar)``."""
    # Computed by the class and never read - kept, with its result dropped.
    _twu_solve_mw(tb)
    tc_n_alkane = tb / (
        0.533272
        + 0.343831e-3 * tb
        + 2.526167e-7 * tb**2
        - 1.65848e-10 * tb**3
        + 4.60774e24 * tb**-13
    )
    phi = 1.0 - tb / tc_n_alkane
    sg_n_alkane = 0.843593 - 0.128624 * phi - 3.36159 * phi**3 - 13749 * phi**12
    pc_n_alkane = (
        0.318317 + 0.099334 * phi**0.5 + 2.89698 * phi + 3.0054 * phi * phi + 8.65163 * phi**4
    ) ** 2
    vc_n_alkane = (0.82055 + 0.715468 * phi + 2.21266 * phi**3 + 13411.1 * phi**14) ** -8

    delta_st = math.exp(5.0 * (sg_n_alkane - density)) - 1.0
    f_t = delta_st * (-0.270159 * tb**-0.5 + (0.0398285 - 0.706691 * tb**-0.5) * delta_st)
    tc = tc_n_alkane * ((1 + 2 * f_t) / (1 - 2 * f_t)) ** 2

    delta_sp = math.exp(0.5 * (sg_n_alkane - density)) - 1.0
    delta_sv = math.exp(4.0 * (sg_n_alkane**2 - density**2)) - 1.0
    f_v = delta_sv * (0.347776 * tb**-0.5 + (-0.182421 + 2.24890 * tb**-0.5) * delta_sv)
    vc = vc_n_alkane * ((1 + 2 * f_v) / (1 - 2 * f_v)) ** 2
    f_p = delta_sp * (
        (2.53262 - 34.4321 * tb**-0.5 - 0.00230193 * tb)
        + (-11.4277 + 187.934 * tb**-0.5 + 0.00414963 * tb) * delta_sp
    )
    pc = (
        pc_n_alkane * (tc / tc_n_alkane) * (vc_n_alkane / vc) * ((1 + 2 * f_p) / (1 - 2 * f_p)) ** 2
    )
    # MPa to bar, the class's own conversion.
    return tc, pc * 10.0


def _lee_kesler_tc_pc(tb: float, density: float) -> tuple[float, float]:
    tc = (
        189.8
        + 450.6 * density
        + (0.4244 + 0.1174 * density) * tb
        + (0.1441 - 1.0069 * density) * 1e5 / tb
    )
    log_pc = (
        3.3864
        - 0.0566 / density
        - (0.43639 + 4.1216 / density + 0.21343 / density / density) * 1e-3 * tb
        + (0.47579 + 1.182 / density + 0.15302 / density / density) * 1e-6 * tb * tb
        - (2.4505 + 9.9099 / density / density) * 1e-10 * tb**3
    )
    return tc, math.exp(log_pc) * 10.0


def _riazi_daubert_tc_pc(molar_mass: float, density: float) -> tuple[float, float]:
    tc = (
        5.0
        / 9.0
        * 554.4
        * math.exp(-1.3478e-4 * molar_mass - 0.61641 * density)
        * molar_mass**0.2998
        * density**1.0555
    )
    pc = (
        0.068947
        * 4.5203e4
        * math.exp(-1.8078e-3 * molar_mass - 0.3084 * density)
        * molar_mass**-0.8063
        * density**1.6015
    )
    return tc, pc


def _evaluate(model: str, molar_mass: float, density: float, supplied_tb: float | None) -> dict:
    """One model over one cut, in g/mol and g/cm3."""
    tb_override = supplied_tb if (supplied_tb is not None and supplied_tb > 0.0) else None

    if model in (
        "pedersen_srk",
        "pedersen_srk_heavy_oil",
        "pedersen_pr",
        "pedersen_pr2",
        "pedersen_pr_heavy_oil",
    ):
        coefs = _pedersen_coefs(model, molar_mass)
        tc = (
            coefs[0][0] * density
            + coefs[0][1] * math.log(molar_mass)
            + coefs[0][2] * molar_mass
            + coefs[0][3] / molar_mass
        )
        pc_bar = math.exp(
            0.01325
            + coefs[1][0]
            + coefs[1][1] * density ** coefs[1][4]
            + coefs[1][2] / molar_mass
            + coefs[1][3] / molar_mass**2
        )
        exponent = (
            coefs[2][0]
            + coefs[2][1] * molar_mass
            + coefs[2][2] * density
            + coefs[2][3] * molar_mass**2
        )
        if model == "pedersen_pr2":
            tb = (
                tb_override
                if tb_override is not None
                else (
                    1928.3
                    - 1.695e5
                    * molar_mass**-0.03522
                    * density**3.266
                    * math.exp(
                        -4.922e-3 * molar_mass - 4.7685 * density + 3.462e-3 * molar_mass * density
                    )
                )
                / 1.8
            )
        else:
            tb = tb_override if tb_override is not None else _srk_tb(molar_mass, density)
        return {
            "tc": tc,
            "pc_bar": pc_bar,
            "tb": tb,
            "acentric_factor": _acentric_edmister(tc, tb, pc_bar),
            "attraction_exponent": exponent,
        }

    if model == "riazi_daubert":
        if molar_mass > RIAZI_DAUBERT_SWITCH:
            tc = (
                SRK_OIL[0][0] * density
                + SRK_OIL[0][1] * math.log(molar_mass)
                + SRK_OIL[0][2] * molar_mass
                + SRK_OIL[0][3] / molar_mass
            )
            pc_bar = math.exp(
                0.01325
                + SRK_OIL[1][0]
                + SRK_OIL[1][1] * density ** SRK_OIL[1][4]
                + SRK_OIL[1][2] / molar_mass
                + SRK_OIL[1][3] / molar_mass**2
            )
        else:
            tc, pc_bar = _riazi_daubert_tc_pc(molar_mass, density)
        tb = (
            tb_override
            if tb_override is not None
            else 97.58 * molar_mass**0.3323 * density**0.04609
        )
        return {
            "tc": tc,
            "pc_bar": pc_bar,
            "tb": tb,
            "acentric_factor": _acentric_kesler_lee(tc, tb, pc_bar, density),
            "attraction_exponent": None,
        }

    if model == "lee_kesler":
        tb = tb_override if tb_override is not None else _base_tb(molar_mass, density)
        tc, pc_bar = _lee_kesler_tc_pc(tb, density)
        return {
            "tc": tc,
            "pc_bar": pc_bar,
            "tb": tb,
            "acentric_factor": _acentric_kesler_lee(tc, tb, pc_bar, density),
            "attraction_exponent": None,
        }

    if model == "twu":
        tb = tb_override if tb_override is not None else _base_tb(molar_mass, density)
        tc, pc_bar = _twu_tc_pc(density, tb)
        return {
            "tc": tc,
            "pc_bar": pc_bar,
            "tb": tb,
            "acentric_factor": _acentric_edmister(tc, tb, pc_bar),
            "attraction_exponent": None,
        }

    if model == "cavett":
        tb = (
            tb_override
            if tb_override is not None
            else 97.58 * molar_mass**0.3323 * density**0.04609
        )
        tc_base, pc_base = _lee_kesler_tc_pc(tb, density)
        api = 141.5 / density - 131.5
        if api < 30.0:
            tc = tc_base * (1.0 + 0.002 * (30.0 - api))
            pc_bar = pc_base * (1.0 + 0.001 * (30.0 - api))
        else:
            tc, pc_bar = tc_base, pc_base
        if tb / tc >= 1.0:
            acentric = _acentric_kesler_lee(tc, tb, pc_bar, density)
        else:
            acentric = min(max(_acentric_edmister(tc, tb, pc_bar), 0.0), 1.5)
        return {
            "tc": tc,
            "pc_bar": pc_bar,
            "tb": tb,
            "acentric_factor": acentric,
            "attraction_exponent": None,
        }

    # `standing`, and the default `getModel` falls back to for anything unrecognised.
    tc, pc_bar = _riazi_daubert_tc_pc(molar_mass, density)
    tb = tb_override if tb_override is not None else _base_tb(molar_mass, density)
    return {
        "tc": tc,
        "pc_bar": pc_bar,
        "tb": tb,
        "acentric_factor": _acentric_kesler_lee(tc, tb, pc_bar, density),
        "attraction_exponent": None,
    }


def tbp_cut_properties(
    molar_mass: Q,
    density: Q,
    model: str | None = None,
    boiling_point: Q | None = None,
) -> TbpCutPropertiesResult:
    """A TBP cut's critical properties, by any of NeqSim's ten models.

    Args:
        model: which of NeqSim's ten models evaluates the cut. ``None`` takes the class's
            own default, ``pedersen_srk``.
        molar_mass: the cut's molar mass.
        density: the cut's normal liquid density at 15 C.
        boiling_point: the cut's normal boiling point, when it is known. Every model returns
            it instead of correlating one; ``None`` runs the correlation.

    Raises:
        OutOfRangeError: if ``molar_mass`` or ``density`` is not positive, or if a supplied
            ``boiling_point`` is not positive.

    Example:
        >>> import azoth
        >>> q = azoth.ureg.Quantity
        >>> r = tbp_cut_properties(
        ...     "lee_kesler", q(0.5, "kg/mol"), q(880.0, "kg/m**3")
        ... )
        >>> round(r.tc.magnitude, 6)
        906.211219
    """
    spec = _spec_for(CALC_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    values = {
        "molar_mass": input_to_si(spec, "molar_mass", molar_mass),
        "density": input_to_si(spec, "density", density),
        # Absent is a skipped check, not a pass - `apply_checks` emits RANGE_CHECK_SKIPPED.
        "boiling_point": (
            None if boiling_point is None else input_to_si(spec, "boiling_point", boiling_point)
        ),
    }
    apply_checks(checks.on_input, values.get, warnings)

    # The boundary: g/mol and g/cm3, which is what every correlation above is written in.
    molar_mass_g = values["molar_mass"] * 1.0e3
    specific_gravity = values["density"] / 1.0e3

    cut = _evaluate(
        model if model is not None else DEFAULT_MODEL,
        molar_mass_g,
        specific_gravity,
        values["boiling_point"],
    )

    watson_k = (1.8 * cut["tb"]) ** (1.0 / 3.0) / specific_gravity

    apply_checks(
        checks.derived,
        lambda name: (
            cut["acentric_factor"]
            if name == "acentric_factor"
            else cut["attraction_exponent"]
            if name == "attraction_exponent"
            else None
        ),
        warnings,
    )

    return TbpCutPropertiesResult(
        tc=from_si(cut["tc"], "K"),
        pc=from_si(cut["pc_bar"] * 1.0e5, "Pa"),
        boiling_temperature=from_si(cut["tb"], "K"),
        acentric_factor=cut["acentric_factor"],
        attraction_exponent=cut["attraction_exponent"],
        watson_k=watson_k,
        warnings=tuple(warnings),
    )
