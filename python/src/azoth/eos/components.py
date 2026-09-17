"""The component databank, and looking a substance up by name.

`data/components/components.csv` ships with the library, generated from NeqSim
(Equinor/NTNU, Apache-2.0) by `tools/gen_databank.py`; `NOTICE` carries the
attribution. It is how the `eos` namespace gets `Tc`, `Pc` and `omega` without every
calculation taking them as caller arguments.

    >>> from azoth.eos import component, from_names
    >>> methane = component("methane")
    >>> methane.Tc
    <Quantity(190.56, 'kelvin')>
    >>> fluid = from_names(["methane", "n-butane"])
    >>> len(fluid)
    2

# What is in it, and what is not

173 substances: NeqSim's `COMP.csv` filtered to the types a cubic equation of state
can describe. The ions are excluded because a cubic has no notion of one and NeqSim
fills their critical properties with a shared default - 29 rows carrying the same
`Pc`, `omega` and `Vc` is not a coincidence, and shipping them would ship
plausible-looking wrong numbers. `COMP_EXT.csv` is not vendored at all: 86 MB of
heavy fluids this library cannot characterise.

# Looking up is a call the caller makes

:class:`~azoth.eos.mixture.Component` still carries **no name**, and no calculation
gained a name parameter. A caller who wants a substance by name calls
:func:`component` and passes the result, which keeps the lookup visible in their own
code rather than hidden inside a flash. That was the property `mixture.py` was
protecting when it said a name would tempt something to look a value up; the answer
is not to refuse names, it is to make the lookup an explicit step.

# Provenance

Every row cites the same source. It is deliberately not a per-value citation: the
provenance here is institutional - a project, a version, a file - and a per-row URL
for 173 rows would be a citation-shaped thing that is not a citation.
"""

from __future__ import annotations

import csv
import io
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, replace
from functools import cache

from azoth import keycard
from azoth._data import find
from azoth.core.errors import InvalidInputError, PropertyUnavailableError
from azoth.core.units import Q, ureg
from azoth.eos.cubic import CUBICS, Cubic
from azoth.eos.mixture import Component, Mixture
from azoth.eos.reference.molar_enthalpy_entropy import IdealGasModel

COMPONENTS_CSV = "data/components/components.csv"
KIJ_CSV = "data/components/kij.csv"
UNIFAC_COMP_CSV = "data/components/UNIFACcomp.csv"
UNIFAC_GROUP_CSV = "data/components/UNIFACGroupParam.csv"
UNIFAC_INTER_CSV = "data/components/UNIFACInterParam.csv"
MBWR32_CSV = "data/components/mbwr32.csv"

#: Columns the loader reads, in order. Named rather than positional because this
#: file's shape is the generator's contract, and a column inserted in the middle
#: should fail loudly rather than shift every value one field left.
COLUMNS = (
    "name",
    "cas",
    "formula",
    "molar_mass_kg_per_mol",
    "tc_k",
    "pc_pa",
    "acentric_factor",
    "critical_volume_m3_per_mol",
    "liquid_density_kg_per_m3",
    "cpa",
    "cpb",
    "cpc",
    "cpd",
    "cpe",
    "schwartzentruber1",
    "schwartzentruber2",
    "schwartzentruber3",
    "mc1",
    "mc2",
    "mc3",
    "mcpr1",
    "mcpr2",
    "mcpr3",
    "umrcpa_mc1",
    "umrcpa_mc2",
    "umrcpa_mc3",
    "umrcpa_mc4",
    "umrcpa_mc5",
    "antoine_type",
    "antoinea",
    "antoineb",
    "antoinec",
    "antoined",
    "antoinee",
    "dipole_moment_debye",
    "viscosity_correction_factor",
    "citation",
)


@dataclass(frozen=True, slots=True)
class DatabankEntry:
    """One substance, as the file records it.

    Wider than :class:`~azoth.eos.mixture.Component` on purpose: the file carries
    the CAS number, the formula and the liquid density, which are worth reporting
    and are not part of what a cubic takes.

    **Those wider fields are optional, and a keycard is why.** A keycard supplies the
    parameters a *model* needs - ``Tc``, ``Pc`` and ``omega`` - and nothing else,
    because nothing else is required to run a cubic. A substance a keycard adds
    therefore has no CAS number here, and ``None`` says so rather than a blank string
    standing in for one. Filling those in from a similar substance would be inventing
    data, which is the failure mode this whole module is arranged against.
    """

    name: str
    cas: str | None
    formula: str | None
    Tc: Q
    Pc: Q
    omega: float
    molar_mass: Q | None
    critical_volume: Q | None
    liquid_density: Q | None
    #: The five coefficients of NeqSim's `Cp` polynomial, in J/(mol*K**n), or ``None``
    #: for a substance a keycard supplied. A keycard gives the parameters a *cubic*
    #: needs - `Tc`, `Pc` and `omega` - and this is not one of them, so a substance it
    #: adds has no enthalpy until its coefficients are supplied too. ``None`` says that
    #: rather than a row of zeros standing in for a polynomial.
    cp: tuple[float, float, float, float, float] | None
    #: Fitted parameters for the alpha correlations that need them, keyed by the
    #: correlation's short name. Empty for a substance a keycard supplied, which has
    #: no such columns.
    alpha_params: dict[str, tuple[float, ...]]
    #: The five Antoine coefficients ``A``-``E``, NeqSim's internal scale. A
    #: keycard-supplied substance has none, so the tuple is all zero.
    antoine: tuple[float, float, float, float, float]
    #: The raw `AntoineVapPresLiqType` label, the messy upstream vocabulary
    #: :func:`form_from_type` cleans.
    antoine_type: str
    #: Dipole moment in debye, and NeqSim's viscosity correction factor. A
    #: keycard-supplied substance has none, so both are zero.
    dipole_moment_debye: float
    viscosity_correction_factor: float
    citation: str | None
    #: Where these values came from: the vendored databank, or the keycard in force.
    #: Not part of a citation - it is the *provenance of the lookup*, which a caller
    #: needs when a result turns out to depend on which file was in play.
    source: str = "databank"

    def component(self, alpha: str | None = None) -> Component:
        """This entry as the thing a calculation takes.

        ``alpha`` names the alpha correlation, so the entry can attach the fitted
        parameters that correlation reads. Delft (1998) is the one correlation whose
        parameter is a flag rather than a fitted value: methane gets the fitted
        cubic, everything else the Soave form.
        """
        params: tuple[float, ...] = ()
        if alpha is not None:
            if alpha == "delft1998":
                params = (1.0,) if self.name == "methane" else ()
            else:
                params = self.alpha_params.get(alpha, ())
        return Component(
            Tc=self.Tc,
            Pc=self.Pc,
            omega=self.omega,
            molar_mass=self.molar_mass,
            alpha_params=params,
        )

    def antoine_form(self) -> str:
        """The cleaned Antoine form this entry's coefficients belong to, one of
        ``"pow10"``, ``"pow10kpa"``, ``"exp"`` or ``"wagner"``."""
        return form_from_type(self.antoine_type)

    def __repr__(self) -> str:
        return f"DatabankEntry({self.name!r}, Tc={self.Tc}, Pc={self.Pc}, omega={self.omega})"


def _rows(text: str) -> list[dict[str, str]]:
    body = [
        line for line in text.splitlines() if line.strip() and not line.lstrip().startswith("#")
    ]
    return [dict(row) for row in csv.DictReader(io.StringIO("\n".join(body)))]


@cache
def _table() -> dict[str, DatabankEntry]:
    """Every component, keyed by lower-case name.

    Cached because the file is immutable within a process and every lookup would
    otherwise re-read and re-parse it.
    """
    entries: dict[str, DatabankEntry] = {}
    for row in _rows(find(COMPONENTS_CSV).read_text(encoding="utf-8")):
        missing = [column for column in COLUMNS if column not in row]
        if missing:
            raise InvalidInputError(
                COMPONENTS_CSV, f"a row is missing {missing}; the file is malformed"
            )
        name = row["name"]
        entries[name] = DatabankEntry(
            name=name,
            cas=row["cas"],
            formula=row["formula"],
            Tc=ureg.Quantity(float(row["tc_k"]), "K"),
            Pc=ureg.Quantity(float(row["pc_pa"]), "Pa"),
            omega=float(row["acentric_factor"]),
            molar_mass=ureg.Quantity(float(row["molar_mass_kg_per_mol"]), "kg/mol"),
            critical_volume=ureg.Quantity(float(row["critical_volume_m3_per_mol"]), "m**3/mol"),
            liquid_density=ureg.Quantity(float(row["liquid_density_kg_per_m3"]), "kg/m**3"),
            cp=(
                float(row["cpa"]),
                float(row["cpb"]),
                float(row["cpc"]),
                float(row["cpd"]),
                float(row["cpe"]),
            ),
            alpha_params={
                "schwartzentruber": _three(
                    row, "schwartzentruber1", "schwartzentruber2", "schwartzentruber3"
                ),
                "mollerup": _three(
                    row, "schwartzentruber1", "schwartzentruber2", "schwartzentruber3"
                ),
                "matcop": _three(row, "mc1", "mc2", "mc3"),
                "matcop_pr": _three(row, "mcpr1", "mcpr2", "mcpr3"),
                # `matcop_prumr` reads no databank columns: NeqSim's term 17 is
                # constructed without parameters, so it reduces to the base Soave form.
                "matcop_5prumr": (
                    float(row["umrcpa_mc1"]),
                    float(row["umrcpa_mc2"]),
                    float(row["umrcpa_mc3"]),
                    float(row["umrcpa_mc4"]),
                    float(row["umrcpa_mc5"]),
                ),
            },
            antoine=(
                float(row["antoinea"]),
                float(row["antoineb"]),
                float(row["antoinec"]),
                float(row["antoined"]),
                float(row["antoinee"]),
            ),
            antoine_type=row["antoine_type"],
            dipole_moment_debye=float(row["dipole_moment_debye"]),
            viscosity_correction_factor=float(row["viscosity_correction_factor"]),
            citation=row["citation"],
        )
    return entries


def _three(row: Mapping[str, str], a: str, b: str, c: str) -> tuple[float, float, float]:
    return (float(row[a]), float(row[b]), float(row[c]))


def form_from_type(label: str) -> str:
    """Map NeqSim's raw ``AntoineVapPresLiqType`` label onto the clean form name.

    ``exp`` and ``log`` are one formula under two names, and ``loglog``/``log10``
    have no branch in NeqSim's dispatch, so they fall through to Wagner - a defect
    this reproduces rather than silently repairs.
    """
    if label == "pow10":
        return "pow10"
    if label == "pow10KPa":
        return "pow10kpa"
    if label in ("exp", "log"):
        return "exp"
    return "wagner"


@cache
def _kij() -> dict[tuple[str, str], float]:
    """Interaction parameters, keyed by the ordered pair of names."""
    pairs: dict[tuple[str, str], float] = {}
    for row in _rows(find(KIJ_CSV).read_text(encoding="utf-8")):
        pair = (row["component_a"], row["component_b"])
        value = float(row["kij_pr"])
        pairs[pair] = value
        pairs[(pair[1], pair[0])] = value
    return pairs


@cache
def _nrtl() -> dict[tuple[str, str], tuple[float, float]]:
    """NRTL `(alpha, gij)`, keyed by ordered pair and stored both ways round.

    `gij` is the Kelvin energy for the ordered key `(first, second)`, so the reversed
    key carries the reversed energy; `alpha` is symmetric and appears under both keys.
    """
    pairs: dict[tuple[str, str], tuple[float, float]] = {}
    for row in _rows(find(KIJ_CSV).read_text(encoding="utf-8")):
        a = row["component_a"]
        b = row["component_b"]
        alpha = float(row["nrtlalpha"])
        pairs[(a, b)] = (alpha, float(row["nrtlgij"]))
        pairs[(b, a)] = (alpha, float(row["nrtlgji"]))
    return pairs


def available(*, card: keycard.Keycard | None = None) -> tuple[str, ...]:
    """Every name available, sorted: the databank plus whatever a card adds.

    `card` is the card this call reads; with none passed the library reads the data
    it ships.
    """
    extra = set(card.components) if card is not None else set()
    return tuple(sorted(set(_table()) | extra))


#: NeqSim's Wilke-Chang association parameters, keyed by solvent name. The six
#: uppercase keys (`MEG`, `DEG`, `TEG`, `MDEA`, `MEA`, `DEA`) are never matched by
#: the lowercased lookup below, so they fall through to 1.0 exactly as NeqSim's own
#: `getAssociationParameter` does.
_WILKE_CHANG_PHI = {
    "water": 2.26,
    "h2o": 2.26,
    "d2o": 2.26,
    "methanol": 1.9,
    "ethanol": 1.5,
    "1-propanol": 1.2,
    "2-propanol": 1.2,
    "1-butanol": 1.0,
    "n-butanol": 1.0,
    "MEG": 1.5,
    "DEG": 1.4,
    "TEG": 1.3,
    "MDEA": 1.5,
    "MEA": 1.7,
    "DEA": 1.5,
    "acetic acid": 1.3,
    "formic acid": 1.6,
}


def wilke_chang_phi(name: str) -> float:
    """The Wilke-Chang association parameter for a solvent, by name.

    A non-associated solvent - most hydrocarbons - is 1.0. The lookup is lowercased,
    mirroring NeqSim, which is why the six uppercase keys in the table never match.
    """
    return _WILKE_CHANG_PHI.get(name.strip().lower(), 1.0)


def entry(name: str, *, card: keycard.Keycard | None = None) -> DatabankEntry:
    """One substance's full record, with any keycard override already applied.

    A keycard wins over the databank, by name, parameter by parameter: a card that
    overrides only ``omega`` keeps the shipped ``Tc`` and ``Pc``. That is the point
    of naming parameters rather than replacing records whole - a user correcting one
    value should not have to restate the others, and should not silently lose them
    if they do not.

    Raises:
        PropertyUnavailableError: if the name is in neither the databank nor the
            keycard, or is in the keycard without every parameter a cubic reads. An
            error rather than a default: substituting a similar substance would
            produce a plausible answer for the wrong fluid, and the caller would
            have no way to see it.
    """
    key = name.strip().lower()
    base = _table().get(key)
    override = card.component(key) if card is not None else None

    if override is None:
        if base is None:
            raise PropertyUnavailableError(
                name,
                "critical constants",
                f"not in the component databank ({len(_table())} substances) and not in "
                f"the loaded keycard. `available()` lists them; `component()` takes a "
                f"Tc, Pc and omega directly for anything else.",
            )
        return base

    if base is None:
        missing = sorted(keycard.CUBIC_PARAMETERS - set(override))
        if missing:
            raise PropertyUnavailableError(
                name,
                "critical constants",
                f"the keycard supplies {sorted(override)} but a cubic needs "
                f"{sorted(keycard.CUBIC_PARAMETERS)}; {missing} is missing. A "
                f"partial component is refused rather than completed from a similar "
                f"substance, which would be inventing data.",
            )
        return DatabankEntry(
            name=name.strip(),
            cas=None,
            formula=None,
            Tc=override["Tc"],
            Pc=override["Pc"],
            omega=_as_float(override["omega"], key),
            molar_mass=None,
            critical_volume=None,
            liquid_density=None,
            # The card may also supply the polynomial, in which case the substance has
            # an enthalpy path; without it, it is a cubic only.
            cp=_cp(override),
            alpha_params={},
            antoine=(0.0, 0.0, 0.0, 0.0, 0.0),
            antoine_type="",
            dipole_moment_debye=0.0,
            viscosity_correction_factor=0.0,
            citation=None,
            source="keycard",
        )

    return replace(
        base,
        Tc=override.get("Tc", base.Tc),
        Pc=override.get("Pc", base.Pc),
        omega=_as_float(override["omega"], key) if "omega" in override else base.omega,
        cp=_cp(override) or base.cp,
        source="keycard",
    )


def _as_float(value: Q, name: str) -> float:
    """A dimensionless quantity as a plain float, refusing a dimensioned one.

    The check is `dimensionality` rather than `check("[dimensionless]")`: pint has no
    `[dimensionless]` dimension, so the obvious spelling raises a `ValueError` about
    the registry rather than returning a verdict about the value.
    """
    if value.dimensionality != ureg.dimensionless.dimensionality:
        raise InvalidInputError(
            f"components.{name}.omega",
            f"the acentric factor is a pure number, but got {value}",
        )
    return float(value.to("dimensionless").magnitude)


def _cp(override: Mapping[str, Q]) -> tuple[float, float, float, float, float] | None:
    """The five heat-capacity coefficients if the card states all of them, else `None`.

    Each is already a quantity in its canonical unit, so the magnitude is the value the
    ideal-gas model takes. A partial polynomial is refused at load, so a card that got
    here names all five or none.
    """
    names = ("cp_a", "cp_b", "cp_c", "cp_d", "cp_e")
    if all(name in override for name in names):
        return (
            float(override["cp_a"].magnitude),
            float(override["cp_b"].magnitude),
            float(override["cp_c"].magnitude),
            float(override["cp_d"].magnitude),
            float(override["cp_e"].magnitude),
        )
    return None


def component(name: str, *, card: keycard.Keycard | None = None) -> Component:
    """One substance as a calculation takes it - `Tc`, `Pc` and `omega`.

    Raises:
        PropertyUnavailableError: as :func:`entry`.
    """
    return entry(name, card=card).component()


def kij_for(
    names: tuple[str, ...], *, card: keycard.Keycard | None = None
) -> dict[tuple[int, int], float]:
    """The interaction pairs the databank knows, for a list of components.

    Only pairs where both names are present are returned, and only where a value
    exists - an unlisted pair is zero, which is the ideal-mixture assumption and is
    what `mixture()` already does with an omitted pair.

    **A keycard's value wins over the databank's**, including when the keycard's
    value is exactly zero: overriding a fitted pair back to ideal mixing is a
    deliberate act, and a rule that treated zero as "absent" would silently undo it.

    Keyed by index into `names`, which is the form `mixture()` takes.
    """
    pairs: dict[tuple[int, int], float] = {}
    for i, a in enumerate(names):
        for j in range(i + 1, len(names)):
            from_keycard = card.kij_for(a, names[j]) if card is not None else None
            value = (
                from_keycard
                if from_keycard is not None
                else _kij().get((a.strip().lower(), names[j].strip().lower()))
            )
            if value is not None and value != 0.0:
                pairs[(i, j)] = value
    return pairs


@dataclass(frozen=True, slots=True)
class NrtlParameters:
    """The NRTL parameters of a mixture, both matrices flattened row-major.

    A caller-supplied record the way :class:`Component` is, and carrying **no names**
    for the same reason: :func:`nrtl_parameters` does the lookup, and a model that
    looked up its own inputs would answer from a file the caller never mentioned.
    """

    #: ``alpha[i][j]``, ``N x N`` row-major. Symmetric with a zero diagonal.
    alpha: tuple[float, ...]
    #: ``dij[i][j] = g_ij`` in kelvin, ``N x N`` row-major. Directional.
    dij: tuple[float, ...]


def nrtl_parameters(names: Sequence[str]) -> NrtlParameters:
    """The NRTL `alpha` and `Dij` matrices for a list of components, by name.

    Both matrices come from the same ``INTER.csv`` row, so ``alpha[i][j]`` and
    ``alpha[j][i]`` are the same number while ``dij[i][j]`` and ``dij[j][i]`` are not.
    ``dij[i][j] = g_ij`` is in Kelvin, and the diagonal of each is zero.

    This is the name-to-matrix resolution `eos.nrtl_activity_coefficients` leaves to
    its caller, the way :func:`bwrs_coefficients` resolves the MBWR-32 set.

    A pair the interaction table does not carry is ``0.0``, which is an ideal
    interaction - what NeqSim's NRTL does with an absent row, and why a pair the table
    has never seen returns ``gamma = 1``. A name that is in *neither* table is refused
    instead, because a typo is a mistake rather than a pair the table happens not to
    cover.

    Raises:
        PropertyUnavailableError: if a name is in neither the databank nor the keycard.
    """
    for name in names:
        entry(name)
    table = _nrtl()
    n = len(names)
    alpha = [0.0] * (n * n)
    dij = [0.0] * (n * n)
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            row = table.get((names[i].strip().lower(), names[j].strip().lower()))
            if row is None:
                continue
            alpha[i * n + j] = row[0]
            dij[i * n + j] = row[1]
    return NrtlParameters(alpha=tuple(alpha), dij=tuple(dij))


@cache
def _unifac_group() -> dict[int, tuple[float, float, int]]:
    """UNIFAC group constants, keyed by subgroup number: `(R, Q, main group)`."""
    out: dict[int, tuple[float, float, int]] = {}
    for row in _rows(find(UNIFAC_GROUP_CSV).read_text(encoding="utf-8")):
        out[int(row["secondary"])] = (
            float(row["volumer"]),
            float(row["surfareaq"]),
            int(row["main"]),
        )
    return out


@cache
def _unifac_aij() -> dict[tuple[int, int], float]:
    """UNIFAC main-group interactions `a_mn` (Kelvin), keyed by main-group pair."""
    out: dict[tuple[int, int], float] = {}
    for row in _rows(find(UNIFAC_INTER_CSV).read_text(encoding="utf-8")):
        m = int(row["maingroup"])
        for key, value in row.items():
            if key.startswith("n") and value.strip():
                out[(m, int(key[1:]))] = float(value)
    return out


@cache
def _unifac_members() -> dict[str, tuple[tuple[int, int], ...]]:
    """UNIFAC group memberships, keyed by name: `(subgroup, count)` tuples."""
    out: dict[str, tuple[tuple[int, int], ...]] = {}
    for row in _rows(find(UNIFAC_COMP_CSV).read_text(encoding="utf-8")):
        name = row["name"].strip().lower()
        subs = tuple(
            (int(key[3:]), int(value))
            for key, value in row.items()
            if key.startswith("sub") and value.strip() not in ("", "0")
        )
        out[name] = subs
    return out


@dataclass(frozen=True, slots=True)
class UnifacParameters:
    """The UNIFAC parameters of a mixture, each matrix flattened row-major.

    A caller-supplied record the way :class:`Component` is, and carrying **no names**
    for the same reason: :func:`unifac_parameters` does the lookup, and a model that
    looked up its own inputs would answer from a file the caller never mentioned.
    """

    #: Per-component group counts, ``N x G`` row-major.
    groups: tuple[float, ...]
    #: The volume ``R`` of each group, length ``G``.
    group_r: tuple[float, ...]
    #: The surface area ``Q`` of each group, length ``G``.
    group_q: tuple[float, ...]
    #: The main-group interaction matrix, ``G x G`` row-major, in kelvin.
    aij: tuple[float, ...]


def unifac_parameters(names: Sequence[str]) -> UnifacParameters:
    """The UNIFAC inputs for a list of components, by name.

    ``groups`` is `N x G` (one row per component, one column per group),
    ``group_r``/``group_q`` are length `G`, and ``aij`` is `G x G` (Kelvin). `G` is
    the union of the named components' subgroups, sorted by subgroup number, with
    absent groups counted zero.

    This is the name-to-matrix resolution `eos.unifac_activity_coefficients` leaves to
    its caller, the way :func:`nrtl_parameters` resolves NRTL's matrices.

    Raises:
        PropertyUnavailableError: if a name has no UNIFAC group assignment.
    """
    group = _unifac_group()
    aij_table = _unifac_aij()
    members = _unifac_members()

    union: list[int] = []
    for name in names:
        subs = members.get(name.strip().lower())
        if subs is None:
            raise PropertyUnavailableError(
                name,
                "UNIFAC group assignment",
                "not in UNIFACcomp.csv; a UNIFAC activity coefficient needs a group "
                "decomposition for every component",
            )
        for subgroup, _ in subs:
            if subgroup not in union:
                union.append(subgroup)
    union.sort()
    g = len(union)

    group_r = [0.0] * g
    group_q = [0.0] * g
    aij = [0.0] * (g * g)
    for k, subgroup in enumerate(union):
        r, q, main = group[subgroup]
        group_r[k] = r
        group_q[k] = q
        for m, other in enumerate(union):
            _, _, other_main = group[other]
            aij[k * g + m] = aij_table.get((main, other_main), 0.0)

    groups = [0.0] * (len(names) * g)
    for i, name in enumerate(names):
        for subgroup, count in members[name.strip().lower()]:
            groups[i * g + union.index(subgroup)] = float(count)

    return UnifacParameters(
        groups=tuple(groups),
        group_r=tuple(group_r),
        group_q=tuple(group_q),
        aij=tuple(aij),
    )


@dataclass(frozen=True, slots=True)
class UniquacParameters:
    """The UNIQUAC volume and surface parameters of a mixture.

    A caller-supplied record the way :class:`Component` is, and carrying **no names**
    for the same reason: :func:`uniquac_parameters` does the lookup.
    """

    #: The van der Waals volume parameter ``r_i`` of each component.
    r: tuple[float, ...]
    #: The van der Waals surface-area parameter ``q_i`` of each component.
    q: tuple[float, ...]


def uniquac_parameters(names: Sequence[str]) -> UniquacParameters:
    """The UNIQUAC `r` and `q` for a list of components, by name.

    ``r_i = sum_k nu_ik R_k`` and ``q_i = sum_k nu_ik Q_k``, the group sums
    :func:`unifac_parameters` forms internally - which is what NeqSim's
    ``ComponentGEUnifac.getR``/``getQ`` compute, and what a UNIQUAC ``r``/``q`` means
    when no fitted value exists.

    **Not** NeqSim's ``rUNIQUAQ``/``qUNIQUAQ`` columns, which ``ComponentGEUniquac``
    reads. Those are ``0.0`` for 109 of the 112 components the table carries - only
    water, acetic acid and ``H2S`` have values - so a UNIQUAC built from them divides
    by zero for almost every real mixture.

    Raises:
        PropertyUnavailableError: if a name has no UNIFAC group decomposition.
    """
    params = unifac_parameters(names)
    g = len(params.group_r)
    n = len(names)
    r = [0.0] * n
    q = [0.0] * n
    for i in range(n):
        for k in range(g):
            r[i] += params.groups[i * g + k] * params.group_r[k]
            q[i] += params.groups[i * g + k] * params.group_q[k]
    return UniquacParameters(r=tuple(r), q=tuple(q))


@dataclass(frozen=True, slots=True)
class VanLaarAcidParameters:
    """The Taleb acid identity of each component of a mixture.

    A caller-supplied record the way :class:`Component` is, and carrying **no names**
    for the same reason: :func:`van_laar_acid_parameters` does the lookup.
    """

    #: Per component: ``1`` water, ``2`` nitric acid, ``3`` sulfuric acid, ``0`` for a
    #: species the model does not cover.
    acid_index: tuple[int, ...]


#: The Taleb (1996) acid index of each spelling, after lower-casing and trimming.
#: NeqSim's ``ComponentGEVanLaarAcid.acidIndexOf`` accepts the formulae as well as the
#: names, and so does this, because a caller writing ``HNO3`` means the same substance
#: as one writing ``nitric acid``.
_ACID_INDEX: dict[str, int] = {
    "water": 1,
    "h2o": 1,
    "nitric acid": 2,
    "hno3": 2,
    "sulfuric acid": 3,
    "sulphuric acid": 3,
    "h2so4": 3,
}


def van_laar_acid_parameters(names: Sequence[str]) -> VanLaarAcidParameters:
    """The Taleb (1996) acid identity of each name, in the components' order.

    ``1`` water, ``2`` nitric acid, ``3`` sulfuric acid, ``0`` for anything else.

    A name the *databank* does not carry is refused. A name it does carry but that is
    not one of the three acids resolves to ``0``, which is not an error:
    :func:`azoth.eos.van_laar_acid_activity_coefficients` gives such a component its
    penalty, which is what NeqSim does with a dissolved carrier gas.

    Raises:
        PropertyUnavailableError: if a name is in neither the databank nor the keycard.
    """
    index: list[int] = []
    for name in names:
        entry(name)
        index.append(_ACID_INDEX.get(name.strip().lower(), 0))
    return VanLaarAcidParameters(acid_index=tuple(index))


def _cubic(name: str) -> Cubic:
    """The cubic named by its short name, ``"pr"``, ``"srk"`` or ``"rk"``."""
    try:
        return CUBICS[name]
    except KeyError:
        raise InvalidInputError(
            "eos", f"unknown cubic {name!r}; expected 'pr', 'srk' or 'rk'"
        ) from None


#: The alpha correlations, keyed by short name. Four feed the one Soave alpha form
#: with a different `m`; `twucoon` and `gassem2001` are non-Soave correlations,
#: `danesh` is Soave's form with `m` scaled by 1.21 above the critical temperature,
#: and the parameterized forms read fitted per-component columns.
_ALPHAS: frozenset[str] = frozenset(
    {
        "pr",
        "srk",
        "pr78",
        "twu",
        "twucoon",
        "gassem2001",
        "danesh",
        "schwartzentruber",
        "mollerup",
        "matcop",
        "matcop_pr",
        "matcop_prumr",
        "matcop_5prumr",
        "delft1998",
    }
)


def _alpha_name(cubic_name: str, alpha: str | None) -> str:
    """The alpha correlation named, defaulting to the cubic's own Soave form."""
    if alpha is None:
        return "srk" if cubic_name == "srk" else "pr"
    if alpha not in _ALPHAS:
        raise InvalidInputError(
            "alpha", f"unknown alpha {alpha!r}; expected one of {sorted(_ALPHAS)}"
        )
    return alpha


def from_names(
    names: list[str],
    *,
    card: keycard.Keycard | None = None,
    eos: str = "pr",
    alpha: str | None = None,
) -> Mixture:
    """A :class:`~azoth.eos.mixture.Mixture` from a list of databank names.

    The interaction parameters come from the databank too, so a caller writing
    `from_names(["methane", "n-butane"])` gets the published `kij` rather than a
    silent zero.

    A mixture alone, for a caller who has their own heat-capacity coefficients or wants
    only to flash. :func:`mixture_of` returns the polynomial beside it, and is what a
    calculation takes - a mixture without one cannot produce an enthalpy.

    Raises:
        PropertyUnavailableError: if any name is not in the databank.
        InvalidInputError: if the list is empty, or a pair is malformed.
    """
    cubic = _cubic(eos)
    alpha_name = _alpha_name(eos, alpha)
    resolved = [name.strip().lower() for name in names]
    components = tuple(entry(name, card=card).component(alpha=alpha_name) for name in resolved)
    return mixture(
        components, kij=kij_for(tuple(resolved), card=card), cubic=cubic, alpha=alpha_name
    )


def mixture_of(
    names: list[str],
    *,
    card: keycard.Keycard | None = None,
    eos: str = "pr",
    alpha: str | None = None,
) -> tuple[Mixture, IdealGasModel]:
    """A mixture and its ideal-gas model, from a list of databank names.

    The two come back together because they are one object in practice: a mixture
    without heat-capacity coefficients cannot produce an enthalpy, and building them
    from two separate lookups invites a call that names different components in each.

    This is the one path a calculation takes to its data. Until it existed, a model's
    spec carried ``Tc``, ``Pc``, ``omega``, ``kij`` and ``cp_a`` through ``cp_e`` as
    nine parallel vectors - numbers in a spec file that nothing could check against
    anything. They are in `data/components/` instead, ported from NeqSim.

    Args:
        names: substance names, matched without regard to case or surrounding space.

    Returns:
        ``(the mixture, the ideal-gas model)``, one entry per name in each vector.

    Raises:
        PropertyUnavailableError: if any name is not in the databank. A name it does
            not have is refused rather than approximated.
        InvalidInputError: if the list is empty, or a substance has no heat-capacity
            coefficients - which is the case for one a keycard added, because a card
            supplies the parameters a *cubic* needs and a polynomial is not one.
    """
    if not names:
        raise InvalidInputError("components", "a mixture needs at least one component")
    resolved = [name.strip().lower() for name in names]
    entries = [entry(name, card=card) for name in resolved]

    missing = [e.name for e in entries if e.cp is None]
    if missing:
        raise InvalidInputError(
            "components",
            f"no heat-capacity coefficients for {missing}. The databank carries them for "
            f"every substance it ships; one a keycard adds needs its own, because a cubic "
            f"needs `Tc`, `Pc` and `omega` and an enthalpy needs the polynomial as well",
        )

    return (
        mixture(
            tuple(e.component(alpha=_alpha_name(eos, alpha)) for e in entries),
            kij=kij_for(tuple(resolved), card=card),
            cubic=_cubic(eos),
            alpha=_alpha_name(eos, alpha),
        ),
        IdealGasModel(
            cp_a=tuple(e.cp[0] for e in entries),  # type: ignore[index]
            cp_b=tuple(e.cp[1] for e in entries),  # type: ignore[index]
            cp_c=tuple(e.cp[2] for e in entries),  # type: ignore[index]
            cp_d=tuple(e.cp[3] for e in entries),  # type: ignore[index]
            cp_e=tuple(e.cp[4] for e in entries),  # type: ignore[index]
        ),
    )


def from_model(name: str, *, card: keycard.Keycard | None = None) -> Mixture:
    """A :class:`~azoth.eos.mixture.Mixture` from a model a keycard declares.

    A keycard's ``models`` section names a cubic variant and the substances it is for
    - named choices from closed vocabularies and nothing to execute, so there is no
    code in a keycard and nothing a keycard can do that this library has not already
    implemented.

    Every declaration is checked against what this build runs, and a name outside the
    vocabularies is refused when the *keycard is loaded* rather than when this is
    called - a model that silently fell back to Peng-Robinson would be a wrong answer
    with no symptom.

    Raises:
        PropertyUnavailableError: if no keycard is loaded, or it declares no model of
            that name.
        PropertyUnavailableError: from :func:`from_names`, if a component of the model
            cannot be resolved.
    """
    model = card.model(name) if card is not None else None
    if model is None:
        where = "the loaded keycard" if card is not None else "no keycard is loaded"
        known = sorted(card.models) if card is not None else []
        raise PropertyUnavailableError(
            name,
            "model definition",
            f"not declared in {where}. Declared models: {known}. A model is a "
            f"keycard's `models` section, not something a calculation resolves on "
            f"its own.",
        )
    return from_names(list(model.components), card=card)


@dataclass(frozen=True, slots=True)
class BwrsCoefficients:
    """The MBWR-32 coefficients of one substance.

    The 32 fitted coefficients ``a0``..``a31`` and the critical density ``rhoc`` in
    mol/L, verbatim from NeqSim's ``MBWR32param`` table. Only methane and ethane have
    them.
    """

    #: ``a0``..``a31``, in the native mol/L, MPa convention.
    a: tuple[float, ...]
    #: The critical density, in mol/L.
    rhoc: float


@cache
def _mbwr32_table() -> dict[str, BwrsCoefficients]:
    rows = _rows(find(MBWR32_CSV).read_text(encoding="utf-8"))
    out: dict[str, BwrsCoefficients] = {}
    for row in rows:
        name = row["name"].strip().lower()
        out[name] = BwrsCoefficients(
            a=tuple(float(row[f"a{i}"]) for i in range(32)),
            rhoc=float(row["rhoc"]),
        )
    return out


def bwrs_coefficients(name: str) -> BwrsCoefficients:
    """The MBWR-32 coefficients of a substance, by name.

    Raises:
        PropertyUnavailableError: if the name is not in the MBWR-32 table. Only
            methane and ethane have parameters; anything else is refused rather than
            given an estimated density.
    """
    table = _mbwr32_table()
    key = name.strip().lower()
    try:
        return table[key]
    except KeyError:
        raise PropertyUnavailableError(
            name,
            "MBWR-32 coefficients",
            "only methane and ethane have MBWR-32 parameters in NeqSim's mbwr32param; "
            "anything else would need an estimated critical density",
        ) from None


# Imported at the bottom because `mixture` lives with the types this module builds
# on, and importing it at the top would make the cycle explicit for no gain.
from azoth.eos.mixture import mixture  # noqa: E402

__all__ = [
    "BwrsCoefficients",
    "DatabankEntry",
    "NrtlParameters",
    "UnifacParameters",
    "UniquacParameters",
    "VanLaarAcidParameters",
    "available",
    "bwrs_coefficients",
    "component",
    "entry",
    "from_model",
    "from_names",
    "kij_for",
    "nrtl_parameters",
    "unifac_parameters",
    "uniquac_parameters",
    "van_laar_acid_parameters",
    "wilke_chang_phi",
]
