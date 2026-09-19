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

286 substances: NeqSim's `COMP.csv` filtered to the types a cubic equation of state
can describe. The ions are excluded because a cubic has no notion of one and NeqSim
fills their critical properties with a shared default - the 62 ion rows carry only 26
distinct critical sets between them, the commonest shared by 27, which is not a
coincidence, and shipping them would ship plausible-looking wrong numbers.
`COMP_EXT.csv` is not vendored at all: 86 MB of heavy fluids this library cannot
characterise.

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
for 286 rows would be a citation-shaped thing that is not a citation.
"""

from __future__ import annotations

import csv
import io
from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass, replace
from functools import cache
from typing import NamedTuple

from azoth import keycard
from azoth._data import find
from azoth.core.errors import InvalidInputError, PropertyUnavailableError
from azoth.core.units import Q, ureg
from azoth.eos.cubic import CUBICS, Cubic
from azoth.eos.mixture import Component, Mixture
from azoth.eos.reference._association import family_of
from azoth.eos.reference.molar_enthalpy_entropy import IdealGasModel

COMPONENTS_CSV = "data/components/components.csv"
KIJ_CSV = "data/components/kij.csv"
UNIFAC_COMP_CSV = "data/components/UNIFACcomp.csv"
UNIFAC_GROUP_CSV = "data/components/UNIFACGroupParam.csv"
UNIFAC_INTER_CSV = "data/components/UNIFACInterParam.csv"
UNIFAC_INTER_B_CSV = "data/components/UNIFACInterParamB.csv"
UNIFAC_INTER_C_CSV = "data/components/UNIFACInterParamC.csv"
UNIFAC_COMP_UMRPRU_CSV = "data/components/UNIFACcompUMRPRU.csv"
MBWR32_CSV = "data/components/mbwr32.csv"
PITZER_PARAMETERS_CSV = "data/components/PitzerParameters.csv"
COMPSALT_CSV = "data/components/COMPSALT.csv"

#: The `REFERENCESTATETYPE` value NeqSim treats as the symmetric (Raoult) reference.
#:
#: `ComponentGE.fugcoef` compares the component's column against exactly this string and
#: falls to the Henry's-law branch for anything else - including the literal `0.0` two of
#: the table's rows carry.
SOLVENT = "solvent"

#: The `COMPTYPE` value whose rows carry no critical constants of their own.
#:
#: `COMP.csv` files 62 rows under it, and for 27 of them the four critical columns hold one
#: default set - `Tc = 444.15` as the file states it, which the generator's 273.15 makes
#: `717.3 K`, then `Pc = 290.89 bar`, `omega = 0.344` and `Vc = 99 cm3/mol`. The remaining
#: 35 inherit a neutral parent's numbers (`MDEA+` carries MDEA's) or carry values with no
#: stated source. `DatabankEntry.component_type` is where the tag reaches this side, and
#: :func:`_refuse_ions` refuses rather than reading the filler.
ION = "ion"

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
    "umrcpa_a0",
    "umrcpa_b",
    "umrcpa_assocenergy",
    "umrcpa_assocvolume",
    "umrcpa_associating",
    "umrcpa_racketz",
    "antoine_type",
    "antoinea",
    "antoineb",
    "antoinec",
    "antoined",
    "antoinee",
    "dipole_moment_debye",
    "viscosity_correction_factor",
    "referencestatetype",
    "associationscheme",
    "associationsites",
    "associationenergy",
    "associationboundingvolume_srk",
    "associationboundingvolume_pr",
    "acpa_srk",
    "bcpa_srk",
    "mcpa_srk",
    "acpa_pr",
    "bcpa_pr",
    "mcpa_pr",
    "racketzcpa",
    "volcorrcpa_t",
    "msaft",
    "sigma_saft_m",
    "epsiksaft",
    "citation",
)


#: The site schemes the databank names, each mapped to the number of sites the scheme
#: has. Upstream also writes ``0``, its marker for a component with no scheme at all.
#:
#: A site *count* cannot stand in for this mapping: the table carries ``1A`` on a
#: component with zero sites and both ``2A`` and ``2B`` on components with two, and
#: NeqSim's ``setAssociationScheme`` switches on the name.
SITE_SCHEMES: dict[str, int] = {"1A": 1, "2A": 2, "2B": 2, "4C": 4}

#: The schemes that cannot self-associate, because every site carries the same charge
#: sign and NeqSim's bond test is a product sign. A component with one of these still
#: *cross*-associates with an oppositely-charged partner.
SELF_BONDLESS_SCHEMES = frozenset({"1A", "2A"})

#: The factor from NeqSim's internal scale for a fitted attraction and covolume to SI,
#: `Component.java:526-531`. **The one crossing between a keycard and the table**: a card
#: states the published SI value, the compiled table carries the internal scale, and
#: :func:`_card_association` is where the two meet.
NEQSIM_INTERNAL_TO_SI = 1.0e-5


@dataclass(frozen=True, slots=True)
class UmrCpaParameters:
    """The ``UMRCPA_*`` fitted set, in the scale the table states it in.

    NeqSim's own conversions, from ``Component.java:526-539``: ``UMRCPA_a0`` is in
    ``bar L**2/mol**2`` and becomes internal ``a`` by ``x 1e4``, ``UMRCPA_b`` in ``L/mol``
    by ``x 1e2``, and ``UMRCPA_assocEnergy`` in ``bar L/mol`` - which is 100 J/mol, so the
    ``x 1e2`` there is a unit conversion and not a scale. ``UMRCPA_assocVolume`` is the
    dimensionless ``kappa_AB`` and is used as it stands.
    """

    #: ``UMRCPA_a0``, in ``bar L**2/mol**2``.
    a0: float
    #: ``UMRCPA_b``, in ``L/mol``.
    b: float
    #: ``UMRCPA_assocEnergy``, in ``bar L/mol``.
    energy: float
    #: ``UMRCPA_assocVolume``, the dimensionless ``kappa_AB``.
    volume: float
    #: ``UMRCPA_racketZ``, the Rackett compressibility of the volume correction.
    racket_z: float

    @property
    def attraction(self) -> float:
        """The fitted attraction in SI, ``Pa*m**6/mol**2``."""
        return self.a0 * 1.0e4 * NEQSIM_INTERNAL_TO_SI

    @property
    def covolume(self) -> float:
        """The fitted covolume in SI, ``m**3/mol``."""
        return self.b * 1.0e2 * NEQSIM_INTERNAL_TO_SI

    @property
    def energy_j_per_mol(self) -> float:
        """The association energy in J/mol."""
        return self.energy * 1.0e2


@dataclass(frozen=True, slots=True)
class AssociationParameters:
    """One component's association parameters, as the databank records them.

    The two fitted cubic sets are in NeqSim's internal scale - ``a`` in
    ``Pa*m**6/mol**2 * 1e5`` and ``b`` in ``m**3/mol * 1e5`` - because the table carries
    what the source table carries. The conversion belongs with the model that reads
    them. Water is the check: :attr:`b_srk` is 1.4515, i.e. ``1.4515e-5 m**3/mol``,
    against the cubic's own ``2.11e-5``.

    **Those fitted values are why a CPA mixture is not determined by ``Tc`` and
    ``Pc``.** A component that associates carries its own attraction and covolume, and
    a model that quietly used the cubic's would be wrong by tens of percent.
    """

    #: The scheme's name, one of :data:`SITE_SCHEMES`.
    scheme: str
    #: The site count the table states. Carried beside the scheme because the two
    #: disagree upstream, and a silent resolution would hide that.
    sites: int
    #: The association energy ``eps``, in J/mol.
    energy: float
    #: ``kappa_AB`` for the SRK family.
    volume_srk: float
    #: The fitted attraction for SRK-CPA, in NeqSim's internal scale.
    a_srk: float
    #: The fitted covolume for SRK-CPA, in NeqSim's internal scale.
    b_srk: float
    #: The SRK alpha correlation's ``m``.
    m_srk: float
    #: ``kappa_AB`` for the PR family.
    volume_pr: float
    #: The fitted attraction for PR-CPA, in NeqSim's internal scale.
    a_pr: float
    #: The fitted covolume for PR-CPA, in NeqSim's internal scale.
    b_pr: float
    #: The PR alpha correlation's ``m``.
    m_pr: float
    #: ``racketZCPA``, the Rackett compressibility NeqSim's CPA volume correction reads.
    racket_z: float
    #: ``volcorrCPA_T``, the CPA volume-translation coefficient.
    #:
    #: NeqSim's ``SystemSrkCPA`` calls ``useVolumeCorrection(true)``, so its root carries
    #: a translation this library does not yet apply. It is **not** the cause of the CPA
    #: root gap - for water the shift is 2.5e-05 internal, two orders of magnitude short.
    volume_correction: float
    #: The ``UMRCPA_*`` set, which replaces the families' fitted values where it is present.
    #:
    #: A **third** fitted set beside the SRK and PR families, not a cubic family of its
    #: own: NeqSim's ``Component.java:534-539`` overrides ``aCPA``/``bCPA``/
    #: ``associationVolume``/``associationEnergy`` from these columns for a
    #: ``ComponentUMRCPA`` whose row carries them.
    umr_cpa: UmrCpaParameters | None = None

    @property
    def self_bonds(self) -> bool:
        """Whether this component's own sites bond with each other."""
        return self.scheme not in SELF_BONDLESS_SCHEMES

    def has_fitted_set(self, family: str) -> bool:
        """Whether this record carries a usable fitted set at one cubic family.

        **NeqSim's own guard, and load-bearing.** ``ComponentSrkCPA`` substitutes the fitted
        values only ``if (Math.abs(aCPA) > 1e-6)``, in its internal units - and the table has
        rows where that matters in both directions: CO2 names the ``2A`` scheme with every
        fitted value zero, and H2S and benzene carry an SRK set and a PR set that is entirely
        zero. Substituting unconditionally gives those a covolume of zero, which makes the
        reduced pressure ``NaN`` rather than a different answer.
        """
        return abs(self.a_srk if family == "srk" else self.a_pr) > 1.0e-6

    def attraction(self, family: str) -> float:
        """The fitted attraction ``a0`` at one family, in SI - ``Pa*m**6/mol**2``."""
        internal = self.a_srk if family == "srk" else self.a_pr
        return internal * NEQSIM_INTERNAL_TO_SI

    def covolume(self, family: str) -> float:
        """The fitted covolume at one family, in SI - ``m**3/mol``.

        **Not the cubic's own.** NeqSim's ``ComponentSrkCPA`` substitutes this for the ``b``
        a cubic derives from ``Tc`` and ``Pc``, and the difference is large: water's is
        ``1.4515e-5`` against ``0.08664 R Tc/Pc``'s ``2.11e-5``.
        """
        internal = self.b_srk if family == "srk" else self.b_pr
        return internal * NEQSIM_INTERNAL_TO_SI

    def alpha_m(self, family: str) -> float:
        """The fitted Soave alpha coefficient ``m`` at one family.

        NeqSim's ``ComponentSrkCPA`` calls ``getAttractiveTerm().setm(mCPA)``, so the alpha
        function is Soave's form with a *fitted* coefficient rather than the one
        ``0.480 + 1.574 omega - 0.176 omega**2`` would give.
        """
        return self.m_srk if family == "srk" else self.m_pr


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
    #: NeqSim's `REFERENCESTATETYPE`, which decides the branch
    #: ``ComponentGE.fugcoef`` takes: :data:`SOLVENT` gives ``gamma_i P0_i / P`` and
    #: anything else a Henry's-law coefficient. A keycard-supplied substance is named
    #: ``solvent``, because a card states what a cubic reads and a cubic has no
    #: reference state.
    reference_state: str
    #: The substance class, exactly as NeqSim's ``COMPTYPE`` column states it. Read for
    #: one thing, and it is a refusal rather than a term: no cubic is built over an
    #: ``"ion"``, whose critical columns hold a default rather than a measurement.
    #: Never computed with.
    component_type: str
    #: The association parameters, or ``None`` for a component the table gives no
    #: scheme - which is 181 of the 348 compiled rows, and is the table's own marker
    #: rather than a missing value.
    association: AssociationParameters | None
    #: PC-SAFT's segment number ``m``, dimensionless. **Zero means the table carries no
    #: PC-SAFT set**, which is how 106 of the 348 rows spell absence - the column is never
    #: blank - so a model must refuse a zero rather than solve for a fluid with no
    #: segments. The rule ``ComponentSrkCPA``'s ``|aCPA| > 1e-6`` guard already follows.
    m_saft: float
    #: PC-SAFT's segment diameter ``sigma``, in **metres**. NeqSim's `COMP.csv` carries it
    #: in ångström and divides by ``1e10`` as it reads it (``Component.java:548``); this
    #: table is compiled SI, so the conversion happened once in the generator.
    sigma_saft: float
    #: PC-SAFT's segment energy over Boltzmann's constant, ``epsilon/k``, in K.
    epsik_saft: float
    #: SAFT-VR-Mie's repulsive exponent ``lambda_r``, dimensionless. **It is an absence
    #: marker too, on exactly the rows :attr:`m_mie` is**: the table carries the standard
    #: ``12`` on 336 of its 348 rows and a fitted exponent on the 12 with a set, so
    #: ``lambda_r_mie == 12.0`` and ``m_mie == 0.0`` hold on the same rows.
    lambda_r_mie: float
    #: SAFT-VR-Mie's attractive exponent ``lambda_a``, dimensionless, ``6`` on 347 rows.
    #: Not a marker: the fitted rows carry ``6`` as well, apart from ``co2`` at ``5.055``.
    lambda_a_mie: float
    #: SAFT-VR-Mie's segment number ``m``, dimensionless. **Zero means the table carries
    #: no SAFT-VR-Mie set**, which is how 336 of the 348 rows spell absence - 12 carry
    #: one, the light alkanes with ``co2``, nitrogen and water. ``sigma_mie`` and
    #: ``epsik_mie`` are zero on exactly those rows, so a model refuses rather than
    #: solving for a fluid with no segments.
    m_mie: float
    #: SAFT-VR-Mie's segment diameter ``sigma``, in metres.
    sigma_mie: float
    #: SAFT-VR-Mie's segment energy over Boltzmann's constant, ``epsilon/k``, in K.
    epsik_mie: float
    #: The ionic charge, in units of the elementary charge; zero for a neutral. **Zero
    #: does not mean "not an ion"**: four rows typed ``"ion"`` carry it - ``h+pzcoo-`` is
    #: a zwitterion and ``caco3``, ``nacl`` and ``cacl2`` are neutral salts filed with
    #: the ions. :attr:`component_type` is the marker; this is the value.
    ionic_charge: float
    #: The Deshmukh-Mather ion diameter, in ångström. **The table's absence marker is
    #: ``0.0``, not a blank**, which is what 18 of the 62 ``"ion"`` rows carry. A model
    #: that needs one must refuse a zero: a zero diameter is the infinite-dilution limit
    #: of the term it belongs to, so a plausible number would come out of it.
    deshmukh_mather_diameter: float
    #: ``DIELECTRICPARAMETER1..5``, the dielectric mixing rule's coefficients. Carried
    #: exactly as NeqSim stores them. The unit is not established - the manifest records
    #: ``neqsim-internal`` - so no conversion is applied and none is guessed.
    dielectric: tuple[float, float, float, float, float]
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
            # Carried whether or not the mixture runs it: an associating equation of state
            # reads it and a cubic ignores it, and whether a *phase model* associates is
            # the model's decision rather than the substance's - `SystemNRTL` builds a
            # classical phase over the same methanol and water that `SystemSrkCPA` builds
            # an associating one over.
            association=self.association,
        )

    def antoine_form(self) -> str:
        """The cleaned Antoine form this entry's coefficients belong to, one of
        ``"dippr101"``, ``"pow10"``, ``"pow10kpa"``, ``"exp"`` or ``"wagner"``.

        **The empty string means the row has no correlation at all** - its label is
        ``none``, upstream's marker for unavailable data - and a caller must refuse rather
        than evaluate. See :func:`form_from_type`.
        """
        return form_from_type(self.antoine_type, self.antoine[4])

    def __repr__(self) -> str:
        return f"DatabankEntry({self.name!r}, Tc={self.Tc}, Pc={self.Pc}, omega={self.omega})"


def _rows(text: str) -> list[dict[str, str]]:
    body = [
        line for line in text.splitlines() if line.strip() and not line.lstrip().startswith("#")
    ]
    return [dict(row) for row in csv.DictReader(io.StringIO("\n".join(body)))]


def _refuse_ions(names: list[str]) -> None:
    """Refuse a cubic over an ion, naming every one.

    An ion's critical columns in NeqSim's table are a default rather than a measurement,
    so a cubic built from them returns plausible-looking wrong numbers - which is the
    failure this module is arranged against. An ion belongs to the electrolyte models,
    which read its charge and diameter instead.
    """
    if not names:
        return
    verb = "is" if len(names) == 1 else "are"
    raise InvalidInputError(
        "components",
        f"{', '.join(names)} {verb} an ion, and a cubic has no notion of one. NeqSim's table "
        f"fills its ion rows' critical columns with a default - one shared set on 27 of the "
        f"62, a neutral parent's numbers on most of the rest - so a cubic built from them "
        f"returns plausible-looking wrong numbers, and a card that adds an ion has no better "
        f"constants to give. An ion belongs to the electrolyte models, which read its charge "
        f"and diameter instead.",
    )


def _absent_is_zero(text: str) -> float:
    """One parameter, where an empty cell means zero.

    **Both upstream tables use a blank and a ``0`` interchangeably for "no parameter"**,
    so a blank infers nothing the row does not already state. Four cells in the compiled
    tables are blank: ``associationboundingvolume_pr`` on ``h2so4``, ``hno3`` and
    ``c2h4-``, and ``nrtlgij`` on the ``k+``/``h2s`` pair, which carries ``KIJPR = 1``
    with no NRTL energy beside it. A *malformed* value still raises - the distinction is
    between an absent parameter and a broken one, which ``float`` keeps.
    """
    return float(text) if text.strip() else 0.0


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
            reference_state=row["referencestatetype"].strip(),
            component_type=row["comptype"].strip().lower(),
            association=_association(row),
            m_saft=float(row["msaft"]),
            sigma_saft=float(row["sigma_saft_m"]),
            epsik_saft=float(row["epsiksaft"]),
            lambda_r_mie=float(row["lambdarsaftvrmie"]),
            lambda_a_mie=float(row["lambdaasaftvrmie"]),
            m_mie=float(row["msaftvrmie"]),
            sigma_mie=float(row["sigma_saft_vr_mie_m"]),
            epsik_mie=float(row["epsiksaftvrmie"]),
            ionic_charge=float(row["ioniccharge"]),
            deshmukh_mather_diameter=float(row["deshmationicdiameter"]),
            dielectric=(
                float(row["dielectricparameter1"]),
                float(row["dielectricparameter2"]),
                float(row["dielectricparameter3"]),
                float(row["dielectricparameter4"]),
                float(row["dielectricparameter5"]),
            ),
            citation=row["citation"],
        )
    return entries


def _three(row: Mapping[str, str], a: str, b: str, c: str) -> tuple[float, float, float]:
    return (float(row[a]), float(row[b]), float(row[c]))


def _association(row: Mapping[str, str]) -> AssociationParameters | None:
    """One row's association parameters, or ``None`` where it carries no scheme.

    **An empty cell means zero here.** The upstream table uses a blank and a ``0``
    interchangeably for "no parameter": the three blanks in the compiled table are
    ``associationboundingvolume_pr`` on ``h2so4``, ``hno3`` and ``c2h4-``, whose every
    sibling association cell is already ``0.0``. See :func:`_absent_is_zero` for the
    fourth. So a blank infers nothing the row does not already state.
    """
    scheme = row["associationscheme"].strip()
    if scheme not in SITE_SCHEMES:
        return None

    def value(column: str) -> float:
        raw = row[column].strip()
        return 0.0 if not raw else float(raw)

    return AssociationParameters(
        scheme=scheme,
        sites=int(value("associationsites")),
        energy=value("associationenergy"),
        volume_srk=value("associationboundingvolume_srk"),
        a_srk=value("acpa_srk"),
        b_srk=value("bcpa_srk"),
        m_srk=value("mcpa_srk"),
        volume_pr=value("associationboundingvolume_pr"),
        a_pr=value("acpa_pr"),
        b_pr=value("bcpa_pr"),
        m_pr=value("mcpa_pr"),
        racket_z=value("racketzcpa"),
        volume_correction=value("volcorrcpa_t"),
        umr_cpa=_parse_umr_cpa(value),
    )


def _parse_umr_cpa(
    value: Callable[[str], float],
) -> UmrCpaParameters | None:
    """One row's ``UMRCPA_*`` set, or ``None`` where the row does not carry one.

    NeqSim's guard is two-part, ``Component.java:534``: ``UMRCPA_associating`` is one
    **and** ``|UMRCPA_a0| > 1e-20``. Seven of the 286 rows pass it - water, methanol,
    ethanol and the four glycols - and a row that fails it keeps the PR family's fitted
    values, which is what ``ComponentUMRCPA`` reads first for being a ``ComponentPR``.
    """
    if abs(value("umrcpa_associating") - 1.0) > 0.0 or abs(value("umrcpa_a0")) <= 1.0e-20:
        return None
    return UmrCpaParameters(
        a0=value("umrcpa_a0"),
        b=value("umrcpa_b"),
        energy=value("umrcpa_assocenergy"),
        volume=value("umrcpa_assocvolume"),
        racket_z=value("umrcpa_racketz"),
    )


def _card_association(
    base: AssociationParameters | None, stated: keycard.Association
) -> AssociationParameters:
    """A card's association, over the table's, parameter by parameter.

    The rule a card follows for the cubic's parameters, applied to the association's: a
    parameter the card names is the card's, and one it does not is the table's. A
    substance the table gives no scheme therefore gets the card's set and zeros beside it,
    which is NeqSim's own ``|aCPA| > 1e-6`` guard reached from the other side - a cubic
    whose attraction no fitted value replaces.

    **The site count is the table's where the scheme is, and the scheme's where the scheme
    is not.** A record carries a count *beside* its scheme because the two disagree
    upstream - the table's ``1A`` rows carry a count of zero - and a card that restated the
    scheme would otherwise resolve that disagreement silently, which is what the field
    exists to keep visible. A card naming a *different* scheme is stating a different
    molecule and has no count to inherit.
    """
    stated_values = stated.parameters

    def number(name: str, table: float) -> float:
        quantity = stated_values.get(name)
        return table if quantity is None else float(quantity.magnitude)

    def fitted(name: str, table: float) -> float:
        quantity = stated_values.get(name)
        return table if quantity is None else float(quantity.magnitude) / NEQSIM_INTERNAL_TO_SI

    if base is not None and base.scheme == stated.scheme:
        sites = base.sites
    else:
        sites = SITE_SCHEMES[stated.scheme]

    return AssociationParameters(
        scheme=stated.scheme,
        sites=sites,
        energy=number("energy", 0.0 if base is None else base.energy),
        volume_srk=number("volume_srk", 0.0 if base is None else base.volume_srk),
        a_srk=fitted("a_srk", 0.0 if base is None else base.a_srk),
        b_srk=fitted("b_srk", 0.0 if base is None else base.b_srk),
        m_srk=number("m_srk", 0.0 if base is None else base.m_srk),
        volume_pr=number("volume_pr", 0.0 if base is None else base.volume_pr),
        a_pr=fitted("a_pr", 0.0 if base is None else base.a_pr),
        b_pr=fitted("b_pr", 0.0 if base is None else base.b_pr),
        m_pr=number("m_pr", 0.0 if base is None else base.m_pr),
        racket_z=0.0 if base is None else base.racket_z,
        volume_correction=0.0 if base is None else base.volume_correction,
        # A card states the parameters a *cubic* reads and carries no `UMRCPA_*` set, so
        # the table's survives a card untouched - the whole point of naming one parameter
        # rather than restating a record.
        umr_cpa=None if base is None else base.umr_cpa,
    )


def form_from_type(label: str, e: float) -> str:
    """Map NeqSim's raw ``AntoineVapPresLiqType`` label onto the clean form name.

    **``e`` is part of the question, not only the label.** Twenty rows in NeqSim's
    ``COMP.csv`` carry DIPPR-101 coefficients - ``exp(A + B/T + C ln T + D T**E)`` -
    under the label ``log``, which names the two-term exponential instead. Reading the
    label alone therefore picks the wrong correlation for them, by thirty to ninety
    orders of magnitude; ``i-pentane`` comes out at 8.3e38 bar.

    The rule is NeqSim's own, ``Component.usesDipprVaporPressureCorrelation``: a
    non-zero ``e`` decides, except that ``pow10`` and ``pow10KPa`` keep precedence
    because those coefficients are log10-based and would not survive the exponential
    form.

    ``exp`` and ``log`` are otherwise one formula under two names, and
    ``loglog``/``log10`` have no branch in NeqSim's dispatch, so they fall through to
    Wagner - a defect this reproduces rather than silently repairs.

    **``none`` is not a form, and it is empty here rather than a fall-through.** It is
    the marker upstream added in ``83b64e5`` (PR #3775) for a row whose correlation is
    *unavailable*: 313 of its 389 rows now carry it, with all five coefficients zero.
    Sent to Wagner those give ``exp(0) * Pc = Pc``, so an unavailable correlation would
    come back as the component's **critical pressure** - measured, ``1.82e6 Pa`` for
    ``nc12``, whose vapour pressure there is about ``42 Pa``. A caller that gets ``""``
    back must refuse, as :func:`azoth.eos.reference._ge_phase.saturation` does.
    """
    if label == "none":
        return ""
    if abs(e) > 1e-12 and label not in ("pow10", "pow10KPa"):
        return "dippr101"
    if label == "pow10":
        return "pow10"
    if label == "pow10KPa":
        return "pow10kpa"
    if label in ("exp", "log"):
        return "exp"
    return "wagner"


@cache
def _kij() -> dict[tuple[str, str], tuple[float, float]]:
    """Interaction parameters, keyed by the ordered pair of names.

    `(srk, pr)`, because **the cubic chooses its column and the two are different
    fits**. NeqSim selects on the phase class - `PhasePrEos` reads `KIJPR` and every
    other phase reads `KIJSRK` - and 194 of the 957 in-scope pairs differ: benzene/methane
    is `0` in one and `0.0209` in the other.
    """
    pairs: dict[tuple[str, str], tuple[float, float]] = {}
    for row in _rows(find(KIJ_CSV).read_text(encoding="utf-8")):
        pair = (row["component_a"], row["component_b"])
        columns = (_absent_is_zero(row["kijsrk"]), _absent_is_zero(row["kij_pr"]))
        pairs[pair] = columns
        pairs[(pair[1], pair[0])] = columns
    return pairs


@cache
def _cpa_kij() -> dict[tuple[str, str], tuple[float, float]]:
    """The associating interaction columns, keyed by the ordered pair.

    `(srk, pr)`, which are two fits rather than one converted. A separate reader from
    :func:`_kij` because it is a separate column of the same file: NeqSim's
    `CPAMixingRuleHandler` reads these and a classical rule reads `KIJSRK` or
    `KIJPR` depending on its cubic, and on water/methanol they differ by a factor of two.
    """
    pairs: dict[tuple[str, str], tuple[float, float]] = {}
    for row in _rows(find(KIJ_CSV).read_text(encoding="utf-8")):
        a, b = row["component_a"], row["component_b"]
        columns = (_absent_is_zero(row["cpakij_srk"]), _absent_is_zero(row["cpakij_pr"]))
        pairs[(a, b)] = columns
        pairs[(b, a)] = columns
    return pairs


@cache
@cache
def pcsaft_kij_for(names: tuple[str, ...]) -> dict[tuple[int, int], float]:
    """The *PC-SAFT* interaction pairs for a list of components, by index.

    The sibling of :func:`kij_for` and :func:`cpa_kij_for` for the third column, with the
    same rules: an unlisted pair is zero, which is the ideal-mixture default NeqSim
    substitutes.

    **`KIJPCSAFT` is a fit of its own.** Methane/n-butane is ``0.022`` here against the
    ``0.01289789`` its SRK and PR columns share, and 42 of the 957 in-scope pairs carry one
    at all - so this is a sparse column rather than a convention for one of the others.

    **NeqSim reads it in the only configuration that runs.** Its default mixing rule leaves
    the matrix null and a *mixture* throws a NullPointerException before reaching a ``k_ij``
    - a pure fluid never asks for a pair - so a caller must set the classic rule, and the
    branch that serves it reads this column. There is no usable NeqSim PC-SAFT that ignores
    it, which is why no card states one yet: a keycard's pair record carries the cubic's
    ``kij`` and the two CPA columns.
    """
    stored = _pcsaft_kij()
    pairs: dict[tuple[int, int], float] = {}
    for i, a in enumerate(names):
        for j in range(i + 1, len(names)):
            value = stored.get((a.strip().lower(), names[j].strip().lower()))
            if value is not None and value != 0.0:
                pairs[(i, j)] = value
    return pairs


@cache
def _pcsaft_kij() -> dict[tuple[str, str], float]:
    """`KIJPCSAFT`, keyed by the ordered pair and stored both ways round.

    Separate from :func:`_kij` and :func:`_cpa_kij` because it is a third column of the
    same file rather than a third convention for one of them.
    """
    pairs: dict[tuple[str, str], float] = {}
    for row in _rows(find(KIJ_CSV).read_text(encoding="utf-8")):
        a, b = row["component_a"], row["component_b"]
        value = _absent_is_zero(row["kijpcsaft"])
        pairs[(a, b)] = value
        pairs[(b, a)] = value
    return pairs


def _nrtl() -> dict[tuple[str, str], tuple[float, float]]:
    """NRTL `(alpha, gij)`, keyed by ordered pair and stored both ways round.

    `gij` is the Kelvin energy for the ordered key `(first, second)`, so the reversed
    key carries the reversed energy; `alpha` is symmetric and appears under both keys.
    """
    pairs: dict[tuple[str, str], tuple[float, float]] = {}
    for row in _rows(find(KIJ_CSV).read_text(encoding="utf-8")):
        a = row["component_a"]
        b = row["component_b"]
        alpha = _absent_is_zero(row["nrtlalpha"])
        pairs[(a, b)] = (alpha, _absent_is_zero(row["nrtlgij"]))
        pairs[(b, a)] = (alpha, _absent_is_zero(row["nrtlgji"]))
    return pairs


def _column(row: Mapping[str, str], name: str) -> float:
    """One column of a row, where a blank is the table's own absence marker."""
    return _absent_is_zero(row[name])


@dataclass(frozen=True, slots=True)
class PitzerRecord:
    """One ion pair's Pitzer parameters, as `PitzerParameters.csv` states them.

    The three binary coefficients are **temperature-dependent**, and the form is NeqSim's
    own (``PhasePitzer``, the comment above its ``beta0T1`` field):

    .. code-block:: text

        beta0(T) = beta0_25 + t1 (1/T - 1/Tr) + t2 ln(T/Tr),   Tr = 298.15 K

    and likewise for ``beta1`` and ``Cphi``. The ``*_25`` names are the values at ``Tr``
    and the ``*_t`` pairs are ``(t1, t2)``, so a caller evaluates the correlation rather
    than reading a number that is only right at one temperature.
    """

    ion1: str
    ion2: str
    beta0_25: float
    beta1_25: float
    cphi_25: float
    beta0_t: tuple[float, float]
    beta1_t: tuple[float, float]
    cphi_t: tuple[float, float]
    #: ``beta2``, nonzero on exactly four rows - ``ca++``, ``mg++``, ``sr++`` and
    #: ``fe++`` against ``so4--``, all negative. So it is a parameter most pairs do not
    #: have rather than a column a model can read unconditionally.
    beta2_25: float
    #: ``theta``, the same-sign ion mixing parameter (Harvie and Weare, 1984). **Zero on
    #: every row of the shipped table**: all 30 rows are a cation against an anion, so
    #: the same-sign term has no parameters here and the zero is an absence rather than a
    #: fitted ideal solution.
    theta: float
    #: ``psi``, the ternary parameter, stated against the pair whose interaction it
    #: modifies. **Zero on every row too**, for :attr:`theta`'s reason: a ternary
    #: parameter needs a same-sign pair to hang off and the table carries none.
    psi_common_ion: float
    #: The pair's fitted validity range, in K. **A caller outside it is extrapolating**.
    t_min: float
    t_max: float
    reference: str


@dataclass(frozen=True, slots=True)
class SaltRecord:
    """One salt's record, as `COMPSALT.csv` states it."""

    name: str
    cation: str
    anion: str
    cation_stoichiometry: float
    anion_stoichiometry: float
    #: The five solubility-product coefficients. **Their form is not established here** -
    #: the manifest records the columns as ``neqsim-internal`` - so they are carried as
    #: the file states them and a model that needs the correlation must read NeqSim's
    #: ``ChemicalReactionOperations`` rather than assume a polynomial.
    ksp: tuple[float, float, float, float, float]
    #: ``Vdelta``, the molar volume change on dissolution. Unit unestablished, as above.
    volume_delta: float
    #: ``waterstoc``, how many waters of hydration the dissolution carries. Zero on every
    #: row of the shipped table.
    water_stoichiometry: float


@cache
def pitzer_parameters() -> tuple[PitzerRecord, ...]:
    """Every row of the Pitzer pair table, in file order."""
    out = []
    for row in _rows(find(PITZER_PARAMETERS_CSV).read_text(encoding="utf-8")):
        ion1 = row["ion1"].strip().lower()
        if not ion1:
            continue

        out.append(
            PitzerRecord(
                ion1=ion1,
                ion2=row["ion2"].strip().lower(),
                beta0_25=_column(row, "beta0_25"),
                beta1_25=_column(row, "beta1_25"),
                cphi_25=_column(row, "cphi_25"),
                beta0_t=(_column(row, "beta0_t1"), _column(row, "beta0_t2")),
                beta1_t=(_column(row, "beta1_t1"), _column(row, "beta1_t2")),
                cphi_t=(_column(row, "cphi_t1"), _column(row, "cphi_t2")),
                beta2_25=_column(row, "beta2_25"),
                theta=_column(row, "theta"),
                psi_common_ion=_column(row, "psi_common_ion"),
                t_min=_column(row, "tmin"),
                t_max=_column(row, "tmax"),
                reference=row["reference"].strip(),
            )
        )
    return tuple(out)


def pitzer_pair(first: str, second: str) -> PitzerRecord | None:
    """The record for an ion pair, either order round, or ``None``.

    Both orders because the table stores a pair once and the interaction is symmetric. A
    pair the table does not carry is ``None`` rather than a zero: a model that evaluated
    unstated parameters as zero would return an ideal-solution answer wearing a Pitzer
    model's name.
    """
    a, b = first.strip().lower(), second.strip().lower()
    for record in pitzer_parameters():
        if (record.ion1, record.ion2) in ((a, b), (b, a)):
            return record
    return None


@cache
def salts() -> tuple[SaltRecord, ...]:
    """Every row of the salt table, in file order."""
    out = []
    for row in _rows(find(COMPSALT_CSV).read_text(encoding="utf-8")):
        name = row["saltname"].strip()
        if not name:
            continue

        out.append(
            SaltRecord(
                name=name,
                cation=row["ion1"].strip().lower(),
                anion=row["ion2"].strip().lower(),
                cation_stoichiometry=_column(row, "stoc1"),
                anion_stoichiometry=_column(row, "stoc2"),
                ksp=(
                    _column(row, "kspwater"),
                    _column(row, "kspwater2"),
                    _column(row, "kspwater3"),
                    _column(row, "kspwater4"),
                    _column(row, "kspwater5"),
                ),
                volume_delta=_column(row, "vdelta"),
                water_stoichiometry=_column(row, "waterstoc"),
            )
        )
    return tuple(out)


def salt(name: str) -> SaltRecord | None:
    """The record for a salt by name, matched without regard to case or surrounding space."""
    key = name.strip().lower()
    for record in salts():
        if record.name.lower() == key:
            return record
    return None


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
    # The card's entry for the substance, and the parameters out of it: `ion` is a class
    # rather than a parameter, so it is read separately and the map below is a map.
    stated_component = card.component(key) if card is not None else None
    override = None if stated_component is None else stated_component.parameters
    stated = card.association_for(key) if card is not None else None
    # `None` is "the card says nothing", which is not `False`: a card silent on the
    # question leaves the databank's class in force.
    ion = card.is_ion(key) if card is not None else None

    if override is None and stated is None:
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
        if override is None:
            raise PropertyUnavailableError(
                name,
                "critical constants",
                f"the keycard states an association for it but no "
                f"{sorted(keycard.CUBIC_PARAMETERS)}, and a substance the databank does not "
                f"have needs every parameter a cubic reads. An association is added to a "
                f"fluid, not to a name.",
            )
        # **An ion is exempt from the cubic's parameters**, and that is the flag's whole
        # point: an ion has no meaningful constants, so requiring them would make a user
        # invent the filler `mixture_of` refuses for exactly that reason.
        missing = [] if ion else sorted(keycard.CUBIC_PARAMETERS - set(override))
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
            # An ion has no meaningful `Tc`, `Pc` or `omega`, so a card need not state
            # them and the entry carries zeros. That is not a silent default: `mixture_of`
            # refuses a cubic over an ion, so nothing ever reads them.
            Tc=override["Tc"] if not ion else ureg.Quantity(0.0, "K"),
            Pc=override["Pc"] if not ion else ureg.Quantity(0.0, "Pa"),
            omega=0.0 if ion else _as_float(override["omega"], key),
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
            reference_state=SOLVENT,
            # The card's own statement, or `other` - the class that permits a cubic, which
            # is what a card-added substance is unless it says otherwise.
            component_type=ION if ion else "other",
            ionic_charge=_card_charge(override),
            # The card states metres and the databank holds ångström; this is the crossing.
            deshmukh_mather_diameter=_card_diameter(override),
            dielectric=_card_dielectric(override) or (0.0, 0.0, 0.0, 0.0, 0.0),
            # A card-added substance has no table row to inherit an association from, so
            # the card's is the whole of it - or none, if it states none.
            association=None if stated is None else _card_association(None, stated),
            # A card carries no PC-SAFT set, and zero is how this table spells absence:
            # a card states the parameters a cubic reads, and a model needing `m`, `sigma`
            # and `epsilon/k` refuses rather than inventing them.
            m_saft=0.0,
            sigma_saft=0.0,
            epsik_saft=0.0,
            lambda_r_mie=0.0,
            lambda_a_mie=0.0,
            m_mie=0.0,
            sigma_mie=0.0,
            epsik_mie=0.0,
            citation=None,
            source="keycard",
        )

    return replace(
        base,
        Tc=base.Tc if override is None else override.get("Tc", base.Tc),
        Pc=base.Pc if override is None else override.get("Pc", base.Pc),
        omega=(
            base.omega
            if override is None or "omega" not in override
            else _as_float(override["omega"], key)
        ),
        cp=base.cp if override is None else (_cp(override) or base.cp),
        association=(
            base.association if stated is None else _card_association(base.association, stated)
        ),
        # **The class is the table's and a card never changes it.** It is what decides
        # whether a cubic may be built at all, so a card able to clear it could hand a
        # cubic the filler the refusal exists to keep out of one.
        component_type=base.component_type,
        ionic_charge=(
            base.ionic_charge if override is None else _card_charge(override, base.ionic_charge)
        ),
        deshmukh_mather_diameter=(
            base.deshmukh_mather_diameter
            if override is None
            else _card_diameter(override, base.deshmukh_mather_diameter)
        ),
        dielectric=(
            base.dielectric if override is None else (_card_dielectric(override) or base.dielectric)
        ),
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


def _card_charge(override: Mapping[str, Q], shipped: float = 0.0) -> float:
    """The charge number the card states, or the databank's.

    A charge number is dimensionless - Pitzer's ionic strength is `1/2 sum m z^2` and
    NeqSim's own column holds the same number - so the canonical unit is
    ``dimensionless`` and the magnitude is the value a model takes.
    """
    if "ionic_charge" not in override:
        return shipped
    return float(override["ionic_charge"].to("dimensionless").magnitude)


def _card_diameter(override: Mapping[str, Q], shipped_angstrom: float = 0.0) -> float:
    """The ion diameter the card states, **in ångström whatever the card wrote**.

    A card states it in the card's own length units and the databank holds ångström -
    NeqSim's `ComponentDesmukhMather` multiplies its column by `1e-10` at the point of
    use - so this is the one crossing, and it is here rather than at every read.
    """
    if "deshmukh_mather_diameter" not in override:
        return shipped_angstrom
    return float(override["deshmukh_mather_diameter"].to("angstrom").magnitude)


def _card_dielectric(
    override: Mapping[str, Q],
) -> tuple[float, float, float, float, float] | None:
    """The five dielectric coefficients if the card states all of them, else `None`.

    Each is already a quantity in its own unit - the polynomial's `d0` is dimensionless,
    `d1` is a temperature, and `d2`..`d4` are inverse temperatures - so the magnitudes are
    the values the correlation takes.
    """
    names = ("dielectric_1", "dielectric_2", "dielectric_3", "dielectric_4", "dielectric_5")
    if all(name in override for name in names):
        return (
            float(override["dielectric_1"].magnitude),
            float(override["dielectric_2"].to("K").magnitude),
            float(override["dielectric_3"].to("1/K").magnitude),
            float(override["dielectric_4"].to("1/K**2").magnitude),
            float(override["dielectric_5"].to("1/K**3").magnitude),
        )
    return None


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
    names: tuple[str, ...], *, eos: str = "pr", card: keycard.Keycard | None = None
) -> dict[tuple[int, int], float]:
    """The interaction pairs the databank knows, for a list of components.

    Only pairs where both names are present are returned, and only where a value
    exists - an unlisted pair is zero, which is the ideal-mixture assumption and is
    what `mixture()` already does with an omitted pair.

    **A keycard's value wins over the databank's**, including when the keycard's
    value is exactly zero: overriding a fitted pair back to ideal mixing is a
    deliberate act, and a rule that treated zero as "absent" would silently undo it.

    Keyed by index into `names`, which is the form `mixture()` takes.

    ``eos`` selects the column, as it does in NeqSim: only Peng-Robinson reads
    ``KIJPR`` and every other cubic reads ``KIJSRK``.
    """
    # NeqSim's own rule: exactly its `PhasePrEos` reads `KIJPR`, and every other
    # phase - SRK, RK, and even `PhaseUMRCPA`, which extends `PhasePrEos`
    # without being it - reads `KIJSRK`.
    column = 1 if eos.strip().lower() == "pr" else 0
    pairs: dict[tuple[int, int], float] = {}
    for i, a in enumerate(names):
        for j in range(i + 1, len(names)):
            from_keycard = card.kij_for(a, names[j]) if card is not None else None
            if from_keycard is not None:
                # A card's zero is a value, not an absence - see the docstring.
                pairs[(i, j)] = from_keycard
                continue
            stored = _kij().get((a.strip().lower(), names[j].strip().lower()))
            if stored is not None and stored[column] != 0.0:
                pairs[(i, j)] = stored[column]
    return pairs


#: The associating interaction columns, and the cubic families each belongs to. The two
#: are separate fits, and a family is chosen by the cubic the mixture runs - Soave's shape
#: reads `cpakij_srk` and Peng-Robinson's `cpakij_pr`.
CPA_FAMILIES: tuple[str, ...] = ("srk", "pr")


def cpa_kij_for(
    names: tuple[str, ...], family: str, *, card: keycard.Keycard | None = None
) -> dict[tuple[int, int], float]:
    """The *associating* interaction pairs for a list of components, by index.

    The sibling of :func:`kij_for` for the column an associating model reads, with the
    same rules: an unlisted pair is zero, and a keycard's value wins over the databank's
    including when it is exactly zero.

    ``family`` is one of :data:`CPA_FAMILIES` and is required rather than defaulted,
    because the two columns are separate fits and a default would be a silent choice
    between two different fluids.

    Raises:
        InvalidInputError: if `family` is not one of :data:`CPA_FAMILIES`.
    """
    if family not in CPA_FAMILIES:
        raise InvalidInputError(
            "family",
            f"{family!r} is not an associating cubic family; expected one of {list(CPA_FAMILIES)}",
        )
    column = CPA_FAMILIES.index(family)
    pairs: dict[tuple[int, int], float] = {}
    for i, a in enumerate(names):
        for j in range(i + 1, len(names)):
            from_keycard = card.cpa_kij_for(a, names[j], family) if card is not None else None
            if from_keycard is not None:
                pairs[(i, j)] = from_keycard
                continue
            stored = _cpa_kij().get((a.strip().lower(), names[j].strip().lower()))
            if stored is not None and stored[column] != 0.0:
                pairs[(i, j)] = stored[column]
    return pairs


@dataclass(frozen=True, slots=True)
class GeNrtlPhaseParameters:
    """The parameters of an NRTL activity-coefficient *phase*.

    The NRTL matrices :class:`NrtlParameters` carries, beside what a phase needs and an
    activity coefficient does not: each component's pure-liquid vapour pressure, because
    the phase's fugacity coefficient is ``gamma_i P0_i / P``.

    **Flat rather than a nested record per component.** The transport carries a record's
    own fields flattened, one list per field, so a list of records would have nothing to
    cross in. The vapour-pressure columns are therefore parallel to the NRTL matrices
    and to each other, five coefficients per component in ``antoine_coefficients``.
    """

    #: ``alpha[i][j]``, ``N x N`` row-major. Symmetric with a zero diagonal.
    alpha: tuple[float, ...]
    #: ``dij[i][j] = g_ij`` in kelvin, ``N x N`` row-major. Directional.
    dij: tuple[float, ...]
    #: NeqSim's Antoine label for each component, in order.
    antoine_type: tuple[str, ...]
    #: The five coefficients ``A``-``E`` of each component, component-major.
    antoine_coefficients: tuple[float, ...]
    #: Critical temperature of each component, in K.
    antoine_tc: tuple[float, ...]
    #: Critical pressure of each component, in Pa.
    antoine_pc: tuple[float, ...]


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
    return _unifac_interaction(UNIFAC_INTER_CSV)


def _unifac_interaction(path: str) -> dict[tuple[int, int], float]:
    """One main-group interaction table, keyed by main-group pair.

    The three of them - `a`, and UNIFAC-PSRK's `b` and `c` - have the same shape, so
    they are read the same way rather than three times over.
    """
    out: dict[tuple[int, int], float] = {}
    for row in _rows(find(path).read_text(encoding="utf-8")):
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


@dataclass(frozen=True, slots=True)
class UnifacPsrkParameters:
    """The UNIFAC-PSRK parameters of a mixture, each matrix flattened row-major.

    The same basis as :class:`UnifacParameters`, with the interaction split into the
    three terms UNIFAC-PSRK fits separately: ``a_mn(T) = a_mn + b_mn T + c_mn T**2``.
    """

    #: Per-component group counts, ``N x G`` row-major.
    groups: tuple[float, ...]
    #: The volume ``R`` of each group, length ``G``.
    group_r: tuple[float, ...]
    #: The surface area ``Q`` of each group, length ``G``.
    group_q: tuple[float, ...]
    #: The constant term of the interaction, ``G x G`` row-major, in kelvin.
    aij: tuple[float, ...]
    #: The linear term, ``G x G`` row-major, in kelvin per kelvin.
    bij: tuple[float, ...]
    #: The quadratic term, ``G x G`` row-major, in kelvin per kelvin squared.
    cij: tuple[float, ...]


def unifac_psrk_parameters(names: Sequence[str]) -> UnifacPsrkParameters:
    """The UNIFAC-PSRK inputs for a list of components, by name.

    ``a``, ``b`` and ``c`` are resolved from NeqSim's ``UNIFACInterParam``,
    ``UNIFACInterParamB`` and ``UNIFACInterParamC`` over the same group union, so the
    three matrices line up column for column. The interaction is
    ``a_mn(T) = a_mn + b_mn T + c_mn T**2``, NeqSim's
    ``ComponentGEUnifacPSRK.calcaij``.

    A pair the tables do not carry is zero in all three terms, so a mixture whose pairs
    all have ``b = c = 0`` reduces this to plain :func:`unifac_parameters` exactly.

    Raises:
        PropertyUnavailableError: if a name has no UNIFAC group assignment.
    """
    group = _unifac_group()
    members = _unifac_members()
    a_table = _unifac_aij()
    b_table = _unifac_interaction(UNIFAC_INTER_B_CSV)
    c_table = _unifac_interaction(UNIFAC_INTER_C_CSV)

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
    bij = [0.0] * (g * g)
    cij = [0.0] * (g * g)
    for k, subgroup in enumerate(union):
        r, q, main = group[subgroup]
        group_r[k] = r
        group_q[k] = q
        for m, other in enumerate(union):
            _, _, other_main = group[other]
            aij[k * g + m] = a_table.get((main, other_main), 0.0)
            bij[k * g + m] = b_table.get((main, other_main), 0.0)
            cij[k * g + m] = c_table.get((main, other_main), 0.0)

    groups = [0.0] * (len(names) * g)
    for i, name in enumerate(names):
        for subgroup, count in members[name.strip().lower()]:
            groups[i * g + union.index(subgroup)] = float(count)

    return UnifacPsrkParameters(
        groups=tuple(groups),
        group_r=tuple(group_r),
        group_q=tuple(group_q),
        aij=tuple(aij),
        bij=tuple(bij),
        cij=tuple(cij),
    )


@dataclass(frozen=True, slots=True)
class UnifacUmrpruParameters:
    """The UNIFAC-UMR-PRU parameters of a mixture, each matrix flattened row-major.

    The interaction is evaluated about 298.15 K rather than about zero:
    ``a_mn(T) = a_mn + b_mn (T - 298.15) + c_mn (T - 298.15)**2``.
    """

    #: Per-component group counts, ``N x G`` row-major.
    groups: tuple[float, ...]
    #: The volume ``R`` of each group, length ``G``.
    group_r: tuple[float, ...]
    #: The surface area ``Q`` of each group, length ``G``.
    group_q: tuple[float, ...]
    #: The constant term of the interaction, ``G x G`` row-major, in kelvin.
    aij: tuple[float, ...]
    #: The linear term, ``G x G`` row-major, in kelvin per kelvin.
    bij: tuple[float, ...]
    #: The quadratic term, ``G x G`` row-major, in kelvin per kelvin squared.
    cij: tuple[float, ...]


#: The UMR-PRU interaction tables of each parameter set, keyed by the name the spec's
#: enum declares. Repo-relative, because that is how every data file is addressed here.
_UMRPRU_SETS: dict[str, tuple[str, str, str]] = {
    "umr": (
        "data/components/UNIFACInterParamA_UMR.csv",
        "data/components/UNIFACInterParamB_UMR.csv",
        "data/components/UNIFACInterParamC_UMR.csv",
    ),
    "umrmc": (
        "data/components/UNIFACInterParamA_UMRMC.csv",
        "data/components/UNIFACInterParamB_UMRMC.csv",
        "data/components/UNIFACInterParamC_UMRMC.csv",
    ),
}


@cache
def _umrpru_members() -> dict[str, tuple[tuple[int, int], ...]]:
    """UMR-PRU group memberships, from the 139-subgroup decomposition."""
    out: dict[str, tuple[tuple[int, int], ...]] = {}
    for row in _rows(find(UNIFAC_COMP_UMRPRU_CSV).read_text(encoding="utf-8")):
        name = row["name"].strip().lower()
        subs = tuple(
            (int(key[3:]), int(value))
            for key, value in row.items()
            if key.startswith("sub") and value.strip() not in ("", "0")
        )
        out[name] = subs
    return out


def unifac_umrpru_parameters(names: Sequence[str], parameters: str) -> UnifacUmrpruParameters:
    """The UNIFAC-UMR-PRU inputs for a list of components, by name and parameter set.

    ``parameters`` names which interaction set to read, ``"umr"`` or ``"umrmc"``.
    NeqSim decides this from ``getComponent(0).getAttractiveTermNumber()`` - the
    ``_umrmc`` tables when it is 13, 19 or 22 - and a component here carries no such
    field, so the caller states which equation of state they are pairing the model with.

    The group decomposition is NeqSim's ``UNIFACcompUMRPRU``, which carries 139
    subgroups rather than the 133 of ``UNIFACcomp``; the group constants are the shared
    table, so only the decomposition and the interaction differ.

    Raises:
        PropertyUnavailableError: if a name has no UMR-PRU group assignment.
        InvalidInputError: if ``parameters`` names no set.
    """
    try:
        a_path, b_path, c_path = _UMRPRU_SETS[parameters]
    except KeyError:
        raise InvalidInputError(
            "parameters",
            f"unknown UMR-PRU parameter set {parameters!r}; expected one of {sorted(_UMRPRU_SETS)}",
        ) from None

    group = _unifac_group()
    members = _umrpru_members()
    a_table = _unifac_interaction(a_path)
    b_table = _unifac_interaction(b_path)
    c_table = _unifac_interaction(c_path)

    union: list[int] = []
    for name in names:
        subs = members.get(name.strip().lower())
        if subs is None:
            raise PropertyUnavailableError(
                name,
                "UNIFAC-UMR-PRU group assignment",
                "not in UNIFACcompUMRPRU.csv; a UNIFAC activity coefficient needs a "
                "group decomposition for every component",
            )
        for subgroup, _ in subs:
            if subgroup not in union:
                union.append(subgroup)
    union.sort()
    g = len(union)

    group_r = [0.0] * g
    group_q = [0.0] * g
    aij = [0.0] * (g * g)
    bij = [0.0] * (g * g)
    cij = [0.0] * (g * g)
    for k, subgroup in enumerate(union):
        r, q, main = group[subgroup]
        group_r[k] = r
        group_q[k] = q
        for m, other in enumerate(union):
            _, _, other_main = group[other]
            aij[k * g + m] = a_table.get((main, other_main), 0.0)
            bij[k * g + m] = b_table.get((main, other_main), 0.0)
            cij[k * g + m] = c_table.get((main, other_main), 0.0)

    groups = [0.0] * (len(names) * g)
    for i, name in enumerate(names):
        for subgroup, count in members[name.strip().lower()]:
            groups[i * g + union.index(subgroup)] = float(count)

    return UnifacUmrpruParameters(
        groups=tuple(groups),
        group_r=tuple(group_r),
        group_q=tuple(group_q),
        aij=tuple(aij),
        bij=tuple(bij),
        cij=tuple(cij),
    )


class _PhaseAntoine(NamedTuple):
    """The per-component vapour-pressure columns every GE phase carries, in parallel."""

    kinds: tuple[str, ...]
    coefficients: tuple[float, ...]
    tcs: tuple[float, ...]
    pcs: tuple[float, ...]


def _phase_antoine(names: Sequence[str], card: keycard.Keycard | None) -> _PhaseAntoine:
    """The vapour-pressure records a phase reads, or a refusal.

    Shared by every activity-coefficient phase, because the two things it checks are the
    same for all of them: a component tagged anything but ``SOLVENT`` takes a Henry's-law
    coefficient in ``ComponentGE.fugcoef``, which this library does not implement, and a
    component with no correlation has no ``P0`` to compose.
    """
    kinds: list[str] = []
    coefficients: list[float] = []
    tcs: list[float] = []
    pcs: list[float] = []
    for name in names:
        record = entry(name, card=card)
        if record.reference_state != SOLVENT:
            raise InvalidInputError(
                "components",
                f"`{record.name}` is tagged `referenceStateType = "
                f"{record.reference_state}` in NeqSim's component database, so "
                f"`ComponentGE.fugcoef` gives it a Henry's-law fugacity coefficient "
                f"rather than `gamma_i P0_i / P`. That branch is not ported, and this "
                f"phase is the Raoult one. The substances the databank tags `solvent` - "
                f"water, the alcohols, the glycols - are the ones it describes.",
            )
        if record.antoine == (0.0, 0.0, 0.0, 0.0, 0.0):
            raise PropertyUnavailableError(
                name,
                "Antoine vapour-pressure coefficients",
                "the phase's fugacity coefficient is `gamma_i P0_i / P`, so every "
                "component needs a correlation; a keycard supplies the parameters a "
                "cubic reads and not these",
            )
        kinds.append(record.antoine_form())
        coefficients.extend(record.antoine)
        tcs.append(record.Tc.to_base_units().magnitude)
        pcs.append(record.Pc.to_base_units().magnitude)
    return _PhaseAntoine(tuple(kinds), tuple(coefficients), tuple(tcs), tuple(pcs))


@dataclass(frozen=True, slots=True)
class GeVanLaarAcidPhaseParameters:
    """The parameters of a Van Laar acid activity-coefficient *phase*.

    :class:`VanLaarAcidParameters`' acid identity, beside each component's Antoine columns
    and the critical constants they need. The columns are carried for *every* component
    but read only for the ones with ``acid_index == 0``: the three modelled species take
    their ``P0`` from :func:`azoth.eos.nitric_sulfuric_acid_vapor_pressure` instead, which
    is the whole point of the phase.
    """

    #: Per component: ``1`` water, ``2`` nitric acid, ``3`` sulfuric acid, ``0`` for a
    #: species the model does not cover.
    acid_index: tuple[int, ...]
    #: NeqSim's Antoine label for each component, in order. Read only where
    #: ``acid_index`` is zero.
    antoine_type: tuple[str, ...]
    #: The five coefficients ``A``-``E`` of each component, component-major.
    antoine_coefficients: tuple[float, ...]
    #: Critical temperature of each component, in K.
    antoine_tc: tuple[float, ...]
    #: Critical pressure of each component, in Pa.
    antoine_pc: tuple[float, ...]


def ge_van_laar_acid_phase_parameters(
    names: Sequence[str], *, card: keycard.Keycard | None = None
) -> GeVanLaarAcidPhaseParameters:
    """The parameters of a Van Laar acid phase for a list of components, by name.

    **This resolver does not refuse a Henry's-law solute, and its four siblings do.** The
    refusal there exists because ``ComponentGE.fugcoef`` branches on ``referenceStateType``,
    so computing the Raoult expression for a ``solute`` component would be the wrong branch
    with no symptom. ``ComponentGEVanLaarAcid`` *overrides* ``fugcoef`` to ignore the tag -
    and it has to, because both acids are tagged ``solute``. What this refuses is a
    component the model does not cover that also has no Antoine correlation, because then
    its ``P0`` has no source.

    Raises:
        PropertyUnavailableError: if a name is in neither the databank nor the keycard, or
            if an uncovered component carries no Antoine correlation.
    """
    acid = van_laar_acid_parameters(names)

    kinds: list[str] = []
    coefficients: list[float] = []
    tcs: list[float] = []
    pcs: list[float] = []
    for name, index in zip(names, acid.acid_index, strict=True):
        record = entry(name, card=card)
        if index == 0 and record.antoine == (0.0, 0.0, 0.0, 0.0, 0.0):
            raise PropertyUnavailableError(
                name,
                "Antoine vapour-pressure coefficients",
                "a species the acid model does not cover takes its `P0` from Antoine, and "
                "the databank carries none for this one",
            )
        kinds.append(record.antoine_form())
        coefficients.extend(record.antoine)
        tcs.append(record.Tc.to_base_units().magnitude)
        pcs.append(record.Pc.to_base_units().magnitude)

    return GeVanLaarAcidPhaseParameters(
        acid_index=acid.acid_index,
        antoine_type=tuple(kinds),
        antoine_coefficients=tuple(coefficients),
        antoine_tc=tuple(tcs),
        antoine_pc=tuple(pcs),
    )


@dataclass(frozen=True, slots=True)
class GeUniquacPhaseParameters:
    """The parameters of a UNIQUAC activity-coefficient *phase*.

    The volume and surface parameters :class:`UniquacParameters` carries, beside what a
    phase needs and an activity coefficient does not: each component's pure-liquid vapour
    pressure, because the phase's fugacity coefficient is ``gamma_i P0_i / P``.

    ``aij`` is **not** here. It is the caller's - no upstream table carries a UNIQUAC
    interaction matrix - and it arrives as its own argument, the way it does for
    :func:`azoth.eos.uniquac_activity_coefficients`.
    """

    #: The van der Waals volume parameter ``r_i`` of each component.
    r: tuple[float, ...]
    #: The van der Waals surface-area parameter ``q_i`` of each component.
    q: tuple[float, ...]
    #: NeqSim's Antoine label for each component, in order.
    antoine_type: tuple[str, ...]
    #: The five coefficients ``A``-``E`` of each component, component-major.
    antoine_coefficients: tuple[float, ...]
    #: Critical temperature of each component, in K.
    antoine_tc: tuple[float, ...]
    #: Critical pressure of each component, in Pa.
    antoine_pc: tuple[float, ...]


def ge_uniquac_phase_parameters(
    names: Sequence[str], *, card: keycard.Keycard | None = None
) -> GeUniquacPhaseParameters:
    """The parameters of a UNIQUAC phase for a list of components, by name.

    ``r`` and ``q`` come from :func:`uniquac_parameters`, so the phase and
    ``eos.uniquac_activity_coefficients`` cannot disagree about them; this adds each
    component's Antoine vapour-pressure correlation beside them.

    Raises:
        PropertyUnavailableError: if a name has no UNIFAC group assignment, or if the
            databank carries it without an Antoine correlation.
        InvalidInputError: if a name is a Henry's-law solute, which this phase does not
            implement.
    """
    uniquac = uniquac_parameters(names)
    antoine = _phase_antoine(names, card)
    return GeUniquacPhaseParameters(
        r=uniquac.r,
        q=uniquac.q,
        antoine_type=antoine.kinds,
        antoine_coefficients=antoine.coefficients,
        antoine_tc=antoine.tcs,
        antoine_pc=antoine.pcs,
    )


@dataclass(frozen=True, slots=True)
class GeWilsonPhaseParameters:
    """The parameters of a Wilson activity-coefficient *phase*.

    Only the vapour-pressure columns. The Wilson correlation reads each component's molar
    mass and critical temperature, which are what :class:`~azoth.eos.mixture.Mixture`
    carries and :func:`mixture_of` already resolves, so the phase takes that mixture
    beside this record rather than a second copy of the same two vectors.
    """

    #: NeqSim's Antoine label for each component, in order.
    antoine_type: tuple[str, ...]
    #: The five coefficients ``A``-``E`` of each component, component-major.
    antoine_coefficients: tuple[float, ...]
    #: Critical temperature of each component, in K.
    antoine_tc: tuple[float, ...]
    #: Critical pressure of each component, in Pa.
    antoine_pc: tuple[float, ...]


def ge_wilson_phase_parameters(
    names: Sequence[str], *, card: keycard.Keycard | None = None
) -> GeWilsonPhaseParameters:
    """The parameters of a Wilson phase for a list of components, by name.

    Raises:
        PropertyUnavailableError: if a name is in neither the databank nor the keycard, or
            if the databank carries it without an Antoine correlation.
        InvalidInputError: if a name is a Henry's-law solute, which this phase does not
            implement.
    """
    antoine = _phase_antoine(names, card)
    return GeWilsonPhaseParameters(
        antoine_type=antoine.kinds,
        antoine_coefficients=antoine.coefficients,
        antoine_tc=antoine.tcs,
        antoine_pc=antoine.pcs,
    )


@dataclass(frozen=True, slots=True)
class GeUnifacPhaseParameters:
    """The parameters of a UNIFAC activity-coefficient *phase*.

    The group tables :class:`UnifacParameters` carries, beside what a phase needs and an
    activity coefficient does not: each component's pure-liquid vapour pressure, because
    the phase's fugacity coefficient is ``gamma_i P0_i / P``.
    """

    #: Per-component group counts, ``N x G`` row-major.
    groups: tuple[float, ...]
    #: The volume ``R`` of each group, length ``G``.
    group_r: tuple[float, ...]
    #: The surface area ``Q`` of each group, length ``G``.
    group_q: tuple[float, ...]
    #: The main-group interaction matrix, ``G x G`` row-major, in kelvin.
    aij: tuple[float, ...]
    #: NeqSim's Antoine label for each component, in order.
    antoine_type: tuple[str, ...]
    #: The five coefficients ``A``-``E`` of each component, component-major.
    antoine_coefficients: tuple[float, ...]
    #: Critical temperature of each component, in K.
    antoine_tc: tuple[float, ...]
    #: Critical pressure of each component, in Pa.
    antoine_pc: tuple[float, ...]


def ge_unifac_phase_parameters(
    names: Sequence[str], *, card: keycard.Keycard | None = None
) -> GeUnifacPhaseParameters:
    """The parameters of a UNIFAC phase for a list of components, by name.

    The group tables come from :func:`unifac_parameters`, so the phase and
    ``eos.unifac_activity_coefficients`` cannot disagree about them; this adds each
    component's Antoine vapour-pressure correlation beside them.

    Raises:
        PropertyUnavailableError: if a name has no UNIFAC group assignment, or if the
            databank carries it without an Antoine correlation.
        InvalidInputError: if a name is a Henry's-law solute, which this phase does not
            implement.
    """
    unifac = unifac_parameters(names)
    antoine = _phase_antoine(names, card)
    return GeUnifacPhaseParameters(
        groups=unifac.groups,
        group_r=unifac.group_r,
        group_q=unifac.group_q,
        aij=unifac.aij,
        antoine_type=antoine.kinds,
        antoine_coefficients=antoine.coefficients,
        antoine_tc=antoine.tcs,
        antoine_pc=antoine.pcs,
    )


def ge_nrtl_phase_parameters(
    names: Sequence[str], *, card: keycard.Keycard | None = None
) -> GeNrtlPhaseParameters:
    """The parameters of an NRTL phase for a list of components, by name.

    The NRTL matrices come from :func:`nrtl_parameters`, so the phase and
    ``eos.nrtl_activity_coefficients`` cannot disagree about them; this adds each
    component's Antoine vapour-pressure correlation beside them.

    Raises:
        PropertyUnavailableError: if a name is in neither the databank nor the keycard,
            or if the databank carries it without an Antoine correlation - which is what
            a keycard-added substance is, since a card states the parameters a cubic
            reads and not these.
    """
    nrtl = nrtl_parameters(names)
    antoine = _phase_antoine(names, card)

    return GeNrtlPhaseParameters(
        alpha=nrtl.alpha,
        dij=nrtl.dij,
        antoine_type=antoine.kinds,
        antoine_coefficients=antoine.coefficients,
        antoine_tc=antoine.tcs,
        antoine_pc=antoine.pcs,
    )


def _interaction_pairs(
    resolved: list[str],
    eos: str,
    *,
    associating: bool,
    card: keycard.Keycard | None,
) -> dict[tuple[int, int], float]:
    """:func:`kij_for` or :func:`cpa_kij_for`, whichever the model reads.

    One function rather than a ternary at each of the two callers, because the choice is
    the *same* choice in both and a second copy of it is a second place to get the family
    wrong.
    """
    if associating:
        return cpa_kij_for(tuple(resolved), family_of(eos), card=card)
    return kij_for(tuple(resolved), eos=eos, card=card)


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
    associating: bool = False,
) -> Mixture:
    """A :class:`~azoth.eos.mixture.Mixture` from a list of databank names.

    The interaction parameters come from the databank too, so a caller writing
    `from_names(["methane", "n-butane"])` gets the published `kij` rather than a
    silent zero.

    A mixture alone, for a caller who has their own heat-capacity coefficients or wants
    only to flash. :func:`mixture_of` returns the polynomial beside it, and is what a
    calculation takes - a mixture without one cannot produce an enthalpy.

    ``associating`` selects the *associating* interaction column and turns the Wertheim
    contribution on, which is the whole of what makes one of these a CPA fluid: the
    attraction is mixed with `cpakij_SRK`/`cpakij_PR` rather than `KIJPR`, and the fitted
    `a` and `b` replace the cubic's. It defaults to off because the same substances are an
    associating fluid under one phase model and a classical one under another.

    Raises:
        PropertyUnavailableError: if any name is not in the databank.
        InvalidInputError: if the list is empty, or a pair is malformed.
    """
    cubic = _cubic(eos)
    alpha_name = _alpha_name(eos, alpha)
    resolved = [name.strip().lower() for name in names]
    entries = [entry(name, card=card) for name in resolved]
    _refuse_ions([e.name for e in entries if e.component_type == ION])
    components = tuple(e.component(alpha=alpha_name) for e in entries)
    return mixture(
        components,
        kij=_interaction_pairs(resolved, eos, associating=associating, card=card),
        cubic=cubic,
        alpha=alpha_name,
        associating=associating,
    )


def umr_cpa_mixture_of(
    names: list[str],
    *,
    card: keycard.Keycard | None = None,
) -> tuple[Mixture, IdealGasModel]:
    """The UMR-CPA fluid: Peng-Robinson, the UMR mixing rule and the association.

    NeqSim's ``SystemUMRCPAEoS``, and the only mixture in this library whose attraction is
    mixed by a universal rule rather than an interaction matrix. Three choices are made
    here and nowhere else:

    * the mixing rule is ``"umr"`` over the ``_umrmc`` UNIFAC tables, which is what
      ``SystemUMRCPAEoS``'s construction and ``SystemThermo``'s ``"UMR-CPA"`` selector
      both name, and the interaction matrix is **zero** - the rule reads no ``kij`` column;
    * the alpha is the five-parameter Mathias-Copeman term seeded with ``UMRCPA_MC1..5``,
      which is ``ComponentUMRCPA.setAttractiveTerm``'s term 22;
    * a component whose row carries no ``UMRCPA_MC`` set is **refused**. NeqSim falls back
      per component - term 19 seeded with ``MCPR1..3`` for a non-associating one, term 1
      with ``mCPA`` for an associating one - and this library's alpha is one per mixture,
      so a mixture of a component that carries the set and one that does not could not be
      expressed at all.

    Raises:
        PropertyUnavailableError: if a name is not in the databank, has no
            ``UMRCPA_MC1..5`` set, or has no ``UNIFACcompUMRPRU`` group decomposition.
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

    alpha = "matcop_5prumr"
    for e in entries:
        coefficients = e.alpha_params.get(alpha, ())
        if not any(abs(c) > 1.0e-20 for c in coefficients):
            raise PropertyUnavailableError(
                e.name,
                "UMRCPA_MC1..5",
                "carries no UMR-CPA Mathias-Copeman set; the UMR-CPA model's attractive "
                "term is that set, and NeqSim's per-component fallback to term 19 or term "
                "1 is not expressible as one mixture-level alpha",
            )

    tables = unifac_umrpru_parameters(resolved, "umrmc")
    fluid = mixture(
        tuple(e.component(alpha=alpha) for e in entries),
        cubic=_cubic("pr"),
        alpha=alpha,
        associating=True,
        mixing_rule="umr",
        umr=tables,
    )
    return (
        fluid,
        IdealGasModel(
            cp_a=tuple(e.cp[0] for e in entries),  # type: ignore[index]
            cp_b=tuple(e.cp[1] for e in entries),  # type: ignore[index]
            cp_c=tuple(e.cp[2] for e in entries),  # type: ignore[index]
            cp_d=tuple(e.cp[3] for e in entries),  # type: ignore[index]
            cp_e=tuple(e.cp[4] for e in entries),  # type: ignore[index]
        ),
    )


def mixture_of(
    names: list[str],
    *,
    card: keycard.Keycard | None = None,
    eos: str = "pr",
    alpha: str | None = None,
    associating: bool = False,
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
    _refuse_ions([e.name for e in entries if e.component_type == ION])

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
            kij=_interaction_pairs(resolved, eos, associating=associating, card=card),
            cubic=_cubic(eos),
            alpha=_alpha_name(eos, alpha),
            associating=associating,
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
    "SELF_BONDLESS_SCHEMES",
    "SITE_SCHEMES",
    "SOLVENT",
    "AssociationParameters",
    "BwrsCoefficients",
    "DatabankEntry",
    "GeNrtlPhaseParameters",
    "GeUnifacPhaseParameters",
    "GeUniquacPhaseParameters",
    "GeVanLaarAcidPhaseParameters",
    "GeWilsonPhaseParameters",
    "NrtlParameters",
    "PitzerRecord",
    "SaltRecord",
    "UnifacParameters",
    "UnifacPsrkParameters",
    "UnifacUmrpruParameters",
    "UniquacParameters",
    "VanLaarAcidParameters",
    "available",
    "bwrs_coefficients",
    "component",
    "entry",
    "from_model",
    "from_names",
    "ge_nrtl_phase_parameters",
    "ge_unifac_phase_parameters",
    "ge_uniquac_phase_parameters",
    "ge_van_laar_acid_phase_parameters",
    "ge_wilson_phase_parameters",
    "kij_for",
    "nrtl_parameters",
    "pitzer_pair",
    "pitzer_parameters",
    "salt",
    "salts",
    "unifac_parameters",
    "unifac_psrk_parameters",
    "unifac_umrpru_parameters",
    "uniquac_parameters",
    "van_laar_acid_parameters",
    "wilke_chang_phi",
]
