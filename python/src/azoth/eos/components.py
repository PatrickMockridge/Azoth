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
from dataclasses import dataclass, replace
from functools import cache

from azoth import keycard
from azoth._data import find
from azoth.core.errors import InvalidInputError, PropertyUnavailableError
from azoth.core.units import Q, ureg
from azoth.eos.mixture import Component, Mixture
from azoth.eos.reference.molar_enthalpy_entropy import IdealGasModel

COMPONENTS_CSV = "data/components/components.csv"
KIJ_CSV = "data/components/kij.csv"

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
    citation: str | None
    #: Where these values came from: the vendored databank, or the keycard in force.
    #: Not part of a citation - it is the *provenance of the lookup*, which a caller
    #: needs when a result turns out to depend on which file was in play.
    source: str = "databank"

    def component(self) -> Component:
        """This entry as the thing a calculation takes."""
        return Component(Tc=self.Tc, Pc=self.Pc, omega=self.omega)

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
            citation=row["citation"],
        )
    return entries


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


def available(*, card: keycard.Keycard | None = None) -> tuple[str, ...]:
    """Every name available, sorted: the databank plus whatever a card adds.

    `card` is the card this call reads; the loaded one is consulted when none is
    passed. See :func:`azoth.keycard.in_force`.
    """
    card = keycard.in_force(card)
    extra = set(card.components) if card is not None else set()
    return tuple(sorted(set(_table()) | extra))


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
    card = keycard.in_force(card)
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
        missing = sorted(set(keycard.COMPONENT_PARAMETERS) - set(override))
        if missing:
            raise PropertyUnavailableError(
                name,
                "critical constants",
                f"the keycard supplies {sorted(override)} but a cubic needs "
                f"{sorted(keycard.COMPONENT_PARAMETERS)}; {missing} is missing. A "
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
            # A keycard supplies the parameters a *cubic* needs, and a polynomial is not
            # one of them. `None` says so rather than a row of zeros standing in for a
            # heat capacity - `mixture_of` is where a caller finds out, by name.
            cp=None,
            citation=None,
            source="keycard",
        )

    return replace(
        base,
        Tc=override.get("Tc", base.Tc),
        Pc=override.get("Pc", base.Pc),
        omega=_as_float(override["omega"], key) if "omega" in override else base.omega,
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
    card = keycard.in_force(card)
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


def from_names(names: list[str], *, card: keycard.Keycard | None = None) -> Mixture:
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
    resolved = [name.strip().lower() for name in names]
    components = tuple(component(name, card=card) for name in resolved)
    return mixture(components, kij=kij_for(tuple(resolved), card=card))


def mixture_of(
    names: list[str], *, card: keycard.Keycard | None = None
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
            tuple(e.component() for e in entries),
            kij=kij_for(tuple(resolved), card=card),
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
    card = keycard.in_force(card)
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


# Imported at the bottom because `mixture` lives with the types this module builds
# on, and importing it at the top would make the cycle explicit for no gain.
from azoth.eos.mixture import mixture  # noqa: E402

__all__ = [
    "DatabankEntry",
    "available",
    "component",
    "entry",
    "from_model",
    "from_names",
    "kij_for",
]
