"""The keycard: the file a user supplies to extend or override what azoth ships.

Sections: ``keyholder``, ``components``, ``associations``, ``kij``, ``fluids``,
``fittings``, ``coefficients`` and ``models``. A name the databank already has replaces
its critical constants; a name it does not have adds one.

    >>> import azoth
    >>> card = azoth.keycard.load("keycard.toml")
    >>> azoth.eos.component("methane", card=card)   # the card's values, not the databank's

**The card is a value.** :func:`load` returns what the file says and stores nothing;
every call that reads a card is handed one, and a call handed none reads the data
this library ships. There is no `current()` and no `clear()`. The library implements;
the engineer decides, and the card is where that decision is recorded.

This module does not validate against the JSON Schema - that is
``tools/check_user_data.py``, and ``test_keycard_loader.py`` holds the two to each
other.
"""

from __future__ import annotations

import tomllib
from collections.abc import Mapping
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import pint

from azoth.core._units_gen import UNIT_VOCABULARY as _UNIT_VOCABULARY
from azoth.core.errors import KeycardError
from azoth.core.units import Q, ureg

#: The format version this module reads. A keycard declaring anything else is refused
#: rather than guessed at.
SCHEMA_VERSION = 2

#: The component parameters this build reads, and the unit each is converted to.
#: A parameter outside this table is refused, because a value nothing reads is a value
#: silently ignored. The unit is checked and not assumed: a `Tc` in Celsius converted
#: as though it were kelvin is a plausible wrong number with no symptom.
COMPONENT_PARAMETERS: Mapping[str, str] = {
    "Tc": "K",
    "Pc": "Pa",
    "omega": "dimensionless",
    "cp_a": "J/(mol*K)",
    "cp_b": "J/(mol*K**2)",
    "cp_c": "J/(mol*K**3)",
    "cp_d": "J/(mol*K**4)",
    "cp_e": "J/(mol*K**5)",
}

#: The subset of `COMPONENT_PARAMETERS` a cubic reads, and therefore what a card must
#: supply to add a substance. The rest are what an *enthalpy* needs, optional for a
#: caller who only flashes.
CUBIC_PARAMETERS: frozenset[str] = frozenset({"Tc", "Pc", "omega"})

#: The association parameters this build reads, and the unit each is converted to,
#: declared in `specs/schema/association.schema.json`.
#:
#: **The fitted values are stated in SI, not in the source table's internal scale.** A
#: holder's regression is in ``Pa*m**6/mol**2`` and ``m**3/mol`` - water's attraction is
#: ``0.12277``, which is the published CPA value - while
#: ``data/components/components.csv`` carries ``12277.0``, which is NeqSim's own scale and
#: belongs to the table that wrote it. ``azoth.eos.components`` is the one crossing.
ASSOCIATION_PARAMETERS: Mapping[str, str] = {
    "energy": "J/mol",
    "volume_srk": "dimensionless",
    "a_srk": "Pa*m**6/mol**2",
    "b_srk": "m**3/mol",
    "m_srk": "dimensionless",
    "volume_pr": "dimensionless",
    "a_pr": "Pa*m**6/mol**2",
    "b_pr": "m**3/mol",
    "m_pr": "dimensionless",
}

#: The association site schemes this build implements, and the schema's enum.
#:
#: The names are the databank's, and they are not interchangeable with the site *count*:
#: ``2A`` and ``2B`` both carry two sites and bond differently, and ``1A`` and ``2A`` carry
#: no bonding pair at all - so a card states the name and the kernel reads the name.
ASSOCIATION_SCHEMES: tuple[str, ...] = ("1A", "2A", "2B", "4C")

#: The units a keycard may declare, compiled from `specs/vocabulary/vocabulary.toml`
#: into `azoth.core._units_gen`. The schema's enum is generated from the same table;
#: a unit outside it fails at load rather than at first use.
UNIT_VOCABULARY: tuple[str, ...] = _UNIT_VOCABULARY

#: The model vocabularies this build implements, each equal to the schema's enum.
MODEL_KINDS = ("cubic_eos",)
MODEL_SHAPES = ("peng_robinson",)
MODEL_ALPHAS = ("peng_robinson",)
MODEL_MIXING_RULES = ("classical_kij",)

#: Every key a model definition may carry, which is the schema's
#: `models.additionalProperties.properties` plus `required`.
MODEL_KEYS = frozenset({"kind", "shape", "alpha", "mixing_rule", "components", "alpha_parameters"})


@dataclass(frozen=True, slots=True)
class Model:
    """One declarative model definition: named choices, and nothing to execute."""

    name: str
    kind: str
    shape: str
    alpha: str
    mixing_rule: str
    components: tuple[str, ...]

    def __repr__(self) -> str:
        return f"Model({self.name!r}, {self.shape}, {self.alpha}, {list(self.components)})"


@dataclass(frozen=True, slots=True)
class Association:
    """One substance's association, as a card states it.

    ``parameters`` is keyed by the names in :data:`ASSOCIATION_PARAMETERS` and holds each
    value in that parameter's canonical unit. ``scheme`` is separate because it is a name
    from a closed vocabulary where every parameter beside it is a number in a unit.
    """

    scheme: str
    parameters: Mapping[str, Q]

    def __repr__(self) -> str:
        return f"Association({self.scheme!r}, {sorted(self.parameters)})"


@dataclass(frozen=True, slots=True)
class Keycard:
    """A loaded keycard. Immutable, and the whole of what a file said."""

    keyholder: str | None
    licence: str | None
    components: Mapping[str, Mapping[str, Q]] = field(default_factory=dict)
    associations: Mapping[str, Association] = field(default_factory=dict)
    kij: Mapping[tuple[str, str], float] = field(default_factory=dict)
    cpa_kij: Mapping[tuple[str, str, str], float] = field(default_factory=dict)
    coefficients: Mapping[str, Mapping[str, Q]] = field(default_factory=dict)
    conventions: Mapping[str, Mapping[str, str]] = field(default_factory=dict)
    models: Mapping[str, Model] = field(default_factory=dict)
    fittings: tuple[Mapping[str, Any], ...] = ()
    fluids: Mapping[str, tuple[Mapping[str, Any], ...]] = field(default_factory=dict)
    path: Path | None = None

    def component(self, name: str) -> Mapping[str, Q] | None:
        """The parameters this keycard gives a substance, or ``None``."""
        return self.components.get(name.strip().lower())

    def association_for(self, name: str) -> Association | None:
        """The association this card states for a substance, or ``None``."""
        return self.associations.get(name.strip().lower())

    def names(self) -> set[str]:
        """Every substance this card states anything about.

        The union of the two sections that name components, because a card may state a
        scheme for a substance whose critical constants it does not touch - and a caller
        comparing the two readers has to ask about the same set on both sides.
        """
        return set(self.components) | set(self.associations)

    def model(self, name: str) -> Model | None:
        """A model definition by name, or ``None``."""
        return self.models.get(name)

    def coefficient(self, calc_id: str, name: str) -> Q | None:
        """A coefficient this keycard supplies for one calculation's argument."""
        return self.coefficients.get(calc_id, {}).get(name)

    def kij_for(self, a: str, b: str) -> float | None:
        """An interaction parameter for a pair, in either order, or ``None``."""
        first, second = sorted((a.strip().lower(), b.strip().lower()))
        return self.kij.get((first, second))

    def cpa_kij_for(self, a: str, b: str, family: str) -> float | None:
        """The *associating* interaction parameter for a pair at one cubic family.

        ``family`` is ``"srk"`` or ``"pr"``, and it is part of the key because the two
        are separate fits rather than one converted: water/methanol is ``-0.153`` in
        both here, but the column is the source's own and nothing derives one from the
        other.
        """
        first, second = sorted((a.strip().lower(), b.strip().lower()))
        return self.cpa_kij.get((first, second, family))

    def __repr__(self) -> str:
        holder = self.keyholder or "unnamed"
        return (
            f"Keycard({holder!r}, {len(self.components)} component(s), "
            f"{len(self.kij)} kij pair(s), {len(self.coefficients)} coefficient(s), "
            f"{len(self.models)} model(s))"
        )


#: There is no module-level card here and no accessor for one: a card is a value a
#: caller passes, so nothing in this library reads a global. `test_keycard_loader.py`
#: asserts that no `Keycard` binding exists at module level and that there is no
#: `current` or `clear` to reach for.


def load(path: str | Path) -> Keycard:
    """Read a keycard from a TOML file, returning it and setting nothing.

    ``card`` is passed to the calls that should read it, as
    ``azoth.eos.component("methane", card=card)``; a call without one reads the data
    this library ships.

    Raises:
        KeycardError: if the file cannot be read, does not parse, or does not
            describe a keycard this build can use.
    """
    where = Path(path)
    try:
        text = where.read_text(encoding="utf-8")
    except OSError as exc:
        raise KeycardError(str(where), f"cannot be read: {exc.strerror or exc}") from exc
    try:
        document = tomllib.loads(text)
    except tomllib.TOMLDecodeError as exc:
        raise KeycardError(str(where), f"does not parse as TOML: {exc}") from exc
    return use(document, path=where)


def use(document: Any, *, path: Path | None = None) -> Keycard:
    """Build a keycard from an already-parsed document.

    The same validation as :func:`load`, separated so a caller holding a mapping - a
    notebook, a test, a config service - does not have to write it to disk first.
    Like :func:`load` it returns the card and sets nothing.
    """
    where = str(path) if path is not None else "<mapping>"

    if not isinstance(document, Mapping):
        raise KeycardError(where, "must be a mapping with a schema_version key")

    version = document.get("schema_version")
    if version != SCHEMA_VERSION:
        raise KeycardError(
            where,
            f"declares schema_version {version!r}; this build reads {SCHEMA_VERSION}. "
            f"Refusing rather than guessing: version 1 had no `components`, "
            f"`coefficients` or `models`, so one of those in a version-1 file would "
            f"be silently dropped.",
        )

    known = {
        "schema_version",
        "keyholder",
        "components",
        "associations",
        "kij",
        "fluids",
        "fittings",
        "coefficients",
        "models",
    }
    unknown = sorted(set(document) - known)
    if unknown:
        raise KeycardError(
            where,
            f"has section(s) {unknown}, which this format does not define. Known: "
            f"{sorted(known)}. Refused rather than ignored: a section nothing reads "
            f"is data that looks in use and is not.",
        )

    kij, cpa_kij = _kij(document.get("kij"), where)
    return Keycard(
        keyholder=_keyholder(document, where),
        licence=_licence(document),
        components=_components(document.get("components"), where),
        associations=_associations(document.get("associations"), where),
        kij=kij,
        cpa_kij=cpa_kij,
        coefficients=_coefficients(document.get("coefficients"), where, "coefficients"),
        conventions=_conventions(document.get("coefficients")),
        models=_models(document.get("models"), where),
        fittings=tuple(document.get("fittings") or ()),
        fluids={name: tuple(rows) for name, rows in (document.get("fluids") or {}).items()},
        path=path,
    )


def _keyholder(document: Mapping[str, Any], where: str) -> str | None:
    holder = document.get("keyholder")
    if holder is None:
        return None
    if not isinstance(holder, Mapping) or not str(holder.get("name", "")).strip():
        raise KeycardError(where, "`keyholder` needs a non-empty `name`")
    return str(holder["name"]).strip()


def _licence(document: Mapping[str, Any]) -> str | None:
    """The licence reference, if the keyholder named one.

    No `where` parameter, unlike its neighbours: it cannot fail. An absent or
    malformed `keyholder` is caught by `_keyholder` and an absent `licence` is simply
    absent, so there is nothing here to attach a location to.
    """
    holder = document.get("keyholder")
    if not isinstance(holder, Mapping):
        return None
    licence = holder.get("licence")
    return str(licence) if licence is not None else None


def _quantity(value: Any, unit: Any, expected: str, where: str, field_name: str) -> Q:
    """A number and a unit, converted to the unit this build reads it in.

    The conversion is the point. A coefficient declared in the wrong unit is a
    factor with no symptom - a `Pc` in bar read as pascal is 1e-5, and the answer
    that comes out is wrong and entirely reasonable-looking - so anything that
    cannot be converted raises rather than being passed through.
    """
    if not isinstance(value, int | float) or isinstance(value, bool):
        raise KeycardError(where, f"`{field_name}` must be a number, got {value!r}")
    if unit not in UNIT_VOCABULARY:
        raise KeycardError(
            where,
            f"`{field_name}` declares the unit {unit!r}, which is not in the "
            f"vocabulary this format uses. Accepted: {list(UNIT_VOCABULARY)}. "
            f"`kelvin` is the one people reach for; the vocabulary spells it `K`.",
        )
    try:
        return ureg.Quantity(float(value), unit).to(expected)
    except pint.errors.PintError as exc:
        raise KeycardError(
            where,
            f"`{field_name}` is given as {value!r} {unit!r}, which cannot be read as "
            f"{expected}: {exc}",
        ) from exc


def _components(raw: Any, where: str) -> dict[str, dict[str, Q]]:
    if raw is None:
        return {}
    if not isinstance(raw, Mapping):
        raise KeycardError(where, "`components` must be a mapping of names to parameters")

    out: dict[str, dict[str, Q]] = {}
    for name, parameters in raw.items():
        field_name = f"components.{name}"
        if not isinstance(parameters, Mapping) or not parameters:
            raise KeycardError(where, f"`{field_name}` must carry at least one parameter")
        resolved: dict[str, Q] = {}
        for parameter, body in parameters.items():
            if parameter not in COMPONENT_PARAMETERS:
                raise KeycardError(
                    where,
                    f"`{field_name}.{parameter}` is not a parameter this build reads. "
                    f"Accepted: {sorted(COMPONENT_PARAMETERS)}. Refused rather than "
                    f"stored: a value nothing reads is data that looks in use and is not.",
                )
            if not isinstance(body, Mapping):
                raise KeycardError(
                    where, f"`{field_name}.{parameter}` must be a mapping with value and unit"
                )
            resolved[parameter] = _quantity(
                body.get("value"),
                body.get("unit"),
                COMPONENT_PARAMETERS[parameter],
                where,
                f"{field_name}.{parameter}",
            )
        present = [p for p in resolved if p.startswith("cp_")]
        if present and len(present) != 5:
            raise KeycardError(
                where,
                f"`{field_name}` needs all five of `cp_a` through `cp_e`, or none of them",
            )
        out[str(name).strip().lower()] = resolved
    return out


def _associations(raw: Any, where: str) -> dict[str, Association]:
    """The association each substance the card states one for, by lower-cased name.

    A section of its own rather than a key under ``components``, because a card may state a
    scheme for a substance whose critical constants it does not touch - and because the
    scheme is a name from a closed vocabulary where every other component parameter is a
    number in a unit, so the two cannot sit in one map without one of them being written as
    the other.
    """
    if raw is None:
        return {}
    if not isinstance(raw, Mapping):
        raise KeycardError(where, "`associations` must be a mapping of names to parameters")

    out: dict[str, Association] = {}
    for name, body in raw.items():
        field_name = f"associations.{name}"
        if not isinstance(body, Mapping) or not body:
            raise KeycardError(where, f"`{field_name}` must carry a scheme")
        scheme = body.get("scheme")
        if scheme not in ASSOCIATION_SCHEMES:
            raise KeycardError(
                where,
                f"`{field_name}.scheme` is {scheme!r}, which this build does not "
                f"implement. Accepted: {list(ASSOCIATION_SCHEMES)}. Refused rather than "
                f"accepted-and-ignored: the schemes are not interchangeable with each "
                f"other or with a site count - `2A` and `2B` both carry two sites and "
                f"bond differently - so a scheme evaluated as another is a different "
                f"fluid with no symptom.",
            )
        resolved: dict[str, Q] = {}
        for parameter, value in body.items():
            if parameter == "scheme":
                continue
            parameter_field = f"{field_name}.{parameter}"
            if parameter not in ASSOCIATION_PARAMETERS:
                raise KeycardError(
                    where,
                    f"`{parameter_field}` is not a parameter this build reads. Accepted: "
                    f"{sorted(ASSOCIATION_PARAMETERS)}. Refused rather than stored: a "
                    f"value nothing reads is data that looks in use and is not.",
                )
            if not isinstance(value, Mapping):
                raise KeycardError(
                    where, f"`{parameter_field}` must be a mapping with value and unit"
                )
            resolved[parameter] = _quantity(
                value.get("value"),
                value.get("unit"),
                ASSOCIATION_PARAMETERS[parameter],
                where,
                parameter_field,
            )
        out[str(name).strip().lower()] = Association(scheme=str(scheme), parameters=resolved)
    return out


def _kij(
    raw: Any, where: str
) -> tuple[dict[tuple[str, str], float], dict[tuple[str, str, str], float]]:
    """The interaction rows, as the classical parameter and the associating pair.

    One walk producing two maps, because it is one document with two columns: `value` is
    ``KIJPR`` and the two ``cpa_value_*`` are ``cpakij_SRK``/``cpakij_PR``. They are
    separate columns with separate fits - water/methanol is ``-0.153`` against
    ``-0.0789`` - so neither may stand in for the other, and the family is part of the
    second map's key rather than a pair carrying one value.
    """
    classical: dict[tuple[str, str], float] = {}
    associating: dict[tuple[str, str, str], float] = {}
    if raw is None:
        return classical, associating
    if not isinstance(raw, list):
        raise KeycardError(where, "`kij` must be a list of rows")

    for index, row in enumerate(raw):
        field_name = f"kij[{index}]"
        if not isinstance(row, Mapping):
            raise KeycardError(where, f"`{field_name}` must be a mapping")
        a, b, value = row.get("component_a"), row.get("component_b"), row.get("value")
        if not isinstance(a, str) or not a.strip():
            raise KeycardError(where, f"`{field_name}` needs a non-empty `component_a`")
        if not isinstance(b, str) or not b.strip():
            raise KeycardError(where, f"`{field_name}` needs a non-empty `component_b`")
        if not isinstance(value, int | float) or isinstance(value, bool):
            raise KeycardError(where, f"`{field_name}` needs a numeric `value`")
        if a.strip().lower() == b.strip().lower():
            raise KeycardError(
                where,
                f"`{field_name}` pairs {a!r} with itself. A component does not interact "
                f"with itself, and a self-pair would silently rescale its attraction.",
            )
        first, second = sorted((a.strip().lower(), b.strip().lower()))
        classical[(first, second)] = float(value)
        for key, family in (("cpa_value_srk", "srk"), ("cpa_value_pr", "pr")):
            stated = row.get(key)
            if stated is None:
                continue
            if not isinstance(stated, int | float) or isinstance(stated, bool):
                raise KeycardError(where, f"`{field_name}` needs a numeric `{key}`")
            associating[(first, second, family)] = float(stated)
    return classical, associating


def _coefficients(raw: Any, where: str, section: str) -> dict[str, dict[str, Q]]:
    if raw is None:
        return {}
    if not isinstance(raw, Mapping):
        raise KeycardError(where, f"`{section}` must map calculation ids to arguments")

    out: dict[str, dict[str, Q]] = {}
    for calc_id, arguments in raw.items():
        if not isinstance(arguments, Mapping) or not arguments:
            raise KeycardError(where, f"`{section}.{calc_id}` must carry at least one argument")
        resolved: dict[str, Q] = {}
        for name, body in arguments.items():
            field_name = f"{section}.{calc_id}.{name}"
            if not isinstance(body, Mapping):
                raise KeycardError(where, f"`{field_name}` must be a mapping with value and unit")
            unit = body.get("unit")
            if unit not in UNIT_VOCABULARY:
                raise KeycardError(
                    where,
                    f"`{field_name}` declares the unit {unit!r}, which is not in the "
                    f"vocabulary this format uses. Accepted: {list(UNIT_VOCABULARY)}. "
                    f"A unit is required rather than optional: a coefficient checked "
                    f"against the spec's is what catches a value in the wrong one, and "
                    f"a wrong unit is a factor with no symptom.",
                )
            try:
                resolved[name] = _coefficient_value(body["value"], unit, where, field_name)
            except KeyError as exc:
                raise KeycardError(where, f"`{field_name}` cannot be read: {exc}") from exc
            except (TypeError, ValueError, pint.errors.PintError) as exc:
                raise KeycardError(where, f"`{field_name}` cannot be read: {exc}") from exc
        out[str(calc_id)] = resolved
    return out


def _coefficient_value(raw: Any, unit: str, where: str, field_name: str) -> Any:
    """One coefficient: a number, a vector or a matrix, all in the unit beside it.

    A vector or a matrix comes back as a **list of quantities, one per entry**, which
    is the convention every declared vector and matrix input already follows -
    `eos.uniquac_activity_coefficients` converts its `aij` entry by entry, and
    `azoth.core.units.to_si_shaped` reads that shape. Wrapping the whole matrix in one
    quantity would need NumPy, and a coefficient is not worth a dependency.

    **A matrix is not a convenience.** NeqSim carries no UNIQUAC interaction table -
    `PhaseGEUniquac` is handed the SRK `kij` as its `aij` - so a caller who wants
    UNIQUAC has to supply that matrix or nobody can. The keycard is where a value the
    holder is accountable for belongs, and until this it could carry only scalars.

    The unit is one unit for the whole value: a matrix whose rows carried different
    units is a table, not a coefficient.
    """
    if isinstance(raw, bool) or not isinstance(raw, (int, float, list)):
        raise KeycardError(
            where,
            f"`{field_name}` is {type(raw).__name__} {raw!r}. A value is a number, a "
            f"non-empty vector or a non-empty rectangular matrix of numbers.",
        )
    if not isinstance(raw, list):
        return ureg.Quantity(raw, unit)
    if not raw:
        raise KeycardError(
            where,
            f"`{field_name}` is an empty list. A value is a number, a non-empty "
            f"vector or a non-empty rectangular matrix of numbers.",
        )
    if all(isinstance(entry, (int, float)) and not isinstance(entry, bool) for entry in raw):
        return [ureg.Quantity(entry, unit) for entry in raw]

    width: int | None = None
    rows: list[list[Any]] = []
    for i, row in enumerate(raw):
        if not isinstance(row, list) or not row:
            raise KeycardError(
                where,
                f"`{field_name}` row {i} is {row!r}. A value is a number, a non-empty "
                f"vector or a non-empty rectangular matrix of numbers.",
            )
        if width is None:
            width = len(row)
        elif len(row) != width:
            raise KeycardError(
                where,
                f"`{field_name}` row {i} has {len(row)} entries while row 0 has "
                f"{width}; a matrix is rectangular.",
            )
        if any(isinstance(entry, bool) or not isinstance(entry, (int, float)) for entry in row):
            raise KeycardError(
                where,
                f"`{field_name}` row {i} is {row!r}, which is not a row of numbers.",
            )
        rows.append([ureg.Quantity(entry, unit) for entry in row])
    return rows


def _conventions(raw: Any) -> dict[str, dict[str, str]]:
    """The `convention` string on each coefficient, where one was supplied.

    Cannot fail either, and so takes no `where`: a `convention` is free text whose
    only job is to travel with a value, and one that is the wrong type is simply not
    a convention rather than an error to report.
    """
    out: dict[str, dict[str, str]] = {}
    if not isinstance(raw, Mapping):
        return out
    for calc_id, arguments in raw.items():
        if not isinstance(arguments, Mapping):
            continue
        for name, body in arguments.items():
            if isinstance(body, Mapping) and isinstance(body.get("convention"), str):
                out.setdefault(str(calc_id), {})[str(name)] = str(body["convention"])
    return out


def _models(raw: Any, where: str) -> dict[str, Model]:
    if raw is None:
        return {}
    if not isinstance(raw, Mapping):
        raise KeycardError(where, "`models` must be a mapping of names to definitions")

    out: dict[str, Model] = {}
    for name, body in raw.items():
        field_name = f"models.{name}"
        if not isinstance(body, Mapping):
            raise KeycardError(where, f"`{field_name}` must be a mapping")

        for key, accepted in (
            ("kind", MODEL_KINDS),
            ("shape", MODEL_SHAPES),
            ("alpha", MODEL_ALPHAS),
            ("mixing_rule", MODEL_MIXING_RULES),
        ):
            chosen = body.get(key)
            if chosen not in accepted:
                raise KeycardError(
                    where,
                    f"`{field_name}.{key}` is {chosen!r}, which this build does not "
                    f"implement. Accepted: {list(accepted)}. Refused rather than "
                    f"accepted-and-ignored: a model declaring a variant and evaluated "
                    f"as another is a wrong answer with no symptom.",
                )

        components = body.get("components")
        if not isinstance(components, list) or not components:
            raise KeycardError(where, f"`{field_name}` needs a non-empty `components` list")

        parameters = body.get("alpha_parameters") or {}
        if parameters:
            raise KeycardError(
                where,
                f"`{field_name}.alpha_parameters` gives {sorted(parameters)}, but the "
                f"alpha function this build implements (`peng_robinson`) takes none. "
                f"A fitted alpha is a different equation, not a parameterisation of "
                f"this one.",
            )

        unknown = sorted(set(body) - MODEL_KEYS)
        if unknown:
            raise KeycardError(
                where,
                f"`{field_name}` carries {unknown}, which no reader here knows. "
                f"Accepted: {sorted(MODEL_KEYS)}. Refused rather than stored and "
                f"ignored: a key a caller can write and never see again is a promise "
                f"the format does not keep.",
            )

        out[str(name)] = Model(
            name=str(name),
            kind=str(body["kind"]),
            shape=str(body["shape"]),
            alpha=str(body["alpha"]),
            mixing_rule=str(body["mixing_rule"]),
            components=tuple(str(component) for component in components),
        )
    return out


def coefficient_value(calc_id: str, name: str, given: Any, *, card: Keycard | None = None) -> Any:
    """A coefficient for a calculation: what the caller passed, else the card's.

    **The card hands back the argument exactly as a caller would have passed it**, so
    nothing downstream can tell the two apart. For an input the spec declares
    dimensioned that is a quantity in the spec's own unit - `eos.uniquac_activity_coefficients`
    takes its `aij` as nested quantities, and a card supplying SI magnitudes there would
    be a second representation of the same input. For one it declares `dimensionless`
    it is a bare number, which is what a dimensionless input is everywhere in this
    library.

    An explicit argument always wins and the card is not consulted.

    The dimension is checked on the way: a coefficient declared in bar where the spec
    says pascal is caught here, because the conversion is what catches it. See
    :func:`azoth.core.units.to_si_shaped` for what a vector or matrix is made of.

    Raises:
        InvalidInputError: if the caller passed nothing and no card supplies it, or if
            the card's value carries the wrong dimensions. A missing coefficient is an
            error rather than a default: a plausible value nobody chose is a wrong answer
            with no symptom.
    """
    from azoth.core.errors import InvalidInputError
    from azoth.core.units import to_si_shaped

    if given is not None:
        return given

    supplied = card.coefficient(calc_id, name) if card is not None else None
    if supplied is None:
        holder = f"the keycard from {card.keyholder!r}" if card and card.keyholder else "no keycard"
        raise InvalidInputError(
            name,
            f"`{calc_id}` needs a value and none was passed; {holder} supplies one. "
            f"Either pass it, or add a `coefficients.{calc_id}.{name}` entry.",
        )

    declaration = _declaration(calc_id, name)
    unit = declaration["unit"]
    interval = bool(declaration.get("interval", False))
    # The conversion *is* the dimension check, so this is not a wasted call: it is what
    # refuses a `Tc` in pascals. Its result is used only where the spec is dimensionless,
    # because there the SI base magnitude is the number a caller passes.
    si = to_si_shaped(supplied, unit, name, interval=interval)
    if unit == "dimensionless":
        return si
    return _in_spec_unit(supplied, unit)


def _declaration(calc_id: str, name: str) -> Mapping[str, Any]:
    """One input's declaration, from whichever registry holds the id.

    **A coefficient may belong to a model, and this used to look only among the
    calculations.** `azoth._registry_gen` holds the calcs and `azoth._models_gen` the
    models; `eos.uniquac_activity_coefficients` is the latter, and asking the calc
    registry alone refused its id as unknown. Nothing noticed because the only caller
    was `hydraulics.orifice_flow`, which is a calc - the bug was waiting for the first
    model to want a card-supplied argument, which is exactly what the matrix work added.

    Raises:
        InvalidInputError: if the id is in neither registry, since there is then no
            declared unit to check a card's value against.
    """
    from azoth import _models_gen
    from azoth._registry_gen import BY_ID
    from azoth.core.errors import InvalidInputError

    document = BY_ID.get(calc_id) or _models_gen.MODEL_BY_ID.get(calc_id)
    if document is None:
        raise InvalidInputError(
            "coefficients",
            f"`{calc_id}` is in neither registry, so there is no spec declaring what "
            f"its `{name}` is and no unit to check a card's value against",
        )
    declaration: Mapping[str, Any] = document["inputs"][name]
    return declaration


def _in_spec_unit(value: Any, unit: str) -> Any:
    """A card's quantity, or nested quantities, restated in the spec's own unit.

    Recursive for the same reason `to_si_shaped` is: a vector and a matrix are lists of
    quantities, one per entry.
    """
    from azoth.core.units import unit_for

    if isinstance(value, list):
        return [_in_spec_unit(entry, unit) for entry in value]
    return value.to(unit_for(unit))


__all__ = [
    "ASSOCIATION_PARAMETERS",
    "ASSOCIATION_SCHEMES",
    "COMPONENT_PARAMETERS",
    "MODEL_ALPHAS",
    "MODEL_KINDS",
    "MODEL_MIXING_RULES",
    "MODEL_SHAPES",
    "SCHEMA_VERSION",
    "Association",
    "Keycard",
    "Model",
    "coefficient_value",
    "load",
    "use",
]
