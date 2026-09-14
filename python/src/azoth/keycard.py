"""The keycard: the file a user supplies to extend or override what azoth ships.

Sections: ``keyholder``, ``components``, ``kij``, ``fluids``, ``fittings``,
``coefficients`` and ``models``. A name the databank already has replaces its
critical constants; a name it does not have adds one.

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
class Keycard:
    """A loaded keycard. Immutable, and the whole of what a file said."""

    keyholder: str | None
    licence: str | None
    components: Mapping[str, Mapping[str, Q]] = field(default_factory=dict)
    kij: Mapping[tuple[str, str], float] = field(default_factory=dict)
    coefficients: Mapping[str, Mapping[str, Q]] = field(default_factory=dict)
    conventions: Mapping[str, Mapping[str, str]] = field(default_factory=dict)
    models: Mapping[str, Model] = field(default_factory=dict)
    fittings: tuple[Mapping[str, Any], ...] = ()
    fluids: Mapping[str, tuple[Mapping[str, Any], ...]] = field(default_factory=dict)
    path: Path | None = None

    def component(self, name: str) -> Mapping[str, Q] | None:
        """The parameters this keycard gives a substance, or ``None``."""
        return self.components.get(name.strip().lower())

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

    return Keycard(
        keyholder=_keyholder(document, where),
        licence=_licence(document),
        components=_components(document.get("components"), where),
        kij=_kij(document.get("kij"), where),
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


def _kij(raw: Any, where: str) -> dict[tuple[str, str], float]:
    if raw is None:
        return {}
    if not isinstance(raw, list):
        raise KeycardError(where, "`kij` must be a list of rows")

    out: dict[tuple[str, str], float] = {}
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
        out[(first, second)] = float(value)
    return out


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
                resolved[name] = ureg.Quantity(float(body["value"]), unit)
            except (KeyError, TypeError, ValueError, pint.errors.PintError) as exc:
                raise KeycardError(where, f"`{field_name}` cannot be read: {exc}") from exc
        out[str(calc_id)] = resolved
    return out


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

    An explicit argument always wins and the card is not consulted. The value is
    converted against the *spec's* declared unit for that input, so a coefficient
    declared in bar and read as pascal is caught here.

    Raises:
        InvalidInputError: if the caller passed nothing and no card supplies it.
            A missing coefficient is an error rather than a default: a plausible
            value nobody chose is a wrong answer with no symptom.
    """
    from azoth._registry_gen import spec  # local: the registry is generated, not core
    from azoth.core.errors import InvalidInputError
    from azoth.core.units import input_to_si

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
    return input_to_si(spec(calc_id), name, supplied)


__all__ = [
    "COMPONENT_PARAMETERS",
    "MODEL_ALPHAS",
    "MODEL_KINDS",
    "MODEL_MIXING_RULES",
    "MODEL_SHAPES",
    "SCHEMA_VERSION",
    "Keycard",
    "Model",
    "coefficient_value",
    "load",
    "use",
]
