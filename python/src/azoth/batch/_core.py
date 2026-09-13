"""Running one calculation over arrays, driven by the spec.

# What this file is

The whole batch layer in one place: validate the inputs against the calc's spec,
convert them once, hand them to whichever backend is selected, and hand the columns to
a per-calc builder that names them.

It is spec-driven for the reason the rest of the project is. The input names and units
come from the spec, so a batch call cannot disagree with the scalar call about what an
argument is called or what it is in, and adding a calculation does not mean restating
either.

# Units: the batch API speaks in the spec's canonical units

`rho=[998.0, 1000.0]` means 998 and 1000 kg/m**3, because that is what the spec declares
for `rho`. Inputs are plain numbers in those units and outputs come back as SI base
magnitudes with a unit map.

That is a deliberate departure from the scalar API, which takes and returns `pint`
quantities. Converting a sequence of quantities means one `pint` operation per element,
which is the same order of cost as the per-element boundary crossing the batch API
exists to remove. So the conversion here is **one factor applied to a whole array**,
taken from `pint` rather than restated - see `_si_per_canonical`.

The cost is that a caller who passes metres where the spec says millimetres gets a
thousand-fold error with nothing to catch it. That is the same error the scalar API's
`UnitMismatchError` exists to prevent, and this is the one place it is genuinely given
up. It is stated here, in the result's docstring, and on the batch page.
"""

from __future__ import annotations

import dataclasses
import importlib
from array import array
from collections.abc import Callable, Iterable, Mapping, Sequence
from enum import Enum
from functools import cache
from typing import Any, get_args, get_origin, get_type_hints

from azoth._dispatch import select
from azoth._registry_gen import BY_ID
from azoth.core.errors import InvalidInputError
from azoth.core.result import FlowRegime
from azoth.core.units import unit_for, ureg
from azoth.core.warnings import Warning

#: A calc's columns, its unit map, and its per-element warnings.
type Columns = dict[str, Any]
type Builder = Callable[[Columns, dict[str, str], tuple[tuple[Warning, ...], ...]], Any]


@cache
def _si_per_canonical(spec_unit: str) -> float:
    """How many SI base units one of a spec's canonical units is worth.

    Taken from `pint` rather than written down, so `mm` is 0.001 because the units
    library says so and not because someone typed it. Restating what the units library
    knows is what produced the `mm` bug this repository already had once.

    Cached because it is looked up once per input per batch rather than once per
    element, and the whole point of this module is that the per-element work is a
    multiply.
    """
    return float(ureg.Quantity(1.0, unit_for(spec_unit)).to_base_units().magnitude)


@cache
def _base_unit(spec_unit: str) -> str:
    """The SI base unit one of a spec's canonical units reduces to.

    The inverse half of `_si_per_canonical`, and needed because the reference backend has
    to rebuild `pint` quantities from the SI magnitudes `_to_si` produced. Rebuilding them
    with the *declared* unit instead is the double conversion that made a 50 mm orifice
    come out at 7.7e-09 m**3/s against the extension's 7.7e-03 - a factor of 10**6, through
    the square on `d` - on the Python path only, with every test still green.

    Both spellings come from `pint` rather than being written down, so they cannot disagree
    about what `mm` is: `_si_per_canonical("mm") == 0.001` and `_base_unit("mm") == "meter"`.
    """
    return str(ureg.Quantity(1.0, unit_for(spec_unit)).to_base_units().units)


def _to_si(values: Sequence[float], spec_unit: str) -> list[float]:
    """Convert a whole input array to SI base magnitudes.

    The multiply is skipped when the factor is exactly 1, which is the common case: every
    unit in the vocabulary except `mm` is its own SI base unit.
    """
    factor = _si_per_canonical(spec_unit)
    if factor == 1.0:
        return list(values)
    return [value * factor for value in values]


def _check_inputs(calc_id: str, spec: Mapping[str, Any], arrays: Mapping[str, Any]) -> int:
    """Validate the supplied arrays, returning the batch length.

    Raises:
        InvalidInputError: for a missing input, an unknown one, or arrays of differing
            lengths. All three are caller errors with no sensible default: a missing
            array cannot be guessed, an unknown name is a typo, and differing lengths
            have no element-wise meaning.

    Only the third is reachable through a public wrapper, and that is the intended order:
    the wrappers declare every input as a keyword-only parameter, so a missing array or a
    misspelt name is a `TypeError` at the call site - and a `mypy` error before that -
    rather than something this function has to be the one to notice. The first two are
    kept because `run` is reachable directly, and because a check that cannot fire is
    still the right shape for the day a wrapper is generated instead of written.
    """
    declared = spec["inputs"]
    required = {name for name, d in declared.items() if not d.get("optional", False)}

    missing = sorted(required - set(arrays))
    if missing:
        raise InvalidInputError(
            calc_id,
            f"batch call is missing required input(s) {missing}; the spec "
            f"declares {sorted(declared)}",
        )
    unknown = sorted(set(arrays) - set(declared))
    if unknown:
        raise InvalidInputError(
            calc_id, f"batch call passes {unknown}, which the spec does not declare"
        )

    lengths = {name: len(values) for name, values in arrays.items()}
    distinct = set(lengths.values())
    if len(distinct) > 1:
        raise InvalidInputError(
            calc_id,
            f"batch inputs have differing lengths ({lengths}); every array must be the "
            f"same length, because this API does not broadcast - a length-one array "
            f"where a length-N array was meant produces plausible answers rather than "
            f"an error",
        )
    return distinct.pop() if distinct else 0


def batchable() -> tuple[str, ...]:
    """Every calc id the batch API can run, sorted.

    The complement of what `_require_batchable` rejects, exposed so a test can assert
    that both backends carry an arm for every one of these and nothing else. That
    assertion is what makes "a new calc gets a batch variant" fail loudly instead of
    silently: without it, adding a calc leaves `batch_run` raising `NotImplementedError`
    at call time, which is a state nothing in the build would notice.
    """
    return tuple(
        calc_id
        for calc_id, spec in sorted(BY_ID.items())
        if not any(declared.get("unit") is None for declared in spec["inputs"].values())
    )


def _column_kinds(spec: Mapping[str, Any]) -> dict[str, str]:
    """Which output fields a batch call carries, and which kind each one is.

    Derived from the result dataclass's own type, not from the values it happens to hold.
    That distinction is the whole point: `darcy_weisbach`'s `re` is `float | None`, and
    sniffing the first element would make its column *numeric* when viscosity was supplied
    and *labels* when it was not. A field whose kind depends on which arguments the caller
    passed is not a field a caller can write against, and the extension sends a `NaN`
    column either way, so the two backends would also have disagreed.

    The rules, in order:

    * drop `warnings`, which is not an output;
    * drop anything the caller supplied - a friction factor comes back as an argument, not
      as an answer;
    * drop anything that is a per-element *structure* (`crane_k_factors`'s `components`, a
      tuple of per-fitting breakdowns). A batch column is one sequence of scalars; a
      sequence of sequences has no column shape, and flattening it would need a length
      column the caller did not ask for that differs per element. That calc has no batch
      form for a separate reason - see `_require_batchable` - so this is a rule kept honest
      by a test rather than one doing work today;
    * everything whose type is an enum is a **label** column; everything else is
      **numeric**.

    Derived rather than listed per calc so a calc gaining an output gains it here too, and
    so this side cannot disagree with the Rust arm about which columns exist: if they do
    disagree, the cross-language batch test says so.
    """
    from azoth.core.result import RESULT_TYPES

    result_type: Any = RESULT_TYPES[spec["id"]]
    hints = get_type_hints(result_type)
    inputs = set(spec["inputs"])
    kinds: dict[str, str] = {}
    for field in dataclasses.fields(result_type):
        hint = hints[field.name]
        if field.name == "warnings" or field.name in inputs:
            continue
        if get_origin(hint) in (tuple, list):
            continue
        kinds[field.name] = "labels" if _is_enum_hint(hint) else "values"
    return kinds


def _is_enum_hint(hint: object) -> bool:
    """Whether a field's type is (or is an optional) enum.

    Any `Enum`, not one of them by name. This checked `FlowRegime` specifically until
    `eos.pr_z_factor`'s `root_structure` became the first enum output outside the
    hydraulics namespace - at which point a field that was correctly a label column
    in the scalar API came back as a numeric one, and the batch arm's label went into
    it as a string a caller could not read as a number. The paragraph above already
    stated the general rule; the code now implements it.

    `FlowRegime | None` is still the motivating case: an absent regime is not any of
    the three real ones, so it survives as `None` in the label column. So does an
    absent `RootStructure`.
    """
    return any(
        isinstance(argument, type) and issubclass(argument, Enum)
        for argument in (hint, *get_args(hint))
        if argument is not None
    )


def _require_batchable(calc_id: str, spec: Mapping[str, Any]) -> None:
    """Reject a calc whose inputs are not all numbers.

    `crane_k_factors` takes `fittings`: a list of registry ids, with no unit, because a
    fitting is a category and not a magnitude. There is no batch shape for it - N
    elements each carrying their own list of names is a sequence of sequences, and one
    shared list would compute the same answer N times, which is not what a batch API is
    for. So it is excluded rather than half-supported, and a test asserts the excluded set
    is exactly this, so the exclusion stays a decision instead of becoming an oversight.

    Derived from the spec, which is what makes it mechanical: a declared input with no
    `unit` is not a magnitude. Every dimensionless input does declare `dimensionless`.
    """
    categorical = sorted(
        name for name, declared in spec["inputs"].items() if declared.get("unit") is None
    )
    if categorical:
        raise InvalidInputError(
            calc_id,
            f"{calc_id} takes {categorical}, which the spec declares without a unit - a "
            f"category, not a magnitude. The batch API carries numeric columns only, so "
            f"this calculation has no batch form. Call it once per element instead.",
        )


def _label(value: Any) -> str | None:
    """A value as an enum label, or `None` if it is absent."""
    if value is None:
        return None
    if isinstance(value, FlowRegime):
        return str(value.value)
    return str(value)


def _as_column(values: Iterable[Any], kind: str) -> Any:
    """One output's values, as either a numeric array or a label tuple.

    `kind` comes from the field's declared type (`_column_kinds`), never from the values,
    so a column's kind cannot change because a caller omitted an optional argument.
    """
    if kind == "labels":
        return tuple(_label(value) for value in values)
    # A dimensioned value arrives as a `pint` quantity from the reference and as a plain
    # magnitude from the extension, so both are reduced to the magnitude here. *Base*
    # units, because that is what the extension sends: `uom`'s `.value` is always the SI
    # base magnitude, so reducing only one side would put the two a factor away for any
    # unit that is not already a base one.
    return array("d", (_magnitude(value) for value in values))


def _magnitude(value: Any) -> float:
    """One output value as an SI base magnitude.

    `NaN` for an absent one, which is what the extension sends for `darcy_weisbach`'s `re`
    when viscosity was omitted. `None` would have to become a label column, and a numeric
    field that turns into a label column because an optional argument was left out is a
    field no caller can write against - and would have made the two backends disagree.
    `NaN` rather than zero because a Reynolds number of zero is a physical claim, and the
    input that would let this calc make one was withheld.
    """
    if value is None:
        return float("nan")
    if hasattr(value, "magnitude"):
        return float(value.to_base_units().magnitude)
    return float(value)


def _run_reference(calc_id: str, spec: Mapping[str, Any], si: Mapping[str, list[float]]) -> Columns:
    """The pure-Python path: the scalar reference, once per element.

    Rebuilds quantities for the dimensioned inputs, which costs a `pint` operation per
    element - accepted rather than avoided, because this is the reference implementation
    and its job is to be the thing the fast path is checked against. A caller who wants
    speed selects the Rust backend.
    """
    namespace, _, function_name = calc_id.rpartition(".")
    module = importlib.import_module(f"azoth.{namespace}.reference.{function_name}")
    function = getattr(module, function_name)

    declared = spec["inputs"]
    kinds = _column_kinds(spec)
    collected: dict[str, list[Any]] = {name: [] for name in kinds}
    warnings: list[tuple[Warning, ...]] = []

    for index in range(len(next(iter(si.values()), []))):
        kwargs: dict[str, Any] = {}
        for name, declaration in declared.items():
            if name not in si:
                continue
            unit = declaration.get("unit")
            # Rebuilt in the *base* unit, because `run` already put `si` in base
            # magnitudes. Using the declared unit here is a second conversion, and it is
            # the one that made a 50 mm orifice disagree with the extension by 10**6 on
            # this backend alone.
            kwargs[name] = (
                ureg.Quantity(si[name][index], _base_unit(unit))
                if unit is not None and unit != "dimensionless"
                else si[name][index]
            )
        result = function(**kwargs)
        for field in kinds:
            collected[field].append(getattr(result, field))
        warnings.append(tuple(result.warnings))

    columns: Columns = {}
    units: dict[str, str] = {}
    for field, kind in kinds.items():
        columns[field] = _as_column(collected[field], kind)
        if kind == "values":
            units[field] = _canonical_unit(spec, field)
    return {"columns": columns, "units": units, "warnings": tuple(warnings)}


def _canonical_unit(spec: Mapping[str, Any], field: str) -> str:
    """The spec's declared unit for an output field.

    Falls back to the SI base unit when the spec declares none, which happens for an
    output the spec marks as a derived diagnostic rather than a declared output.
    """
    declaration = spec["outputs"].get(field)
    if declaration is None or declaration.get("unit") is None:
        return "dimensionless"
    return str(declaration["unit"])


def _run_rust(calc_id: str, si: Mapping[str, list[float]]) -> Columns:
    """The Rust path: one call, N evaluations, arrays in and arrays out."""
    from azoth import _core
    from azoth._rust_bridge import _warnings

    result = _core.batch_run(calc_id, {name: list(values) for name, values in si.items()})
    columns: Columns = {}
    units: dict[str, str] = {}
    for column in result.columns:
        if column.labels is not None:
            columns[column.name] = tuple(column.labels)
        else:
            columns[column.name] = array("d", column.values or [])
            units[column.name] = column.unit
    # The transported warnings become real ones through the same adapter the scalar
    # bridge uses, rather than a second copy of the same three-field construction. Two
    # copies of a cross-language conversion is one more place for the two to disagree,
    # and this one is already tested by every scalar path.
    warnings = tuple(_warnings(element) for element in result.warnings)
    return {"columns": columns, "units": units, "warnings": warnings}


def run(calc_id: str, arrays: Mapping[str, Any], build: Builder) -> Any:
    """Run one calculation over arrays of inputs.

    The public entry is a per-calc function in :mod:`azoth.batch.hydraulics` or
    :mod:`azoth.batch.thermal`, which owns the argument names and the result type. This
    function is the machinery they share.

    Raises:
        InvalidInputError: for a missing, unknown or mismatched input, as
            `_check_inputs` describes.
        OutOfRangeError: from whichever element violates a hard bound. **Fail-fast** -
            the whole call raises rather than returning the elements that worked with a
            silent gap where one did not. An exception from the extension is converted to
            the same class the reference raises, so a caller cannot tell which backend
            answered by what they caught.
    """
    spec = BY_ID[calc_id]
    _require_batchable(calc_id, spec)
    length = _check_inputs(calc_id, spec, arrays)

    declared = spec["inputs"]
    si: dict[str, list[float]] = {}
    for name, values in arrays.items():
        declaration = declared[name]
        unit = declaration.get("unit")
        if unit is None or unit == "dimensionless":
            si[name] = list(values)
        else:
            si[name] = _to_si(values, unit)

    raw = _run_rust(calc_id, si) if select() == "rust" else _run_reference(calc_id, spec, si)
    if length != len(raw["warnings"]):  # pragma: no cover - defensive
        raise InvalidInputError(
            calc_id,
            f"asked for {length} element(s) but the backend returned {len(raw['warnings'])}",
        )
    return build(raw["columns"], raw["units"], raw["warnings"])


def sequence(value: Sequence[float], field: str) -> list[float]:
    """Coerce one input argument to a list of numbers, rejecting a bare scalar.

    A scalar where an array belongs is rejected rather than broadcast: a length-one
    array where a length-N array was meant would produce plausible answers rather than an
    error, and this API has no way to tell that case from the intended one.
    """
    if isinstance(value, int | float) or hasattr(value, "magnitude"):
        raise InvalidInputError(
            field,
            f"expected a sequence of numbers for {field!r}, got a single value. The "
            f"batch API does not broadcast; pass a list, even for one element.",
        )
    return [float(v) for v in value]
