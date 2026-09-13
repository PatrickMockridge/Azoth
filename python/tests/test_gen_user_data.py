"""The user data generator: does what it writes match what the loaders read.

``tools/gen_user_data.py`` compiles a checked ``keycard.yaml`` into the CSVs
under ``data/``. The dangerous property of a generator like that is not that it
crashes - it is that it writes something *nearly* right: a column renamed, a value
reformatted, a header line missing. Every one of those produces a file that looks
like data, loads without complaint in one language, and silently disagrees with
the other.

So these tests check three things, and none of them is "the tool runs":

* **The format is the committed format.** The shipped fittings registry is
  regenerated row for row from its own parsed contents and must come back
  byte-identical, so the generator reproduces the file the loaders already read
  rather than a new dialect of it.
* **The columns are the loaders' contract.** The generator's column tuples are
  compared against the header rows of the committed files, which are what both
  implementations key off.
* **Numbers keep their text.** ``0.0000200`` and ``2e-05`` are the same float and
  not the same statement about precision, so a value the user wrote is carried
  through as written.
"""

from __future__ import annotations

import csv
import importlib.util
import sys
from pathlib import Path
from typing import Any

import pytest
import yaml

REPO_ROOT = Path(__file__).resolve().parents[2]
GENERATOR = REPO_ROOT / "tools" / "gen_user_data.py"
FITTINGS_CSV = REPO_ROOT / "data" / "fittings" / "crane_k_factors.csv"


def _load_generator() -> Any:
    """Import `tools/gen_user_data.py` by path.

    `tools/` is not on the import path and is not part of the installed package -
    it is build tooling, alongside `check_user_data.py` and `gen_registry.py`,
    which the other tool tests reach the same way.
    """
    sys.path.insert(0, str(REPO_ROOT / "tools"))
    spec = importlib.util.spec_from_file_location("gen_user_data", GENERATOR)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


gen = _load_generator()


def data_rows(path: Path) -> list[dict[str, str]]:
    """The data rows of a shipped CSV, with the comment banner stripped.

    Exactly what both loaders do - filter `#` lines, then parse - so a test using
    this reads the file the same way the implementations do.
    """
    body = [
        line
        for line in path.read_text(encoding="utf-8").splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ]
    return list(csv.DictReader(body))


def test_the_generator_reproduces_the_shipped_registry_byte_for_byte() -> None:
    """Row for row, the committed file comes back out of the generator unchanged.

    The strongest available check that the generator writes the format the loaders
    already read, rather than a plausible new one. It would catch a renamed column,
    a reordered field, a changed line terminator, a lost header, or a number
    rendered differently - none of which any single-language test would notice,
    because a CSV that both sides parse consistently can still parse *wrongly*
    into the same wrong values.

    What it does not prove: that the banner is right, since the banner is what the
    generator emits and this compares against a file the generator already wrote.
    The banner is checked separately, by `test_the_banner_stops_claiming_...`.
    """
    rows = data_rows(FITTINGS_CSV)
    assert rows, "the shipped registry has no data rows; the reader is wrong, not the file"
    # `id` is the user-facing name; `fitting_id` is the column. The generator is
    # what bridges them, which is one of the things this test pins.
    as_user_rows = [{**row, "id": row["fitting_id"]} for row in rows]

    assert gen.render_fittings(as_user_rows) == FITTINGS_CSV.read_text(encoding="utf-8"), (
        "the generator no longer reproduces data/fittings/crane_k_factors.csv. "
        "Either the file was edited by hand, or the generator's format drifted."
    )


def test_the_generator_writes_the_columns_the_loaders_read() -> None:
    """The column tuples are the cross-language contract, so they are compared.

    Rust's `Fitting` and Python's `PropertyPoint` both key off these exact header
    names. A generator with its own idea of them produces a file the loaders parse
    into empty or misaligned fields - and `csv.DictReader` does not fail on a
    missing column, it returns `None`, which becomes an empty string and then a
    citation that says nothing.
    """
    fittings_header = next(iter(data_rows(FITTINGS_CSV)))
    assert tuple(fittings_header) == gen.FITTING_COLUMNS, (
        "the generator's fitting columns differ from the shipped file's header"
    )

    for fluid in ("water", "air"):
        header = next(iter(data_rows(REPO_ROOT / "data" / "fluids" / f"{fluid}.csv")))
        assert tuple(header) == gen.FLUID_COLUMNS, f"{fluid}.csv header differs"


def test_a_number_keeps_the_text_the_user_wrote() -> None:
    """`0.0000200` survives as `0.0000200`, not as `2e-05`.

    The same float, and a different statement: the first says three significant
    figures, which is information an engineer reading a property table needs and a
    tool has no business discarding. This is why the loader keeps scalar text.
    """
    document = gen.yaml.load("t: 0.0000200\nn: 30.0\nb: 2e-05\n", Loader=gen.TextLoader)
    assert document == {"t": "0.0000200", "n": "30.0", "b": "2e-05"}

    rendered = gen.render_fluid(
        "water",
        [
            {
                "temperature_c": "0.0",
                "density_kg_m3": "999.8",
                "dynamic_viscosity_pa_s": "0.0000200",
                "citation": "x",
                "verify_status": "unverified",
            }
        ],
    )
    assert "0.0000200" in rendered, "the significant figures were reformatted away"
    assert "2e-05" not in rendered


def test_the_banner_stops_claiming_every_value_is_a_placeholder() -> None:
    """The loudest line in the file is a claim about the file, so it is derived.

    `ERROR: EVERY VALUE ... IS A DUMMY` is true of the shipped registry and false
    the moment a user supplies real numbers. A banner that stayed put would be
    worse than absent - it would teach readers to skip the one line that matters
    most, which is the failure mode this whole mechanism is built to avoid.
    """
    dummy = {
        "id": "x",
        "family": "bend",
        "name": "x",
        "n_ld": "1",
        "f_t_basis": "f_t",
        "citation": "DUMMY",
        "verify_status": "estimated_dummy",
        "source_ref": "",
        "source_locator": "",
    }
    assert "EVERY VALUE IN THIS FILE IS AN ESTIMATED DUMMY" in gen.render_fittings([dummy])

    real = {**dummy, "verify_status": "verified", "citation": "read it", "source_ref": "https://x"}
    rendered = gen.render_fittings([real])
    assert "EVERY VALUE IN THIS FILE IS AN ESTIMATED DUMMY" not in rendered, (
        "the placeholder banner survived a file that is not all placeholders"
    )
    assert "NOT EVERY VALUE IN THIS FILE IS A PLACEHOLDER" in rendered
    # And it must tell the user the thing that actually matters to them.
    assert "do not commit" in rendered.lower()


def test_the_banner_does_not_claim_the_tool_wrote_a_file_it_did_not() -> None:
    """The shipped water table is real data that never went near this tool.

    A banner saying "generated from keycard.yaml" would be false about it, and
    the water and air tables are legitimately committed - they are published facts,
    not a licensed table. So the non-placeholder banner has to be true of both
    cases, which means it names the restriction conditionally rather than
    asserting it.
    """
    water = data_rows(REPO_ROOT / "data" / "fluids" / "water.csv")
    rendered = gen.render_fluid("water", water)

    assert "GENERATED FILE" not in rendered
    assert "IF YOU RAN" in rendered
    # The rows themselves survive, which is the point of rendering them.
    assert water[0]["citation"] in rendered
    assert water[0]["verify_status"] in rendered


def test_a_fluid_that_no_implementation_can_read_is_refused() -> None:
    """`data/fluids/example_fluid.csv` would be a file nothing loads.

    Both languages hardcode the fluid list - Rust through `include_str!` and a
    `match`, Python through `_BUILTINS` - so generating a table for an unregistered
    name produces a file that looks like data in use and is read by nothing. That
    is the failure this project is organised against, so it is refused with the
    list of places to register it rather than written out.
    """
    document = {
        "fittings": None,
        "fluids": {"example_fluid": [{"temperature_c": "0"}, {"temperature_c": "1"}]},
    }
    with pytest.raises(gen.UnregisteredFluid) as excinfo:
        gen.plan({k: v for k, v in document.items() if v is not None})
    assert excinfo.value.name == "example_fluid"
    assert set(excinfo.value.known) >= {"water", "air"}, "built-ins not discovered"


def test_the_generator_agrees_with_the_checker_about_the_template() -> None:
    """The template generates, and does so through the checker's own rules.

    A file that passes `check_user_data.py` and is refused by the generator would
    mean two definitions of what a valid row is - the drift this project removes
    wherever it appears. The template is the file both are documented against, so
    it is the one to run them both over.
    """
    template = REPO_ROOT / "keycard.example.yaml"
    document = gen.yaml.load(template.read_text(encoding="utf-8"), Loader=gen.TextLoader)
    report = gen.validate(document, template)
    assert not report.errors, (
        f"the template passes tools/check_user_data.py but the generator rejects it: "
        f"{report.errors}. Two definitions of a valid row, which is the drift this "
        f"project removes wherever it appears."
    )
    assert report.by_status, "no rows were counted, so validation did not run"
    # `unstated` is the expected label: a keycard row carries no verification
    # status, which is the point of the change.


def test_the_template_can_be_generated_once_the_fluids_are_registered() -> None:
    """With the fluids section removed, the template compiles end to end.

    The template deliberately ships a fluid that is not built in, so that it can
    show the shape without pretending the repository carries a third table. Taking
    that section out leaves a document the tool can write, which is what makes the
    rest of this file's claims reachable rather than hypothetical.
    """
    template = REPO_ROOT / "keycard.example.yaml"
    document = gen.yaml.load(template.read_text(encoding="utf-8"), Loader=gen.TextLoader)
    document.pop("fluids")

    outputs = gen.plan(document)
    assert len(outputs) == 1
    relative, text = outputs[0]
    assert relative == Path("data/fittings/crane_k_factors.csv")

    # It is a valid registry: strip the banner and parse it back.
    body = [
        line for line in text.splitlines() if line.strip() and not line.lstrip().startswith("#")
    ]
    rows = list(csv.DictReader(body))
    assert rows and rows[0]["fitting_id"] == "90_elbow"
    assert rows[0]["n_ld"] == "30.0", "the template's own text was reformatted"


def test_the_shipped_registry_is_still_all_placeholders() -> None:
    """The generator must not have quietly upgraded the repository's own data.

    Nothing in this file writes to `data/` - every test works on rendered text -
    but the failure this guards is the one the generator makes *possible*: run it
    against a real data file and the placeholders are gone. A test asserting they
    are still there is cheap, and it fails at the moment that matters.
    """
    statuses = {row["verify_status"] for row in data_rows(FITTINGS_CSV)}
    assert statuses == {"estimated_dummy"}, (
        f"data/fittings/crane_k_factors.csv is no longer all placeholders ({statuses}). "
        f"If you generated licensed data into it, restore it with "
        f"`git checkout -- data/` - see docs/src/copyright.md."
    )

    for fluid in ("water", "air"):
        statuses = {
            row["verify_status"]
            for row in data_rows(REPO_ROOT / "data" / "fluids" / f"{fluid}.csv")
        }
        assert statuses == {"unverified"}, f"{fluid}.csv statuses changed: {statuses}"


def test_an_unknown_schema_version_is_refused(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    """A newer document shape is refused rather than guessed at.

    `TextLoader` carries every scalar as text, so a version arrives as a string
    rather than a number. Reading it as a float or an int would accept `1.5` or
    reject `"2"`, and either would be a silent decision about a format this tool
    does not know - so it is compared as text, and anything else is refused.
    """
    path = tmp_path / "keycard.yaml"
    path.write_text(
        yaml.safe_dump(
            {
                "schema_version": 99,
                "fittings": [
                    {
                        "id": "x",
                        "family": "bend",
                        "name": "x",
                        "n_ld": 1.0,
                        "f_t_basis": "f_t",
                        "citation": "DUMMY",
                    }
                ],
            }
        ),
        encoding="utf-8",
    )
    assert gen.main([str(path)]) == 1
    assert "schema_version" in capsys.readouterr().err


def test_a_document_that_fails_its_checks_writes_nothing(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    """Validation comes before writing, and a failure stops the write.

    Checking and writing are separate steps so that they cannot become one step
    that does both when one of them fails - which is the same reason
    `check_user_data.py` does not write. The failure that motivates it: a partial
    or unchecked file landing in `data/` and being embedded by the next build.
    """
    path = tmp_path / "keycard.yaml"
    path.write_text(
        yaml.safe_dump(
            {
                "schema_version": 2,
                "fittings": [
                    {
                        "id": "x",
                        "family": "bend",
                        "name": "x",
                        "n_ld": -1.0,  # non-positive, which the checker refuses
                        "f_t_basis": "f_t",
                        "citation": "DUMMY",
                    }
                ],
            }
        ),
        encoding="utf-8",
    )
    assert gen.main([str(path)]) == 1
    err = capsys.readouterr().err
    assert "non-positive equivalent length" in err
    assert "Nothing was written" in err
