# Worked notebooks

Five Jupyter notebooks using the published `azoth-engine` package. They read no repository file:
each one installs, imports `azoth`, and runs. Their outputs are committed so GitHub renders the
tables and figures, and CI re-executes them so they cannot stop running quietly.

| notebook | what it shows |
|---|---|
| [`01-simple-calculations.ipynb`](01-simple-calculations.ipynb) | pint quantities, Reynolds number, three friction factors, Darcy-Weisbach, pump power, conduction, a refused `degC` difference, and the provenance block. |
| [`02-intermediate-thermodynamics.ipynb`](02-intermediate-thermodynamics.ipynb) | a Peng-Robinson mixture, TP flash, vapour fraction across temperature, the P-T phase envelope, saturation pressure, enthalpy and transport properties at the flash roots. |
| [`03-advanced-calculations.ipynb`](03-advanced-calculations.ipynb) | an LMTD correction factor, packed-column sizing and its operating window, a liquid property table from a flash, and a check that did not run. |
| [`04-process-simulation.ipynb`](04-process-simulation.ipynb) | the unit-operation palette, a flowsheet document opened with the embedded palette, its stream table, an active recycle, an edit-and-run, and a unit operation called directly. |
| [`05-data-science.ipynb`](05-data-science.ipynb) | the batch API over arrays, pandas tables, a compressibility sweep, a conduction sweep, and ISO 6976 over three gas blends. |

## Run them

Python 3.12 or later:

```bash
python -m venv .venv
. .venv/bin/activate
pip install azoth-engine jupyterlab matplotlib pandas
jupyter lab
```

The `notebook` extra in `pyproject.toml` names the same tools; it ships with the next release.

To execute all five headlessly and check the committed shape:

```bash
python tools/check_notebooks.py
```

`--write` refreshes the committed outputs. The tool is refused under CI; the CI job executes
the notebooks without rewriting them.

## The one wheel-specific trap

`run_flowsheet` and `validate` default to an on-disk `specs/unit_ops` palette that the wheel does
not ship. From an installed package, use `azoth.process.Session`, which opens the palette compiled
into the extension, or pass a real `palette_dir`. Notebook 04 demonstrates `Session`.
