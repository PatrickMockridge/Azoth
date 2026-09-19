//! The data tables, exposed to Python so the two sides can be compared.
//!
//! It exposes three things:
//!
//! * the parsed rows, so every field can be compared row by row;
//! * the **raw embedded text**, so the comparison can be made on bytes rather than on
//!   values. Parsed values agreeing is a weaker claim: two different files can parse
//!   to the same values, and a stale copy bundled into a wheel would look precisely
//!   like that - an agreement that proves nothing about which file was read;
//! * the repo-relative path of each file, so Python knows which file on disk is
//!   supposed to be the same one, rather than hardcoding the mapping in a test.
//!
//! Not a data API. Nothing a calculation needs comes from here - the calcs read their
//! own tables inside the crate, and the Python reference reads the files. This is
//! introspection in the same sense `warning_codes` and `result_fields` are: it exists
//! so a test can assert a cross-language claim instead of asserting it in prose.

use azoth_eos::association::SiteScheme;
use azoth_eos::databank;
use azoth_hydraulics::{fittings, fluids};
use pyo3::prelude::*;

/// One file this build embeds, with its bytes.
#[pyclass(frozen, skip_from_py_object, module = "azoth._core", name = "DataFile")]
#[derive(Debug, Clone, PartialEq)]
pub struct PyDataFile {
    /// The name the file is addressed by, e.g. `water` or `fittings`.
    #[pyo3(get)]
    pub name: String,
    /// Repo-relative path, e.g. `data/fluids/water.csv`.
    #[pyo3(get)]
    pub path: String,
    /// The exact text this build embedded. UTF-8, because `include_str!` requires it.
    #[pyo3(get)]
    pub text: String,
}

#[pymethods]
impl PyDataFile {
    fn __repr__(&self) -> String {
        format!("DataFile({}, {} bytes)", self.path, self.text.len())
    }
}

/// One row of the fittings registry, transported.
#[pyclass(
    frozen,
    skip_from_py_object,
    module = "azoth._core",
    name = "FittingRow"
)]
#[derive(Debug, Clone, PartialEq)]
pub struct PyFittingRow {
    /// Stable identifier.
    #[pyo3(get)]
    pub id: String,
    /// Coarse category.
    #[pyo3(get)]
    pub family: String,
    /// Human-readable name.
    #[pyo3(get)]
    pub name: String,
    /// Equivalent length ratio, `L_eq / D`.
    #[pyo3(get)]
    pub n_ld: f64,
    /// What the coefficient is multiplied by, conventionally `f_t`.
    #[pyo3(get)]
    pub f_t_basis: String,
    /// Where the value came from, or a statement that it came from nowhere.
    #[pyo3(get)]
    pub citation: String,
    /// `verified`, `unverified` or `estimated_dummy`.
    #[pyo3(get)]
    pub verify_status: String,
}

/// One row of a fluid property table, transported.
///
/// Carries the fluid's name on every row rather than once per table, so a comparison
/// in Python can flatten every row of every fluid into one list and still know which
/// fluid each came from. A table-per-call shape would need the test to iterate fluids
/// and tables separately, which is one more place to get the loops wrong.
#[pyclass(frozen, skip_from_py_object, module = "azoth._core", name = "FluidRow")]
#[derive(Debug, Clone, PartialEq)]
pub struct PyFluidRow {
    /// The fluid this row belongs to.
    #[pyo3(get)]
    pub fluid: String,
    /// Temperature in degrees Celsius.
    #[pyo3(get)]
    pub temperature_c: f64,
    /// Density in kg/m^3.
    #[pyo3(get)]
    pub density_kg_m3: f64,
    /// Dynamic viscosity in Pa*s.
    #[pyo3(get)]
    pub dynamic_viscosity_pa_s: f64,
    /// Where the value came from.
    #[pyo3(get)]
    pub citation: String,
    /// `verified`, `unverified` or `estimated_dummy`.
    #[pyo3(get)]
    pub verify_status: String,
}

/// One row of the component databank, transported.
///
/// The columns are the ones a calculation reads. The file carries more - the CAS
/// number, the formula, the liquid density - and those are deliberately not here: a
/// cross-language test that compared them would be testing a field neither
/// implementation acts on, and the fields that matter would still need their own case.
#[pyclass(
    frozen,
    skip_from_py_object,
    module = "azoth._core",
    name = "ComponentRow"
)]
#[derive(Debug, Clone, PartialEq)]
pub struct PyComponentRow {
    /// Lower-case name, which is the key the table is addressed by.
    #[pyo3(get)]
    pub name: String,
    /// Critical temperature, in K.
    #[pyo3(get)]
    pub tc_k: f64,
    /// Critical pressure, in Pa.
    #[pyo3(get)]
    pub pc_pa: f64,
    /// Acentric factor, dimensionless.
    #[pyo3(get)]
    pub acentric_factor: f64,
    /// The five `Cp` coefficients, in J/(mol*K**n).
    ///
    /// Optional because the row shape is also what an *overlay* reports, and a
    /// substance a keycard adds has no polynomial: a keycard supplies the parameters a
    /// cubic reads, and this is not one of them. Every row of the shipped table has
    /// all five, so a `None` here can only come from an overlay.
    #[pyo3(get)]
    pub cp_a: Option<f64>,
    #[pyo3(get)]
    pub cp_b: Option<f64>,
    #[pyo3(get)]
    pub cp_c: Option<f64>,
    #[pyo3(get)]
    pub cp_d: Option<f64>,
    #[pyo3(get)]
    pub cp_e: Option<f64>,
    /// NeqSim's `COMPTYPE` for the row: the substance class.
    ///
    /// A string rather than a flag because the table carries ten of them, and `ion` is
    /// the one this crate acts on - it is what `mixture_of` refuses a cubic over. The
    /// others are reported so a comparison sees the whole field rather than the one value
    /// currently read.
    #[pyo3(get)]
    pub component_type: String,
    /// The ionic charge as a **charge number**, in units of the elementary charge.
    #[pyo3(get)]
    pub ionic_charge: f64,
    /// The Deshmukh-Mather ion diameter, in ångström, as the table stores it.
    #[pyo3(get)]
    pub deshmukh_mather_diameter: f64,
    /// The five dielectric-constant coefficients `d0`..`d4`.
    #[pyo3(get)]
    pub dielectric: [f64; 5],
    /// The association site scheme's name - `"1A"`, `"2A"`, `"2B"` or `"4C"` - or the
    /// empty string for a component the table gives no scheme.
    ///
    /// A string rather than a number because the scheme *name* is what NeqSim switches
    /// on: the site count cannot reconstruct it, since the table carries `1A` at zero
    /// sites and both `2A` and `2B` at two.
    #[pyo3(get)]
    pub association_scheme: String,
    /// The site count the table states, zero for a component with no scheme.
    #[pyo3(get)]
    pub association_sites: u32,
    /// The association energy `eps`, in J/mol.
    #[pyo3(get)]
    pub association_energy: f64,
    /// `kappa_AB` for the SRK family.
    #[pyo3(get)]
    pub association_volume_srk: f64,
    /// The fitted SRK-CPA attraction, in NeqSim's internal scale.
    #[pyo3(get)]
    pub association_a_srk: f64,
    /// The fitted SRK-CPA covolume, in NeqSim's internal scale.
    #[pyo3(get)]
    pub association_b_srk: f64,
    /// The SRK alpha correlation's `m`.
    #[pyo3(get)]
    pub association_m_srk: f64,
    /// `kappa_AB` for the PR family.
    #[pyo3(get)]
    pub association_volume_pr: f64,
    /// The fitted PR-CPA attraction, in NeqSim's internal scale.
    #[pyo3(get)]
    pub association_a_pr: f64,
    /// The fitted PR-CPA covolume, in NeqSim's internal scale.
    #[pyo3(get)]
    pub association_b_pr: f64,
    /// The PR alpha correlation's `m`.
    #[pyo3(get)]
    pub association_m_pr: f64,
    /// `racketZCPA`, the Rackett compressibility NeqSim's CPA volume correction reads.
    #[pyo3(get)]
    pub association_racket_z: f64,
    /// `volcorrCPA_T`, the CPA volume-translation coefficient.
    #[pyo3(get)]
    pub association_volume_correction: f64,
}

/// One row of the interaction table, transported.
#[pyclass(frozen, skip_from_py_object, module = "azoth._core", name = "KijRow")]
#[derive(Debug, Clone, PartialEq)]
pub struct PyKijRow {
    /// The pair, with the lower-sorting name first.
    #[pyo3(get)]
    pub component_a: String,
    #[pyo3(get)]
    pub component_b: String,
    /// The Peng-Robinson binary interaction parameter, `KIJPR`.
    #[pyo3(get)]
    pub kij_pr: f64,
    /// The Soave-Redlich-Kwong binary interaction parameter, `KIJSRK`, which every cubic
    /// but Peng-Robinson reads.
    #[pyo3(get)]
    pub kij_srk: f64,
}

/// Every data file this build embeds, with its bytes.
///
/// # Errors
/// Returns an error if the embedded text is malformed, which is a build-time
/// invariant rather than a user condition.
#[pyfunction]
#[must_use]
pub fn data_files() -> Vec<PyDataFile> {
    let mut out = vec![
        PyDataFile {
            name: "fittings".to_string(),
            path: fittings::embedded_path().to_string(),
            text: fittings::embedded_csv().to_string(),
        },
        PyDataFile {
            name: "components".to_string(),
            path: databank::COMPONENTS_PATH.to_string(),
            text: databank::embedded_components().to_string(),
        },
        PyDataFile {
            name: "kij".to_string(),
            path: databank::KIJ_PATH.to_string(),
            text: databank::embedded_kij().to_string(),
        },
        PyDataFile {
            name: "unifaccomp".to_string(),
            path: databank::UNIFAC_COMP_PATH.to_string(),
            text: databank::embedded_unifac_comp().to_string(),
        },
        PyDataFile {
            name: "unifacgroupparam".to_string(),
            path: databank::UNIFAC_GROUP_PATH.to_string(),
            text: databank::embedded_unifac_group().to_string(),
        },
        PyDataFile {
            name: "unifacinterparam".to_string(),
            path: databank::UNIFAC_INTER_PATH.to_string(),
            text: databank::embedded_unifac_inter().to_string(),
        },
        PyDataFile {
            name: "mbwr32".to_string(),
            path: databank::MBWR32_PATH.to_string(),
            text: databank::embedded_mbwr32().to_string(),
        },
    ];
    for fluid in fluids::available_fluids() {
        if let (Some(path), Some(text)) =
            (fluids::embedded_path(fluid), fluids::embedded_csv(fluid))
        {
            out.push(PyDataFile {
                name: fluid.to_string(),
                path: path.to_string(),
                text: text.to_string(),
            });
        }
    }
    out
}

/// Every row of the fittings registry, as this crate parsed it.
///
/// # Errors
/// Returns an error if the embedded registry is malformed.
#[pyfunction]
pub fn fittings_rows(py: Python<'_>) -> PyResult<Vec<PyFittingRow>> {
    let rows = fittings::registry().map_err(|e| crate::errors::to_pyerr(py, e))?;
    Ok(rows
        .iter()
        .map(|row| PyFittingRow {
            id: row.id.clone(),
            family: row.family.clone(),
            name: row.name.clone(),
            n_ld: row.n_ld,
            f_t_basis: row.f_t_basis.clone(),
            citation: row.citation.clone(),
            verify_status: row.status.as_str().to_string(),
        })
        .collect())
}

/// Every row of the component databank, as the `azoth-eos` crate parsed it.
#[pyfunction]
#[must_use]
pub fn component_rows() -> Vec<PyComponentRow> {
    databank::all_entries().into_iter().map(row_of).collect()
}

/// One entry in the table's own row shape.
///
/// One definition, used by the table and by an overlay's rows as well: the two are
/// compared field by field in `python/tests/test_data_agreement.py`, and a second
/// builder would be a second place the shape can differ.
pub(crate) fn row_of(entry: &databank::Entry) -> PyComponentRow {
    PyComponentRow {
        name: entry.name.clone(),
        tc_k: entry.tc,
        pc_pa: entry.pc,
        acentric_factor: entry.omega,
        cp_a: entry.cp.map(|cp| cp[0]),
        cp_b: entry.cp.map(|cp| cp[1]),
        cp_c: entry.cp.map(|cp| cp[2]),
        cp_d: entry.cp.map(|cp| cp[3]),
        cp_e: entry.cp.map(|cp| cp[4]),
        component_type: entry.class.clone(),
        ionic_charge: entry.ionic_charge,
        deshmukh_mather_diameter: entry.deshmukh_mather_diameter,
        dielectric: entry.dielectric,
        association_scheme: entry
            .association
            .as_ref()
            .map(|a| scheme_name(a.scheme).to_string())
            .unwrap_or_default(),
        association_sites: entry.association.as_ref().map_or(0, |a| a.sites),
        association_energy: entry.association.as_ref().map_or(0.0, |a| a.energy),
        association_volume_srk: entry.association.as_ref().map_or(0.0, |a| a.volume_srk),
        association_a_srk: entry.association.as_ref().map_or(0.0, |a| a.a_srk),
        association_b_srk: entry.association.as_ref().map_or(0.0, |a| a.b_srk),
        association_m_srk: entry.association.as_ref().map_or(0.0, |a| a.m_srk),
        association_volume_pr: entry.association.as_ref().map_or(0.0, |a| a.volume_pr),
        association_a_pr: entry.association.as_ref().map_or(0.0, |a| a.a_pr),
        association_b_pr: entry.association.as_ref().map_or(0.0, |a| a.b_pr),
        association_m_pr: entry.association.as_ref().map_or(0.0, |a| a.m_pr),
        association_racket_z: entry.association.as_ref().map_or(0.0, |a| a.racket_z),
        association_volume_correction: entry
            .association
            .as_ref()
            .map_or(0.0, |a| a.volume_correction),
    }
}

/// The databank name of a site scheme, the inverse of `SiteScheme::from_databank_name`.
///
/// The two are written out in both directions rather than one deriving from the other,
/// because a round trip through a wrong name would agree with itself.
///
/// [`SiteScheme::NonAssociating`] maps to the empty string, which is also what a row with
/// no scheme reports - and it is unreachable here besides: the loader answers `None` to
/// the table's `0`, so a row never carries it.
fn scheme_name(scheme: SiteScheme) -> &'static str {
    match scheme {
        SiteScheme::NonAssociating => "",
        SiteScheme::OneA => "1A",
        SiteScheme::TwoA => "2A",
        SiteScheme::TwoB => "2B",
        SiteScheme::FourC => "4C",
    }
}

/// Every row of the interaction table, as the `azoth-eos` crate parsed it.
#[pyfunction]
#[must_use]
pub fn kij_rows() -> Vec<PyKijRow> {
    databank::all_kij()
        .into_iter()
        .map(|(component_a, component_b, kij_pr, kij_srk)| PyKijRow {
            component_a,
            component_b,
            kij_pr,
            kij_srk,
        })
        .collect()
}

/// Every row of a built-in fluid's table, as this crate parsed it.
///
/// # Errors
/// Returns an error for a fluid this build does not carry, rather than an empty list:
/// an unknown fluid and a fluid with no rows are different problems, and returning
/// `[]` for both would make a typo look like an empty table.
#[pyfunction]
pub fn fluid_rows(py: Python<'_>, name: &str) -> PyResult<Vec<PyFluidRow>> {
    let table = fluids::provider_for(name).map_err(|e| crate::errors::to_pyerr(py, e))?;
    Ok(table
        .points
        .iter()
        .map(|point| PyFluidRow {
            fluid: table.name.clone(),
            temperature_c: point.temperature_c,
            density_kg_m3: point.density_kg_m3,
            dynamic_viscosity_pa_s: point.dynamic_viscosity_pa_s,
            citation: point.citation.clone(),
            verify_status: point.status.as_str().to_string(),
        })
        .collect())
}
