//! The shell-and-tube material table, vendored from NeqSim's design database.
//!
//! `ShellAndTubeDesignCalculator.loadMaterialFromDB` reads the `HeatExchangerTubeMaterials`
//! table by exact grade; this is that table, compiled to `data/heat_exchanger/tube_materials.csv`
//! and embedded. **`nace_compliant` stays the file's text**: the class tests it against `Yes`,
//! which no row holds, and the port reproduces that defect rather than repairing it.

use std::sync::LazyLock;

/// Where the compiled table lives, as the extension reports it.
pub const TABLE_PATH: &str = "data/heat_exchanger/tube_materials.csv";

/// The compiled table.
const TABLE_CSV: &str = include_str!("../../../data/heat_exchanger/tube_materials.csv");

/// One material grade, as the class reads it.
#[derive(Debug, Clone, PartialEq)]
pub struct TubeMaterial {
    /// The grade key, matched by exact equality.
    pub material_grade: String,
    /// The specification string the class prefixes the grade with.
    pub specification: String,
    /// The material family.
    pub material_type: String,
    /// The minimum yield strength, MPa.
    pub smys_mpa: f64,
    /// The allowable stress, MPa.
    pub allowable_stress_mpa: f64,
    /// The thermal conductivity, W/(m K).
    pub thermal_conductivity_wpm_k: f64,
    /// The NACE text, `true` or `false` in the file.
    pub nace_compliant: String,
}

impl TubeMaterial {
    /// The class's own test: `"Yes".equalsIgnoreCase(text)`, which no row satisfies.
    #[must_use]
    pub fn is_nace_compliant(&self) -> bool {
        self.nace_compliant.eq_ignore_ascii_case("Yes")
    }
}

/// The compiled table, parsed once.
static TABLE: LazyLock<Vec<TubeMaterial>> = LazyLock::new(|| {
    TABLE_CSV
        .lines()
        .skip(1)
        .filter(|line| !line.trim().is_empty())
        .map(|line| {
            let fields: Vec<&str> = line.split(',').collect();
            TubeMaterial {
                material_grade: fields[0].to_string(),
                specification: fields[1].to_string(),
                material_type: fields[2].to_string(),
                smys_mpa: fields[3].parse().unwrap_or(0.0),
                allowable_stress_mpa: fields[4].parse().unwrap_or(0.0),
                thermal_conductivity_wpm_k: fields[5].parse().unwrap_or(0.0),
                nace_compliant: fields[6].to_string(),
            }
        })
        .collect()
});

/// The embedded table, for the extension's `data_files()`.
#[must_use]
pub fn embedded_csv() -> &'static str {
    TABLE_CSV
}

/// The embedded table's path.
#[must_use]
pub fn embedded_path() -> &'static str {
    TABLE_PATH
}

/// Looks a grade up by exact equality; the first row wins, as a SQL lookup would.
#[must_use]
pub fn lookup(grade: &str) -> Option<&'static TubeMaterial> {
    TABLE
        .iter()
        .find(|material| material.material_grade == grade)
}

#[cfg(test)]
mod tests {
    use super::lookup;

    #[test]
    fn every_row_reads_and_the_nace_text_is_never_yes() {
        let carbon = lookup("SA-179").expect("the default tube grade");
        assert_eq!(carbon.allowable_stress_mpa, 92.0);
        assert!(!carbon.is_nace_compliant());
        let stainless = lookup("SA-213-TP316").expect("a stainless grade");
        assert_eq!(stainless.thermal_conductivity_wpm_k, 16.2);
        // The file's `NACECompliant` is `true`, and the class's test is `Yes`.
        assert_eq!(stainless.nace_compliant, "true");
        assert!(!stainless.is_nace_compliant());
        assert!(lookup("SA-999").is_none());
    }
}
