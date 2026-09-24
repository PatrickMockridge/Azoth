//! ISO 6976's per-component constants, read from `data/standards/iso6976.csv`.
//!
//! The CSV is embedded with `include_str!` rather than read at runtime, for the reason
//! `azoth-hydraulics`' fittings registry gives: it removes any question of which file was
//! loaded, and the Python side opens the same path from the repository, so the two
//! implementations read one file rather than two copies that are supposed to match. The two
//! *are* checked against each other - `python/tests/test_data_agreement.py` compares the
//! parsed rows and the embedded bytes.
//!
//! **The table is the standard's, and it is not the component databank.** It carries
//! exactly the components ISO 6976 defines, which is narrower: 40 of the standard's 56 rows
//! name a substance the databank has, and a gas naming one of the other sixteen has no
//! calorific value here. The manifest dispositions all 56 with the reason.

use azoth_core::{AzothError, Result};

/// How close a stated reference temperature has to be to one the table carries.
///
/// **A tolerance rather than equality, because the caller states kelvin.** The standard's
/// temperatures are stated in degrees Celsius and the library's convention is kelvin, so a
/// caller writes `288.7` for the 60 °F reference - and `288.7 - 273.15` is
/// `15.550000000000011`, which `==` would refuse. NeqSim compares exactly, but its values
/// are its own literals; a caller's are not, and a refusal over the last bit of a
/// subtraction would be a port that cannot be used.
const TOLERANCE: f64 = 1.0e-9;

/// The embedded table. Path is relative to this source file.
const ISO6976_CSV: &str = include_str!("../../../data/standards/iso6976.csv");

/// One component's row: the standard's tabulated properties, in SI.
///
/// The reference temperatures are the standard's own, and each field's name carries the one
/// it belongs to: `compression_factor_0c` is `Z` at 0 °C, `superior_calorific_value_j_per_mol_60f`
/// is `Hsup` at 60 °F. A caller selects the pair its reference temperatures name - the table
/// has no column for an arbitrary temperature, and the standard defines none.
#[derive(Debug, Clone, PartialEq)]
pub struct Iso6976Row {
    /// The component's name, as the component databank spells it.
    pub name: String,
    /// Molar mass, kg/mol.
    pub molar_mass: f64,
    /// The compression factor at 0 °C.
    pub compression_factor_0c: f64,
    /// The compression factor at 15 °C - and at 15.55 °C, whose column the table does not carry.
    pub compression_factor_15c: f64,
    /// The compression factor at 20 °C.
    pub compression_factor_20c: f64,
    /// The summation factor at 0 °C, whose z-weighted sum is squared and subtracted from one.
    pub summation_factor_0c: f64,
    /// The summation factor at 15 °C.
    pub summation_factor_15c: f64,
    /// The summation factor at 20 °C.
    pub summation_factor_20c: f64,
    /// Superior (gross) molar calorific value at 0 °C, J/mol.
    pub superior_calorific_value_0c: f64,
    /// Superior calorific value at 15 °C, J/mol.
    pub superior_calorific_value_15c: f64,
    /// Superior calorific value at 20 °C, J/mol.
    pub superior_calorific_value_20c: f64,
    /// Superior calorific value at 25 °C, J/mol.
    pub superior_calorific_value_25c: f64,
    /// Superior calorific value at 60 °F, J/mol.
    pub superior_calorific_value_60f: f64,
    /// Inferior (net) molar calorific value at 0 °C, J/mol.
    pub inferior_calorific_value_0c: f64,
    /// Inferior calorific value at 15 °C, J/mol.
    pub inferior_calorific_value_15c: f64,
    /// Inferior calorific value at 20 °C, J/mol.
    pub inferior_calorific_value_20c: f64,
    /// Inferior calorific value at 25 °C, J/mol.
    pub inferior_calorific_value_25c: f64,
    /// Inferior calorific value at 60 °F, J/mol.
    pub inferior_calorific_value_60f: f64,
    /// How many carbon atoms the molecule carries, which the standard's tables index by.
    pub carbon_count: u32,
    /// Where the row came from.
    pub citation: String,
}

impl Iso6976Row {
    /// The compression factor at a reference temperature in degrees Celsius.
    ///
    /// # Errors
    /// [`AzothError::InvalidInput`] for a temperature the standard's table has no column
    /// for, which is every one outside `{0, 15, 15.55, 20}`.
    pub fn compression_factor_at(&self, celsius: f64) -> Result<f64> {
        if (celsius - 0.0).abs() < TOLERANCE {
            Ok(self.compression_factor_0c)
        } else if (celsius - 15.0).abs() < TOLERANCE || (celsius - 15.55).abs() < TOLERANCE {
            Ok(self.compression_factor_15c)
        } else if (celsius - 20.0).abs() < TOLERANCE {
            Ok(self.compression_factor_20c)
        } else {
            Err(AzothError::invalid_input(
                "volumetric_reference_temperature",
                format!(
                    "the standard's table carries a compression factor at 0, 15 and 20 C, and \
                     {celsius} C is not one of them (15.55 C is the 60 F reference and reads \
                     the 15 C column)"
                ),
            ))
        }
    }

    /// The summation factor at a reference temperature in degrees Celsius.
    ///
    /// # Errors
    /// [`AzothError::InvalidInput`] for a temperature the standard's table has no column
    /// for, which is every one outside `{0, 15, 15.55, 20}`.
    pub fn summation_factor_at(&self, celsius: f64) -> Result<f64> {
        if (celsius - 0.0).abs() < TOLERANCE {
            Ok(self.summation_factor_0c)
        } else if (celsius - 15.0).abs() < TOLERANCE || (celsius - 15.55).abs() < TOLERANCE {
            Ok(self.summation_factor_15c)
        } else if (celsius - 20.0).abs() < TOLERANCE {
            Ok(self.summation_factor_20c)
        } else {
            Err(AzothError::invalid_input(
                "volumetric_reference_temperature",
                format!(
                    "the standard's table carries a summation factor at 0, 15 and 20 C, and \
                     {celsius} C is not one of them (15.55 C is the 60 F reference and reads \
                     the 15 C column)"
                ),
            ))
        }
    }

    /// A calorific value at a reference temperature in degrees Celsius.
    ///
    /// # Errors
    /// [`AzothError::InvalidInput`] for a temperature the standard's table has no column
    /// for, which is every one outside `{0, 15, 15.55, 20, 25}`.
    pub fn calorific_value_at(&self, celsius: f64, superior: bool) -> Result<f64> {
        let column = if (celsius - 0.0).abs() < TOLERANCE {
            (
                self.superior_calorific_value_0c,
                self.inferior_calorific_value_0c,
            )
        } else if (celsius - 15.0).abs() < TOLERANCE {
            (
                self.superior_calorific_value_15c,
                self.inferior_calorific_value_15c,
            )
        } else if (celsius - 15.55).abs() < TOLERANCE {
            (
                self.superior_calorific_value_60f,
                self.inferior_calorific_value_60f,
            )
        } else if (celsius - 20.0).abs() < TOLERANCE {
            (
                self.superior_calorific_value_20c,
                self.inferior_calorific_value_20c,
            )
        } else if (celsius - 25.0).abs() < TOLERANCE {
            (
                self.superior_calorific_value_25c,
                self.inferior_calorific_value_25c,
            )
        } else {
            return Err(AzothError::invalid_input(
                "energy_reference_temperature",
                format!(
                    "the standard's table carries a calorific value at 0, 15, 15.55, 20 and \
                     25 C, and {celsius} C is not one of them"
                ),
            ));
        };
        Ok(if superior { column.0 } else { column.1 })
    }
}

/// The whole table, parsed.
///
/// # Errors
/// [`AzothError::InvalidInput`] if the embedded table is malformed, which is a build
/// problem rather than a caller's: the file is generated by `tools/gen_databank.py` and
/// checked by `--check`.
pub fn table() -> Result<Vec<Iso6976Row>> {
    let mut reader = csv::Reader::from_reader(ISO6976_CSV.as_bytes());
    let headers = reader
        .headers()
        .map_err(|e| AzothError::invalid_input("iso6976", format!("no header row: {e}")))?
        .clone();
    let index = |name: &str| -> Result<usize> {
        headers.iter().position(|h| h == name).ok_or_else(|| {
            AzothError::invalid_input(
                "iso6976",
                format!("the embedded table has no `{name}` column"),
            )
        })
    };
    let mut rows = Vec::new();
    for record in reader.records() {
        let record = record
            .map_err(|e| AzothError::invalid_input("iso6976", format!("malformed CSV row: {e}")))?;
        let number = |name: &str| -> Result<f64> {
            let raw = record.get(index(name)?).unwrap_or("").trim();
            raw.parse().map_err(|e| {
                AzothError::invalid_input("iso6976", format!("`{name}` is {raw:?}: {e}"))
            })
        };
        let text = |name: &str| -> Result<String> {
            Ok(record
                .get(index(name)?)
                .unwrap_or("")
                .trim()
                .trim_matches('"')
                .to_string())
        };
        let carbon = record.get(index("carbon_count")?).unwrap_or("0").trim();
        rows.push(Iso6976Row {
            name: text("name")?,
            molar_mass: number("molar_mass_kg_per_mol")?,
            compression_factor_0c: number("compression_factor_0c")?,
            compression_factor_15c: number("compression_factor_15c")?,
            compression_factor_20c: number("compression_factor_20c")?,
            summation_factor_0c: number("summation_factor_0c")?,
            summation_factor_15c: number("summation_factor_15c")?,
            summation_factor_20c: number("summation_factor_20c")?,
            superior_calorific_value_0c: number("superior_calorific_value_j_per_mol_0c")?,
            superior_calorific_value_15c: number("superior_calorific_value_j_per_mol_15c")?,
            superior_calorific_value_20c: number("superior_calorific_value_j_per_mol_20c")?,
            superior_calorific_value_25c: number("superior_calorific_value_j_per_mol_25c")?,
            superior_calorific_value_60f: number("superior_calorific_value_j_per_mol_60f")?,
            inferior_calorific_value_0c: number("inferior_calorific_value_j_per_mol_0c")?,
            inferior_calorific_value_15c: number("inferior_calorific_value_j_per_mol_15c")?,
            inferior_calorific_value_20c: number("inferior_calorific_value_j_per_mol_20c")?,
            inferior_calorific_value_25c: number("inferior_calorific_value_j_per_mol_25c")?,
            inferior_calorific_value_60f: number("inferior_calorific_value_j_per_mol_60f")?,
            carbon_count: carbon.parse().map_err(|e| {
                AzothError::invalid_input("iso6976", format!("`carbon_count` is {carbon:?}: {e}"))
            })?,
            citation: text("citation")?,
        });
    }
    Ok(rows)
}

/// The row for one component, matched by name.
///
/// # Errors
/// [`AzothError::InvalidInput`] for a name the standard's table does not carry, **including
/// the sixteen rows the compiled table drops**: a gas containing propylene or carbon
/// monoxide has no ISO 6976 calorific value in this library rather than a partial one, and
/// the message says which side of that the name is on.
pub fn row(name: &str) -> Result<Iso6976Row> {
    let wanted = name.trim().to_lowercase();
    table()?
        .into_iter()
        .find(|row| row.name == wanted)
        .ok_or_else(|| {
            AzothError::invalid_input(
                "components",
                format!(
                    "{name} has no row in ISO 6976's table. The standard defines 56 \
                     substances and this library compiles the 40 the component databank also \
                     carries; a gas naming one of the others has no calorific value here."
                ),
            )
        })
}

/// The embedded bytes, so a test can compare them against the file Python reads.
#[must_use]
pub fn embedded_csv() -> &'static str {
    ISO6976_CSV
}
