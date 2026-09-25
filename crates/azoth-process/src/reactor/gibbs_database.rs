//! `GibbsReactor`'s own species database, and the per-species arithmetic that reads it.
//!
//! **This is not the databank's formation properties.** `GIBBSENERGYOFFORMATION` is `used`
//! from P10 and its values are not these; `GibbsReactor` carries its own two files and its
//! own correlations, and the whole point of this tier is that a port composed from the
//! databank's properties would answer with different numbers than the class it ports.
//!
//! The two files are embedded rather than read at runtime, for the reason
//! `azoth-standards`' ISO 6976 constants give: it removes any question of which file was
//! loaded. `tools/gen_databank.py` writes them and `databank/manifest.toml` dispositions
//! every column of both.
//!
//! **Three facts about the source table are load-bearing and are reproduced rather than
//! tidied**, each measured against the file:
//!
//! * The eight counts are `O, N, C, H, S, Na, Ar, Z` and **`Z` is a charge**, not a count:
//!   `OH-` is `-1` and `SO4--` is `-2`.
//! * The class names only seven of them. `elementNames` is
//!   `{"O","N","C","H","S","Ar","Z"}`, so its `Ar` balance reads the file's **`Na`**
//!   column and its `Z` balance reads the file's **`Ar`** column - which is zero on all 79
//!   rows, so that balance row is never active and its Lagrange multiplier is never solved
//!   for. The charge is read into `elements[7]` and no loop reaches it: **the solve does
//!   not balance charge**.
//! * `calculateJ`'s javadoc and its code disagree on the sign of three terms. The code is
//!   what is ported; see [`GibbsSpecies::j_term`].

use std::collections::BTreeMap;

use azoth_core::{AzothError, Result};

/// `GibbsReactor.java:427`'s `elementNames`, which has **seven** entries for an
/// eight-column table.
///
/// Carried as the class writes it, including the shift: index 5 is named `Ar` and reads
/// the file's `Na` column, index 6 is named `Z` and reads the file's `Ar` column.
pub const ELEMENT_NAMES: [&str; 7] = ["O", "N", "C", "H", "S", "Ar", "Z"];

/// `REFERENCE_TEMPERATURE`, the 298.15 K the formation properties are stated at.
pub const REFERENCE_TEMPERATURE: f64 = 298.15;

/// `R_KJ`, the gas constant the class's Gibbs arithmetic uses, kJ/(mol·K).
///
/// **A fourth gas constant in this workspace**, and not interchangeable with the others:
/// ISO 6976's table uses `8.314510`, the Sm³ conversions `8.3144621`, and
/// `KineticReaction` `8.31446`. This one is CODATA's 2018 value and it is the class's.
pub const R_KJ: f64 = 8.314_462_618e-3;

/// The four heat-capacity coefficients the element subtraction reads, per element, in the
/// order `O, N, C, H, S` - `calculateCorrectedHeatCapacityCoeffs`' own five arrays.
const ELEMENT_CP: [[f64; 4]; 5] = [
    [12.73, 7.60e-3, -3.58e-6, 6.56e-10],
    [14.4415, -7.85e-4, 4.04e-6, -1.44e-9],
    [8.43, 0.0, 0.0, 0.0],
    [14.544, -9.60e-4, 2.00e-6, -4.35e-10],
    [17.815, 0.001, 0.000, 0.000],
];

/// The compiled species table. Path is relative to this source file.
const SPECIES_CSV: &str = include_str!("../../../../data/reactors/gibbs_reactor.csv");

/// The compiled coefficient table, joined to `SPECIES_CSV` by lowercased molecule name.
const COEFFICIENTS_CSV: &str = include_str!("../../../../data/reactors/gibbs_reactor_coeffs.csv");

/// One species' order-5 Gibbs and enthalpy polynomials.
///
/// Both are evaluated as `c0*T^5 + c1*T^4 + c2*T^3 + c3*T^2 + c4*T + c5` and returned
/// with no unit conversion, so both are kJ/mol. Measured on CO2 at 298.15 K the Gibbs
/// polynomial gives `-394.2203` against the same row's tabulated `-394.4` kJ/mol, so the
/// two are independent fits of one quantity and not a redundant copy of it.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Polynomials {
    /// The Gibbs energy polynomial's six coefficients, highest power first.
    pub gibbs: [f64; 6],
    /// The enthalpy polynomial's six coefficients, highest power first.
    pub enthalpy: [f64; 6],
}

/// One row of the class's database: an element vector, three formation properties, and
/// optionally the two polynomials.
#[derive(Debug, Clone, PartialEq)]
pub struct GibbsSpecies {
    /// The file's own spelling, which is the join key lowercased.
    pub name: String,
    /// The file's eight counts, in its own column order: O, N, C, H, S, Na, Ar, charge.
    ///
    /// **Eight, and the class reads seven of them under seven names.** See
    /// [`ELEMENT_NAMES`] and [`GibbsSpecies::balance_elements`].
    pub elements: [f64; 8],
    /// The standard enthalpy of formation at 298.15 K, kJ/mol.
    pub hf298: f64,
    /// The standard Gibbs energy of formation at 298.15 K, kJ/mol.
    pub gf298: f64,
    /// The standard molar entropy at 298.15 K, J/(mol·K).
    pub sf298: f64,
    /// The two polynomials, or `None` for the 34 of the 49 species the coefficient file
    /// does not carry - which is the class's `Double.isNaN` test for its fallback branch.
    pub polynomials: Option<Polynomials>,
}

impl GibbsSpecies {
    /// The eight counts as the class's *element balance* reads them: seven entries, under
    /// `O, N, C, H, S, Ar, Z`.
    ///
    /// **This drops the charge and shifts the last two names**, which is the class's own
    /// arithmetic rather than a tidying of it - see the module doc. `balance[5]` is the
    /// file's `Na` count and `balance[6]` is the file's `Ar` count.
    #[must_use]
    pub fn balance_elements(&self) -> [f64; 7] {
        let mut out = [0.0; 7];
        out.copy_from_slice(&self.elements[..7]);
        out
    }

    /// `calculateCorrectedHeatCapacityCoeffs`: the component's Cp polynomial less the
    /// contribution of its constituent elements, `[dA, dB, dC, dD]`.
    ///
    /// `cp` is the component's **own** `A, B, C, D` - NeqSim's `getCpA()` through
    /// `getCpD()`, which for azoth is the component databank's `cp_a` through `cp_d`. It
    /// is an argument rather than a field because it is a property of the fluid's
    /// component and not of this database's row.
    ///
    /// Only the first five counts are read, so the `Na`/`Ar`/charge columns do not enter.
    #[must_use]
    pub fn corrected_heat_capacity_coefficients(&self, cp: [f64; 4]) -> [f64; 4] {
        let mut corrected = cp;
        for (element, coefficients) in ELEMENT_CP.iter().enumerate() {
            for (power, coefficient) in coefficients.iter().enumerate() {
                corrected[power] -= self.elements[element] * coefficient;
            }
        }
        corrected
    }

    /// `calculateJ`: the corrected formation enthalpy, kJ/mol.
    ///
    /// **The code, not the javadoc.** The comment says
    /// `J = ΔH°f - ΔA·TR - ΔB/2·TR² - ΔC/3·TR³ - ΔD/4·TR⁴`, which would subtract every
    /// term; the code is `deltaHf298 - (dA*TR - dB/2*TR² - dC/3*TR³ - dD/4*TR⁴)/1000`,
    /// which subtracts the first and adds the other three. The two differ for any species
    /// whose corrected `B`, `C` or `D` is non-zero, which is most of them, so this is not
    /// a comment that can be followed. The arithmetic is ported as written and the
    /// divergence is recorded here and in the spec.
    #[must_use]
    pub fn j_term(&self, cp: [f64; 4]) -> f64 {
        let [da, db, dc, dd] = self.corrected_heat_capacity_coefficients(cp);
        let tr = REFERENCE_TEMPERATURE;
        self.hf298
            - (da * tr - db / 2.0 * tr.powi(2) - dc / 3.0 * tr.powi(3) - dd / 4.0 * tr.powi(4))
                / 1000.0
    }

    /// `calculateI`, the reference-temperature term of the fallback Gibbs branch.
    #[must_use]
    pub fn i_term(&self, cp: [f64; 4]) -> f64 {
        let [da, db, dc, dd] = self.corrected_heat_capacity_coefficients(cp);
        let tr = REFERENCE_TEMPERATURE;
        (1.0 / R_KJ)
            * (-(self.j_term(cp) / tr)
                + (da * tr.ln() + db / 2.0 * tr + dc / 6.0 * tr.powi(2) + dd / 12.0 * tr.powi(3))
                    / 1000.0)
    }

    /// `calculateGibbsEnergy`, kJ/mol.
    ///
    /// **The polynomial branch is taken if the species has *any* coefficient**, so a
    /// species the coefficient file carries never reaches the fallback - including
    /// `oxygen`, `hydrogen` and `nitrogen`, whose six coefficients are all exactly zero and
    /// whose Gibbs energy this therefore returns as exactly `0.0` at every temperature.
    /// The fallback branch is the class's order-3 integral from 298.15 K.
    #[must_use]
    pub fn gibbs_energy(&self, temperature: f64, cp: [f64; 4]) -> f64 {
        if let Some(polynomials) = &self.polynomials {
            return evaluate(&polynomials.gibbs, temperature);
        }
        let [da, db, dc, dd] = self.corrected_heat_capacity_coefficients(cp);
        let t = temperature;
        let delta_gf_rt_ref = self.gf298 / (R_KJ * REFERENCE_TEMPERATURE);
        let delta_gf_rt = delta_gf_rt_ref
            + self.i_term(cp)
            + (1.0 / R_KJ)
                * (self.j_term(cp) / t
                    + (-da * t.ln() - db / 2.0 * t - dc / 6.0 * t.powi(2) - dd / 12.0 * t.powi(3))
                        / 1000.0);
        delta_gf_rt * R_KJ * t
    }

    /// `calculateEnthalpy`, kJ/mol, on the same two branches.
    #[must_use]
    pub fn enthalpy(&self, temperature: f64, cp: [f64; 4]) -> f64 {
        if let Some(polynomials) = &self.polynomials {
            return evaluate(&polynomials.enthalpy, temperature);
        }
        let [da, db, dc, dd] = self.corrected_heat_capacity_coefficients(cp);
        let t = temperature;
        let tr = REFERENCE_TEMPERATURE;
        self.hf298
            + (da * (t - tr)
                + db / 2.0 * (t.powi(2) - tr.powi(2))
                + dc / 3.0 * (t.powi(3) - tr.powi(3))
                + dd / 4.0 * (t.powi(4) - tr.powi(4)))
                / 1000.0
    }

    /// `calculateEntropy`, **kJ/(mol·K)**.
    ///
    /// `cp0` is the component's ideal-gas heat capacity at this temperature - NeqSim's
    /// `getCp0(T)`, which azoth has as `azoth_eos::ideal_gas_cp`. The class's javadoc says
    /// the return is J/(mol·K) while dividing by 1000, and the division is what the
    /// callers see: the method hands back kJ/(mol·K) under a doc comment that says
    /// otherwise. Ported as written.
    #[must_use]
    pub fn entropy(&self, temperature: f64, cp0: f64) -> f64 {
        (self.sf298 + cp0 * (temperature / REFERENCE_TEMPERATURE).ln()) / 1000.0
    }
}

/// `c0*T^5 + c1*T^4 + c2*T^3 + c3*T^2 + c4*T + c5`.
fn evaluate(coefficients: &[f64; 6], temperature: f64) -> f64 {
    let mut total = 0.0;
    for (index, coefficient) in coefficients.iter().enumerate() {
        total += coefficient * temperature.powi(5 - index as i32);
    }
    total
}

/// The class's species database, the two compiled tables joined on the lowercased name.
#[derive(Debug, Clone)]
pub struct GibbsDatabase {
    species: Vec<GibbsSpecies>,
    by_name: BTreeMap<String, usize>,
}

impl GibbsDatabase {
    /// The shipped database.
    ///
    /// # Errors
    /// [`AzothError::InvalidInput`] if either embedded table is malformed, which is a
    /// build defect rather than a caller condition - the tables are compiled and
    /// `gen_databank.py --check` holds them.
    pub fn shipped() -> Result<Self> {
        Self::from_csv(SPECIES_CSV, COEFFICIENTS_CSV)
    }

    /// Parse the two tables and join them.
    ///
    /// **The join is on the lowercased name on both sides**, because the class lowercases
    /// both: `componentMap` is keyed by `molecule.toLowerCase()` and `extraCoeffMap` by
    /// `parts[0].trim().toLowerCase()`. The two files do not agree on case - the
    /// coefficient file writes `co2` where the species file writes `CO2` - so a
    /// case-sensitive join would silently keep 8 of the 15 coefficient rows.
    ///
    /// # Errors
    /// [`AzothError::InvalidInput`] on a short row, an unparseable number, or a species
    /// whose row does not carry eight counts.
    pub fn from_csv(species_csv: &str, coefficients_csv: &str) -> Result<Self> {
        let mut coefficients: BTreeMap<String, Polynomials> = BTreeMap::new();
        for record in parse(coefficients_csv)? {
            let name = record[0].to_lowercase();
            let mut values = [0.0; 12];
            for (index, value) in record[1..13].iter().enumerate() {
                values[index] = number(value, &name)?;
            }
            coefficients.insert(
                name,
                Polynomials {
                    gibbs: values[..6].try_into().expect("six"),
                    enthalpy: values[6..].try_into().expect("six"),
                },
            );
        }

        let mut species = Vec::new();
        let mut by_name = BTreeMap::new();
        for record in parse(species_csv)? {
            let name = record[0].to_string();
            let mut elements = [0.0; 8];
            for (index, value) in record[1..9].iter().enumerate() {
                elements[index] = number(value, &name)?;
            }
            // The species is the file's spelling, so the lookup lowercases the query.
            let polynomials = coefficients.get(&name.to_lowercase()).copied();
            let hf298 = number(&record[9], &name)?;
            let gf298 = number(&record[10], &name)?;
            let sf298 = number(&record[11], &name)?;
            by_name.insert(name.to_lowercase(), species.len());
            species.push(GibbsSpecies {
                name,
                elements,
                hf298,
                gf298,
                sf298,
                polynomials,
            });
        }
        if species.is_empty() {
            return Err(AzothError::invalid_input(
                "gibbs_database",
                "the species table compiled to no rows",
            ));
        }
        Ok(Self { species, by_name })
    }

    /// Every species, in the table's own order.
    #[must_use]
    pub fn species(&self) -> &[GibbsSpecies] {
        &self.species
    }

    /// The species a fluid's component name refers to.
    ///
    /// The lookup is `componentMap.get(componentName.toLowerCase())`, so it is
    /// case-insensitive and a component azoth's databank carries but this table does not
    /// resolves to `None` - which is the class's `comp == null` and is a component the
    /// solve skips rather than refuses.
    #[must_use]
    pub fn get(&self, component: &str) -> Option<&GibbsSpecies> {
        self.by_name
            .get(&component.to_lowercase())
            .map(|index| &self.species[*index])
    }

    /// How many species the database carries.
    #[must_use]
    pub fn len(&self) -> usize {
        self.species.len()
    }

    /// Whether the database is empty, which it never is.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.species.is_empty()
    }
}

/// One table's records, each a vector of fields, header included.
///
/// **`has_headers(false)`, because `records()` skips the header by default and the header
/// is a record here.** Left at the default the iterator yields the 49 data rows of the
/// species table and this function's own `remove(0)` then eats the first of them - which is
/// `CO2`, the row a test looks for first. Reading the header as a record and removing it
/// explicitly is the only reading in which the two steps cannot both be applied.
///
/// The compiled tables are comma-separated and are generated, so a malformed one is a
/// build defect.
fn parse(csv: &str) -> Result<Vec<Vec<String>>> {
    let mut records = Vec::new();
    let mut reader = csv::ReaderBuilder::new()
        .has_headers(false)
        .from_reader(csv.as_bytes());
    for record in reader.records() {
        let record = record.map_err(|error| {
            AzothError::invalid_input(
                "gibbs_database",
                format!("the table could not be read: {error}"),
            )
        })?;
        records.push(record.iter().map(str::to_string).collect());
    }
    if records.is_empty() {
        return Err(AzothError::invalid_input(
            "gibbs_database",
            "the table has no header record",
        ));
    }
    records.remove(0);
    Ok(records)
}

/// One field as a number, naming the row when it is not.
fn number(field: &str, row: &str) -> Result<f64> {
    field.trim().parse::<f64>().map_err(|_| {
        AzothError::invalid_input(
            "gibbs_database",
            format!("the row for {row:?} carries {field:?}, which is not a number"),
        )
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_shipped_table_is_the_size_the_manifest_records() {
        let database = GibbsDatabase::shipped().expect("the shipped tables load");
        assert_eq!(database.len(), 49, "the compiled species table");
        let with_polynomials = database
            .species()
            .iter()
            .filter(|s| s.polynomials.is_some())
            .count();
        assert_eq!(
            with_polynomials, 15,
            "the compiled coefficient table joins to 15"
        );
    }

    /// **The two things the source table gets wrong are reproduced, and this is where
    /// they are pinned.** `argon` carries its atom count in the column the file labels
    /// `Na`; `OH-` carries a charge of `-1` in the eighth; and the class's names for the
    /// last two balances are shifted by one against the file's columns.
    #[test]
    fn the_element_columns_say_what_the_class_reads_them_as() {
        let database = GibbsDatabase::shipped().expect("the shipped tables load");

        let argon = database.get("argon").expect("argon is a row");
        assert_eq!(argon.elements, [0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0]);
        // The class's `Ar` balance - index 5 - is the file's `Na` column, and that is
        // where argon's count is.
        assert_eq!(argon.balance_elements()[5], 1.0);
        // The class's `Z` balance - index 6 - is the file's `Ar` column, which is zero
        // on every row, so that balance row is never an active element.
        assert_eq!(argon.balance_elements()[6], 0.0);

        let hydroxide = database.get("OH-").expect("OH- is a row");
        assert_eq!(hydroxide.elements[7], -1.0, "the eighth column is a charge");
        assert!(
            database.species().iter().all(|s| s.elements[6] == 0.0),
            "the file's `Ar` column is zero on every row"
        );
        assert_eq!(
            database
                .species()
                .iter()
                .filter(|s| s.elements[7] != 0.0)
                .count(),
            9,
            "the nine ions are the only rows carrying a charge"
        );
    }

    /// The species azoth's component table does not carry are not compiled, and the three
    /// the coefficient file names that nothing can reach are not either.
    #[test]
    fn the_unreachable_rows_are_not_in_the_table() {
        let database = GibbsDatabase::shipped().expect("the shipped tables load");
        for name in ["FeS", "S2", "TEA", "HCOO-", "SO3", "NO"] {
            assert!(
                database.get(name).is_none(),
                "{name} should not be compiled"
            );
        }
        // A species NeqSim's own component table carries but azoth's does not.
        assert!(database.get("NO").is_none());
        // And the lookup is case-insensitive, as the class's is.
        assert!(database.get("CO2").is_some());
        assert!(database.get("co2").is_some());
    }

    /// The polynomial branch is exact arithmetic on the file's coefficients, and the
    /// zero-coefficient species return exactly zero - which is the class's behaviour and
    /// not an accident of the port.
    #[test]
    fn the_polynomial_branch_is_the_file_s_coefficients() {
        let database = GibbsDatabase::shipped().expect("the shipped tables load");
        let cp = [0.0; 4];

        let oxygen = database.get("oxygen").expect("oxygen is a row");
        assert!(
            oxygen.polynomials.is_some(),
            "oxygen has a row of zero coefficients"
        );
        assert_eq!(oxygen.gibbs_energy(900.0, cp), 0.0);
        assert_eq!(oxygen.enthalpy(900.0, cp), 0.0);

        // CO2's Gibbs polynomial at the reference temperature is close to its tabulated
        // standard Gibbs energy - the two are independent fits, so this is a band and not
        // an equality. The band is what makes a column misalignment visible.
        let co2 = database.get("CO2").expect("CO2 is a row");
        let at_reference = co2.gibbs_energy(REFERENCE_TEMPERATURE, cp);
        assert!(
            (at_reference - (-394.4)).abs() < 1.0,
            "CO2's polynomial gives {at_reference}, and its tabulated standard Gibbs energy is -394.4"
        );
    }

    /// A species with no polynomial takes the fallback branch, and the branch is the
    /// class's `I` and `J` terms.
    #[test]
    fn the_fallback_branch_uses_the_corrected_heat_capacity() {
        let database = GibbsDatabase::shipped().expect("the shipped tables load");
        let formic_acid = database.get("formic acid").expect("formic acid is a row");
        assert!(formic_acid.polynomials.is_none());

        // Methane's Cp polynomial as NeqSim ships it, which is what the class reads
        // through `getCpA()`..`getCpD()`.
        let cp = [37.978_352, -0.074_618_15, 0.000_301_881, -2.83e-7];
        // HCOOH is O=2, N=0, C=1, H=2, S=0, so the subtraction is
        // A - 2*A_O - 0*A_N - 1*A_C - 2*A_H - 0*A_S.
        assert_eq!(formic_acid.elements[..5], [2.0, 0.0, 1.0, 2.0, 0.0]);
        let corrected = formic_acid.corrected_heat_capacity_coefficients(cp);
        let expected = cp[0] - 2.0 * 12.73 - 8.43 - 2.0 * 14.544;
        assert!(
            (corrected[0] - expected).abs() < 1e-12,
            "expected {expected}, got {}",
            corrected[0]
        );

        // The fallback reproduces the reference state exactly: at 298.15 K the integral
        // terms vanish and the Gibbs energy is the tabulated one.
        //
        // **This does not pin `J`'s sign, and it is worth being clear about that.**
        // `J` enters `I` and the temperature bracket with opposite signs, so the two
        // cancel identically at `T = REFERENCE_TEMPERATURE` whatever `J` is - the
        // assertion holds for every sign convention, including the one the class's
        // javadoc states and the class's code does not. What it does pin is that the two
        // terms are mutually consistent and that the reference state is reproduced at
        // all. `J`'s sign is pinned by the NeqSim capture, which reads a temperature
        // away from 298.15.
        let at_reference = formic_acid.gibbs_energy(REFERENCE_TEMPERATURE, cp);
        assert!(
            (at_reference - formic_acid.gf298).abs() < 1e-9,
            "at the reference temperature the fallback returns gf298, and it returned {at_reference}"
        );
    }
}
