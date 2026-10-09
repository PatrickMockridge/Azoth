//! `characterization.pedersen_plus_split` - a plus fraction divided into carbon-number cuts.
//!
//! Spec: `specs/models/characterization/pedersen_plus_split.toml`. Oracle:
//! `validation/neqsim/captures/plus_fraction_probe.tsv`.
//!
//! # Two damped Newton solves, not one
//!
//! The cuts' abundances are `z_i = exp(a + b*CN)` and their gravities `rho_i = c + d*ln(CN)`.
//! `PedersenPlusModelSolver.solve` fixes `(a, b)` against the plus fraction's mole fraction and
//! mean molar mass, writes the result back to the model, and only then fixes `(c, d)` against the
//! tabulated anchor gravity and the plus fraction's own. The second solve reads the first's
//! answer, so the two are sequential.
//!
//! Four things in that solver are reproduced rather than repaired, and each is a number:
//!
//! * the two loops damp by different constants - `iter/(iter + 50)` for the abundances and
//!   `iter/(iter + 5)` for the gravities;
//! * neither stops when its residual falls below the tolerance alone: `iter < 3` forces three
//!   steps, and the residual that ends the loop is the one measured **before** the final step;
//! * every abundance iteration also runs an inner fixed point on the gravity slope alone, with
//!   the intercept held at `0.80`;
//! * the gravity solve starts from the constructor's own `(0.80, 0.0408709562)` rather than from
//!   what the abundance phase left on the model.
//!
//! Three Jacobian entries are not the derivatives they stand for and are not corrected here: the
//! abundance row's `(1,0)` is the class's own `mTot1*zSum - mTot1*zSum`, identically zero; its
//! `(1,1)` divides by the sum of squares of `z` where the derivative divides by the square of the
//! sum; and the gravity Jacobian's `(1,1)` is zeroed and never set, while its `(1,0)` carries the
//! opposite sign to the derivative.
//!
//! # The linear solve is Jama's
//!
//! `Jac.solve(fvec)` is a 2x2 LU decomposition with partial pivoting, and the pivot test matters
//! here: the gravity Jacobian's `(1,0)` sits near one against an `(0,0)` of exactly one, so a
//! last-bit difference in it chooses a different elimination order. [`solve_2x2`] transcribes
//! `Jama.Matrix.LUDecomposition` operation for operation rather than replacing it with the closed
//! form, which agrees to twelve places and then stops agreeing.
//!
//! # A range ends at C80's table, not at C200's
//!
//! `PVTsimMolarMass` runs C6 to C200 and `PVTsimDensities` stops at C80, and the class indexes
//! both as `[CN - 6]`. The gravity table is read at exactly one row - the anchor below - so the
//! short table is enough, and it is stored here at its own length so that the asymmetry is
//! visible rather than assumed away.

use std::sync::OnceLock;

use azoth_core::units::{MassDensity, MolarMass, kilograms_per_cubic_meter, kilograms_per_mole};
use azoth_core::{AzothError, Result, apply_checks};

use crate::model_gen;
use crate::results::PedersenPlusSplitResult;

/// The PVTsim tables, compiled from the pinned `PlusFractionModel.java` by `tools/gen_databank.py`.
const PLUS_FRACTION_CSV: &str = include_str!("../../../data/characterization/plus_fraction.csv");

/// The carbon number the tables begin at, which is the class's own `[i - 6]` index base.
const TABLE_FIRST_CARBON_NUMBER: usize = 6;

/// The residual 2-norm both loops stop below, and the inner fixed point's own threshold.
const RESIDUAL_TOLERANCE: f64 = 1.0e-6;

/// The cap on each Newton loop, which the class applies to the abundances and the gravities alike.
const MAX_ITERATIONS: u32 = 200;

/// A loop is not allowed to stop before this many steps, whatever its residual.
const MIN_ITERATIONS: u32 = 3;

/// The abundance step's damping: `iter/(iter + 50)`.
const ABUNDANCE_DAMPING: f64 = 50.0;

/// The gravity step's damping: `iter/(iter + 5)`, a tenth of the other's.
const GRAVITY_DAMPING: f64 = 5.0;

/// The cap on the inner fixed point, five times the outer loops'.
const DENSITY_SLOPE_MAX_ITERATIONS: u32 = 1000;

/// `coefs[0]` when the abundance solve starts, seeded rather than solved.
const Z_SEED_INTERCEPT: f64 = 0.1;

/// `coefs[2]`'s constructor value, which the abundance phase holds fixed and the gravity solve
/// starts from.
const DENSITY_SEED_INTERCEPT: f64 = 0.80;

/// `coefs[3]`'s constructor value, which the gravity solve starts from **regardless** of what the
/// abundance phase's inner fixed point left behind.
const DENSITY_SEED_SLOPE: f64 = 0.040_870_956_2;

/// The two PVTsim tables, as `PedersenPlusModel` declares them.
struct PlusFractionTables {
    /// One molar mass per carbon number from C6, in g/mol. This *is* the cut's molar mass.
    molar_mass: Vec<f64>,
    /// One specific gravity per carbon number from C6, in g/cm3. **Shorter than `molar_mass`**:
    /// NeqSim's `PVTsimDensities` ends at C80 against a molar mass table that reaches C200.
    densities: Vec<f64>,
}

impl PlusFractionTables {
    /// The specific gravity the gravity line is anchored to.
    fn anchor(&self, first_carbon_number: usize) -> Result<f64> {
        let index = first_carbon_number
            .checked_sub(TABLE_FIRST_CARBON_NUMBER)
            .filter(|index| *index < self.densities.len());
        index.map_or_else(
            || {
                Err(AzothError::invalid_input(
                    "first_carbon_number",
                    format!(
                        "the density table runs C{} to C{}; the gravity line is anchored to the \
                         first cut's own row, so a plus fraction outside that span has no anchor",
                        TABLE_FIRST_CARBON_NUMBER,
                        TABLE_FIRST_CARBON_NUMBER + self.densities.len() - 1
                    ),
                ))
            },
            |index| Ok(self.densities[index]),
        )
    }

    /// The cut's tabulated molar mass in g/mol.
    fn cut_molar_mass(&self, carbon_number: usize) -> Result<f64> {
        let index = carbon_number
            .checked_sub(TABLE_FIRST_CARBON_NUMBER)
            .filter(|index| *index < self.molar_mass.len());
        index.map_or_else(
            || {
                Err(AzothError::invalid_input(
                    "last_carbon_number",
                    format!(
                        "the molar-mass table runs C{} to C{}",
                        TABLE_FIRST_CARBON_NUMBER,
                        TABLE_FIRST_CARBON_NUMBER + self.molar_mass.len() - 1
                    ),
                ))
            },
            |index| Ok(self.molar_mass[index]),
        )
    }
}

/// The tables, parsed once from the embedded CSV.
fn tables() -> &'static PlusFractionTables {
    static TABLES: OnceLock<PlusFractionTables> = OnceLock::new();
    TABLES.get_or_init(|| parse_tables().expect("the embedded PVTsim plus-fraction tables should parse"))
}

fn parse_tables() -> std::result::Result<PlusFractionTables, String> {
    let mut reader = csv::ReaderBuilder::new()
        .has_headers(true)
        .from_reader(PLUS_FRACTION_CSV.as_bytes());
    let mut molar_mass = Vec::new();
    let mut densities = Vec::new();
    let mut previous: Option<usize> = None;
    for (row, record) in reader.records().enumerate() {
        let record = record.map_err(|error| format!("row {}: {error}", row + 2))?;
        let field = |index: usize| record.get(index).unwrap_or("").trim().to_string();
        let carbon_number: usize = field(0)
            .parse()
            .map_err(|error| format!("row {}: carbon number: {error}", row + 2))?;
        let expected = previous.map_or(TABLE_FIRST_CARBON_NUMBER, |last| last + 1);
        if carbon_number != expected {
            return Err(format!("row {}: C{carbon_number} after C{expected}", row + 2));
        }
        previous = Some(carbon_number);
        molar_mass.push(
            field(1)
                .parse::<f64>()
                .map_err(|error| format!("row {}: molar mass: {error}", row + 2))?,
        );
        // Past C80 the density cell is empty, which is how the shorter table is carried.
        let density = field(2);
        if !density.is_empty() {
            if densities.len() != carbon_number - TABLE_FIRST_CARBON_NUMBER {
                return Err(format!(
                    "row {}: C{carbon_number}'s density is not contiguous with the rows above it",
                    row + 2
                ));
            }
            densities.push(
                density
                    .parse::<f64>()
                    .map_err(|error| format!("row {}: density: {error}", row + 2))?,
            );
        }
    }
    Ok(PlusFractionTables {
        molar_mass,
        densities,
    })
}

/// `c + d*ln(CN)`, the gravity line, in g/cm3.
fn gravity_line(c: f64, d: f64, carbon_number: usize) -> f64 {
    c + d * (carbon_number as f64).ln()
}

/// The state the two solves carry between steps: the model's four coefficients and the three
/// numbers the plus fraction was described by.
struct Split {
    /// `coefs[0]`, the abundance line's intercept.
    a: f64,
    /// `coefs[1]`, the abundance line's slope per carbon number.
    b: f64,
    /// `coefs[2]`, the gravity line's intercept.
    c: f64,
    /// `coefs[3]`, the gravity line's slope per natural logarithm of the carbon number.
    d: f64,
    /// The first carbon number, inclusive.
    first: usize,
    /// One past the last carbon number.
    last: usize,
    /// The plus fraction's mole fraction.
    z_plus: f64,
    /// The plus fraction's molar mass, mol/kg.
    m_plus: f64,
    /// The plus fraction's specific gravity, g/cm3.
    dens_plus: f64,
}

impl Split {
    /// `exp(a + b*CN)` for one carbon number.
    fn abundance(&self, carbon_number: usize) -> f64 {
        (self.a + self.b * carbon_number as f64).exp()
    }
}

/// `setfvecAB`: the abundance residuals, on the **kg/mol** scale.
fn abundance_residuals(split: &Split, table: &PlusFractionTables) -> Result<[f64; 2]> {
    let mut z_sum = 0.0;
    let mut m_sum = 0.0;
    for i in split.first..split.last {
        let z = split.abundance(i);
        z_sum += z;
        m_sum += z * (table.cut_molar_mass(i)? / 1000.0);
    }
    Ok([z_sum - split.z_plus, m_sum / z_sum - split.m_plus])
}

/// `setJacAB`, including the two entries that are not derivatives.
fn abundance_jacobian(split: &Split, table: &PlusFractionTables) -> Result<[[f64; 2]; 2]> {
    let mut n_tot = 0.0;
    let mut n_tot2 = 0.0;
    for i in split.first..split.last {
        let z = split.abundance(i);
        n_tot += z;
        n_tot2 += (i as f64) * z;
    }
    let mut m_tot1 = 0.0;
    let mut m_tot2 = 0.0;
    let mut z_squares = 0.0;
    let mut z_sum = 0.0;
    let mut z_times_carbon = 0.0;
    for i in split.first..split.last {
        let z = split.abundance(i);
        let m = table.cut_molar_mass(i)? / 1000.0;
        m_tot1 += m * z;
        m_tot2 += (i as f64) * m * z;
        z_squares += z.powf(2.0);
        z_sum += z;
        z_times_carbon += (i as f64) * z;
    }
    // `(1,0)` is the class's own `mTot1*zSum - mTot1*zSum`, identically zero, so the pair is
    // solved as a triangular system. `(1,1)` divides by `z_squares`, the sum of squares of `z`,
    // where the derivative divides by the square of `z_sum`.
    Ok([
        [n_tot, n_tot2],
        [0.0, (m_tot2 * z_sum - m_tot1 * z_times_carbon) / z_squares],
    ])
}

/// `setCoefs(double[])`'s inner `do`-`while`: the gravity slope alone, with `c` held fixed.
///
/// The update is `d * rho_plus / <rho>`, a rescaling rather than a step, and it reads the old
/// slope twice. The class runs it inside every abundance iteration.
fn density_slope_fixed_point(
    split: &Split,
    c: f64,
    mut d: f64,
    table: &PlusFractionTables,
) -> Result<f64> {
    let mut iterations = 0_u32;
    loop {
        iterations += 1;
        let mut m_sum = 0.0;
        let mut dens_sum = 0.0;
        for i in split.first..split.last {
            let z = split.abundance(i);
            let m = table.cut_molar_mass(i)? / 1000.0;
            m_sum += z * m;
            dens_sum += z * m / gravity_line(c, d, i);
        }
        dens_sum = m_sum / dens_sum;
        d += 1.0 * (split.dens_plus - dens_sum) / dens_sum * d;
        if !((split.dens_plus - dens_sum).abs() > RESIDUAL_TOLERANCE
            && iterations < DENSITY_SLOPE_MAX_ITERATIONS)
        {
            return Ok(d);
        }
    }
}

/// `setfvecCD`: the gravity residuals, on the **g/mol** scale the class uses here.
fn gravity_residuals(
    split: &Split,
    anchor: f64,
    table: &PlusFractionTables,
) -> Result<[f64; 2]> {
    let mut temp = 0.0;
    let mut temp2 = 0.0;
    for i in split.first..split.last {
        let weight = split.abundance(i) * table.cut_molar_mass(i)?;
        temp += weight;
        temp2 += weight / gravity_line(split.c, split.d, i);
    }
    // The anchor is the carbon number **below** the first cut against the first cut's own row.
    Ok([
        gravity_line(split.c, split.d, split.first - 1) - anchor,
        temp / temp2 - split.dens_plus,
    ])
}

/// `setJacCD`, with its `(1,1)` left at zero and its `(1,0)` at the negative of the derivative.
fn gravity_jacobian(split: &Split, table: &PlusFractionTables) -> Result<[[f64; 2]; 2]> {
    let mut temp = 0.0;
    let mut temp2 = 0.0;
    let mut temp3 = 0.0;
    for i in split.first..split.last {
        let weight = split.abundance(i) * table.cut_molar_mass(i)?;
        let rho = gravity_line(split.c, split.d, i);
        temp += weight;
        temp2 += weight / rho;
        temp3 -= weight / rho.powf(2.0);
    }
    let derivative = temp / (temp2 * temp2) * temp3;
    Ok([[1.0, ((split.first - 1) as f64).ln()], [derivative, 0.0]])
}

/// `Jama.Matrix.LUDecomposition`'s solve, transcribed for a 2x2 system. `None` when singular.
///
/// Written out rather than closed-formed because the pivot test is on a quantity that sits near
/// one, and a closed form would silently commit to the other elimination order.
fn solve_2x2(matrix: [[f64; 2]; 2], rhs: [f64; 2]) -> Option<[f64; 2]> {
    let mut lu = matrix;
    let mut pivot = [0_usize, 1];
    for j in 0..2 {
        let mut column = [lu[0][j], lu[1][j]];
        for i in 0..2 {
            let mut sum = 0.0;
            for k in 0..i.min(j) {
                sum += lu[i][k] * column[k];
            }
            column[i] -= sum;
            lu[i][j] = column[i];
        }
        let mut p = j;
        for i in (j + 1)..2 {
            if column[i].abs() > column[p].abs() {
                p = i;
            }
        }
        if p != j {
            lu.swap(p, j);
            pivot.swap(p, j);
        }
        if lu[j][j] != 0.0 {
            for i in (j + 1)..2 {
                lu[i][j] /= lu[j][j];
            }
        }
    }
    if lu[0][0] == 0.0 || lu[1][1] == 0.0 {
        return None;
    }
    let mut x = [rhs[pivot[0]], rhs[pivot[1]]];
    for k in 0..2 {
        for i in (k + 1)..2 {
            x[i] -= x[k] * lu[i][k];
        }
    }
    for k in (0..2).rev() {
        x[k] /= lu[k][k];
        for i in 0..k {
            x[i] -= x[k] * lu[i][k];
        }
    }
    Some(x)
}

/// A plus fraction divided into carbon-number cuts by Pedersen's two solves.
///
/// # Errors
/// * [`AzothError::OutOfRange`] if an input is outside its declared range, which is where a
///   carbon number below C6, above C80 for the first cut, or above C200 for the last is refused.
/// * [`AzothError::InvalidInput`] if the plus fraction is lighter than the first cut the table
///   gives it, if its range is empty, or if a Jacobian comes out singular.
///
/// # Example
/// ```
/// use azoth_core::units::{kilograms_per_cubic_meter, kilograms_per_mole};
/// use azoth_characterization::pedersen_plus_split;
///
/// let split = pedersen_plus_split(
///     kilograms_per_mole(0.4),
///     kilograms_per_cubic_meter(850.0),
///     0.1,
///     20,
///     80,
/// )?;
/// assert_eq!(split.cut_z.len(), 60);
/// assert!((split.z_intercept - -2.48637183497003).abs() < 1e-12);
/// assert!((split.cut_molar_mass[0].value - 0.275).abs() < 1e-15);
/// # Ok::<(), azoth_core::AzothError>(())
/// ```
pub fn pedersen_plus_split(
    molar_mass: MolarMass,
    density: MassDensity,
    mole_fraction: f64,
    first_carbon_number: usize,
    last_carbon_number: usize,
) -> Result<PedersenPlusSplitResult> {
    let spec = &model_gen::PEDERSEN_PLUS_SPLIT_SPEC;
    let mut warnings = Vec::new();

    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "molar_mass" => Some(molar_mass.value),
            "density" => Some(density.value),
            "mole_fraction" => Some(mole_fraction),
            "first_carbon_number" => Some(first_carbon_number as f64),
            "last_carbon_number" => Some(last_carbon_number as f64),
            _ => None,
        },
        &mut warnings,
    )?;

    if last_carbon_number <= first_carbon_number {
        return Err(AzothError::invalid_input(
            "last_carbon_number",
            format!(
                "{last_carbon_number} leaves no cut above {first_carbon_number}: the range is \
                 half-open, so at least one carbon number has to lie between them"
            ),
        ));
    }

    let table = tables();
    // `characterizePlusFraction`'s own test: the table's first cut against the plus fraction's
    // molar mass in g/mol. NeqSim logs and returns without a split; this refuses, because a
    // caller who asked for a split that was not made has to be told.
    let first_cut_gmol = table.cut_molar_mass(first_carbon_number)?;
    if first_cut_gmol > molar_mass.value * 1000.0 {
        return Err(AzothError::invalid_input(
            "molar_mass",
            format!(
                "the plus fraction is lighter than the first cut the table gives it, at {} g/mol \
                 against C{first_carbon_number}'s {first_cut_gmol} g/mol",
                molar_mass.value * 1000.0
            ),
        ));
    }
    let anchor = table.anchor(first_carbon_number)?;

    let mut split = Split {
        // The class seeds the abundance intercept, and takes the slope from the mole fraction.
        a: Z_SEED_INTERCEPT,
        b: mole_fraction.ln() / first_carbon_number as f64,
        c: DENSITY_SEED_INTERCEPT,
        d: DENSITY_SEED_SLOPE,
        first: first_carbon_number,
        last: last_carbon_number,
        z_plus: mole_fraction,
        m_plus: molar_mass.value,
        dens_plus: density.value / 1000.0,
    };

    // The abundance solve. Each step is damped, published to the model, and followed by the
    // inner fixed point on the gravity slope that the class runs from inside `setCoefs`.
    let mut seed = [split.a, split.b];
    let mut iterations = 0_u32;
    loop {
        iterations += 1;
        let f = abundance_residuals(&split, table)?;
        let jacobian = abundance_jacobian(&split, table)?;
        let dx = solve_2x2(jacobian, f).ok_or_else(|| singular("molar_mass"))?;
        let scale = f64::from(iterations) / (f64::from(iterations) + ABUNDANCE_DAMPING);
        seed[0] -= dx[0] * scale;
        seed[1] -= dx[1] * scale;
        split.a = seed[0];
        split.b = seed[1];
        split.d = density_slope_fixed_point(&split, split.c, split.d, table)?;
        let norm = (f[0] * f[0] + f[1] * f[1]).sqrt();
        if !((norm > RESIDUAL_TOLERANCE || iterations < MIN_ITERATIONS)
            && iterations < MAX_ITERATIONS)
        {
            break;
        }
    }

    // The gravity solve, from the solver's own seed rather than from what the walk above left.
    let mut seed = [DENSITY_SEED_INTERCEPT, DENSITY_SEED_SLOPE];
    let mut iterations = 0_u32;
    loop {
        iterations += 1;
        let f = gravity_residuals(&split, anchor, table)?;
        let jacobian = gravity_jacobian(&split, table)?;
        let dx = solve_2x2(jacobian, f).ok_or_else(|| singular("density"))?;
        let scale = f64::from(iterations) / (f64::from(iterations) + GRAVITY_DAMPING);
        seed[0] -= dx[0] * scale;
        seed[1] -= dx[1] * scale;
        split.c = seed[0];
        split.d = seed[1];
        let norm = (f[0] * f[0] + f[1] * f[1]).sqrt();
        if !((norm > RESIDUAL_TOLERANCE || iterations < MIN_ITERATIONS)
            && iterations < MAX_ITERATIONS)
        {
            break;
        }
    }

    let count = last_carbon_number - first_carbon_number;
    let mut cut_z = Vec::with_capacity(count);
    let mut cut_molar_mass = Vec::with_capacity(count);
    let mut cut_density = Vec::with_capacity(count);
    for i in first_carbon_number..last_carbon_number {
        cut_z.push(split.abundance(i));
        cut_molar_mass.push(kilograms_per_mole(table.cut_molar_mass(i)? / 1000.0));
        cut_density.push(kilograms_per_cubic_meter(
            gravity_line(split.c, split.d, i) * 1000.0,
        ));
    }

    Ok(PedersenPlusSplitResult {
        cut_z,
        cut_molar_mass,
        cut_density,
        z_intercept: split.a,
        z_slope: split.b,
        density_intercept: split.c,
        density_slope: split.d,
        warnings,
    })
}

/// NeqSim's `LUDecomposition` throws `Matrix is singular.` and the caller does not catch it; this
/// refuses instead of returning a division by zero, naming the solve that failed.
fn singular(field: &str) -> AzothError {
    AzothError::invalid_input(
        field,
        "the solve's Jacobian is singular at the plus fraction's own numbers, so this split has \
         no step to take",
    )
}
