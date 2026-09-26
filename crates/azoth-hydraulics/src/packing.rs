//! NeqSim's packing specifications: the built-in table, the compiled CSV and the lookup.
//!
//! `PackingSpecificationLibrary` and `PackingSpecification`, which `PackingHydraulicsCalculator`
//! takes its geometry from. The table is 22 built-ins - thirteen random and nine structured,
//! written into the class - **and the CSV is loaded after them, so a CSV row wins a name
//! collision**: measured, `Pall-Ring-50` resolves to the file's *plastic* row (111.1 m**2/m**3,
//! 0.919, factor 180) rather than the built-in's metal one (120.0, 0.96, 66), because both
//! normalize to `pallring50` and the later registration replaces the earlier.
//!
//! The file is `data/packing/packing.csv`, compiled from NeqSim's `designdata/Packing.csv` by
//! `tools/gen_databank.py`'s sibling for this table; `databank/manifest.toml` records every one of
//! its columns.

use std::collections::HashMap;
use std::sync::LazyLock;

/// A packing's geometry and material constants.
///
/// The class validates every one of these positivity: a non-positive area, void fraction, packing
/// factor, critical surface tension or Billet constant is refused rather than carried, which is
/// why the constructor is fallible.
#[derive(Debug, Clone, PartialEq)]
pub struct PackingSpecification {
    /// The display name, which for a random packing the library derives from the raw name and the
    /// nominal size - `Pall-Ring-50`, not `pallring`.
    pub name: String,
    /// `random` or `structured`.
    pub category: String,
    /// The material, which the critical surface tension is read from.
    pub material: String,
    /// The nominal size in millimetres; zero for a structured packing, whose rows carry none.
    pub nominal_size_mm: f64,
    /// The specific surface area, in m**2/m**3.
    pub specific_surface_area: f64,
    /// The void fraction.
    pub void_fraction: f64,
    /// The packing factor, in 1/m.
    pub packing_factor: f64,
    /// The critical surface tension of the material, in N/m.
    pub critical_surface_tension: f64,
    /// The Billet-Schultes liquid constant, which the `BILLET_SCHULTES_1999` multiplier reads.
    pub billet_liquid_constant: f64,
    /// The Billet-Schultes gas constant.
    pub billet_gas_constant: f64,
}

/// NeqSim's own alias normalisation: lowercased with every non-alphanumeric removed.
///
/// `Pall-Ring-50`, `pallring 50` and `PALL_RING_50` are one key.
#[must_use]
pub fn normalize(name: &str) -> String {
    name.to_lowercase()
        .chars()
        .filter(char::is_ascii_alphanumeric)
        .collect()
}

/// The critical surface tension a material's packing takes, in N/m.
///
/// `PackingSpecificationLibrary.criticalSurfaceTension`: a plastic packing 0.033, a ceramic one
/// 0.061, and anything else 0.075.
#[must_use]
pub fn critical_surface_tension(material: &str) -> f64 {
    let normalized = material.to_lowercase();
    if normalized.contains("plastic") {
        return 0.033;
    }
    if normalized.contains("ceramic") {
        return 0.061;
    }
    0.075
}

/// The display name a CSV row's raw name, category and size make.
///
/// `PackingSpecificationLibrary.displayName`: a structured packing keeps the raw name, and a
/// random one is named after its family - `pallring` becomes `Pall-Ring-N`, `rashig`/`raschig`/
/// `rachig` become `Raschig-Ring-N`, and anything else keeps its raw name with the rounded size
/// appended.
#[must_use]
pub fn display_name(raw_name: &str, category: &str, size_mm: f64) -> String {
    if category.eq_ignore_ascii_case("structured") {
        return raw_name.to_string();
    }
    let normalized = normalize(raw_name);
    let rounded = size_mm.round() as i64;
    if normalized.contains("pallring") {
        return format!("Pall-Ring-{rounded}");
    }
    if normalized.contains("rashig")
        || normalized.contains("raschig")
        || normalized.contains("rachig")
    {
        return format!("Raschig-Ring-{rounded}");
    }
    format!("{raw_name}-{rounded}")
}

/// A built-in random packing's row: name, material, size in mm, area, void fraction, packing
/// factor, and the two Billet constants.
type RandomRow = (&'static str, &'static str, f64, f64, f64, f64, f64, f64);

/// A built-in structured packing's row, which carries no nominal size.
type StructuredRow = (&'static str, &'static str, f64, f64, f64, f64, f64);

/// NeqSim's built-in random packings, `registerBuiltIns`' thirteen.
///
/// `registerBuiltIns`' thirteen, whose values the class's own comment calls "built-in
/// Kister/Onda engineering estimates" and "vendor-style" ones.
const RANDOM: &[RandomRow] = &[
    ("Pall-Ring-25", "metal", 25.0, 210.0, 0.94, 157.0, 1.0, 0.40),
    ("Pall-Ring-38", "metal", 38.0, 164.0, 0.95, 92.0, 1.0, 0.40),
    ("Pall-Ring-50", "metal", 50.0, 120.0, 0.96, 66.0, 1.0, 0.40),
    (
        "Raschig-Ring-25",
        "ceramic",
        25.0,
        190.0,
        0.68,
        580.0,
        1.0,
        0.40,
    ),
    (
        "Raschig-Ring-50",
        "ceramic",
        50.0,
        95.0,
        0.74,
        155.0,
        1.0,
        0.40,
    ),
    ("IMTP-25", "metal", 25.0, 226.0, 0.97, 134.0, 1.05, 0.42),
    ("IMTP-40", "metal", 40.0, 151.0, 0.97, 79.0, 1.05, 0.42),
    ("IMTP-50", "metal", 50.0, 102.0, 0.98, 56.0, 1.05, 0.42),
    ("IMTP-70", "metal", 70.0, 72.0, 0.98, 36.0, 1.05, 0.42),
    (
        "Berl-Saddle-25",
        "ceramic",
        25.0,
        260.0,
        0.68,
        360.0,
        0.95,
        0.38,
    ),
    (
        "Berl-Saddle-38",
        "ceramic",
        38.0,
        165.0,
        0.70,
        220.0,
        0.95,
        0.38,
    ),
    (
        "Berl-Saddle-50",
        "ceramic",
        50.0,
        105.0,
        0.72,
        150.0,
        0.95,
        0.38,
    ),
    (
        "Intalox-Saddle-25",
        "ceramic",
        25.0,
        255.0,
        0.78,
        200.0,
        1.0,
        0.40,
    ),
];

/// NeqSim's built-in structured packings, `registerBuiltIns`' nine.
///
/// `registerBuiltIns`' nine, whose nominal size is zero - the structured constructor takes none.
const STRUCTURED: &[StructuredRow] = &[
    ("Mellapak-125Y", "metal", 125.0, 0.99, 33.0, 1.2, 0.45),
    ("Mellapak-250Y", "metal", 250.0, 0.98, 66.0, 1.2, 0.45),
    ("Mellapak-350Y", "metal", 350.0, 0.97, 105.0, 1.2, 0.45),
    ("Mellapak-500Y", "metal", 500.0, 0.96, 180.0, 1.2, 0.45),
    ("Flexipac-1Y", "metal", 135.0, 0.99, 36.0, 1.15, 0.44),
    ("Flexipac-2Y", "metal", 220.0, 0.98, 60.0, 1.15, 0.44),
    ("Flexipac-3Y", "metal", 340.0, 0.97, 100.0, 1.15, 0.44),
    ("Sulzer-BX", "metal", 500.0, 0.90, 140.0, 1.3, 0.50),
    ("Sulzer-CY", "metal", 750.0, 0.88, 220.0, 1.3, 0.50),
];

/// The compiled table, which is `data/packing/packing.csv`.
const PACKING_CSV: &str = include_str!("../../../data/packing/packing.csv");

/// The registry, built once: the built-ins first, then the file.
static PACKINGS: LazyLock<HashMap<String, PackingSpecification>> = LazyLock::new(|| {
    let mut table = HashMap::new();
    for (name, material, size, area, void, factor, liquid, gas) in RANDOM {
        register(
            &mut table,
            PackingSpecification {
                name: (*name).to_string(),
                category: "random".to_string(),
                material: (*material).to_string(),
                nominal_size_mm: *size,
                specific_surface_area: *area,
                void_fraction: *void,
                packing_factor: *factor,
                critical_surface_tension: critical_surface_tension(material),
                billet_liquid_constant: *liquid,
                billet_gas_constant: *gas,
            },
        );
    }
    for (name, material, area, void, factor, liquid, gas) in STRUCTURED {
        register(
            &mut table,
            PackingSpecification {
                name: (*name).to_string(),
                category: "structured".to_string(),
                material: (*material).to_string(),
                nominal_size_mm: 0.0,
                specific_surface_area: *area,
                void_fraction: *void,
                packing_factor: *factor,
                critical_surface_tension: critical_surface_tension(material),
                billet_liquid_constant: *liquid,
                billet_gas_constant: *gas,
            },
        );
    }
    load_csv(&mut table);
    table
});

/// Register one specification under its normalized name.
fn register(
    table: &mut HashMap<String, PackingSpecification>,
    specification: PackingSpecification,
) {
    table.insert(normalize(&specification.name), specification);
}

/// Load the compiled CSV **after** the built-ins, which is what lets a file row replace one.
///
/// `loadCsvPackings`, including its two derivations: the Billet liquid constant is the `Cp`
/// column where it is positive and one where it is not, and the gas constant is `Ch` divided by
/// six, floored at 0.4. A row whose area, void fraction or factor is not positive is skipped
/// rather than carried, which is the class's own `parseCsvLine`.
fn load_csv(table: &mut HashMap<String, PackingSpecification>) {
    for line in PACKING_CSV.lines().skip(1) {
        let fields: Vec<&str> = line.split(',').collect();
        if fields.len() < 9 {
            continue;
        }
        let Ok(size) = fields[3].trim().parse::<f64>() else {
            continue;
        };
        let (Ok(area), Ok(void), Ok(factor)) = (
            fields[4].trim().parse::<f64>(),
            fields[5].trim().parse::<f64>(),
            fields[6].trim().parse::<f64>(),
        ) else {
            continue;
        };
        if area <= 0.0 || void <= 0.0 || factor <= 0.0 {
            continue;
        }
        let cp = fields[7].trim().parse::<f64>().unwrap_or(1.0);
        let ch = fields[8].trim().parse::<f64>().unwrap_or(0.4);
        let category = fields[1].trim();
        let material = fields[2].trim();
        register(
            table,
            PackingSpecification {
                name: display_name(fields[0].trim(), category, size),
                category: category.to_string(),
                material: material.to_string(),
                nominal_size_mm: size,
                specific_surface_area: area,
                void_fraction: void,
                packing_factor: factor,
                critical_surface_tension: critical_surface_tension(material),
                billet_liquid_constant: if cp > 0.0 { cp } else { 1.0 },
                billet_gas_constant: if ch > 0.0 { ch / 6.0 } else { 0.4 },
            },
        );
    }
}

/// A packing by name or alias, or `None`.
#[must_use]
pub fn packing(name: &str) -> Option<&'static PackingSpecification> {
    PACKINGS.get(&normalize(name))
}

/// A packing by name, falling back to `Pall-Ring-50`, which is the class's own default.
///
/// # Panics
/// Never in practice: the fallback's own row is a built-in, so it is registered by construction.
#[must_use]
pub fn packing_or_default(name: &str) -> &'static PackingSpecification {
    packing(name)
        .or_else(|| packing("Pall-Ring-50"))
        .expect("Pall-Ring-50 is a built-in packing")
}

/// Every registered packing's name, sorted - the class's own `getPackingNames` reports them in
/// registration order, and a sorted list is the same set without depending on a `HashMap`'s.
#[must_use]
pub fn packing_names() -> Vec<&'static str> {
    let mut names: Vec<&'static str> = PACKINGS.values().map(|p| p.name.as_str()).collect();
    names.sort_unstable();
    names
}
