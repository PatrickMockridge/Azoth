//! NeqSim's packing specifications: the built-in table, the compiled CSV and the lookup.
//!
//! `PackingSpecificationLibrary` and `PackingSpecification`, which `PackingHydraulicsCalculator`
//! takes its geometry from. The table is 22 built-ins - thirteen random and nine structured,
//! written into the class - registered under their name and a `name-material` alias. The CSV then
//! adds **material-qualified** entries and never replaces a built-in (upstream PR #4405, #4380):
//! `Pall-Ring-50` is the built-in metal estimate (120.0 m**2/m**3, 0.96, factor 66),
//! `Pall-Ring-50-metal` the same, and `Pall-Ring-50-plastic` the file's plastic row (111.1 m**2/m**3,
//! 0.919, factor 180). Before the fix the file was keyed by display name alone, so the plastic row's
//! later line overwrote the built-in's geometry.
//!
//! A row whose category is not `random` or `structured` is not a packing at all - the file's
//! `placeholder` rows are stand-ins - and an unknown name is refused rather than defaulted.
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
///
/// `display_names` is `PackingSpecificationLibrary.DISPLAY_NAMES` - the names `register` was called
/// with, in registration order, so the aliases are keys without being listed.
struct Registry {
    by_key: HashMap<String, PackingSpecification>,
    display_names: Vec<String>,
}

static PACKINGS: LazyLock<Registry> = LazyLock::new(|| {
    let mut registry = Registry {
        by_key: HashMap::new(),
        display_names: Vec::new(),
    };
    let built_in = |name: &str,
                    material: &str,
                    category: &str,
                    size: f64,
                    area: f64,
                    void: f64,
                    factor: f64,
                    liquid: f64,
                    gas: f64| {
        PackingSpecification {
            name: name.to_string(),
            category: category.to_string(),
            material: material.to_string(),
            nominal_size_mm: size,
            specific_surface_area: area,
            void_fraction: void,
            packing_factor: factor,
            critical_surface_tension: critical_surface_tension(material),
            billet_liquid_constant: liquid,
            billet_gas_constant: gas,
        }
    };
    for (name, material, size, area, void, factor, liquid, gas) in RANDOM {
        let specification = built_in(
            name, material, "random", *size, *area, *void, *factor, *liquid, *gas,
        );
        register(&mut registry, specification.clone());
        register_alias(&mut registry, &format!("{name}-{material}"), &specification);
    }
    for (name, material, area, void, factor, liquid, gas) in STRUCTURED {
        let specification = built_in(
            name,
            material,
            "structured",
            0.0,
            *area,
            *void,
            *factor,
            *liquid,
            *gas,
        );
        register(&mut registry, specification.clone());
        register_alias(&mut registry, &format!("{name}-{material}"), &specification);
    }
    load_csv(&mut registry);
    registry
});

/// Register one specification under its normalized name, and list it as a display name.
fn register(registry: &mut Registry, specification: PackingSpecification) {
    let name = specification.name.clone();
    register_alias(registry, &name, &specification);
    registry.display_names.push(name);
}

/// Register an extra alias for an existing specification.
fn register_alias(registry: &mut Registry, alias: &str, specification: &PackingSpecification) {
    registry
        .by_key
        .insert(normalize(alias), specification.clone());
}

/// Whether a name or alias resolves.
fn contains(registry: &Registry, name: &str) -> bool {
    registry.by_key.contains_key(&normalize(name))
}

/// Load the compiled CSV **after** the built-ins, **without replacing any of them**.
///
/// `loadCsvPackings` as upstream PR #4405 left it: the display name carries the material -
/// `Pall-Ring-50-plastic` - and a row whose qualified name is already registered is skipped, so
/// the file's metal 50 mm ring (`112.6`, `0.951`, factor `85`) does not displace the built-in
/// metal estimate (`120.0`, `0.96`, `66`). A qualified name that is new is also registered under
/// its unqualified form, unless that name already exists - which is how `Pall-Ring-35` and
/// `Pall-Ring-15` are reachable and `Pall-Ring-50` stays the built-in.
///
/// The two derivations are the class's: the Billet liquid constant is the `Cp` column where it is
/// positive and one where it is not, and the gas constant is `Ch` divided by six, floored at 0.4.
/// A row whose category is neither `random` nor `structured` - the file's four `placeholder`
/// rows - is not a packing, and a row whose area, void fraction or factor is not positive is
/// skipped rather than carried; both are `parseCsvLine`'s own checks.
fn load_csv(registry: &mut Registry) {
    for line in PACKING_CSV.lines().skip(1) {
        let fields: Vec<&str> = line.split(',').collect();
        if fields.len() < 9 {
            continue;
        }
        let category = fields[1].trim();
        if !category.eq_ignore_ascii_case("random") && !category.eq_ignore_ascii_case("structured")
        {
            continue;
        }
        let material = fields[2].trim().to_lowercase();
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
        let name = format!(
            "{}-{material}",
            display_name(fields[0].trim(), category, size)
        );
        if contains(registry, &name) {
            continue;
        }
        let specification = PackingSpecification {
            name: name.clone(),
            category: category.to_string(),
            material: material.clone(),
            nominal_size_mm: size,
            specific_surface_area: area,
            void_fraction: void,
            packing_factor: factor,
            critical_surface_tension: critical_surface_tension(&material),
            billet_liquid_constant: if cp > 0.0 { cp } else { 1.0 },
            billet_gas_constant: if ch > 0.0 { ch / 6.0 } else { 0.4 },
        };
        register(registry, specification.clone());
        let suffix = format!("-{material}");
        let legacy = name.strip_suffix(&suffix).unwrap_or(&name);
        if !contains(registry, legacy) {
            register_alias(registry, legacy, &specification);
        }
    }
}

/// A packing by name or alias, or `None`.
#[must_use]
pub fn packing(name: &str) -> Option<&'static PackingSpecification> {
    PACKINGS.by_key.get(&normalize(name))
}

/// A packing by name, or `None` - **the class's own refusal**, since upstream PR #4405.
///
/// NeqSim's `getOrDefault` used to answer an unknown preset with `Pall-Ring-50`; #4380 replaced
/// it with `require`, whose `IllegalArgumentException` this reports as `None` and the caller
/// turns into an `InvalidInput`. A default would make a misspelt packing a silent geometry
/// change, which is the defect the fix is about.
#[must_use]
pub fn required_packing(name: &str) -> Option<&'static PackingSpecification> {
    packing(name)
}

/// Every registered packing's display name, in registration order - `getPackingNames`' own list.
///
/// Aliases are lookup keys and are not listed; `Pall-Ring-50-plastic` is a display name because
/// the file registered it, while a built-in's `Pall-Ring-50-metal` is an alias.
#[must_use]
pub fn packing_names() -> Vec<String> {
    PACKINGS.display_names.clone()
}
