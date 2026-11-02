//! `process.separator_performance` - NeqSim's droplet-based separation calculator.
//!
//! **The standard path of `SeparatorPerformanceCalculator`**, with the three classes it reads:
//! the Rosin-Rammler/log-normal `DropletSizeDistribution`, the `GradeEfficiencyCurve` family,
//! and `DropletSettlingCalculator`'s Schiller-Naumann settling. The enhanced chain in the
//! class - flow regime, inlet device, vessel geometry and the internals database - is a
//! different solver and is not here.
//!
//! **The class's quirks are reproduced rather than corrected**: the three-phase gas
//! entrainment writes the whole `1 - overall` to both liquid fractions, the gas carry-under
//! writes one fraction to both liquid phases, and the liquid-liquid stage uses a fixed 300 s
//! residence time. The efficiency integrals are the midpoint sums over 50 log-spaced classes
//! without renormalising, and `erf`/`inverseErf` are the class's own approximations.

use azoth_core::{AzothError, Result};

/// `DropletSettlingCalculator.G`.
const G: f64 = 9.81;
/// `DropletSettlingCalculator.MAX_ITERATIONS`.
const MAX_ITERATIONS: usize = 50;
/// `DropletSettlingCalculator.TOLERANCE`.
const TOLERANCE: f64 = 1.0e-6;
/// `DropletSizeDistribution.numberOfClasses`.
const NUMBER_OF_CLASSES: usize = 50;
/// `calcLiquidLiquidSeparation`'s fixed residence time, s.
const LIQUID_LIQUID_RESIDENCE_TIME: f64 = 300.0;
/// `estimateLiquidResidenceTime`'s velocity floor, m/s.
const MIN_LIQUID_VELOCITY: f64 = 0.0001;
/// The near-dry level floor in `estimateLiquidResidenceTime`.
const MIN_PHASE_FRACTION: f64 = 0.01;
/// `DropletSettlingCalculator.calcCriticalDiameter`'s `Double.MAX_VALUE`: the class's sentinel
/// for a section with no residence time or no density difference. The declared bounds keep the
/// branch unreachable; the name is what says so.
// numerics-ok: the class writes `Double.MAX_VALUE` itself.
const NO_SETTLING_CUT: f64 = f64::MAX;

/// `GradeEfficiencyCurve.wireMeshDefault`'s d_50, m.
const WIRE_MESH_CUT: f64 = 5.0e-6;
/// `GradeEfficiencyCurve.wireMeshDefault`'s sharpness.
const WIRE_MESH_SHARPNESS: f64 = 2.5;
/// `GradeEfficiencyCurve.wireMeshDefault`'s maximum efficiency.
const WIRE_MESH_MAX: f64 = 0.998;
/// `GradeEfficiencyCurve.vanePackDefault`'s d_50, m.
const VANE_PACK_CUT: f64 = 12.0e-6;
/// `GradeEfficiencyCurve.vanePackDefault`'s sharpness.
const VANE_PACK_SHARPNESS: f64 = 2.0;
/// `GradeEfficiencyCurve.vanePackDefault`'s maximum efficiency.
const VANE_PACK_MAX: f64 = 0.995;
/// `GradeEfficiencyCurve.axialCycloneDefault`'s d_50, m.
const AXIAL_CYCLONE_CUT: f64 = 3.0e-6;
/// `GradeEfficiencyCurve.axialCycloneDefault`'s sharpness.
const AXIAL_CYCLONE_SHARPNESS: f64 = 4.0;
/// `GradeEfficiencyCurve.axialCycloneDefault`'s maximum efficiency.
const AXIAL_CYCLONE_MAX: f64 = 0.998;
/// `GradeEfficiencyCurve.platePack`'s sharpness.
const PLATE_PACK_SHARPNESS: f64 = 2.0;

/// The class's `DistributionType`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DistributionKind {
    /// Rosin-Rammler (Weibull).
    RosinRammler,
    /// Log-normal.
    LogNormal,
}

/// The class's `DistributionType` with its parameters.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DropletDistribution {
    /// Which closed form.
    pub kind: DistributionKind,
    /// `d_63.2` for Rosin-Rammler, `d_50` for log-normal.
    pub characteristic_diameter: f64,
    /// The Rosin-Rammler `q` or the log-normal geometric standard deviation.
    pub spread: f64,
}

impl DropletDistribution {
    /// `cumulativeFraction`: the volume fraction smaller than `d`.
    fn cumulative(&self, d: f64) -> f64 {
        if d <= 0.0 {
            return 0.0;
        }
        match self.kind {
            DistributionKind::RosinRammler => {
                1.0 - (-(d / self.characteristic_diameter).powf(self.spread)).exp()
            }
            DistributionKind::LogNormal => {
                let z = (d.ln() - self.characteristic_diameter.ln()) / self.spread;
                0.5 * (1.0 + erf(z / std::f64::consts::SQRT_2))
            }
        }
    }

    /// `inverseCDF`: the diameter at cumulative fraction `f`.
    fn inverse(&self, f: f64) -> f64 {
        if f <= 0.0 {
            return 0.0;
        }
        if f >= 1.0 {
            return self.characteristic_diameter * 10.0;
        }
        match self.kind {
            DistributionKind::RosinRammler => {
                self.characteristic_diameter * (-(1.0 - f).ln()).powf(1.0 / self.spread)
            }
            DistributionKind::LogNormal => (self.characteristic_diameter.ln()
                + self.spread * std::f64::consts::SQRT_2 * inverse_erf(2.0 * f - 1.0))
            .exp(),
        }
    }

    /// `getDiscreteClasses`: `(lower, midpoint, volume fraction)` over the 0.1%-99.9% span.
    fn classes(&self) -> [(f64, f64, f64); NUMBER_OF_CLASSES] {
        let mut d_min = self.inverse(0.001);
        if d_min <= 0.0 {
            d_min = self.characteristic_diameter * 0.01;
        }
        let d_max = self.inverse(0.999);
        let log_min = d_min.ln();
        let log_max = d_max.ln();
        let step = (log_max - log_min) / NUMBER_OF_CLASSES as f64;
        let mut classes = [(0.0, 0.0, 0.0); NUMBER_OF_CLASSES];
        for (index, class) in classes.iter_mut().enumerate() {
            let lower = (log_min + index as f64 * step).exp();
            let upper = (log_min + (index + 1) as f64 * step).exp();
            let middle = (log_min + (index as f64 + 0.5) * step).exp();
            *class = (
                lower,
                middle,
                self.cumulative(upper) - self.cumulative(lower),
            );
        }
        classes
    }
}

/// The class's `InternalsType`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Internals {
    /// The gravity settling section itself, never a configured internals curve.
    Gravity,
    /// Wire mesh demister.
    WireMesh,
    /// Vane pack.
    VanePack,
    /// Axial cyclone.
    AxialCyclone,
    /// Plate pack.
    PlatePack,
    /// A user lookup table.
    Custom,
}

/// One row of a lookup table, in SI.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LookupPoint {
    /// The droplet diameter, m.
    pub diameter: f64,
    /// The capture efficiency.
    pub efficiency: f64,
}

/// One configured internals curve, before its defaults are resolved.
#[derive(Debug, Clone, PartialEq)]
pub struct CurveSettings {
    /// Which correlation.
    pub internals: Internals,
    /// `d_50`; absent takes the type's default.
    pub cut_diameter: Option<f64>,
    /// The sigmoid's sharpness; absent takes the type's default.
    pub sharpness: Option<f64>,
    /// The sigmoid's `eta_max`; absent takes the type's default.
    pub max_efficiency: Option<f64>,
    /// A custom curve's lookup table.
    pub points: Option<Vec<LookupPoint>>,
}

/// The class's `GradeEfficiencyCurve`, resolved.
#[derive(Debug, Clone, PartialEq)]
struct GradeCurve {
    internals: Internals,
    cut_diameter: f64,
    sharpness: f64,
    max_efficiency: f64,
    min_efficiency: f64,
    points: Vec<LookupPoint>,
}

impl GradeCurve {
    /// `GradeEfficiencyCurve.gravity`.
    fn gravity(cut_diameter: f64) -> Self {
        Self {
            internals: Internals::Gravity,
            cut_diameter,
            sharpness: 2.0,
            max_efficiency: 1.0,
            min_efficiency: 0.0,
            points: Vec::new(),
        }
    }

    /// Resolve one declared curve, with the class's per-type defaults.
    fn resolve(role: &str, settings: &CurveSettings) -> Result<Self> {
        if settings.internals == Internals::Gravity {
            return Err(AzothError::invalid_input(
                role,
                "gravity is the settling section, not a configured internals curve",
            ));
        }
        let required = |field: &str| -> Result<f64> {
            settings.cut_diameter.ok_or_else(|| {
                AzothError::invalid_input(
                    field,
                    format!("{role} = {} needs a cut diameter", name(settings.internals)),
                )
            })
        };
        let required_max = |field: &str| -> Result<f64> {
            settings.max_efficiency.ok_or_else(|| {
                AzothError::invalid_input(
                    field,
                    format!(
                        "{role} = {} needs a maximum efficiency",
                        name(settings.internals)
                    ),
                )
            })
        };
        let (cut_diameter, sharpness, max_efficiency) = match settings.internals {
            Internals::Gravity => (0.0, 2.0, 1.0),
            Internals::WireMesh => (
                settings.cut_diameter.unwrap_or(WIRE_MESH_CUT),
                settings.sharpness.unwrap_or(WIRE_MESH_SHARPNESS),
                settings.max_efficiency.unwrap_or(WIRE_MESH_MAX),
            ),
            Internals::VanePack => (
                settings.cut_diameter.unwrap_or(VANE_PACK_CUT),
                settings.sharpness.unwrap_or(VANE_PACK_SHARPNESS),
                settings.max_efficiency.unwrap_or(VANE_PACK_MAX),
            ),
            Internals::AxialCyclone => (
                settings.cut_diameter.unwrap_or(AXIAL_CYCLONE_CUT),
                settings.sharpness.unwrap_or(AXIAL_CYCLONE_SHARPNESS),
                settings.max_efficiency.unwrap_or(AXIAL_CYCLONE_MAX),
            ),
            Internals::PlatePack => (
                required("cut_diameter")?,
                settings.sharpness.unwrap_or(PLATE_PACK_SHARPNESS),
                required_max("max_efficiency")?,
            ),
            Internals::Custom => (0.0, 2.0, 1.0),
        };
        let points = if settings.internals == Internals::Custom {
            let points = settings.points.clone().ok_or_else(|| {
                AzothError::invalid_input(
                    "points",
                    format!("{role} = custom needs a lookup table of at least two points"),
                )
            })?;
            if points.len() < 2 {
                return Err(AzothError::invalid_input(
                    "points",
                    format!(
                        "{role} = custom needs at least two points, not {}",
                        points.len()
                    ),
                ));
            }
            points
                .into_iter()
                .map(|point| LookupPoint {
                    diameter: point.diameter,
                    efficiency: point.efficiency.clamp(0.0, 1.0),
                })
                .collect()
        } else {
            if settings.points.is_some() {
                return Err(AzothError::invalid_input(
                    "points",
                    format!(
                        "{role} = {} does not read a lookup table",
                        name(settings.internals)
                    ),
                ));
            }
            Vec::new()
        };
        Ok(Self {
            internals: settings.internals,
            cut_diameter,
            sharpness,
            max_efficiency,
            min_efficiency: 0.0,
            points,
        })
    }

    /// `getEfficiency`: the fractional capture at one diameter.
    fn efficiency(&self, diameter: f64) -> f64 {
        if diameter <= 0.0 {
            return self.min_efficiency;
        }
        if self.internals == Internals::Custom {
            return interpolate_custom(&self.points, diameter);
        }
        let ratio = diameter / self.cut_diameter;
        let eta = if self.internals == Internals::Gravity {
            ratio * ratio
        } else {
            self.max_efficiency * (1.0 - (-0.693 * ratio.powf(self.sharpness)).exp())
        };
        eta.min(self.max_efficiency).max(self.min_efficiency)
    }

    /// `calcOverallEfficiency`: the class's midpoint sum, not renormalised.
    fn overall(&self, distribution: &DropletDistribution) -> f64 {
        sum_classes(&distribution.classes(), |middle| self.efficiency(middle))
    }
}

/// One efficiency's sum over a class list, clamped to `[0, 1]`.
fn sum_classes(
    classes: &[(f64, f64, f64); NUMBER_OF_CLASSES],
    capture: impl Fn(f64) -> f64,
) -> f64 {
    let total: f64 = classes
        .iter()
        .map(|(_, middle, fraction)| capture(*middle) * fraction)
        .sum();
    // numerics-ok: `calcOverallEfficiency` is `Math.max(0.0, Math.min(1.0, sum))`.
    total.clamp(0.0, 1.0)
}

/// `calcCombinedEfficiency`: the two stages in series, per size class.
fn combined(distribution: &DropletDistribution, first: &GradeCurve, second: &GradeCurve) -> f64 {
    sum_classes(&distribution.classes(), |middle| {
        1.0 - (1.0 - first.efficiency(middle)) * (1.0 - second.efficiency(middle))
    })
}

/// The class's `interpolateCustom`.
fn interpolate_custom(points: &[LookupPoint], diameter: f64) -> f64 {
    if points.len() < 2 {
        return 0.0;
    }
    if diameter <= points[0].diameter {
        return points[0].efficiency;
    }
    let last = points[points.len() - 1];
    if diameter >= last.diameter {
        return last.efficiency;
    }
    for window in points.windows(2) {
        if diameter >= window[0].diameter && diameter < window[1].diameter {
            let fraction =
                (diameter - window[0].diameter) / (window[1].diameter - window[0].diameter);
            return window[0].efficiency + fraction * (window[1].efficiency - window[0].efficiency);
        }
    }
    last.efficiency
}

/// The class's `erf`: Abramowitz and Stegun 7.1.26.
fn erf(x: f64) -> f64 {
    let negative = x < 0.0;
    let abs_x = x.abs();
    let t = 1.0 / (1.0 + 0.3275911 * abs_x);
    let poly = t
        * (0.254829592
            + t * (-0.284496736 + t * (1.421413741 + t * (-1.453152027 + t * 1.061405429))));
    let result = 1.0 - poly * (-abs_x * abs_x).exp();
    if negative { -result } else { result }
}

/// The class's `inverseErf`: Winitzki's rational approximation.
fn inverse_erf(x: f64) -> f64 {
    if x <= -1.0 {
        return -6.0;
    }
    if x >= 1.0 {
        return 6.0;
    }
    let a = 0.147;
    let ln_term = (1.0 - x * x).ln();
    let pi_inv = 2.0 / (std::f64::consts::PI * a);
    let half_ln = ln_term / 2.0;
    let sqrt_term = ((pi_inv + half_ln).powi(2) - ln_term / a).sqrt();
    let result = (sqrt_term - pi_inv - half_ln).sqrt();
    if x >= 0.0 { result } else { -result }
}

/// `calcDragCoefficient`: Schiller-Naumann below `Re = 1000`, Newton above.
fn drag_coefficient(reynolds: f64) -> f64 {
    if reynolds <= 0.0 {
        return 1.0e10;
    }
    if reynolds < 1000.0 {
        (24.0 / reynolds) * (1.0 + 0.15 * reynolds.powf(0.687))
    } else {
        0.44
    }
}

/// `calcReynolds`.
fn reynolds(diameter: f64, velocity: f64, density: f64, viscosity: f64) -> f64 {
    if viscosity <= 0.0 {
        return 0.0;
    }
    density * velocity.abs() * diameter / viscosity
}

/// `calcTerminalVelocity`: positive when the dispersed phase is heavier.
fn terminal_velocity(diameter: f64, continuous: f64, dispersed: f64, viscosity: f64) -> f64 {
    if diameter <= 0.0 || continuous <= 0.0 || viscosity <= 0.0 {
        return 0.0;
    }
    let delta = (dispersed - continuous).abs();
    if delta < 1.0e-6 {
        return 0.0;
    }
    let mut velocity = diameter * diameter * delta * G / (18.0 * viscosity);
    for _ in 0..MAX_ITERATIONS {
        let re = reynolds(diameter, velocity, continuous, viscosity);
        let coefficient = drag_coefficient(re);
        let updated = (4.0 * G * diameter * delta / (3.0 * coefficient * continuous)).sqrt();
        if (updated - velocity).abs() / (velocity + 1.0e-30) < TOLERANCE {
            velocity = updated;
            break;
        }
        velocity = 0.5 * (velocity + updated);
    }
    if dispersed > continuous {
        velocity
    } else {
        -velocity
    }
}

/// `calcCriticalDiameter`: the smallest droplet that settles the height in the time.
fn critical_diameter(
    available_height: f64,
    residence_time: f64,
    continuous: f64,
    dispersed: f64,
    viscosity: f64,
) -> f64 {
    if residence_time <= 0.0 || available_height <= 0.0 {
        return NO_SETTLING_CUT;
    }
    let required = available_height / residence_time;
    let delta = (dispersed - continuous).abs();
    if delta < 1.0e-6 {
        return NO_SETTLING_CUT;
    }
    let mut diameter = (18.0 * viscosity * required / (delta * G)).sqrt();
    for _ in 0..MAX_ITERATIONS {
        let velocity = terminal_velocity(diameter, continuous, dispersed, viscosity).abs();
        if velocity <= 0.0 {
            break;
        }
        let ratio = required / velocity;
        let updated = diameter * ratio.sqrt();
        if (updated - diameter).abs() / (diameter + 1.0e-30) < TOLERANCE {
            diameter = updated;
            break;
        }
        diameter = 0.5 * (diameter + updated);
    }
    diameter
}

/// `estimateLiquidResidenceTime`.
fn estimate_liquid_residence_time(
    _diameter: f64,
    length: f64,
    _orientation: Orientation,
    liquid_level_fraction: f64,
    gas_velocity: f64,
) -> f64 {
    let gas_fraction = (1.0 - liquid_level_fraction).max(MIN_PHASE_FRACTION);
    let liquid_fraction = liquid_level_fraction.max(MIN_PHASE_FRACTION);
    let liquid_velocity = gas_velocity * gas_fraction / liquid_fraction;
    length / liquid_velocity.max(MIN_LIQUID_VELOCITY)
}

/// The vessel's orientation, `Separator`'s `orientation` as the class reads it.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Orientation {
    /// Horizontal, with the gas cap above the liquid.
    Horizontal,
    /// Vertical, with the gas above the liquid surface.
    Vertical,
}

/// The class's `calculate` arguments plus its configured distributions and curves.
#[derive(Debug, Clone, PartialEq)]
pub struct PerformanceSettings {
    /// Gas density, kg/m3.
    pub gas_density: f64,
    /// Oil density, kg/m3; zero is absent.
    pub oil_density: f64,
    /// Water density, kg/m3; zero is absent.
    pub water_density: f64,
    /// Gas viscosity, Pa s.
    pub gas_viscosity: f64,
    /// Oil viscosity, Pa s.
    pub oil_viscosity: f64,
    /// Water viscosity, Pa s.
    pub water_viscosity: f64,
    /// Superficial gas velocity, m/s.
    pub gas_velocity: f64,
    /// Vessel internal diameter, m.
    pub vessel_diameter: f64,
    /// Vessel length, m.
    pub vessel_length: f64,
    /// Vessel orientation.
    pub orientation: Orientation,
    /// Liquid level fraction.
    pub liquid_level_fraction: f64,
    /// The entrained liquid's distribution; absent leaves the gas-liquid branch unrun.
    pub gas_liquid: Option<DropletDistribution>,
    /// The entrained gas's distribution.
    pub gas_bubble: Option<DropletDistribution>,
    /// The dispersed oil's distribution in water.
    pub oil_in_water: Option<DropletDistribution>,
    /// The dispersed water's distribution in oil.
    pub water_in_oil: Option<DropletDistribution>,
    /// The gas outlet's internals.
    pub mist_eliminator: Option<CurveSettings>,
    /// The liquid-liquid coalescer.
    pub coalescer: Option<CurveSettings>,
    /// Whether the gravity stage is part of the gas-liquid efficiency.
    pub include_gravity_section: bool,
    /// The multiplier on `oil_in_gas` and `water_in_gas`.
    pub liquid_in_gas_factor: f64,
    /// The multiplier on `gas_in_oil` and `gas_in_water`.
    pub gas_carry_under_factor: f64,
    /// The multiplier on `oil_in_water` and `water_in_oil`.
    pub liquid_liquid_factor: f64,
}

/// The six fractions and the five efficiencies the standard path reports.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PerformanceOutcome {
    /// Oil carried out with the gas.
    pub oil_in_gas: f64,
    /// Water carried out with the gas.
    pub water_in_gas: f64,
    /// Gas carried under into the oil.
    pub gas_in_oil: f64,
    /// Gas carried under into the water.
    pub gas_in_water: f64,
    /// Oil left in the water.
    pub oil_in_water: f64,
    /// Water left in the oil.
    pub water_in_oil: f64,
    /// The gravity stage's integrated efficiency.
    pub gravity_section_efficiency: f64,
    /// The mist eliminator's integrated efficiency.
    pub mist_eliminator_efficiency: f64,
    /// The two stages in series.
    pub overall_gas_liquid_efficiency: f64,
    /// The oil-in-water gravity efficiency.
    pub liquid_liquid_gravity_efficiency: f64,
    /// The gas section's critical diameter, m.
    pub gravity_cut_diameter: f64,
}

impl Default for PerformanceOutcome {
    fn default() -> Self {
        Self {
            oil_in_gas: 0.0,
            water_in_gas: 0.0,
            gas_in_oil: 0.0,
            gas_in_water: 0.0,
            oil_in_water: 0.0,
            water_in_oil: 0.0,
            gravity_section_efficiency: 0.0,
            mist_eliminator_efficiency: 0.0,
            overall_gas_liquid_efficiency: 0.0,
            liquid_liquid_gravity_efficiency: 0.0,
            gravity_cut_diameter: 0.0,
        }
    }
}

/// The class's standard path.
///
/// # Errors
/// [`AzothError::InvalidInput`] for a curve missing what its type requires, or a lookup table
/// on a type that does not read one. The declared bounds keep the class's degenerate cut
/// (`f64::MAX`) unreachable.
pub fn separator_performance(settings: &PerformanceSettings) -> Result<PerformanceOutcome> {
    let mist = settings
        .mist_eliminator
        .as_ref()
        .map(|curve| GradeCurve::resolve("mist_eliminator", curve))
        .transpose()?;
    let coalescer = settings
        .coalescer
        .as_ref()
        .map(|curve| GradeCurve::resolve("oil_water_coalescer", curve))
        .transpose()?;
    let mut outcome = PerformanceOutcome::default();
    let has_liquid = settings.oil_density > 0.0 || settings.water_density > 0.0;
    let liquid_density = if settings.oil_density > 0.0 {
        settings.oil_density
    } else {
        settings.water_density
    };
    let liquid_viscosity = if settings.oil_density > 0.0 {
        settings.oil_viscosity
    } else {
        settings.water_viscosity
    };

    if let Some(distribution) = &settings.gas_liquid {
        if has_liquid {
            let (height, residence) = gas_section(settings);
            outcome.gravity_cut_diameter = critical_diameter(
                height,
                residence,
                settings.gas_density,
                liquid_density,
                settings.gas_viscosity,
            );
            let gravity = GradeCurve::gravity(outcome.gravity_cut_diameter);
            if settings.include_gravity_section {
                outcome.gravity_section_efficiency = gravity.overall(distribution);
            }
            if let Some(mist) = &mist {
                outcome.mist_eliminator_efficiency = mist.overall(distribution);
            }
            outcome.overall_gas_liquid_efficiency = match (settings.include_gravity_section, &mist)
            {
                (true, Some(mist)) => combined(distribution, &gravity, mist),
                (true, None) => outcome.gravity_section_efficiency,
                (false, Some(_)) => outcome.mist_eliminator_efficiency,
                (false, None) => 0.0,
            };
            let entrainment = 1.0 - outcome.overall_gas_liquid_efficiency;
            if settings.oil_density > 0.0 && settings.water_density > 0.0 {
                outcome.oil_in_gas = entrainment;
                outcome.water_in_gas = entrainment;
            } else if settings.oil_density > 0.0 {
                outcome.oil_in_gas = entrainment;
            } else {
                outcome.water_in_gas = entrainment;
            }
        }
    }

    if let Some(distribution) = &settings.gas_bubble {
        if has_liquid {
            let height = liquid_height(settings);
            let residence = if height <= 0.0 {
                0.0
            } else {
                estimate_liquid_residence_time(
                    settings.vessel_diameter,
                    settings.vessel_length,
                    settings.orientation,
                    settings.liquid_level_fraction,
                    settings.gas_velocity,
                )
            };
            let cut = critical_diameter(
                height,
                residence,
                liquid_density,
                settings.gas_density,
                liquid_viscosity,
            );
            let gravity = GradeCurve::gravity(cut);
            let efficiency = gravity.overall(distribution);
            outcome.gas_in_oil = 1.0 - efficiency;
            outcome.gas_in_water = outcome.gas_in_oil;
        }
    }

    if settings.oil_density > 0.0 && settings.water_density > 0.0 {
        let settling_height = liquid_height(settings) / 2.0;
        if let Some(distribution) = &settings.water_in_oil {
            let cut = critical_diameter(
                settling_height,
                LIQUID_LIQUID_RESIDENCE_TIME,
                settings.oil_density,
                settings.water_density,
                settings.oil_viscosity,
            );
            let gravity = GradeCurve::gravity(cut);
            let efficiency = gravity.overall(distribution);
            outcome.water_in_oil = match &coalescer {
                Some(coalescer) => (1.0 - efficiency) * (1.0 - coalescer.overall(distribution)),
                None => 1.0 - efficiency,
            };
        }
        if let Some(distribution) = &settings.oil_in_water {
            let cut = critical_diameter(
                settling_height,
                LIQUID_LIQUID_RESIDENCE_TIME,
                settings.water_density,
                settings.oil_density,
                settings.water_viscosity,
            );
            let gravity = GradeCurve::gravity(cut);
            let efficiency = gravity.overall(distribution);
            outcome.liquid_liquid_gravity_efficiency = efficiency;
            outcome.oil_in_water = match &coalescer {
                Some(coalescer) => (1.0 - efficiency) * (1.0 - coalescer.overall(distribution)),
                None => 1.0 - efficiency,
            };
        }
    }

    outcome.oil_in_gas = clamp01(outcome.oil_in_gas * settings.liquid_in_gas_factor);
    outcome.water_in_gas = clamp01(outcome.water_in_gas * settings.liquid_in_gas_factor);
    outcome.gas_in_oil = clamp01(outcome.gas_in_oil * settings.gas_carry_under_factor);
    outcome.gas_in_water = clamp01(outcome.gas_in_water * settings.gas_carry_under_factor);
    outcome.oil_in_water = clamp01(outcome.oil_in_water * settings.liquid_liquid_factor);
    outcome.water_in_oil = clamp01(outcome.water_in_oil * settings.liquid_liquid_factor);
    Ok(outcome)
}

/// The gas section's available height and residence time.
fn gas_section(settings: &PerformanceSettings) -> (f64, f64) {
    let gas_fraction = 1.0 - settings.liquid_level_fraction;
    let height = match settings.orientation {
        Orientation::Vertical => settings.vessel_length * gas_fraction,
        Orientation::Horizontal => settings.vessel_diameter * gas_fraction,
    };
    let residence = if gas_fraction > 0.0 && settings.gas_velocity > 0.0 {
        settings.vessel_length / settings.gas_velocity
    } else {
        0.0
    };
    (height, residence)
}

/// The liquid's height: the full length vertical, the diameter horizontal.
fn liquid_height(settings: &PerformanceSettings) -> f64 {
    match settings.orientation {
        Orientation::Vertical => settings.vessel_length * settings.liquid_level_fraction,
        Orientation::Horizontal => settings.vessel_diameter * settings.liquid_level_fraction,
    }
}

/// `clamp01`, the class's own.
fn clamp01(value: f64) -> f64 {
    value.clamp(0.0, 1.0)
}

/// The spec's spelling for one internals type.
fn name(internals: Internals) -> &'static str {
    match internals {
        Internals::Gravity => "gravity",
        Internals::WireMesh => "wire_mesh",
        Internals::VanePack => "vane_pack",
        Internals::AxialCyclone => "axial_cyclone",
        Internals::PlatePack => "plate_pack",
        Internals::Custom => "custom",
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn distribution(cut: f64) -> DropletDistribution {
        DropletDistribution {
            kind: DistributionKind::RosinRammler,
            characteristic_diameter: cut,
            spread: 2.6,
        }
    }

    #[test]
    fn the_class_default_wire_mesh_captures_a_100_micron_distribution() {
        let mut settings = base();
        settings.gas_liquid = Some(distribution(100.0e-6));
        settings.mist_eliminator = Some(CurveSettings {
            internals: Internals::WireMesh,
            cut_diameter: None,
            sharpness: None,
            max_efficiency: None,
            points: None,
        });
        let outcome = separator_performance(&settings).expect("the default curve");
        assert!(outcome.overall_gas_liquid_efficiency > 0.99);
        assert!(outcome.oil_in_gas < 0.01);
    }

    #[test]
    fn a_three_phase_state_writes_the_whole_entrainment_to_both_liquids() {
        let mut settings = base();
        settings.water_density = 1000.0;
        settings.water_viscosity = 1.0e-3;
        settings.gas_liquid = Some(distribution(100.0e-6));
        settings.mist_eliminator = Some(CurveSettings {
            internals: Internals::WireMesh,
            cut_diameter: None,
            sharpness: None,
            max_efficiency: None,
            points: None,
        });
        let outcome = separator_performance(&settings).expect("the default curve");
        assert_eq!(outcome.oil_in_gas, outcome.water_in_gas);
    }

    #[test]
    fn a_custom_curve_needs_its_table() {
        let mut settings = base();
        settings.mist_eliminator = Some(CurveSettings {
            internals: Internals::Custom,
            cut_diameter: None,
            sharpness: None,
            max_efficiency: None,
            points: None,
        });
        assert!(separator_performance(&settings).is_err());
    }

    fn base() -> PerformanceSettings {
        PerformanceSettings {
            gas_density: 50.0,
            oil_density: 800.0,
            water_density: 0.0,
            gas_viscosity: 1.5e-5,
            oil_viscosity: 2.0e-3,
            water_viscosity: 0.0,
            gas_velocity: 1.0,
            vessel_diameter: 2.0,
            vessel_length: 8.0,
            orientation: Orientation::Horizontal,
            liquid_level_fraction: 0.5,
            gas_liquid: Some(distribution(100.0e-6)),
            gas_bubble: None,
            oil_in_water: None,
            water_in_oil: None,
            mist_eliminator: None,
            coalescer: None,
            include_gravity_section: true,
            liquid_in_gas_factor: 1.0,
            gas_carry_under_factor: 1.0,
            liquid_liquid_factor: 1.0,
        }
    }
}
