//! `unit_ops.shortcut_distillation_column` - Fenske-Underwood-Gilliland on a flash's K-values.

use azoth_core::units::{Power, Pressure, watts};
use azoth_core::{AzothError, Result};
use azoth_eos::{Phase, pt_flash};

use crate::stream::Stream;

/// The closed-form column's answer: the two products and the eight scalars the class exposes.
///
/// A struct rather than a tuple because the scalars are the answer as much as the streams
/// are: `unit_ops.distillation_column`'s profile is a solved tray-by-tray state, and this
/// class's answer is seven numbers and a feed tray. Naming them keeps the model's declared
/// outputs a rearrangement of one value.
#[derive(Debug, Clone)]
pub struct ShortcutColumn {
    /// The distillate outlet.
    pub distillate: Stream,
    /// The bottoms outlet.
    pub bottoms: Stream,
    /// `alpha_LK/HK`, the K-value ratio the whole calculation is written in.
    pub relative_volatility: f64,
    /// Fenske's minimum stages at total reflux.
    pub minimum_stages: f64,
    /// Underwood's minimum reflux.
    pub minimum_reflux_ratio: f64,
    /// The stages Molokanov's fit gives at the actual reflux.
    pub actual_stages: f64,
    /// `minimum_reflux_ratio * reflux_ratio_multiplier`.
    pub actual_reflux_ratio: f64,
    /// The feed stage counted from the top.
    pub feed_tray_number: i64,
    /// The condenser duty the class reports.
    pub condenser_duty: Power,
    /// The reboiler duty the class reports.
    pub reboiler_duty: Power,
}

/// The class's K-value fallback for a single-phase feed.
///
/// `(Pc/P) exp(5.373 (1 + omega) (1 - Tc/T))`. **The 5.373 is the class's**, not the 5.37 of
/// Wilson's own estimate, and it is kept as written.
const WILSON_EXPONENT: f64 = 5.373;

/// The class's latent heat stand-in, J/mol, in `computeDuties`.
const AVERAGE_LATENT_HEAT: f64 = 30000.0;

/// The class's own split-fraction bounds, from `boundedSplitFraction`.
const SPLIT_FLOOR: f64 = 1.0e-12;

/// The class's own Underwood stopping rule.
const UNDERWOOD_TOLERANCE: f64 = 1.0e-10;
/// The second Underwood stopping rule: a bracket this narrow stops the bisection.
const UNDERWOOD_BRACKET: f64 = 1.0e-12;
/// The bisection's iteration cap, which the class's own loop uses.
const UNDERWOOD_ITERATIONS: usize = 200;

/// Solve a shortcut distillation column.
///
/// The arithmetic is `ShortcutDistillationColumn.run`'s, in its own order: flash the feed,
/// take `K_i = y_i/x_i` or the Wilson fallback, form `alpha_i = K_i / K_HK`, then Fenske,
/// Underwood, Gilliland through Molokanov, Kirkbride, the duties and the two products.
///
/// **Two of the class's answers are refused rather than reproduced.** `alpha_LK/HK <= 1` is
/// the class's own `solved = false`; a `reflux_ratio_multiplier` at or below one leaves
/// Gilliland at `Infinity` stages, or takes its silent `Y = 0.5` fallback. And a flash whose
/// K-values straddled one is [`Phase::Trivial`], where NeqSim reads a phase type this
/// library's model does not prove.
///
/// # Errors
/// [`azoth_core::AzothError::InvalidInput`] for a key that is not a component, a relative
/// volatility at or below one, a reflux multiplier at or below one, and the trivial-solution
/// case; whatever the flash refuses.
#[allow(clippy::too_many_arguments)] // one parameter per declared input
pub fn shortcut_distillation_column(
    feed: &Stream,
    light_key: &str,
    heavy_key: &str,
    light_key_recovery_distillate: f64,
    heavy_key_recovery_bottoms: f64,
    reflux_ratio_multiplier: f64,
    condenser_pressure: Option<Pressure>,
    reboiler_pressure: Option<Pressure>,
) -> Result<ShortcutColumn> {
    let lk = index_of(&feed.components, light_key, "light_key")?;
    let hk = index_of(&feed.components, heavy_key, "heavy_key")?;
    if reflux_ratio_multiplier <= 1.0 {
        return Err(AzothError::invalid_input(
            "reflux_ratio_multiplier",
            format!(
                "a reflux multiplier of {reflux_ratio_multiplier} is at or below the minimum: at one \
                 the class returns an infinite stage count and at less than one it takes its silent \
                 Y = 0.5 fallback, and neither is a column"
            ),
        ));
    }

    let (mixture, _ideal_gas) = feed.mixture()?;
    let n = feed.z.len();
    let flash = pt_flash(&mixture, feed.t, feed.p, &feed.z)?;

    // `ShortcutDistillationColumn` asks whether both a gas and a liquid phase are there and
    // divides one by the other; where either is absent it has no ratio and estimates one.
    let k_values: Vec<f64> = match flash.phase {
        Phase::TwoPhase => (0..n)
            .map(|i| {
                if flash.x[i] > 1.0e-20 {
                    flash.y[i] / flash.x[i]
                } else {
                    1.0e10
                }
            })
            .collect(),
        Phase::AllVapour | Phase::AllLiquid => {
            let components = mixture.components();
            (0..n)
                .map(|i| {
                    let c = &components[i];
                    (c.pc.value / feed.p.value)
                        * (WILSON_EXPONENT * (1.0 + c.omega) * (1.0 - c.tc.value / feed.t.value))
                            .exp()
                })
                .collect()
        }
        Phase::Trivial => {
            return Err(AzothError::invalid_input(
                "feed_t / feed_p",
                format!(
                    "the flash at {} K and {} Pa is a trivial solution - every K-value straddled \
                     one, so nothing here proves whether the feed is a gas or a liquid. The \
                     class reads its phase type, which decides both the Wilson fallback and the \
                     feed quality, and this model does not decide it",
                    feed.t.value, feed.p.value
                ),
            ));
        }
    };

    let k_hk = k_values[hk];
    let relative_volatility = k_values[lk] / k_hk;
    if relative_volatility <= 1.0 {
        return Err(AzothError::invalid_input(
            "light_key / heavy_key",
            format!(
                "{light_key} is not more volatile than {heavy_key}: alpha_LK/HK is \
                 {relative_volatility}, and the class logs an error and leaves every answer at \
                 its field initialiser"
            ),
        ));
    }
    let alpha: Vec<f64> = k_values.iter().map(|k| k / k_hk).collect();

    // The class's `computeFeedQuality`: one phase is 0 for a gas and 1 for anything else, two
    // phases are the liquid's share of the moles, which on a molar flash is `1 - beta`.
    let q = match flash.phase {
        Phase::TwoPhase => 1.0 - flash.beta.unwrap_or(0.0),
        Phase::AllVapour => 0.0,
        Phase::AllLiquid => 1.0,
        Phase::Trivial => unreachable!("refused above"),
    };

    // Fenske, on the recoveries clamped exactly as `boundedSplitFraction` clamps them.
    let x_lk_d = bounded(light_key_recovery_distillate);
    let x_hk_d = bounded(1.0 - heavy_key_recovery_bottoms);
    let x_lk_b = bounded(1.0 - light_key_recovery_distillate);
    let x_hk_b = bounded(heavy_key_recovery_bottoms);
    let minimum_stages = ((x_lk_d / x_hk_d) * (x_hk_b / x_lk_b)).ln() / relative_volatility.ln();

    // Underwood's root, between the heavy key's relative volatility and the light key's.
    let theta = underwood_root(&alpha, &feed.z, q, lk);

    let fractions = estimate_distillate_fractions(
        &alpha,
        lk,
        hk,
        relative_volatility,
        minimum_stages,
        light_key_recovery_distillate,
        heavy_key_recovery_bottoms,
    );

    // Underwood's second equation, on the distillate composition the split fractions imply.
    let mut x_d: Vec<f64> = (0..n).map(|i| feed.z[i] * fractions[i]).collect();
    let total_d: f64 = x_d.iter().sum();
    if total_d > 0.0 {
        for value in &mut x_d {
            *value /= total_d;
        }
    }
    let mut r_min_plus_one = 0.0;
    for i in 0..n {
        if x_d[i] > 1.0e-15 && (alpha[i] - theta).abs() > 1.0e-10 {
            r_min_plus_one += alpha[i] * x_d[i] / (alpha[i] - theta);
        }
    }
    let minimum_reflux_ratio = (r_min_plus_one - 1.0).max(0.0);

    // Gilliland through Molokanov's fit, and the stages it inverts to.
    let actual_reflux_ratio = minimum_reflux_ratio * reflux_ratio_multiplier;
    let x_gilliland = (actual_reflux_ratio - minimum_reflux_ratio) / (actual_reflux_ratio + 1.0);
    let y_gilliland = if (0.0..=1.0).contains(&x_gilliland) {
        1.0 - (((1.0 + 54.4 * x_gilliland) / (11.0 + 117.2 * x_gilliland))
            * ((x_gilliland - 1.0) / x_gilliland.sqrt()))
        .exp()
    } else {
        // The class's `Y = 0.5` fallback. Unreachable for a multiplier above one, and kept
        // because that is what the class does with an `X` outside the interval.
        0.5
    };
    let actual_stages = (y_gilliland + minimum_stages) / (1.0 - y_gilliland);

    // Kirkbride, on the class's own argument - see the model spec's assumptions for how it
    // differs from the correlation as usually quoted.
    let z_lk_f = feed.z[lk];
    let z_hk_f = feed.z[hk];
    let b_over_d = bottoms_to_distillate(&feed.z, &fractions);
    let x_hk_distillate = x_hk_d * z_hk_f;
    let x_lk_bottoms = x_lk_b * z_lk_f;
    let kirkbride_ratio = if x_lk_bottoms > 1.0e-15 && b_over_d > 1.0e-15 {
        ((z_hk_f / z_lk_f) * (x_lk_bottoms / x_hk_distillate) * (b_over_d * b_over_d)).powf(0.206)
    } else {
        0.0
    };
    let feed_tray_number = (actual_stages / (1.0 + kirkbride_ratio)).round() as i64 + 1;

    // The duties, which are estimates: a hard-coded latent heat and one per cent of the
    // feed's own enthalpy, in W because the feed's flow is mol/s.
    let total_distillate_fraction = (0..n)
        .map(|i| feed.z[i] * fractions[i])
        .sum::<f64>()
        .clamp(0.001, 0.999);
    let distillate_flow = feed.n * total_distillate_fraction;
    let vapour_flow = distillate_flow * (actual_reflux_ratio + 1.0);
    let condenser_duty = -vapour_flow * AVERAGE_LATENT_HEAT;
    let reboiler_duty = -condenser_duty + (feed.n * feed.h.value) * 0.01;

    // The products: the feed's fluid with the split's moles at the stated pressure, at the
    // feed's temperature.
    let mut distillate_moles = vec![0.0; n];
    let mut bottoms_moles = vec![0.0; n];
    for i in 0..n {
        let total = feed.n * feed.z[i];
        distillate_moles[i] = (total * fractions[i]).max(0.0);
        bottoms_moles[i] = (total - distillate_moles[i]).max(0.0);
    }
    let distillate = product(feed, &distillate_moles, condenser_pressure)?;
    let bottoms = product(feed, &bottoms_moles, reboiler_pressure)?;

    Ok(ShortcutColumn {
        distillate,
        bottoms,
        relative_volatility,
        minimum_stages,
        minimum_reflux_ratio,
        actual_stages,
        actual_reflux_ratio,
        feed_tray_number,
        condenser_duty: watts(condenser_duty),
        reboiler_duty: watts(reboiler_duty),
    })
}

/// One product stream: the feed's fluid at the split's composition and the stated pressure.
///
/// The class sets the pressure only when the stated one is positive, so an unstated or
/// non-positive pressure leaves the product at the feed's.
fn product(feed: &Stream, moles: &[f64], stated: Option<Pressure>) -> Result<Stream> {
    let total: f64 = moles.iter().sum();
    if total <= 0.0 {
        return Err(AzothError::invalid_input(
            "light_key_recovery_distillate / heavy_key_recovery_bottoms",
            "this split leaves one product empty, and a stream with no moles has no \
             composition",
        ));
    }
    let z: Vec<f64> = moles.iter().map(|m| m / total).collect();
    let p = match stated {
        Some(p) if p.value > 0.0 => p,
        _ => feed.p,
    };
    Stream::from_pt(feed.components.clone(), z, total, p, feed.t)
}

/// The class's `boundedSplitFraction`: a finite fraction pulled inside `[1e-12, 1 - 1e-12]`,
/// and a non-finite one replaced by `0.5`.
fn bounded(split: f64) -> f64 {
    if !split.is_finite() {
        return 0.5;
    }
    split.clamp(SPLIT_FLOOR, 1.0 - SPLIT_FLOOR)
}

/// The fraction of each feed component the class sends to the distillate.
///
/// The two keys take their recoveries; anything more volatile than the light key takes
/// `0.999` and anything less volatile than the heavy key `0.001`; and a component between the
/// keys takes a Fenske-like expression in `n_min`. Every one is then bounded.
fn estimate_distillate_fractions(
    alpha: &[f64],
    lk: usize,
    hk: usize,
    relative_volatility: f64,
    minimum_stages: f64,
    light_key_recovery_distillate: f64,
    heavy_key_recovery_bottoms: f64,
) -> Vec<f64> {
    (0..alpha.len())
        .map(|i| {
            let fraction = if i == lk {
                light_key_recovery_distillate
            } else if i == hk {
                1.0 - heavy_key_recovery_bottoms
            } else if alpha[i] > relative_volatility {
                0.999
            } else if alpha[i] < 1.0 {
                0.001
            } else {
                1.0 / (1.0
                    + (alpha[i] / relative_volatility).powf(-minimum_stages)
                        * (1.0 - light_key_recovery_distillate)
                        / light_key_recovery_distillate)
            };
            bounded(fraction)
        })
        .collect()
}

/// The class's `computeBottomsToDistillateRatio`, on the same split fractions.
fn bottoms_to_distillate(z_feed: &[f64], fractions: &[f64]) -> f64 {
    let total_d: f64 = (0..z_feed.len()).map(|i| z_feed[i] * fractions[i]).sum();
    let total_b: f64 = (0..z_feed.len())
        .map(|i| z_feed[i] * (1.0 - fractions[i]))
        .sum();
    if total_d > 1.0e-15 {
        total_b / total_d
    } else {
        1.0
    }
}

/// The class's `solveUnderwood`: bisection on `[1 + 1e-6, alpha_LK/HK - 1e-6]`.
///
/// The bracket's upper end is the **light key's** relative volatility and not the array's
/// maximum, which is a different number whenever a component is more volatile than the light
/// key - and in the captured four-component row two of them are.
fn underwood_root(alpha: &[f64], z_feed: &[f64], q: f64, light_key: usize) -> f64 {
    let mut low = 1.0 + 1.0e-6;
    let mut high = alpha[light_key] - 1.0e-6;
    let mut theta = 0.5 * (low + high);
    for _ in 0..UNDERWOOD_ITERATIONS {
        let mid = 0.5 * (low + high);
        let f_mid = underwood_function(alpha, z_feed, mid, q);
        if f_mid.abs() < UNDERWOOD_TOLERANCE || (high - low) < UNDERWOOD_BRACKET {
            theta = mid;
            break;
        }
        let f_low = underwood_function(alpha, z_feed, low, q);
        if f_low * f_mid < 0.0 {
            high = mid;
        } else {
            low = mid;
        }
        theta = mid;
    }
    theta
}

/// The class's `underwoodFunction`: `sum_i alpha_i z_i / (alpha_i - theta) - (1 - q)`.
fn underwood_function(alpha: &[f64], z_feed: &[f64], theta: f64, q: f64) -> f64 {
    let mut sum = 0.0;
    for i in 0..alpha.len() {
        if (alpha[i] - theta).abs() > 1.0e-10 {
            sum += alpha[i] * z_feed[i] / (alpha[i] - theta);
        }
    }
    sum - (1.0 - q)
}

/// The position of a named component, or a refusal naming what the caller asked for.
fn index_of(components: &[String], name: &str, which: &str) -> Result<usize> {
    components.iter().position(|c| c == name).ok_or_else(|| {
        AzothError::invalid_input(
            which,
            format!("{name} is not one of this fluid's components"),
        )
    })
}
